#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <dirent.h>

#include "../lib/RTKLIB/src/rtklib.h"
#include "satinsmap.h"
//#include "INS_GNSS.h"

/* constants -----------------------------------------------------------------*/
#define MAXVEL      0.1           /* max velocity for using non-holonomic constraint */
#define MAXGYRO     (10.0*D2R)    /* max rotation speed for using non-holonomic constraint */
#define VARVEL      SQR(0.05)     /* initial variance of receiver vel ((m/s)^2) */
#define MINZC       15            /* min count for zero velocity update once */


/* global variables ----------------------------------------------------------*/
FILE *fp_lane;       /* Lane coordinate file pointer */
FILE *fimu;          /* Imu datafile pointer  */
FILE *imu_tactical; /* Imu datafile pointer */
lane_t lane;
imuraw_t imu_obs_global={{0}}; 
pva_t pva_global={{0}};  
pva_t pvagnss={{0}};
insgnss_opt_t insgnssopt={0};
sol_t *solw;         /* gnss solution structure buffer */
obsb_t obsw={0};          /* observation data buffer */
ins_states_t *insw;   /* ins states buffer */
int zvu_counter = 0;
int gnss_w_counter = 0;
int ins_w_counter = 0;
const double Omge[9]={0,OMGE,0,-OMGE,0,0,0,0,0}; /* (5.18) */

char *outpath1[] = {"../out/"};  
FILE *out_PVA;
FILE *out_clock_file;
FILE *out_IMU_bias_file;
FILE *out_tropo_file;
FILE *out_amb_file;
FILE *out_KF_SD_file;
FILE *out_KF_state_error;
FILE *out_raw_fimu;
FILE *out_KF_residuals; 

/* global states index -------------------------------------------------------*/
int IA = 0, NA = 0;   /* index and number of attitude states */
int IV = 0, NV = 0;   /* index and number of velocity states */
int IP = 0, NP = 0;   /* index and number of position states */
int iba = 0, nba = 0; /* index and number of accl bias states */
int ibg = 0, nbg = 0; /* index and number of gyro bias states */
int irc = 0, nrc = 0; /* index and number of receiver clock state */
int irr = 0, nrr = 0; /* index and number of receiver clock drift state */
int IT = 0, NT = 0;   /* index and number of tropo state */
int IN = 0, NN = 0;   /* index and number of ambiguities state */

/* extract substrings from string ---------------------------------------------
* extract unsigned/signed bits from input string and modify and return res
* args   : unsigned char *buff I byte data
*          int    pos    I      position from start of data
*          int    len    I      length
* return : res extracted string
*-----------------------------------------------------------------------------*/
extern void getstr1(const char *input, int offset, int len, char *res)
{
        strncpy(res, input + offset, len); 
}

/* Summation of matrix diagonal ------------------------------------------------
* compute summation of the square root of the diagonal elements of a squared matrix
* args   : double *A       I   matrix (nxn)
*          int  n       I   number of rows and columns of matrix A
           int ini, end  I from (ini) index to (end) index in matrix A
* return : summation
*-----------------------------------------------------------------------------*/
extern double Sumdiag(const double *A, int n, int ini, int end) 
{
    double sum=0.0;
    int i; 

    for (i=ini;i<end;i++) sum+=sqrt(A[i+i*n]);

    return sum;
}

/* 3D distance ----------------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
* return : geometric distance (m)
*-----------------------------------------------------------------------------*/
extern double dist(const double *rs, const double *rr)
{
    double r,e[3];
    int i;

    for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
    return r=norm(e,3);
}

/* Minimum value --------------------------------------------------------------
* find the minimum value of a vector
* args   : double a	I	array containing values
*	 : int	n	I	size of array
*        : double min	O	minimum value
* return : position of minimum value within the array
*-----------------------------------------------------------------------------*/
extern int min(double *a, int n, double *min) {
  int c, index;

  *min = a[0];
  index = 0;

  for (c=1;c<n;c++) {
     if (a[c]<*min){
       index = c;
       *min = a[c];
     }
  }
  return index;
}

/* Minimum resdiuals value --------------------------------------------------------------
* find the minimum residuals (close to zero) value of a vector
* args   : double a	I	array containing values
*	 : int	n	I	size of array
*        : double min	O	minimum value
* return : position of minimum value within the array
*-----------------------------------------------------------------------------*/
extern int resmin(double *a, int n, double *min) {
  int c, index;

  *min = a[0];
  index = 0;

  for (c=1;c<n;c++) {
     if (abs(a[c])<abs(*min)){
       index = c;
       *min = a[c];
     }
  }
  return index;
}

/* Geodetic to Local topocentric & Local to Geodetic system	--------------------
% DG2LG  Converts lat,lon,h to local geodetic coordinates.
%   Local origin at lat,lon.  If astronomic lat,h input,
%   then output is in local astronomic system.  Vectorized.
% Version: 2017-10-11
% Usage:   [dx,dy]=dg2lg(lat,h,lat,lon)
% Input:   pos - vector of latitude and longitude differences (rad)
%          lat0  - vector of lats of local system origins (rad)
%          h0    - vector of hts of local system origins
% Output:  e(1-x)   -  vector of y (E) coordinates in local system
%          e(2-y)   -  vector of x (N) coordinates in local system

% Adapted from Copyright (c) 2011, Michael R. Craymer
*-----------------------------------------------------------------------------*/
extern d2lgs(double lat0, double h0, double* pos, double* e){
double v, r, dlat, dlon;
double a,e2;

dlat=pos[0];
dlon=pos[1];
a=RE_WGS84;
//e2=FE_WGS84;

v=a/sqrt(1-e2*pow(sin(lat0), 2) );
r=v*(1-e2)/(1-e2*(pow(sin(lat0),2)));
e[1]=dlat*(r+h0); //x(N) local system
e[0]=dlon*cos(lat0)*(v+h0);//y(E) local system
}

/* Azimuth between two points *************************************************/
double azt(double *enu1, double *enu2){
  double DE,DN,delta;
  double azm,q; //azimuth and its quadrant

  DE=enu2[0]-enu1[0]; DN=enu2[1]-enu1[1];

  delta = ( DE / DN );
  //                                angle in radians
  if (delta<0){ // 2nd or 4th Q
      if ( DN < 0){
          // 2nd quadrant
          azm = 180*(PI/180) + atan(delta); //adding because the value is negative
          q=2;
        }else{
          // 4th quadrant
          azm = 360*(PI/180) + atan(delta);//adding because the value is negative
          q=4;
        }
   }else if ( DE < 0){   // 1st or 3rd Q
          // 3rd quadrant
          azm = 180*(PI/180) + atan(delta);
          q=3;
        }
      else{ // 1st quadrant
          azm = atan(delta);
          q=1;
      }

  return azm; //in radians
}

/* Lane lever-arm correction ---------------------------------------------------
* Account the lever-arm into the lanes according to the navigaiton direction
* args:	IO double *r[3]	vector in xyz with the lane point [m]
*       I double *recbody[3]  receiver body-frame offsets {x,y,z(height)}
*-----------------------------------------------------------------------------*/
void antlevarm(double *enu, double az, double *recbody){
  double de,dn;
  double ang,c;
  /* correction factor in body frame xy resultant*/
  c=sqrt(recbody[0]*recbody[0]+recbody[1]*recbody[1]);

  // *** ESPECIAL CASE! antenna at the right side of platform origin:
  /* Consider when az=0=360, 90, 180, 270 */
  if (0.0<az && az<(90*D2R)) {
    ang=((90*D2R)-az);
    de=sin(ang)*c;
    dn=cos(ang)*c;
    //printf("I\nang:%lf\n",ang*R2D);
  //  printf("de:%lf\n", de);
  //  printf("dn:%lf\n", dn);
    enu[0]=enu[0]+de;
    enu[1]=enu[1]-dn;
  }else if ((90*D2R)<az&&az<(180*D2R)) {
    ang=(az-PI/2);
    de=sin(ang)*c;
    dn=cos(ang)*c;
  //  printf("II\nang:%lf\n",ang*R2D);
  //  printf("de:%lf\n", de);
  //  printf("dn:%lf\n", dn);
    enu[0]=enu[0]-de;
    enu[1]=enu[1]-dn;
  }else if ((180*D2R)<az && az<(270*D2R)) {
    ang=((270*D2R)-az);
    de=sin(ang)*c;
    dn=cos(ang)*c;
  //  printf("III\nang:%lf\n",ang*R2D);
  //  printf("de:%lf\n", de);
  //  printf("dn:%lf\n", dn);
    enu[0]=enu[0]-de;
    enu[1]=enu[1]+dn;
  }else{
    ang=az-(270*D2R);
    de=sin(ang)*c;
    dn=cos(ang)*c;
  //  printf("IV\nang:%lf\n",ang*R2D);
  //  printf("de:%lf\n", de);
  //  printf("dn:%lf\n", dn);
    enu[0]=enu[0]+de;
    enu[1]=enu[1]+dn;
  }

  // Adding lever-arm to the lanes (which is on the ground)
  enu[2]=enu[2]+recbody[2];
}

/* 3D distance interpolation ---------------------------------------------------
* Creates points between two other points
* args   : FILE *fp       I   file containing the points to be interpolated {xyz}
*          double  d      I   space beteween points in meters
*         int n           I   number of points
* return : file output containing the list of itnerpolated points
*-----------------------------------------------------------------------------*/
void lineInterp(FILE *fp, double d, int n)
{
  double D[3], p[3*n], p1[3], p2[3], pn[3], nrm=0.0, daux=0.0;
  char str[50], fileout[150];
  int  i,j;
  double pos0[3],pos1[3],pos2[3],posn[3],enu1[3],enu2[3],enun[3],enu[3],pos[3];
  //pos0[0]=45.94312323*D2R; pos0[1]=-66.64106835*D2R; pos0[2]=52.7309; //BMO corner 1
  //pos0[0]=45.9431226252*D2R; pos0[1]=-66.6410698327*D2R; pos0[2]=54.5210901409; //BMO corner 1 from #5 with kinematic height
  pos0[0]=45.9431226252*D2R; pos0[1]=-66.6410698327*D2R; pos0[2]=54.4000901409; //BMO corner 1 from #5 with static height
  double r0[3],r[3],dr[3],dpos[3];
  double recbody[3]; /* Antena offsets in body-frame {x(nav.direc),y(right),z(down)} in meters */
  double az;
  //recbody[0]=0;recbody[1]=0.113;recbody[2]=1.629; /* exp4*/
  //recbody[0]=0;recbody[1]=0.113;recbody[2]=1.595; /* exp5 static height*/
  recbody[0]=0;recbody[1]=0.113;recbody[2]=1.716; /* exp5 kin. height*/
  FILE *fpout;

  pos2ecef(pos0, r0);

  strcat(fileout,outpath1[0]);
  strcat(fileout,"LineInterp.txt");
  fpout=fopen("../out/LineInterp_rtklib5_static_height_rtkplot.txt","w");

  for (i = 0; i < n; i++) {
    fgets(str, 50, fp);
    puts(str);
    sscanf(str, "%lf %lf %lf", &p[i*3], &p[i*3+1], &p[i*3+2]);
  }

  for (i = 0; i < n; i++) {
    if(i==n-1){
      for (j=0; j<3; j++) {p1[j]=p[i*3+j];p2[j]=p[(0)*3+j];}
    }else{
      for (j=0; j<3; j++) {p1[j]=p[i*3+j];p2[j]=p[(i+1)*3+j];}
     }

    for (j = 0; j < 3; j++) D[j]=p2[j]-p1[j];
    nrm=norm(D,3);

    for (j = 0; j < 3; j++) dr[j]=p1[j]-r0[j];
    ecef2enu(pos0, dr, enu1);

    for (j = 0; j < 3; j++) dr[j]=p2[j]-r0[j];
    ecef2enu(pos0, dr, enu2);
    az=azt(enu1, enu2);

    printf("\n Azimuth1: %lf\n",az*R2D);

    /* Antenna lever arm correction - apply it on the lanes instead
    antlevarm(enu1,az,recbody);
    enu2ecef(pos0,enu1,r);
    for (j = 0; j < 3; j++) p1[j]=r[j]+r0[j];

    antlevarm(enu2,az,recbody);
    enu2ecef(pos0,enu2,r);
    for (j = 0; j < 3; j++) p2[j]=r[j]+r0[j];   */

    fprintf(fpout, "%11.5lf %11.5lf %11.5lf\n",p1[0],p1[1],p1[2]);

    daux=0.1;
    while(daux<=nrm){
    for (j=0; j < 3; j++) pn[j]=p1[j]+daux*(D[j]/nrm);
    fprintf(fpout, "%11.5lf %11.5lf %11.5lf\n",pn[0],pn[1],pn[2]);
    daux+=d;
	 }
   fprintf(fpout, "%11.5lf %11.5lf %11.5lf\n",p2[0],p2[1],p2[2]);
  fprintf(fpout, "\n");
 }
   printf("Finished 3D line interpolation\n");
   fclose(fpout);
}

/* C Program to List Files in Directory from: http://www.sanfoundry.com/c-program-list-files-directory/
 * It needs the dirent.h library ---------------------------------------------*/
extern int listfiles(void)
{
    DIR *d;
    struct dirent *dir;
    char* dirpath = "/media/emerson/My_files/IAG_IASPEI_2017/IAG_IASPEI_2017/tropo/";
    char* cpychar[70], buff[150];
    d = opendir(dirpath);
    FILE* fp;
    FILE* res;
    int i, j;
    res = fopen("/media/emerson/My_files/IAG_IASPEI_2017/IAG_IASPEI_2017/troporesults.txt", "w");
 i=0;
 if (d){
        while ((dir = readdir(d)) != NULL){ i++; }
 }
 i=i-2; /* Accounting for . and .. */
 printf("Number of files: %d\n", i);

 rewinddir(d);

    if (d)
    {

        while ((dir = readdir(d)) != NULL)
        {
	    j=0;
            printf("%s\n", dir->d_name);
 	    strcpy(cpychar,dirpath);
	    strcat(cpychar,dir->d_name);
	    //printf("Path: %s \n",cpychar);
	    fp=fopen(cpychar, "r");
	    if(fp==NULL){printf("Error in file opening: %s\n", dir->d_name);}
	    /* Send*/
	    while ( fgets(buff, 150, fp)!= NULL ){
		   j++;
		}
 	    //printf("Size of file: %d\n",j);
	    //puts(buff);
            fprintf(res,"%s ",dir->d_name);
            fprintf(res, buff);

        }
        memset(cpychar,0,sizeof(cpychar));
        closedir(d);
        fclose(fp);
    }
    fclose(res);
    return(0);
}

/* Beginning of rnx2rtkp -----------------------------------------------------*/
#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     8                   /* max number of input files */

static const char rcsid[]="$Id: rnx2rtkp.c,v 1.1 2008/07/17 21:55:16 ttaka Exp $";

/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:moving-base,",
"                 5:fixed,6:ppp-kinematic,7:ppp-static) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
" -y level  output soltion status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]"
};


/* show message --------------------------------------------------------------*/
extern int showmsg(char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}

static void dumpobs(obs_t *obs)
{
    gtime_t time={0};
    int i;
    char str[64];
    printf("obs : n=%d\n",obs->n);
    for (i=0;i<2;i++) { //obs->n
        time2str(obs->data[i].time,str,3);
        printf("%s : %2d %2d %13.3f %13.3f %13.3f %13.3f  %d %d\n",str,obs->data[i].sat,
               obs->data[i].rcv,obs->data[i].L[0],obs->data[i].L[1],
               obs->data[i].P[0],obs->data[i].P[1],obs->data[i].LLI[0],obs->data[i].LLI[1]);

        assert(1<=obs->data[i].sat&&obs->data[i].sat<=32);
        assert(timediff(obs->data[i].time,time)>=-DTTOL);

        time=obs->data[i].time;
    }
}

/* Return a skew matrix (W -3x3) af a vector (vec - 3x1)  --------------------*/
extern void vec2skew (double *vec, double *W){
  W[0]= 0.0;    W[1]=-vec[2]; W[2]= vec[1];
  W[3]= vec[2]; W[4]= 0.0;    W[5]=-vec[0];
  W[6]=-vec[1]; W[7]= vec[0]; W[8]= 0.0;
}

/* matl2c */
/* matl2c ----------------------------------------------------------------------
* description: converts a matrix A (l by c) into A (l by c) in column-major
order, or tranpose A
* GNSS solutions
* args :
*       double*   A_l  I  matrix in line-order (l by c)
*       int       l    I  number of lines
*       int       c    I  number of columns
*       double*   A_c  IO matrix converted into column-order (l by c)
*-----------------------------------------------------------------------------*/
extern void matl2c(double *A_l, int l, int c, double *A_c){
  int i,j;

  for (i=0;i<l;i++) {
    for (j=0;j<c;j++) {
      A_c[i*c+j]=A_l[j*l+i];
    }
  }

}

/* Check if GNSS and IMU measurements are synchronized
   return 1 if yes, 0 if is not, status=-1 if GNSS time behing IMU ----------------*/
int gnssimusync(float tgps, um7pack_t *imu){
 float timu;

 timu=time2gpst(imu->time,&imu->gpsw);
 timu=imu->sec;

 printf("TGPS: %f, TIMU: %f \nTrue GNSS and IMU diff: %lf\n",tgps, timu, fabs(tgps-timu));

  if (  fabs(tgps-timu)>0.41  ) {

    if (tgps-timu > 0) { /* gnss ahead */
      imu->status=0;
      imu->count=0;
      if (tgps-timu > 0.5) {
        printf("GNSS time ahead imu time!!\n");
        /* Read next imu to catch up, however it would skip one imu meas? */
      }
      return 0;
    }else{ /* gnss back */
      printf("GNSS time behind imu time!! INS stays in the same epoch\n");
      imufileback();
      imu->status=-1;
      imu->count=0;
      return 0;
    }
  }else{ /* synchronized */
    return 1;
  }

}

/* Loosley-coupled INS-(MAP or GNSS) algorithm */
extern void ins_LC (double* gnss_xyz_ini_pos, double* gnss_xyz_ini_cov, double* gnss_enu_vel, double ini_pos_time,
um7pack_t *imu, pva_t *pvap, imuraw_t *imuobsp )
{
  double *z,*H,*R,*xp,*Pp, *Q, *PHI,*G,*u, dt;
  double *xk, *x_, *Pk, *P_, *v, *v_, *Q_,*F;
  double *F1,*Q1,*K1,*I1,*xp1;
  double *P_aux,*Q_aux;
  double E[9],Ea[9],I[9], Cbn_aux[9],l[3]; /* Lever-arm from INS to GNSS antenna in b-frame {x(nav.direc),y(right),z(down)} in meters */
  int i,j,nx,n,nv, info,iter;

  /* Number of states (nx) and measureprintf("\n\n **********GNSS/INS **************************************\n\n" );ments (n)  */
  nx=15; /* ins 15-states {de,dv,dr,dba,dbg} system model in local-nav-frame */
  n=  6; /* Gnss-ins differences in position and velocity */
  nv=n;

  xp=mat(nx,1); Pp=zeros(nx,nx); Q=zeros(nx,nx);
  PHI=mat(nx,nx); G=zeros(nx,nx-3);  u=mat(nx-3,1);
  z=mat(nv,1); H=zeros(nv,nx); R=zeros(nv,nv);
  x_=mat(nx,1); xk=mat(nx,1);
  P_=zeros(nx,nx); Pk=zeros(nx,nx);
  v=mat(nv,1); v_=mat(nv,1);
  Q_=mat(nx,nx); F=mat(nx,nx);
  P_aux=mat(nx,nx); Q_aux=mat(nx,nx);
  /* For the Filter*/
  F1=mat(nx,n);Q1=mat(n,nx);K1=mat(nx,n);I1=eye(nx);xp1=mat(nx,1);

  /* Initialize the ins error-state and covariances xp, Pp  (Groves 2013, page 594)*/
  if (iter==0) {
    for (i=0;i<nx;i++) xp[i]=0.0;
  }
//  posvelP0(Pp,nx,rtk);
//  attbiasP0(&imuobsp, &pvap, Pp, nx);

  /* Generic intialization for Loosley coupled GNSS/INS  */
  initialize_P0(Pp,nx);

  //inssysmatrix(PHI, G, nx, &pvap, &imuobsp, dt);
  /* Transition matrix using  (14.73 - page 589) */
  inssysmatrix1(PHI, nx, &pvap, &imuobsp, dt);
  /* System noise */
  sysnoise1(Q, dt, nx);

  /* Kalman filter system propagation */

  /* Propagate state estimates using (3.14) noting that all states are zero
  due to closed-loop correction.*/
   for (i=0;i<nx;i++) x_[i]=0.0;
   //matmul("NN",nx,1,nx,1.0,PHI,xp,0.0,x_); /* x_= T.xp */

  /* Propagate state error covariance using (3.46, page 100) */

  for (i = 0; i < nx; i++) {
    for (j = 0; j < nx; j++) {
      P_aux[i*nx+j]= Pp[i*nx+j]+0.5*Q[i*nx+j]; /* ...*(Pp + 0.5*Q)*... */
    }
  }
  matmul("NT",nx,nx,nx,1.0,P_aux,PHI,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
  matmul("NN",nx,nx,nx,1.0,PHI,Q_aux,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

  for (i = 0; i < nx; i++) {
    for (j = 0; j < nv; j++) {
      P_[i*nx+j]= Q_[i*nx+j]+0.5*Q[i*nx+j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
    }
  }

  /* Decide whether to perform or not the measurement update  ? */


  /* Kalman filter or EKF measurement update */
  measvec(z, gnss_xyz_ini_pos, &pvap, &imuobsp, l);
  measmatrixH (H, nv, nx, &pvap);
  measnoiseR (R, nv, gnss_xyz_ini_cov);
  //measnoiseR (R, nv, rtk);

  /* Innovation term */
  matmul("NN",nv,1,nx,1.0,H,x_,0.0,v_);  /* v=z-H*x_ */
  for (i=0;i<nv;i++) v[i]=z[i]-v_[i];

  for (i = 0; i < nv; i++) {
  //   printf("v[%d]:%lf\n",i,z[i]);
  }

  /* Update  */
  /*  how about some iterations??  */
    matcpy(Q1,R,n,n);
    matcpy(xp1,x_,nx,1);
    matmul("NT",nx,n,nx,1.0,P_,H,0.0,F1);      /* Q=H'*P*H+R  */
    matmul("NN",n,n,nx,1.0,H,F1,1.0,Q1);
    if (!(info=matinv(Q1,n))) {
        matmul("NN",nx,n,n,1.0,F1,Q1,0.0,K1);   /* K=P*H*Q^-1  */
        matmul("NN",nx,1,n,1.0,K1,v,1.0,xp1);   /* xp=x+K*v  */
        matmul("NT",nx,nx,n,-1.0,K1,H,1.0,I1);  /* Pp=(I-K*H')*P  */
        matmul("NN",nx,nx,nx,1.0,I1,P_,0.0,P_);
    }

      for (i=0;i<nx;i++) x_[i]=xp1[i];

      /* Updated state estimates and state estimation erros covariance matrix*/
      for (i=0;i<nx;i++) xp[i]=x_[i];
      for (i=0;i<nx*nx;i++) Pp[i]=P_[i];

  /* Correct inertial navigation solution */

  /* Attitude - */
  /* Shin (2001) Cbn(+)=(I+E)*Cbn(-) */
  /* Groves(2013) Cbn(+)=(I3-[delta_att[3] ^])*Cbn(-)*/
  I[0]=I[4]=I[8]=1.0;
  I[1]=I[2]=I[3]=I[5]=I[6]=I[7]=0.0;
  vec2skew (x_, E);
  for (i=0;i<9;i++) Ea[9]=I[i]-E[i];
  matmul("NN",3,3,3,1.0,Ea,pvap->Cbn,0.0,Cbn_aux);
  for (i=0;i<9;i++) {
    pva_global.Cbn[i]=Cbn_aux[i];
  }

  Cbn2euler(imu,pvap->Cbn);
  for (i=0;i<3;i++) {pva_global.A[i]=pvap->A[i]=imu->aea[i];}

  /* Velocity - v(+)=v(-)-dv */
  for(i=3;i<6;i++) {pva_global.v[i-3]=pvap->v[i-3]=pvap->v[i-3]-x_[i];}

  /* Position - r(+)=r(-)-dr */
  for(i=6;i<9;i++) {pva_global.r[i-6]=pvap->r[i-6]=pvap->r[i-6]-x_[i];}

  /* Update Bias estimates Groves(2013) */

  /* Accelerometers bias - ba(+)=ba(-)+dba */
  for(i=9;i<12;i++) {imu_obs_global.ba[i-9]=imuobsp->ba[i-9]=imuobsp->ba[i-9]+x_[i];}

  /* Accelerometers bias - bg(+)=bg(-)+dbg */
  for(i=12;i<15;i++) {imu_obs_global.bg[i-12]=imuobsp->bg[i-12]=imuobsp->bg[i-12]+x_[i];}


  /* Next cycle */
  /* Update */
  imu_obs_global=*imuobsp;
  pva_global=*pvap;
  for (i=0;i<3;i++) pvagnss.r[i]=gnss_xyz_ini_pos[i];
  pvagnss.sec=ini_pos_time;


/*  printf("\nGNSS and INS are synchronized\n");
printf("\nFINAL IMU and PVA  AFTER *******************************{\n\n");
printf("acc.: %lf, %lf, %lf\n",imu_obs_global.fb[0],imu_obs_global.fb[1],imu_obs_global.fb[2] );
printf("Gyro.: %lf, %lf, %lf\n",imu_obs_global.wibb[0],imu_obs_global.wibb[1],imu_obs_global.wibb[2] );
printf("Ba.: %lf, %lf, %lf\n",imu_obs_global.ba[0],imu_obs_global.ba[1],imu_obs_global.ba[2] );
printf("Bgo.: %lf, %lf, %lf\n",imu_obs_global.bg[0],imu_obs_global.bg[1],imu_obs_global.bg[2] );
printf("P: %lf, %lf, %lf\n",pvap->r[0],pvap->r[1],pvap->r[2] );
printf("V: %lf, %lf, %lf\n",pvap->v[0],pvap->v[1],pvap->v[2] );
printf("A: %lf, %lf, %lf\n",pva_global.A[0],pva_global.A[1],pva_global.A[2] );
printf("*********************************************************}\n\n"); */

free(z);free(H);free(R);
free(xp);free(Pp); free(Q);
free(PHI);free(G);free(u);
free(x_);free(xk);free(P_);free(Pk);
free(v);free(Q_);free(v_);
free(P_aux);free(Q_aux);
free(F);
free(F1); free(Q1); free(K1); free(I1); free(xp1);

}

/* INS/GNSS Loosley coupled (LS) integration --------------------------------------------------------
* description: Compute PVA solution and ins biases by integrating INS and
* GNSS solutions
* args :
*       double*   gnss_xyz_ini_pos  I gnss position solution in {x,y,z}
*       double*   gnss_xyz_ini_cov  I gnss position variance/covariance (m^2)
* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or {c_ee,c_nn,c_uu,c_en,c_nu,c_ue}
*       double*   gnss_enu_vel      I gnss velocity in {e,n,up}
*       float     ini_pos_time      I position time in week and tow in gps time
* Ins/gnss model based on Grove(2013), from section 14.2.4 and Shin (2001)
*-----------------------------------------------------------------------------*/
extern void insgnssLC (double* gnss_xyz_ini_pos, double* gnss_xyz_ini_cov,\
   double* gnss_enu_vel, double ini_pos_time, um7pack_t *imu, pva_t *pvap, imuraw_t *imuobsp)
{
  /*imuraw_t imuobsp={0};
  um7pack_t imu={0};
  pva_t pvap={0};*/
  double *z,*H,*R,*xp,*Pp, *Q, *PHI,*G,*u, dt;
  double *xk, *x_, *Pk, *P_, *v, *v_, *Q_,*F;
  double *F1,*Q1,*K1,*I1,*xp1;
  double *P_aux,*Q_aux;
  double E[9],Ea[9],I[9], Cbn_aux[9],l[3]; /* Lever-arm from INS to GNSS antenna in b-frame {x(nav.direc),y(right),z(down)} in meters */
  int i,j,nx,n,nv, info,iter;
  char* imuf, datastring[150];
  char c=0;
  double conf_not_ini_time, pos[3];

//  printf("\n\n ************************** IT STARTS HERE *****************************************************{\n");

  /* GNSS antenna lever arm to the INS body-frame  */
  l[0]=0.009;l[1]=0.1653;l[2]=0.0;

  /* Update the latest solution determined in the KF
  imuobsp=imu_obs_global;
  pvap=pva_global;*/

  /* Read imu current measurements
  while(imu->status!=1){
    if(!fgets(datastring, 150, fimu)) return;
    parseimudata(datastring,&imu);
  }

  while(gnssimusync(ini_pos_time, &imu)!=1){ //Run ins-nav until catch up with GNSS time (if GNSS ahead)
    /* If GNSS behind IMU time
   if (imu->status==-1) break;
    /* Read imu current obs and compute an ins solution
    while( imu->status != 1 ){
    //imuf=fgets(datastring, 150, fimu);
    if( !fgets(datastring, 150, fimu)
     || ferror( fimu ) || feof( fimu ) ) break;
    parseimudata(datastring,&imu);
  }
}
*/
  /* Time span between measurements (s) ?? */
  /* Time span from previous measurement */
  conf_not_ini_time=imuobsp->sec;
  if (!imuobsp->sec) {
      imuobsp->sec= imu->sec;
      dt = imu->sec-imuobsp->sec; /* In fraction of seconds */
      pvagnss.sec= ini_pos_time;
    }else{
    dt= imu->sec-imuobsp->sec; /* In fraction of seconds */
    //printf("Time span: %lf \n",dt);
    dt=ini_pos_time-pvagnss.sec; /*Using Gnss time-span for the KF*/
  }

  /* Compute inertial navigation solution */
  /* Updating pvap and imuobsp with current ins solution and measurement  */
//  printf("\n *** INSNAV **************************************\n");
  if(insnav1(gnss_xyz_ini_pos, gnss_enu_vel, ini_pos_time, imu, pvap, imuobsp)!=1){
    printf("INS navigation computation error! \n");
  }

 /* Loop for refinement or feedback?, what condition ???*/
 // for (iter = 0; iter < 1; iter++) {

  /* Compute gnss-inertial integrated solution */
  if(gnssimusync(ini_pos_time, imu)==1){
  //  printf("\n\n **********GNSS/INS UPDATE ********************************\n\n" );
  pvap->sec=ini_pos_time;

    /* Number of states (nx) and measureprintf("\n\n **********GNSS/INS **************************************\n\n" );ments (n)  */
    nx=15; /* ins 15-states {de,dv,dr,dba,dbg} system model in local-nav-frame */
    n=  6; /* Gnss-ins differences in position and velocity */
    nv=n;

    xp=mat(nx,1); Pp=zeros(nx,nx); Q=zeros(nx,nx);
    PHI=mat(nx,nx); G=zeros(nx,nx-3);  u=mat(nx-3,1);
    z=mat(nv,1); H=zeros(nv,nx); R=zeros(nv,nv);
    x_=mat(nx,1); xk=mat(nx,1);
    P_=zeros(nx,nx); Pk=zeros(nx,nx);
    v=mat(nv,1); v_=mat(nv,1);
    Q_=mat(nx,nx); F=mat(nx,nx);
    P_aux=mat(nx,nx); Q_aux=mat(nx,nx);
    /* For the Filter*/
    F1=mat(nx,n);Q1=mat(n,nx);K1=mat(nx,n);I1=eye(nx);xp1=mat(nx,1);

    /* Initialize the ins error-state and covariances xp, Pp  (Groves 2013, page 594)*/
    if (iter==0) {
      for (i=0;i<nx;i++) xp[i]=0.0;
    }
  //  posvelP0(Pp,nx,rtk);
  //  attbiasP0(&imuobsp, pvap, Pp, nx);

    /* Generic intialization for Loosley coupled GNSS/INS  */
    initialize_P0(Pp,nx);

    //inssysmatrix(PHI, G, nx, &pvap, &imuobsp, dt);
    /* Transition matrix using  (14.73 - page 589) */
    inssysmatrix1(PHI, nx, pvap, imuobsp, dt);
    /* System noise */
    sysnoise1(Q, dt, nx);

    /* Kalman filter system propagation */

    /* Propagate state estimates using (3.14) noting that all states are zero
    due to closed-loop correction.*/
     for (i=0;i<nx;i++) x_[i]=0.0;
     //matmul("NN",nx,1,nx,1.0,PHI,xp,0.0,x_); /* x_= T.xp */

    /* Propagate state error covariance using (3.46, page 100) */

    for (i = 0; i < nx; i++) {
      for (j = 0; j < nx; j++) {
        P_aux[i*nx+j]= Pp[i*nx+j]+0.5*Q[i*nx+j]; /* ...*(Pp + 0.5*Q)*... */
      }
    }
    matmul("NT",nx,nx,nx,1.0,P_aux,PHI,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
    matmul("NN",nx,nx,nx,1.0,PHI,Q_aux,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

    for (i = 0; i < nx; i++) {
      for (j = 0; j < nv; j++) {
        P_[i*nx+j]= Q_[i*nx+j]+0.5*Q[i*nx+j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
      }
    }

    /* Decide whether to perform or not the measurement update  ? */


    /* Kalman filter or EKF measurement update */
    measvec(z, gnss_xyz_ini_pos, pvap, imuobsp, l);
    measmatrixH (H, nv, nx, pvap);
    measnoiseR (R, nv, gnss_xyz_ini_cov);
    //measnoiseR (R, nv, rtk);

    /* Innovation term */
    matmul("NN",nv,1,nx,1.0,H,x_,0.0,v_);  /* v=z-H*x_ */
    for (i=0;i<nv;i++) v[i]=z[i]-v_[i];

    for (i = 0; i < nv; i++) {
    //   printf("v[%d]:%lf\n",i,v[i]);
    }

    /* Update  */

    /*  how about some iterations??  */

//    printf("Before filter\n");
  //  if ((info=filter(x_,P_,HT,v,R,nx,nv))) { //status (0:ok,<0:error)

      matcpy(Q1,R,n,n);
      matcpy(xp1,x_,nx,1);
      matmul("NT",nx,n,nx,1.0,P_,H,0.0,F1);      /* Q=H'*P*H+R  */
      matmul("NN",n,n,nx,1.0,H,F1,1.0,Q1);
      if (!(info=matinv(Q1,n))) {
          matmul("NN",nx,n,n,1.0,F1,Q1,0.0,K1);   /* K=P*H*Q^-1  */
          matmul("NN",nx,1,n,1.0,K1,v,1.0,xp1);   /* xp=x+K*v  */
          matmul("NT",nx,nx,n,-1.0,K1,H,1.0,I1);  /* Pp=(I-K*H')*P  */
          matmul("NN",nx,nx,nx,1.0,I1,P_,0.0,P_);
      }

        for (i=0;i<nx;i++) x_[i]=xp1[i];

        /* Updated state estimates and state estimation erros covariance matrix*/
        for (i=0;i<nx;i++) xp[i]=x_[i];
        for (i=0;i<nx*nx;i++) Pp[i]=P_[i];

    /* Correct inertial navigation solution */

    /* Attitude - */
    /* Shin (2001) Cbn(+)=(I+E)*Cbn(-) */
    /* Groves(2013) Cbn(+)=(I3-[delta_att[3] ^])*Cbn(-)*/
    I[0]=I[4]=I[8]=1.0;
    I[1]=I[2]=I[3]=I[5]=I[6]=I[7]=0.0;
    vec2skew (x_, E);
    for (i=0;i<9;i++) Ea[9]=I[i]-E[i];
    matmul("NN",3,3,3,1.0,Ea,pvap->Cbn,0.0,Cbn_aux);
    for (i=0;i<9;i++) {
      pva_global.Cbn[i]=Cbn_aux[i];
    }

    for (i = 0; i < 9; i++) {
    //  printf("E[%d]: %lf\n",i, E[i]);
    }

    Cbn2euler(imu,pvap->Cbn);
    for (i=0;i<3;i++) {pvap->A[i]=imu->aea[i];}

    /* Velocity - v(+)=v(-)-dv */
    for(i=3;i<6;i++) {pvap->v[i-3]=pvap->v[i-3]-x_[i];}

    /* Position - r(+)=r(-)-dr */
    for(i=6;i<9;i++) {pvap->r[i-6]=pvap->r[i-6]-x_[i];}

    /* Update Bias estimates Groves(2013) */

    /* Accelerometers bias - ba(+)=ba(-)+dba */
    for(i=9;i<12;i++) {imuobsp->ba[i-9]=imuobsp->ba[i-9]+x_[i];}

    /* Accelerometers bias - bg(+)=bg(-)+dbg */
    for(i=12;i<15;i++) {imuobsp->bg[i-12]=imuobsp->bg[i-12]+x_[i];}


  /* GLOBAL FILES - IMU bias and LC Covariances output */
    fprintf(out_IMU_bias_file,"%f %lf %lf %lf %lf %lf %lf\n", ini_pos_time, \
      imuobsp->ba[0],imuobsp->ba[1],imuobsp->ba[2],\
      imuobsp->bg[0],imuobsp->bg[1], imuobsp->bg[2] );

      /* Generate KF uncertainty output record */
      fprintf(out_KF_SD_file,"%f\t", ini_pos_time);
      for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) {
          if (i==j) {
            fprintf(out_KF_SD_file,"%lf\t", sqrt(Pp[i*nx+j]) );
          }
        }
      }
      fprintf(out_KF_SD_file,"\n");


    /* Next cycle */
    /* Update
    imu_obs_global=imuobsp;
    pva_global=pvap;*/
    for (i=0;i<3;i++) pvagnss.r[i]=gnss_xyz_ini_pos[i];
    pvagnss.sec=ini_pos_time;


  /*  printf("\nGNSS and INS are synchronized\n");
  printf("\nFINAL IMU and PVA  AFTER *******************************{\n\n");
  printf("acc.: %lf, %lf, %lf\n",imu_obs_global.fb[0],imu_obs_global.fb[1],imu_obs_global.fb[2] );
  printf("Gyro.: %lf, %lf, %lf\n",imu_obs_global.wibb[0],imu_obs_global.wibb[1],imu_obs_global.wibb[2] );
  printf("Ba.: %lf, %lf, %lf\n",imu_obs_global.ba[0],imu_obs_global.ba[1],imu_obs_global.ba[2] );
  printf("Bgo.: %lf, %lf, %lf\n",imu_obs_global.bg[0],imu_obs_global.bg[1],imu_obs_global.bg[2] );
  printf("P: %lf, %lf, %lf\n",pva_global.r[0],pva_global.r[1],pva_global.r[2] );
  printf("V: %lf, %lf, %lf\n",pva_global.v[0],pva_global.v[1],pva_global.v[2] );
  printf("A: %lf, %lf, %lf\n",pva_global.A[0],pva_global.A[1],pva_global.A[2] );
  printf("*********************************************************}\n\n");// */

  free(z);free(H);free(R);
  free(xp);free(Pp); free(Q);
  free(PHI);free(G);free(u);
  free(x_);free(xk);free(P_);free(Pk);
  free(v);free(Q_);free(v_);
  free(P_aux);free(Q_aux);
  free(F);
  free(F1); free(Q1); free(K1); free(I1); free(xp1);


  }else{
    printf("\n GNSS and INS are not synchronized\n");
    /* Update
    imu_obs_global=imuobsp;
    pva_global=pvap;*/
    for (i=0;i<3;i++) pvagnss.r[i]=gnss_xyz_ini_pos[i];
    pvagnss.sec=ini_pos_time;
  }

 //} /* iteration end */

}

/* Determine heading from the map or Ground-truth trajectory -----------------*/
extern double headfrommap(double* closest_lane_pos){
  double next_lane_pos[3], origin_llh_pos[3], origin_xyz_pos[3];
  double dr[3], enu_closest[3], enu_next[3], az=0.0;
  int j;

  //Origin - BMO corner 1
  origin_llh_pos[0]= 45.94312323*D2R;
  origin_llh_pos[1]=-66.64106835*D2R;
  origin_llh_pos[2]= 52.7309;

  pos2ecef(origin_llh_pos,origin_xyz_pos);

  for (j=0;j<3;j++) next_lane_pos[j]=lane.buffer[(lane.currsearch+1)*3+j];

  /* Transform to local NED system  */
  for (j=0;j<3;j++) dr[j]=closest_lane_pos[j]-origin_xyz_pos[j];
  ecef2enu(origin_llh_pos, dr, enu_closest);

  for (j=0;j<3;j++) dr[j]=next_lane_pos[j]-origin_xyz_pos[j];
  ecef2enu(origin_llh_pos, dr, enu_next);

  return az=azt(enu_closest, enu_next);
}

/* Determine velocity using the previous and current closest map point position*/
extern void velfrommap(double *xyz_curr_pos, double *xyz_prev_pos, double *ned_vel, float dt){
  double vxyz[3], venu[3];
  int i;

  /* cartesian velocity */
  if (dt<=0.0) {
    for (i=0;i<3;i++) vxyz[i]=(xyz_curr_pos[i]-xyz_prev_pos[i])/(1.0);
  }else{
    for (i=0;i<3;i++) vxyz[i]=(xyz_curr_pos[i]-xyz_prev_pos[i])/(dt);
  }
  /* loca enu velocity */
    ecef2enu(xyz_prev_pos,vxyz,venu);
    ned_vel[2]=-venu[2]; /* Down */
    /* E,N,D velocity */
    ned_vel[0]=venu[1];
    ned_vel[1]=venu[0];
}

/* INS/MAP integration ---------------------------------------------------------
* description: Compute PVA solution and ins biases by integrating INS and
* the Ground-truth trajectory
* args :
*       double*   pva_t PVA_prev_sol  I system last PVA solution
*       float     prev_time      I position time in week and tow in gps time
*       imuraw_t* imu_curr_meas   I IMU measurement structure
*-----------------------------------------------------------------------------*/
extern void insmap (){
  int j;
  double xyz_ini_pos[3], xyz_ini_cov[3], ned_ini_vel[3], pos[3], head_angle=0.0;
  double xyz_prev_clst_pos[3];
  float ini_pos_time;
  char datastring[150];
  imuraw_t imu_obs_prev={0};
  um7pack_t imu_curr_meas={0};
  pva_t PVA_prev_sol={0};

  //FILE* facc, *fgyr;
  FILE* fp, *fv, *fatt;
  //facc = fopen("../out/insmap_test_acceleration.txt","a");
  //fgyr = fopen("../out/insmap_test_gyro.txt","a");
  fatt = fopen("../out/insmap_test_attitude.txt","a");
  fp = fopen("../out/insmap_test_pos.txt","a");
  fv = fopen("../out/insmap_test_vel.txt","a");


    /* Update the with previous solution determined in the navigation solution */
    imu_obs_prev=imu_obs_global;
    PVA_prev_sol=pva_global;
    for (j=0;j<3;j++) xyz_prev_clst_pos[j]=lane.buffer[j]; /* Previous closest point */

    /* Read imu current obs and compute an ins solution */
    while(imu_curr_meas.status!=1){
      if(!fgets(datastring, 150, fimu)) return;
      parseimudata(datastring,&imu_curr_meas);
      //printf("HERE, status: %d\n",imu_curr_meas.status);
    }
    ini_pos_time=imu_curr_meas.sec;
      /* Fetching data to the first buffer structure (lane.buffer)	*/
    //inibuff();

    /* First position solution     */
    if (norm(pva_global.r,3)==0.0) {
      /* Use first position of the buffer */
      for (j=0;j<3;j++) xyz_ini_pos[j]=lane.buffer[j];
      imu_obs_prev.sec=ini_pos_time;
    }else{
      /* use previous solution */
      pos2ecef(PVA_prev_sol.r, xyz_ini_pos);
    }

    /* Fixing time
    if (imu_curr_meas.sec-imu_obs_prev.sec==0.0) {
      ini_pos_time=imu_curr_meas.sec=imu_curr_meas.sec+0.01;
    }  */

    /* Lane covariance information */
    for (j=0;j<3;j++) xyz_ini_cov[j]=0.01; //Use some value

    /* Closest point computation  			*/
    //clstpt(&xyz_ini_pos);

    /* Closest point position in the lane is the initial INS position    */
    for (j=0;j<3;j++) xyz_ini_pos[j]=lane.buffer[lane.currsearch*3+j];

    /* Constrain velocity due to GNSS outage  */

     /* If previous solution was a GNSS/INS+MAP, and use previous until after
     a few seconds(how many seconds?)
     for (j=0;j<3;j++) ned_ini_vel[j]=pvagnss.v[j]=PVA_prev_sol.v[j]; */

     /* Otherwise, previous velocity will not be reliable, deduce it from map */
     if (norm(pva_global.r,3)==0.0){// || imu_curr_meas.sec<405026.3) { /* Static constrain exp1:405026.3 exp2:220149.0 */
       for (j=0;j<3;j++) ned_ini_vel[j]=pvagnss.v[j]=0.0;
     }else{
       /* Deduce from map */
       velfrommap(xyz_ini_pos,xyz_prev_clst_pos,&ned_ini_vel, ini_pos_time-imu_obs_prev.sec);
       for (j=0;j<3;j++) pvagnss.v[j]=ned_ini_vel[j];
     }

     //printf("Times: curr %f, prev: %f \n",ini_pos_time,imu_obs_prev.sec );

     //printf("VELFROM MAP: %lf, %lf, %lf\n%lf, %lf, %lf\nTime differences:%f",ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2],
      //  xyz_prev_clst_pos[0],xyz_prev_clst_pos[1],xyz_prev_clst_pos[2],ini_pos_time-imu_obs_prev.sec);

    /* Heading (body-frame) or Azimuth (local-nav-frame)? constrain using lanes  */
    head_angle=headfrommap(xyz_ini_pos);

   //printf("Azimuth: %lf\n",head_angle*R2D);
   PVA_prev_sol.A[2]=head_angle;

    /* Compute inertial navigation solution */
    /* Updating PVA_prev_sol and imu_obs_prev with current IMU measurement */
    if(insnav1(xyz_ini_pos, ned_ini_vel, ini_pos_time, &imu_curr_meas, &PVA_prev_sol, &imu_obs_prev)!=1){
          printf("INS navigation computation error! \n");
    }else{

    //  printf("IMU measurement time: %f, status: %d, count: %d\n",imu_curr_meas.sec,imu_curr_meas.status, imu_curr_meas.count);
      /*  printf("acc.: %lf, %lf, %lf\n",imu_curr_meas.a[0],imu_curr_meas.a[1],imu_curr_meas.a[2] );
        printf("Gyro.: %lf, %lf, %lf\n",imu_curr_meas.g[0],imu_curr_meas.g[1],imu_curr_meas.g[2] );
        printf("Mag.: %lf, %lf, %lf\n",imu_curr_meas.m[0],imu_curr_meas.m[1],imu_curr_meas.m[2] );   */
      //  printf("P: %lf, %lf, %lf\n",imu_curr_meas.r[0]*R2D,imu_curr_meas.r[1]*R2D,imu_curr_meas.r[2] );
      //  printf("V: %lf, %lf, %lf\n",imu_curr_meas.v[0],imu_curr_meas.v[1],imu_curr_meas.v[2] );
      //  printf("A: %lf, %lf, %lf\n",imu_curr_meas.aea[0],imu_curr_meas.aea[1],imu_curr_meas.aea[2] );

      /* Call Loosley-coupled algorithm */
      //printf("\n *********** INS/MAP UPDATE *******************************\n" );
    //  ins_LC (xyz_ini_pos, xyz_ini_cov, ned_ini_vel, ini_pos_time, &imu_curr_meas, &PVA_prev_sol, &imu_obs_prev);


      /* Update the with latest solution determined in the navigation solution */
      imu_obs_global=imu_obs_prev;
      pva_global=PVA_prev_sol;
    /*  Printing imu measurements */
    //  fclose(fvp);
    //fprintf(facc,"%f %lf %lf %lf\n",imu_curr_meas.sec, imu_curr_meas.a[0], imu_curr_meas.a[1], imu_curr_meas.a[2]); // accelerations
    //fprintf(fgyr,"%f %lf %lf %lf\n",imu_obs_global.sec, imu_curr_meas.g[0]*R2D, imu_curr_meas.g[1]*R2D, imu_curr_meas.g[2]*R2D); // rate-rotations
    fprintf(fatt,"%f %lf %lf %lf\n",imu_obs_global.sec, pva_global.A[0]*R2D, pva_global.A[1]*R2D, pva_global.A[2]*R2D); // Euler-rotations
    //fprintf(facc,"%f %f %f %f \n",imu.sec, imu.m[0], imu.m[1], imu.m[2]); // magnetometers

    /* Printing velocity and/or positions into a file*/
    ecef2pos(pvagnss.r, pos);
    fprintf(fp,"%f %lf %lf %lf\n",imu_obs_global.sec, pva_global.r[0]*R2D, pva_global.r[1]*R2D, pva_global.r[2]); // in llh
    fprintf(fv,"%f %lf %lf %lf\n",imu_obs_global.sec, pva_global.v[0], pva_global.v[1], pva_global.v[2]); // in NED
  }

  /* Update the with latest solution determined in the navigation solution */
  imu_obs_global=imu_obs_prev;
  pva_global=PVA_prev_sol;

  //fclose(facc); fclose(fgyr);
  fclose(fp); fclose(fv); fclose(fatt);
}

/* From Shin 2001, page 42   */
//void IMU_meas_interpolation(double t_imu_prev, double t_imu_curr, double t_gps, \
  um7pack_t *imu_curr, imuraw_t *imu_prev){
void IMU_meas_interpolation(double t_imu_prev, double t_imu_curr, double t_gps, \
  ins_states_t *imu_curr, imuraw_t *imu_prev){
    int i;
    double imu_t_gps_a[3], imu_t_gps_g[3];

    for (i=0;i<3;i++) imu_t_gps_a[i]=imu_prev->fb0[i]+\
    ((imu_curr->data.fb0[i]-imu_prev->fb0[i])/(t_imu_curr-t_imu_prev))*(t_gps-t_imu_prev);

    for (i=0;i<3;i++) imu_t_gps_g[i]=imu_prev->wibb0[i]+\
    ((imu_curr->data.wibb0[i]-imu_prev->wibb0[i])/(t_imu_curr-t_imu_prev))*(t_gps-t_imu_prev);

    /* Updating imu_curr */
    for (i=0;i<3;i++) imu_curr->data.fb0[i]= imu_t_gps_a[i];
    for (i=0;i<3;i++) imu_curr->data.wibb0[i]= imu_t_gps_g[i];
    imu_curr->time=imu_curr->data.sec=t_gps;
  }

/* From plotdata.cpp in ../app/rtkplot */
// update Multipath ------------------------------------------------------------
void UpdateMp(const obs_t *Obs, const nav_t *Nav)
{
    obsd_t *data;
    double lam1,lam2,I,C,*Mp, B;
    int i,j,k,f1,f2,sat,sys,per,per_=-1,n;

   /*
    for (i=0;i<NFREQ+NEXOBS;i++) {
        delete [] Mp[i]; Mp[i]=NULL;
    }*/

    if (n<=0) return;


    /*for (i=0;i<NFREQ+NEXOBS;i++) {
        Mp[i]=new double[Obs.n];
    }*/

    Mp=mat(n,NFREQ+NEXOBS);

    //ReadWaitStart();
    //ShowLegend(NULL);

    for (i=0;i<n;i++) {
        data=Obs->data+i;
        sys=satsys(data->sat,NULL);

        for (j=0;j<NFREQ+NEXOBS;j++) {
            Mp[i*(NFREQ+NEXOBS)+j]=0.0;

            code2obs(data->code[j],&f1);

            if (sys==SYS_CMP) {
                if      (f1==5) f1=2; /* B2 */
                else if (f1==4) f1=3; /* B3 */
            }
            if      (sys==SYS_GAL) f2=f1==1?3:1; /* E1/E5a */
            else if (sys==SYS_SBS) f2=f1==1?3:1; /* L1/L5 */
            else if (sys==SYS_CMP) f2=f1==1?2:1; /* B1/B2 */
            else                   f2=f1==1?2:1; /* L1/L2 */

            lam1=satwavelen(data->sat,f1-1,&Nav);
            lam2=satwavelen(data->sat,f2-1,&Nav);
            if (lam1==0.0||lam2==0.0) continue;

            if (data->P[j]!=0.0&&data->L[j]!=0.0&&data->L[f2-1]) {
                C=SQR(lam1)/(SQR(lam1)-SQR(lam2));
                I=lam1*data->L[j]-lam2*data->L[f2-1];
                Mp[i*(NFREQ+NEXOBS)+j]=data->P[j]-lam1*data->L[j]+2.0*C*I;
            }
        }
    }
    for (sat=1;sat<=MAXSAT;sat++) for (i=0;i<NFREQ+NEXOBS;i++) {
        sys=satsys(sat,NULL);

        for (j=k=n=0,B=0.0;j<Obs->n;j++) {
            if (Obs->data[j].sat!=sat) continue;

            code2obs(Obs->data[j].code[i],&f1);

            if (sys==SYS_CMP) {
                if      (f1==5) f1=2; /* B2 */
                else if (f1==4) f1=3; /* B3 */
            }
            if      (sys==SYS_GAL) f2=f1==1?3:1;
            else if (sys==SYS_CMP) f2=f1==1?2:1;
            else                   f2=f1==1?2:1;

            if ((Obs->data[j].LLI[i]&1)||(Obs->data[j].LLI[f2-1]&1)||
                fabs(Mp[i*(NFREQ+NEXOBS)+j]-B)>0.1) {  //THRES_SLIP

                for (;k<j;k++) if (Obs->data[k].sat==sat) Mp[i*(NFREQ+NEXOBS)+k]-=B;
                B=Mp[i*(NFREQ+NEXOBS)+j]; n=1; k=j;
            }
            else {
                if (n==0) k=j;
                B+=(Mp[i*(NFREQ+NEXOBS)+j]-B)/++n;
            }
        }
        if (n>0) {
            for (;k<j;k++) if (Obs->data[k].sat==sat) Mp[i*(NFREQ+NEXOBS)+k]-=B;
        }
        per=sat*100/MAXSAT;
        if (per!=per_) {
            //ShowMsg(s.sprintf("updating multipath... (%d%%)",(per_=per)));
            //Application->ProcessMessages();
        }
    }
    //ReadWaitEnd();
}
/* transpose matrix-----------------------------------------------------------
 * args   : double  *A      I   matrix
 *          int      n      I   rows of transpose matrix
 *          int      m      I   cols of transpose matrix
 *          double  *At     O   transpose of matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void matt(const double *A,int n,int m,double *At)
{
    int i,j; for (i=0;i<n;i++) for (j=0;j<m;j++) At[i+j*n]=A[j+i*m];
}

/* Convert row order matrix into column order ----------------------------
 * args    :  double *A   I  input matrix in row order to convert
 *            int m,n     I   rows and cols of input matrix
 *            double *B   O  output matrix in row order
 * return  : none
 * ----------------------------------------------------------------------*/
extern void row_to_column_order(const double *A,int m,int n, double *B)
{
    int i,j;
    for (i=0;i<m;i++) {
        for (j=0;j<n;j++) B[j*m+i]=A[i*n+j];
    }
}

/* asign a block matrix to another matrix by giving index -----------------
 * args    :  double *A   IO  input matrix for asign in elements
 *            int m,n     I   rows and cols of input matrix
 *            double *B   I   asign matrix
 *            int p,q     I   rows and cols of asign matrix
 *            int isr,isc I   start row and col to asign matrix
 * return  : none
 * ----------------------------------------------------------------------*/
extern void asi_blk_mat(double *A,int m,int n,const double *B,int p ,int q,
                        int isr,int isc)
{
    int i,j; for (i=isr;i<isr+p;i++) {
        for (j=isc;j<isc+q;j++) A[i+j*m]=B[(i-isr)+(j-isc)*p];
    }
}
/* correction direction cosine matrix by attitude errors---------------------
 * args   :  double *dx  I   attitude errors
 *           double *C   IO  direction cosine matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void corratt(const double *dx,double *C)
{
    int i;
    double T[9],*I=eye(3);

    skewsym3(dx,T);
    for (i=0;i<9;i++) I[i]-=T[i];

    matcpy(T,C,3,3);
    matmul("NN",3,3,3,1.0,I,T,0.0,C);
    free(I);
}

/* correction imu accl. and gyro. measurements-----------------------------
 * args     :    I  double *accl  imu accl. measurements
 *               I  double *gyro  imu gyro. measurements
 *               I  double *Ma    non-orthogonal between sensor axes and
 *                                body frame for accelerometers
 *               I  double *Mg    non-orthogonal between sensor axes and
 *                                body frame for gyroscopes
 *               I  double *ba    accelemeter-bias
 *               I  double *bg    gyro-bias
 *               I  double *Gg    g-dependent bias for a gyro triad
 *               O  double *cor_accl
 *                         *cor_gyro  corrected imu accl. and gyro. measurements
 * return : none
 * -----------------------------------------------------------------------*/
extern void ins_errmodel2(const double *accl,const double *gyro,const double *Ma,
                          const double *Mg,const double *ba,const double *bg,
                          const double *Gg,double *cor_accl,double *cor_gyro)
{
    int i,j;
    double Mai[9],Mgi[9],*I=eye(3),T[9]={0},Gf[3]={0};

    for (i=0;i<3;i++) for (j=0;j<3;j++) Mai[i+j*3]=I[i+j*3]+Ma[i+j*3];
    for (i=0;i<3;i++) for (j=0;j<3;j++) Mgi[i+j*3]=I[i+j*3]+Mg[i+j*3];

    if (!matinv(Mai,3)&&!matinv(Mgi,3)) {
        matmul("NN",3,1,3,1.0,Mai,accl,0.0,T);
        matmul("NN",3,1,3,1.0,Mgi,gyro,0.0,T+3);
    }
    if (cor_accl) {
        for (i=0;i<3;i++) cor_accl[i]=T[i]-ba[i];
        matmul("NN",3,1,3,1.0,Gg,accl,0.0,Gf);
    }
    if (cor_gyro) {
        for (i=0;i<3;i++) cor_gyro[i]=T[3+i]-bg[i]-Gf[i];
    }
    free(I);
}

/* get the acceleration of body in ecef-frame by input acceleration measurements
 * args     :  double *fib    I  input acceleration measurements
 *             double *Cbe    I  body-frame to ecef-frame convert matrix
 *             double *re     I  position of imu-body in ecef-frame
 *             double *ve     I  velocity of imu-body in ecef-frame
 *             double *ae     O  acceleration of imu-body in ecef-frame
 * return :none
 * --------------------------------------------------------------------------*/
extern void getaccl(const double *fib,const double *Cbe,const double *re,
                    const double *ve,double *ae)
{
    int i; double fe[3],ge[3],cori[3];

    matmul3v("N",Cbe,fib,fe);
    Gravity_ECEF(re,ge);
    matmul3v("N",Omge,ve,cori);

    for (i=0;i<3;i++) {
        ae[i]=fe[i]+ge[i]-2.0*cori[i];
    }
}

/* close loop for non-holonomic constraint-----------------------------------*/
extern void clp(ins_states_t *ins,const insgnss_opt_t *opt,const double *x)
{
    int i;
    double *I=eye(3),fibc[3],omgbc[3],ang[3];

    /* close-loop attitude correction */
    corratt(x,ins->Cbe);

    /* close-loop velocity and position correction */
    ins->ve[0]-=x[xiV()+0];
    ins->ve[1]-=x[xiV()+1];
    ins->ve[2]-=x[xiV()+2]; 

    ins->re[0]-=x[xiP()+0];
    ins->re[1]-=x[xiP()+1];
    ins->re[2]-=x[xiP()+2]; 

    /* close-loop accl and gyro bias */
    ins->data.ba[0]+=x[xiBa()+0];
    ins->data.ba[1]+=x[xiBa()+1];
    ins->data.ba[2]+=x[xiBa()+2];

    ins->data.bg[0]+=x[xiBg()+0];
    ins->data.bg[1]+=x[xiBg()+1];
    ins->data.bg[2]+=x[xiBg()+2];

    /* close-loop residual scale factors of gyro and accl */
    //for (i=isg;i<isg+nsg;i++) ins->Mg[i-isg+(i-isg)*3]+=x[i];

    //for (i=isa;i<isa+nsa;i++) ins->Ma[i-isa+(i-isa)*3]+=x[i];

    /* correction imu accl and gyro measurements data 
    ins_errmodel2(ins->data.fb0,ins->data.wibb0,
                  ins->data.Ma,ins->data.Mg,
                  ins->data.ba,ins->data.bg,ins->data.Gg,
                  fibc,omgbc);
    */
    /* correction imu-body accelerometer */
    getaccl(fibc,ins->Cbe,ins->re,ins->ve,ins->data.fbe);

    free(I);
}

/* close loop for non-holonomic constraint-----------------------------------*/
void clp1(pva_t *PVA_sol, const double *x)
{
    int i;
    double Skew_x_est_new[9]={0.0}, I[9]={1,0,0,0,1,0,0,0,1};
    double I_less_Skew[9], est_C_b_e_new[9], est_IMU_bias_new[6];
    double est_v_eb_e_new[3], est_r_eb_e_new[3];
    double r[3];

  //  pos2ecef(PVA_sol->r,r);
    for (i=0;i<3;i++) r[i]= PVA_sol->re[i]; /* xyz_ecef */

    /* Correct attitude, velocity, and position using (14.7-9) */

    Skew_symmetric(x, Skew_x_est_new);
    for (i=0;i<9;i++) I_less_Skew[i]= I[i]-Skew_x_est_new[i];
    matmul_row("NN",3,3,3,1.0,I_less_Skew,PVA_sol->Cbe,0.0,est_C_b_e_new);


    for (i=3;i<6;i++) est_v_eb_e_new[i-3] = PVA_sol->ve[i-3] - x[i];
    for (i=6;i<9;i++) est_r_eb_e_new[i-6] = r[i-6] - x[i];

    /* Update IMU bias and GNSS receiver clock estimates */
    for (i=0;i<6;i++) est_IMU_bias_new[i] = PVA_sol->out_IMU_bias_est[i] + x[i+9];

    //PVA_sol->clock_offset_drift[0] = x[15];
    //PVA_sol->clock_offset_drift[1] = x[16];

    //ecef2pos(r,PVA_sol->r);
    for (i=0;i<3;i++) PVA_sol->ve[i]=est_v_eb_e_new[i];
    for (i=0;i<9;i++) PVA_sol->Cbe[i]=est_C_b_e_new[i];
    for (i=0;i<6;i++) PVA_sol->out_IMU_bias_est[i]=est_IMU_bias_new[i];

    /* Update the Navigation Solution */

    double est_L_b=0.0, est_lambda_b=0.0, est_h_b=0.0, est_v_eb_n[3]={0.0};
    double est_C_b_n[9]={0.0},est_C_b_n_T[9]={0.0};

    /* Convert navigation solution to NED  */
    ECEF_to_NED(est_r_eb_e_new,est_v_eb_e_new,est_C_b_e_new,\
      &est_L_b,&est_lambda_b,&est_h_b,\
      est_v_eb_n,est_C_b_n);

    /* Transposing Cbn */
    row_to_column_order(est_C_b_n,3,3,est_C_b_n_T);

    PVA_sol->r[0] = est_L_b;
    PVA_sol->r[1] = est_lambda_b;
    PVA_sol->r[2] = est_h_b;
    for (i=0;i<3;i++) PVA_sol->v[i] = est_v_eb_n[i];
    for (i=0;i<9;i++) PVA_sol->Cbn[i]=est_C_b_n[i];
    CTM_to_Euler(PVA_sol->A, est_C_b_n_T);

}

/* zero velocity update for ins navigation -----------------------------------
 * args    :  insstate_t *ins  IO  ins state
 *            insopt_t *opt    I   ins options
 *            imud_t *imu      I   imu measurement data
 *            int flag         I   static flag (1: static, 0: motion)
 * return  : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int zvu1(pva_t *PVA_sol, const um7pack_t *imu, int nx)
{
    int info=0,i, j;
    static int nz=0;
    double *x,*H,*R,*v, *P, I[9]={-1,0,0,0,-1,0,0,0,-1};

    trace(3,"zvu:");
    printf("zvu\n");

    //flag&=nz++>MINZC?nz=0,true:false;

    //if (!flag) return info;

    x=zeros(1,nx); H=zeros(3,nx);
    R=zeros(3,3); v=zeros(3,1); P=zeros(nx,nx);

    for (i = 0; i < nx; i++) x[i]=0.000000000000001;

    /* sensitive matrix */
    asi_blk_mat(H,3,nx,I,3,3,0,3);

    printf("H_matrix in zvu\n");
     for (i = 0; i < 3; i++) {
       for (j = 0; j < nx; j++) {
         printf("%lf ",H[j*3+i]);
       }
       printf("\n");
     }

    /* variance matrix */
    R[0]=R[4]=R[8]=sqrt(0.05);

    v[0]=PVA_sol->ve[0];
    v[1]=PVA_sol->ve[1];
    v[2]=PVA_sol->ve[2]; /* residual vector */

    for (i = 0; i < nx; i++) {
      for (j = 0; j < nx; j++) {
        P[j*nx+i] = PVA_sol->P[i*nx+j]; 
      }
    }  

    printf("P_matrix in zvu\n");
     for (i = 0; i < nx; i++) {
       for (j = 0; j < nx; j++) {
         (i==j?printf("%lf ",PVA_sol->P[i*nx+j]):0.0);
       }
       printf("\n"); 
     }


    if (norm(PVA_sol->v,3)<0.1&&norm(imu->g,3)<(10.0*D2R)) {

        /* ekf filter */
        info=filter(x,P,H,v,R,nx,3);

        printf("\nParameter vector zvu after\n");
         for (i = 0; i < 17; i++) {
           printf("%.10lf\n", x[i]);
         }

        /* solution fail */
        if (info) {
            trace(2,"zero velocity update filter error\n");
            info=0;
        }
        else {
            /* solution ok */
            //ins->stat=INSS_ZVU;
            info=1;
            clp1(PVA_sol,x);
            /* Full weight matrix  */
            for (i=0; i<nx; i++) {
              for (j=0;j<nx;j++) {
                PVA_sol->P[i*nx+j]=P[j*nx+i];
              }
            }
            trace(3,"zero velocity update ok\n");
        }
    }
    free(x); free(H);
    free(R); free(v); free(P);
    return info;
}

/* convert dcm to roll-pitch-yaw---------------------------------------------
 * args  :  double *C   I  direction cosine matrix
 *          double *rpy O  roll,pitch,yaw {rad}
 * note: C=Rz*Ry*Rx
 * return: none
 * --------------------------------------------------------------------------*/
extern void c2rpy(const double *C,double *rpy)
{
    rpy[0]=atan2(-C[5],C[8]);
    rpy[1]=asin ( C[2]);
    rpy[2]=atan2(-C[1],C[0]);
}

/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jacobian_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
}


/* measurement sensitive-matrix for non-holonomic----------------------------*/
static int bldnhc1(const um7pack_t *imu, const double *Cbe,
                  const double *ve,int nx,double *v,double *H,double *R)
{
  int i, j, nv,IA,IV;
  double C[9],CT[9], T[9], TT[9], vb[3],r[2],S[9];
  double T1[9], Ceb[9], I_x_Cbe[2*3]={0.0}, I[6]={0,-1,0,0,0,-1};
  double I_x_Cbe_T[2*3]={0.0};

  trace(3,"bldnhc:\n");

  /* Row to column order Or matrix transpose */
  row_to_column_order(Cbe,3,3, Ceb);

  /* velocity in body-frame */
  matmul_row("NN",3,1,3,1.0, Ceb, ve, 0.0, vb);

  matmul_row("NN",2,3,3,1.0,I,Ceb, 0.0, I_x_Cbe);

  /* Row to column order Or matrix transpose */
  row_to_column_order(I_x_Cbe,2,3,I_x_Cbe_T);

  Skew_symmetric(ve, CT);
  matmul_row("NN",3,3,3,1.0,Ceb,CT,0.0,TT);

  /* Passing T to column order */
  row_to_column_order(TT,3,3, T);

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      C[i*3+j]=Cbe[i*3+j]; /* Cbe(row) == Ceb(column) */
    }
  }

  //#define UPD_IN_EULER   0                /* update attitude in euler angles space when coupled measurement */

  /*
  #if UPD_IN_EULER */
  if (0) {
  jacobian_prot_pang(Cbe,S);
  matcpy(T1,T,3,3);
  matmul("NN",3,3,3,1.0,T1,S,0.0,T);
  }
  /* #endif
*/
printf("C:\n");
for (i=0;i<3;i++) {
  for (j=0;j<3;j++) {
    printf("%lf ", C[i*3+j]);
  }
  printf("\n");
}
printf("T:\n");
for (i=0;i<3;i++) {
  for (j=0;j<3;j++) {
    printf("%lf ", T[i*3+j]);
  }
  printf("\n");
}

/* build residual vector */
    for (nv=0,i=1;i<3;i++) {

        /* check velocity measurement */
        if (fabs(vb[i])>0.5) { //MAXVEL=0.5
            trace(2,"too large velocity measurement\n");
            continue;
        }
        /* check gyro measurement */
        if (fabs(norm(imu->g,3))>30.0*D2R) {
            trace(2,"too large vehicle turn\n");
            continue;
        }
        /* Include the first block when lever arm is known */
        //H[0+nv*nx]=T[i]; H[0+1+nv*nx]=T[i+3]; H[0+2+nv*nx]=T[i+6];
        H[3+nv*nx]=C[i]; H[3+1+nv*nx]=C[i+3]; H[3+2+nv*nx]=C[i+6];
        //H[3+nv*nx]=C[i]; H[3+3+nv*nx]=C[i+3]; H[3+6+nv*nx]=C[i+6];

        v[nv  ]=vb[i];
        r[nv++]=SQR(0.05);//VARVEL;
    }
    for (i=0;i<nv;i++) R[i+i*nv]=r[i];
    return nv;
}

/* using non-holonomic constraint for ins navigation---------------------------
 * args    :  insstate_t* ins  IO  ins state
 *            insopt_t* opt    I   ins options
 *            imud_t* imu      I   imu measurement data
 * return  : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int nhc1(pva_t *PVA_sol, const um7pack_t *imu, int nx)
{
    int info=0,nv,i,j;
    double *H,*HT,*v,*R,*x, *P;

    trace(3,"nhc:\n");
    printf("nhc:\n");

    H=zeros(2,nx); HT=zeros(2,nx); R=zeros(2,2);
    v=zeros(2,1); x=zeros(1,nx); P=zeros(nx,nx);

    for (i = 0; i < 17; i++) x[i]=0.000000000000001;

    printf("P_matrix\n");
     for (i = 0; i < nx; i++) {
       for (j = 0; j < nx; j++) {
         P[j*nx+i]=PVA_sol->P[i*nx+j];
         (i==j?printf("%lf ",P[j*nx+i]):0.0);
       }
       printf("\n");
     }

    nv=bldnhc1(imu,PVA_sol->Cbe,PVA_sol->ve,nx,v,H,R);

    printf("H_matrix row\n");
     for (i = 0; i < 2; i++) {
       for (j = 0; j < nx; j++) {
         printf("%lf ",H[i*nx+j]);
       }
       printf("\n");
     }

     /* Passing T to column order */
     row_to_column_order(H,2,nx, HT);

    printf("H_matrix column\n");
     for (i = 0; i < 2; i++) {
       for (j = 0; j < nx; j++) {
         printf("%lf ",HT[i*2+j]);
       }
       printf("\n");
     }
printf("\n");


    if (nv>0) {

        /* kalman filter */
        info=filter(x,P,H,v,R,nx,nv);

        printf("\nParameter vector nhc after\n");
         for (i = 0; i < 17; i++) {
           printf("%.10lf\n", x[i]);
         }

        /*  check ok? */
        if (info) {
            trace(2,"non-holonomic constraint filter fail\n");
            info=0;
        }
        else {
            /* solution ok */
            //ins->stat=INSS_NHC;
            info=1;
            clp1(PVA_sol,x);
            /* Full weight matrix  */
            for (i=0; i<nx; i++) {
              for (j=0;j<nx;j++) {
                PVA_sol->P[i*nx+j]=P[j*nx+i];
              }
            }
            trace(3,"use non-holonomic constraint ok\n");
        }
    }
    free(H); free(HT); free(v);
    free(R); free(x);
    return info;
}


/* Core function -------------------------------------------------------------
* Description: Receive raw GNSS and INS data and determine a PVA solution
* args:
* rtk_t *rtk  IO
* const obsd_t *obs I
* int n I
* const nav_t *nav  I
* return:
* obs.: this function is ran inside rtkpos function of RTKlib
------------------------------------------------------------------------------*/
extern void core(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav){
  double xyz_ini_pos[3], xyz_ini_cov[3], gnss_vel_cov[6];
  double pos[3], head_angle=0.0, vxyz[3], venu_fromgnss[3];
  double xyz_prev_clst_pos[3], llh_pos[3], aux, gan[3], wiee[3];
  double gnss_xyz_ini_cov[6], gnss_ned_cov[3], gnss_vel_ned_cov[6];
  double Qvenu[9]={0.0}, Qposenu[9]={0.0}, P[9]={0.0}, P_vel[9];
  double Qenu[9]={0};
  double ned_ini_vel[3]={0};
  double ini_pos_time;
  double enu_prev[3]={0.0}, enu_curr[3]={0.0}, azmt=0.0;
  char datastring[150];
  imuraw_t imu_obs_prev={0};
  um7pack_t imu_curr_meas={{0}};
  pva_t PVA_prev_sol={{0}};
  int i, j, gnss_vel=0;
  int TC_or_LC=1; /* TC:1 and LC=0 */
  int Tact_or_Low_IMU=1; /* Tactical=1 and Low grade=0*/
  //FILE *ppp_llh;
  char str[100], check=1;
  double G = 9.80665;
  /* ZUPT detection, based on   GREJNER-BRZEZINSKA et al. (2002) */
  // horizontal velocity and gyro components tolerance leves based on static INS data:
  double vn0, ve0, vn0std, ve0std, gyrx0, gyry0, gyrx0std, gyry0std;
   vn0 = -0.002743902895411; ve0 = -0.002219817510341;
   vn0std = 0.001958782893669; ve0std = 0.001549122618107;
   gyrx0 = -0.000152716; gyry0 = 0.000386451;
   gyrx0std = 0.000483409; gyry0std = 0.001271191;

  printf("\n *****************  CORE BEGINS ***********************\n");

  /* Update with previous solution */
  imu_obs_prev=imu_obs_global;
  PVA_prev_sol=pva_global;
  PVA_prev_sol.t_s=pvagnss.sec; /* last state/GNSS estimation time */

  /* Initialize time from GNSS */
  ini_pos_time=time2gpst(rtk->sol.time,NULL);

  /* Initialize position from GNSS SPP */
  for (j=0;j<3;j++) xyz_ini_pos[j]=rtk->sol.rr[j];

  /* Velocity from GNSS */
  for (j=0;j<3;j++) vxyz[j]=(xyz_ini_pos[j]-pvagnss.r[j])/\
  (ini_pos_time-pvagnss.sec);
  ecef2enu(pvagnss.r,vxyz,venu_fromgnss); /* Here is ENU! */
  /* To NED */
  pvagnss.v[0]=venu_fromgnss[1];
  pvagnss.v[1]=venu_fromgnss[0];
  pvagnss.v[2]=-venu_fromgnss[2];

  printf("NED_ini_vel.: t:%lf - %lf %lf %lf\n", ini_pos_time, rtk->sol.rr[3],rtk->sol.rr[4],rtk->sol.rr[5]);


  /* Velocity from GNSS in case there is no doppler observable */
  if (fabs(rtk->sol.rr[3]) <= 0.000001 && fabs(rtk->sol.rr[4]) <= 0.000000001 && \
  fabs(rtk->sol.rr[5]) <= 0.000000001) {
    /* NO VELOCITY */
    for (j=0;j<3;j++) ned_ini_vel[j]=venu_fromgnss[j];
    printf("NO_velocity: t:%lf - %lf %lf %lf\n", ini_pos_time, ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2]);
    gnss_vel=0;
  }else{
    ecef2pos(xyz_ini_pos,llh_pos);
    ecef2enu(llh_pos,rtk->sol.rr+3,ned_ini_vel); /* Here is ENU! */
    printf("Yes_velocity: t:%lf - %lf %lf %lf\n", ini_pos_time, ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2]);
    gnss_vel=1;
  }

//  ecef2pos(xyz_ini_pos,llh_pos);
//  ecef2enu(llh_pos,rtk->sol.rr+3,ned_ini_vel); /* Here is ENU! */

  /* ENU to NED velocity from gnss solution */
     aux=ned_ini_vel[1];
     ned_ini_vel[1]=ned_ini_vel[0];
     ned_ini_vel[0]=aux;
     ned_ini_vel[2]=-ned_ini_vel[2];

  /* Position covariance and velocity */
  for (j = 0; j < 6; j++) gnss_xyz_ini_cov[j]=rtk->sol.qr[j];
  P[0]     =rtk->sol.qr[0]; /* xx or ee */
  P[4]     =rtk->sol.qr[1]; /* yy or nn */
  P[8]     =rtk->sol.qr[2]; /* zz or uu */
  P[1]=P[3]=rtk->sol.qr[3]; /* xy or en */
  P[5]=P[7]=rtk->sol.qr[4]; /* yz or nu */
  P[2]=P[6]=rtk->sol.qr[5]; /* zx or ue */
  covenu(llh_pos,P,Qposenu);
  /* Using it where? */
  gnss_ned_cov[0]=Qposenu[4];
  gnss_ned_cov[1]=Qposenu[0];
  gnss_ned_cov[2]=Qposenu[8];

  /* Velocity NED covariance from GNSS */
  P_vel[0]     =rtk->sol.qrv[0]; /* xx or ee */
  P_vel[4]     =rtk->sol.qrv[1]; /* yy or nn */
  P_vel[8]     =rtk->sol.qrv[2]; /* zz or uu */
  P_vel[1]=P_vel[3]=rtk->sol.qrv[3]; /* xy or en */
  P_vel[5]=P_vel[7]=rtk->sol.qrv[4]; /* yz or nu */
  P_vel[2]=P_vel[6]=rtk->sol.qrv[5]; /* zx or ue */
  covenu(llh_pos,P_vel,Qvenu);

  /* Velocity covariance in NED */
  gnss_vel_ned_cov[0] = Qvenu[4];
  gnss_vel_ned_cov[1] = Qvenu[0];
  gnss_vel_ned_cov[2] = Qvenu[8];

  /*  Map information ********
  // Fetching data to the first buffer structure (lane.buffer)
  //inibuff();
  // Closest point computation
  //clstpt(xyz_ini_pos);
  for (j=0;j<3;j++) xyz_prev_clst_pos[j]=lane.buffer[j];
  // Heading (body-frame) or Azimuth (local-nav-frame)? constrain using lanes
  head_angle=headfrommap(xyz_prev_clst_pos);
  printf("Heading from map: %lf deg\n",head_angle*R2D);
 */

 /* Processing window - integrates when GNSS and INS mea. are closer by 0.1s
 The do while takes care when INS or GNSS is too ahead from each other         */
 do {

   printf("NEXT INS MEASUREMENT READING  **********************************\n");


   memset(&imu_curr_meas, 0, sizeof(um7pack_t));

    /*  Low-cost imu  */
   /* Read imu current measurements
   while(imu_curr_meas.status!=1){
     if(!fgets(datastring, 150, fimu)) {printf("END OF FILE?\n");check=NULL;break;}
     puts(datastring);
     parseimudata(datastring,&imu_curr_meas);
   }

   printf("TIMS FROM INS: %lf, %f\n", imu_curr_meas.sec, imu_curr_meas.internal_time);

   /* Turn-on initial biases             //exp4
 imu_curr_meas.a[0]=imu_curr_meas.a[0]-0.4453528;//-0.455972;//*G;
 imu_curr_meas.a[1]=imu_curr_meas.a[1]-1.162248;//-1.02027;//-0.972083*G;
 imu_curr_meas.a[2]=imu_curr_meas.a[2];//+8.0764127;//+8.1495;//-1.64328*G;

   // Fixing imu time
   imu_curr_meas.sec=imu_curr_meas.sec+ \
   ((double)imu_curr_meas.internal_time-floor((double)imu_curr_meas.internal_time));   
*/
   if (imu_obs_prev.sec <= 0.0) {
     rewind(imu_tactical); 
   }   

   /* Tactical IMU reading  */
   // For March experiments: 2,1,0 - For previous: 2, 0, 1
     check=fgets(str, 100, imu_tactical);
     sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &imu_curr_meas.sec, &imu_curr_meas.a[2],\
     &imu_curr_meas.a[1],&imu_curr_meas.a[0], &imu_curr_meas.g[2],&imu_curr_meas.g[1],\
     &imu_curr_meas.g[0]);

     if (check==NULL) {
         /* end of file */
          printf("END OF INS FILE: %s\n", check);
          //rewind(imu_tactical);
         return;
       }else{imu_curr_meas.status=1;} 
        
     
     /* Axis mirroring - For MArch experiments  */
     //imu_curr_meas.a[1]=-imu_curr_meas.a[1];
     //imu_curr_meas.a[2]=-imu_curr_meas.a[2];
     //imu_curr_meas.g[1]=-imu_curr_meas.g[1]; 
     //imu_curr_meas.g[2]=-imu_curr_meas.g[2]; 

    /* raw acc. to m/s*s and rate ve locity from degrees to radians  */
    for (j= 0;j<3;j++) imu_curr_meas.a[j]=imu_curr_meas.a[j]*G;
    for (j=0;j<3;j++) imu_curr_meas.g[j]=imu_curr_meas.g[j]*D2R;

    /* IMU was mounted with Z axis towards up */
    //imu_curr_meas.a[2] = -imu_curr_meas.a[2];

    /* Turn-on initial biases             //city collection  
  imu_curr_meas.a[0]=imu_curr_meas.a[0]-0.4453528;
  imu_curr_meas.a[1]=imu_curr_meas.a[1]-1.162248;;
  imu_curr_meas.a[2]=imu_curr_meas.a[2];*/ 


   /* Output raw INS */
   fprintf(out_raw_fimu, "%lf %lf %lf %lf %lf %lf %lf\n",imu_curr_meas.sec, \
   imu_curr_meas.a[0],\
   imu_curr_meas.a[1],imu_curr_meas.a[2], imu_curr_meas.g[0],imu_curr_meas.g[1],
   imu_curr_meas.g[2]);

   /* If INS time is ahead of GNSS exit INS loop */
   printf("GPS AND IMU: %lf %lf\n", ini_pos_time, imu_curr_meas.sec);
  if (ini_pos_time-imu_curr_meas.sec < -0.01) { // for tactical, -1.65 for consurmer
    /* Back INS file reader */
    printf("Coarse.GNSS time behind imu time!! INS stays in the same epoch\n");
    rewind(imu_tactical);
    //imufileback();

    if(imu_obs_prev.sec > 0.0){
     /* Not the first solution, use previous sol.*/
   }else{
     /* The first solution */
     for (j=0;j<3;j++) PVA_prev_sol.re[j]=xyz_ini_pos[j];
     ecef2pos(PVA_prev_sol.re,PVA_prev_sol.r);
     for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
     if (TC_or_LC) {
       /* Tightly*/
     }
   }

   for (i = 0; i < 17; i++) {
     for (j=0;j<17;j++) (i==j?PVA_prev_sol.P[i*17+j]=0.01:0.0);
   }

    /* Attitude will not be actually the previous since it hasn't
    been initialized yet !! */
    PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
    PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;
    PVA_prev_sol.sec=imu_curr_meas.sec;
    PVA_prev_sol.time=imu_curr_meas.time;
    PVA_prev_sol.Nav_or_KF=0; /* Use nav. since it is not integrated */

    break; 
  }else{

    /* Interpolate INS measurement to match GNSS time  */
    if (imu_obs_prev.sec>0.0) {
      if (ini_pos_time < imu_curr_meas.sec && ini_pos_time > imu_obs_prev.sec) {
        printf("INTERPOLATE INS MEAS TO GNSS TIME \n");
        /* imu_curr_meas is modified  */
        IMU_meas_interpolation(imu_obs_prev.sec, imu_curr_meas.sec, ini_pos_time, \
        &imu_curr_meas, &imu_obs_prev);
        /* Return imu file pointer to process the imu_curr in the next epoch */
       //imufileback();
     }
   }

  printf("GNSS.INSPREV.INS.TIMES.DIFF: %lf, %lf, %lf, %lf\n",ini_pos_time,imu_obs_prev.sec,\
  imu_curr_meas.sec, ini_pos_time-imu_curr_meas.sec );

  /* Attitude Initialization (Groves, 2013)  */

  /* Earth rotation vector in e-frame	*/
  wiee[0]=0;wiee[1]=0;wiee[2]=OMGE;

  /* local apparent gravity vector */
  appgrav(llh_pos, gan, wiee);

  /* Coarse alignment or use of previous solution? */
  /* Levelling and gyrocompassing, update imu->aea[] vector */
  printf("\n Coarse.obs:               %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", \
    imu_curr_meas.sec, imu_obs_prev.fb[0],imu_obs_prev.fb[1],imu_obs_prev.fb[2], \
  imu_curr_meas.a[0],imu_curr_meas.a[1],imu_curr_meas.a[2]);

  /* Azimuth initialization with gnss */
  double delta_xyz_prev[3], delta_xyz_prev_curr[3];
  for (j=0;j<3;j++) delta_xyz_prev[j]=0.0; //difference to itself, since it is the origin
  for (j=0;j<3;j++) delta_xyz_prev_curr[j]=xyz_ini_pos[j]-pvagnss.r[j];//PVA_prev_sol.r[j];
  ecef2enu(llh_pos, delta_xyz_prev,enu_prev);
  ecef2enu(llh_pos, delta_xyz_prev_curr, enu_curr);
  azmt=azt(enu_prev, enu_curr);  /*azimuth in radians */

  if(imu_obs_prev.sec > 0.0){
    /* Not the first solution */

    if (PVA_prev_sol.Nav_or_KF) {
      /* If coming from an integration */
      //coarseAlign(imu_obs_prev.fb, imu_obs_prev.wibb, PVA_prev_sol.v, gan,\
         PVA_prev_sol.A); // it modifies PVA_prev_sol.A[] vector */

      /* Azimuth from GNSS */
      //PVA_prev_sol.A[2]=azmt;

    }else{
      /* From navigation -> uses previous PVA_old
      coarseAlign(imu_obs_prev.fb, imu_obs_prev.wibb, PVA_prev_sol.v, gan,\
         PVA_prev_sol.A); // it modifies PVA_prev_sol.A[] vector */
    }

  //  printf("Coarse.NOT the First sol.: %lf, %lf, %lf, %lf \n", \
  //    imu_curr_meas.sec, PVA_prev_sol.A[0],PVA_prev_sol.A[1], PVA_prev_sol.A[2]);
  }else{
    /* The first solution*/
      if (norm(PVA_prev_sol.v,3)>0.0) {
        /* It was already initialized by GNSS*/
      }else{
        /*Use current gnss velocity as an approximation */
        for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
      }
      /* Initialize velocity with GNSS for coarseAlignment */
      for (j=0;j<3;j++) imu_curr_meas.v[j]=PVA_prev_sol.v[j];
      coarseAlign(imu_curr_meas.a, imu_curr_meas.g, PVA_prev_sol.v, gan, \
        imu_curr_meas.aea);
      for (j=0;j<3;j++) PVA_prev_sol.A[j]=imu_curr_meas.aea[j]; //Attitude in A_nb frame
      printf("Coarse.First sol.:      %lf, %lf, %lf, %lf \n", \
      imu_curr_meas.sec, PVA_prev_sol.A[0],PVA_prev_sol.A[1],PVA_prev_sol.A[2]);

      /* Azimuth from GNSS */
      //PVA_prev_sol.A[2]=azmt;
    }

    //PVA_prev_sol.A[2]=imu_curr_meas.aea[2]=head_angle;

/* Tightly or Loosley-coupled integration */
/* If enough GNSS measurements uses Tightly otherwise Loosley-coupled*/
  if (TC_or_LC) {

    /* Tightly coupled INS/GNSS */
    if(imu_obs_prev.sec > 0.0){
      printf("Coarse.************************************ First NAV. solution\n");
      TC_INS_GNSS_core(rtk, obs, n, nav, &imu_curr_meas, (double)imu_curr_meas.sec\
      - imu_obs_prev.sec, gnss_ned_cov, gnss_vel_ned_cov, &PVA_prev_sol, Tact_or_Low_IMU, TC_or_LC);
    }else{
      for (j=0;j<3;j++) PVA_prev_sol.re[j]=xyz_ini_pos[j];
      ecef2pos(PVA_prev_sol.re,PVA_prev_sol.r);
      for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
      for (j=0;j<3;j++) PVA_prev_sol.A[j]=imu_curr_meas.aea[j];
      PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
      PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;
      PVA_prev_sol.sec=imu_curr_meas.sec;
      PVA_prev_sol.time=imu_curr_meas.time;
    }

  }else{

    /* Loosley-coupled INS/GNSS */
    if(imu_obs_prev.sec > 0.0){
      //insgnssLC
      LC_INS_GNSS_core(xyz_ini_pos, gnss_xyz_ini_cov, gnss_ned_cov, ned_ini_vel, gnss_vel_ned_cov, ini_pos_time, \
      &imu_curr_meas, (double)imu_curr_meas.sec - imu_obs_prev.sec,\
      &PVA_prev_sol, &imu_obs_prev, Tact_or_Low_IMU, TC_or_LC);
    }else{
      for (j=0;j<3;j++) PVA_prev_sol.re[j]=xyz_ini_pos[j];
      ecef2pos(PVA_prev_sol.re,PVA_prev_sol.r);
      for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
      for (j=0;j<3;j++) PVA_prev_sol.A[j]=imu_curr_meas.aea[j];
      PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
      PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;
      PVA_prev_sol.sec=imu_curr_meas.sec;
      PVA_prev_sol.time=imu_curr_meas.time;
    }
  }
  printf("PVA.v: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.v[0],PVA_prev_sol.ve[1],PVA_prev_sol.ve[2]);
  printf("PVA.A: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.A[0],PVA_prev_sol.A[1],PVA_prev_sol.A[2]);
  printf("PVA.r: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.re[0],PVA_prev_sol.re[1],PVA_prev_sol.re[2]);

  if(imu_obs_prev.sec > 0.0 )   //&& PVA_prev_sol.Nav_or_KF==1){
  {
  /* Non-holonomic constraints  */
  nhc1(&PVA_prev_sol, &imu_curr_meas, 17);

  }

  printf("PVA.v: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.ve[0],PVA_prev_sol.ve[1],PVA_prev_sol.ve[2]);
  printf("PVA.A: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.A[0],PVA_prev_sol.A[1],PVA_prev_sol.A[2]);
  printf("PVA.r: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.re[0],PVA_prev_sol.re[1],PVA_prev_sol.re[2]);

  

  /* ZUPT detection, based on   GREJNER-BRZEZINSKA et al. (2002) */
   /* Condition */
   if ( (fabs(PVA_prev_sol.v[0]) - vn0) <= 3*vn0std && (fabs(PVA_prev_sol.v[1]) - ve0) <= 3*ve0std ) {
     printf("Static.zupt.vel: %lf, Dif|v-v0|= %lf <= %lf and %lf < %lf \n",PVA_prev_sol.sec, \
   (fabs(PVA_prev_sol.v[0]) - vn0), 3*vn0std, (fabs(PVA_prev_sol.v[1]) - ve0), 3*ve0std );
   //for (j=0;j<3;j++) PVA_prev_sol.v[j]=0.0;
   /* Zero velocity update */
   if (zvu_counter++>10) {
     printf("ZVU UPDATE: counter: %d", zvu_counter);
     zvu1(&PVA_prev_sol, &imu_curr_meas, 17);
     zvu_counter=0;
   }
   }else{
     zvu_counter=0;
   }

  /* Output PVA solution */
  if (PVA_prev_sol.sec<0.0) {
  }else{
    fprintf(out_PVA,"%lf %.12lf %.12lf %lf %lf %lf %lf %lf %lf %lf %d\n",\
    PVA_prev_sol.sec, PVA_prev_sol.r[0]*R2D, PVA_prev_sol.r[1]*R2D, PVA_prev_sol.r[2],\
    PVA_prev_sol.v[0], PVA_prev_sol.v[1], PVA_prev_sol.v[2],
    PVA_prev_sol.A[0]*R2D,PVA_prev_sol.A[1]*R2D,PVA_prev_sol.A[2]*R2D, PVA_prev_sol.Nav_or_KF );
  }



/* GYRO static detection
   if ( (fabs(imu_curr_meas.g[0]) - gyrx0 <= 3*gyrx0std) && (fabs(imu_curr_meas.g[1]) - gyry0) <= 3*gyry0std ) {
     printf("Static.zupt.gyr: %lf, Dif|g-g0|=%lf <= %lf and %lf < %lf \n",PVA_prev_sol.sec, \
   (fabs(imu_curr_meas.g[0]) - gyrx0), 3*gyrx0std, (fabs(imu_curr_meas.g[1]) - gyry0), 3*gyry0std );
   }
*/


  /* Update measurements */
  for (j=0;j<3;j++){
  imu_obs_prev.fb[j]=imu_curr_meas.a[j];
  imu_obs_prev.wibb[j]=imu_curr_meas.g[j];
  }
  imu_obs_prev.sec=imu_curr_meas.sec;
  imu_obs_prev.time=imu_curr_meas.time;


} // end if else If INS time is ahead of GNSS exit INS loop condition

   /* Loop conditions */
   if (ini_pos_time-imu_curr_meas.sec<0.0001) {
     printf("UPDATE INS AND GNSS \n");
     break;
   }else{
     if (ini_pos_time-imu_curr_meas.sec<-0.0001) {
       printf("UPDATE INS AND GNSS \n");
       break;
     }else{
       if (check==NULL) {
         /* end of file */
          printf("END OF INS FILE: %s\n", check);
         return;
       }else{printf("KEEP UPDATING INS \n");}
     }
   }

 } while(ini_pos_time-imu_curr_meas.sec > 0.0001 || \
   ini_pos_time-imu_curr_meas.sec > -0.0001 );

 printf("OUT OF LOOP\n \
 IMU_TIME_AND_GNSS_TIMES: %f %f\n", imu_curr_meas.sec, ini_pos_time);

  // printf("Vel NED GNSS: %lf, %lf, %lf\n",ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2]);
  // printf("Vel NED INS: %lf, %lf, %lf\n",PVA_prev_sol.v[0],PVA_prev_sol.v[1],PVA_prev_sol.v[2]);

   /* Initialize position from GNSS SPP */
   if (ini_pos_time-imu_curr_meas.sec < -0.01) { //-0.01
   /* INS time ahead of GPS, don't pass anything to the satellite positioning */
    }else{
   /* INS and GNSS times are synchronized */
   /* State and covariances to satellite positioning */
   //for (j=0;j<3;j++) rtk->sol.rr[j]=PVA_prev_sol.re[j];
   if (PVA_prev_sol.out_errors[6]<=0.0 || PVA_prev_sol.out_errors[7]<=0.0 || \
    PVA_prev_sol.out_errors[8]<=0.0 ) {
   }else{
       //for (j=0;j<3;j++) rtk->sol.qr[j]=(1.0/PVA_prev_sol.out_errors[j+6]);
     }
   }
   /* Making GNSS position and velocity be the INS initialization on the next
    processing, instead of the Navigation-derived one */
printf("NED_ini_vel.: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2]);
printf("NED_ini_vel.: t:%lf - %lf %lf %lf\n", PVA_prev_sol.sec, PVA_prev_sol.v[0],PVA_prev_sol.v[1],PVA_prev_sol.v[2]);
if (gnss_vel) {
  /* It came with velocities */
  /*
  if (norm(ned_ini_vel,3) > 0.5) { // ned_ini_vel or PVA_prev_sol. ?
      printf("IS MOVING\n");
    }else{printf("IS STATIC\n");
      for (j=0;j<3;j++) ned_ini_vel[j]=0.0;
    }*/

    if ( (fabs(ned_ini_vel[0]) - vn0) <= 3*vn0std && (fabs(ned_ini_vel[1]) - ve0) <= 3*ve0std ) {
      printf("Static.zupt.ned_ini: %lf, Dif|v-v0|= %lf <= %lf and %lf < %lf \n",PVA_prev_sol.sec, \
    (fabs(PVA_prev_sol.v[0]) - vn0), 3*vn0std, (fabs(PVA_prev_sol.v[1]) - ve0), 3*ve0std );
    //for (j=0;j<3;j++) ned_ini_vel[j]=0.0;

    /* Zero velocity update */
    zvu1(&PVA_prev_sol, &imu_curr_meas, 17);

    }


}else{
  /* It didn't come with velocities */
  for (j=0;j<3;j++) ned_ini_vel[j]=PVA_prev_sol.v[j];
}


/*
      if ( (fabs(PVA_prev_sol.v[0]) - vn0) <= 3*vn0std*10 && (fabs(PVA_prev_sol.v[1]) - ve0) <= 3*ve0std*10 ) {
        printf("Static.zupt.ned_ini: %lf, Dif|v-v0|= %lf <= %lf and %lf < %lf \n",PVA_prev_sol.sec, \
      (fabs(PVA_prev_sol.v[0]) - vn0), 3*vn0std, (fabs(PVA_prev_sol.v[1]) - ve0), 3*ve0std );
      for (j=0;j<3;j++) ned_ini_vel[j]=0.0;
      }
*/


/*
if ( norm(ned_ini_vel,3) <= 0.000001 && norm(PVA_prev_sol.v,3) <= 0.000001) {
  for (j=0;j<3;j++) ned_ini_vel[j]=PVA_prev_sol.v[j]=0.0;
}else{
  if (norm(ned_ini_vel,3) > 0.0005 && PVA_prev_sol.Nav_or_KF==1) {
    /* Use it for next
  }else{
    if (norm(PVA_prev_sol.v,3) > 0.005 && PVA_prev_sol.Nav_or_KF==1) {
      for (j=0;j<3;j++) ned_ini_vel[j]=PVA_prev_sol.v[j];
    }
  }
}
*/

    if (PVA_prev_sol.Nav_or_KF) {
      /* Integrated solution */
      /* Zero velocity update */
      //zvu(&PVA_prev_sol, &imu_curr_meas, 17);

      for (j=0;j<3;j++) PVA_prev_sol.re[j]=xyz_ini_pos[j];
      ecef2pos(PVA_prev_sol.re,PVA_prev_sol.r);
      for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
      PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
      PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;

      PVA_prev_sol.t_s = time2gpst(rtk->sol.time,NULL); /* time of week in (GPST) */
      for (j=0;j<3;j++) pvagnss.r[j]=xyz_ini_pos[j];
      pvagnss.sec=PVA_prev_sol.t_s; 

    }else{
      /* Navigation solution */
      //for (j=0;j<3;j++) PVA_prev_sol.re[j]=xyz_ini_pos[j];
      //for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
      PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
      PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;
      PVA_prev_sol.t_s = time2gpst(rtk->sol.time,NULL); /* time of week in (GPST) */
      for (j=0;j<3;j++) pvagnss.r[j]=xyz_ini_pos[j];
      pvagnss.sec=PVA_prev_sol.t_s;
    }

    printf("AZIMUTH: %lf: %lf\n", time2gpst(rtk->sol.time,NULL), azmt*R2D);

   /* Update */

   imu_obs_global=imu_obs_prev;
   pva_global=PVA_prev_sol;


 printf("\n *****************  CORE ENDS ***********************\n");
}

/* input imu measurement data-------------------------------------------------*/
static int inputimu(ins_states_t *ins, int week){
  char check, str[150];
  um7pack_t imu_curr_meas={0};
  int i, j;

  if(insgnssopt.Tact_or_Low){
    /* Tactical KVH input */
    // For March experiments: 2,1,0 - For previous: 2, 0, 1
    check=fgets(str, 150, imu_tactical);
    sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &ins->time, &ins->data.fb0[2],\
    &ins->data.fb0[1],&ins->data.fb0[0], &ins->data.wibb0[2],&ins->data.wibb0[1],\
    &ins->data.wibb0[0]);

    ins->data.time=gpst2time(week, ins->time);

    /* raw acc. to m/s*s and rate ve locity from degrees to radians  */
    for (j= 0;j<3;j++) ins->data.fb0[j]=ins->data.fb0[j]*Gcte;
    for (j=0;j<3;j++) ins->data.wibb0[j]=ins->data.wibb0[j]*D2R;

    if (check==NULL) {
      /* end of file */
      printf("END OF INS FILE: %s\n", check);
      //rewind(imu_tactical);
      return 0;
    }else{
      imu_curr_meas.status=1;
      return 1;
      }

  }else{
    /* Low cost IMU um7 input*/
   while(imu_curr_meas.status!=1){
     if(!fgets(str, 150, fimu)) {
       printf("END OF FILE?\n");
       check=NULL;
       return 0;
       }
      parseimudata(str,&imu_curr_meas); 

    // Fixing imu time
   imu_curr_meas.sec=imu_curr_meas.sec+ \
   ((double)imu_curr_meas.internal_time\
   -floor((double)imu_curr_meas.internal_time));
   }
   ins->time=imu_curr_meas.sec;
   for (i=0;i<3;i++) ins->data.fb0[i]=imu_curr_meas.a[i];
   for (i=0;i<3;i++) ins->data.wibb0[i]=imu_curr_meas.g[i];
   ins->data.time=gpst2time(week, imu_curr_meas.sec);

   return 1; 
  } 
    
}

/* Stores gnss measurement and observation to structure buffer by time */
void gnssbuffer(sol_t *sol, obsd_t *data, int n){
  int i, stat=0;
  //printf("GNSS_meas: iter: %d window: %d\n", gnss_w_counter, insgnssopt.gnssw);

  if(gnss_w_counter >= insgnssopt.gnssw ){
     for(i=0;i<insgnssopt.gnssw-1;i++){
      solw[i]=solw[i+1];
    }
     for (i=0;i<n;i++) obsw.data0[i]=obsw.data1[i];
     obsw.n0=obsw.n1;
     for (i=0;i<n;i++) obsw.data1[i]=obsw.data2[i];
     obsw.n1=obsw.n2;
     for (i=0;i<n;i++) obsw.data2[i]=data[i];
     obsw.n2=n;
    solw[insgnssopt.gnssw-1]=*sol;
  }else{
    solw[gnss_w_counter]=*sol;
    if(gnss_w_counter==0) {for (i=0;i<n;i++) obsw.data0[i]=data[i];
       obsw.n0=n;}
    else if(gnss_w_counter==1) {for (i=0;i<n;i++) obsw.data1[i]=data[i];
    obsw.n1=n;}
    else {for (i=0;i<n;i++) obsw.data2[i]=data[i]; 
    obsw.n2=n;}
  }
}

void insbufferpass(ins_states_t *insin, ins_states_t *insout){
  int i,j;

  matcpy(insin->x,insout->x,insout->nx,1);
  matcpy(insin->P,insout->P,insout->nx,insout->nx);
  matcpy(insin->P0,insout->P0,insout->nx,insout->nx);

}


/* Stores ins measurement to structure buffer by time */
void insbuffer(ins_states_t *insc){
  int i;
  //printf("INS_meas: iter: %d window: %d\n", ins_w_counter, insgnssopt.insw);
  if(ins_w_counter >= insgnssopt.insw ){
    for(i=0;i<insgnssopt.insw-1;i++){
      insw[i]=insw[i+1];
      insbufferpass(insw+i, insw+i+1);
    }
    insw[insgnssopt.insw-1]=*insc; 
    insbufferpass(insc, insw+insgnssopt.insw-1);
  }else{
    ins_buffinit(insw+ins_w_counter, insc->nx);
    insw[ins_w_counter]=*insc;
    insbufferpass(insc, insw+ins_w_counter);
  }
}

/* conversion matrix of ned frame to ecef frame--------------------------------
 * conversion matrix of ned to ecef
 * args   : double *pos     I   position {lat,lon,height} (rad/m)
 *          double *Cne     O   convertion matrix between frame
 * return : none
 * --------------------------------------------------------------------------*/
extern void ned2xyz(const double *pos,double *Cne)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);

    Cne[0]=-sinp*cosl; Cne[3]=-sinl; Cne[6]=-cosp*cosl;
    Cne[1]=-sinp*sinl; Cne[4]= cosl; Cne[7]=-cosp*sinl;
    Cne[2]=      cosp; Cne[5]= 0.0;  Cne[8]=-sinp;
}
  
/* estimate heading using the navigation frame velocity-----------------------
 * args   :  double *vel  I  velocity in navigation frame (ned-frame)
 * return :  heading (rad) 
 * --------------------------------------------------------------------------*/
extern double vel2head(const double *vel)
{
    return atan2(vel[1],fabs(vel[0])<1E-4?1E-4:vel[0]);
}     

/* ned frame to body frame---------------------------------------------------
 * args    : double rpy      I  attitude {roll,picth,yaw} (rad)
 *           double Cnb      O  matrix of ned to body frame
 * return  : none
 * -------------------------------------------------------------------------*/
extern void rpy2dcm(const double *rpy,double *Cnb)
{
    double sin_phi=sin(rpy[0]),cos_phi=cos(rpy[0]),
           sin_theta=sin(rpy[1]),cos_theta=cos(rpy[1]),
           sin_psi=sin(rpy[2]),cos_psi=cos(rpy[2]);
    Cnb[0]= cos_theta*cos_psi;
    Cnb[3]= cos_theta*sin_psi;
    Cnb[6]=-sin_theta;
    Cnb[1]=-cos_phi*sin_psi+sin_phi*sin_theta*cos_psi;
    Cnb[4]= cos_phi*cos_psi+sin_phi*sin_theta*sin_psi;
    Cnb[7]= sin_phi*cos_theta;
    Cnb[2]= sin_phi*sin_psi+cos_phi*sin_theta*cos_psi;
    Cnb[5]=-sin_phi*cos_psi+cos_phi*sin_theta*sin_psi;
    Cnb[8]= cos_phi*cos_theta;
}

/* multiply 3d matries -------------------------------------------------------*/
extern void matmul3(const char *tr, const double *A, const double *B, double *C)
{
    matmul(tr,3,3,3,1.0,A,B,0.0,C);
}
/* multiply 3d matrix and vector ---------------------------------------------*/
extern void matmul3v(const char *tr, const double *A, const double *b, double *c)
{
    char t[]="NN";
    t[0]=tr[0];
    matmul(t,3,1,3,1.0,A,b,0.0,c);
}
/* 3d skew symmetric matrix --------------------------------------------------*/ 
extern void skewsym3(const double *ang, double *C)
{
    C[0]=0.0;     C[3]=-ang[2]; C[6]=ang[1];
    C[1]=ang[2];  C[4]=0.0;     C[7]=-ang[0];
    C[2]=-ang[1]; C[5]=ang[0];  C[8]=0.0;
}
/* 3d skew symmetric matrix --------------------------------------------------*/ 
extern void skewsym3x(double x,double y,double z,double *C)
{
    C[0]=0.0; C[3]=-z;  C[6]=y;
    C[1]=z;   C[4]=0.0; C[7]=-x;
    C[2]=-y;  C[5]=x;   C[8]=0.0;
}
/* set all matrix elements to zero ------------------------------------------*/
extern void setzero(double *A,int n,int m)
{
    int i,j;
    for (i=0;i<n;i++) for (j=0;j<m;j++) A[i*m+j]=0.0;
}
/* D=A*B*C---------------------------------------------------------------------*/
extern void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D)
{
    char tr_[8];
    double *T=mat(n,q);
    matmul(tr,n,q,p,1.0,A,B,0.0,T);
    sprintf(tr_,"N%c",tr[2]);
    matmul(tr_,n,m,q,1.0,T,C,0.0,D); free(T);
}

/* gnss antenna position/velecity transform to ins position/velecity---------
 * args  :  double *pos    I  position of gnss antenna (ecef)
 *          double *vel    I  velecity of gnss antenna (ecef)
 *          double *Cbe    I  transform matrix of body-frame to ecef-frame
 *          double *lever  I  arm lever
 *          imud_t *imu    I  imu measurement data
 *          double *posi   O  ins position (ecef) (NULL :no output)
 *          double *veli   O  ins velecity (ecef) (NULL :no output)
 * return : none
 * -------------------------------------------------------------------------*/
extern void gapv2ipv(const double *pos,const double *vel,const double *Cbe,
                     const double *lever,const imuraw_t *imu,double *posi,
                     double *veli)
{
    int i; double T[9],TT[9];

    if (posi) {
        matcpy(posi,pos,3,1);
        matmul("NN",3,1,3,-1.0,Cbe,lever,1.0,posi);
    }
    if (veli) {
        skewsym3(imu->wibb0,T);
        matmul("NN",3,1,3,1.0,T,lever,0.0,TT);
        matmul("NN",3,1,3,1.0,Cbe,TT,0.0,T);

        matmul33("NNN",Omge,Cbe,lever,3,3,3,1,TT);
        for (i=0;i<3;i++) veli[i]=vel[i]-T[i]+TT[i];
    }
}

/* give gps antenna position/velocity to initialize ins states when static ---
 * args    :  gtime_t time     I  time of antenna position/velocity
 *            double *rr       I  gps antenna position (ecef)
 *            double *vr       I  gps antenna velocity (ecef)
 *            ins_states_t *insc  O  initialed ins states
 * return  : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int coarse_align(gtime_t time,const double *rr,const double *vr,
                     ins_states_t *insc)
{
    double llh[3],vn[3],C[9],rpy[3]={0}, sineyaw, coseyaw; 
    int i;

    trace(3,"coarse_align: time=%s\n",time_str(time,4));

    /* velocity in navigation frame */
    ecef2pos(rr,llh);
    ned2xyz(llh,C);
    matmul("TN",3,1,3,1.0,C,vr,0.0,vn);

    matcpy(insc->rn,llh,1,3);
    matcpy(insc->vn,vn ,1,3);

    /* Levelling from Groves (2013) */
    rpy[0]= atan2(-insc->data.fb0[1],-insc->data.fb0[2]); /* roll - phi angle */  
    rpy[1]= atan (insc->data.fb0[0]/sqrt(insc->data.fb0[1]*insc->data.fb0[1]+insc->data.fb0[2]*insc->data.fb0[2])); /* pitch - theta angle */

    /* Gyrocompasing by Groves(2013,page 199)  */
    /* yaw - psi angle */
    sineyaw=-insc->data.wibb0[1]*cos(rpy[0])+\
    insc->data.wibb0[2]*sin(rpy[0]);
    coseyaw= insc->data.wibb0[0]*cos(rpy[1])+\
    insc->data.wibb0[1]*sin(rpy[0])*sin(rpy[1])+\
    insc->data.wibb0[2]*cos(rpy[0])*sin(rpy[1]);
    rpy[2]=atan2(sineyaw,coseyaw);  /* yaw */

    for(i=0;i<3;i++) insc->pan[i]=rpy[i];

    rpy2dcm(rpy,C);
    /* Row to column order Or matrix transpose */
    row_to_column_order(C,3,3,insc->pCbn);

    ned2xyz(llh,C);
    matmul("NN",3,3,3,1.0,C,insc->pCbn,0.0,insc->pCbe);

    /* find closest imu measurement index */
    if (&insc->data) {
      if (fabs(timediff(time,insc->data.time))>0.005) { 
        //DTTOL:tolerance of time difference
          printf("Time tolerance limit error: %lf \n", \
          fabs(timediff(time,insc->data.time)));
         return 0;
      }
    } 

    /* initial ins position */
    gapv2ipv(rr,vr,insc->pCbe,insgnssopt.lever,&insc->data, \
             insc->pre,insc->pve);
    return 1;
}

/* give gps antenna position/velocity to initial ins states-----------------
 * args    :  gtime_t time     I  time of antenna position/velocity
 *            double *rr       I  gps antenna position (ecef)
 *            double *vr       I  gps antenna velocity (ecef)
 *            ins_states_t *insc  O  initialed ins states
 * return  : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int ant2inins(gtime_t time,const double *rr,const double *vr,
                     ins_states_t *insc)
{
    double llh[3],vn[3],C[9],rpy[3]={0};
    int i;

    trace(3,"ant2inins: time=%s\n",time_str(time,4));

    /* velocity in navigation frame */
    ecef2pos(rr,llh);
    ned2xyz(llh,C);
    matmul("TN",3,1,3,1.0,C,vr,0.0,vn);

    matcpy(insc->prn,llh,1,3);
    matcpy(insc->pvn,vn ,1,3);

    rpy[2]=vel2head(vn); /* yaw */

    for(i=0;i<3;i++) insc->pan[i]=rpy[i];

    rpy2dcm(rpy,C);
    /* Row to column order Or matrix transpose */
    row_to_column_order(C,3,3,insc->pCbn);

    ned2xyz(llh,C);
    matmul("NN",3,3,3,1.0,C,insc->pCbn,0.0,insc->pCbe);

    /* find closest imu measurement index */
    if (&insc->data) {

      if (fabs(timediff(time,insc->data.time))>0.005) { 
        // DTTOL:tolerance of time difference
         printf("Time tolerance limit error: %lf \n",\
          fabs(timediff(time,insc->data.time)));
         return 0;
      } 
      /* align imu measurement data */
      if (norm(insc->data.wibb0,3)>(10*D2R)) {
        printf("Car rotating: %lf \n", norm(insc->data.wibb0,3));
        return 0;} 
        //MAXROT: 10*D2R max rotation of vehicle velocity matching alignment */
    }
    /* initial ins position */
    gapv2ipv(rr,vr,insc->pCbe,insgnssopt.lever,&insc->data,\ 
             insc->pre,insc->pve);
    return 1;
}


/* use rtk solutions to initial ins states ------------------------------------
 * args   :  sol_t *solw      I  solution buffer
 *           insstates_t *ins IO ins states
 *           insgnssopt_t ingssopt I ins/gnss options
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
int init_inspva(sol_t *sol, ins_states_t *insc){
  int i,j, k, n=insgnssopt.gnssw;
  double dt[n-1],rr[3],vr[3];

  for (i=0;i<n-1;i++){
    dt[i]=timediff(sol[i+1].time,sol[i].time);
  }
  /* check time continuity */
  if (norm(dt,n-1)>SQRT(n-1)+0.1
    ||norm(dt,n-1)==0.0) {
    return 0;
  }
  k=(n-1)/2;
  matcpy(rr,sol[k+1].rr+0,1,3);
  matcpy(vr,sol[k+1].rr+3,1,3);
  if (norm(vr,3)==0.0) {
    for (i=0;i<3;i++) {
      vr[i]=(sol[k+1].rr[i]-sol[k].rr[i])/dt[k];
    }
  }
  if (norm(vr,3)<5.0) { 
    // 5.0: min velocity for ins velocity match alignment
    printf("Velocity error: %lf\n", norm(vr,3));

    if (!coarse_align(sol[k+1].time,rr,vr,insc)) {
      printf("coarse_align error\n");
      return 0;
    }else{
      return 1;
    }
  }

  /* initial ins state use single positioning */
  if (!ant2inins(sol[k+1].time,rr,vr,insc)) {
      printf("ant2inins error\n");
     return 0;
  }
  return 1;             
                                     
}

/* Update ins states and measurements */
void insupdt(ins_states_t *ins){
  int i,j;
   /* States update */
   ins->ptime = ins->time;
   for(i=0;i<3;i++){
    ins->pvn[i]=ins->vn[i]; ins->pve[i]=ins->ve[i]; 
    ins->prn[i]=ins->rn[i]; ins->pre[i]=ins->re[i];
    ins->pan[i]=ins->an[i]; ins->pae[i]=ins->ae[i];
   }
   for(i=0;i<9;i++){
    ins->pCbn[i]=ins->Cbn[i]; ins->pCbe[i]=ins->Cbe[i]; 
   }
   for(i=0;i<6;i++) ins->pdtr[i]=ins->dtr[i];
   ins->pdtrr=ins->dtrr;  
 
   /* Measurements update*/
   ins->pdata=ins->data;
   
}

/* Interpolate INS measurement to match GNSS time  -------------------------
*  return: 0 obs doesn't match, 1: ok                                      */
int interp_ins2gpst(double gnss_time, ins_states_t *ins)
{
  if (ins->ptime>0.0) {
   if (gnss_time <= ins->time && gnss_time > ins->ptime) {
      printf("INTERPOLATE INS MEAS TO GNSS TIME \n"); 
      /* insc is modified  */
      IMU_meas_interpolation(ins->ptime, ins->time, gnss_time, \
      ins, &ins->pdata);
      /* Return imu file pointer to process the imu_curr in the next epoch */
      //imufileback(); //only for low-cost imu
      return 1;
    }else return 0;
  }else return 0; 

}

/* initialize ins control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void insinit(ins_states_t *ins, insgnss_opt_t *insopt, prcopt_t *opt, int nsat) 
{
    int i;

    trace(3,"insinit :\n");

    ins->nb=ins->nx=insgnssopt.mode<1?15:ppptcnx(opt);
    //ins->nb=ins->nx=insgnssopt.mode<1?15:(xnRx(opt)+nsat);
   // ins->nb=opt->mode<=PMODE_FIXED?NR(opt):0; //what is its use??  
    ins->dt=0.0;
    ins->x=zeros(ins->nx,1);
    ins->P=zeros(ins->nx,ins->nx);
    ins->P0=zeros(ins->nx,ins->nx);
   // ins->xa=zeros(ins->nb,1);
    ins->F =eye  (ins->nx); 
   // ins->Pa=zeros(ins->nb,ins->nb);

    /* initialize parameter indices */
    initPNindex(opt);

    getP0(insopt, ins->P, ins->nx);
    //getP0(insopt, ins->Pa, ins->nx);    
 
}
/* initialize buffer */
extern void ins_buffinit(ins_states_t *ins, int nx) 
{
  int i;
    trace(3,"insinit :\n");

    ins->nx=nx;
    ins->x=zeros(nx,1);    
    ins->P=zeros(nx,nx);
    ins->P0=zeros(nx,nx);
   // ins->xa=zeros(nx,1);
    ins->F =eye  (nx); 
   // ins->Pa=zeros(nx,nx);
 
}   

/* free ins control ------------------------------------------------------------
* free memory for ins control struct
* args   : ins_states_t    *ins    IO  ins control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void insfree(ins_states_t *ins)
{
    trace(3,"insfree :\n"); 

    ins->nx=ins->nb=0;ins->x =NULL;
    free(ins->x );  ins->P =NULL;
    free(ins->P );   ins->P0 =NULL;
    free(ins->P0);  ins->F = NULL;
   // free(ins->xa); ins->xa=NULL; printf("Here 3\n"); 
    free(ins->F );   
   // free(ins->Pa); ins->Pa=NULL; printf("Here 5\n");
}

extern void print_ins_pva(ins_states_t *ins){
  int i,j;

  printf("ins->re and ins->pre:\n");
  for (i=0;i<3;i++) printf("%lf ", ins->re[i]);
  printf("\n");
  for (i=0;i<3;i++) printf("%lf ", ins->pre[i]);
  printf("\n");

  printf("ins->ve and ins->pve:\n");
  for (i=0;i<3;i++) printf("%lf ", ins->ve[i]);
  printf("\n");
  for (i=0;i<3;i++) printf("%lf ", ins->pve[i]);
  printf("\n");

  printf("ins->Cbe and ins->pCbe:\n");
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++){
      printf("%lf ", ins->Cbe[i*3+j]);
    }
    printf("\n");
  }
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++){
      printf("%lf ", ins->pCbe[i*3+j]); 
    }
    printf("\n");
  }

  printf("x:\n");
  for (j = 0; j < ins->nx; j++) printf("%lf,",ins->x[j]); printf("\n");

  printf("P:\n");
  for (j = 0; j < ins->nx; j++) printf("%lf,",ins->P[j*ins->nx+j]); 
  printf("\n");

   printf("P0:\n");
  for (j = 0; j < ins->nx; j++) printf("%lf,",ins->P0[j*ins->nx+j]); 
  printf("\n");
}
extern void memset_ins_pva(ins_states_t *ins){
  int i,j;

  for (i=0;i<3;i++){
    ins->re[i]=ins->pre[i]=0.0;
    ins->ve[i]=ins->pve[i]=0.0;
    ins->rn[i]=ins->prn[i]=0.0;
    ins->vn[i]=ins->pvn[i]=0.0;
  }
   for (i=0;i<9;i++){
    ins->Cbe[i]=ins->pCbe[i]=0.0;
    ins->Cbn[i]=ins->pCbn[i]=0.0;
  }

}
/* measurement sensitive-matrix for non-holonomic----------------------------*/
static int bldnhc(const insgnss_opt_t *opt,const imuraw_t *imu,const double *Cbe,
                  const double *ve,int nx,double *v,double *H,double *R)
{
    int i,nv,IA,IV;
    double C[9],T[9],vb[3],r[2],S[9];
    double T1[9];

    trace(3,"bldnhc:\n");

    IA=xiA(); IV=xiV();

    /* velocity in body-frame */
    matmul("TN",3,1,3,1.0,Cbe,ve,0.0,vb);

    skewsym3(ve,C);
    matmul("TN",3,3,3,1.0,Cbe,C,0.0,T);
    matt(Cbe,3,3,C);

#if UPD_IN_EULER
    jacobian_prot_pang(Cbe,S);
    matcpy(T1,T,3,3);
    matmul("NN",3,3,3,1.0,T1,S,0.0,T);
#endif
    /* build residual vector */
    for (nv=0,i=1;i<3;i++) {

        /* check velocity measurement */
        if (fabs(vb[i])>MAXVEL) {
            trace(2,"too large velocity measurement\n");
            continue;
        }
        /* check gyro measurement */
        if (fabs(norm(imu->wibb,3))>30.0*D2R) {
            trace(2,"too large vehicle turn\n");
            continue;
        }
        H[IA+nv*nx]=T[i]; H[IA+1+nv*nx]=T[i+3]; H[IA+2+nv*nx]=T[i+6];
        H[IV+nv*nx]=C[i]; H[IV+1+nv*nx]=C[i+3]; H[IV+2+nv*nx]=C[i+6];
        
        v[nv  ]=vb[i];
        r[nv++]=VARVEL;
    }
    for (i=0;i<nv;i++) R[i+i*nv]=r[i];
    return nv;
}
/* using non-holonomic constraint for ins navigation---------------------------
 * args    :  insstate_t* ins  IO  ins state
 *            insopt_t* opt    I   ins options
 *            imud_t* imu      I   imu measurement data
 * return  : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int nhc(ins_states_t *ins,const insgnss_opt_t *opt)
{
    const imuraw_t *imu=&ins->data;
    int nx=ins->nx,info=0,nv, i;
    double *H,*v,*R,*x;  

    trace(3,"nhc:\n");
    printf("nhc:\n");

    H=zeros(2,nx); R=zeros(2,2);
    v=zeros(2,1); x=zeros(1,nx);
    for (i = 0; i < nx; i++) x[i]=1E-17;

    nv=bldnhc(opt,imu,ins->Cbe,ins->ve,nx,v,H,R);
    if (nv>0) {

        /* kalman filter */
        info=filter(x,ins->P,H,v,R,nx,nv);
        printf("nhc.x:\n");
        for (i=0; i < nx; i++){
          printf("%lf ", x[i]);
        }

        /*  check ok? */
        if (info) {
            trace(2,"non-holonomic constraint filter fail\n"); 
            printf("non-holonomic constraint filter fail\n");
            info=0;
        }
        else {
            /* solution ok */
            //ins->stat=INSS_NHC; 
            info=1;
            clp(ins,opt,x);
            trace(3,"use non-holonomic constraint ok\n");
            printf("use non-holonomic constraint ok\n");
        }
    }
    free(H); free(v);
    free(R); free(x);
    return info;
}
/* zero velocity update for ins navigation -----------------------------------
 * args    :  insstate_t *ins  IO  ins state
 *            insopt_t *opt    I   ins options
 *            imud_t *imu      I   imu measurement data
 *            int flag         I   static flag (1: static, 0: motion)
 * return  : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int zvu(ins_states_t *ins,const insgnss_opt_t *opt,int flag)
{
    imuraw_t *imu=&ins->data;
    int nx=ins->nx,info=0, i;
    static int nz=0;
    double *x,*H,*R,*v,I[9]={-1,0,0,0,-1,0,0,0,-1};

    trace(3,"zvu:\n");
    printf("zvu:\n");

   // flag&=nz++>MINZC?nz=0,true:false;

    if (!flag) return info;

    x=zeros(1,nx); H=zeros(3,nx);
    R=zeros(3,3); v=zeros(3,1);
    for (i = 0; i < nx; i++) x[i]=1E-17;

    /* sensitive matrix */
    asi_blk_mat(H,3,nx,I,3,3,0,3);

    /* variance matrix */
    R[0]=R[4]=R[8]=VARVEL;

    v[0]=ins->ve[0];
    v[1]=ins->ve[1];
    v[2]=ins->ve[2]; /* residual vector */

    if (norm(v,3)<MAXVEL&&norm(imu->wibb,3)<MAXGYRO) { 

        /* ekf filter */
        info=filter(x,ins->P,H,v,R,nx,3);

        printf("zvu.x:\n");
        for (i=0; i < nx; i++){
          printf("%lf ", x[i]);
        }

        /* solution fail */
        if (info) {
            trace(2,"zero velocity update filter error\n");
            printf("zero velocity update filter error\n");
            info=0;
        }
        else {
            /* solution ok */
            //ins->stat=INSS_ZVU;
            info=1;
            clp(ins,opt,x);
            trace(3,"zero velocity update ok\n");
            printf("zero velocity update ok\n");
        }
    }
    free(x); free(H);
    free(R); free(v);
    return info;
}
/* ZUPT detection, based on   GREJNER-BRZEZINSKA et al. (2002) 
  If static: ++zvu_counter else zvu=0 if not ---------------------------------------*/
void detstc(ins_states_t *ins){
  // horizontal velocity and gyro components tolerance leves based on static INS data:
  double vn0, ve0, vn0std, ve0std, gyrx0, gyry0, gyrx0std, gyry0std;
   vn0 = -0.002743902895411; ve0 = -0.002219817510341;
   vn0std = 0.001958782893669; ve0std = 0.001549122618107;
   gyrx0 = -0.000152716; gyry0 = 0.000386451;
   gyrx0std = 0.000483409; gyry0std = 0.001271191;
   printf("static detection \n");

     if ( (fabs(ins->vn[0]) - vn0) <= 3*vn0std && (fabs(ins->vn[1]) - ve0) <= 3*ve0std ) {
       zvu_counter++;
       printf("static detection: yes: %d\n", zvu_counter);
      }else {zvu_counter=0;printf("static detection: no: %d\n", zvu_counter);}
}

/* Output imu raw data to file */
void outputrawimu(ins_states_t *insc){
  fprintf(out_raw_fimu, "%lf %lf %lf %lf %lf %lf %lf\n",insc->data.sec, \
     insc->data.fb0[0], insc->data.fb0[1],insc->data.fb0[2], insc->data.wibb0[0],insc->data.wibb0[1],
     insc->data.wibb0[2]);
}

/* Output ins/gnss solution record to respective files */
void outputinsgnsssol(ins_states_t *insc, insgnss_opt_t *opt, prcopt_t *gnssopt, 
 int n, obsd_t *obs){
  int i,j, sat;

  /* Output PVA solution     */ 
  if (insc->ptime>0.0) {
    fprintf(out_PVA,"%lf %.12lf %.12lf %lf %lf %lf %lf %lf %lf %lf %d\n",\
    insc->time, insc->rn[0]*R2D, insc->rn[1]*R2D, insc->rn[2],\
    insc->vn[0], insc->vn[1], insc->vn[2],
    insc->an[0]*R2D,insc->an[1]*R2D,insc->an[2]*R2D, opt->Nav_or_KF);  
  }

  /* Generate IMU bias output record */
  fprintf(out_IMU_bias_file, "%lf %lf %lf %lf %.10lf %.10lf %.10lf %d\n", insc->time,
        insc->data.ba[0], insc->data.ba[1], insc->data.ba[2],
        insc->data.bg[0], insc->data.bg[1], insc->data.bg[2], 
        opt->Nav_or_KF);
  
  /* Generate clock output record */
  if (xnRc(gnssopt)>1){
    fprintf(out_clock_file, "%lf %lf %lf %lf %d\n", insc->time, insc->dtr[0], 
    insc->dtr[1], insc->dtrr, opt->Nav_or_KF);
  }else{
    fprintf(out_clock_file, "%lf %lf %lf %d\n", insc->time, insc->dtr[0], 
    insc->dtrr, opt->Nav_or_KF);
  }  
  /* Generate Tropospheric delay  output record */
  fprintf(out_tropo_file, "%lf %lf %d\n", insc->time, insc->x[xiTr(gnssopt)],
  opt->Nav_or_KF);
 
  /* Generate Ambiguities output record */
  for (i=0;i<n&&i<MAXOBS;i++) {
    sat=obs[i].sat;
    fprintf(out_amb_file, "%lf %d %lf %d\n", insc->time, sat, insc->x[xiBs(gnssopt, sat)],
    opt->Nav_or_KF);  
  }   
 
  /* Generate KF uncertainty output record */
  fprintf(out_KF_SD_file, "%lf ", insc->time);
   for (i = 0; i < insc->nx; i++){
     for (j = 0; j < insc->nx; j++){
       (i==j?fprintf(out_KF_SD_file, "%lf ", SQRT(fabs(insc->P[i * insc->nx + j])) ):0);
      }
    }
  fprintf(out_KF_SD_file, "%d\n", opt->Nav_or_KF); 
} 

/* GNSS covariance to KF weights from gnss solution */
void pvclkCovfromgnss(rtk_t *rtk, ins_states_t *ins){
  int nx=ins->nx;

  /* Velocity full covariance */
  ins->P[IV*nx+IV]=rtk->sol.qrv[0];                                 /* xx or ee */
  ins->P[(IV+1)*nx+(IV+1)]=rtk->sol.qrv[1];                         /* yy or nn */
  ins->P[(IV+2)*nx+(IV+2)]=rtk->sol.qrv[2];                         /* zz or uu */
  ins->P[IV*nx+IV+1]=ins->P[(IV+1)*nx+IV]=rtk->sol.qrv[3];          /* xy or en */
  ins->P[(IV+1)*nx+(IV+2)]=ins->P[(IV+2)*nx+(IV+1)]=rtk->sol.qrv[4];/* yz or nu */
  ins->P[IV*nx+IV+2]=ins->P[(IV+2)*nx+IV]=rtk->sol.qrv[5];          /* zx or ue */

  /* Position full covariance */
  ins->P[IP*nx+IP]=rtk->sol.qr[0];                                 /* xx or ee */
  ins->P[(IP+1)*nx+(IP+1)]=rtk->sol.qr[1];                         /* yy or nn */
  ins->P[(IP+2)*nx+(IP+2)]=rtk->sol.qr[2];                         /* zz or uu */
  ins->P[IP*nx+IP+1]=ins->P[(IP+1)*nx+IP]=rtk->sol.qr[3];          /* xy or en */
  ins->P[(IP+1)*nx+(IP+2)]=ins->P[(IP+2)*nx+(IP+1)]=rtk->sol.qr[4];/* yz or nu */
  ins->P[IP*nx+IP+2]=ins->P[(IP+2)*nx+IP]=rtk->sol.qr[5];          /* zx or ue */

}

/* Core function -------------------------------------------------------------
* Description: Receive raw GNSS and INS data and determine a PVA solution
* args:
* rtk_t *rtk  IO
* const obsd_t *obs I
* int n I
* const nav_t *nav  I
* return:
* obs.: this function is ran inside rtkpos function of RTKlib
------------------------------------------------------------------------------*/
extern void core1(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav){ 
  int i, j, week, flag, core_count=0;
  double gnss_time, rr[3], ve[3];
  ins_states_t insc={0};
  prcopt_t *opt = &rtk->opt; 
 
 
  printf("\n *****************  CORE BEGINS *******************: %lf\n", time2gpst(rtk->sol.time,&week));

  /* Save gnss position and velocity */
  matcpy(rr,rtk->sol.rr,1,3);
  matcpy(ve,rtk->sol.rr+3,1,3);

  /* initialize ins/gnss parameter default uncertainty */ 
  //ig_paruncinit(&insgnssopt); 
  kf_par_unc_init(&insgnssopt);

  /* initialize ins state */
  insinit(&insc, &insgnssopt, opt, n);
  
  /* Initialize time from GNSS */
  gnss_time=time2gpst(rtk->sol.time,&week);

  /* Feed gnss solution and measurement buffers */
  gnssbuffer(&rtk->sol, obs, n);

  /* Ins and integration loop
  Processing window - integrates when GNSS and INS mea. are closer by 0.1s
  The do while takes care when INS or GNSS is too ahead from each other         */  
  do {
    if (insgnssopt.ins_ini){  
      printf("\n **** Ins Loop starts ****: %d \n", ins_w_counter);
    } 
    
    if(core_count>0){
      memset_ins_pva(&insc);    
    }

    /* Initialize ins states with previous state from ins buffer */
      if(ins_w_counter>insgnssopt.insw-1){
        /* After filling ins buffer */
        //printf("After filling ins buffer: %d\n", insgnssopt.insw-1);
        //print_ins_pva(insw+insgnssopt.insw-1);
       insc=insw[insgnssopt.insw-1];
       }else {
         /* First epoch */
        if(ins_w_counter>1) {
          //printf("Before filling ins buffer: %d\n", ins_w_counter-1);
          //print_ins_pva(insw+ins_w_counter-1);
          insc=insw[ins_w_counter-1];
        }
       }  

    /* input ins */ 
    if(!inputimu(&insc, week)) {printf(" ** End of imu file **\n"); break;}

    /* Output raw INS */
    if(insc.pdata.sec > 0.0 ){
      outputrawimu(&insc);
    }

    if (ins_w_counter >= 1){
      insc.dt = insc.time - insc.ptime;
    }

    printf("Ins time: %lf %lf %lf\n", insc.dt,insc.time,insc.ptime); 

    /* If INS time is ahead of GNSS exit INS loop */ 
    if (gnss_time-insc.time < -0.01) {
      // for tactical, -1.65 for consumer
      printf("\n ** Ins time ahead gps time **\n", insc.time,gnss_time); 
   
      /* Back INS file reader up */
      rewind(imu_tactical);

      insc.stat=-1;
      break;  
    }else{ 
      /* INS Navigation and/or INS/GNSS Integration */

      /* correction imu accl and gyro measurements data */
      ins_errmodel2(&insc.data.fb0,&insc.data.wibb0,
                  &insc.data.Ma,&insc.data.Mg,
                  &insc.data.ba,&insc.data.bg,&insc.data.Gg,
                  &insc.data.fb,&insc.data.wibb);

      /* initial ins states */
      // Here it's where the PVA initialization with the alignment is done 
      if (gnss_w_counter>2 && insgnssopt.ins_ini!=1){ 
        //150 means 1s of ins data, thus perform initialization only in the beginning 
        if(init_inspva(solw, &insc)){
          printf(" ** Ins initialization ok: %lf **\n", insc.time);
          insgnssopt.ins_ini=1;
        }else{
          printf(" ** Ins initialization error: %lf **\n", insc.time);
          insupdt(&insc);
          /* Add current ins measurement to buffer */  
          insbuffer(&insc);
          ins_w_counter++;
          insgnssopt.ins_ini=0;   
          continue;
        } 
      }      
   
      /* Check if ins and gps observations match for integration */
      flag=interp_ins2gpst(gnss_time, &insc);
      // HOW ABOUT INCLUDING THIS TO MODIFY FLAG TOO?: 
      /*if (fabs(GNSS_measurements->sec - INS_measurements->sec) < 0.0001 && GNSS_measurements[0].gdop[0] < 2.5 && no_GNSS_meas >= 4 &&
      norm(TC_KF_config->init_pos_unc_ned, 2) < 5.0)  */

      /* Integration */
      if (insgnssopt.mode){ // Add condition: if less than four satellites mode=0;
         /* Tightly-coupled ins/gnss */
         TC_INS_GNSS_core1(rtk, obs, n, nav, &insc, &insgnssopt, flag);
      }else{   
         /* Loosley-coupled ins/gnss */
         LC_INS_GNSS_core1(rtk, obs, n, nav, &insc, &insgnssopt, flag); 
      }  
      
      /* Non-holonomic constraints */
      if(insc.pdata.sec > 0.0 ){ 
        printf("nhc update\n");
        nhc(&insc,&insgnssopt); 
      }

      detstc(&insc);  
      printf("ZVU counter: %d\n", zvu_counter);
      /* Zero velocity update */
      if (zvu_counter>10) {
        printf("ZVU UPDATE: counter: %d", zvu_counter);  
        /* Zero-velocity constraints */
        zvu(&insc,&insgnssopt,1); 
      }  

      /* Output PVA, clock, imu bias solution     */
      if(insc.ptime>0.0){ 
        outputinsgnsssol(&insc, &insgnssopt, &rtk->opt, n, obs);  
      }

    } // end If INS time ahead of GNSS condition 

    /* Re-initializing ins with gnss solution when integration occurs */
     if (insgnssopt.Nav_or_KF==1){
      for(i=0;i<3;i++) insc.re[i]=rr[i];
      for(i=0;i<3;i++) insc.ve[i]=ve[i]; 
      /* Clock solution */
      insc.dtr[0]=rtk->sol.dtr[0]*CLIGHT; 
      update_ins_state_n(&insc);  
      /* Gnss solution covariance to ins */
      pvclkCovfromgnss(rtk, &insc); 
    }
 
    /* Zeroing closed-loop states */
     for (i=0;i<xnCl();i++) insc.x[i]=0.0; 

     /* Update ins states and measurements */  
     insupdt(&insc);

     /* Add current ins measurement to buffer */
     insbuffer(&insc);

     /* Global ins counter */ 
     ins_w_counter++;
     core_count++;
     
     printf("Ins counter updated: %d\n", ins_w_counter);
     printf("\n **** Ins Loop ends **** \n");

     /* Loop conditions */  
     if (gnss_time-insc.time<0.0001) { 
      printf("UPDATE INS AND GNSS: gpst: %lf, imut: %lf, dt diff: %lf\n", gnss_time, insc.time, gnss_time-insc.time);
      break;
     }else{
       if (gnss_time-insc.time<-0.0001) {
       printf("UPDATE INS AND GNSS: gpst: %lf, imut: %lf, dt diff: %lf\n", gnss_time, insc.time, gnss_time-insc.time);
       break;
       }else{
         printf("KEEP UPDATING INS: gpst: %lf, imut: %lf, dt diff: %lf\n", gnss_time, insc.time, gnss_time-insc.time);
        }
      }

  } while(gnss_time-insc.time > 0.0001 || \ 
   gnss_time-insc.time > -0.0001 );

  printf("\n ** Out of loop **\ngpst: %lf, imut: %lf, dt diff: %lf\n", gnss_time, insc.time, gnss_time-insc.time); 
 
   /* Update */  
 
   /* Global gnss counters */ 
   gnss_w_counter++; 

   /* Free memory */ 
   insfree(&insc);  

 printf("\n *****************  CORE ENDS ***********************\n");
}

int main(void){
/* Variables declaration =====================================================*/
prcopt_t prcopt=prcopt_default;
solopt_t solopt=solopt_default;
filopt_t filopt={""};
gtime_t ts={0},te={0};
double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3];
int i,j,k,n,ret;
char *infile[MAXFILE],*outfile=""; 
//char file[] = "../data/Lanes_XYZ_rtklib";//_corr";
//char file[] = "../out/LineInterp_rtklib3.txt";//_corr";
char file[] = "../out/LineInterp_rtklib5_kin_static_heights_using_corners_exp4.txt";//leverarm accounted
//fp_lane=fopen(file,"r");//lane_t lane;
//char imufilename[]="../data/LOG__010.SBF_SBF_ASCIIIn.txt"; //MEM-IMU file path exp1
//char imufilename[]="../data/LOG__033.SBF_SBF_ASCIIIn.txt"; //MEM-IMU file path exp2
//char imufilename[]="../data/LOG__038.SBF_SBF_ASCIIIn.txt"; //IMU file path exp3
//char imufilename[]="../data/LOG__040.SBF_SBF_ASCIIIn.txt"; //IMU file path exp4
//char imufilename[]="/home/emerson/Desktop/Connected_folders/data/imu_ascii_new.txt";
FILE *res;
//  char residualsfname[]="../out/exp1_PPP.stat"; //Residuals file
char residualsfname[]="../out/PPP_car_back.pos.stat"; //Residuals file
char tracefname[]="../out/trace.txt"; //trace file
int l=0,c;   
int argc; // Size of file or options? 
   
/* Global structures initialization */ 
insgnssopt.mode = 1; /* Tightly=1, or Loosley=0 coupled solution */
insgnssopt.Tact_or_Low = 1;       /* Type of inertial, tact=1, low=0 */
insgnssopt.scalePN = 0;             /* use extended Process noise model */ 
insgnssopt.gnssw = 3; 
insgnssopt.insw = 10; 
insgnssopt.exphi = 0;            /* use precise system propagate matrix for ekf */
/* ins sthocastic process noises: */
insgnssopt.baproopt=INS_RANDOM_WALK;
insgnssopt.bgproopt=INS_RANDOM_WALK;
insgnssopt.saproopt=INS_RANDOM_WALK;
insgnssopt.sgproopt=INS_RANDOM_WALK;
solw=(sol_t*)malloc(sizeof(sol_t)*insgnssopt.gnssw);  /* gnss solution structure window size allocation */
insw=(ins_states_t*)malloc(sizeof(ins_states_t)*insgnssopt.insw);   /* ins states window size allocation */
 
strcpy(filopt.trace,tracefname); 
   
/* Global TC_KF_INS_GNSS output files          */
out_PVA=fopen("../out/out_PVA.txt","w");
out_clock_file=fopen("../out/out_clock_file.txt","w");
out_IMU_bias_file=fopen("../out/out_IMU_bias.txt","w");
out_tropo_file=fopen("../out/out_tropo_bias.txt","w");
out_amb_file=fopen("../out/out_amb_bias.txt","w");
out_KF_state_error=fopen("../out/out_KF_state_error.txt","w");
out_KF_SD_file=fopen("../out/out_KF_SD.txt","w");
out_raw_fimu=fopen("../out/out_raw_imu.txt","w");
out_KF_residuals=fopen("../out/out_KF_residuals.txt","w");
//imu_tactical=fopen("../data/imu_ascii_new_1.txt", "r");
imu_tactical=fopen("../data/imu_ascii_new_timesync2.txt", "r");
fimu=fopen("../data/LOG__040.SBF_SBF_ASCIIIn.txt","r"); 
 
rewind(imu_tactical);                                                 



/* Declarations from rnx2rtkp source code program */
//clk93stream.rtcm3  CLK930600.rtcm3
// Base station: /home/emerson/Desktop/rtk_simul/PNW2_140448.16o

/* PPP-Kinematic  1st experimet
char *argv[] = {"./rnx2rtkp", "../data/SEPT2640.17O", "../data/igs19674.*", "../data/SEPT2640.17N", "-o", "../out/PPP_phs2_exp1.pos", "-k", "../config/opts3.conf"};

char *comlin = "./rnx2rtkp ../data/SEPT2640.17O ../data/igs19674.*  ../data/SEPT2640.17N -o ../out/PPP_phs2_exp1.pos -k ../config/opts3.conf";
  */ // Command line

/* PPP-Kinematic  2nd experimet
  char *argv[] = {"./rnx2rtkp", "../data/LOG__033.18O", "../data/igs20012.sp3", "../data/brdc1350.18n", "-o", "../out/PPP_test2.pos", "-k", "../config/opts3.conf"};

  char *comlin = "./rnx2rtkp ../data/LOG__033.18O ../data/igs20012.sp3  ../data/brdc1350.18n -o ../out/PPP_test2.pos -k ../config/opts3.conf";
     */// Command line

/* PPP-Kinematic  3rd experimet
  char *argv[] = {"./rnx2rtkp", "../data/LOG__038.18O", "../data/LOG__038_gps.nav", "-o", "../out/SPP_test3.pos", "-k", "../config/opts4.conf"};
  char *comlin = "./rnx2rtkp ../data/LOG__038.18O ../data/LOG__038_gps.nav -o ../out/SPP_test3.pos -k ../config/opts4.conf";
         // Command line */

/* PPP-Kinematic  4th experimet - GPS AND GLONASS */
 //char *argv[] = {"./rnx2rtkp", "../data/LOG__040.18o", "../data/BRDC00IGS_R_20181930000_01D_MN.rnx", "../data/grm20094.clk", "../data/grm20094.sp3", "-o", "../out/PPP_bmo.pos", "-k", "../config/opts3.conf", "-x", "5"};
 
 //char *comlin = "./rnx2rtkp ../data/LOG__040.18o ../data/BRDC00IGS_R_20181930000_01D_MN.rnx ../data/igs20094.* -o ../out/PPP_mod_exp4.pos -k ../config/opts3.conf";
// Command line 


/* PPP-Kinematic  Kinematic Positioning dataset
 char *argv[] = {"./rnx2rtkp", "../data/CAR_2890.18O", "../data/CAR_2890.18N", "../data/igs20232.*", "-o", "../out/PPP_car_back.pos", "-k", "../config/opts3.conf", "-x 5"};
 char *comlin = "./rnx2rtkp ../data/CAR_2890.18O ../data/CAR_2890.18N ../data/igs20232.* -o ../out/PPP_car_back.pos -k ../config/opts3.conf -x 5";
  */    // Command line 

/* PPP-Kinematic  Kinematic Positioning dataset  GPS+GLONASS */
//char *argv[] = {"./rnx2rtkp", "../data/CAR_2890.18O", "../data/BRDC00IGS_R_20182890000_01D_MN.rnx", "../data/grm20232.clk","../data/grm20232.sp3", "-o", "../out/PPP_car_back.pos", "-k", "../config/opts3.conf", "-x", "5"};

/* PPP-Kinematic  Kinematic Positioning dataset  March 21, 2019 GPS */
char *argv[] = {"./rnx2rtkp", "../data/APS_center.19O", "../data/APS_center.19N", "../data/igs20452.clk","../data/igs20452.sp3", "-o", "../out/PPP_march21.pos", "-k", "../config/opts3.conf", "-x", "5"};
 
/* PPP-AR Kinematic
char *argv[] = {"./rnx2rtkp", "../data/SEPT2640.17O", "../data/grg19674.*", "../data/SEPT2640.17N", "-o", "../out/exp1_PPP_amb_mod_constr.pos", "-k", "../config/opts3.conf"};

char *comlin = "./rnx2rtkp ../data/SEPT2640.17O ../data/grg19674.*  ../data/SEPT2640.17N -o ../out/exp1_PPP_amb_mod_constr.pos -k ../config/opts3.conf";
  */  // Command line

/* PPK - Post Processed Kinematic
   char *argv[] = {"./rnx2rtkp", "../data/2_APS_2644.17O", "../data/UNB3264.17o", "../data/1_APS_2640.17N", "../data/igs19674.sp3", "-o", "../out/PPK_exp1.pos", "-k", "../config/relative.conf"};

   char *comlin = "./rnx2rtkp ../data/2_APS_2644.17O ../data/UNB3264.17o ../data/1_APS_2640.17N ../data/igs19674.sp3 -o ../out/PPK_exp1.pos -k ../config/relative.conf";
// Command line*/

/* Static data
 char *argv[] = {"./rnx2rtkp", "/home/emerson/Desktop/rtk_simul/unbj_constrain_test/unbj0600.16o", "/home/emerson/Desktop/rtk_simul/unbj_constrain_test/brdc0600.16n", "/home/emerson/Desktop/rtk_simul/unbj_constrain_test/igs18861.sp3", "-o", "/home/emerson/Desktop/rtk_simul/unbj_constrain_test/sol_unbj1.pos", "-k", "/home/emerson/Desktop/rtk_simul/opts3.conf"};

   char *comlin = "./rnx2rtkp /home/emerson/Desktop/rtk_simul/unbj_constrain_test/unbj0600.16o /home/emerson/Desktop/rtk_simul/unbj_constrain_test/brdc0600.16n /home/emerson/Desktop/rtk_simul/unbj_constrain_test/igs18861.sp3 -o /home/emerson/Desktop/rtk_simul/unbj_constrain_test/sol_unbj1.pos -k /home/emerson/Desktop/rtk_simul/opts3.conf"; // Command line
*/

 //argc= sizeof(comlin) / sizeof(char);  
    argc=11; 

    prcopt.mode  =PMODE_KINEMA;
  //  prcopt.navsys=SYS_GPS|SYS_GLO;   
    prcopt.refpos=1;
    prcopt.glomodear=1;
    solopt.timef=0;
    sprintf(solopt.prog ,"%s ver.%s",PROGNAME,VER_RTKLIB);
    sprintf(filopt.trace,"%s.trace",PROGNAME);
/* Declarations from rnx2rtkp */

/* Reference lane size measurement */
/*for(;;){
  c =fgetc(fp_lane);
    if( c == EOF ){break;
    }else if ( c == '\n') {
      l++;
    }
}
   fseek(fp_lane, 0L, SEEK_END);
   lane.size = ftell(fp_lane);
   lane.sizeline = lane.size/l;
   rewind(fp_lane);
*/

/* load options from configuration file */
//printf("%d\n",argc);

    for (i=1;i<argc;i++) {
    printf("input: %s \n", argv[i]);
        if (!strcmp(argv[i],"-k")&&i+1<argc) {
            resetsysopts();
            if (!loadopts(argv[++i],sysopts)) return -1;
            getsysopts(&prcopt,&solopt,&filopt);
        }
    }
    for (i=1,n=0;i<argc;i++) {
        if      (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-ts")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",es,es+1,es+2);
            sscanf(argv[++i],"%lf:%lf:%lf",es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",ee,ee+1,ee+2);
            sscanf(argv[++i],"%lf:%lf:%lf",ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-ti")&&i+1<argc) tint=atof(argv[++i]);
        else if (!strcmp(argv[i],"-k")&&i+1<argc) {++i; continue;}
        else if (!strcmp(argv[i],"-p")&&i+1<argc) prcopt.mode=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-f")&&i+1<argc) prcopt.nf=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-m")&&i+1<argc) prcopt.elmin=atof(argv[++i])*D2R;
        else if (!strcmp(argv[i],"-v")&&i+1<argc) prcopt.thresar[0]=atof(argv[++i]);
        else if (!strcmp(argv[i],"-s")&&i+1<argc) strcpy(solopt.sep,argv[++i]);
        else if (!strcmp(argv[i],"-d")&&i+1<argc) solopt.timeu=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-b")) prcopt.soltype=1;
        else if (!strcmp(argv[i],"-c")) prcopt.soltype=2;
        else if (!strcmp(argv[i],"-i")) prcopt.modear=2;
        else if (!strcmp(argv[i],"-h")) prcopt.modear=3;
        else if (!strcmp(argv[i],"-t")) solopt.timef=1;
        else if (!strcmp(argv[i],"-u")) solopt.times=TIMES_UTC;
        else if (!strcmp(argv[i],"-e")) solopt.posf=SOLF_XYZ;
        else if (!strcmp(argv[i],"-a")) solopt.posf=SOLF_ENU;
        else if (!strcmp(argv[i],"-n")) solopt.posf=SOLF_NMEA;
        else if (!strcmp(argv[i],"-g")) solopt.degf=1;
        else if (!strcmp(argv[i],"-r")&&i+3<argc) {
            prcopt.refpos=0;
            for (j=0;j<3;j++) prcopt.rb[j]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i],"-l")&&i+3<argc) {
            prcopt.refpos=0;
            for (j=0;j<3;j++) pos[j]=atof(argv[++i]);
            for (j=0;j<2;j++) pos[j]*=D2R;
            pos2ecef(pos,prcopt.rb);
        }
        else if (!strcmp(argv[i],"-y")&&i+1<argc) solopt.sstat=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-x")&&i+1<argc) solopt.trace=atoi(argv[++i]);
        else if (*argv[i]=='-') printhelp();
        else if (n<MAXFILE) infile[n++]=argv[i];
    }                                                     //end of for
    if (n<=0) {
        showmsg("error : no input file");
        return -2;
    }

/* Start rnx2rtkp processing  ----------------------------------- --*/ 
/* Processing ------------------------------------------------------*/
 ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
 if (!ret) fprintf(stderr,"%40s\r","");
 
  for (i=0;i<insgnssopt.insw;i++) insfree(insw+i); 
  free(solw); free(insw); 

 /* ins navigation only */
 //imu_tactical_navigation(imu_tactical);

 /* Closing global files          */
  //fclose(fp_lane);
  //fclose (fimu);
  fclose(out_PVA);
  fclose(out_clock_file);
  fclose(out_IMU_bias_file); 
  fclose(out_tropo_file);
  fclose(out_amb_file);
  fclose(out_KF_SD_file);
  fclose(out_raw_fimu);
  fclose(imu_tactical);
  fclose(out_KF_state_error); 
  fclose(out_KF_residuals);  
  fclose(fimu);                  

/*
  char posfile[]="../out/out_PVA.txt";
  imuposplot(posfile);
  char gyrofile[]="../out/out_PVA.txt";   
  imueulerplot(gyrofile);                            */

/* Other function call after processing ------------------------------*/


/* INS/GNSS plots  */ 
char out_raw_fimu[]="../out/out_raw_imu.txt"; 
imuaccplot(out_raw_fimu);
imugyroplot(out_raw_fimu);
char gyrofile[]="../out/out_PVA.txt";  
imueulerplot(gyrofile);
char velfile[]="../out/out_PVA.txt"; 
imuvelplot(velfile);
char posfile[]="../out/out_PVA.txt";
imuposplot(posfile);
char imu_bias[]="../out/out_IMU_bias.txt";        
imuaccbiasplot(imu_bias); 
imugyrobiasplot(imu_bias);
char imu_KF_stds[]="../out/out_KF_SD.txt"; 
KF_att_stds_plot(imu_KF_stds);
KF_vel_stds_plot(imu_KF_stds);
KF_pos_stds_plot(imu_KF_stds);
char imu_KF_clk[]="../out/out_clock_file.txt"; 
KF_clock_plot(imu_KF_clk);
/*char KF_states[]="../out/out_KF_state_error.txt";
KF_state_errors_plot_att(KF_states);
KF_state_errors_plot_vel(KF_states);
KF_state_errors_plot_pos(KF_states);
KF_state_errors_plot_accb(KF_states);
KF_state_errors_plot_gyrb(KF_states);
KF_state_errors_plot_clk(KF_states);  */
char imu_KF_res[]="../out/out_KF_residuals.txt";
KF_residuals_plot(imu_KF_res);
char tropo_file[]="../out/out_tropo_bias.txt";
tropo_plot(tropo_file);
char amb_file[]="../out/out_amb_bias.txt";
amb_plot(amb_file);                                  


/* Plot IMU ra measurements
char imu_raw_meas[]="../out/out_raw_imu.txt";
imuaccplot(imu_raw_meas);
imugyroplot(imu_raw_meas);
*/

/* Adjusting IMU tactical time
FILE *new=fopen("../data/imu_ascii_new_1.txt", "w");
char str[100];
char timeini[]="2018 10 16 19 5 35.600";
gtime_t t_ini;
double gpst;
int week=2023;

str2time(timeini, 0, 23, &t_ini);
gpst=time2gpst(t_ini, &week);
printf("Time of week: %lf\n",gpst );


double t_0=gpst+622.0;
//t_0 = 244535.296;  // Tactical clock corect time experiment
double t_end, t=0.0, t_beggin=0.0;
double a[3], g[3];
i=0;
while ( fgets(str, 100, imu_tactical)!= NULL ){
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,%lf", &hh, &m,&s, &ht);
  sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &t, &a[0], &a[1], &a[2], &g[0], &g[1], &g[2]);
  if(i==0){
      t_end = t_0;
      t_beggin=t;
  }else{
    t_end = t_0 + (t-t_beggin);
  }
  //printf("%lf %lf %lf\n", posp[0]*R2D, posp[1]*R2D, posp[2]);
  //t=(double)(hh+(m/60.0)+(s/3600.0));
  fprintf(new, "%.3lf %lf %lf %lf %lf %lf %lf\n", t_end, a[0], a[1], a[2], g[0], g[1], g[2]);
  i++;
}

fclose(new);*/



/*Processing imu_tactical only */
//imuaccplot("/home/emerson/Desktop/Connected_folders/data/imu_ascii.txt");
//imugyroplot("/home/emerson/Desktop/Connected_folders/data/imu_ascii.txt");
//imu_tactical_navigation(imu_tactical);


/* IMU-MAP test
  while ( 1 ) {
    if(ferror(fimu)||feof(fimu)) break;
    insmap ();
}*/

/* IMU PLOTS
char accfile[]="../out/ins_acceleration.txt";
char gyrofile[]="../out/ins_gyro.txt";
char posfile[]="../out/insmap_test_pos.txt";
char velfile[]="../out/insmap_test_vel.txt";
   imuaccplot(accfile);
   imugyroplot(gyrofile);
//   imueulerplot();
   imuUpplot(posfile);
   imuvelplot(velfile);*/
/*
   char raw_imu[]="../out/out_raw_imu.txt";
   imuaccplot(raw_imu);
   imugyroplot(raw_imu);
   char gyrofile[]="../out/out_PVA.txt";
   imueulerplot(gyrofile);
   char velfile[]="../out/out_PVA.txt";
   imuvelplot(velfile);
   char posfile[]="../out/out_PVA.txt";
   imuposplot(posfile);
/**/

   // amb_and_var_plot ();

  /* Printing residuals for the grad seminar */
  //resprint4();
  //resinnovbounds();
  //resinnovqchisq();
  //resinnovqchisqtger();  

  //sppstdnsat();

  //return ret;
  /* Calling interpolation function
  FILE *fp=fopen("../data/Bmo_corners_static_from_rtkplot","r"); // static (Antenna centered) exp5
  //FILE *fp=fopen("../data/Bmo_corners_xyz_rtklib","r"); // Corner centered
//  FILE *fp=fopen("../data/Bmo_cartesian_corners_derived_with_height","r"); // Kinematic (antenna centered) exp5
  lineInterp(fp,0.1,4);
  fclose(fp);
*/


 /* For Map-matching plots
 double r[3], enu[3],r0[3], dr[3], pos0[3];
 double lat0=45.94312323*D2R, lon0=-66.64106835*D2R, h0=52.7309; // BMO corner 1
 double x,y,z;
 char str[150];

 //FILE *flaenu=fopen("../data/Lanes_enu_rtklib_corr.txt","w");
 FILE *flaenu=fopen("../data/Lanes_enu_rtklib3.txt","w");

 pos0[0]=lat0;pos0[1]=lon0;pos0[2]=h0; // datum
 pos2ecef(pos0, r0);

 printf("%lf,%lf,%lf\n",r0[0],r0[1],r0[2]);

  rewind(fp_lane);
  while ( fgets(str, 150, fp_lane)!= NULL ){
    //puts(str);
    sscanf(str, "%lf %lf %lf", &r[0], &r[1], &r[2]);

    for (j = 0; j < 3; j++) dr[j]=r[j]-r0[j];
    ecef2enu(pos0, dr, enu);

    //sscanf(str, "%*4d%*[/]%*2d%*[/]%*4d %*2d%*[:]%*2d%*[:]%*12lf   %13lf  %14lf\n", &y, &x);
    fprintf(flaenu, "%lf %lf %lf\n", enu[0], enu[1],enu[2]);
   }
   fclose(flaenu);
*/

//  resprint4();

 /* Height and 3D distance comparison plot  */
//  heightplot();
 //coordcompplot ();
 //coordcompcalcplot ();

//mapmatchplot();
//mapmatchplot1 ();

/*Print misclosure and residuals values from SPP vs MM  using Rescomp file */
//resprint2();
//resprint1();
//resprint3();

 /* IMU functions call	tests
 pos[0]= 1761566.127;
 pos[1]= -4078757.743;
 pos[2]= 4560891.627;
 insnav(fimu, pos);*/
 //fclose(fimu);

 /* Residuals plots functions
 res=fopen(residualsfname,"r");
 resprint(res);
 fclose(res);*/


/* Height comparison file preparation
FILE *f1,*f2,*f3,*f4,*f1out,*f2out,*f3out;
char str[150];
double t,s,ht,ep[6],posr[3],r[3],posp[3],enup[3],enu[3],dr[3],r0[3],pos0[3];
pos0[0]=45.94312323*D2R, pos0[1]=-66.64106835*D2R, pos0[2]=52.7309; // BMO corner 1
pos0[0]=(45+(56.0/60.0)+(35.241079/3600.0))*D2R;  //BMO corner with lever arm
pos0[1]=-(66+(38.0/60.0)+(27.85004/3600.0))*D2R;
pos0[2]=54.360;
int ww,hh,m;
gtime_t taux;

//f1=fopen("../out/Relative_ALLOFF_ANTENNASON_AMB_ON_wirthout_height.pos","r");
//f2=fopen("../out/PPPonly2Ddist.txt","r");
//f3=fopen("../out/PPP_only_exp1 (copy).pos","r");
//f3=fopen("../out/PPP_cand_enu.txt","r");
//f4=fopen("../out/PPK_enu.txt","r");
//f1out=fopen("../data/Magnet_kin_height.txt","w");
f2out=fopen("../data/Lanes_5_llh.csv","w");
//f2out=fopen("../out/PPK2Ddist.txt","w");
//f3out=fopen("../data/NrCan_kin_height1.txt","w");
//f3out=fopen("../data/GAPS_kin_height.txt","w");
//f3out=fopen("../data/exp1_PPP_back_height.pos","w");
//f3out=fopen("../data/exp1_PPP_back_height.pos","w");
pos2ecef(pos0, r0);
//printf("%lf,%lf,%lf\n",r0[0],r0[1],r0[2]);
double avge=0.0,avgn=0.0,avgup=0.0;
double avppponly[3],avgpppmod[3],avgrtk[3];
avppponly[0]=0.395488;avppponly[1]=0.436018;avppponly[2]=1.315044;
avgpppmod[0]=0.000212;avgpppmod[1]=0.002623;avgpppmod[2]=0.002917;
avgrtk[0]=0.179752;avgrtk[1]=0.277061;avgrtk[2]=0.271385;


/* Opening Reference Lane file and coverting cartesia coordinates to geodetic

while ( fgets(str, 200, fp_lane)!= NULL ){
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,%lf", &hh, &m,&s, &ht);
  sscanf(str, "%lf %lf %lf", &r[0], &r[1], &r[2]);
  ecef2pos(r, posp);
  //printf("%lf %lf %lf\n", posp[0]*R2D, posp[1]*R2D, posp[2]);
  //t=(double)(hh+(m/60.0)+(s/3600.0));
  fprintf(f2out, "%11.9lf,  %11.9lf, %11.9lf\n", posp[0]*R2D, posp[1]*R2D, posp[2]);
}
*/

/*
i=0;
f2=fopen("../out/PPPonly2Ddist.txt","r");
while ( fgets(str, 200, f2)!= NULL ){
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,%lf", &hh, &m,&s, &ht);
  sscanf(str, "%*lf %lf %lf %lf", &enup[0], &enup[1], &enup[2]);
  avge+=(enup[0]-avppponly[0])*(enup[0]-avppponly[0]);
  avgn+=(enup[1]-avppponly[1])*(enup[1]-avppponly[1]);
  avgup+=(enup[2]-avppponly[2])*(enup[2]-avppponly[2]);
  i++;
  //t=(double)(hh+(m/60.0)+(s/3600.0));
  //fprintf(f1out, "%lf %lf\n", t,ht);
}

printf("PPP-only std in E, N, UP: %lf, %lf, %lf\n",sqrt(avge/(i-1)),sqrt(avgn/(i-1)),sqrt(avgup/(i-1)));
fclose(f2);

i=0;
avge=0.0,avgn=0.0,avgup=0.0;
f2=fopen("../out/PPPmod2Ddist.txt","r");
while ( fgets(str, 200, f2)!= NULL ){
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,%lf", &hh, &m,&s, &ht);
  sscanf(str, "%*lf %lf %lf %lf", &enup[0], &enup[1], &enup[2]);
  avge+=(enup[0]-avgpppmod[0])*(enup[0]-avgpppmod[0]);
  avgn+=(enup[1]-avgpppmod[1])*(enup[1]-avgpppmod[1]);
  avgup+=(enup[2]-avgpppmod[2])*(enup[2]-avgpppmod[2]);
  i++;
  //t=(double)(hh+(m/60.0)+(s/3600.0));
  //fprintf(f1out, "%lf %lf\n", t,ht);
}

//printf("PPP-mod averages in E, N, UP: %lf, %lf, %lf\n",avge/i,avgn/i,avgup/i);
printf("PPP-mod std in E, N, UP: %lf, %lf, %lf\n",sqrt(avge/(i-1)),sqrt(avgn/(i-1)),sqrt(avgup/(i-1)));
fclose(f2);

i=0;
avge=0.0,avgn=0.0,avgup=0.0;
f2=fopen("../out/PPK2Ddist.txt","r");
while ( fgets(str, 200, f2)!= NULL ){
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,%lf", &hh, &m,&s, &ht);
  sscanf(str, "%*lf %lf %lf %lf", &enup[0], &enup[1], &enup[2]);
//  avge+=enup[0];
//  avgn+=enup[1];
//  avgup+=enup[2];
avge+=(enup[0]-avgrtk[0])*(enup[0]-avgrtk[0]);
avgn+=(enup[1]-avgrtk[1])*(enup[1]-avgrtk[1]);
avgup+=(enup[2]-avgrtk[2])*(enup[2]-avgrtk[2]);

  i++;
  //t=(double)(hh+(m/60.0)+(s/3600.0));
  //fprintf(f1out, "%lf %lf\n", t,ht);
}
//printf("RTK averages in E, N, UP: %lf, %lf, %lf\n",avge/i,avgn/i,avgup/i);
printf("RTK std in E, N, UP: %lf, %lf, %lf\n",sqrt(avge/(i-1)),sqrt(avgn/(i-1)),sqrt(avgup/(i-1)));
fclose(f2);


if(f1 == NULL) {
     perror("Error opening file");
     return(-1);
  }


while ( fgets(str1, 150, f1)!= NULL ){
            //2017/09/21 16:23:36.400   45.943136265  -66.641063914    56.5241
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf   %lf  %lf    %lf", &hh, &m,&s, &posr[0], &posr[1], &posr[2]);

           //2017/09/21 16:23:36.400,  45.943124507, -66.641093833,   52.0725,  6,  9, 23.6700, 14.2659, 34.0400, 11.4957,  3.8339, -5.5948,  0.00,   0.0
           //2018/07/12 13:40:06.000,  45.943122547, -66.641069572,   54.3898,  1,  9,  0.0065,  0.0063,  0.0137,  0.0022,  0.0025, -0.0051,  0.00,   5.4
  sscanf(str1, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,   %lf,  %lf,    %lf,", &hh, &m,&s, &posr[0], &posr[1], &posr[2]);
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf%*c%*lf%*c%*lf%*c%lf", &hh, &m,&s, &ht);

  t=(double)(hh+(m/60.0)+(s/3600.0));

  posr[0]=posr[0]*D2R; posr[1]=posr[1]*D2R;
  pos2ecef(posr, r);
 // printf("Read position: %lf %lf %lf\n",r[0],r[1], r[2]);
  // Fetching data to the first buffer structure (lane.buffer)
   inibuff();
  // Closest point computation
  clstpt(&r);
  // Closest point positions to local
	for (j=0;j<3;j++) posp[j]=lane.buffer[lane.currsearch*3+j];
  //printf("clst: %lf %lf %lf\n", posp[0], posp[1], posp[2]);
  for (j = 0; j < 3; j++) dr[j]=posp[j]-r0[j];
  ecef2enu(pos0, dr, enup);
  ecef2pos(posr, posp);
  //printf("%lf\n", posp[2]);

  // When xyz directly
  //sscanf(str, "%lf %lf %lf", &r[0], &r[1], &r[2]);

  for (j = 0; j < 3; j++) dr[j]=r[j]-r0[j];
  ecef2enu(pos0, dr, enu);


  //printf("rtk: %lf %lf %lf,clst: %lf %lf %lf\n", enu[0],enu[1],enu[2],enup[0],enup[1],enup[2]);
  // en and up differences
  //fprintf(f2out, "%lf %lf %lf\n", t, sqrt( pow(enup[0]-enu[0],2) + pow(enup[1]-enu[1],2)),enup[2]-enu[2]);

  // e, n, and up differences
fprintf(f2out, "%lf %lf %lf %lf\n", t, fabs(enup[0]-enu[0]), fabs(enup[1]-enu[1]), fabs(enup[2]-enu[2]));

  // e, n, and up
  //fprintf(f2out, "%lf %lf %lf\n", enu[0], enu[1], enu[2]);

  //fprintf(f2out, "%lf %lf\n",t,posr[2]);

}

/*
while ( fgets(str, 200, f3)!= NULL ){
  // NRCan csv file readings
  //sscanf(str, "%*lf%*c%*lf%*c%lf%*c%lf", &ht,&t);
  // GAPS pos file reading
  //sscanf(str, " %*4d   %*1d  %*2d  %2d  %2d  %lf       %*d       %*d    %*lf      %*d       %*d    %*lf        %lf", &hh, &m,&s, &ht);
  // Magnet tools
  //2017/09/21 16:23:36.400,53.559
  //sscanf(str, " %*4d/%*d/%*d  %2d:%2d:%lf,%lf", &hh, &m,&s, &ht);
  // RTKLIB - 2017/09/21 16:23:36.400   45.943136265  -66.641063914    56.5241
sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf   %*lf  %*lf    %lf", &hh, &m,&s, &ht);
//sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf%*c%*lf%*c%*lf%*c%lf", &hh, &m,&s, &ht);
t=(double)(hh+(m/60.0)+(s/3600.0));
   //printf("%d,%d,%lf,%lf\n",hh,m,s,ht);
  //printf("%lf,%lf\n",t, ht);
//  fprintf(f3out, "%lf %lf\n", t,ht);
}
*/
//fclose(f1);fclose(f3);fclose(f4);
//fclose(f1out);
//fclose(f2out);
//fclose(f3out);    */

 //heightplot();
 //resinnovbounds ();
 //mapmatchplot1 ();
 //mapmatchplot ();
 //nsat();

 //imuposplot();
 //imuaccplot();
 //imuvelplot();
 //imugyroplot();


 printf("\n\n SUCCESSFULLY EXECUTED!  \n\n");
 return;
}
