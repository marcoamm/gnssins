#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <dirent.h>

#include "rtklib.h"
#include "satinsmap.h"

/* global variables ----------------------------------------------------------*/

FILE *fp_lane;       /* Lane coordinate file pointer */
FILE *fimu;          /* Imu datafile pointer  */
FILE *imu_tactical; /* Imu datafile pointer */
lane_t lane;
imuraw_t imu_obs_global={{0}};
pva_t pva_global={{0}};
pva_t pvagnss={{0}};

char *outpath1[] = {"../out/"};
FILE *out_PVA;
FILE *out_clock_file;
FILE *out_IMU_bias_file;
FILE *out_KF_SD_file;
FILE *out_raw_fimu;


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
  double xyz_ini_pos[3], xyz_ini_cov[3], pos[3], head_angle=0.0;
  double xyz_prev_clst_pos[3];
  double ned_ini_vel[3]={0};
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
void IMU_meas_interpolation(float t_imu_prev, float t_imu_curr, float t_gps, \
  um7pack_t *imu_curr, imuraw_t *imu_prev){
    int i;
    double imu_t_gps_a[3], imu_t_gps_g[3];

    for (i=0;i<3;i++) imu_t_gps_a[i]=imu_prev->fb[i]+\
    ((imu_curr->a[i]-imu_prev->fb[i])/(t_imu_curr-t_imu_prev))*(t_gps-t_imu_prev);

    for (i=0;i<3;i++) imu_t_gps_g[i]=imu_prev->wibb[i]+\
    ((imu_curr->g[i]-imu_prev->wibb[i])/(t_imu_curr-t_imu_prev))*(t_gps-t_imu_prev);

    /* Updating imu_curr */
    for (i=0;i<3;i++) imu_curr->a[i]= imu_t_gps_a[i];
    for (i=0;i<3;i++) imu_curr->g[i]= imu_t_gps_g[i];
    imu_curr->sec=t_gps;
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
  double pos[3], head_angle=0.0, vxyz[3];
  double xyz_prev_clst_pos[3], llh_pos[3], aux, gan[3], wiee[3];
  double gnss_xyz_ini_cov[6];
  double Qenu[9]={0};
  double ned_ini_vel[3]={0};
  float ini_pos_time;
  char datastring[150];
  imuraw_t imu_obs_prev={0};
  um7pack_t imu_curr_meas={0};
  pva_t PVA_prev_sol={{0}};
  int j;
  int TC_or_LC=1; /* TC:1 and LC=0 */
  int Tact_or_Low_IMU=1; /* Tactical=1 and Low grade=0*/
  //FILE *ppp_llh;
  char str[100];
  double G = 9.80665;

  printf("\n *****************  CORE BEGINS ***********************\n");

  /* Update with previous solution */
  imu_obs_prev=imu_obs_global;
  PVA_prev_sol=pva_global;

  /* Initialize time from GNSS */
  ini_pos_time=time2gpst(rtk->sol.time,NULL);

  /* Initialize position from GNSS SPP */
  for (j=0;j<3;j++) xyz_ini_pos[j]=rtk->sol.rr[j];

  /* Velocity from GNSS in case there is no doppler observable */
  if (rtk->sol.rr[3] <= 0.0 || rtk->sol.rr[4] <= 0.0) {
    /* NO VELOCITY */
    for (j=0;j<3;j++) vxyz[j]=(xyz_ini_pos[j]-pvagnss.r[j])/(ini_pos_time-pvagnss.sec);
    ecef2enu(pvagnss.r,vxyz,ned_ini_vel); /* Here is ENU! */
    printf("Estimated velocity: %lf %lf %lf\n",ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2]);
  }else{
    ecef2pos(xyz_ini_pos,llh_pos);
    ecef2enu(llh_pos,rtk->sol.rr+3,ned_ini_vel); /* Here is ENU! */
  } 

  /* ENU to NED velocity from gnss solution */
     aux=ned_ini_vel[1];
     ned_ini_vel[1]=ned_ini_vel[0];
     ned_ini_vel[0]=aux;
     ned_ini_vel[2]=-ned_ini_vel[2];

  /* Position covariance */
  for (j = 0; j < 6; j++) gnss_xyz_ini_cov[j]=rtk->sol.qr[j];
  for (j = 0; j < 6; j++) gnss_vel_cov[j]=rtk->sol.qrv[j];
  covenu(llh_pos, gnss_vel_cov, Qenu);

  /* Velocity covariance in NED */
  gnss_vel_cov[0] = Qenu[1];
  gnss_vel_cov[1] = Qenu[0];
  gnss_vel_cov[2] = -Qenu[2];

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

  /* Initialize velocity with GNSS for coarseAlignment */
  for (j=0;j<3;j++) imu_curr_meas.v[j]=ned_ini_vel[j];

 /* Processing window - integrates when GNSS and INS mea. are closer by 0.1s
 The do while takes care when INS or GNSS is too ahead from each other         */
 do {

 /*  Low-cost imu
   memset(&imu_curr_meas, 0, sizeof(um7pack_t));
   /* Read imu current measurements
   while(imu_curr_meas.status!=1){
     if(!fgets(datastring, 150, fimu)) return;
     parseimudata(datastring,&imu_curr_meas);
   }
   // Fixing imu time
   imu_curr_meas.sec=imu_curr_meas.sec+ \
   (imu_curr_meas.internal_time-(float)floor(imu_curr_meas.internal_time));
*/

  /* Tactical IMU reading */
   fgets(str, 100, imu_tactical);
   sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &imu_curr_meas.sec, &imu_curr_meas.a[2],\
   &imu_curr_meas.a[1],&imu_curr_meas.a[0], &imu_curr_meas.g[2],&imu_curr_meas.g[1],\
   &imu_curr_meas.g[0]);
   imu_curr_meas.status=1;

   /* raw acc. to m/s*s and rate velocity from degrees to radians */
   for (j=0;j<3;j++) imu_curr_meas.a[j]=imu_curr_meas.a[j]*G;
   for (j=0;j<3;j++) imu_curr_meas.g[j]=imu_curr_meas.g[j]*D2R;


   /* Output raw INS */
   fprintf(out_raw_fimu, "%f %lf %lf %lf %lf %lf %lf\n",imu_curr_meas.sec, \
   imu_curr_meas.a[0],\
   imu_curr_meas.a[1],imu_curr_meas.a[2], imu_curr_meas.g[0],imu_curr_meas.g[1],
   imu_curr_meas.g[2]);
   
   printf("GNSS.INS.INSPREV.TIMES.DIFF: %.3f, %.3lf, %lf\n",roundf(ini_pos_time),\
   imu_curr_meas.sec, ini_pos_time-imu_curr_meas.sec );

  /* Interpolate INS measurement to match GNSS time
  if (imu_obs_prev.sec>0.0) {
    if (ini_pos_time < imu_curr_meas.sec && ini_pos_time > imu_obs_prev.sec) {
      printf("INTERPOLATE INS MEAS TO GNSS TIME \n");
      /* imu_curr_meas is modified
      IMU_meas_interpolation(imu_obs_prev.sec, imu_curr_meas.sec, ini_pos_time, \
        &imu_curr_meas, &imu_obs_prev);
      /*Return imu file pointer to process the imu_curr in the next epoch
       //imufileback();
    }
  }
 */
  /* Stationary-condition detection --------------------*/
   /* Horizontal velocity threshold for land-vehicles:
   low-cost imu < 0.5 m/s per axis for aviation-grade IMU < 0.0075 m/s per axis*/
   /* Stationary detection walking --------------------
    Use a window of 0.5 second and 1.5 m/s-2 threshold is suitable
    gn(llh_pos[0],llh_pos[2])
    if (fabs( (norm(imu_curr_meas.a,3)-gn(llh_pos[0],llh_pos[2])) )<1.626) {
    //  printf("Is Stationary \t");
     }else{ //printf("Is Moving \t");
    }
  */

  /* Horizontal condition based on velocities derived from INS measurements /
    0.5590 corresponds to the horizontal resultant threshold
  if (norm(PVA_prev_sol.v,2) > 0.55901) {
    printf("IS MOVING\n");
  }else{printf("IS STATIC\n");}  */

  /* Stationary-condition */
  if (norm(ned_ini_vel,3) > 0.5) {
  //  printf("IS MOVING\n");
  }else{//printf("IS STATIC\n");
    //for (j=0;j<3;j++) PVA_prev_sol.v[j]=0.0;
  }

  /* Attitude Initialization (Groves, 2013)  */

  /* Earth rotation vector in e-frame	*/
  wiee[0]=0;wiee[1]=0;wiee[2]=OMGE;

  /* local apparent gravity vector */
  appgrav(llh_pos, gan, wiee);

  /* Coarse alignment or use of previous solution? */
  /* Levelling and gyrocompassing, update imu->aea[] vector */
  coarseAlign(&imu_curr_meas, gan);

  /* Levelled period - 7 initial minutes
  if ( fabs((ini_pos_time-394469.000)/60.0) < 7) {
    PVA_prev_sol.A[0]=imu_curr_meas.aea[0]=0.0;
    PVA_prev_sol.A[1]=imu_curr_meas.aea[1]=0.0;
  }
  */
  //PVA_prev_sol.A[2]=imu_curr_meas.aea[2]=head_angle;

 /* If INS time is ahead of GNSS */
  if (ini_pos_time-imu_curr_meas.sec < -0.01) {
    /* Back INS file reader */
    printf("GNSS time behind imu time!! INS stays in the same epoch\n");
    rewind(imu_tactical);
    //imufileback();
    break;
  }else{

/* Tightly or Loosley-coupled integration */
/* If enough GNSS measurements uses Tightly otherwise Loosley-coupled*/
  if (TC_or_LC) {

    /* Tightly coupled INS/GNSS */
    if(imu_obs_prev.sec > 0.0){
      TC_INS_GNSS_core(rtk, obs, n, nav, &imu_curr_meas, (double)imu_curr_meas.sec\
      - imu_obs_prev.sec, &PVA_prev_sol, Tact_or_Low_IMU);
    }else{
      for (j=0;j<3;j++) PVA_prev_sol.r[j]=xyz_ini_pos[j];
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
      LC_INS_GNSS_core(xyz_ini_pos, gnss_xyz_ini_cov, ned_ini_vel, gnss_vel_cov, ini_pos_time, \
      &imu_curr_meas, (double)imu_curr_meas.sec - imu_obs_prev.sec,\
      &PVA_prev_sol, &imu_obs_prev, Tact_or_Low_IMU);
    }else{
      for (j=0;j<3;j++) PVA_prev_sol.r[j]=xyz_ini_pos[j];
      for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
      for (j=0;j<3;j++) PVA_prev_sol.A[j]=imu_curr_meas.aea[j];
      PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
      PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;
      PVA_prev_sol.sec=imu_curr_meas.sec;
      PVA_prev_sol.time=imu_curr_meas.time;
    }
  }

  /* Update measurements */
  for (j=0;j<3;j++){
  imu_obs_prev.fb[j]=imu_curr_meas.a[j];
  imu_obs_prev.wibb[j]=imu_curr_meas.g[j];
  }
  imu_obs_prev.sec=imu_curr_meas.sec;
  imu_obs_prev.time=imu_curr_meas.time;

  ecef2pos(PVA_prev_sol.r,llh_pos);

  /* Output PVA solution */
  fprintf(out_PVA,"%f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",\
  PVA_prev_sol.sec, llh_pos[0]*R2D, llh_pos[1]*R2D, llh_pos[2],\
  PVA_prev_sol.v[0], PVA_prev_sol.v[1], PVA_prev_sol.v[2],
  PVA_prev_sol.A[0]*R2D,PVA_prev_sol.A[1]*R2D,PVA_prev_sol.A[2]*R2D );
}

   /* Loop conditions */
   if (ini_pos_time-imu_curr_meas.sec<0.0001) {
     printf("UPDATE INS AND GNSS \n");
     break;
   }else{
     if (ini_pos_time-imu_curr_meas.sec<-0.0001) {
       printf("UPDATE INS AND GNSS \n");
       break;
     }else{
       printf("KEEP UPDATING INS \n");
     }
   }

 } while(ini_pos_time-imu_curr_meas.sec > 0.0001 || \
   ini_pos_time-imu_curr_meas.sec > -0.0001 );


 printf("OUT OF LOOP\n \
 IMU_TIME_AND_GNSS_TIMES: %f %f\n", imu_curr_meas.sec, ini_pos_time);


  // printf("Vel NED GNSS: %lf, %lf, %lf\n",ned_ini_vel[0],ned_ini_vel[1],ned_ini_vel[2]);
  // printf("Vel NED INS: %lf, %lf, %lf\n",PVA_prev_sol.v[0],PVA_prev_sol.v[1],PVA_prev_sol.v[2]);

   /* Making GNSS position and velocity be the INS initialization on the next
    processing, instead of the Navigation-derived one */
   for (j=0;j<3;j++) PVA_prev_sol.r[j]=xyz_ini_pos[j];
   for (j=0;j<3;j++) PVA_prev_sol.v[j]=ned_ini_vel[j];
   PVA_prev_sol.clock_offset_drift[0]=rtk->sol.dtr[0]*CLIGHT;
   PVA_prev_sol.clock_offset_drift[1]=rtk->sol.dtrr;

   /* Update *//* observation data doppler frequency (Hz) */
   PVA_prev_sol.t_s = time2gpst(rtk->sol.time,NULL); /* time of week in (GPST) */
   imu_obs_global=imu_obs_prev;
   pva_global=PVA_prev_sol;

   for (j=0;j<3;j++) pvagnss.r[j]=xyz_ini_pos[j];
   pvagnss.sec=PVA_prev_sol.t_s;

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
char imufilename[]="../data/LOG__040.SBF_SBF_ASCIIIn.txt"; //IMU file path exp4
//char imufilename[]="/home/emerson/Desktop/Connected_folders/data/imu_ascii_new.txt";
//fimu=fopen(imufilename,"r");
FILE *res;
//  char residualsfname[]="../out/exp1_PPP.stat"; //Residuals file
char residualsfname[]="../out/PPP_car_back.pos.stat"; //Residuals file
int l=0,c;
int argc; // Size of file or options?

/* Global TC_KF_INS_GNSS output files*/
out_PVA=fopen("../out/out_PVA.txt","w");
out_clock_file=fopen("../out/out_clock_file.txt","w");
out_IMU_bias_file=fopen("../out/out_IMU_bias.txt","w");
out_KF_SD_file=fopen("../out/out_KF_SD.txt","w");
out_raw_fimu=fopen("../out/out_raw_imu.txt","w");
imu_tactical=fopen("../data/imu_ascii_new.txt", "r");

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

/* PPP-Kinematic  4th experimet
 char *argv[] = {"./rnx2rtkp", "../data/LOG__040.18o", "../data/LOG__040_gps.18n", "../data/igs20094.*", "-o", "../out/PPP_mod_exp4.pos", "-k", "../config/opts3.conf"};
 char *comlin = "./rnx2rtkp ../data/LOG__040.18o ../data/LOG__040_gps.18n ../data/igs20094.* -o ../out/PPP_mod_exp4.pos -k ../config/opts3.conf";
// Command line
*/

/* PPP-Kinematic  Kinematic Positioning dataset */
 char *argv[] = {"./rnx2rtkp", "../data/CAR_2890.18O", "../data/CAR_2890.18N", "../data/igs20232.*", "-o", "../out/PPP_car_back.pos", "-k", "../config/opts3.conf"};
 char *comlin = "./rnx2rtkp ../data/CAR_2890.18O ../data/CAR_2890.18N ../data/igs20232.* -o ../out/PPP_car_back.pos -k ../config/opts3.conf";
      // Command line

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
 argc= sizeof(comlin) / sizeof(char);

    prcopt.mode  =PMODE_KINEMA;
    prcopt.navsys=SYS_GPS|SYS_GLO;
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

 /* Closing global files */
  //fclose(fp_lane);
  //fclose (fimu);
  fclose(imu_tactical);
  fclose(out_PVA);
  fclose(out_clock_file);
  fclose(out_IMU_bias_file);
  fclose(out_KF_SD_file);
  fclose(out_raw_fimu);


/* Other function call after processing ------------------------------*/

/*
FILE *new=fopen("/home/emerson/Desktop/Connected_folders/data/imu_ascii_new.txt", "w");
char str[100];
float t_0=242844.99572, t_end, t;
double a[3], g[3];
i=0;
while ( fgets(str, 100, imu_tactical)!= NULL ){
  //sscanf(str, "%*4d%*c%*2d%*c%*2d %2d%*c%2d%*c%lf,%lf", &hh, &m,&s, &ht);
  sscanf(str, "%f %lf %lf %lf %lf %lf %lf", &t, &a[0], &a[1], &a[2], &g[0], &g[1], &g[2]);
  if(i==0){
      t_end = t_0;
  }else{
    t_end = t_0 + t;
  }
  //printf("%lf %lf %lf\n", posp[0]*R2D, posp[1]*R2D, posp[2]);
  //t=(double)(hh+(m/60.0)+(s/3600.0));
  fprintf(new, "%6.5f %lf %lf %lf %lf %lf %lf\n", t_end, a[0], a[1], a[2], g[0], g[1], g[2]);
  i++;
}
*/

/*Processing imu_tactical only */
//imuaccplot("/home/emerson/Desktop/Connected_folders/data/imu_ascii.txt");
//imugyroplot("/home/emerson/Desktop/Connected_folders/data/imu_ascii.txt");
//imu_tactical_navigation(imu_tactical);
/**/
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
