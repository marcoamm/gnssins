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
int dz_counter = 0;
res_t resid={0};
static_info_t staticInfo={0};
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

/* Local apparent gravity ----------------------------------------------------
* description: computes the gravity vector in the navigation frame (n-frame)
* args	:	double *pos[3]	I	{llh} vector [m] in e-frame
*	 	double *gan[3]	O	apparent gravity vector in n-frame [m/s*s]
* Reference: Chatfiled (1997, pag. 80)
*----------------------------------------------------------------------------*/
void appgrav (double* pos, double* gan, double* wiee){
 double ge[3], r[3];
 double auxvec[3], aux1vec[3], gae[3], Ren[9];
 int i;

 /* llh to xyz	*/
 pos2ecef(pos, r);

 /* apparent gravity vector in e-frame
 from normal gravity vector at the location and altitude	*/
 ge[0]=0;ge[1]=0;ge[2]=gn(pos[0], pos[2]);

 cross3(wiee, r, auxvec);
 cross3(wiee, auxvec, aux1vec);

 for (i=0; i<3;i++){gae[i]=ge[i]-aux1vec[i];}

 /* Rotation matrix from e to n-frame	*/
 Ce2n(pos[0],pos[1],Ren);

 /* Local apparent gravity vector in n-frame	*/
 matmul_row("NN", 3, 1, 3, 1.0, Ren, gae, 0.0, gan);

 for (i=0; i<3;i++) gan[i]=-gan[i];

 /* Forcing gan back to ge, as in Shin (2001) page 24 */
 for (i=0; i<3;i++) gan[i]=ge[i];


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
/* input imu measurement data --------------------------------------------------*/
/* Return imu file pointer one complete set of observation  */
extern void imufileback(){
  int j;
  char c;

  for (j = 0; j <= 3 ; j++) {  // Returns 2 sets of obs
    while (c != '\n'){
      fseek(imu_tactical, -1L, SEEK_CUR);
      c = fgetc(imu_tactical);
      fseek(imu_tactical, -1L, SEEK_CUR);
    }
    c=0;
  }
  /* to match the begining of the line after the \n */
  fseek(imu_tactical, 1, SEEK_CUR);
}
/* Parse and organize the MEMs-IMU received packets ----------------------------
* description: IMU
* args   : char *strline       I   char with one line containing one sensor IMU readings
* One epoch message type:
489629.200,1946,COM3,36,$PCHRS,0,101.633,3.32,0.33,-0.49,*7E
489629.200,1946,COM3,42,$PCHRS,1,101.633,0.3717,0.1030,-0.7725,*77
489629.200,1946,COM3,42,$PCHRS,2,101.633,0.4922,0.2708,-1.6283,*7D, where
secofweek sync, GPS week, port, ?, header filed, count, sensor time, x,y,z,checksum
    O  imu.status -1: meassage not complete, 0: error reading, 1: ok
*-----------------------------------------------------------------------------*/
void parseimudata(char* strline, um7pack_t* imu)
{
   char *token;
   int i=1, sscanstat;
   float sensor_x,sensor_y,sensor_z;
   double imucurrtime;
   double G = 9.80665;

  sscanstat=sscanf(strline, "%lf,%d,%*4s,%*2d,%*6s,%d,%f,%*f,%*f,%*f,%*s",
  &imu->sec, &imu->gpsw, &imu->count, &imu->internal_time);

  /* Checking if received message is complete and synchronized	*/
  if (sscanstat < 4) {
    imu->status=-1;
    return;
  }

  sscanstat=sscanf(strline, "%lf,%d,%*4s,%*2d,%*6s,%d,%*f,%f,%f,%f,%3c",
  &imu->sec, &imu->gpsw, &imu->count, &sensor_x, &sensor_y, &sensor_z, &imu->checksum);
  /*printf("%6.1f, %d, %d, %f, %f, %f, %c \n",
  imu->sec, imu->gpsw, imu->count, sensor_x, sensor_y, sensor_z, imu->checksum);*/
  imu->time=gpst2time(imu->gpsw,imu->sec);

 /* Check if message was read successfully	*/
 if (sscanstat < 7) {
   imu->status=0;
   return;
 }

 if (imu->count == 0){
  // printf("Gyro data: %d \n", imu->count);
   imu->g[0]=sensor_x*D2R; imu->g[1]=sensor_y*D2R;imu->g[2]=sensor_z*D2R; /* rad/s */
   imu->tgyr=imu->internal_time;
 }else if (imu->count == 1){
       	//  printf("Accelerometer data: %d \n", imu->count);
	  imu->a[0]=sensor_x*G;imu->a[1]=sensor_y*G;imu->a[2]=sensor_z*G; /* in gravities to m/s2 */
    imu->tacc=imu->internal_time;
       }else {
	 // printf("Magnetometer data: %d \n", imu->count);
	  imu->m[0]=sensor_x;imu->m[1]=sensor_y;imu->m[2]=sensor_z; /* uT unitless */
    imu->tmag=imu->internal_time;
       }

   /* Checking if data is complete for processing	*/
   /* Using time of each data measurement */
   if (imu->tgyr==0.0 || imu->tacc==0.0 ) {
     /* Not complete yet */
     imu->status=0;
     return;
   }

   if (imu->tgyr==imu->tacc) {
     /* Data complete and synchronized */
     imu->status=1;
   }else{
     /* Data not synchronized */
     imu->status=0;
   }
   //if ( imu->count==2 ) imu->status = 1;

}
/* Imu input data --------------  */
static int inputimu(prcopt_t *opt, ins_states_t *ins, int week){
  char check, str[150];
  um7pack_t imu_curr_meas={0};
  int i, j;

  if(insgnssopt.Tact_or_Low){
    /* Tactical KVH input */
    // For Oct/March experiments: 2,1,0
    check=fgets(str, 150, imu_tactical);
    sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &ins->time, &ins->data.fb0[2],\
    &ins->data.fb0[1],&ins->data.fb0[0], &ins->data.wibb0[2],&ins->data.wibb0[1],\
    &ins->data.wibb0[0]);
    
    ins->time=ins->time+16; /* Accounting for leap seconds */
    ins->data.sec=ins->time;

    ins->data.time=gpst2time(week, ins->time);

    /* raw acc. to m/s*s and rate ve locity from degrees to radians  */
    for (j= 0;j<3;j++) ins->data.fb0[j]=ins->data.fb0[j]*Gcte;
    for (j=0;j<3;j++) ins->data.wibb0[j]=ins->data.wibb0[j]*D2R;

    printf("Acfilt: %lf %lf %lf - %lf %lf %lf\n", ins->pdata.fb0[0],ins->pdata.fb0[1],\
    ins->pdata.fb0[2], ins->data.fb0[0],ins->data.fb0[1],ins->data.fb0[2]);

    //Filter acceleration data low pass t = t-1*a + t*b
    //ins->data.fb0[0] = ins->pdata.fb0[0] * opt->acfilt[0] + ins->data.fb0[0] * opt->acfilt[1];
    //ins->data.fb0[1] = ins->pdata.fb0[1] * opt->acfilt[0] + ins->data.fb0[1] * opt->acfilt[1];
    //ins->data.fb0[2] = ins->pdata.fb0[2] * opt->acfilt[0] + ins->data.fb0[2] * opt->acfilt[1];

    printf("IMU.raw.read: %lf %lf %lf %lf %lf %lf %lf - check: %d\n", ins->time, ins->data.fb0[2],\
    ins->data.fb0[1],ins->data.fb0[0], ins->data.wibb0[2],ins->data.wibb0[1],\
    ins->data.wibb0[0], check);

    /* Turn-on bias on x and y axis */
    //  ins->data.fb0[0]=ins->data.fb0[0]-0.869565;//-0.850372;
    //  ins->data.fb0[1]=ins->data.fb0[1]+0.193414;//+0.200865;

    if (check==NULL && check ==0) {
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
     if(!fgets(str, 150, imu_tactical)) {
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
   ins->data.sec=ins->time=imu_curr_meas.sec;
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
      //insw[i].data.sec=insw[i+1].data.sec;
      //insw[i].data=insw[i+1].data;
      insw[i]=insw[i+1];
      insbufferpass(insw+i, insw+i+1);
    }
    //insw[insgnssopt.insw-1].data.sec=insc->data.sec;
    //insw[insgnssopt.insw-1].data=insc->data; 
    insw[insgnssopt.insw-1]=*insc; 
    insbufferpass(insc, insw+insgnssopt.insw-1);
  }else{
    ins_buffinit(insw+ins_w_counter, insc->nx);
    //insw[insgnssopt.insw-1].data.sec=insc->data.sec;
    //insw[insgnssopt.insw-1].data=insc->data;
    insw[ins_w_counter]=*insc;
    insbufferpass(insc, insw+ins_w_counter);
  }
}

/* Stores last 10 values on a vector by sliding them*/
extern void windowSlider(double *vec, int size, double value){
  int i;
    for(i=0;i<size-1;i++){ 
      vec[i]=vec[i+1];
    }
    vec[size]=value; 
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


/* Coarse alignment -----------------------------------------------------------
* description: Coarse alignment of roll, pitch and yaw
* args   : um7pack_t* imu       IO    imu structure
*	   double* gn	I	local gravity
           double acc[3]	I	INS accelerations from previous epoch {xyz}
*          double v[3]	I	GNSS/or INS/GNSS previous state velocities {NED}
           gan[3]  I  gravity vector
           euler_angles[3]  IO  attitude angles in _nb frame {roll, pitch, yaw}
Reference: Shin (2005, pag.42) and Groves (2013, pages 198 to 199)
*-----------------------------------------------------------------------------*/
//void coarseAlign(um7pack_t* imu, double* gan)
void coarseAlign(double *acc, double *gyr, double *vel, double *gan, double *euler_angles)
{
 double r1[9], r2[9], aux[3], sign=1.0;
 double sineyaw, coseyaw;

 /* Is it a static or kinematic alignment? Use GNSS velocity or any velocity in imu->v */
 if ( norm(vel,3) < 0.5){
   printf("LEVELLING.STATIC, %lf, %lf, %lf\n", vel[0], vel[1],vel[2]);
                                    //static alignment
  /* Then we can solve for the roll(x(phi)) and pitch(y(theta)) angles using accelerometers */
  //if(acc[2]<0.0){sign=-1.0;}
  //accea[0]= sign*(asin((acc[1])/norm(gan,3)));
  //accea[1]=-sign*(asin((acc[0])/norm(gan,3)));

  /* Levelling from Groves (2013) */
  /* roll - phi angle */
  euler_angles[0]= atan2(-acc[1],-acc[2]);

  /* pitch - theta angle */
  euler_angles[1]= atan (acc[0]/sqrt(acc[1]*acc[1]+acc[2]*acc[2]));

  /* Gyrocompasing by Groves(2013,page 199)  */
  /* yaw - psi angle */
  sineyaw=-gyr[1]*cos(euler_angles[0])+gyr[2]*sin(euler_angles[0]);
  coseyaw= gyr[0]*cos(euler_angles[1])+gyr[1]*sin(euler_angles[0])*sin(euler_angles[1])+\
  gyr[2]*cos(euler_angles[0])*sin(euler_angles[1]);
  euler_angles[2]=atan2(sineyaw,coseyaw);

}else{                            //kinematic alignment
 /* Roll can be initialized to zero with an uncertainty of +-5 degrees,
 in most cases, on the road */
 printf("LEVELLING.KINE, %lf, %lf, %lf\n", vel[0], vel[1],vel[2] );

 euler_angles[0]=0.0;

 /* Levelling from Groves (2013) */
 /* roll - phi angle */
 //euler_angles[0]= atan2(-acc[1],-acc[2]); /* CORRECT IT LATER*/

 /* The pitch(y(theta)) and heading(z(psi)) angles using velocities */
 euler_angles[1]= atan(vel[2]/sqrt(vel[0]*vel[0]+vel[1]*vel[1]));
 euler_angles[2]= atan(vel[2]/vel[1]);
 /* Can be improved using an EKF with LHU model, see Shin (2005) chapter 3	*/
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
   for(i=0;i<3;i++){
    ins->pdata.fb0[i]=ins->data.fb0[i];
    ins->pdata.fb[i]=ins->data.fb[i];
    ins->pdata.wibb0[i]=ins->data.wibb0[i];
    ins->pdata.wibb[i]=ins->data.wibb[i];
   }
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

    ins->nb=ins->nx=ppptcnx(opt);
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
    int nx=ins->nx,info=0,nv, i, j;
    double *H,*v,*R,*x;  

    trace(3,"nhc:\n");
    printf("nhc:\n");

    H=zeros(2,nx); R=zeros(2,2);
    v=zeros(2,1); x=zeros(1,nx);
    for (i = 0; i < nx; i++) x[i]=1E-17;

    nv=bldnhc(opt,imu,ins->Cbe,ins->ve,nx,v,H,R);

    // i = xiV()+1;
    // for (i = xiV()+1; i < xiV() + 3; i++){
    //   R[i-(xiV()+1)+(i-(xiV()+1))*nv] = SQRT(ins->P[i + i * nx]);
    // }

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
    int nx=ins->nx,info=0, i, j;
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
    // i = xiV();
    // for (i = xiV(); i < xiV() + 3; i++){
    //   R[(i-xiV())+(i-xiV())*3] = SQRT(ins->P[i + i * nx]);
    // }        

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

/* Quasi_stationary IMU calibration  ************************************************/
extern int finealign(ins_states_t *ins,const insgnss_opt_t *opt,int flag)
{
    imuraw_t *imu=&ins->data;
    int nx=ins->nx,info=0, i;
    static int nz=0;
    double *x,*H,*R,*v,I[9]={-1,0,0,0,-1,0,0,0,-1};

    trace(3,"fine alignment:\n");
    printf("fine alignment:\n");

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
  double vn0, ve0, vn0std, ve0std, gyrx0, gyry0, gyrx0std, gyry0std,vnveRes;
   vn0 = -0.002743902895411; ve0 = -0.002219817510341;
   vn0std = 0.001958782893669; ve0std = 0.001549122618107;
   gyrx0 = -0.000152716; gyry0 = 0.000386451;
   gyrx0std = 0.000483409; gyry0std = 0.001271191;
   vnveRes = 10*(sqrt((vn0std*vn0std-ve0std*ve0std)));
   printf("static detection \n");

   //  if ( (fabs(ins->vn[0]) - vn0) <= 3*vn0std && (fabs(ins->vn[1]) - ve0) <= 3*ve0std ) {
    if ( norm(ins->vn, 2) < 1.0 ) {
       zvu_counter++;
       printf("static detection: 1 %lf\n", ins->time);
      }else {zvu_counter=0;printf("static detection: 0 %lf\n", ins->time);} 
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

   /* Generate IMU bias output record */
   fprintf(out_IMU_bias_file, "%lf %lf %lf %lf %.10lf %.10lf %.10lf %d\n", insc->time,
        insc->data.ba[0], insc->data.ba[1], insc->data.ba[2],
        insc->data.bg[0], insc->data.bg[1], insc->data.bg[2], 
        opt->Nav_or_KF);

   /* Generate KF uncertainty output record */
   if (opt->Nav_or_KF ){  
    fprintf(out_KF_SD_file, "%lf ", insc->time);
     for (i = 0; i < insc->nx; i++){
      for (j = 0; j < insc->nx; j++){
        (i==j?fprintf(out_KF_SD_file, "%lf ", SQRT(fabs(insc->P[i * insc->nx + j])) ):0);
       }
     }
   fprintf(out_KF_SD_file, "%d\n", opt->Nav_or_KF); 
   }  
 }

  if (opt->Nav_or_KF){
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
 
 }

 } 

/* GNSS covariance to KF weights from gnss solution */
void pvclkCovfromgnss(rtk_t *rtk, ins_states_t *ins){
  int nx=ins->nx;

  printf("GNSS vel cov: %lf %lf %lf\n",rtk->sol.qrv[0], rtk->sol.qrv[1], rtk->sol.qrv[2] );
  printf("GNSS pos cov: %lf %lf %lf\n",rtk->sol.qr[0], rtk->sol.qr[1], rtk->sol.qr[2] );

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

/* GNSS solution quality check, return 1 if it passess */  
int gnssQC(rtk_t *rtk, int nsat){
  int i,j;
  double gnss_xyz_cov[6], P[9], pos[3], Qposenu[9], gnss_ned_cov[3];

  printf("GNSS quality check:\n");

  /* Position covariance and velocity */
  for (j = 0; j < 6; j++) gnss_xyz_cov[j]=rtk->sol.qr[j];
  P[0]     =rtk->sol.qr[0]; /* xx or ee */
  P[4]     =rtk->sol.qr[1]; /* yy or nn */
  P[8]     =rtk->sol.qr[2]; /* zz or uu */
  P[1]=P[3]=rtk->sol.qr[3]; /* xy or en */
  P[5]=P[7]=rtk->sol.qr[4]; /* yz or nu */
  P[2]=P[6]=rtk->sol.qr[5]; /* zx or ue */
  
  covenu(pos,P,Qposenu);
  /* Using it where? */
  gnss_ned_cov[0]=Qposenu[4];
  gnss_ned_cov[1]=Qposenu[0];
  gnss_ned_cov[2]=Qposenu[8];

  printf("GNSS quality check: %lf, %lf Nsat: %d\n",norm(gnss_ned_cov,2), rtk->sol.gdop[0], nsat);

  if (norm(gnss_ned_cov,2)<20.0){ //for SPP a reasonable value is 15 m, for others 5 m
    if (rtk->sol.gdop[0]<4.0) 
    {
      printf("GNSS quality check ok\n"); 
      return 1;      
    }
    return 1;
  }
  printf("GNSS quality check fail\n");
  return 0;               
}

/* Stationary detection based on gnss and rotation test from imu measurements ******************************
  if check=1, gnss velocity only, if check=0, imu rotation check only
  return: 1: static, 0: no static */
void statRotat(ins_states_t *ins, double gnss_time, int check){
  int i;

  printf("bias static and rotation detection \n");

  if(gnss_w_counter >= insgnssopt.gnssw ) i=insgnssopt.gnssw-1;
  else i=gnss_w_counter;

  printf("GNSS ACCEL: %lf %lf %lf\n", gnss_time, norm(solw[i].rr+3,2), norm(ins->data.fb0 ,2)); // for debugging


  if (check){
    /* Static detection by GNSS */
      /* Velocity threshold: as low as 0.0075 m/s per axis (norm=0.0129) for aviagtion-grade 
      and 0.5 m/s (norm=0.8660) for consumer-grade sensors  
      Using 0.3 m/s (norm=0.519615)*/
    if (norm(solw[i].rr+3,3)<0.519615){
    /* Static */
    if (staticInfo.static_counter<10) {
      staticInfo.vel_gnss[staticInfo.static_counter]=1;
      staticInfo.gnss_time[staticInfo.static_counter]=gnss_time;
    }else {
      for (i = 0; i < 9; i++) {
        staticInfo.vel_gnss[i]=staticInfo.vel_gnss[i+1];
        staticInfo.gnss_time[i]=staticInfo.gnss_time[i+1];
      }
      staticInfo.vel_gnss[9]=1;
      staticInfo.gnss_time[9]=gnss_time;
    }
  }else{
    /* Moving */
    if (staticInfo.static_counter<10) {
      staticInfo.vel_gnss[staticInfo.static_counter]=0;
      staticInfo.gnss_time[staticInfo.static_counter]=gnss_time;
    }else {
      for (i = 0; i < 9; i++) {
        staticInfo.vel_gnss[i]=staticInfo.vel_gnss[i+1];
        staticInfo.gnss_time[i]=staticInfo.gnss_time[i+1];
      }
      staticInfo.vel_gnss[9]=0;
      staticInfo.gnss_time[9]=gnss_time;
    }
  }
  staticInfo.static_counter++;
  }else{  /* Gyro turning test */
    if (norm(ins->data.wibb0,3)<MAXGYRO){
    /* Straight */
    if (staticInfo.static_counter<10) {
      staticInfo.gyros[staticInfo.static_counter]=1;
      staticInfo.ins_time[staticInfo.static_counter]=ins->time;
    }else {
      for (i = 0; i < 9; i++) {
        staticInfo.gyros[i]=staticInfo.vel_gnss[i+1];
        staticInfo.ins_time[i]=staticInfo.ins_time[i+1];
      }
      staticInfo.gyros[9]=1;
      staticInfo.ins_time[9]=ins->time;
    }
  }else{
    /* Turning */
    if (staticInfo.static_counter<10) {
      staticInfo.gyros[staticInfo.static_counter]=0;
      staticInfo.ins_time[staticInfo.static_counter]=ins->time;
    }else {
      for (i = 0; i < 9; i++) {
        staticInfo.gyros[i]=staticInfo.vel_gnss[i+1];
        staticInfo.ins_time[i]=staticInfo.ins_time[i+1];
      }
      staticInfo.gyros[9]=0;
      staticInfo.ins_time[9]=ins->time;
    }
  }
  staticInfo.gyro_counter++;   
 }
 
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
  int i, j, week, flag, core_count=0;
  double gnss_time, rr[3], ve[3];
  ins_states_t insc={{{0}}};
  prcopt_t *opt = &rtk->opt; 
 
 
  printf("\n *****************  CORE BEGINS *******************: %lf\n", time2gpst(rtk->sol.time,&week));

  /* Save gnss position and velocity */ 
  matcpy(rr,rtk->sol.rr,1,3);
  matcpy(ve,rtk->sol.rr+3,1,3);

  /* initialize ins/gnss parameter default uncertainty */ 
  //ig_paruncinit(&insgnssopt); 
  kf_par_unc_init(&insgnssopt);
  
    printf("Insc.pdata: %lf %lf %lf - %lf %lf %lf\n", insc.pdata.fb0[0],insc.pdata.fb0[1],\
    insc.pdata.fb0[2], insc.data.fb0[0],insc.data.fb0[1],insc.data.fb0[2]);

  /* initialize ins state */
  insinit(&insc, &insgnssopt, opt, n);
  
  /* Initialize time from GNSS */
  gnss_time=time2gpst(rtk->sol.time,&week);

  /* Feed gnss solution and measurement buffers */
  gnssbuffer(&rtk->sol, obs, n);

  /* Check if imu file is not at the end */
  if(!insgnssopt.ins_EOF) return;

  /* Static check with GNSS */
  statRotat(&insc, gnss_time, 1);

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
       insc.data=insw[insgnssopt.insw-1].data;
       insc.data.sec=insw[insgnssopt.insw-1].data.sec; 
       }else {
         /* First epoch */
        if(ins_w_counter>1) {
          //printf("Before filling ins buffer: %d\n", ins_w_counter-1);
          //print_ins_pva(insw+ins_w_counter-1);
          insc=insw[ins_w_counter-1];
          insc.data=insw[ins_w_counter-1].data;
          insc.data.sec=insw[ins_w_counter-1].data.sec;
        }
       }  

    /* input ins */ 
    if(!inputimu(opt, &insc, week)) {insgnssopt.ins_EOF=0;printf(" ** End of imu file **\n"); return;}

    printf("Insc.pdata1: %lf %lf %lf - %lf %lf %lf\n", insc.pdata.fb0[0],insc.pdata.fb0[1],\
    insc.pdata.fb0[2], insc.data.fb0[0],insc.data.fb0[1],insc.data.fb0[2]);

    /* Propagation time */
    insc.proptime = gnss_time-0.5;

    /* Output raw INS */
    if(insc.pdata.sec > 0.0 ){ 
      outputrawimu(&insc);
    }
    /* Ins time propagation */
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

      /* GNSS solution quality check */
      if(flag) flag=gnssQC(rtk, n); 

      //if (flag)  

      /* Integration */ 
      printf("GNSS time and PROP time: %lf %lf\n", gnss_time, insc.proptime );
      TC_INS_GNSS_core1(rtk, obs, n, nav, &insc, &insgnssopt, flag);

      
      /* Non-holonomic constraints */
      if(insc.pdata.sec > 0.0 ){ 
        printf("nhc update\n");
        nhc(&insc,&insgnssopt);  
      }

      /* Static and rotation detection - It modifies staticInfo structure */
      printf("Test: %d %lf\n",insgnssopt.Nav_or_KF, fabs(gnss_time-insc.ptctime));
      statRotat(&insc, gnss_time, 0);


      printf("Static by gnss: %lf %d\n", gnss_time, (staticInfo.static_counter>10?staticInfo.vel_gnss[9]:staticInfo.vel_gnss[staticInfo.static_counter]));
      printf("Straight by Rotat: %lf %d\n", insc.time, (staticInfo.static_counter>10?staticInfo.gyros[9]:staticInfo.gyros[staticInfo.static_counter]));

      detstc(&insc);  
      
      /* Zero velocity update */  
      //if (zvu_counter>10) {  
      printf("Stat condition: %d \n", staticInfo.static_counter>10?staticInfo.vel_gnss[9]:staticInfo.vel_gnss[staticInfo.static_counter]);
      if((staticInfo.static_counter>10?staticInfo.vel_gnss[9]:staticInfo.vel_gnss[staticInfo.static_counter]) ){

        /* If straight and static do a Fine alignment */
        if ((staticInfo.static_counter>10?staticInfo.gyros[9]:staticInfo.gyros[staticInfo.static_counter])){
          //finealign(&insc,&insgnssopt,1);
        }
        
        /* Bias estimation */
        printf("ZVU UPDATE: 1 %lf\n", insc.time);
        /* Zero-velocity constraints */
        zvu(&insc,&insgnssopt,1); 
        zvu_counter=0;  
      }else{printf("ZVU UPDATE: 0 %lf\n", insc.time);}

      /* Output PVA, clock, imu bias solution     */ 
      if(insc.ptime>0.0){  
        outputinsgnsssol(&insc, &insgnssopt, &rtk->opt, n, obs);  
      }

    } // end If INS time ahead of GNSS condition 

    /* Re-initializing ins with gnss solution when integration occurs */
     if (insgnssopt.Nav_or_KF==1){
      for(i=0;i<3;i++) insc.re[i]=rtk->sol.rr[i];
      for(i=0;i<3;i++) insc.ve[i]=rtk->sol.rr[i+3];  
      /* Clock solution */
      //insc.dtr[0]=rtk->sol.dtr[0]*CLIGHT;  
      update_ins_state_n(&insc);  
      /* Gnss solution covariance to ins */ 
      //pvclkCovfromgnss(rtk, &insc);    
    }
 
    /* Zeroing closed-loop states */
     for (i=0;i<xnCl();i++) insc.x[i]=0.0;   

     /* Update ins states and measurements */
     insc.ptctime=gnss_time;    
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
FILE *res;

char residualsfname[]="../out/PPP_car_back.pos.stat"; //Residuals file
char tracefname[]="../out/trace.txt"; //trace file
int l=0,c;   
   
/* Global structures initialization */ 
insgnssopt.Tact_or_Low = 1;       /* Type of inertial, tact=1, low=0 */
insgnssopt.scalePN = 0;             /* use extended Process noise model */ 
insgnssopt.gnssw = 3; 
insgnssopt.insw = 10;  
insgnssopt.exphi = 0;            /* use precise system propagate matrix for ekf */
/* ins sthocastic process noises: */
insgnssopt.baproopt=INS_GAUSS_MARKOV;
insgnssopt.bgproopt=INS_GAUSS_MARKOV; 
insgnssopt.saproopt=INS_GAUSS_MARKOV;
insgnssopt.sgproopt=INS_GAUSS_MARKOV;
solw=(sol_t*)malloc(sizeof(sol_t)*insgnssopt.gnssw);  /* gnss solution structure window size allocation */
insw=(ins_states_t*)malloc(sizeof(ins_states_t)*insgnssopt.insw);   /* ins states window size allocation */
insgnssopt.ins_EOF=1;

/* Residuals structure */
resid.nv_w=10;
resid.data=(res_epoch_t*)malloc(sizeof(res_epoch_t)*resid.nv_w);  /* residuals structure window size allocation */

strcpy(filopt.trace,tracefname); 
   
/* Global TC_KF_INS_GNSS output files    */     
out_PVA=fopen("../out/out_PVA.txt","w"); 
out_clock_file=fopen("../out/out_clock_file.txt","w");  
out_IMU_bias_file=fopen("../out/out_IMU_bias.txt","w");
out_tropo_file=fopen("../out/out_tropo_bias.txt","w");
out_amb_file=fopen("../out/out_amb_bias.txt","w");
out_KF_state_error=fopen("../out/out_KF_state_error.txt","w");
out_KF_SD_file=fopen("../out/out_KF_SD.txt","w");
out_raw_fimu=fopen("../out/out_raw_imu.txt","w");
out_KF_residuals=fopen("../out/out_KF_residuals.txt","w");
imu_tactical=fopen("../data/19032019/imu_ascii.txt", "r");

rewind(imu_tactical);                                                   

/* PPP-Kinematic  Kinematic Positioning dataset  GPS+GLONASS */  
//char *argv[] = {"./rnx2rtkp", "../data/16102018/CAR_2890.18O", "../data/16102018/BRDC00IGS_R_20182890000_01D_MN.nav", "../data/16102018/grm20232.clk","../data/16102018/grm20232.sp3", "-o", "../out/PPP.pos", "-k", "../config/opts3.conf", "-x", "5"};

/* PPP-Kinematic  Kinematic Positioning dataset  March 21, 2019 GPS */
//char *argv[] = {"./rnx2rtkp", "../data/APS_center.19O", "../data/BRDC00WRD_S_20190780000_01D_MN.rnx", "../data/igs20452.clk","../data/igs20452.sp3", "-o", "../out/PPP.pos", "-k", "../config/opts3.conf", "-x", "5"};

char *argv[11] = {"./rnx2rtkp", "../data/19032019/observations.rnx", "../data/19032019/navigation.nav", "../data/19032019/orbit.sp3","../data/19032019/clock.clk", "-o", "../out/PPP_march21.pos", "-k", "../config/opts3.conf", "-x", "5"};

/* PPP-AR Kinematic   
char *argv[] = {"./rnx2rtkp", "../data/SEPT2640.17O", "../data/grg19674.*", "../data/SEPT2640.17N", "-o", "../out/exp1_PPP_amb_mod_constr.pos", "-k", "../config/opts3.conf"};*/

 //argc= sizeof(comlin) / sizeof(char);  
   int argc = 11;

    prcopt.mode  =PMODE_KINEMA;
  //  prcopt.navsys=SYS_GPS|SYS_GLO;   
    prcopt.refpos=1;
    prcopt.glomodear=1;
    solopt.timef=0;
    sprintf(solopt.prog ,"%s ver.%s",PROGNAME,VER_RTKLIB);
    sprintf(filopt.trace,"%s.trace",PROGNAME);

/* Declarations from rnx2rtkp */

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
  
   /* Start rnx2rtkp processing  ---------*/ 
  ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
  if (!ret) fprintf(stderr,"%40s\r","");
 
  for (i=0;i<insgnssopt.insw;i++) insfree(insw+i);  
  free(resid.data);
  free(solw); free(insw); 

 /* ins navigation only */
 //imu_tactical_navigation(imu_tactical); 

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

  char posfile[]="../out/out_PVA.txt"; 
  imuposplot(posfile);                

/* INS/GNSS plots  */ 
//  char out_raw_fimu[]="../out/out_raw_imu.txt";      
//  imuaccplot(out_raw_fimu);
// imugyroplot(out_raw_fimu);
// char gyrofile[]="../out/out_PVA.txt";   
// imueulerplot(gyrofile);
// char velfile[]="../out/out_PVA.txt"; 
// imuvelplot(velfile);
// char posfile[]="../out/out_PVA.txt"; 
// imuposplot(posfile);
// char imu_bias[]="../out/out_IMU_bias.txt";        
// imuaccbiasplot(imu_bias);  
// imugyrobiasplot(imu_bias);
// char imu_KF_stds[]="../out/out_KF_SD.txt"; 
// KF_att_stds_plot(imu_KF_stds);
// KF_vel_stds_plot(imu_KF_stds);
// KF_pos_stds_plot(imu_KF_stds);
// char imu_KF_clk[]="../out/out_clock_file.txt"; 
// KF_clock_plot(imu_KF_clk);
// char imu_KF_res[]="../out/out_KF_residuals.txt"; 
// KF_residuals_plot(imu_KF_res);
// char tropo_file[]="../out/out_tropo_bias.txt";
// tropo_plot(tropo_file);
//  char amb_file[]="../out/out_amb_bias.txt"; 
//  amb_plot(amb_file);                                  

 printf("\n\n SUCCESSFULLY EXECUTED!  \n\n");
 return;
}
