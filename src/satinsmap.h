/*------------------------------------------------------------------------------
* satinsmap.h : satinsmap constants, types and function prototypes
*
*          Copyright (C) 2015-2018 by Cavalheri, E. P. All rights reserved.
*
* version : $Revision: 1.1 $ $Date:  $
* history : 2017/10/09 1.0  satinsmap ver.1.0.0
*-----------------------------------------------------------------------------*/
#ifndef SATINSMAP_H
#define SATINSMAP_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <pthread.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif

/* constants -----------------------------------------------------------------*/

#define	SPC		0.1	/* Lane points spacing in meters (m)	*/
#define WBS		10*SPC	/* Whole buffer size in meters (m)	*/
#define BUFFSIZE	25 /* (WBS/SPC/2) Buffer search half size of vector
 (if buffsize/2=WBS(m)/SPC(m)*s), to convert into vector positions*/

/* Earth parameters ----------------------------------------------------------*/
#define e_2 ((2*FE_WGS84)-(FE_WGS84*FE_WGS84))
#define RN(lat) (RE_WGS84/sqrt(1-e_2*(sin(lat)*sin(lat)))) /* Prime vertical radii */
#define RM(lat) ( RE_WGS84*(1-e_2) / pow( 1-e_2*( sin(lat)*sin(lat) ),(3/2) ) ) /* Merdidian radius of curvature */
#define reeS(lat) (RN(lat)*sqrt(cos(lat)*cos(lat)+(1-e_2)*(1-e_2)*sin(lat)*sin(lat))) /* Geocentric radius at the surface*/
/* Normal gravity by Schwars and Wei (2000), p.30 from Shin (2001) */
#define	a1gn	9.7803267715
#define	a2gn	0.0052790414
#define	a3gn	0.0000232718
#define	a4gn	-0.0000030876910891
#define	a5gn	0.0000000043977311
#define	a6gn	0.0000000000007211
#define gn(lat,h) ( a1gn*(1+a2gn*(sin(lat)*sin(lat))+a3gn*(sin(lat)*sin(lat)*sin(lat)*sin(lat))) + h*(a4gn+a5gn*(sin(lat)*sin(lat))) + a6gn*h*h   ) /* Normal gravity on the ellipsoidal surface WGS84 - Defense Mapping Agency */

/* coordinate rotation matrices (Jekeli, 2001; Shin, 2001)--------------------*/
/* Earth-Centered Inertial to Earth-Centered Fixed (i to e-frame)
   where t is the angle (g.gg) and X is the output Rotation matrix	*/
#define Ci2e(t,X) do { \
 (X)[0]=cos(t); (X)[1]=sin(t); (X)[3]=-sin(t);(X)[4]=cos(t);\
 (X)[2]=(X)[5]=(X)[6]=(X)[7]=(X)[8]=0;\
} while (0) //In this case t=wie*t, earth rotation and the time since start of navigation

/* Earth-Centered Fixed to navigation frame (NED) (e to n-frame)	*/
#define Ce2n(lat,lon,X) do {\
 (X)[0]=-cos(lon)*sin(lat); (X)[1]=-sin(lon)*sin(lat); (X)[2]= cos(lat);\
 (X)[3]=-sin(lon);          (X)[4]=cos(lon);           (X)[5]=0;\
 (X)[6]=-cos(lon)*cos(lat); (X)[7]=-sin(lon)*cos(lat); (X)[8]=-sin(lat);\
} while (0)
/* n-frame to e-frame	*/
#define Cn2e(lat,lon,X) do {\
 (X)[0]=-cos(lon)*sin(lat); (X)[1]=-sin(lon); (X)[2]= -cos(lon)*cos(lat);\
 (X)[3]=-sin(lon)*sin(lat); (X)[4]=cos(lon);  (X)[5]=-sin(lon)*cos(lat);\
 (X)[6]=cos(lat);           (X)[7]=0;         (X)[8]=-sin(lat);\
} while (0)

/* Body-frame (x(nav dir.),y(dext. sys.), down(vertical)) to navigation frame (NED) (e to n-frame)	*/
#define Cb2n(phi,theta,psi,X) do {\
 (X)[0]=cos(theta)*cos(psi);(X)[1]=-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);(X)[2]=sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);\
 (X)[3]=cos(theta)*sin(psi);(X)[4]=cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);(X)[5]=-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi); \
 (X)[6]=-sin(theta);(X)[7]=sin(phi)*cos(theta);(X)[8]=cos(phi)*cos(theta); \
} while (0)
/* n-frame to b-frame */
#define Cn2b(phi,theta,psi,X) do {\
 (X)[0]=cos(theta)*cos(psi);(X)[1]=cos(theta)*sin(psi);(X)[2]=-sin(theta);\
 (X)[3]=-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);(X)[4]=cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);(X)[5]=sin(phi)*cos(theta); \
 (X)[6]=sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);(X)[7]=-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);(X)[8]=cos(phi)*cos(theta); \
} while (0)

/* R1, R2 and R3 rotation matrices	*/
#define R1(phi,X) do {\
 (X)[0]=1; (X)[1]=0;        (X)[2]=0;\
 (X)[3]=0; (X)[4]=cos(phi); (X)[5]=-sin(phi);\
 (X)[6]=0; (X)[7]=sin(phi); (X)[8]=cos(phi);\
} while (0)

#define R2(theta,X) do {\
 (X)[0]=cos(theta);  (X)[1]=0; (X)[2]=sin(theta);\
 (X)[3]=0;           (X)[4]=1; (X)[5]=0;\
 (X)[6]=-sin(theta); (X)[7]=0; (X)[8]=cos(theta);\
} while (0)

#define R3(psi,X) do {\
 (X)[0]=cos(psi); (X)[1]=-sin(psi);(X)[2]=0;\
 (X)[3]=sin(psi); (X)[4]=cos(psi); (X)[5]=0;\
 (X)[6]=0;        (X)[7]=0;        (X)[8]=0;\
} while (0)

/* type definitions ----------------------------------------------------------*/

typedef struct {        /* Reference Lane record */
    long int	size;	/* Lane file size in bytes, or lines */
    int sizeline; /* Size of the lines in file  */
    int sizebuffer; /* Amount of useful data in buffer  */
    double	buffer [BUFFSIZE*2*3];	/* Buffer reference lane points in xyz coodinates [m] */
    int	currsearch;	/* Closest point position in the search buffer in bytes */
} lane_t;

typedef struct {        /* UM7 sensor package records */
    int count; 		/* Flag for type of data, 0:gyro,1:accelerometer,2:magnetometer */
    float internal_time; /* amount of time in seconds since the sensor was on */
    double sec;		/*  If GPS is connected, this is synchronized to UTC time of day
			(in seconds of the week) 		*/
    int gpsw;		/* GPS week	*/
    gtime_t time;	/* GPST time	*/
    double Ba[6]; /* Accelerometer bias and scale factors {bax,bay,baz,sax,say,saz} */
    double Bg[6]; /* Gyroscope bias and scale factors {bgx,bgy,bgz,sgx,sgy,sgz} */
    double g[3];	/* Gyro data for xyz in B-frame-components (rad/sec)	*/
    double aea[3];	/* Attitude euler angles for xyz{roll,pitch,yaw} in I-frame-components (rad)	*/
    double a[3];	/* Accelerometer data for xyz in B-frame-components in gravities (m/s^2)	*/
    double m[3];	/* Magnetometer data for xyz in B-frame-components in unit-form (unitless)		*/
    double v[3],r[3];	/* Velocity (m/s) and position (m) estimates in I-frame	*/
    char checksum;	/* hex representation of the single-byte checksum of the packet		*/
    int status;		/* If it received a complete package=1, else =0	*/
    float tacc,tgyr,tmag; /* Accelerometer, gyroscopes, and magnetometer reading times */
} um7pack_t;

typedef struct {        /* IMU sensor measurements */
    float sec;		/* amount of time in seconds since the sensor was on		*/
    gtime_t time;	/* GPS UTC time	*/
    double wibb[3];	/* rate-gyro b-frame {wx,wy,wz} (deg/sec)	*/
    double stdw[3];	/* rate-gyro stds */
    double fb[3];	/* Specific-force b-frame {fx,fy,fz} (g or m/s^2)	*/
    double stdf [3];	/* specific-force stds */
    double ba[6]; /* Accelerometer bias and drift {bax,bay,baz,dbax,dbay,dbaz} */
    double stdba[3];	/* Acc. bias stds */
    double bg[6]; /* Gyroscope bias and drift {bgx,bgy,bgz,dbgx,dbgy,dbgz} */
    double stdbg[3];	/* Gyro. bias stds */
    double sma[6]; /* Accelerometer scale factors and cross-coupling {sax,say,saz,max,may,maz} */
    double stdsma[6];	/* Acc. scale factors bias stds */
    double smg[6]; /* Gyroscope scale factors and cross-coupling {sgx,sgy,sgz,mgx,mgy,mgz} */
    double stdsmg[6];	/* Gyro. scale factors bias stds */
} imuraw_t;

typedef struct {        /* Position, velocity and attitude structure (PVA) */
    double sec;		/* amount of time in seconds since the sensor was on		*/
    float t_s;    /* State-time estimation */
    gtime_t time;	/* GPS UTC time	*/
    double clock_offset_drift[2]; /* Esimated receiver clock offset and drift */
    double r[3];	/* Position in n-frame {n,e,d} or {lat,lon,h} (m or rads)	*/
    double v[3]; /* Velocity in n-frame {vn,ve,vd} (m/s) */
    double A[3]; /* Attitude vector xyz {roll,pitch,yaw} (rad) */
    double Cbn[9]; /* Attitude matrix {rad}*/
    double P[9]; /* Attitude [3], velocity [3] and position [3] uncertainties (m^2 / (m/s)^2 rad)*/
    double out_errors[17];  /*(weights) Attitude, velocity, position, acc. and gyro bias, clock offset and drift errors */
    double out_IMU_bias_est[6]; /* Estimated IMU acc. and gyro. bias*/
} pva_t;

/* global variables ----------------------------------------------------------*/
extern lane_t lane;
extern imuraw_t imu_obs_global;
extern pva_t pva_global;
extern pva_t pvagnss;
extern FILE *fp_lane;
extern FILE *fimu;
extern FILE *out_raw_fimu;
extern FILE *out_PVA;
extern FILE *out_clock_file;
extern FILE *out_IMU_bias_file;
extern FILE *out_KF_SD_file;
extern FILE *imu_tactical;
//fp_lane=fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/data/Lanes1_XYZ_2016_ITRF08.txt","r");
       /* Lane coordinate file pointer */

/* function declaration ------------------------------------------------------*/

/* utilities functions */
extern void getstr1(const char *input, int offset, int len, char *res);
extern double dist(const double *rs, const double *rr);
extern int min(double *a, int n, double *min);
extern int resmin(double *a, int n, double *min);
extern int listfiles(void);
extern int showmsg(char *format, ...);
extern void settspan(gtime_t ts, gtime_t te);
extern void settime(gtime_t time);
extern void vec2skew (double *vec, double *W);

/* map-matching functions	*/
void extern inibuff ();
extern int match (rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
  double *mmcand);

/* System positioning */
extern void core(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);

/* imu-mems functions --------------------------------------------------------*/
extern void inssysmatrix(double *PHI, double *G, int nx, pva_t *pva,
  imuraw_t *imuobs, double dt);
extern void inssysmatrix1(double *PHI, int nx, pva_t *pva, imuraw_t *imuobsc,
  double dt);
extern void sysnoise(double *Q, double dt, int nx);
extern void sysnoise1(double *Q, double dt, int nx);
extern void attbiasP0(imuraw_t *imu, pva_t *pvap, double *Pp, int nx);
extern void posvelP0 (double *P, int nx, rtk_t *rtk);
extern void initialize_P0(double *P_matrix, int n);
extern void measvec(double *z, double* xyz_ini_pos, pva_t *pvap, imuraw_t *imu, double *l);
extern void measmatrixH (double *H, int n, int nx, pva_t *pva);
extern void measnoiseR (double *R, int nv, double* gnss_xyz_ini_cov);
extern void imufileback();
extern int insnav (rtk_t *rtk, um7pack_t *imu, pva_t *xp, imuraw_t *imuobsp);
extern int insnav1(double* xyz_ini_pos, double* gnss_enu_vel, float ini_pos_time, um7pack_t *imu, pva_t *pvap, imuraw_t *imuobsp);
extern void insgnssLC (double* gnss_xyz_ini_pos, double* gnss_xyz_ini_cov,\
   double* gnss_enu_vel, double ini_pos_time, um7pack_t *imu, pva_t *pvap, imuraw_t *imuobsp);
extern void insmap ();
extern void ins_LC (double* gnss_xyz_ini_pos, double* gnss_xyz_ini_cov, double* gnss_enu_vel, double ini_pos_time, um7pack_t *imu, pva_t *pvap, imuraw_t *imuobsp);

/* plot functions ------------------------------------------------------------*/
extern void mapmatchplot ();
extern void resprint(FILE* res);
extern void imuaccplot(char* filename);
extern void imugyroplot(char* filename);
extern void imumagplot(char* filename);
extern void imueulerplot(char* filename);
extern void imuposplot(char* filename);
extern void imuUpplot(char* filename);
extern void imuvelplot(char* filename);
extern void imuaccbiasplot(char* filename);
extern void imugyrobiasplot(char* filename);
extern void KF_att_stds_plot(char* filename);
extern void KF_vel_stds_plot(char* filename);
extern void KF_pos_stds_plot(char* filename);
extern void imugyrobiasplot(char* filename);
extern void KF_clock_plot(char* filename);
extern void nsat ();

/* geodetic and positioning functions ----------------------------------------*/
extern d2lgs(double lat, double h, double* pos, double* e);
extern void undiffppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);

#ifdef __cplusplus
}
#endif
#endif /* SATINSMAP_H */
