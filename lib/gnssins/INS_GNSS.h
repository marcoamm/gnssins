/*------------------------------------------------------------------------------
* INS_GNSS.h : INS_GNSS constants, types and function prototypes
*
*          Copyright (C) 2015-2018 by Cavalheri E.P. and Mendonca M.A.
*
* version : $Revision: 1.1 $ $Date:  $
* history : 2018/11/01 1.0  INS_GNSS ver.1.0.0
*-----------------------------------------------------------------------------*/
#ifndef INS_GNSS_H
#define INS_GNSS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "rtklib.h"

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

/* constants  */
#define GME         3.986004415E+14 /* earth gravitational constant */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */

#define MIN_NSAT_SOL 4              /* min satellite number for solution */

#define MAXDT 3600.0 /* max time difference for ins-gnss coupled */
#define MAXVAR 1E10  /* max variance for reset covariance matrix */
#define STD_POS 2.5  /* pos-std-variance of gnss positioning (m)*/
#define STD_VEL 0.1  /* vel-std-variance of gnss positioning (m/s) */

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(10.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define MAXINOP 1000.0      /* max innovations for updates ins states */
#define MAXINOV 100.0       /* max innovations for updates ins states */
#define MAXVEL 5.0          /* max velocity for refine attitude */
#define MINVEL 3.0          /* min velocity for refine attitude */
#define MAXGYRO (5.0 * D2R) /* max gyro measurement for refine attitude */
#define MAXANG 3.0          /* max angle for refine attitude */
#define MAXUPDTIMEINT 60.0  /* max update time internal is same as correlation time */
#define MAXSYNDIFF 1.0      /* max time difference of ins and gnss time synchronization */
#define CORRETIME 360.0     /* correlation time for gauss-markov process */
#define MAXROT (10.0 * D2R) /* max rotation of vehicle when velocity matching alignment */
#define MAXSOLS 5           /* max number of solutions for reboot lc  */
#define MAXDIFF 10.0        /* max time difference between solution */
#define MAXVARDIS (10.0)    /* max variance of disable estimated state */
#define REBOOT 1            /* ins loosely coupled reboot if always update fail */
#define USE_MEAS_COV 0      /* use measurement covariance {xx,yy,zz,xy,xz,yz} for update, not just {xx,yy,zz} */
#define CHKNUMERIC 1        /* check numeric for given value */
#define NOINTERP 0          /* no interpolate ins position/velocity when gnss measurement if need */
#define COR_IN_ROV 0        /* correction attitude in rotation vector,otherwise in euler angles */

#define NNAC 3 /* number of accl process noise */
#define NNGY 3 /* number of gyro process noise */
#define NNBA 3 /* number of accl bias process noise */
#define NNBG 3 /* number of gyro bias process noise */
#define NNPX (NNBA + NNBG + NNAC + NNGY)
#define INAC 0                    /* index of accl process noise */
#define INGY NNAC                 /* index of gyro process noise */
#define INBA (NNAC + NNGY)        /* index of accl bias process noise */
#define INBG (NNAC + NNGY + NNBA) /* index of gyro bias process noise */

#define CORRETIME 360.0 /* correlation time for gauss-markov process */
#define ORDERS 10       /* orders of approximate exponential of matrix */

/* Earth parameters */
#define RE_GRS80    6378137.0           /* earth semimajor axis (GRS80) (m) */
#define FE_GRS80    (1.0/298.257222101) /* earth flattening (GRS80) */
#define mu      3.986004418E14 /* WGS84 Earth gravitational constant (m^3 s^-2) */
#define J_2     1.082627E-3 /* WGS84 Earth's second gravitational constant */
#define e_2	((2*FE_GRS80)-(FE_GRS80*FE_GRS80)) 				/* second eccentricity */
#define RN(lat) (RE_GRS80/sqrt(1-e_2*(sin(lat)*sin(lat)))) 			/* Prime vertical radii */
#define RM(lat) ( RE_GRS80*(1-e_2) / pow( 1-e_2*( sin(lat)*sin(lat) ),(3/2) ) ) /* Merdidian radius of curvature */
#define reeS(lat) (RN(lat)*sqrt(cos(lat)*cos(lat)+(1-e_2)*(1-e_2)*sin(lat)*sin(lat))) /* Geocentric radius at the surface*/
/* Normal gravity by Schwars and Wei (2000), p.30 from Shin (2001) */
#define	a1gn	9.7803267715
#define	a2gn	0.0052790414
#define	a3gn	0.0000232718
#define	a4gn	-0.0000030876910891
#define	a5gn	0.0000000043977311
#define	a6gn	0.0000000000007211
#define gn(lat,h) ( a1gn*(1+a2gn*(sin(lat)*sin(lat))+a3gn*(sin(lat)*sin(lat)*sin(lat)*sin(lat))) + h*(a4gn+a5gn*(sin(lat)*sin(lat))) + a6gn*h*h   ) /* Normal gravity on the ellipsoidal surface WGS84 - Defense Mapping Agency */

/* Coordinate rotation matrices (Jekeli, 2001; Shin, 2001) */

/* Earth-Centered Inertial to Earth-Centered Fixed (i to e-frame)
   where t is the angle (g.gg) and X is the output Rotation matrix	*/
#define Ci2e(t,X) do { \
 (X)[0]=cos(t); (X)[1]=sin(t); (X)[3]=-sin(t);(X)[4]=cos(t);\
 (X)[2]=(X)[5]=(X)[6]=(X)[7]=0;\
 (X)[8]=1;\
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

typedef struct {        /* initialization of pos., vel. and attitude rrors */
  double  delta_r_eb_n[3];     /* position error resolved along NED (m) */
  double  delta_v_eb_n[3];     /* velocity error resolved along NED (m/s) */
  double  delta_eul_nb_n[3];   /* attitude error as NED Euler angles (rad) */
} initialization_errors;

typedef struct {        /* IMU errors */
  double  delta_r_eb_n[3];     /* position error resolved along NED (m) */
  double   b_a[3];              /* Accelerometer biases (m/s^2) */
  double   b_g[3];              /* Gyro biases (rad/s) */
  double  M_a[9];              /* Accelerometer scale factor and cross coupling errors */
  double  M_g[9];              /* Gyro scale factor and cross coupling errors */
  double  G_g[9];              /* Gyro g-dependent biases (rad-sec/m) */
  double   accel_noise_root_PSD; /* Accelerometer noise root PSD (m s^-1.5) */
  double   gyro_noise_root_PSD; /* Gyro noise root PSD (rad s^-0.5) */
  float   accel_quant_level;   /* Accelerometer quantization level (m/s^2) */
  float   gyro_quant_level;    /* Gyro quantization level (rad/s) */
} IMU_errors;

typedef struct {      /* GNSS configuration */
  double epoch_interval;      /* Interval between GNSS epochs (s) */
  double init_est_r_ea_e[3];     /* Initial estimated position (m; ECEF) */
  int no_sat;              /* Number of satellites in constellation */
  double r_os;              /* Orbital radius of satellites (m) */
  float inclination;       /* Inclination angle of satellites (deg) */
  float const_delta_lambda; /* Longitude offset of constellation (deg) */
  float const_delta_t;     /* Timing offset of constellation (s) */
  float mask_angle;          /* Mask angle (deg) */
  float SIS_err_SD;        /* Signal in space error SD (m) */
  float zenith_iono_err_SD; /* Zenith ionosphere error SD (m) */
  float zenith_trop_err_SD; /* Zenith troposphere error SD (m) */
  float code_track_err_SD; /* Code tracking error SD (m) */
  float rate_track_err_SD; /* Range rate tracking error SD (m/s) */
  float rx_clock_offset;   /* Receiver clock offset at time=0 (m) */
  float rx_clock_drift;    /* Receiver clock drift at time=0 (m/s) */
} GNSS_config;

typedef struct {      /* Tighlty coupled Kalman Filter configuration */
 float init_att_unc[3];           /* Initial attitude uncertainty per axis (rad) */
 float init_vel_unc[3];           /* Initial velocity uncertainty per axis (m/s) */
 float init_pos_unc[3];           /* Initial position uncertainty per axis (m) */
 double init_pos_unc_ned[3];           /* Initial position uncertainty per axis (m) in NED */
 double init_b_a_unc;           /* Initial accel. bias uncertainty (m/s^2) */
 double init_b_g_unc;           /* Initial gyro. bias uncertainty (rad/s) */
 float init_clock_offset_unc;  /* Initial clock offset uncertainty per axis (m) */
 float init_clock_drift_unc;   /* Initial clock drift uncertainty per axis (m/s) */
 double gyro_noise_PSD;         /* Gyro noise PSD (rad^2/s) */
 double accel_noise_PSD;        /* Accelerometer noise PSD (m^2 s^-3) */
 double accel_bias_PSD;         /* Accelerometer bias random walk PSD (m^2 s^-5) */
 double gyro_bias_PSD;          /* Gyro bias random walk PSD (rad^2 s^-3) */
 float clock_freq_PSD;         /* Receiver clock frequency-drift PSD (m^2/s^3) */
 float clock_phase_PSD;        /* Receiver clock phase-drift PSD (m^2/s) */
 float pseudo_range_SD;        /* Pseudo-range measurement noise SD (m) */
 float range_rate_SD;          /* Pseudo-range rate measurement noise SD (m/s) */
} TC_KF_config;

typedef struct {      /* Loosely coupled Kalman Filter configuration */
 float init_att_unc[3];           /* Initial attitude uncertainty per axis (rad) */
 float init_vel_unc[3];           /* Initial velocity uncertainty per axis (m/s) */
 float init_pos_unc[3];           /* Initial position uncertainty per axis (m) */
 double init_pos_unc_ned[3];           /* Initial position uncertainty per axis (m) in NED */
 double init_b_a_unc;           /* Initial accel. bias uncertainty (m/s^2) */
 double init_b_g_unc;           /* Initial gyro. bias uncertainty (rad/s) */
 double gyro_noise_PSD;         /* Gyro noise PSD (rad^2/s) */
 double accel_noise_PSD;        /* Accelerometer noise PSD (m^2 s^-3) */
 double accel_bias_PSD;         /* Accelerometer bias random walk PSD (m^2 s^-5) */
 double gyro_bias_PSD;          /* Gyro bias random walk PSD (rad^2 s^-3) */
 float pos_meas_SD;        /* Position measurement noise SD per axis (m) */
 float vel_meas_SD;          /* Velocity measurement noise SD per axis (m/s) */
} LC_KF_config;

typedef struct {      /* Position velocity and attitude solution structure */
  double time;             /* time (sec) */
  double latitude;        /* latitude (rad) */
  double longitude;       /* longitude (rad) */
  double height;           /* height (m) */
  double ned_velocity[3];  /* north, east, down velocity (m/s) */
  double euler_angles[3];  /* roll, pitch, yaw angle of body w.r.t NED (rad) */
  double ve[3],re[3], ae[3]; /*states in ecef */
  double C_b_e[9];        /* Tranformation matrix from body-to-ECEF frame  */
  double C_b_n[9];
  int Nav_or_KF;          /* Navigation sol:0 KF Integrated sol:1   */
  double P[17*17];        /* Full TC KF weight matrix (it accomodates for the LC[15*15]) */
} PVAT_solution;


typedef struct {      /* Position velocity and attitude solution structure */
  unsigned char sat;            /* Satellite number */
  gtime_t time;       /* receiver sampling time (GPST) */
  double sec;		/* amount of time in seconds since the sensor was onor in week time	*/
  double tt;          /* time difference between current and previous (s) */
  double P[NFREQ+NEXOBS];             /* Pseudo-range measurements (m)  */
  double L[NFREQ+NEXOBS];    /* Pseudo-range rate measurements (m/s) or cycles*/
  double D[NFREQ+NEXOBS];     /* observation data doppler frequency (Hz) */
  double Sat_r_eb_e[3];    /* Satellite ECEF position (m) */
  double Sat_v_eb_e[3];    /* Satellite ECEF velocity (m/s)time (sec) */
  double gdop[4];         /* DOPs {GDOP,PDOP,HDOP,VDOP} */
} GNSS_measurements;

typedef struct {        /* IMU sensor measurements */
    double sec;		/* amount of time in seconds since the sensor was on or in week time*/
    gtime_t time;	/* GPS UTC time	*/
    double tt;          /* time difference between current and previous (s) */
    double omega_ib_b[3];	/* rate-gyro b-frame {wx,wy,wz} (rad/sec)	*/
    double omega_stds[3];	/* rate-gyro stds */
    double f_ib_b[3];	/* Specific-force b-frame {fx,fy,fz} (m/s^2)	*/
    double f_stds [3];	/* specific-force stds */
} INS_measurements;

/* global variables ----------------------------------------------------------*/

/* function declaration ------------------------------------------------------*/

/* utilities functions */

#ifdef __cplusplus
}
#endif
#endif /* INS_GNSS_H */
