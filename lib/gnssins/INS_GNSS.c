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
#include "INS_GNSS.h"
#include "../../src/satinsmap.h"

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

/* struct declaration -------------------------------------------------------*/


/* Functions declaration  -------------------------------------*/

/* number of estimated insppptc states ------------------------------------------------*/
extern int ppptcnx(const prcopt_t *opt)
{
  printf("ppprcnx: XnX: %d\n", xnX(opt));
  return xnX(opt);
}
/* initial ins-gnss coupled ekf estimated states and it covariance-----------
 * args  :  insopt_t *opt    I  ins options
 *          insstate_t *ins  IO ins states
 * return : none
 * note   : it also can initial ins tightly coupled
 * --------------------------------------------------------------------------*/
extern void initPNindex(const prcopt_t *opt)
{
  trace(3, "iniPNindex:\n");

  /* initial global states index and numbers */
  IA = xiA();
  NA = xnA();
  IV = xiV();
  NV = xnV();
  IP = xiP();
  NP = xnP();
  iba = xiBa();
  nba = xnBa();
  ibg = xiBg();
  nbg = xnBg();
  irc = xiRc();
  nrc = xnRc(opt);
  irr = xiRr(opt);
  nrr = xnRr();
  IT = xiTr(opt);
  NT = xnT(opt);
  IN = xiBs(opt, 1);
  NN = xnB();
}

/* Utility ----------------------------------------------------------*/
static inline int8_t sign(double val)
{
  /* Corresponding sign() function in Matlab
  From: https://forum.arduino.cc/index.php?topic=37804.0 solution*/
  if (val < 0.0)
    return -1;
  if (val > 0.0)
    return 1;
  return 0;
}

/* precise tropospheric model ------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, const double *azel,
                       const prcopt_t *opt, const double *x, double *dtdx,
                       double *var)
{
  const double zazel[] = {0.0, PI / 2.0};
  double zhd, m_h, m_w, cotz, grad_n, grad_e;

  /* zenith hydrostatic delay */
  zhd = tropmodel(time, pos, zazel, 0.0);

  /* mapping function */
  m_h = tropmapf(time, pos, azel, &m_w);

  if ((opt->tropopt == TROPOPT_ESTG || opt->tropopt == TROPOPT_CORG) && azel[1] > 0.0)
  {

    /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
    cotz = 1.0 / tan(azel[1]);
    grad_n = m_w * cotz * cos(azel[0]);
    grad_e = m_w * cotz * sin(azel[0]);
    m_w += grad_n * x[1] + grad_e * x[2];
    dtdx[1] = grad_n * (x[0] - zhd);
    dtdx[2] = grad_e * (x[0] - zhd);
  }
  dtdx[0] = m_w;
  *var = SQR(0.01);
  return m_h * zhd + m_w * (x[0] - zhd);
}

void Initialize_INS_GNSS_LCKF_consumer_grade_IMU(initialization_errors *initialization_errors,
                                                 IMU_errors *IMU_errors, GNSS_config *GNSS_config, LC_KF_config *LC_KF_config)
{
  int i;
  /*%SCRIPT Tightly coupled INS/GNSS demo:
%   Profile_1 (60s artificial car motion with two 90 deg turns)
%   Consumer-grade IMU
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% Created 12/4/12 by Paul Groves

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details     */

  /* Constants  */
  double deg_to_rad = 0.01745329252;
  double rad_to_deg = 1 / deg_to_rad;
  double micro_g_to_meters_per_second_squared = 9.80665E-6;

  /* CONFIGURATION */
  /* Input truth motion profile filename
input_profile_name = 'Profile_2.csv'; */
  /* Output motion profile and error filenames
output_profile_name = 'INS_GNSS_Demo_9_Profile.csv';
output_errors_name = 'INS_GNSS_Demo_9_Errors.csv';  */

  /* Attitude initialization error (deg, converted to rad; @N,E,D) */
  initialization_errors->delta_eul_nb_n[0] = -0.5 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[1] = 0.4 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[2] = 2 * deg_to_rad; /* rad */

  /* Attitude initialization error (deg, converted to rad; @N,E,D) From um7 datasheet - static values */
  initialization_errors->delta_eul_nb_n[0] = 1 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[1] = 1 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[2] = 3 * deg_to_rad; /* rad */

  /* Accelerometer biases (micro-g, converted to m/s^2; body axes) */
  IMU_errors->b_a[0] = 9000 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[1] = -13000 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[2] = 8000 * micro_g_to_meters_per_second_squared;
  /* Gyro biases (deg/hour, converted to rad/sec; body axes) */
  IMU_errors->b_g[0] = -180 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = 260 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = -160 * deg_to_rad / 3600;
  /* Gyro biases (deg/s, converted to rad/sec; body axes)  From um7 */
  IMU_errors->b_a[0] = 0.75 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[1] = 0.75 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[2] = 1.5 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_g[0] = 20 * deg_to_rad;
  IMU_errors->b_g[0] = 20 * deg_to_rad;
  IMU_errors->b_g[0] = 20 * deg_to_rad;
  /* Accelerometer scale factor and cross coupling errors (ppm, converted to
 unitless; body axes) */
  IMU_errors->M_a[0] = 50000 * 1E-6;
  IMU_errors->M_a[1] = -15000 * 1E-6;
  IMU_errors->M_a[2] = 10000 * 1E-6;
  IMU_errors->M_a[3] = -7500 * 1E-6;
  IMU_errors->M_a[4] = -60000 * 1E-6;
  IMU_errors->M_a[5] = 12500 * 1E-6;
  IMU_errors->M_a[6] = -12500 * 1E-6;
  IMU_errors->M_a[7] = 5000 * 1E-6;
  IMU_errors->M_a[8] = 20000 * 1E-6;
  /* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
  IMU_errors->M_g[0] = 40000 * 1E-6;
  IMU_errors->M_a[1] = -14000 * 1E-6;
  IMU_errors->M_g[2] = 12500 * 1E-6;
  IMU_errors->M_g[3] = 0 * 1E-6;
  IMU_errors->M_a[4] = -30000 * 1E-6;
  IMU_errors->M_g[5] = -7500 * 1E-6;
  IMU_errors->M_g[6] = 0 * 1E-6;
  IMU_errors->M_a[7] = 0 * 1E-6;
  IMU_errors->M_g[8] = -17500 * 1E-6;

  /* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
  IMU_errors->G_g[0] = 90 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[1] = -110 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[2] = -60 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[3] = -50 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[4] = 190 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[5] = -160 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[6] = 30 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[7] = 110 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[8] = -130 * deg_to_rad / (3600 * 9.80665);

  /* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) */
  IMU_errors->accel_noise_root_PSD = 1000 * micro_g_to_meters_per_second_squared;
  /* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) FROM um7 datasheet */
  IMU_errors->accel_noise_root_PSD = 400 * micro_g_to_meters_per_second_squared;
  /* Gyro noise root PSD (deg per root hour, converted to rad s^-0.5) */
  IMU_errors->gyro_noise_root_PSD = 1 * deg_to_rad / 60.0;
  /* Gyro noise root PSD (deg per second root Hz, converted to rad s^-0.5)  From um7 datasheet */
  IMU_errors->gyro_noise_root_PSD = 0.005 * deg_to_rad;
  /* Accelerometer quantization level (m/s^2) */
  IMU_errors->accel_quant_level = 1E-1;
  /* Gyro quantization level (rad/s) */
  IMU_errors->gyro_quant_level = 2E-3;

  /* Interval between GNSS epochs (s) */
  GNSS_config->epoch_interval = 0.5;

  /* Initial estimated position (m; ECEF) */
  GNSS_config->init_est_r_ea_e[0] = 0;
  GNSS_config->init_est_r_ea_e[1] = 0;
  GNSS_config->init_est_r_ea_e[2] = 0;

  /* Number of satellites in constellation */
  GNSS_config->no_sat = 30;
  /* Orbital radius of satellites (m) */
  GNSS_config->r_os = 2.656175E7;
  /* Inclination angle of satellites (deg) */
  GNSS_config->inclination = 55;
  /* Longitude offset of constellation (deg) */
  GNSS_config->const_delta_lambda = 0;
  /* Timing offset of constellation (s) */
  GNSS_config->const_delta_t = 0;

  /* Mask angle (deg) */
  GNSS_config->mask_angle = 10;
  /* Signal in space error SD (m) *Give residual where corrections are applied */
  GNSS_config->SIS_err_SD = 1;
  /* Zenith ionosphere error SD (m) *Give residual where corrections are applied */
  GNSS_config->zenith_iono_err_SD = 2;
  /* Zenith troposphere error SD (m) *Give residual where corrections are applied*/
  GNSS_config->zenith_trop_err_SD = 0.2;
  /* Code tracking error SD (m) *Can extend to account for multipath */
  GNSS_config->code_track_err_SD = 1;
  /* Range rate tracking error SD (m/s) *Can extend to account for multipath */
  GNSS_config->rate_track_err_SD = 0.02;
  /* Receiver clock offset at time=0 (m); */
  GNSS_config->rx_clock_offset = 10000;
  /* Receiver clock drift at time=0 (m/s); */
  GNSS_config->rx_clock_drift = 100;

  /* Initial attitude uncertainty per axis (deg, converted to rad) */
  for (i = 0; i < 3; i++)
    LC_KF_config->init_att_unc[i] = D2R * 2;
  /* Initial velocity uncertainty per axis (m/s) */
  for (i = 0; i < 3; i++)
    LC_KF_config->init_vel_unc[i] = 0.1;
  /* Initial position uncertainty per axis (m) */
  for (i = 0; i < 3; i++)
    LC_KF_config->init_pos_unc[i] = 10;
  /* Initial accelerometer bias uncertainty per instrument (micro-g, converted
/* to m/s^2) */
  LC_KF_config->init_b_a_unc = 10000 * micro_g_to_meters_per_second_squared;
  /* Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)*/
  LC_KF_config->init_b_g_unc = 200 * deg_to_rad / 3600;

  /* Gyro noise PSD (deg^2 per hour, converted to rad^2/s) */
  LC_KF_config->gyro_noise_PSD = 0.01 * 0.01;
  /* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
  LC_KF_config->accel_noise_PSD = 0.2 * 0.2;
  /* % NOTE: A large noise PSD is modeled to account for the scale-factor and
% cross-coupling errors that are not directly included in the Kalman filter
% model*/
  /* Accelerometer bias random walk PSD (m^2 s^-5) */
  LC_KF_config->accel_bias_PSD = 1.0E-5;
  /* Gyro bias random walk PSD (rad^2 s^-3) */
  LC_KF_config->gyro_bias_PSD = 4.0E-11;

  /* Pseudo-range measurement noise SD (m) */
  LC_KF_config->pos_meas_SD = 2.5;
  /* Pseudo-range rate measurement noise SD (m/s) */
  LC_KF_config->vel_meas_SD = 0.1;
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Initialize_INS_GNSS_TCKF_consumer_grade_IMU(initialization_errors *initialization_errors,
                                                 IMU_errors *IMU_errors, GNSS_config *GNSS_config, TC_KF_config *TC_KF_config)
{
  int i;
  /*%SCRIPT Tightly coupled INS/GNSS demo:
%   Profile_1 (60s artificial car motion with two 90 deg turns)
%   Consumer-grade IMU
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% Created 12/4/12 by Paul Groves

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details     */

  /* Constants  */
  double deg_to_rad = 0.01745329252;
  double rad_to_deg = 1 / deg_to_rad;
  double micro_g_to_meters_per_second_squared = 9.80665E-6;

  /* CONFIGURATION */
  /* Input truth motion profile filename
input_profile_name = 'Profile_2.csv'; */
  /* Output motion profile and error filenames
output_profile_name = 'INS_GNSS_Demo_9_Profile.csv';
output_errors_name = 'INS_GNSS_Demo_9_Errors.csv';  */

  /* Attitude initialization error (deg, converted to rad; @N,E,D) */
  initialization_errors->delta_eul_nb_n[0] = -0.5 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[1] = 0.4 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[2] = 2 * deg_to_rad; /* rad */

  /* Attitude initialization error (deg, converted to rad; @N,E,D) From um7 datasheet - static values */
  initialization_errors->delta_eul_nb_n[0] = 1 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[1] = 1 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[2] = 3 * deg_to_rad; /* rad */

  /* Accelerometer biases (micro-g, converted to m/s^2; body axes) */
  IMU_errors->b_a[0] = 9000 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[1] = -13000 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[2] = 8000 * micro_g_to_meters_per_second_squared;
  /* Gyro biases (deg/hour, converted to rad/sec; body axes) */
  IMU_errors->b_g[0] = -180 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = 260 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = -160 * deg_to_rad / 3600;
  /* Gyro biases (deg/s, converted to rad/sec; body axes)  From um7 */
  IMU_errors->b_a[0] = 0.75 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[1] = 0.75 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[2] = 1.5 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_g[0] = 20 * deg_to_rad;
  IMU_errors->b_g[0] = 20 * deg_to_rad;
  IMU_errors->b_g[0] = 20 * deg_to_rad;
  /* Accelerometer scale factor and cross coupling errors (ppm, converted to
 unitless; body axes) */
  IMU_errors->M_a[0] = 50000 * 1E-6;
  IMU_errors->M_a[1] = -15000 * 1E-6;
  IMU_errors->M_a[2] = 10000 * 1E-6;
  IMU_errors->M_a[3] = -7500 * 1E-6;
  IMU_errors->M_a[4] = -60000 * 1E-6;
  IMU_errors->M_a[5] = 12500 * 1E-6;
  IMU_errors->M_a[6] = -12500 * 1E-6;
  IMU_errors->M_a[7] = 5000 * 1E-6;
  IMU_errors->M_a[8] = 20000 * 1E-6;
  /* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
  IMU_errors->M_g[0] = 40000 * 1E-6;
  IMU_errors->M_a[1] = -14000 * 1E-6;
  IMU_errors->M_g[2] = 12500 * 1E-6;
  IMU_errors->M_g[3] = 0 * 1E-6;
  IMU_errors->M_a[4] = -30000 * 1E-6;
  IMU_errors->M_g[5] = -7500 * 1E-6;
  IMU_errors->M_g[6] = 0 * 1E-6;
  IMU_errors->M_a[7] = 0 * 1E-6;
  IMU_errors->M_g[8] = -17500 * 1E-6;

  /* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
  IMU_errors->G_g[0] = 90 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[1] = -110 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[2] = -60 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[3] = -50 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[4] = 190 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[5] = -160 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[6] = 30 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[7] = 110 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[8] = -130 * deg_to_rad / (3600 * 9.80665);

  /* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) */
  IMU_errors->accel_noise_root_PSD = 1000 * micro_g_to_meters_per_second_squared;
  /* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) FROM um7 datasheet */
  IMU_errors->accel_noise_root_PSD = 400 * micro_g_to_meters_per_second_squared;
  /* Gyro noise root PSD (deg per root hour, converted to rad s^-0.5) */
  IMU_errors->gyro_noise_root_PSD = 1 * deg_to_rad / 60.0;
  /* Gyro noise root PSD (deg per second root Hz, converted to rad s^-0.5)  From um7 datasheet */
  IMU_errors->gyro_noise_root_PSD = 0.005 * deg_to_rad;
  /* Accelerometer quantization level (m/s^2) */
  IMU_errors->accel_quant_level = 1E-1;
  /* Gyro quantization level (rad/s) */
  IMU_errors->gyro_quant_level = 2E-3;

  /* Interval between GNSS epochs (s) */
  GNSS_config->epoch_interval = 0.5;

  /* Initial estimated position (m; ECEF) */
  GNSS_config->init_est_r_ea_e[0] = 0;
  GNSS_config->init_est_r_ea_e[1] = 0;
  GNSS_config->init_est_r_ea_e[2] = 0;

  /* Number of satellites in constellation */
  GNSS_config->no_sat = 30;
  /* Orbital radius of satellites (m) */
  GNSS_config->r_os = 2.656175E7;
  /* Inclination angle of satellites (deg) */
  GNSS_config->inclination = 55;
  /* Longitude offset of constellation (deg) */
  GNSS_config->const_delta_lambda = 0;
  /* Timing offset of constellation (s) */
  GNSS_config->const_delta_t = 0;

  /* Mask angle (deg) */
  GNSS_config->mask_angle = 10;
  /* Signal in space error SD (m) *Give residual where corrections are applied */
  GNSS_config->SIS_err_SD = 1;
  /* Zenith ionosphere error SD (m) *Give residual where corrections are applied */
  GNSS_config->zenith_iono_err_SD = 2;
  /* Zenith troposphere error SD (m) *Give residual where corrections are applied*/
  GNSS_config->zenith_trop_err_SD = 0.2;
  /* Code tracking error SD (m) *Can extend to account for multipath */
  GNSS_config->code_track_err_SD = 1;
  /* Range rate tracking error SD (m/s) *Can extend to account for multipath */
  GNSS_config->rate_track_err_SD = 0.02;
  /* Receiver clock offset at time=0 (m); */
  GNSS_config->rx_clock_offset = 10000;
  /* Receiver clock drift at time=0 (m/s); */
  GNSS_config->rx_clock_drift = 100;

  /* Initial attitude uncertainty per axis (deg, converted to rad) */
  for (i = 0; i < 3; i++)
    TC_KF_config->init_att_unc[i] = D2R * 2.0;
  /* Initial velocity uncertainty per axis (m/s) */
  for (i = 0; i < 3; i++)
    TC_KF_config->init_vel_unc[i] = 0.1;
  /* Initial position uncertainty per axis (m) */
  for (i = 0; i < 3; i++)
    TC_KF_config->init_pos_unc[i] = 10.0;
  /* Initial accelerometer bias uncertainty per instrument (micro-g, converted
/* to m/s^2) */
  TC_KF_config->init_b_a_unc = 10000.0 * micro_g_to_meters_per_second_squared;
  /* Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)*/
  TC_KF_config->init_b_g_unc = 200.0 * deg_to_rad / 3600.0;
  /* Initial clock offset uncertainty per axis (m) */
  TC_KF_config->init_clock_offset_unc = 10.0;
  /* Initial clock drift uncertainty per axis (m/s) */
  TC_KF_config->init_clock_drift_unc = 0.1;

  /* Gyro noise PSD (deg^2 per hour, converted to rad^2/s) */
  TC_KF_config->gyro_noise_PSD = 0.01 * 0.01;
  /* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
  TC_KF_config->accel_noise_PSD = 0.2 * 0.2;
  /* % NOTE: A large noise PSD is modeled to account for the scale-factor and
% cross-coupling errors that are not directly included in the Kalman filter
% model*/
  /* Accelerometer bias random walk PSD (m^2 s^-5) */
  TC_KF_config->accel_bias_PSD = 1.0E-5;
  /* Gyro bias random walk PSD (rad^2 s^-3) */
  TC_KF_config->gyro_bias_PSD = 4.0E-11;
  /* Receiver clock frequency-drift PSD (m^2/s^3) */
  TC_KF_config->clock_freq_PSD = 1;
  /* Receiver clock phase-drift PSD (m^2/s) */
  TC_KF_config->clock_phase_PSD = 1;

  /* Pseudo-range measurement noise SD (m) */
  TC_KF_config->pseudo_range_SD = 2.5;
  /* Pseudo-range rate measurement noise SD (m/s) */
  TC_KF_config->range_rate_SD = 0.1;
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Initialize_INS_GNSS_LCKF_tactical_grade_IMU(initialization_errors *initialization_errors,
                                                 IMU_errors *IMU_errors, GNSS_config *GNSS_config, LC_KF_config *LC_KF_config)
{
  int i;
  /*%INS_GNSS_Demo_9
%SCRIPT Tightly coupled INS/GNSS demo:
%   Profile_2 (175s car)
%   Tactical-grade IMU
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% Created 27/5/12 by Paul Groves

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details   */

  /* Constants  */
  double deg_to_rad = 0.01745329252;
  double rad_to_deg = 1 / deg_to_rad;
  double micro_g_to_meters_per_second_squared = 9.80665E-6;

  /* CONFIGURATION */
  /* Input truth motion profile filename
input_profile_name = 'Profile_2.csv'; */
  /* Output motion profile and error filenames
output_profile_name = 'INS_GNSS_Demo_9_Profile.csv';
output_errors_name = 'INS_GNSS_Demo_9_Errors.csv';  */

  /* Attitude initialization error (deg, converted to rad; @N,E,D) */
  initialization_errors->delta_eul_nb_n[0] = -0.05 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[1] = 0.04 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[2] = 1 * deg_to_rad; /* rad */

  /* Accelerometer biases (micro-g, converted to m/s^2; body axes) */
  IMU_errors->b_a[0] = 900 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[1] = -1300 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[2] = 800 * micro_g_to_meters_per_second_squared;
  /* Gyro biases (deg/hour, converted to rad/sec; body axes) */
  IMU_errors->b_g[0] = -9 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = 13 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = -8 * deg_to_rad / 3600;
  /* Accelerometer scale factor and cross coupling errors (ppm, converted to
 unitless; body axes) */
  IMU_errors->M_a[0] = 500 * 1E-6;
  IMU_errors->M_a[1] = -300 * 1E-6;
  IMU_errors->M_a[2] = 200 * 1E-6;
  IMU_errors->M_a[3] = -150 * 1E-6;
  IMU_errors->M_a[4] = -600 * 1E-6;
  IMU_errors->M_a[5] = 250 * 1E-6;
  IMU_errors->M_a[6] = -250 * 1E-6;
  IMU_errors->M_a[7] = 100 * 1E-6;
  IMU_errors->M_a[8] = 450 * 1E-6;
  /* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
  IMU_errors->M_g[0] = 400 * 1E-6;
  IMU_errors->M_a[1] = -300 * 1E-6;
  IMU_errors->M_g[2] = 250 * 1E-6;
  IMU_errors->M_g[3] = 0 * 1E-6;
  IMU_errors->M_a[4] = -300 * 1E-6;
  IMU_errors->M_g[5] = -150 * 1E-6;
  IMU_errors->M_g[6] = 0 * 1E-6;
  IMU_errors->M_a[7] = 0 * 1E-6;
  IMU_errors->M_g[8] = -350 * 1E-6;

  /* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
  IMU_errors->G_g[0] = 0.9 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[1] = -1.1 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[2] = -0.6 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[3] = -0.5 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[4] = 1.9 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[5] = -1.6 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[6] = 0.3 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[7] = 1.1 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[8] = -1.3 * deg_to_rad / (3600 * 9.80665);

  /* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) */
  IMU_errors->accel_noise_root_PSD = 100 * micro_g_to_meters_per_second_squared;
  /* Gyro noise root PSD (deg per root hour, converted to rad s^-0.5) */
  IMU_errors->gyro_noise_root_PSD = 0.01 * deg_to_rad / 60.0;
  /* Accelerometer quantization level (m/s^2) */
  IMU_errors->accel_quant_level = 1E-2;
  /* Gyro quantization level (rad/s) */
  IMU_errors->gyro_quant_level = 2E-4;

  /* Interval between GNSS epochs (s) */
  GNSS_config->epoch_interval = 0.5;

  /* Initial estimated position (m; ECEF) */
  GNSS_config->init_est_r_ea_e[0] = 0;
  GNSS_config->init_est_r_ea_e[1] = 0;
  GNSS_config->init_est_r_ea_e[2] = 0;

  /* Number of satellites in constellation */
  GNSS_config->no_sat = 30;
  /* Orbital radius of satellites (m) */
  GNSS_config->r_os = 2.656175E7;
  /* Inclination angle of satellites (deg) */
  GNSS_config->inclination = 55;
  /* Longitude offset of constellation (deg) */
  GNSS_config->const_delta_lambda = 0;
  /* Timing offset of constellation (s) */
  GNSS_config->const_delta_t = 0;

  /* Mask angle (deg) */
  GNSS_config->mask_angle = 10;
  /* Signal in space error SD (m) *Give residual where corrections are applied */
  GNSS_config->SIS_err_SD = 1;
  /* Zenith ionosphere error SD (m) *Give residual where corrections are applied */
  GNSS_config->zenith_iono_err_SD = 2;
  /* Zenith troposphere error SD (m) *Give residual where corrections are applied*/
  GNSS_config->zenith_trop_err_SD = 0.2;
  /* Code tracking error SD (m) *Can extend to account for multipath */
  GNSS_config->code_track_err_SD = 1;
  /* Range rate tracking error SD (m/s) *Can extend to account for multipath */
  GNSS_config->rate_track_err_SD = 0.02;
  /* Receiver clock offset at time=0 (m); */
  GNSS_config->rx_clock_offset = 10000;
  /* Receiver clock drift at time=0 (m/s); */
  GNSS_config->rx_clock_drift = 100;

  /* Initial attitude uncertainty per axis (deg, converted to rad) */
  for (i = 0; i < 3; i++)
    LC_KF_config->init_att_unc[i] = D2R * 1;
  /* Initial velocity uncertainty per axis (m/s) */
  for (i = 0; i < 3; i++)
    LC_KF_config->init_vel_unc[i] = 0.1;
  /* Initial position uncertainty per axis (m) */
  for (i = 0; i < 3; i++)
    LC_KF_config->init_pos_unc[i] = 10;
  /* Initial accelerometer bias uncertainty per instrument (micro-g, converted
/* to m/s^2) */
  LC_KF_config->init_b_a_unc = 1000 * micro_g_to_meters_per_second_squared;
  /* Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)*/
  LC_KF_config->init_b_g_unc = 10 * deg_to_rad / 3600;

  /* Gyro noise PSD (deg^2 per hour, converted to rad^2/s) */
  LC_KF_config->gyro_noise_PSD = pow((0.02 * deg_to_rad / 60), 2);
  /* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
  LC_KF_config->accel_noise_PSD = pow((200 * micro_g_to_meters_per_second_squared), 2);
  /* Accelerometer bias random walk PSD (m^2 s^-5) */
  LC_KF_config->accel_bias_PSD = 1.0E-7;
  /* Gyro bias random walk PSD (rad^2 s^-3) */
  LC_KF_config->gyro_bias_PSD = 2.0E-12;

  /* Position measurement noise SD per axis (m)  */
  LC_KF_config->pos_meas_SD = 2.5;
  /* Velocity measurement noise SD per axis (m/s)  */
  LC_KF_config->vel_meas_SD = 0.1;
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Initialize_INS_GNSS_TCKF_tactical_grade_IMU(initialization_errors *initialization_errors,
                                                 IMU_errors *IMU_errors, GNSS_config *GNSS_config, TC_KF_config *TC_KF_config)
{
  int i;
  /*%INS_GNSS_Demo_9
%SCRIPT Tightly coupled INS/GNSS demo:
%   Profile_2 (175s car)
%   Tactical-grade IMU
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% Created 27/5/12 by Paul Groves

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details   */

  /* Constants  */
  double deg_to_rad = 0.01745329252;
  double rad_to_deg = 1 / deg_to_rad;
  double micro_g_to_meters_per_second_squared = 9.80665E-6;

  /* CONFIGURATION */
  /* Input truth motion profile filename
input_profile_name = 'Profile_2.csv'; */
  /* Output motion profile and error filenames
output_profile_name = 'INS_GNSS_Demo_9_Profile.csv';
output_errors_name = 'INS_GNSS_Demo_9_Errors.csv';  */

  /* Attitude initialization error (deg, converted to rad; @N,E,D) */
  initialization_errors->delta_eul_nb_n[0] = -0.05 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[1] = 0.04 * deg_to_rad;
  initialization_errors->delta_eul_nb_n[2] = 1 * deg_to_rad; /* rad */

  /* Accelerometer biases (micro-g, converted to m/s^2; body axes) */
  IMU_errors->b_a[0] = 900 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[1] = -1300 * micro_g_to_meters_per_second_squared;
  IMU_errors->b_a[2] = 800 * micro_g_to_meters_per_second_squared;
  /* Gyro biases (deg/hour, converted to rad/sec; body axes) */
  IMU_errors->b_g[0] = -9 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = 13 * deg_to_rad / 3600;
  IMU_errors->b_g[0] = -8 * deg_to_rad / 3600;
  /* Accelerometer scale factor and cross coupling errors (ppm, converted to
 unitless; body axes) */
  IMU_errors->M_a[0] = 500 * 1E-6;
  IMU_errors->M_a[1] = -300 * 1E-6;
  IMU_errors->M_a[2] = 200 * 1E-6;
  IMU_errors->M_a[3] = -150 * 1E-6;
  IMU_errors->M_a[4] = -600 * 1E-6;
  IMU_errors->M_a[5] = 250 * 1E-6;
  IMU_errors->M_a[6] = -250 * 1E-6;
  IMU_errors->M_a[7] = 100 * 1E-6;
  IMU_errors->M_a[8] = 450 * 1E-6;
  /* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
  IMU_errors->M_g[0] = 400 * 1E-6;
  IMU_errors->M_a[1] = -300 * 1E-6;
  IMU_errors->M_g[2] = 250 * 1E-6;
  IMU_errors->M_g[3] = 0 * 1E-6;
  IMU_errors->M_a[4] = -300 * 1E-6;
  IMU_errors->M_g[5] = -150 * 1E-6;
  IMU_errors->M_g[6] = 0 * 1E-6;
  IMU_errors->M_a[7] = 0 * 1E-6;
  IMU_errors->M_g[8] = -350 * 1E-6;

  /* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
  IMU_errors->G_g[0] = 0.9 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[1] = -1.1 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[2] = -0.6 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[3] = -0.5 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[4] = 1.9 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[5] = -1.6 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[6] = 0.3 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[7] = 1.1 * deg_to_rad / (3600 * 9.80665);
  IMU_errors->G_g[8] = -1.3 * deg_to_rad / (3600 * 9.80665);

  /* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) */
  IMU_errors->accel_noise_root_PSD = 100 * micro_g_to_meters_per_second_squared;
  /* Gyro noise root PSD (deg per root hour, converted to rad s^-0.5) */
  IMU_errors->gyro_noise_root_PSD = 0.01 * deg_to_rad / 60.0;
  /* Accelerometer quantization level (m/s^2) */
  IMU_errors->accel_quant_level = 1E-2;
  /* Gyro quantization level (rad/s) */
  IMU_errors->gyro_quant_level = 2E-4;

  /* Interval between GNSS epochs (s) */
  GNSS_config->epoch_interval = 0.5;

  /* Initial estimated position (m; ECEF) */
  GNSS_config->init_est_r_ea_e[0] = 0;
  GNSS_config->init_est_r_ea_e[1] = 0;
  GNSS_config->init_est_r_ea_e[2] = 0;

  /* Number of satellites in constellation */
  GNSS_config->no_sat = 30;
  /* Orbital radius of satellites (m) */
  GNSS_config->r_os = 2.656175E7;
  /* Inclination angle of satellites (deg) */
  GNSS_config->inclination = 55;
  /* Longitude offset of constellation (deg) */
  GNSS_config->const_delta_lambda = 0;
  /* Timing offset of constellation (s) */
  GNSS_config->const_delta_t = 0;

  /* Mask angle (deg) */
  GNSS_config->mask_angle = 10;
  /* Signal in space error SD (m) *Give residual where corrections are applied */
  GNSS_config->SIS_err_SD = 1;
  /* Zenith ionosphere error SD (m) *Give residual where corrections are applied */
  GNSS_config->zenith_iono_err_SD = 2;
  /* Zenith troposphere error SD (m) *Give residual where corrections are applied*/
  GNSS_config->zenith_trop_err_SD = 0.2;
  /* Code tracking error SD (m) *Can extend to account for multipath */
  GNSS_config->code_track_err_SD = 1;
  /* Range rate tracking error SD (m/s) *Can extend to account for multipath */
  GNSS_config->rate_track_err_SD = 0.02;
  /* Receiver clock offset at time=0 (m); */
  GNSS_config->rx_clock_offset = 10000;
  /* Receiver clock drift at time=0 (m/s); */
  GNSS_config->rx_clock_drift = 100;

  /* Initial attitude uncertainty per axis (deg, converted to rad) */
  for (i = 0; i < 3; i++)
    TC_KF_config->init_att_unc[i] = D2R * 1;
  /* Initial velocity uncertainty per axis (m/s) */
  for (i = 0; i < 3; i++)
    TC_KF_config->init_vel_unc[i] = 0.1;
  /* Initial position uncertainty per axis (m) */
  for (i = 0; i < 3; i++)
    TC_KF_config->init_pos_unc[i] = 10;
  /* Initial accelerometer bias uncertainty per instrument (micro-g, converted
/* to m/s^2) */
  TC_KF_config->init_b_a_unc = 1000 * micro_g_to_meters_per_second_squared;
  /* Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)*/
  TC_KF_config->init_b_g_unc = 10 * deg_to_rad / 3600;
  /* Initial clock offset uncertainty per axis (m) */
  TC_KF_config->init_clock_offset_unc = 10;
  /* Initial clock drift uncertainty per axis (m/s) */
  TC_KF_config->init_clock_drift_unc = 0.1;

  /* Gyro noise PSD (deg^2 per hour, converted to rad^2/s) */
  TC_KF_config->gyro_noise_PSD = pow((0.02 * deg_to_rad / 60), 2);
  /* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
  TC_KF_config->accel_noise_PSD = pow((200 * micro_g_to_meters_per_second_squared), 2);
  /* Accelerometer bias random walk PSD (m^2 s^-5) */
  TC_KF_config->accel_bias_PSD = 1.0E-7;
  /* Gyro bias random walk PSD (rad^2 s^-3) */
  TC_KF_config->gyro_bias_PSD = 2.0E-12;
  /* Receiver clock frequency-drift PSD (m^2/s^3) */
  TC_KF_config->clock_freq_PSD = 1;
  /* Receiver clock phase-drift PSD (m^2/s) */
  TC_KF_config->clock_phase_PSD = 1;

  /* Pseudo-range measurement noise SD (m) */
  TC_KF_config->pseudo_range_SD = 2.5;
  /* Pseudo-range rate measurement noise SD (m/s) */
  TC_KF_config->range_rate_SD = 0.1;
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
extern void Skew_symmetric(double *vec, double *W)
{
  W[0] = 0.0;
  W[1] = -vec[2];
  W[2] = vec[1];
  W[3] = vec[2];
  W[4] = 0.0;
  W[5] = -vec[0];
  W[6] = -vec[1];
  W[7] = vec[0];
  W[8] = 0.0;
}
/* Name of function ------------------------------------------------------------

%Gravitation_ECI - Calculates  acceleration due to gravity resolved about
%ECEF-frame
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 1/4/2012 by Paul Groves
%
% Inputs:
%   r_eb_e  Cartesian position of body frame w.r.t. ECEF frame, resolved
%           about ECEF-frame axes (m)
% Outputs:
%   g       Acceleration due to gravity (m/s^2)

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details
------------------------------------------------------------------------------*/
void Gravity_ECEF(double *r_eb_e, double *g)
{
  double mag_r, z_scale, gamma[3];

  /* Parameters  */
  double R_0 = RE_GRS80;  /* WGS84 Equatorial radius in meters */
  double omega_ie = OMGE; /* Earth rotation rate (rad/s)  */

  /* Begins  */

  /* Calculate distance from center of the Earth  */
  mag_r = norm(r_eb_e, 3);

  /* If the input position is 0,0,0, produce a dummy output */
  if (mag_r == 0)
  {
    g[0] = 0;
    g[1] = 0;
    g[2] = 0;

    /* Calculate gravitational acceleration using (2.142)  */
  }
  else
  {
    z_scale = 5 * pow((r_eb_e[2] / mag_r), 2);
    gamma[0] = -(mu / pow(mag_r, 3)) * (r_eb_e[0] + (1.5 * J_2 * pow((R_0 / mag_r), 2)) *
                                                        ((1 - z_scale) * r_eb_e[0]));
    gamma[1] = -(mu / pow(mag_r, 3)) * (r_eb_e[1] + (1.5 * J_2 * pow((R_0 / mag_r), 2)) *
                                                        ((1 - z_scale) * r_eb_e[1]));
    gamma[2] = -(mu / pow(mag_r, 3)) * (r_eb_e[2] + (1.5 * J_2 * pow((R_0 / mag_r), 2)) *
                                                        ((3 - z_scale) * r_eb_e[2]));

    /* Add centripetal acceleration using (2.133)  */
    g[0] = gamma[0] + omega_ie * omega_ie * r_eb_e[0];
    g[1] = gamma[1] + omega_ie * omega_ie * r_eb_e[1];
    g[2] = gamma[2];

  } //end % if
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Euler_to_CTM(double *delta_eul_nb_n, double *C)
{
  double sin_phi, cos_phi, sin_theta, cos_theta, sin_psi, cos_psi;

  /* Precalculate sines and cosines of the Euler angles */
  sin_phi = sin(delta_eul_nb_n[0]);
  cos_phi = cos(delta_eul_nb_n[0]);
  sin_theta = sin(delta_eul_nb_n[1]);
  cos_theta = cos(delta_eul_nb_n[1]);
  sin_psi = sin(delta_eul_nb_n[2]);
  cos_psi = cos(delta_eul_nb_n[2]);

  /* Calculate coordinate transformation matrix using (2.22) */
  C[0] = cos_theta * cos_psi;
  C[1] = cos_theta * sin_psi;
  C[2] = -sin_theta;
  C[3] = (-cos_phi * sin_psi) + (sin_phi * sin_theta * cos_psi);
  C[4] = (cos_phi * cos_psi) + (sin_phi * sin_theta * sin_psi);
  C[5] = sin_phi * cos_theta;
  C[6] = (sin_phi * sin_psi) + (cos_phi * sin_theta * cos_psi);
  C[7] = (-sin_phi * cos_psi) + (cos_phi * sin_theta * sin_psi);
  C[8] = cos_phi * cos_theta;
}

/* Compute Euler from Cbn ------------------------------------------------------
* description: Compute euler angles from the DCM Cbn
* args   : double *euler  vector with {roll, pitch, yaw} angles    IO
* ref.: (Shin, 2001, page 12) and Groves(prev version of 2013, page 28)
*-----------------------------------------------------------------------------*/
void CTM_to_Euler(double *euler, double *Cbn)
{
  /* Roll (x rotation) - phi
  euler[0]=atan2(Cbn[7],Cbn[8]);
  /* Pitch (y rotation) - theta
//  imu->aea[1]=-1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6])));
  euler[1]=-asin(Cbn[6]);
  /* Yaw (z rotation) - psi
  euler[2]=atan2(Cbn[3],Cbn[0]);*/
  /*
  printf("Euler 0: %lf\n", atan2(Cbn[7],Cbn[8]));
  printf("Euler 1: %lf or %lf\n", -asin(Cbn[6]), -1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6]))) );
  printf("Euler 2: %lf\n",  atan2(Cbn[3],Cbn[0]));*/

  /* Grove for attitude -> Cbn 2.24

 euler[0] = atan2(Cbn[7],Cbn[8]); /* Roll (x rotation) - phi */
  //imu->aea[1]=-1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6])));
  /*euler[1]= -asin(Cbn[6]);   Pitch (y rotation) - theta */

  /* euler[2]= atan2(Cbn[3],Cbn[0]);   Yaw (z rotation) - psi */
  /*
 printf("Euler 0: %lf\n", euler[0]);
 printf("Euler 1: %lf\n", euler[1]);
 printf("Euler 2: %lf\n",  euler[2]);  */

  /* From Groves Matlab code */

  /* Calculate Euler angles using 2.23 */

  euler[0] = atan2(Cbn[5], Cbn[8]); /* Roll (x rotation) - phi */
                                    //imu->aea[1]=-1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6])));
  euler[1] = -asin(Cbn[2]);         /* Pitch (y rotation) - theta */

  euler[2] = atan2(Cbn[1], Cbn[0]); /* Yaw (z rotation) - psi */
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void ECEF_to_NED(double *r_eb_e, double *v_eb_e, double *C_b_e,
                 double *L_b, double *lambda_b, double *h_b, double *v_eb_n,
                 double *C_b_n)
{
  /*%ECEF_to_NED - Converts Cartesian  to curvilinear position, velocity
%resolving axes from ECEF to NED and attitude from ECEF- to NED-referenced
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 2/4/2012 by Paul Groves
%
% Inputs:
%   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
%                 along ECEF-frame axes (m)
%   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
%                 ECEF-frame axes (m/s)
%   C_b_e         body-to-ECEF-frame coordinate transformation matrix
%
% Outputs:
%   L_b           latitude (rad)
%   lambda_b      longitude (rad)
%   h_b           height (m)
%   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
%                 north, east, and down (m/s)
%   C_b_n         body-to-NED coordinate transformation matrix   */

  /* Parameters
   R_0 = 6378137; %WGS84 Equatorial radius in meters
   e = 0.0818191908425; %WGS84 eccentricity   */

  /* Begins  */
  double k1, k2, beta, E, F, P, Q, D, V, G, T;
  double cos_lat, sin_lat, cos_long, sin_long, C_e_n[9];
  double pos[3];

  /* Convert position using Borkowski closed-form exact solution  */
  /* From (2.113)  */

  *lambda_b = atan2(r_eb_e[1], r_eb_e[0]);

  /* From (C.29) and (C.30)  */
  k1 = sqrt(1 - e_2) * fabs(r_eb_e[2]);

  k2 = e_2 * RE_GRS80;
  beta = sqrt((r_eb_e[0] * r_eb_e[0]) + (r_eb_e[1] * r_eb_e[1]));
  E = (k1 - k2) / beta;
  F = (k1 + k2) / beta;

  /* From (C.31)  */
  P = 4 / 3 * ((E * F) + 1);

  /* From (C.32)  */
  Q = 2 * ((E * E) - (F * F));

  /* From (C.33)  */
  D = (P * P * P) + (Q * Q);

  /* From (C.34)  */
  V = pow((sqrt(D) - Q), (1 / 3)) - pow((sqrt(D) + Q), (1 / 3));

  /* From (C.35) */
  G = 0.5 * (sqrt((E * E) + V) + E);

  /* From (C.36) */
  T = sqrt((G * G) + (F - (V * G)) / ((2 * G) - E)) - G;

  /* From (C.37)  */
  *L_b = sign(r_eb_e[2]) * atan((1 - (T * T)) / (2 * T * sqrt(1 - e_2)));

  /* From (C.38) */
  *h_b = (beta - (RE_GRS80 * T)) * cos(*L_b) +
         (r_eb_e[2] - (sign(r_eb_e[2]) * RE_GRS80 * sqrt(1 - e_2))) * sin(*L_b);

  /* using RTKLIB because above computations are not consistent */
  ecef2pos(r_eb_e, pos);
  *L_b = pos[0];
  *lambda_b = pos[1];
  *h_b = pos[2];

  /* Calculate ECEF to NED coordinate transformation matrix using (2.150) */
  cos_lat = cos(*L_b);
  sin_lat = sin(*L_b);
  cos_long = cos(*lambda_b);
  sin_long = sin(*lambda_b);

  C_e_n[0] = -sin_lat * cos_long;
  C_e_n[1] = -sin_lat * sin_long;
  C_e_n[2] = cos_lat;
  C_e_n[3] = -sin_long;
  C_e_n[4] = cos_long;
  C_e_n[5] = 0;
  C_e_n[6] = -cos_lat * cos_long;
  C_e_n[7] = -cos_lat * sin_long;
  C_e_n[8] = -sin_lat;

  /* Transform velocity using (2.73) */
  matmul_row("NN", 3, 1, 3, 1.0, C_e_n, v_eb_e, 0.0, v_eb_n);

  /* Transform attitude using (2.15) */
  matmul_row("NN", 3, 3, 3, 1.0, C_e_n, C_b_e, 0.0, C_b_n);
}

void NED_to_ECEF(double *L_b, double *lambda_b, double *h_b, double *v_eb_n,
                 double *C_b_n, double *r_eb_e, double *v_eb_e, double *C_b_e)
{
  /*NED_to_ECEF - Converts curvilinear to Cartesian position, velocity
%resolving axes from NED to ECEF ECEF_to_NEDand attitude from NED- to ECEF-referenced
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 2/4/2012 by Paul Groves
%
% Inputs:
%   L_b           latitude (rad)
%   lambda_b      longitude (rad)
%   h_b           height (m)
%   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
%                 north, east, and down (m/s)
%   C_b_n         body-to-NED coordinate transformation matrix
%
% Outputs:
%   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
%                 along ECEF-frame axes (m)
%   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
%                 ECEF-frame axes (m/s)
%   C_b_e         body-to-ECEF-frame coordinate transformation matrix

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details

% Parameters
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity  */

  /* Begins  */
  double cos_lat, sin_lat, cos_long, sin_long, C_e_n[9], C_n_e[9], R_E;

  /* Calculate transverse radius of curvature using (2.105)  */
  R_E = RN(*L_b);

  /* Convert position using (2.112)  */
  cos_lat = cos(*L_b);
  sin_lat = sin(*L_b);
  cos_long = cos(*lambda_b);
  sin_long = sin(*lambda_b);
  r_eb_e[0] = (R_E + *h_b) * cos_lat * cos_long;
  r_eb_e[1] = (R_E + *h_b) * cos_lat * sin_long;
  r_eb_e[2] = ((1 - e_2) * R_E + *h_b) * sin_lat;

  /* Calculate ECEF to NED coordinate transformation matrix using (2.150) */
  C_e_n[0] = -sin_lat * cos_long;
  C_e_n[1] = -sin_lat * sin_long;
  C_e_n[2] = cos_lat;
  C_e_n[3] = -sin_long;
  C_e_n[4] = cos_long;
  C_e_n[5] = 0;
  C_e_n[6] = -cos_lat * cos_long;
  C_e_n[7] = -cos_lat * sin_long;
  C_e_n[8] = -sin_lat;

  /*  Calculate ECEF to NED coordinate transformation matrix using (2.150)  */
  C_n_e[0] = -sin_lat * cos_long;
  C_n_e[1] = -sin_long;
  C_n_e[2] = -cos_lat * cos_long;
  C_n_e[3] = -sin_lat * sin_long;
  C_n_e[4] = cos_long;
  C_n_e[5] = -cos_lat * sin_long;
  C_n_e[6] = cos_lat;
  C_n_e[7] = 0.0;
  C_n_e[8] = -sin_lat;

  /* Transform velocity using (2.73)  */
  matmul_row("NN", 3, 1, 3, 1.0, C_n_e, v_eb_n, 0.0, v_eb_e);

  /* Transform attitude using (2.15) */
  matmul_row("NN", 3, 3, 3, 1.0, C_n_e, C_b_n, 0.0, C_b_e);
}

/* function declaration ------------------------------------------------------*/

void pv_ECEF_to_NED(double *r_eb_e, double *v_eb_e, double *L_b,
                    double *lambda_b, double *h_b, double *v_eb_n)
{
  /*pv_ECEF_to_NED - Converts Cartesian  to curvilinear position and velocity
   %resolving axes from ECEF to NED
   %
   % Software for use with "Principles of GNSS, Inertial, and Multisensor
   % Integrated Navigation Systems," Second Edition.
   %
   % This function created 11/4/2012 by Paul Groves
   %
   % Inputs:
   %   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
   %                 along ECEF-frame axes (m)
   %   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
   %                 ECEF-frame axes (m/s)
   %
   % Outputs:
   %   L_b           latitude (rad)
   %   lambda_b      longitude (rad)
   %   h_b           height (m)
   %   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
   %                 north, east, and down (m/s)

   % Copyright 2012, Paul Groves
   % License: BSD; see license.txt for details

   % Parameters
   R_0 = 6378137; %WGS84 Equatorial radius in meters
   e = 0.0818191908425; %WGS84 eccentricity   */

  /* Begins  */
  double k1, k2, beta, E, F, P, Q, D, V, G, T;
  double cos_lat, sin_lat, cos_long, sin_long, C_e_n[9];
  double pos[3];

  /* Convert position using Borkowski closed-form exact solution  */
  /* From (2.113)  */
  *lambda_b = atan2(r_eb_e[1], r_eb_e[0]);

  /* From (C.29) and (C.30)  */
  k1 = sqrt(1 - e_2) * fabs(r_eb_e[2]);
  k2 = e_2 * RE_GRS80;
  beta = sqrt((r_eb_e[0] * r_eb_e[0]) + (r_eb_e[1] * r_eb_e[1]));
  E = (k1 - k2) / beta;
  F = (k1 + k2) / beta;

  /* From (C.31)  */
  P = 4 / 3 * (E * F + 1);

  /* From (C.32)  */
  Q = 2 * (E * E - F * F);

  /* From (C.33)  */
  D = P * P * P + Q * Q;

  /* From (C.34)  */
  V = pow((sqrt(D) - Q), (1 / 3)) - pow((sqrt(D) + Q), (1 / 3));

  /* From (C.35) */
  G = 0.5 * (sqrt(E * E + V) + E);

  /* From (C.36) */
  T = sqrt(G * G + (F - V * G) / (2 * G - E)) - G;

  /* From (C.37)  */
  *L_b = sign(r_eb_e[2]) * atan((1 - T * T) / (2 * T * sqrt(1 - e_2)));

  /* From (C.38) */
  *h_b = (beta - RE_GRS80 * T) * cos(*L_b) +
         (r_eb_e[2] - sign(r_eb_e[2]) * RE_GRS80 * sqrt(1 - e_2)) * sin(*L_b);

  /* Using RTKLIB because computations are not consistent  */
  ecef2pos(r_eb_e, pos);
  *L_b = pos[0];
  *lambda_b = pos[1];
  *h_b = pos[2];

  /* Calculate ECEF to NED coordinate transformation matrix using (2.150) */
  cos_lat = cos(*L_b);
  sin_lat = sin(*L_b);
  cos_long = cos(*lambda_b);
  sin_long = sin(*lambda_b);
  C_e_n[0] = -sin_lat * cos_long;
  C_e_n[1] = -sin_lat * sin_long;
  C_e_n[2] = cos_lat;
  C_e_n[3] = -sin_long;
  C_e_n[4] = cos_long;
  C_e_n[5] = 0;
  C_e_n[6] = -cos_lat * cos_long;
  C_e_n[7] = -cos_lat * sin_long;
  C_e_n[8] = -sin_lat;

  /* Transform velocity using (2.73) */
  matmul_row("NN", 3, 1, 3, 1.0, C_e_n, v_eb_e, 0.0, v_eb_n);
  //matmul("NN",1,3,3,1.0, v_eb_e, C_e_n, 0.0, v_eb_n);

  /* Ends */
}

void pv_NED_to_ECEF(double *L_b, double *lambda_b, double *h_b, double *v_eb_n,
                    double *r_eb_e, double *v_eb_e)
{
  /*%pv_NED_to_ECEF - Converts curvilinear to Cartesian position and velocity
%resolving axes from NED to ECEF
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 11/4/2012 by Paul Groves
%
% Inputs:
%   L_b           latitude (rad)
%   lambda_b      longitude (rad)
%   h_b           height (m)
%   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
%                 north, east, and down (m/s)
%
% Outputs:
%   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
%                 along ECEF-frame axes (m)
%   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
%                 ECEF-frame axes (m/s)

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details

% Parameters
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity  */

  /* Begins  */
  double cos_lat, sin_lat, cos_long, sin_long, C_e_n[9], R_E;

  /* Calculate transverse radius of curvature using (2.105)  */
  R_E = RN(*L_b);

  /* Convert position using (2.112)  */
  cos_lat = cos(*L_b);
  sin_lat = sin(*L_b);
  cos_long = cos(*lambda_b);
  sin_long = sin(*lambda_b);
  r_eb_e[0] = (R_E + *h_b) * cos_lat * cos_long;
  r_eb_e[1] = (R_E + *h_b) * cos_lat * sin_long;
  r_eb_e[2] = ((1 - e_2) * R_E + *h_b) * sin_lat;

  /* Calculate ECEF to NED coordinate transformation matrix using (2.150) */
  C_e_n[0] = -sin_lat * cos_long;
  C_e_n[1] = -sin_lat * sin_long;
  C_e_n[2] = cos_lat;
  C_e_n[3] = -sin_long;
  C_e_n[4] = cos_long;
  C_e_n[5] = 0;
  C_e_n[6] = -cos_lat * cos_long;
  C_e_n[7] = -cos_lat * sin_long;
  C_e_n[8] = -sin_lat;

  /* Transform velocity using (2.73)  */
  matmul("TN", 3, 1, 3, 1.0, C_e_n, v_eb_n, 0.0, v_eb_e);

  /* Ends */
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Calculate_errors_NED(double *est_L_b, double *est_lambda_b,
                          double *est_h_b, double *est_v_eb_n, double *est_C_b_n, double *true_L_b,
                          double *true_lambda_b, double *true_h_b, double *true_v_eb_n,
                          double *true_C_b_n, double *delta_r_eb_n, double *delta_v_eb_n,
                          double *delta_eul_nb_n)
{
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Initialize_NED_attitude(double *C_b_n,
                             initialization_errors *initialization_errors, double *est_C_b_n)
{
  double delta_C_b_n[9];
  int i;

  /* Attitude initialization, using (5.109) and (5.111) */
  for (i = 0; i < 3; i++)
    initialization_errors->delta_eul_nb_n[i] =
        -initialization_errors->delta_eul_nb_n[i];

  Euler_to_CTM(initialization_errors->delta_eul_nb_n, delta_C_b_n);

  /* est_C_b_n = delta_C_b_n * C_b_n; */
  matmul("NN", 3, 3, 3, 1.0, delta_C_b_n, C_b_n, 0.0, est_C_b_n);
}
/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Initialize_TC_P_matrix(TC_KF_config *TC_KF_config, double *P_matrix)
{
  int i, j, npar = 17;
  double att_var[3], vel_var[3], pos_var[3], ba_var, bg_var, dt_off_var, dt_drift_var;

  /* Initialize error covariance matrix */
  // IT ASSUMES THE SAME VALUE FOR EACH AXES, HOWEVER IT MAY DIFFER!!!!!!!!!
  for (i = 0; i < 3; i++)
  {
    att_var[i] = TC_KF_config->init_att_unc[i] * TC_KF_config->init_att_unc[i];
    vel_var[i] = TC_KF_config->init_vel_unc[i] * TC_KF_config->init_vel_unc[i];
    pos_var[i] = TC_KF_config->init_pos_unc[i] * TC_KF_config->init_pos_unc[i];
  }

  ba_var = TC_KF_config->init_b_a_unc * TC_KF_config->init_b_a_unc;
  bg_var = TC_KF_config->init_b_g_unc * TC_KF_config->init_b_g_unc;
  dt_off_var = TC_KF_config->init_clock_offset_unc * TC_KF_config->init_clock_offset_unc;
  dt_drift_var = TC_KF_config->init_clock_drift_unc * TC_KF_config->init_clock_drift_unc;

  /* The 17-states {de,dv,dr,dba,dbg, dt, dtdot} system model order */

  for (i = 0; i < 3; i++)
  {
    printf("att: %lf, vel: %lf, pos: %lf\n", att_var[i], vel_var[i], pos_var[i]);
  }

  /* Attitude  */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      (i == j ? P_matrix[i * npar + j] = att_var[i] : 0.0);
    }
  }

  /* Velocity */
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      (i == j ? P_matrix[i * npar + j] = vel_var[j - 3] : 0.0);
    }
  }

  /* Position */
  for (i = 6; i < 9; i++)
  {
    for (j = 6; j < 9; j++)
    {
      (i == j ? P_matrix[i * npar + j] = pos_var[j - 6] : 0.0);
    }
  }

  /* Acc. bias  */
  for (i = 9; i < 12; i++)
  {
    for (j = 9; j < 12; j++)
    {
      (i == j ? P_matrix[i * npar + j] = ba_var : 0.0);
    }
  }

  /* Gyro. bias  */
  for (i = 12; i < 15; i++)
  {
    for (j = 12; j < 15; j++)
    {
      (i == j ? P_matrix[i * npar + j] = bg_var : 0.0);
    }
  }

  /* Clock offset and drift variances */
  P_matrix[15 * npar + 15] = dt_off_var;
  P_matrix[16 * npar + 16] = dt_drift_var;
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
void Initialize_LC_P_matrix(LC_KF_config *LC_KF_config, double *P_matrix)
{
  int i, j, npar = 15;
  double att_var[3], vel_var[3], pos_var[3], ba_var, bg_var;

  /* Initialize error covariance matrix */
  // IT ASSUMES THE SAME VALUE FOR EACH AXES, HOWEVER IT MAY DIFFER!!!!!!!!!
  for (i = 0; i < 3; i++)
  {
    att_var[i] = LC_KF_config->init_att_unc[i] * LC_KF_config->init_att_unc[i];
    vel_var[i] = LC_KF_config->init_vel_unc[i] * LC_KF_config->init_vel_unc[i];
    pos_var[i] = LC_KF_config->init_pos_unc[i] * LC_KF_config->init_pos_unc[i];
  }

  ba_var = LC_KF_config->init_b_a_unc * LC_KF_config->init_b_a_unc;
  bg_var = LC_KF_config->init_b_g_unc * LC_KF_config->init_b_g_unc;

  /* The 15-states {de,dv,dr,dba,dbg} system model order */

  /* Attitude  */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      (i == j ? P_matrix[i * npar + j] = att_var[j] : 0.0);
    }
  }

  /* Velocity */
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      (i == j ? P_matrix[i * npar + j] = vel_var[j - 3] : 0.0);
    }
  }

  /* Position */
  for (i = 6; i < 9; i++)
  {
    for (j = 6; j < 9; j++)
    {
      (i == j ? P_matrix[i * npar + j] = pos_var[j - 6] : 0.0);
    }
  }

  /* Acc. bias  */
  for (i = 9; i < 12; i++)
  {
    for (j = 9; j < 12; j++)
    {
      (i == j ? P_matrix[i * npar + j] = ba_var : 0.0);
    }
  }

  /* Gyro. bias  */
  for (i = 12; i < 15; i++)
  {
    for (j = 12; j < 15; j++)
    {
      (i == j ? P_matrix[i * npar + j] = bg_var : 0.0);
    }
  }
}

/* Quaternion multiplication ---------------------------------------------------
* description: Multiplies p . q  quaternions
* args   : double *p      I    first quaternion {4x1}
*	         double *q		  I    second quaternion {4x1}
*	         double *pq     IO   resulting quaternion {4x1}
* Groves (2103) from Appendix E - at E.6.3 section, apge E-10
*-----------------------------------------------------------------------------*/
quaternion_mult(double *p, double *q, double *pq)
{
  pq[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  pq[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  pq[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  pq[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

/* Attitude update -------------------------------------------------------
* description: Update the quaternion from gyro corrected measurements
* args   : um7pack_t* imu       IO    current (tk) imu structure
*	   double dt		I     time difference from previous measu. (tk-1)
*	   double* wien[3], wenn[3]	I	rotation vectors
*	   double* C[9]		IO    receives an initial Cbn and outputs the updated one
Reference: Initialy from Shin (2001, pag.22) adapting to Grove (2103) from
 Appendix E - at E.6.3 section
*-----------------------------------------------------------------------------*/
void attitude_update(double *alpha_ib_b, double *omega_ie, float t, double *old_C_b_e, double *new_q_b_e)
{
  double old_q_b_e[4], q_less_plus[4], q_omega[4], q_aux[4], q_aux1[4];
  double mag_alpha, ac = 0.0, as = 0.0;
  int i;

  /* Form quaternions */

  /* Old Cbe quaternion */
  DCM_to_quaternion(old_C_b_e, old_q_b_e);

  /* Reverse quaternion to become old_q_e_b */
  old_q_b_e[0] = -old_q_b_e[0];

  /* Coeficients of q_less_plus quaternion */
  mag_alpha = norm(alpha_ib_b, 3);
  ac = 1 - 0.5 * (mag_alpha * mag_alpha / 4) + ((1 / 24) * (pow(mag_alpha / 2, 4)));
  as = 0.5 - ((1 / 12) * (pow(mag_alpha / 2, 2)));

  /* q_less_plus quaternion */
  q_less_plus[0] = ac;
  for (i = 1; i < 4; i++)
    q_less_plus[i] = as * alpha_ib_b[i - 1];

  /*q_omega quaternion */
  q_omega[0] = 0.0;
  for (i = 1; i < 4; i++)
    q_omega[i] = 0.5 * omega_ie[i - 1] * t;

  quaternion_mult(q_omega, old_q_b_e, q_aux);

  quaternion_mult(old_q_b_e, q_less_plus, q_aux1);

  /* The updated quternion */
  for (i = 0; i < 4; i++)
    new_q_b_e[i] = q_aux1[i] - q_aux[i];

  /* Reverse quaternion */
  new_q_b_e[0] = -new_q_b_e[0];
}

/* Convert rotation matrix to quaternion ----------------- ---------------------
* description: Cbe or Cbeta_alpha to quaternion transformation
* args   : um7pack_t* imu       IO    imu structure
*
Reference: Groves (2013, pag.41)
*-----------------------------------------------------------------------------*/
void DCM_to_quaternion(double *C, double *q)
{
  /* Transformation between DCM Cbn and the quaternion is accomplished by: */
  q[0] = 0.5 * sqrt(1 + C[0] + C[4] + C[8]);

  if (q[0] < 0.001)
  {
    q[1] = 0.5 * sqrt(1 + C[0] + C[4] + C[8]);
    q[0] = 0.25 * (C[5] - C[7]) / (q[1]);
    q[2] = 0.25 * (C[3] - C[1]) / (q[1]);
    q[3] = 0.25 * (C[6] - C[2]) / (q[1]);
  }
  else
  {
    q[1] = 0.25 * (C[5] - C[7]) / (q[0]);
    q[2] = 0.25 * (C[6] - C[2]) / (q[0]);
    q[3] = 0.25 * (C[1] - C[3]) / (q[0]);
  }
}

/* Quaternion to DCM ------------------------------------- ---------------------
* description: transform a quaternion to DCM Cbn amtrix
* args   : double* C[9]       IO    DCM Cbe matrix
* 	       double* q[4]	      I   quaternion
Reference: Groves (2013, pag.41)
*-----------------------------------------------------------------------------*/
void Quaternion_to_DCM(double *q, double *C)
{
  /* The transformation between the quaternion and the DCM Cbn is accomplished by */
  C[0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  C[1] = 2 * (q[1] * q[2] + q[3] * q[0]);
  C[2] = 2 * (q[1] * q[3] - q[2] * q[0]);
  C[3] = 2 * (q[1] * q[2] - q[3] * q[0]);
  C[4] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  C[5] = 2 * (q[2] * q[3] + q[1] * q[0]);
  C[6] = 2 * (q[1] * q[3] + q[2] * q[0]);
  C[7] = 2 * (q[2] * q[3] - q[1] * q[0]);
  C[8] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}

/* Euler to quaternion ----------------------------------- ---------------------
* description: transform an euler vector to quaternion
* args   : double *eul_nb       I    Euler angles vector {3x1}
* 	       double* q_nb[4]	      IO   quaternion {4x1}
Reference: Groves (2013, pag.42)
*-----------------------------------------------------------------------------*/
Euler_to_quaternion(double *eul_nb, double *q_nb)
{
  double phi_half = eul_nb[0] / 2, theta_half = eul_nb[1] / 2,
         yaw_half = eul_nb[2] / 2;

  q_nb[0] = cos(phi_half) * cos(theta_half) * cos(yaw_half) + sin(phi_half) * sin(theta_half) * sin(yaw_half);
  q_nb[1] = sin(phi_half) * cos(theta_half) * cos(yaw_half) - cos(phi_half) * sin(theta_half) * sin(yaw_half);
  q_nb[2] = cos(phi_half) * sin(theta_half) * cos(yaw_half) + sin(phi_half) * cos(theta_half) * sin(yaw_half);
  q_nb[3] = cos(phi_half) * cos(theta_half) * sin(yaw_half) - sin(phi_half) * sin(theta_half) * cos(yaw_half);
}

/* Quaternion to Euler ---------------------------------- ---------------------
* description: transform an quaternion to euler vector
* args   : double *q_nb      I    quaternion {4x1}
* 	       double  *eul_nb      IO   Euler angles vector {3x1}
Reference: Groves (2013, pag.42)
*-----------------------------------------------------------------------------*/
Quaternion_to_euler(double *q, double *eul)
{

  eul[0] = atan2(2 * (q[0] * q[1] + q[2] * q[3]), (1 - 2 * q[1] - 2 * q[2]));
  eul[1] = asin(2 * (q[0] * q[2] - q[1] * q[3]));
  eul[2] = atan2(2 * (q[0] * q[3] + q[1] * q[2]), (1 - 2 * q[2] - 2 * q[3]));

  printf("Quat.euler: %lf\n", 2 * (q[0] * q[2] - q[1] * q[3]));
}

/* Attitude update with quaternion --------------------- ---------------------
* description: quaternion approach for updating the attitude matrix Cbn
* args   : um7pack_t* imu       IO    imu structure
*
Reference: Shin (2001, pag.18)
*-----------------------------------------------------------------------------*/
void Cbn_from_rotations(um7pack_t *imu, double *C)
{
  double q[4], u, qdoti[16], qdot[4];

  /* Building quaternion from a rotation angle vector u{ux,uy,uz} */

  /* QUEST: Is it from the attitude Euler angles or  the rotation rate from gyroscopes measurements */
  u = sqrt(imu->aea[0] * imu->aea[0] + imu->aea[1] * imu->aea[1] + imu->aea[2] * imu->aea[2]);

  q[0] = (imu->g[0] / u) * sin(u / 2);
  q[1] = (imu->g[1] / u) * sin(u / 2);
  q[2] = (imu->g[2] / u) * sin(u / 2);
  q[3] = cos(u / 2);

  /* Normality condition
 if( q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] > 1.000001 || q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] < 0.999999 )
 q[]=q[]/sqrt(qTq); //for each component
*/

  qdoti[0] = qdoti[5] = qdoti[10] = qdoti[15] = 0;
  qdoti[1] = qdoti[11] = imu->g[2];  /*wz*/
  qdoti[4] = qdoti[14] = -imu->g[2]; /*-wz*/
  qdoti[7] = qdoti[8] = imu->g[1];   /*wy*/
  qdoti[13] = qdoti[2] = -imu->g[1]; /*-wy*/
  qdoti[3] = qdoti[6] = imu->g[0];   /*wx*/
  qdoti[9] = qdoti[12] = -imu->g[0]; /*-wx*/

  matmul("NN", 4, 1, 4, 0.5, qdoti, q, 0.0, qdot);

  /* The transformation between the quaternion and the DCM Cbn is accomplished by */
  C[0] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  C[1] = 2 * (q[0] * q[1] - q[2] * q[3]);
  C[2] = 2 * (q[0] * q[2] - q[1] * q[3]);
  C[3] = 2 * (q[0] * q[1] - q[2] * q[3]);
  C[4] = q[1] * q[1] - q[0] * q[0] - q[2] * q[2] + q[3] * q[3];
  C[5] = 2 * (q[1] * q[2] - q[0] * q[3]);
  C[6] = 2 * (q[0] * q[2] - q[1] * q[3]);
  C[7] = 2 * (q[1] * q[2] - q[0] * q[3]);
  C[8] = q[2] * q[2] - q[0] * q[0] - q[1] * q[1] + q[3] * q[3];
}

/* Cbn update from attitude error estimate --------------- ---------------------
* description: Corrects previous Cbe by the estimated attitude error to a new Cbe
* args   : double *est_delta_Psi     I    vector with attitude errors {3x1}
           double *est_C_b_e_old     I    old Cbe transformation matrix {3x3}
           double *est_C_b_e_new     O    new Cbe transformation matrix {3x3}
*
Reference: Grove (2013) Appendix E, section E.6.3 page E-14
*-----------------------------------------------------------------------------*/
void Quaternion_attitude_errror_correction(double *est_delta_Psi,
                                           double *est_C_b_e_old, double *est_C_b_e_new)
{
  double q_new[4], q_old[4], M_delta_Psi[16];

  /* Obtain Cbe old quaternion */
  DCM_to_quaternion(est_C_b_e_old, q_old);

  printf("OLD_QUAT NORM %lf\n", norm(q_old, 4));

  M_delta_Psi[0] = M_delta_Psi[5] = M_delta_Psi[10] = M_delta_Psi[15] = 1.0;
  M_delta_Psi[1] = M_delta_Psi[11] = -est_delta_Psi[0]; /*-wx*/
  M_delta_Psi[4] = M_delta_Psi[14] = est_delta_Psi[0];  /*wx*/
  M_delta_Psi[2] = M_delta_Psi[13] = -est_delta_Psi[1]; /*wy*/
  M_delta_Psi[7] = M_delta_Psi[8] = est_delta_Psi[1];   /*-wy*/
  M_delta_Psi[3] = M_delta_Psi[6] = -est_delta_Psi[2];  /*wz*/
  M_delta_Psi[9] = M_delta_Psi[12] = est_delta_Psi[2];  /*-wz*/

  matmul_row("NN", 4, 1, 4, -0.5, M_delta_Psi, q_old, 0.0, q_new);
  //matmul("NN", 1, 4, 4, -0.5, q_old, M_delta_Psi, 0.0, q_new);

  printf("NEW_QUAT NORM %lf\n", norm(q_new, 4));

  /* Obtain Cbe new from quaternion */
  Quaternion_to_DCM(q_new, est_C_b_e_new);
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
/* TC_KF_Epoch - Implements one cycle of the tightly coupled INS/GNSS
    /* extended Kalman filter plus closed-loop correction of all inertial states
    /*
    /* Software for use with "Principles of GNSS, Inertial, and Multisensor
    /* Integrated Navigation Systems," Second Edition.
    /*
    /* This function created 12/4/2012 by Paul Groves
    /*
    /* Inputs:
    /*   GNSS_measurements     GNSS measurement data:
    /*     Column 1              Pseudo-range measurements (m)
    /*     Column 2              Pseudo-range rate measurements (m/s)
    /*     Columns 3-5           Satellite ECEF position (m)
    /*     Columns 6-8           Satellite ECEF velocity (m/s)
    /*   no_meas               Number of satellites for which measurements are
    /*                         supplied
    /*   tor_s                 propagation interval (s)
    /*   est_C_b_e_old         prior estimated body to ECEF coordinate
    /*                         transformation matrix
    /*   est_v_eb_e_old        prior estimated ECEF user velocity (m/s)
    /*   est_r_eb_e_old        prior estimated ECEF user position (m)
    /*   est_IMU_bias_old      prior estimated IMU biases (body axes)
    /*   est_clock_old         prior Kalman filter state estimates
    /*   P_matrix_old          previous Kalman filter error covariance matrix
    /*   meas_f_ib_b           measured specific force
    /*   est_L_b_old           previous latitude solution
    /*   TC_KF_config
    /*     .gyro_noise_PSD     Gyro noise PSD (rad^2/s)
    /*     .accel_noise_PSD    Accelerometer noise PSD (m^2 s^-3)
    /*     .accel_bias_PSD     Accelerometer bias random walk PSD (m^2 s^-5)
    /*     .gyro_bias_PSD      Gyro bias random walk PSD (rad^2 s^-3)
    /*     .clock_freq_PSD     Receiver clock frequency-drift PSD (m^2/s^3)
    /*     .clock_phase_PSD    Receiver clock phase-drift PSD (m^2/s)
    /*     .pseudo_range_SD    Pseudo-range measurement noise SD (m)
    /*     .range_rate_SD      Pseudo-range rate measurement noise SD (m/s)
    /*
    /* Outputs:
    /*   est_C_b_e_new     updated estimated body to ECEF coordinate
    /*                      transformation matrix
    /*   est_v_eb_e_new    updated estimated ECEF user velocity (m/s)
    /*   est_r_eb_e_new    updated estimated ECEF user position (m)
    /*   est_IMU_bias_new  updated estimated IMU biases
    /*     Rows 1-3          estimated accelerometer biases (m/s^2)
    /*     Rows 4-6          estimated gyro biases (rad/s)
    /*   est_clock_new     updated Kalman filter state estimates
    /*     Row 1             estimated receiver clock offset (m)
    /*     Row 2             estimated receiver clock drift (m/s)
    /*   P_matrix_new      updated Kalman filter error covariance matrix


    /* Copyright 2012, Paul Groves
    /* License: BSD; see license.txt for details  */
void TC_KF_Epoch(GNSS_measurements *GNSS_measurements, int no_meas,
                 const obsd_t *obs, const nav_t *nav,
                 double tor_s, double *est_C_b_e_old, double *est_v_eb_e_old, double *est_r_eb_e_old,
                 double *est_IMU_bias_old, double *est_clock_old, double *P_matrix_old,
                 double *meas_f_ib_b,
                 double est_L_b_old, TC_KF_config *TC_KF_config, double *est_C_b_e_new,
                 double *est_v_eb_e_new, double *est_r_eb_e_new, double *est_IMU_bias_new,
                 double *est_clock_new, double *P_matrix_new)
{

  /* Constants (sone of these could be changed to inputs at a later date) */
  double c = 299792458;   /* Speed of light in m/s */
  double omega_ie = OMGE; /* Earth rotation rate in rad/s */
  double R_0 = RE_WGS84;  /*WGS84 Equatorial radius in meters */
  double e = sqrt(e_2);   /*WGS84 eccentricity                        */

  /* Begins */
  double omega_ie_vec[3], Omega_ie[9] = {0.0};
  double *Phi_matrix, *Phi_transp, *Q_prime_matrix, *x_est_propagated, *x_est_new;
  double *Q_, *Q_aux, *H_matrix_transp;
  double *P_matrix_propagated, *P_matrix, *P_aux, *u_as_e_T, *pred_meas;
  double *H_matrix, *R_matrix, *ones, *delta_r, *delta_z, approx_range, range;
  double range_rate, C_e_I[9], *I_meas, *K_matrix, *K_matrix_inv;
  double meas_f_ib_e[3] = {0.0}, Skew_meas_f_ib_e[9] = {0.0}, est_r_eb_e[3], delta_r_aux[3];
  double geocentric_radius, g[3] = {0.0}, g_est_r_eb_e[9] = {0.0}, *I, *I_par;
  double Omega_x_est_r_eb_old[3], Omega_x_GNSS_meas_r[3], v_plus_Omega_x_est_r[3];
  double v_plus_Omega_x_GNSS_meas_r[3], est_v_plus_Omega_x_est_r[3];
  double C_x_v[3], C_less_v_plus_Omega[3];
  double *F1, *Q1, *K_x_delta_z, *K_x_H, *I_less_k_x_H;
  double Skew_x_est_new[9], I_less_Skew[9];
  double rho, u[3] = {0.0}, rate, vs[3], a[3], e_[3] = {0.0}, cosel, azel[2 * no_meas], pos[3] = {0.0};
  double E[9] = {0.0}, *var;
  int i, j, info;
  double dion = 0.0, vion = 0.0, dtrp = 0.0, vtrp = 0.0, lam_L1;

  printf("INPUT:\n");
  printf("no_meas=%d;\n", no_meas);
  printf("tor_s=%d;\n", no_meas);
  /* Loop measurements */
  for (j = 0; j < no_meas; j++)
  {
    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf;\n", GNSS_measurements[j].P[0],
           GNSS_measurements[j].D[0], GNSS_measurements[j].Sat_r_eb_e[0],
           GNSS_measurements[j].Sat_r_eb_e[1], GNSS_measurements[j].Sat_r_eb_e[2],
           GNSS_measurements[j].Sat_v_eb_e[0], GNSS_measurements[j].Sat_v_eb_e[1],
           GNSS_measurements[j].Sat_v_eb_e[2]);
  }
  printf("est_IMU_bias_old=[");
  for (i = 0; i < 6; i++)
  {
    printf("%lf, \n", est_IMU_bias_old[i]);
  }
  printf("]\n");

  printf("est_clock_old=[");
  for (i = 0; i < 2; i++)
  {
    printf("%lf, \n", est_clock_old[i]);
  }
  printf("]\n");
  printf("P_matrix_old=[");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%.13lf, ", P_matrix_old[i * 17 + j]);
    }
    printf(";\n");
  }
  printf("]\n");

  /* Initialize matrices and vectors */
  Phi_matrix = eye(17);
  Phi_transp = zeros(17, 17);
  Q_prime_matrix = zeros(17, 17);
  x_est_propagated = zeros(17, 1);
  x_est_new = zeros(17, 1);
  P_matrix_propagated = mat(17, 17);
  P_aux = mat(17, 17);
  Q_aux = mat(17, 17);
  P_matrix = mat(17, 17);
  u_as_e_T = zeros(no_meas, 3);
  pred_meas = zeros(no_meas, 2);
  H_matrix = zeros((2 * no_meas), 17);
  ones = mat(no_meas, 1);
  I_meas = eye(no_meas);
  H_matrix_transp = zeros(17, 2 * no_meas);
  R_matrix = zeros((2 * no_meas), (2 * no_meas));
  I_par = eye(17);
  K_matrix_inv = zeros((2 * no_meas), (2 * no_meas));
  K_matrix = mat(17, (2 * no_meas));
  delta_r = mat(no_meas, 3);
  delta_z = mat(2 * no_meas, 1);
  Q_ = mat(17, 17);
  F1 = mat(17, 2 * no_meas);
  Q1 = mat(2 * no_meas, 2 * no_meas);
  K_x_delta_z = mat(17, 1);
  K_x_H = mat(17, 17);
  I_less_k_x_H = mat(17, 17);
  I = eye(3);
  var = mat(no_meas, 2);

  for (i = 0; i < no_meas; i++)
    ones[i] = 1.0;
  for (i = 0; i < 3; i++)
    est_r_eb_e[i] = est_r_eb_e_old[i];

  /* Skew symmetric matrix of Earth rate */
  omega_ie_vec[0] = 0.0;
  omega_ie_vec[1] = 0.0;
  omega_ie_vec[2] = OMGE;
  Skew_symmetric(omega_ie_vec, Omega_ie);

  /* SYSTEM PROPAGATION PHASE */

  /* 1. Determine transition matrix using (14.50) (first-order approx) */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * 17 + j] = Phi_matrix[i * 17 + j] - Omega_ie[i * 3 + j] * tor_s;
    }
  }

  for (i = 0; i < 3; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Phi_matrix[i * 17 + j] = est_C_b_e_old[i * 3 + (j - 12)] * tor_s;
    }
  }

  matmul_row("NN", 3, 1, 3, 1.0, est_C_b_e_old, meas_f_ib_b, 0.0, meas_f_ib_e);
  Skew_symmetric(meas_f_ib_e, Skew_meas_f_ib_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * 17 + j] = -tor_s * Skew_meas_f_ib_e[(i - 3) * 3 + j];
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * 17 + j] = Phi_matrix[i * 17 + j] - 2 * Omega_ie[(i - 3) * 3 + (j - 3)] * tor_s;
    }
  }

  geocentric_radius = R_0 / sqrt(1 - pow((e * sin(est_L_b_old)), 2)) *
                      sqrt(pow(cos(est_L_b_old), 2) + pow((1 - e * e), 2) * pow(sin(est_L_b_old), 2)); /* from (2.137)*/
  Gravity_ECEF(est_r_eb_e_old, g);                                                                     //returns a vector
  matmul_row("NN", 3, 3, 1, 1.0, g, est_r_eb_e, 0.0, g_est_r_eb_e);
  //matmul("NT",3,3,1,1.0,est_r_eb_e_old,g,0.0, g_est_r_eb_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 6; j < 9; j++)
    {
      Phi_matrix[i * 17 + j] = -tor_s * 2 /
                               geocentric_radius * g_est_r_eb_e[(i - 3) * 3 + (j - 6)] / (norm(est_r_eb_e_old, 3));
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Phi_matrix[i * 17 + j] = est_C_b_e_old[(i - 3) * 3 + (j - 9)] * tor_s;
    }
  }

  for (i = 6; i < 9; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * 17 + j] = I[(i - 6) * 3 + (j - 3)] * tor_s;
    }
  }

  Phi_matrix[15 * 17 + 16] = tor_s;

  printf("Phi_matrix\n");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%.13lf ", Phi_matrix[i * 17 + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 2. Determine approximate system noise covariance matrix using (14.82) */
  //Q_prime_matrix(1:3,1:3) = eye(3) * TC_KF_config.gyro_noise_PSD * tor_s;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Q_prime_matrix[i * 17 + j] = (i == j ? I[(i)*3 + (j)] *
                                                 TC_KF_config->gyro_noise_PSD * tor_s
                                           : 0.0);
    }
  }
  //Q_prime_matrix(4:6,4:6) = eye(3) * TC_KF_config.accel_noise_PSD * tor_s;
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Q_prime_matrix[i * 17 + j] = (i == j ? I[(i - 3) * 3 + (j - 3)] *
                                                 TC_KF_config->accel_noise_PSD * tor_s
                                           : 0.0);
    }
  }
  //Q_prime_matrix(10:12,10:12) = eye(3) * TC_KF_config.accel_bias_PSD * tor_s;
  for (i = 9; i < 12; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Q_prime_matrix[i * 17 + j] = (i == j ? I[(i - 9) * 3 + (j - 9)] *
                                                 TC_KF_config->accel_bias_PSD * tor_s
                                           : 0.0);
    }
  }
  //Q_prime_matrix(13:15,13:15) = eye(3) * TC_KF_config.gyro_bias_PSD * tor_s;
  for (i = 12; i < 15; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Q_prime_matrix[i * 17 + j] = (i == j ? I[(i - 12) * 3 + (j - 12)] *
                                                 TC_KF_config->gyro_bias_PSD * tor_s
                                           : 0.0);
    }
  }

  //Q_prime_matrix(16,16) = TC_KF_config.clock_phase_PSD * tor_s;
  Q_prime_matrix[15 * 17 + 15] = TC_KF_config->clock_phase_PSD * tor_s;

  //Q_prime_matrix(17,17) = TC_KF_config->clock_freq_PSD * tor_s;
  Q_prime_matrix[16 * 17 + 16] = TC_KF_config->clock_freq_PSD * tor_s;
  /**/
  printf("Q_matrix\n");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%.13lf ", Q_prime_matrix[i * 17 + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 3. Propagate state estimates using (3.14) noting that only the clock
    /* states are non-zero due to closed-loop correction.

    /* Using values estimated in RTKLIB per epoch */
  for (i = 0; i < 17; i++)
    x_est_propagated[i] = 1E-20;
  x_est_propagated[15] = est_clock_old[0]; // + est_clock_old[1] * tor_s;
  x_est_propagated[16] = est_clock_old[1];
  //x_est_propagated[15] = 0.0; //+ est_clock_old[1] * tor_s;
  //x_est_propagated[16] = 0.0;

  /* 4. Propagate state estimation error covariance matrix using (3.46) */
  //P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *\
    //    Phi_matrix' + 0.5 * Q_prime_matrix;
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      P_aux[i * 17 + j] = P_matrix_old[i * 17 + j] + 0.5 * Q_prime_matrix[i * 17 + j];
    }
  }

  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      Phi_transp[j * 17 + i] = Phi_matrix[i * 17 + j];
    }
  }

  matmul_row("NN", 17, 17, 17, 1.0, P_aux, Phi_transp, 0.0, Q_aux); /*(Pp + 0.5*Q)*PHI^T */
  //matmul("TN",17,17,17,1.0,Phi_matrix,P_aux,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
  matmul_row("NN", 17, 17, 17, 1.0, Phi_matrix, Q_aux, 0.0, Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */
  /*matmul("NN",17,17,17,1.0,Q_aux,Phi_matrix,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      P_matrix_propagated[i * 17 + j] = Q_[i * 17 + j] + 0.5 * Q_prime_matrix[i * 17 + j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
    }
  }

  /* P_matrix_old is already the Propagated one*/
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      // P_matrix_propagated[i*17+j]= P_matrix_old[i*17+j];
    }
  }

  printf("P_matrix_old\n");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%.13lf ", P_matrix_old[i * 17 + j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("P_matrix_propagated\n");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%.13lf ", P_matrix_propagated[i * 17 + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* MEASUREMENT UPDATE PHASE */
  printf("CLOCK OFF. AND DRIFT: %lf, %lf\n", x_est_propagated[15], x_est_propagated[16]);

  ecef2pos(est_r_eb_e_old, pos);
  xyz2enu(pos, E);

  /* Loop measurements */
  for (j = 0; j < no_meas; j++)
  {                                      //no_meas represent the number of visible satellites
    azel[j * 2] = azel[1 + j * 2] = 0.0; /*{azimuth, elevation}*/
    var[j * 2] = var[1 + j * 2] = 0.0;

    /* Predict approx range */
    for (i = 0; i < 3; i++)
      delta_r[j * 3 + i] = GNSS_measurements[j].Sat_r_eb_e[i] - est_r_eb_e_old[i];
    //printf("\ndelta_r: %lf, %lf, %lf\n",delta_r[j*3], delta_r[j*3+1],delta_r[j*3+2]);

    rho = geodist(GNSS_measurements[j].Sat_r_eb_e, est_r_eb_e_old, u);
    satazel(pos, u, azel + j * 2);

    approx_range = norm(delta_r + (j * 3), 3);
    //printf("approx_range: %lf\n",approx_range);

    /* Calculate frame rotation during signal transit time using (8.36) */
    C_e_I[0] = 1;
    C_e_I[1] = omega_ie * approx_range / c;
    C_e_I[2] = 0;
    C_e_I[3] = -omega_ie * approx_range / c;
    C_e_I[4] = 1;
    C_e_I[5] = 0;
    C_e_I[6] = 0;
    C_e_I[7] = 0;
    C_e_I[8] = 1;

    /* Predict pseudo-range using (9.165) */
    matmul_row("NN", 3, 1, 3, 1.0, C_e_I, GNSS_measurements[j].Sat_r_eb_e, 0.0, delta_r_aux);
    //matmul("NN",1,3,3,1.0,GNSS_measurements[j].Sat_r_eb_e,C_e_I,0.0,delta_r_aux);

    /* ionospheric corrections */
    if (!ionocorr(obs[j].time, nav, obs[j].sat, pos, azel + j * 2,
                  IONOOPT_BRDC, &dion, &vion))
    {
      printf("ERROR IONO!\n");
    } //opt.ionoopt

    /* GPS-L1 -> L1/B1 */
    if ((lam_L1 = nav->lam[obs[j].sat - 1][0]) > 0.0)
    {
      dion *= (lam_L1 / lam_carr[0]) * (lam_L1 / lam_carr[0]);
    }

    /* tropospheric corrections */
    if (!tropcorr(obs[j].time, nav, pos, azel + j * 2,
                  TROPOPT_SAAS, &dtrp, &vtrp))
    { //opt.tropopt
      continue;
    }

    printf("IONO AND TROPO: %lf, %lf\n", dion, dtrp);

    for (i = 0; i < 3; i++)
      delta_r[j * 3 + i] = delta_r_aux[i] - est_r_eb_e_old[i];
    range = norm(delta_r + (j * 3), 3);

    printf("Geometric distance - Rtklib: %lf Grove: %lf\n", rho, range);

    range = rho;

    pred_meas[j * 2 + 0] = range + x_est_propagated[15];

    /* Predict line of sight */
    for (i = 0; i < 3; i++)
      u_as_e_T[j * 3 + i] = delta_r[j * 3 + i] / range;
    for (i = 0; i < 3; i++)
      u_as_e_T[j * 3 + i] = u[i];

    /* Predict pseudo-range rate using (9.165) */
    /*range_rate = u_as_e_T[j*3+i] * (C_e_I[] * (GNSS_measurements[j,6:8]' +...
            Omega_ie * GNSS_measurements[j,3:5]') - (est_v_eb_e_old +...
            Omega_ie * est_r_eb_e_old));*/

    /* RTKLIB range-rate */
    /* line-of-sight vector in ecef */
    cosel = cos(azel[1 + j * 2]);
    a[0] = sin(azel[j * 2]) * cosel;
    a[1] = cos(azel[j * 2]) * cosel;
    a[2] = sin(azel[1 + j * 2]);
    matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e_); // This one is from RTKLIB (it is already in row-order)

    /* satellite velocity relative to receiver in ecef */
    for (i = 0; i < 3; i++)
      vs[i] = GNSS_measurements[j].Sat_v_eb_e[i] - est_v_eb_e_old[i];

    /* range rate with earth rotation correction */
    rate = dot(vs, e_, 3) + OMGE / CLIGHT * (GNSS_measurements[j].Sat_v_eb_e[1] * est_r_eb_e_old[0] + GNSS_measurements[j].Sat_r_eb_e[1] * est_v_eb_e_old[0] - GNSS_measurements[j].Sat_v_eb_e[0] * est_r_eb_e_old[1] - GNSS_measurements[j].Sat_r_eb_e[0] * est_v_eb_e_old[1]);

    /* doppler residual
            v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);*/

    matmul_row("NN", 3, 1, 3, 1.0, Omega_ie, est_r_eb_e_old, 0.0, Omega_x_est_r_eb_old);
    //matmul("NN",1,3,3,1.0,est_r_eb_e_old,Omega_ie,0.0, Omega_x_est_r_eb_old);
    matmul_row("NN", 3, 1, 3, 1.0, Omega_ie, GNSS_measurements[j].Sat_r_eb_e, 0.0, Omega_x_GNSS_meas_r);
    //matmul("NN",1,3,3,1.0,GNSS_measurements[j].Sat_r_eb_e,Omega_ie,0.0,Omega_x_GNSS_meas_r);

    for (i = 0; i < 3; i++)
      v_plus_Omega_x_GNSS_meas_r[i] =
          GNSS_measurements[j].Sat_v_eb_e[i] + Omega_x_GNSS_meas_r[i];

    for (i = 0; i < 3; i++)
      est_v_plus_Omega_x_est_r[i] =
          est_v_eb_e_old[i] + Omega_x_est_r_eb_old[i];

    matmul_row("NN", 3, 1, 3, 1.0, C_e_I, v_plus_Omega_x_GNSS_meas_r, 0.0, C_x_v);
    //matmul("NN",1,3,3,1.0,v_plus_Omega_x_GNSS_meas_r,C_e_I,0.0,C_x_v);

    for (i = 0; i < 3; i++)
      C_less_v_plus_Omega[i] =
          C_x_v[i] - est_v_plus_Omega_x_est_r[i];

    matmul_row("NN", 1, 1, 3, 1.0, u_as_e_T + (j * 3), C_less_v_plus_Omega, 0.0, &range_rate);
    //matmul("NN",1,1,3,1.0,C_less_v_plus_Omega,u_as_e_T+(j*3),0.0,&range_rate);

    printf("Geometric distance rate - Rtklib: %lf Grove: %lf\n", rate, range_rate);

    range_rate = rate;

    pred_meas[j * 2 + 1] = range_rate + x_est_propagated[16];

    /* error variance */
    var[j * 2 + 0] = pow(TC_KF_config->pseudo_range_SD / sin(azel[1 + j * 2]), 2);
    var[j * 2 + 1] = pow(TC_KF_config->range_rate_SD / sin(azel[1 + j * 2]), 2);

  } //end for j

  /* 5. Set-up measurement matrix using (14.126) */

  //H_matrix(1:no_meas,7:9) = u_as_e_T(1:no_meas,1:3); - Position
  for (i = 0; i < no_meas; i++)
  {
    for (j = 6; j < 9; j++)
    {
      H_matrix[i * 17 + j] = u_as_e_T[i * 3 + (j - 6)];
    }
  }

  //H_matrix(1:no_meas,16) = ones(no_meas,1);  - Clock offset
  for (i = 0; i < no_meas; i++)
    H_matrix[i * 17 + 15] = 1.0;

  //H_matrix((no_meas + 1):(2 * no_meas),4:6) = u_as_e_T(1:no_meas,1:3); - Velocity
  for (i = no_meas; i < 2 * no_meas; i++)
  {
    for (j = 3; j < 6; j++)
    {
      H_matrix[i * 17 + j] = u_as_e_T[(i - no_meas) * 3 + (j - 3)];
    }
  }

  //H_matrix((no_meas + 1):(2 * no_meas),17) = ones(no_meas,1);  -  Clock drift
  for (i = no_meas; i < 2 * no_meas; i++)
  {
    H_matrix[i * 17 + 16] = 1.0;
  }

  printf("H_matrix\n");
  for (i = 0; i < no_meas * 2; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%.12lf ", H_matrix[i * 17 + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 6. Set-up measurement noise covariance matrix assuming all measurements
    /* are independent and have equal variance for a given measurement type. */

  //R_matrix(1:no_meas,1:no_meas) = eye(no_meas) *TC_KF_config.pseudo_range_SD^2;
  for (i = 0; i < no_meas; i++)
  {
    for (j = 0; j < no_meas; j++)
    {
      R_matrix[i * (2 * no_meas) + j] = I_meas[i * (no_meas) + j] *
                                        pow(TC_KF_config->pseudo_range_SD, 2);
    }
  }
  for (i = 0; i < no_meas; i++)
  {
    for (j = 0; j < no_meas; j++)
    {
      R_matrix[i * (2 * no_meas) + j] = (i == j ? var[i * 2] : 0.0);
    }
  }

  //R_matrix(1:no_meas,(no_meas + 1):(2 * no_meas)) = zeros(no_meas);
  for (i = 0; i < no_meas; i++)
  {
    for (j = no_meas; j < 2 * no_meas; j++)
    {
      R_matrix[i * (2 * no_meas) + j] = 0.0;
    }
  }

  //R_matrix((no_meas + 1):(2 * no_meas),1:no_meas) =  zeros(no_meas);
  for (i = no_meas; i < 2 * no_meas; i++)
  {
    for (j = 0; j < no_meas; j++)
    {
      R_matrix[i * (2 * no_meas) + j] = 0.0;
    }
  }

  //R_matrix((no_meas + 1):(2 * no_meas),(no_meas + 1):(2 * no_meas)) =...
  //  eye(no_meas) * TC_KF_config.range_rate_SD^2;
  for (i = no_meas; i < 2 * no_meas; i++)
  {
    for (j = no_meas; j < 2 * no_meas; j++)
    {
      R_matrix[i * (2 * no_meas) + j] = I_meas[(i - no_meas) * (no_meas) + (j - no_meas)] *
                                        pow(TC_KF_config->range_rate_SD, 2);
    }
  }
  for (i = no_meas; i < 2 * no_meas; i++)
  {
    for (j = no_meas; j < 2 * no_meas; j++)
    {
      R_matrix[i * (2 * no_meas) + j] = (i == j ? var[(i - no_meas) * 2 + 1] : 0.0);
    }
  }

  printf("R_meas_noise\n");
  for (i = 0; i < no_meas * 2; i++)
  {
    for (j = 0; j < no_meas * 2; j++)
    {
      printf("%.12lf ", R_matrix[i * (no_meas * 2) + j]);
    }
    printf("\n");
  }
  printf("\n");

  for (i = 0; i < 2 * no_meas; i++)
  {
    for (j = 0; j < 17; j++)
    {
      H_matrix_transp[j * 2 * no_meas + i] = H_matrix[i * 17 + j];
    }
  }

  /* 7. Calculate Kalman gain using (3.21) */
  //K_matrix = P_matrix_propagated * H_matrix' * inv(H_matrix *...
  //    P_matrix_propagated * H_matrix' + R_matrix);
  matmul_row("NN", 17, 2 * no_meas, 17, 1.0, P_matrix_propagated, H_matrix_transp, 0.0, F1);

  //matmul("TN",2*no_meas,17,17,1.0,H_matrix,P_matrix_propagated,0.0,F1);
  matmul_row("NN", 2 * no_meas, 2 * no_meas, 17, 1.0, H_matrix, F1, 0.0, Q1);
  //matmul("NN",2*no_meas,2*no_meas,17,1.0,F1,H_matrix,0.0,Q1);
  for (i = 0; i < 2 * no_meas; i++)
  {
    for (j = 0; j < 2 * no_meas; j++)
    {
      K_matrix_inv[i * (2 * no_meas) + j] = Q1[i * (2 * no_meas) + j] + R_matrix[i * (2 * no_meas) + j];
    }
  }

  if (!(info = matinv(K_matrix_inv, 2 * no_meas)))
  {
    printf("Invertion status, if 0 > info:error! info:%d\n", info);
    matmul_row("NN", 17, 2 * no_meas, 2 * no_meas, 1.0, F1, K_matrix_inv, 0.0, K_matrix);
    //matmul("NN",2*no_meas,17,2*no_meas,1.0,K_matrix_inv,F1,0.0,K_matrix);
    printf("\n\n****************  INVERTING MATRIX  ***********************\n\n");
  }
  else
  {
    printf("Invertion status, if 0 > info:error! info:%d\n", info);
    printf("\n\n*************** ERROR INVERTING MATRIX*********************\n\n");
  }

  printf("K_gain\n");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < no_meas * 2; j++)
    {
      printf("%.12lf ", K_matrix[i * (no_meas * 2) + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 8. Formulate measurement innovations using (14.119) */

  //delta_z(1:no_meas,1) = GNSS_measurements(1:no_meas,1) - pred_meas(1:no_meas,1);
  //delta_z((no_meas + 1):(2 * no_meas),1) = GNSS_measurements(1:no_meas,2) -...
  //pred_meas(1:no_meas,2);

  /* Clock jump correction  */
  double delta_z_sum = 0.0;
  for (i = 0; i < no_meas; i++)
  {
    delta_z[i] = GNSS_measurements[i].P[0] - pred_meas[i * 2];
    delta_z_sum += delta_z[i];
  }
  if (fabs(delta_z_sum / no_meas) > 150000.0)
  {
    printf("CLOCK JUMP CORR.: %lf >150000.0, CLK: %lf \n",
           fabs(delta_z_sum / no_meas), x_est_propagated[15]);
  }

  for (i = no_meas; i < 2 * no_meas; i++)
  {
    delta_z[i] = GNSS_measurements[i - no_meas].D[0] - pred_meas[(i - no_meas) * 2 + 1];
  }

  printf("Measurement vector:\n");
  for (i = 0; i < no_meas; i++)
    printf("%lf ", GNSS_measurements[i].P[0]);
  for (i = 0; i < no_meas; i++)
    printf("%lf ", GNSS_measurements[i].D[0]);
  printf("\n");

  printf("Predicted measurement vector:\n");
  for (i = 0; i < no_meas; i++)
    printf("%lf ", pred_meas[i * 2]);
  for (i = 0; i < no_meas; i++)
    printf("%lf ", pred_meas[i * 2 + 1]);
  printf("\n");

  printf("Innovation vector:\n");
  for (i = 0; i < 2 * no_meas; i++)
    printf("%lf ", delta_z[i]);

  for (i = 0; i < no_meas; i++)
  {
    fprintf(out_KF_residuals, "%lf %2d %lf %lf\n", GNSS_measurements->sec,
            GNSS_measurements[i].sat, delta_z[i], delta_z[i + no_meas]);
  }

  /* 9. Update state estimates using (3.24) */

  matmul_row("NN", 17, 1, 2 * no_meas, 1.0, K_matrix, delta_z, 0.0, K_x_delta_z);

  for (i = 0; i < 17; i++)
    x_est_new[i] = x_est_propagated[i] + K_x_delta_z[i];

  printf("\nx_new\n");
  for (i = 0; i < 17; i++)
    printf("%lf ", x_est_new[i]);
  printf("\n");
  /*
      fprintf(out_KF_state_error,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
    %lf %lf %lf %lf %lf %lf %lf\n", GNSS_measurements->sec,\
    K_x_delta_z[0],K_x_delta_z[1],K_x_delta_z[2],K_x_delta_z[3],K_x_delta_z[4],K_x_delta_z[5],\
    K_x_delta_z[6],K_x_delta_z[7],K_x_delta_z[8],K_x_delta_z[9],K_x_delta_z[10],\
    K_x_delta_z[11],K_x_delta_z[12],K_x_delta_z[13],K_x_delta_z[14],K_x_delta_z[15],\
    K_x_delta_z[16], K_x_delta_z[17] );*/

  fprintf(out_KF_state_error, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
    %lf %lf %lf %lf %lf %lf %lf\n",
          GNSS_measurements->sec,
          x_est_new[0], x_est_new[1], x_est_new[2], x_est_new[3], x_est_new[4], x_est_new[5],
          x_est_new[6], x_est_new[7], x_est_new[8], x_est_new[9], x_est_new[10],
          x_est_new[11], x_est_new[12], x_est_new[13], x_est_new[14], x_est_new[15],
          x_est_new[16], x_est_new[17]);

  /* 10. Update state estimation error covariance matrix using (3.25) */
  //P_matrix_new = (eye(17) - K_matrix * H_matrix) * P_matrix_propagated;
  matmul_row("NN", 17, 17, 2 * no_meas, 1.0, K_matrix, H_matrix, 0.0, K_x_H);

  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      I_less_k_x_H[i * 17 + j] = I_par[i * 17 + j] - K_x_H[i * 17 + j];
    }
  }

  matmul_row("NN", 17, 17, 17, 1.0, I_less_k_x_H, P_matrix_propagated, 0.0, P_matrix_new);

  /* CLOSED-LOOP CORRECTION */

  /* Correct attitude, velocity, and position using (14.7-9) */

  /* Quaternion attitude correction */
  Quaternion_attitude_errror_correction(x_est_new, est_C_b_e_old, est_C_b_e_new);

  Skew_symmetric(x_est_new, Skew_x_est_new);
  for (i = 0; i < 9; i++)
    I_less_Skew[i] = I[i] - Skew_x_est_new[i];
  matmul_row("NN", 3, 3, 3, 1.0, I_less_Skew, est_C_b_e_old, 0.0, est_C_b_e_new);
  //matmul("NN",3,3,3,1.0,est_C_b_e_old,I_less_Skew,0.0,est_C_b_e_new);

  for (i = 3; i < 6; i++)
    est_v_eb_e_new[i - 3] = est_v_eb_e_old[i - 3] - x_est_new[i];
  for (i = 6; i < 9; i++)
    est_r_eb_e_new[i - 6] = est_r_eb_e_old[i - 6] - x_est_new[i];

  /* Update IMU bias and GNSS receiver clock estimates */
  for (i = 0; i < 6; i++)
    est_IMU_bias_new[i] = est_IMU_bias_old[i] + x_est_new[i + 9];
  est_clock_new[0] = x_est_new[15];
  est_clock_new[1] = x_est_new[16];

  printf("OUTPUT:\n");

  double C_Transp[9], est_C_b_n_T[9], est_L_b = {0.0}, est_lambda_b = {0.0}, est_h_b = {0.0}, est_v_eb_n[3] = {0.0}, est_C_b_n[9] = {0.0};
  double euler_angles[3] = {0.0};

  ECEF_to_NED(est_r_eb_e_new, est_v_eb_e_new, est_C_b_e_new,
              &est_L_b, &est_lambda_b, &est_h_b,
              est_v_eb_n, est_C_b_n);

  printf("est_C_b_n\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%.13lf ", est_C_b_n[i * 3 + j]);
    }
    printf("\n");
  }
  printf("\n");

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      C_Transp[j * 3 + i] = est_C_b_n[i * 3 + j];
    }
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      est_C_b_n_T[i * 3 + j] = C_Transp[i * 3 + j];
    }
  }

  /* */
  CTM_to_Euler(euler_angles, est_C_b_n_T);

  printf("Euler from CTM: %lf, %lf, %lf\n", euler_angles[0],
         euler_angles[1], euler_angles[2]);

  /* Loop measurements */
  printf("est_C_b_e_new\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf ", est_C_b_e_new[i * 3 + j]);
    }
    printf("\n");
  }
  printf("est_v_eb_e_new");
  for (i = 0; i < 3; i++)
  {
    printf("%lf ", est_v_eb_e_new[i]);
  }
  printf("\n");
  printf("est_r_eb_e_new");
  for (i = 0; i < 3; i++)
  {
    printf("%lf ", est_r_eb_e_new[i]);
  }
  printf("\n");
  printf("est_IMU_bias_new:");
  for (i = 0; i < 2; i++)
  {
    printf("%lf ", est_IMU_bias_new[i]);
  }
  printf("\n");
  printf("est_clock_new:");
  for (i = 0; i < 2; i++)
  {
    printf("%lf, \n", est_clock_new[i]);
  }
  printf("\n");

  printf("P_new\n");
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      printf("%lf ", P_matrix_new[i * 17 + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* Ends */

  /* Freeing memory */
  free(Phi_matrix);
  free(Phi_transp);
  free(Q_prime_matrix);
  free(I);
  free(x_est_propagated);
  free(x_est_new);
  free(H_matrix_transp);
  free(P_matrix_propagated);
  free(P_aux);
  free(Q_aux);
  free(P_matrix);
  free(u_as_e_T);
  free(pred_meas);
  free(ones);
  free(I_par);
  free(K_matrix);
  free(delta_r);
  free(delta_z);
  free(Q_);
  free(F1);
  free(Q1);
  free(K_x_delta_z);
  free(I_meas);
  free(R_matrix);
  free(K_matrix_inv);
  free(K_x_H);
  free(I_less_k_x_H);
  free(H_matrix);
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
/* LC_KF_Epoch - Implements one cycle of the loosely coupled INS/GNSS
% Kalman filter plus closed-loop correction of all inertial states
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 12/4/2012 by Paul Groves
%
% Inputs:
%   GNSS_r_eb_e           GNSS estimated ECEF user position (m)
%   GNSS_v_eb_e           GNSS estimated ECEF user velocity (m/s)
%   tor_s                 propagation interval (s)
%   est_C_b_e_old         prior estimated body to ECEF coordinate
%                         transformation matrix
%   est_v_eb_e_old        prior estimated ECEF user velocity (m/s)
%   est_r_eb_e_old        prior estimated ECEF user position (m)
%   est_IMU_bias_old      prior estimated IMU biases (body axes)
%   P_matrix_old          previous Kalman filter error covariance matrix
%   meas_f_ib_b           measured specific force
%   est_L_b_old           previous latitude solution
%   LC_KF_config
%     .gyro_noise_PSD     Gyro noise PSD (rad^2/s)
%     .accel_noise_PSD    Accelerometer noise PSD (m^2 s^-3)
%     .accel_bias_PSD     Accelerometer bias random walk PSD (m^2 s^-5)
%     .gyro_bias_PSD      Gyro bias random walk PSD (rad^2 s^-3)
%     .pos_meas_SD            Position measurement noise SD per axis (m)
%     .vel_meas_SD            Velocity measurement noise SD per axis (m/s)
%
% Outputs:
%   est_C_b_e_new     updated estimated body to ECEF coordinate
%                      transformation matrix
%   est_v_eb_e_new    updated estimated ECEF user velocity (m/s)
%   est_r_eb_e_new    updated estimated ECEF user position (m)
%   est_IMU_bias_new  updated estimated IMU biases
%     Rows 1-3          estimated accelerometer biases (m/s^2)
%     Rows 4-6          estimated gyro biases (rad/s)
%   P_matrix_new      updated Kalman filter error covariance matrix


% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details  */
void LC_KF_Epoch(double time, double *GNSS_r_eb_e, double *GNSS_v_eb_e, int no_meas,
                 double tor_s, double *est_C_b_e_old, double *est_v_eb_e_old, double *est_r_eb_e_old,
                 double *est_IMU_bias_old, double *P_matrix_old, double *meas_f_ib_b,
                 double est_L_b_old, LC_KF_config *LC_KF_config, double *est_C_b_e_new,
                 double *est_v_eb_e_new, double *est_r_eb_e_new, double *est_IMU_bias_new,
                 double *P_matrix_new)
{

  /* Constants (sone of these could be changed to inputs at a later date) */
  double c = 299792458;   /* Speed of light in m/s */
  double omega_ie = OMGE; /* Earth rotation rate in rad/s */
  double R_0 = RE_WGS84;  /*WGS84 Equatorial radius in meters */
  double e = sqrt(e_2);   /*WGS84 eccentricity                        */

  /* Begins */
  double omega_ie_vec[3], Omega_ie[9];
  double *Phi_matrix, *Q_prime_matrix, *x_est_propagated, *x_est_new;
  double *Q_, *Q_aux, *Phi_transp, *H_matrix_transp;
  double *P_matrix_propagated, *P_matrix, *P_aux, *u_as_e_T, *pred_meas;
  double *H_matrix, *R_matrix, *ones, *delta_r, *delta_z, approx_range, range;
  double range_rate, C_e_I[9], *I_meas, *K_matrix, *K_matrix_inv;
  double meas_f_ib_e[3], Skew_meas_f_ib_e[9], est_r_eb_e[3], delta_r_aux[3];
  double geocentric_radius, g[3], g_est_r_eb_e[9], *I, *I_par;
  double Omega_x_est_r_eb_old[3], Omega_x_GNSS_meas_r[3], v_plus_Omega_x_est_r[3];
  double v_plus_Omega_x_GNSS_meas_r[3], est_v_plus_Omega_x_est_r[3];
  double C_x_v[3], C_less_v_plus_Omega[3];
  double *F1, *Q1, *K_x_delta_z, *K_x_H, *I_less_k_x_H;
  double Skew_x_est_new[9], I_less_Skew[9];
  double rho, u[3], rate, vs[3], a[3], e_[3], cosel, azel[2 * no_meas], pos[3];
  double E[9];
  int i, j, info;

  /* Initialize matrices and vectors */
  Phi_matrix = eye(15);
  Phi_transp = zeros(15, 15);
  Q_prime_matrix = zeros(15, 15);
  x_est_propagated = zeros(15, 1);
  x_est_new = zeros(15, 1);
  P_matrix_propagated = mat(15, 15);
  P_aux = mat(15, 15);
  Q_aux = mat(15, 15);
  P_matrix = mat(15, 15);
  u_as_e_T = zeros(no_meas, 3);
  pred_meas = zeros(no_meas, 1);
  H_matrix = zeros((no_meas), 15);
  ones = mat(no_meas, 1);
  I_meas = eye(no_meas);
  H_matrix_transp = zeros(15, no_meas);
  R_matrix = zeros((no_meas), (no_meas));
  I_par = eye(15);
  K_matrix_inv = zeros((no_meas), (no_meas));
  K_matrix = mat(15, (no_meas));
  delta_r = mat(no_meas, 3);
  delta_z = mat(no_meas, 1);
  Q_ = mat(15, 15);
  F1 = mat(15, no_meas);
  Q1 = mat(no_meas, no_meas);
  K_x_delta_z = mat(15, 1);
  K_x_H = mat(15, 15);
  I_less_k_x_H = mat(15, 15);
  I = eye(3);
  /*
    printf("INPUT:\n");
    printf("no_meas=%d;\n",no_meas);
    printf("tor_s=%lf;\n",tor_s);

    printf("est_IMU_bias_old=[");
    for (i = 0; i < 6; i++) {
      printf("%lf, \n", est_IMU_bias_old[i]);
        }
    printf("]\n");

    printf("]\n");
    printf("P_matrix_old=[");
    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        printf("%.13lf, ", P_matrix_old[i*15+j]);
      }
      printf(";\n");
    }
    printf("]\n");
*/
  for (i = 0; i < no_meas; i++)
    ones[i] = 1.0;
  for (i = 0; i < 3; i++)
    est_r_eb_e[i] = est_r_eb_e_old[i];

  /* Skew symmetric matrix of Earth rate */
  omega_ie_vec[0] = 0;
  omega_ie_vec[1] = 0;
  omega_ie_vec[2] = OMGE;
  Skew_symmetric(omega_ie_vec, Omega_ie);

  /* SYSTEM PROPAGATION PHASE */

  /* 1. Determine transition matrix using (14.50) (first-order approx) */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * 15 + j] = Phi_matrix[i * 15 + j] - Omega_ie[i * 3 + j] * tor_s;
    }
  }

  for (i = 0; i < 3; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Phi_matrix[i * 15 + j] = est_C_b_e_old[i * 3 + (j - 12)] * tor_s;
    }
  }

  matmul_row("NN", 3, 1, 3, 1.0, est_C_b_e_old, meas_f_ib_b, 0.0, meas_f_ib_e);
  Skew_symmetric(meas_f_ib_e, Skew_meas_f_ib_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * 15 + j] = -tor_s * Skew_meas_f_ib_e[(i - 3) * 3 + j];
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * 15 + j] = Phi_matrix[i * 15 + j] - 2 * Omega_ie[(i - 3) * 3 + (j - 3)] * tor_s;
    }
  }

  geocentric_radius = R_0 / sqrt(1 - pow((e * sin(est_L_b_old)), 2)) *
                      sqrt(pow(cos(est_L_b_old), 2) + pow((1 - e * e), 2) * pow(sin(est_L_b_old), 2)); /* from (2.137)*/
  Gravity_ECEF(est_r_eb_e_old, g);
  //matmul_row("NT",3,3,1,1.0,g,est_r_eb_e,0.0, g_est_r_eb_e);
  matmul_row("NN", 3, 3, 1, 1.0, g, est_r_eb_e, 0.0, g_est_r_eb_e);
  for (i = 3; i < 6; i++)
  {
    for (j = 6; j < 9; j++)
    {
      Phi_matrix[i * 15 + j] = -tor_s * 2 /
                               geocentric_radius * g_est_r_eb_e[(i - 3) * 3 + (j - 6)] / (norm(est_r_eb_e_old, 3));
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Phi_matrix[i * 15 + j] = est_C_b_e_old[(i - 3) * 3 + (j - 9)] * tor_s;
    }
  }

  for (i = 6; i < 9; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * 15 + j] = I[(i - 6) * 3 + (j - 3)] * tor_s;
    }
  }
  /*
    printf("Phi_matrix\n");
    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        printf("%.13lf ", Phi_matrix[i*15+j]);
      }
      printf("\n");
    }
    printf("\n");
*/

  /* 2. Determine approximate system noise covariance matrix using (14.82) */
  //Q_prime_matrix(1:3,1:3) = eye(3) * LC_KF_config.gyro_noise_PSD * tor_s;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Q_prime_matrix[i * 15 + j] = (i == j ? I[(i)*3 + (j)] *
                                                 LC_KF_config->gyro_noise_PSD * tor_s
                                           : 0.0);
    }
  }
  //Q_prime_matrix(4:6,4:6) = eye(3) * LC_KF_config.accel_noise_PSD * tor_s;
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Q_prime_matrix[i * 15 + j] = (i == j ? I[(i - 3) * 3 + (j - 3)] *
                                                 LC_KF_config->accel_noise_PSD * tor_s
                                           : 0.0);
    }
  }
  //Q_prime_matrix(10:12,10:12) = eye(3) * LC_KF_config.accel_bias_PSD * tor_s;
  for (i = 9; i < 12; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Q_prime_matrix[i * 15 + j] = (i == j ? I[(i - 9) * 3 + (j - 9)] *
                                                 LC_KF_config->accel_bias_PSD * tor_s
                                           : 0.0);
    }
  }
  //Q_prime_matrix(13:15,13:15) = eye(3) * LC_KF_config.gyro_bias_PSD * tor_s;
  for (i = 12; i < 15; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Q_prime_matrix[i * 15 + j] = (i == j ? I[(i - 12) * 3 + (j - 12)] *
                                                 LC_KF_config->gyro_bias_PSD * tor_s
                                           : 0.0);
    }
  }
  /*
    printf("Q_matrix\n");
    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        printf("%.13lf ", Q_prime_matrix[i*15+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
  /* 3. Propagate state estimates using (3.14) noting that all states are zero
       due to closed-loop correction.  */
  for (i = 0; i < 15; i++)
    x_est_propagated[i] = 0.0;

  /* 4. Propagate state estimation error covariance matrix using (3.46) */
  //P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *\
    //    Phi_matrix' + 0.5 * Q_prime_matrix;
  for (i = 0; i < 15; i++)
  {
    for (j = 0; j < 15; j++)
    {
      P_aux[i * 15 + j] = P_matrix_old[i * 15 + j] + 0.5 * Q_prime_matrix[i * 15 + j];
    }
  }

  for (i = 0; i < 15; i++)
  {
    for (j = 0; j < 15; j++)
    {
      Phi_transp[j * 15 + i] = Phi_matrix[i * 15 + j];
    }
  }

  matmul_row("NN", 15, 15, 15, 1.0, P_aux, Phi_transp, 0.0, Q_aux); /* (Pp + 0.5*Q)*PHI^T */
  matmul_row("NN", 15, 15, 15, 1.0, Phi_matrix, Q_aux, 0.0, Q_);    /* PHI*(Pp + 0.5*Q)*PHI^T */

  for (i = 0; i < 15; i++)
  {
    for (j = 0; j < 15; j++)
    {
      P_matrix_propagated[i * 15 + j] = Q_[i * 15 + j] + 0.5 * Q_prime_matrix[i * 15 + j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
    }
  }
  /*
    printf("P_matrix_old\n");
    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        printf("%.13lf ", P_matrix_old[i*15+j]);
      }
      printf("\n");
    }
    printf("\n");

    printf("P_matrix_propagated\n");
    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        printf("%.13lf ", P_matrix_propagated[i*15+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
  /* MEASUREMENT UPDATE PHASE */

  /* 5. Set-up measurement matrix using (14.115)  */

  //H_matrix(1:3,7:9) = -eye(3);
  for (i = 0; i < 3; i++)
  {
    for (j = 6; j < 9; j++)
    {
      H_matrix[i * 15 + j] = (i == j ? -1.0 : 0.0);
    }
  }

  //H_matrix(4:6,4:6) = -eye(3);
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      H_matrix[i * 15 + j] = (i == j ? -1.0 : 0.0);
    }
  }
  /*
   printf("H_matrix\n");
   for (i = 0; i < no_meas; i++) {
     for (j = 0; j < 15; j++) {
       printf("%.12lf ", H_matrix[i*15+j]);
     }
     printf("\n");
   }
   printf("\n");
*/
  /* 6. Set-up measurement noise covariance matrix assuming all components of
    GNSS position and velocity are independent and have equal variance. */

  //R_matrix(1:3,1:3) = eye(3) * LC_KF_config.pos_meas_SD^2;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      R_matrix[i * (no_meas) + j] = I_meas[i * (no_meas) + j] *
                                    pow(LC_KF_config->pos_meas_SD, 2);
    }
  }

  //R_matrix(1:3,4:6) = zeros(3);
  for (i = 0; i < 3; i++)
  {
    for (j = 3; j < 6; j++)
    {
      R_matrix[i * (no_meas) + j] = 0.0;
    }
  }

  //R_matrix(4:6,1:3) = zeros(3);
  for (i = 3; i < 6; i++)
  {
    for (j = 0; j < 3; j++)
    {
      R_matrix[i * (no_meas) + j] = 0.0;
    }
  }

  //R_matrix(4:6,4:6) = eye(3) * LC_KF_config.vel_meas_SD^2;
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      R_matrix[i * (no_meas) + j] = I_meas[(i - 3) * (no_meas) + (j - 3)] *
                                    pow(LC_KF_config->vel_meas_SD, 2);
    }
  }
  /*
    printf("R_meas_noise\n");
    for (i = 0; i < no_meas; i++) {
      for (j = 0; j < no_meas; j++) {
        printf("%.12lf ", R_matrix[i*(no_meas)+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
  for (i = 0; i < no_meas; i++)
  {
    for (j = 0; j < 15; j++)
    {
      H_matrix_transp[j * no_meas + i] = H_matrix[i * 15 + j];
    }
  }

  /* 7. Calculate Kalman gain using (3.21) */
  //K_matrix = P_matrix_propagated * H_matrix' * inv(H_matrix *...
  //    P_matrix_propagated * H_matrix' + R_matrix);
  matmul_row("NN", 15, 6, 15, 1.0, P_matrix_propagated, H_matrix_transp, 0.0, F1);
  matmul_row("NN", 6, 6, 15, 1.0, H_matrix, F1, 0.0, Q1);

  for (i = 0; i < 6; i++)
  {
    for (j = 0; j < 6; j++)
    {
      K_matrix_inv[i * (6) + j] = Q1[i * (6) + j] + R_matrix[i * (6) + j];
    }
  }

  if (!(info = matinv(K_matrix_inv, 6)))
  {
    printf("Invertion status, if 0 > info:error! info:%d\n", info);
    matmul_row("NN", 15, 6, 6, 1.0, F1, K_matrix_inv, 0.0, K_matrix);
    printf("\n\n\n\n********************************** INVERTING MATRIX  **************************************\n\n\n\n\n");
  }
  else
  {
    printf("Invertion status, if 0 > info:error! info:%d\n", info);
    printf("\n\n\n\n********************************** ERROR INVERTING MATRIX**************************************\n\n\n\n\n");
  }
  /*
    printf("K_gain\n");
    for (i = 0; i < 15; i++) {
      for (j = 0; j < no_meas; j++) {
        printf("%.12lf ", K_matrix[i*(no_meas)+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
  /* 8. Formulate measurement innovations using (14.102), noting that zero
    % lever arm is assumed here
  printf("GNSS pos and vel: %lf, %lf, %lf, %lf, %lf, %lf\n",GNSS_r_eb_e[0],GNSS_r_eb_e[1],GNSS_r_eb_e[2],GNSS_v_eb_e[0], \
  GNSS_v_eb_e[1],GNSS_v_eb_e[2]  );

  printf("INS pos and vel: %lf, %lf, %lf, %lf, %lf, %lf\n",est_r_eb_e_old[0],est_r_eb_e_old[1],est_r_eb_e_old[2],est_v_eb_e_old[0], \
est_v_eb_e_old[1],est_v_eb_e_old[2]  );
*/
  for (i = 0; i < 3; i++)
  {
    delta_z[i] = GNSS_r_eb_e[i] - est_r_eb_e_old[i];
  }

  for (i = 3; i < 6; i++)
  {
    delta_z[i] = GNSS_v_eb_e[i - 3] - est_v_eb_e_old[i - 3];
  }
  /* Set a form for position and velocity residuals
    for (i = 0; i < no_meas; i++) {
      fprintf(out_KF_residuals,"%lf %2d %lf %lf\n", time,\
      GNSS_measurements[i].sat, delta_z[i], delta_z[i+no_meas] );
    }
*/

  /* 9. Update state estimates using (3.24) */
  matmul_row("NN", 15, 1, 6, 1.0, K_matrix, delta_z, 0.0, K_x_delta_z);
  for (i = 0; i < 15; i++)
    x_est_new[i] = x_est_propagated[i] + K_x_delta_z[i];

  /* 10. Update state estimation error covariance matrix using (3.25) */
  //P_matrix_new = (eye(15) - K_matrix * H_matrix) * P_matrix_propagated;
  matmul_row("NN", 15, 15, 6, 1.0, K_matrix, H_matrix, 0.0, K_x_H);
  for (i = 0; i < 15; i++)
  {
    for (j = 0; j < 15; j++)
    {
      I_less_k_x_H[i * 15 + j] = I_par[i * 15 + j] - K_x_H[i * 15 + j];
    }
  }

  matmul_row("NN", 15, 15, 15, 1.0, I_less_k_x_H, P_matrix_propagated, 0.0, P_matrix_new);

  fprintf(out_KF_state_error, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
    %lf %lf %lf %lf %lf\n",
          time,
          x_est_new[0], x_est_new[1], x_est_new[2], x_est_new[3], x_est_new[4], x_est_new[5],
          x_est_new[6], x_est_new[7], x_est_new[8], x_est_new[9], x_est_new[10],
          x_est_new[11], x_est_new[12], x_est_new[13], x_est_new[14], x_est_new[15]);

  /* CLOSED-LOOP CORRECTION */

  /*
    for(i=0;i<no_meas;i++) {
       printf("delta_z[%d]:%lf\t",i,delta_z[i]);
    }
    printf("\n");

   for (i=0;i<15;i++) printf("x_est_new.[%d]:%lf\n",i, x_est_new[i]);
*/

  /* Correct attitude, velocity, and position using (14.7-9) */

  /* Quaternion attitude correction
    printf("est_C_b_e_old: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
    est_C_b_e_old[0],est_C_b_e_old[1],est_C_b_e_old[2],est_C_b_e_old[3],\
    est_C_b_e_old[4],est_C_b_e_old[5],est_C_b_e_old[6],est_C_b_e_old[7],est_C_b_e_old[8]);

    Quaternion_attitude_errror_correction(x_est_new, est_C_b_e_old, est_C_b_e_new);
   */

  Skew_symmetric(x_est_new, Skew_x_est_new);

  for (i = 0; i < 9; i++)
    I_less_Skew[i] = I[i] - Skew_x_est_new[i];

  matmul_row("NN", 3, 3, 3, 1.0, I_less_Skew, est_C_b_e_old, 0.0, est_C_b_e_new);

  for (i = 3; i < 6; i++)
    est_v_eb_e_new[i - 3] = est_v_eb_e_old[i - 3] - x_est_new[i];
  for (i = 6; i < 9; i++)
    est_r_eb_e_new[i - 6] = est_r_eb_e_old[i - 6] - x_est_new[i];

  /* Update IMU bias and GNSS receiver clock estimates */
  for (i = 0; i < 6; i++)
    est_IMU_bias_new[i] = est_IMU_bias_old[i] + x_est_new[i + 9];

  /* Ends */

  /* Freeing memory */
  free(Phi_matrix);
  free(Phi_transp);
  free(Q_prime_matrix);
  free(I);
  free(x_est_propagated);
  free(x_est_new);
  free(H_matrix_transp);
  free(P_matrix_propagated);
  free(P_aux);
  free(Q_aux);
  free(P_matrix);
  free(u_as_e_T);
  free(pred_meas);
  free(ones);
  free(I_par);
  free(K_matrix);
  free(delta_r);
  free(delta_z);
  free(Q_);
  free(F1);
  free(Q1);
  free(K_x_delta_z);
  free(I_meas);
  free(R_matrix);
  free(K_matrix_inv);
  free(K_x_H);
  free(I_less_k_x_H);
  free(H_matrix);
}

/* Name of function ------------------------------------------------------------
* Brief description
*Nav_equations_ECEF - Runs precision ECEF-frame inertial navigation
 %equations
 %
 % Software for use with "Principles of GNSS, Inertial, and Multisensor
 % Integrated Navigation Systems," Second Edition.
 %
 % This function created 1/4/2012 by Paul Groves
 %
 % Inputs:
 %   tor_i         time interval between epochs (s)
 %   old_r_eb_e    previous Cartesian position of body frame w.r.t. ECEF
 %                 frame, resolved along ECEF-frame axes (m)
 %   old_C_b_e     previous body-to-ECEF-frame coordinate transformation matrix
 %   old_v_eb_e    previous velocity of body frame w.r.t. ECEF frame, resolved
 %                 along ECEF-frame axes (m/s)
 %   f_ib_b        specific force of body frame w.r.t. ECEF frame, resolved
 %                 along body-frame axes, averaged over time interval (m/s^2)
 %   omega_ib_b    angular rate of body frame w.r.t. ECEF frame, resolved
 %                 about body-frame axes, averaged over time interval (rad/s)
 % Outputs:
 %   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
 %                 along ECEF-frame axes (m)
 %   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
 %                 ECEF-frame axes (m/s)
 %   C_b_e         body-to-ECEF-frame coordinate transformation matrix

 % Copyright 2012, Paul Groves
 % License: BSD; see license.txt for details
*-----------------------------------------------------------------------------*/
void Nav_equations_ECEF(double tor_i,
                        double *old_r_eb_e, double *old_v_eb_e, double *old_C_b_e,
                        double *f_ib_b, double *omega_ib_b, double *r_eb_e,
                        double *v_eb_e, double *C_b_e)
{

  /* parameters  */
  double omega_ie = 7.292115E-5; // Earth rotation rate (rad/s)
  double alpha_ie, mag_alpha;
  double C_Earth[9], alpha_ib_b[3], Alpha_ib_b[9];
  double Alpha_squared[9], second_term[9], first_term[9], C_new_old[9];
  double I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double C_aux[9], first_term_2[9], second_term_2[9], C_b_e_Cbb[9];
  double ave_C_b_e[9], Cbb[9], alpha_ie_vec[3], Alpha_ie[9], last_term[9];
  double f_ib_e[3], omega_ie_vec[3], Omega_ie[9], g[3], Omega_v_eb_e[3];
  double new_q_b_e[4], new_C_b_e[9];
  int i, j;

  printf("INPUT: \n");
  printf("tor_i= %f;\n", tor_i);
  printf("old_r_eb_e=[%lf; %lf; %lf];\n", old_r_eb_e[0], old_r_eb_e[1], old_r_eb_e[2]);
  printf("old_v_eb_e=[%lf; %lf; %lf];\n", old_v_eb_e[0], old_v_eb_e[1], old_v_eb_e[2]);
  printf("f_ib_b=[%lf; %lf; %lf];\n", f_ib_b[0], f_ib_b[1], f_ib_b[2]);
  printf("omega_ib_b=[%lf, %lf, %lf];\n", omega_ib_b[0], omega_ib_b[1], omega_ib_b[2]);
  printf("old_C_b_e=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", old_C_b_e[i * 3 + j]); /* code */
    }
    printf("];\n");
  }

  /* Begins     */

  /* ATTITUDE UPDATE  */
  /* From (2.145) determine the Earth rotation over the update interval
      C_Earth = C_e_i' * old_C_e_i  */
  alpha_ie = omega_ie * tor_i;
  C_Earth[0] = cos(alpha_ie);
  C_Earth[1] = sin(alpha_ie);
  C_Earth[2] = 0.0;
  C_Earth[3] = -sin(alpha_ie);
  C_Earth[4] = cos(alpha_ie);
  C_Earth[5] = 0.0;
  C_Earth[6] = 0.0;
  C_Earth[7] = 0.0;
  C_Earth[8] = 1.0;

  omega_ie_vec[0] = 0;
  omega_ie_vec[1] = 0;
  omega_ie_vec[2] = omega_ie;

  /* Calculate attitude increment, magnitude, and skew-symmetric matrix  */
  for (i = 0; i < 3; i++)
    alpha_ib_b[i] = omega_ib_b[i] * tor_i;
  mag_alpha = norm(alpha_ib_b, 3);
  Skew_symmetric(alpha_ib_b, Alpha_ib_b);

  printf("alpha_ib_b=[%lf, %lf, %lf];\n", alpha_ib_b[0], alpha_ib_b[1], alpha_ib_b[2]);
  printf("Alpha_ib_b=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", Alpha_ib_b[i * 3 + j]); /* code */
    }
    printf("];\n");
  }

  printf("Mag_alpha:%lf\n", mag_alpha);

  /* Attitude Update using Quaternion algebra */
  attitude_update(alpha_ib_b, omega_ie_vec, tor_i, old_C_b_e, new_q_b_e);
  Quaternion_to_DCM(new_q_b_e, new_C_b_e);

  printf("Quaternion new_C_b_e=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", new_C_b_e[i * 3 + j]); /* code */
    }
    printf("];\n");
  }

  /* Obtain coordinate transformation matrix from the new attitude w.r.t. an
      inertial frame to the old using Rodrigues' formula, (5.73)  */
  matmul_row("NN", 3, 3, 3, 1.0, Alpha_ib_b, Alpha_ib_b, 0.0, Alpha_squared);

  for (i = 0; i < 9; i++)
    second_term[i] = (1 - cos(mag_alpha)) /
                     (mag_alpha * mag_alpha) * Alpha_squared[i];
  for (i = 0; i < 9; i++)
    first_term[i] = I[i] + sin(mag_alpha) /
                               mag_alpha * Alpha_ib_b[i];

  if (mag_alpha > 1.E-8)
  {
    for (i = 0; i < 9; i++)
      C_new_old[i] = first_term[i] + second_term[i];
  }
  else
  {
    for (i = 0; i < 9; i++)
      C_new_old[i] = I[i] + Alpha_ib_b[i];
  } // end if mag_alpha

  printf("C_new_old=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", C_new_old[i * 3 + j]); /* code */
    }
    printf("];\n");
  }
  printf("C_Earth=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", C_Earth[i * 3 + j]); /* code */
    }
    printf("];\n");
  }

  /* Update attitude using (5.75)  */
  matmul_row("NN", 3, 3, 3, 1.0, C_Earth, old_C_b_e, 0.0, C_aux);
  //matmul("NN", 3, 3, 3, 1.0, old_C_b_e, C_Earth, 0.0, C_aux);

  matmul_row("NN", 3, 3, 3, 1.0, C_aux, C_new_old, 0.0, C_b_e);
  //matmul("NN", 3, 3, 3, 1.0, C_new_old, C_aux, 0.0, C_b_e);

  /* SPECIFIC FORCE FRAME TRANSFORMATION
     % Calculate the average body-to-ECEF-frame coordinate transformation
     % matrix over the update interval using (5.84) and (5.85)  */
  for (i = 0; i < 9; i++)
    first_term_2[i] = I[i] + second_term[i];
  for (i = 0; i < 9; i++)
    second_term_2[i] = ((1 - (sin(mag_alpha) / mag_alpha)) /
                        (mag_alpha * mag_alpha)) *
                       Alpha_squared[i];
  for (i = 0; i < 9; i++)
    Cbb[i] = first_term_2[i] + second_term_2[i];
  alpha_ie_vec[0] = 0;
  alpha_ie_vec[1] = 0;
  alpha_ie_vec[2] = alpha_ie;
  Skew_symmetric(alpha_ie_vec, Alpha_ie);
  matmul_row("NN", 3, 3, 3, 0.5, Alpha_ie, old_C_b_e, 0.0, last_term);
  //matmul("NN", 3, 3, 3, 0.5, old_C_b_e, Alpha_ie, 0.0, last_term);

  if (mag_alpha > 1.E-8)
  {
    matmul_row("NN", 3, 3, 3, 1.0, old_C_b_e, Cbb, 0.0, C_b_e_Cbb);
    //matmul("NN", 3, 3, 3, 1.0, Cbb, old_C_b_e, 0.0, C_b_e_Cbb);
    for (i = 0; i < 9; i++)
      ave_C_b_e[i] = C_b_e_Cbb[i] - last_term[i];
  }
  else
  {
    for (i = 0; i < 9; i++)
      ave_C_b_e[i] = old_C_b_e[i] - last_term[i];
  } //if mag_alpha

  printf("ave_C_b_e=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", ave_C_b_e[i * 3 + j]); /* code */
    }
    printf("];\n");
  }

  /* Transform specific force to ECEF-frame resolving axes using (5.85) */
  matmul_row("NN", 3, 1, 3, 1.0, ave_C_b_e, f_ib_b, 0.0, f_ib_e);
  //matmul("NN", 1, 3, 3, 1.0, f_ib_b, ave_C_b_e, 0.0, f_ib_e);

  printf("f_ib_e=[%lf, %lf, %lf];\n", f_ib_e[0], f_ib_e[1], f_ib_e[2]);

  /* UPDATE VELOCITY
     % From (5.36), */
  Skew_symmetric(omega_ie_vec, Omega_ie);
  printf("Omega_ie=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", Omega_ie[i * 3 + j]); /* code */
    }
    printf("];\n");
  }
  Gravity_ECEF(old_r_eb_e, g);
  printf("g=[%lf, %lf, %lf];\n", g[0], g[1], g[2]);

  matmul_row("NN", 3, 1, 3, 1.0, Omega_ie, old_v_eb_e, 0.0, Omega_v_eb_e);
  //matmul("NN", 1, 3, 3, 1.0, old_v_eb_e, Omega_ie, 0.0, Omega_v_eb_e);
  for (i = 0; i < 3; i++)
    v_eb_e[i] = old_v_eb_e[i] + tor_i * (f_ib_e[i] + g[i] -
                                         2 * Omega_v_eb_e[i]);

  /* UPDATE CARTESIAN POSITION
     % From (5.38), */
  for (i = 0; i < 3; i++)
    r_eb_e[i] = old_r_eb_e[i] + (v_eb_e[i] + old_v_eb_e[i]) * 0.5 * tor_i;

  printf("OUTPUT: \n");
  printf("P: %lf, %lf, %lf\n", r_eb_e[0], r_eb_e[1], r_eb_e[2]);
  printf("V: %lf, %lf, %lf\n", v_eb_e[0], v_eb_e[1], v_eb_e[2]);
  printf("C_b_e\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf ", C_b_e[i * 3 + j]); /* code */
    }
    printf("\n");
  }
}

/* Propagate Inertial covariance matrix - modify P_matrix_old for Tightly coupled */
propinsstateTC(double *P_matrix_old, const double *est_r_eb_e, const double *est_C_b_e_old,
               const double *meas_f_ib_b, TC_KF_config *TC_KF_config, double tor_s)
{
  double *Phi_matrix, *Phi_transp;
  double *Q_prime_matrix, *Q_, *Q_aux;
  double *P_matrix_propagated, *P_aux;
  double *I, omega_ie_vec[3], Omega_ie[9] = {0.0};
  double meas_f_ib_e[3] = {0.0}, Skew_meas_f_ib_e[9] = {0.0};
  double est_L_b_old = 0.0, llh[3] = {0.0}, geocentric_radius, g[3] = {0.0}, g_est_r_eb_e[3] = {0.0};
  int i, j, n;

  /* Constants */
  double e = sqrt(e_2); /*WGS84 eccentricity                        */

  n = 17; /*tightly*/

  ecef2pos(est_r_eb_e, llh);
  est_L_b_old = llh[0];

  printf("LATITUDE: %lf\n", est_L_b_old * R2D);

  /* Initialize matrices and vectors */
  Phi_matrix = eye(n);
  Phi_transp = zeros(n, n);
  Q_prime_matrix = zeros(n, n);
  P_matrix_propagated = mat(n, n);
  P_aux = mat(n, n);
  Q_aux = mat(n, n);
  Q_ = mat(n, n);
  I = eye(3);

  printf("est_r_eb_e: %lf, %lf, %lf\n", est_r_eb_e[0], est_r_eb_e[1], est_r_eb_e[2]);
  printf("I1 and tor_i: %lf\n", tor_s);
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%.13lf ", I[i * 3 + j] * tor_s);
    }
    printf("\n");
  }
  printf("\n");

  /* Skew symmetric matrix of Earth rate */
  omega_ie_vec[0] = 0.0;
  omega_ie_vec[1] = 0.0;
  omega_ie_vec[2] = OMGE;
  Skew_symmetric(omega_ie_vec, Omega_ie);

  /* SYSTEM PROPAGATION PHASE */

  /* 1. Determine transition matrix using (14.50) (first-order approx) */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * n + j] = Phi_matrix[i * n + j] - Omega_ie[i * 3 + j] * tor_s;
    }
  }

  for (i = 0; i < 3; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Phi_matrix[i * n + j] = est_C_b_e_old[i * 3 + (j - 12)] * tor_s;
    }
  }

  matmul_row("NN", 3, 1, 3, 1.0, est_C_b_e_old, meas_f_ib_b, 0.0, meas_f_ib_e);
  Skew_symmetric(meas_f_ib_e, Skew_meas_f_ib_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * n + j] = -tor_s * Skew_meas_f_ib_e[(i - 3) * 3 + j];
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * n + j] = Phi_matrix[i * n + j] - 2 * Omega_ie[(i - 3) * 3 + (j - 3)] * tor_s;
    }
  }

  geocentric_radius = RE_WGS84 / sqrt(1 - pow((e * sin(est_L_b_old)), 2)) *
                      sqrt(pow(cos(est_L_b_old), 2) + pow((1 - e * e), 2) * pow(sin(est_L_b_old), 2)); /* from (2.137)*/
  Gravity_ECEF(est_r_eb_e, g);                                                                         //returns a vector
  matmul_row("NN", 3, 3, 1, 1.0, g, est_r_eb_e, 0.0, g_est_r_eb_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 6; j < 9; j++)
    {
      Phi_matrix[i * n + j] = -tor_s * 2 /
                              geocentric_radius * g_est_r_eb_e[(i - 3) * 3 + (j - 6)] / (norm(est_r_eb_e, 3));
    }
  }

  printf("Phi_matrix1a\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", Phi_matrix[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  for (i = 3; i < 6; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Phi_matrix[i * n + j] = est_C_b_e_old[(i - 3) * 3 + (j - 9)] * tor_s;
    }
  }

  for (i = 6; i < 9; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * n + j] = I[(i - 6) * 3 + (j - 3)] * tor_s;
    }
  }

  Phi_matrix[15 * n + 16] = tor_s;

  printf("Phi_matrix1\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", Phi_matrix[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 2. Determine approximate system noise covariance matrix using (14.82) */
  //Q_prime_matrix(1:3,1:3) = eye(3) * TC_KF_config.gyro_noise_PSD * tor_s;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i)*3 + (j)] *
                                                (TC_KF_config->gyro_noise_PSD) * tor_s
                                          : 0.0);
    }
  }
  //Q_prime_matrix(4:6,4:6) = eye(3) * TC_KF_config.accel_noise_PSD * tor_s;
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i - 3) * 3 + (j - 3)] *
                                                (TC_KF_config->accel_noise_PSD) * tor_s
                                          : 0.0);
    }
  }
  //Q_prime_matrix(10:12,10:12) = eye(3) * TC_KF_config.accel_bias_PSD * tor_s;
  for (i = 9; i < 12; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i - 9) * 3 + (j - 9)] *
                                                (TC_KF_config->accel_bias_PSD) * tor_s
                                          : 0.0);
    }
  }

  //Q_prime_matrix(13:15,13:15) = eye(3) * TC_KF_config.gyro_bias_PSD * tor_s;
  for (i = 12; i < 15; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i - 12) * 3 + (j - 12)] *
                                                (TC_KF_config->gyro_bias_PSD) * tor_s
                                          : 0.0);
    }
  }

  //Q_prime_matrix(16,16) = TC_KF_config.clock_phase_PSD * tor_s;
  Q_prime_matrix[15 * n + 15] = TC_KF_config->clock_phase_PSD * tor_s;

  //Q_prime_matrix(n,n) = TC_KF_config->clock_freq_PSD * tor_s;
  Q_prime_matrix[16 * n + 16] = TC_KF_config->clock_freq_PSD * tor_s;

  printf("Q_matrix1\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", Q_prime_matrix[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 4. Propagate state estimation error covariance matrix using (3.46) */
  //P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *\
    //    Phi_matrix' + 0.5 * Q_prime_matrix;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      P_aux[i * n + j] = P_matrix_old[i * n + j] + 0.5 * Q_prime_matrix[i * n + j];
    }
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      Phi_transp[j * n + i] = Phi_matrix[i * n + j];
    }
  }

  matmul_row("NN", n, n, n, 1.0, P_aux, Phi_transp, 0.0, Q_aux); /*(Pp + 0.5*Q)*PHI^T */
  //matmul("TN",n,n,n,1.0,Phi_matrix,P_aux,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
  matmul_row("NN", n, n, n, 1.0, Phi_matrix, Q_aux, 0.0, Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */
  /*matmul("NN",n,n,n,1.0,Q_aux,Phi_matrix,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      P_matrix_propagated[i * n + j] = Q_[i * n + j] + 0.5 * Q_prime_matrix[i * n + j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
    }
  }

  printf("P_matrix_old1\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", P_matrix_old[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("P_matrix_propagated1\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", P_matrix_propagated[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* Propagating P matrix */
  for (i = 0; i < n * n; i++)
    P_matrix_old[i] = P_matrix_propagated[i];

  free(Phi_matrix);
  free(Phi_transp);
  free(Q_);
  free(Q_aux);
  free(Q_prime_matrix);
  free(P_matrix_propagated);
  free(P_aux);
  free(I);
}

/* Propagate Inertial covariance matrix - modify P_matrix_old for Loosley coupled */
propinsstateLC(double *P_matrix_old, const double *est_r_eb_e, const double *est_C_b_e_old,
               const double *meas_f_ib_b, LC_KF_config *LC_KF_config, double tor_s)
{
  double *Phi_matrix, *Phi_transp;
  double *Q_prime_matrix, *Q_, *Q_aux;
  double *P_matrix_propagated, *P_aux;
  double I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}, omega_ie_vec[3], Omega_ie[9] = {0.0};
  double meas_f_ib_e[3] = {0.0}, Skew_meas_f_ib_e[9] = {0.0};
  double est_L_b_old = 0.0, llh[3] = {0.0}, geocentric_radius, g[3] = {0.0}, g_est_r_eb_e[3] = {0.0};
  int i, j, n;

  /* Constants */
  double e = sqrt(e_2); /*WGS84 eccentricity                        */

  n = 15; /*loosley*/

  ecef2pos(est_r_eb_e, llh);
  est_L_b_old = llh[0];

  /* Initialize matrices and vectors */
  Phi_matrix = eye(n);
  Phi_transp = zeros(n, n);
  Q_prime_matrix = zeros(n, n);
  P_matrix_propagated = mat(n, n);
  P_aux = mat(n, n);
  Q_aux = mat(n, n);
  Q_ = mat(n, n);

  /* Skew symmetric matrix of Earth rate */
  omega_ie_vec[0] = 0.0;
  omega_ie_vec[1] = 0.0;
  omega_ie_vec[2] = OMGE;
  Skew_symmetric(omega_ie_vec, Omega_ie);

  /* SYSTEM PROPAGATION PHASE */

  /* 1. Determine transition matrix using (14.50) (first-order approx) */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * n + j] = Phi_matrix[i * n + j] - Omega_ie[i * 3 + j] * tor_s;
    }
  }

  for (i = 0; i < 3; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Phi_matrix[i * n + j] = est_C_b_e_old[i * 3 + (j - 12)] * tor_s;
    }
  }

  matmul_row("NN", 3, 1, 3, 1.0, est_C_b_e_old, meas_f_ib_b, 0.0, meas_f_ib_e);
  Skew_symmetric(meas_f_ib_e, Skew_meas_f_ib_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Phi_matrix[i * n + j] = -tor_s * Skew_meas_f_ib_e[(i - 3) * 3 + j];
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * n + j] = Phi_matrix[i * n + j] - 2 * Omega_ie[(i - 3) * 3 + (j - 3)] * tor_s;
    }
  }

  geocentric_radius = RE_WGS84 / sqrt(1 - pow((e * sin(est_L_b_old)), 2)) *
                      sqrt(pow(cos(est_L_b_old), 2) + pow((1 - e * e), 2) * pow(sin(est_L_b_old), 2)); /* from (2.137)*/
  Gravity_ECEF(est_r_eb_e, g);                                                                         //returns a vector
  matmul_row("NN", 3, 3, 1, 1.0, g, est_r_eb_e, 0.0, g_est_r_eb_e);

  for (i = 3; i < 6; i++)
  {
    for (j = 6; j < 9; j++)
    {
      Phi_matrix[i * n + j] = -tor_s * 2 /
                              geocentric_radius * g_est_r_eb_e[(i - 3) * 3 + (j - 6)] / (norm(est_r_eb_e, 3));
    }
  }

  for (i = 3; i < 6; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Phi_matrix[i * n + j] = est_C_b_e_old[(i - 3) * 3 + (j - 9)] * tor_s;
    }
  }

  for (i = 6; i < 9; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Phi_matrix[i * n + j] = I[(i - 6) * 3 + (j - 3)] * tor_s;
    }
  }

  printf("Phi_matrix\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", Phi_matrix[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 2. Determine approximate system noise covariance matrix using (14.82) */
  //Q_prime_matrix(1:3,1:3) = eye(3) * TC_KF_config.gyro_noise_PSD * tor_s;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i)*3 + (j)] *
                                                (LC_KF_config->gyro_noise_PSD) * tor_s
                                          : 0.0);
    }
  }
  //Q_prime_matrix(4:6,4:6) = eye(3) * TC_KF_config.accel_noise_PSD * tor_s;
  for (i = 3; i < 6; i++)
  {
    for (j = 3; j < 6; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i - 3) * 3 + (j - 3)] *
                                                (LC_KF_config->accel_noise_PSD) * tor_s
                                          : 0.0);
    }
  }
  //Q_prime_matrix(10:12,10:12) = eye(3) * TC_KF_config.accel_bias_PSD * tor_s;
  for (i = 9; i < 12; i++)
  {
    for (j = 9; j < 12; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i - 9) * 3 + (j - 9)] *
                                                (LC_KF_config->accel_bias_PSD) * tor_s
                                          : 0.0);
    }
  }

  //Q_prime_matrix(13:15,13:15) = eye(3) * TC_KF_config.gyro_bias_PSD * tor_s;
  for (i = 12; i < 15; i++)
  {
    for (j = 12; j < 15; j++)
    {
      Q_prime_matrix[i * n + j] = (i == j ? I[(i - 12) * 3 + (j - 12)] *
                                                (LC_KF_config->gyro_bias_PSD) * tor_s
                                          : 0.0);
    }
  }

  printf("Q_matrix\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", Q_prime_matrix[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* 4. Propagate state estimation error covariance matrix using (3.46) */
  //P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *\
    //    Phi_matrix' + 0.5 * Q_prime_matrix;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      P_aux[i * n + j] = P_matrix_old[i * n + j] + 0.5 * Q_prime_matrix[i * n + j];
    }
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      Phi_transp[j * n + i] = Phi_matrix[i * n + j];
    }
  }

  matmul_row("NN", n, n, n, 1.0, P_aux, Phi_transp, 0.0, Q_aux); /*(Pp + 0.5*Q)*PHI^T */
  //matmul("TN",n,n,n,1.0,Phi_matrix,P_aux,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
  matmul_row("NN", n, n, n, 1.0, Phi_matrix, Q_aux, 0.0, Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */
  /*matmul("NN",n,n,n,1.0,Q_aux,Phi_matrix,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      P_matrix_propagated[i * n + j] = Q_[i * n + j] + 0.5 * Q_prime_matrix[i * n + j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
    }
  }

  printf("P_matrix_old\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", P_matrix_old[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("P_matrix_propagated\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%.13lf ", P_matrix_propagated[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  /* Propagating P matrix */
  for (i = 0; i < n * n; i++)
    P_matrix_old[i] = P_matrix_propagated[i];

  free(Phi_matrix);
  free(Phi_transp);
  free(Q_);
  free(Q_aux);
  free(Q_prime_matrix);
  free(P_matrix_propagated);
  free(P_aux);
}

/* Name of function ------------------------------------------------------------
* Brief description
* % Inputs:
%   no_epochs    Number of epochs of profile data
%   initialization_errors
%     .delta_r_eb_n     position error resolved along NED (m)
%     .delta_v_eb_n     velocity error resolved along NED (m/s)
%     .delta_eul_nb_n   attitude error as NED Euler angles (rad)
%   IMU_errors
%     .delta_r_eb_n     position error resolved along NED (m)
%     .b_a              Accelerometer biases (m/s^2)
%     .b_g              Gyro biases (rad/s)
%     .M_a              Accelerometer scale factor and cross coupling errors
%     .M_g              Gyro scale factor and cross coupling errors
%     .G_g              Gyro g-dependent biases (rad-sec/m)
%     .accel_noise_root_PSD   Accelerometer noise root PSD (m s^-1.5)
%     .gyro_noise_root_PSD    Gyro noise root PSD (rad s^-0.5)
%     .accel_quant_level      Accelerometer quantization level (m/s^2)
%     .gyro_quant_level       Gyro quantization level (rad/s)
%   GNSS_config
%     .epoch_interval     Interval between GNSS epochs (s)
%     .init_est_r_ea_e    Initial estimated position (m; ECEF)
%     .no_sat             Number of satellites in constellation
%     .r_os               Orbital radius of satellites (m)
%     .inclination        Inclination angle of satellites (deg)
%     .const_delta_lambda Longitude offset of constellation (deg)
%     .const_delta_t      Timing offset of constellation (s)
%     .mask_angle         Mask angle (deg)
%     .SIS_err_SD         Signal in space error SD (m)
%     .zenith_iono_err_SD Zenith ionosphere error SD (m)
%     .zenith_trop_err_SD Zenith troposphere error SD (m)
%     .code_track_err_SD  Code tracking error SD (m)
%     .rate_track_err_SD  Range rate tracking error SD (m/s)
%     .rx_clock_offset    Receiver clock offset at time=0 (m)
%     .rx_clock_drift     Receiver clock drift at time=0 (m/s)
%   LC_KF_config
%     .init_att_unc           Initial attitude uncertainty per axis (rad)
%     .init_vel_unc           Initial velocity uncertainty per axis (m/s)
%     .init_pos_unc           Initial position uncertainty per axis (m)
%     .init_b_a_unc           Initial accel. bias uncertainty (m/s^2)
%     .init_b_g_unc           Initial gyro. bias uncertainty (rad/s)
%     .gyro_noise_PSD         Gyro noise PSD (rad^2/s)
%     .accel_noise_PSD        Accelerometer noise PSD (m^2 s^-3)
%     .accel_bias_PSD         Accelerometer bias random walk PSD (m^2 s^-5)
%     .gyro_bias_PSD          Gyro bias random walk PSD (rad^2 s^-3)
%     .pos_meas_SD            Position measurement noise SD per axis (m)
%     .vel_meas_SD            Velocity measurement noise SD per axis (m/s)
%
%  Outputs:
%   out_profile        Navigation solution as a motion profile array
%   out_errors         Navigation solution error array
%   out_IMU_bias_est   Kalman filter IMU bias estimate array
%   out_clock          GNSS Receiver clock estimate array
%   out_KF_SD          Output Kalman filter state uncertainties
% notes  :
*-----------------------------------------------------------------------------*/
void Loosely_coupled_INS_GNSS(INS_measurements *INS_measurements, int no_meas,
                              PVAT_solution *pvat_old, int no_par, double old_time, double time_last_GNSS, double curr_GNSS_time,
                              double *GNSS_r_eb_e, double *GNSS_v_eb_e,
                              initialization_errors *initialization_errors, IMU_errors *IMU_errors,
                              GNSS_config *GNSS_config, LC_KF_config *LC_KF_config, PVAT_solution *pvat_new,
                              double *out_errors, double *out_IMU_bias_est, double *out_KF_SD,
                              double *meas_f_ib_b, double *meas_omega_ib_b)
{
  double true_L_b, true_lambda_b, true_h_b;
  double old_est_llh[3], old_est_r_eb_e[3], old_est_v_eb_e[3], old_true_C_b_e[9];
  double est_clock[2], est_r_eb_e[3];
  double old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n[3];
  double true_v_eb_n[3], true_C_b_n[9], true_eul_nb[3], old_est_C_b_n[9];
  double old_true_r_eb_e[3], old_true_v_eb_e[3];
  double est_L_b, est_h_b, est_lambda_b, est_v_eb_e[3], est_v_eb_n[3];
  double est_C_b_e[9], est_C_b_n[9];
  double delta_r_eb_n[3], delta_v_eb_n[3], delta_eul_nb_n[3];
  double est_IMU_bias[6] = {0};
  double quant_residuals[6] = {0};
  double *P_matrix, *P_matrix_new;
  double GNSS_epoch, tor_i, time, tor_s;
  int i, j;
  double est_C_b_e_new[9], est_C_b_n_T[9], est_v_eb_e_new[3], est_r_eb_e_new[3];
  double est_IMU_bias_new[6], est_clock_new[2];
  double q_nb[4], llh[3];
  double C_Transp[9] = {0.0};

  printf("\n *****************  LC_INS/GNSS BEGINS ************************\n");
  /**/
  printf("P: %lf, %lf, %lf\n", pvat_old->latitude * R2D, pvat_old->longitude * R2D, pvat_old->height);
  printf("V: %lf, %lf, %lf\n", pvat_old->ned_velocity[0], pvat_old->ned_velocity[1], pvat_old->ned_velocity[2]);
  printf("A: %lf, %lf, %lf\n", pvat_old->euler_angles[0], pvat_old->euler_angles[1], pvat_old->euler_angles[2]);

  /* Initialize true navigation solution */
  time = INS_measurements->sec;
  est_L_b = true_L_b = pvat_old->latitude;
  true_lambda_b = pvat_old->longitude;
  true_h_b = pvat_old->height;
  for (i = 0; i < 3; i++)
    true_v_eb_n[i] = pvat_old->ned_velocity[i];
  for (i = 0; i < 3; i++)
    true_eul_nb[i] = pvat_old->euler_angles[i];
  Euler_to_CTM(true_eul_nb, true_C_b_n);

  printf("true_eul_nb: %lf, %lf, %lf\n", true_eul_nb[0], true_eul_nb[1], true_eul_nb[2]);

  /* Transposing "true_C_b_n" (in _n_b) to actual _b_n frame*/
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      C_Transp[j * 3 + i] = true_C_b_n[i * 3 + j];
    }
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      true_C_b_n[i * 3 + j] = C_Transp[i * 3 + j];
    }
  }

  printf("true_C_b_n from euler\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%.10lf ", true_C_b_n[i * 3 + j]);
    }
    printf("\n");
  }

  NED_to_ECEF(&true_L_b, &true_lambda_b, &true_h_b, true_v_eb_n, true_C_b_n,
              old_true_r_eb_e, old_est_v_eb_e, old_true_C_b_e);

  printf("old_true_C_b_e:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",
         old_true_C_b_e[0], old_true_C_b_e[1], old_true_C_b_e[2],
         old_true_C_b_e[3], old_true_C_b_e[4], old_true_C_b_e[5],
         old_true_C_b_e[6], old_true_C_b_e[7], old_true_C_b_e[8]);

  /* Initialize other matrices */

  /* Initialize Kalman filter P matrix and IMU bias states */
  P_matrix = zeros(no_par, no_par); // in this case 15x15
  P_matrix_new = zeros(no_par, no_par);
  if (out_errors[0] > 0.0 && out_errors[1] > 0.0 && out_errors[2] > 0.0)
  {
    printf("ERRROS ARE NOT ZERO \n");
    for (i = 0; i < no_par; i++)
    {
      for (j = 0; j < no_par; j++)
      {
        P_matrix[i * no_par + j] = pvat_old->P[i * no_par + j];
      }
    }
  }
  else
  {
    printf("INITIALIZING P_matrix TO DEFAULT\n");
    Initialize_LC_P_matrix(LC_KF_config, P_matrix); //IT WILL CHANGE FOR THE UNDIFF/UNCOMB MODEL
  }

  /* Generate IMU bias and clock output records
  out_IMU_bias_est(1,1) = old_time;
  out_IMU_bias_est(1,2:7) = est_IMU_bias';
  out_clock(1,1) = old_time;
  out_clock(1,2:3) = est_clock;  */

  /* Generate KF uncertainty record */

  /* Initialize GNSS model timing */
  GNSS_epoch = 1;

  /* Main loop */

  //for (i = 0; i < no_obs; i++) {
  /* Time interval */
  tor_i = time - old_time;

  if (old_time <= 0.000)
  {
    tor_i = 1.0; //1 second of time interval
    time_last_GNSS = time - tor_i;
  }

  /* Correct IMU errors */
  for (i = 0; i < 6; i++)
    est_IMU_bias[i] = out_IMU_bias_est[i];
  for (i = 0; i < 3; i++)
    meas_f_ib_b[i] = INS_measurements->f_ib_b[i];
  for (i = 0; i < 3; i++)
    meas_omega_ib_b[i] = INS_measurements->omega_ib_b[i];
  for (i = 0; i < 3; i++)
    meas_f_ib_b[i] = meas_f_ib_b[i] - est_IMU_bias[i];
  for (i = 0; i < 3; i++)
    meas_omega_ib_b[i] = meas_omega_ib_b[i] - est_IMU_bias[i + 3];

  /* Update estimated navigation solution */
  printf("\n *****************  NAV_EQUATIONS BEGINS ************************\n");
  Nav_equations_ECEF(tor_i,
                     old_true_r_eb_e, old_est_v_eb_e, old_true_C_b_e, meas_f_ib_b,
                     meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
  printf("\n *****************  NAV_EQUATIONS ENDS **************************\n");

  /* INS covariance propagation */
  propinsstateLC(P_matrix, est_r_eb_e, est_C_b_e, meas_f_ib_b, LC_KF_config, tor_i);

  /* Determine whether to update GNSS simulation and run Kalman filter */

  /* Only process if GNSS and INS timne differences are within 0.4 seconds */
  printf("TimeDIFF: %lf, no_meas: %d, NORM hor.pos.: %lf - time: %lf, time_last_GNSS: %lf\n",
         fabs(curr_GNSS_time - INS_measurements->sec),
         no_meas, norm(LC_KF_config->init_pos_unc_ned, 2), time, time_last_GNSS);

  if (fabs(curr_GNSS_time - INS_measurements->sec) < 0.0001 && norm(LC_KF_config->init_pos_unc_ned, 2) < 3.0)
  {

    printf("LCKF_here\n");
    /* KF time interval */
    tor_s = time - time_last_GNSS;

    /* Run Integration Kalman filter */
    printf("\n *******************  LC_KF_EPOCH BEGINS ************************\n");
    /**/ LC_KF_Epoch(time, GNSS_r_eb_e, GNSS_v_eb_e, no_meas, tor_s, est_C_b_e,
                     est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix,
                     meas_f_ib_b, est_L_b, LC_KF_config, est_C_b_e_new,
                     est_v_eb_e_new, est_r_eb_e_new, est_IMU_bias_new,
                     P_matrix_new);
    printf("\n *******************  LC_KF_EPOCH ENDS **************************\n");
    /*
printf("\n *****************  NAV_EQUATIONS CLOSED-LOOP BEGINS **************\n");

/* Correct IMU errors
for (i=0;i<3;i++) meas_f_ib_b[i]=INS_measurements->f_ib_b[i];
for (i=0;i<3;i++) meas_omega_ib_b[i]=INS_measurements->omega_ib_b[i];
for (i=0;i<3;i++) meas_f_ib_b[i] = meas_f_ib_b[i] - est_IMU_bias_new[i];
for (i=0;i<3;i++) meas_omega_ib_b[i] = meas_omega_ib_b[i] - est_IMU_bias_new[i+3];

for (i=0;i<3;i++) old_est_r_eb_e[i]=est_r_eb_e_new[i];
for (i=0;i<3;i++) old_est_v_eb_e[i]=est_v_eb_e_new[i];
for (i=0;i<9;i++) old_true_C_b_e[i]=est_C_b_e_new[i];

    Nav_equations_ECEF(tor_i,\
        old_est_r_eb_e, old_est_v_eb_e,old_true_C_b_e, meas_f_ib_b,\
        meas_omega_ib_b, est_r_eb_e_new, est_v_eb_e_new, est_C_b_e_new);
printf("\n *****************  NAV_EQUATIONS CLOSED-LOOP ENDS ****************\n");
*/
    time_last_GNSS = time;

    for (i = 0; i < 6; i++)
      out_IMU_bias_est[i] = est_IMU_bias_new[i];

    fprintf(out_IMU_bias_file, "%lf %lf %lf %lf %lf %lf %lf\n", time,
            est_IMU_bias_new[0], est_IMU_bias_new[1], est_IMU_bias_new[2],
            est_IMU_bias_new[3], est_IMU_bias_new[4], est_IMU_bias_new[5]);

    /* Generate KF uncertainty output record */
    fprintf(out_KF_SD_file, "%f\t", time);
    for (i = 0; i < no_par; i++)
    {
      for (j = 0; j < no_par; j++)
      {
        if (i == j)
        {
          out_errors[j] = P_matrix_new[i * no_par + j];
          fprintf(out_KF_SD_file, "%lf\t", sqrt(P_matrix_new[i * no_par + j]));
        }
      }
    }
    fprintf(out_KF_SD_file, "\n");

    /* Full weight matrix */
    for (i = 0; i < 15; i++)
    {
      for (j = 0; j < 15; j++)
      {
        pvat_new->P[i * 15 + j] = P_matrix_new[i * 15 + j];
      }
    }

    pvat_new->Nav_or_KF = 1;
  }
  else
  {
    /* Full weight matrix */
    for (i = 0; i < 15; i++)
    {
      for (j = 0; j < 15; j++)
      {
        pvat_new->P[i * 15 + j] = P_matrix[i * 15 + j];
      }
    }
    for (i = 0; i < 6; i++)
      out_IMU_bias_est[i] = est_IMU_bias[i];
    for (i = 0; i < 3; i++)
      est_r_eb_e_new[i] = est_r_eb_e[i];
    for (i = 0; i < 3; i++)
      est_v_eb_e_new[i] = est_v_eb_e[i];
    for (i = 0; i < 9; i++)
      est_C_b_e_new[i] = est_C_b_e[i];
    pvat_new->Nav_or_KF = 0;
  }

  /* Convert navigation solution to NED  */
  ECEF_to_NED(est_r_eb_e_new, est_v_eb_e_new, est_C_b_e_new,
              &est_L_b, &est_lambda_b, &est_h_b,
              est_v_eb_n, est_C_b_n);

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      C_Transp[j * 3 + i] = est_C_b_n[i * 3 + j];
    }
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      est_C_b_n_T[i * 3 + j] = C_Transp[i * 3 + j];
    }
  }

  double q_b_n[4];
  DCM_to_quaternion(est_C_b_n, q_b_n);
  Quaternion_to_euler(q_b_n, pvat_new->euler_angles);

  printf("Euler from Quaternions: %lf, %lf, %lf\n", pvat_new->euler_angles[0],
         pvat_new->euler_angles[1], pvat_new->euler_angles[2]);

  pvat_new->latitude = est_L_b;
  pvat_new->longitude = est_lambda_b;
  pvat_new->height = est_h_b;
  for (i = 0; i < 3; i++)
    pvat_new->ned_velocity[i] = est_v_eb_n[i];
  CTM_to_Euler(pvat_new->euler_angles, est_C_b_n_T);
  pvat_new->time = time;

  printf("Euler from CTM: %lf, %lf, %lf\n", pvat_new->euler_angles[0],
         pvat_new->euler_angles[1], pvat_new->euler_angles[2]);

  /* Passing Cbe from KF direclty to the next Navigation */
  for (i = 0; i < 9; i++)
    pvat_new->C_b_e[i] = est_C_b_e_new[i];
  for (i = 0; i < 3; i++)
    pvat_new->re[i] = est_r_eb_e_new[i];
  for (i = 0; i < 3; i++)
    pvat_new->ve[i] = est_v_eb_e_new[i];

  /* Free memory */
  free(P_matrix);
  free(P_matrix_new);
  printf("\n *****************  LC_INS/GNSS ENDS **************************\n");
}

/* Name of function ------------------------------------------------------------
* Brief description
* % Inputs:
%   no_epochs    Number of epochs of profile data
%   initialization_errors
%     .delta_r_eb_n     position error resolved along NED (m)
%     .delta_v_eb_n     velocity error resolved along NED (m/s)
%     .delta_eul_nb_n   attitude error as NED Euler angles (rad)
%   IMU_errors
%     .delta_r_eb_n     position error resolved along NED (m)
%     .b_a              Accelerometer biases (m/s^2)
%     .b_g              Gyro biases (rad/s)
%     .M_a              Accelerometer scale factor and cross coupling errors
%     .M_g              Gyro scale factor and cross coupling errors
%     .G_g              Gyro g-dependent biases (rad-sec/m)
%     .accel_noise_root_PSD   Accelerometer noise root PSD (m s^-1.5)
%     .gyro_noise_root_PSD    Gyro noise root PSD (rad s^-0.5)
%     .accel_quant_level      Accelerometer quantization level (m/s^2)
%     .gyro_quant_level       Gyro quantization level (rad/s)
%   GNSS_config
%     .epoch_interval     Interval between GNSS epochs (s)
%     .init_est_r_ea_e    Initial estimated position (m; ECEF)
%     .no_sat             Number of satellites in constellation
%     .r_os               Orbital radius of satellites (m)
%     .inclination        Inclination angle of satellites (deg)
%     .const_delta_lambda Longitude offset of constellation (deg)
%     .const_delta_t      Timing offset of constellation (s)
%     .mask_angle         Mask angle (deg)
%     .SIS_err_SD         Signal in space error SD (m)
%     .zenith_iono_err_SD Zenith ionosphere error SD (m)
%     .zenith_trop_err_SD Zenith troposphere error SD (m)
%     .code_track_err_SD  Code tracking error SD (m)
%     .rate_track_err_SD  Range rate tracking error SD (m/s)
%     .rx_clock_offset    Receiver clock offset at time=0 (m)
%     .rx_clock_drift     Receiver clock drift at time=0 (m/s)
%   TC_KF_config
%     .init_att_unc           Initial attitude uncertainty per axis (rad)
%     .init_vel_unc           Initial velocity uncertainty per axis (m/s)
%     .init_pos_unc           Initial position uncertainty per axis (m)
%     .init_b_a_unc           Initial accel. bias uncertainty (m/s^2)
%     .init_b_g_unc           Initial gyro. bias uncertainty (rad/s)
%     .init_clock_offset_unc  Initial clock offset uncertainty per axis (m)
%     .init_clock_drift_unc   Initial clock drift uncertainty per axis (m/s)
%     .gyro_noise_PSD         Gyro noise PSD (rad^2/s)
%     .accel_noise_PSD        Accelerometer noise PSD (m^2 s^-3)
%     .accel_bias_PSD         Accelerometer bias random walk PSD (m^2 s^-5)
%     .gyro_bias_PSD          Gyro bias random walk PSD (rad^2 s^-3)
%     .clock_freq_PSD         Receiver clock frequency-drift PSD (m^2/s^3)
%     .clock_phase_PSD        Receiver clock phase-drift PSD (m^2/s)
%     .pseudo_range_SD        Pseudo-range measurement noise SD (m)
%     .range_rate_SD          Pseudo-range rate measurement noise SD (m/s)
%
%  Outputs:
%   out_profile        Navigation solution as a motion profile array
%   out_errors         Navigation solution error array
%   out_IMU_bias_est   Kalman filter IMU bias estimate array
%   out_clock          GNSS Receiver clock estimate array
%   out_KF_SD          Output Kalman filter state uncertainties
% notes  :
*-----------------------------------------------------------------------------*/
void Tightly_coupled_INS_GNSS(INS_measurements *INS_measurements,
                              GNSS_measurements *GNSS_measurements,
                              int no_GNSS_meas, const obsd_t *obs, const nav_t *nav, PVAT_solution *pvat_old, int no_par, float old_time,
                              float time_last_GNSS, double *clock_offset_drift,
                              initialization_errors *initialization_errors, IMU_errors *IMU_errors,
                              GNSS_config *GNSS_config, TC_KF_config *TC_KF_config, PVAT_solution *pvat_new,
                              double *out_errors, double *out_IMU_bias_est, double *out_clock,
                              double *out_KF_SD, double *meas_f_ib_b, double *meas_omega_ib_b)
{
  double true_L_b, true_lambda_b, true_h_b;
  double old_est_llh[3], old_est_r_eb_e[3] = {0.0}, old_est_v_eb_e[3] = {0.0},
                         old_true_C_b_e[9] = {0.0};
  double est_clock[2], est_r_eb_e[3];
  double old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n[3];
  double true_v_eb_n[3], true_C_b_n[9] = {0.0}, true_eul_nb[3], old_est_C_b_n[9] = {0.0};
  double old_true_r_eb_e[3], old_true_v_eb_e[3];
  double est_L_b, est_h_b, est_lambda_b, est_v_eb_e[3], est_v_eb_n[3];
  double est_C_b_e[9] = {0.0}, est_C_b_n[9] = {0.0};
  double delta_r_eb_n[3], delta_v_eb_n[3], delta_eul_nb_n[3] = {0.0};
  double est_IMU_bias[6] = {0.0};
  double quant_residuals[6] = {0};
  double *P_matrix, *P_matrix_new;
  double tor_i, time, tor_s;
  int i, j;
  double est_C_b_e_new[9] = {0.0}, est_C_b_n_T[9], est_v_eb_e_new[3] = {0.0}, est_r_eb_e_new[3] = {0.0};
  double est_IMU_bias_new[6] = {0.0}, est_clock_new[2] = {0.0};
  double q_nb[4], llh[3];
  double C_Transp[9] = {0.0}, checkP = 0.0;

  printf("\n *****************  TC_INS/GNSS BEGINS ************************\n");
  /*
  printf("P: %lf, %lf, %lf\n",pvat_old->latitude*R2D,pvat_old->longitude*R2D,pvat_old->height );
  printf("V: %lf, %lf, %lf\n",pvat_old->ned_velocity[0],pvat_old->ned_velocity[1],pvat_old->ned_velocity[2] );
  printf("A: %lf, %lf, %lf\n",pvat_old->euler_angles[0],pvat_old->euler_angles[1],pvat_old->euler_angles[2] );*/

  /* Initialize true navigation solution */
  time = INS_measurements->sec;
  est_L_b = true_L_b = pvat_old->latitude;
  true_lambda_b = pvat_old->longitude;
  true_h_b = pvat_old->height;
  for (i = 0; i < 3; i++)
    true_v_eb_n[i] = pvat_old->ned_velocity[i];
  for (i = 0; i < 3; i++)
    true_eul_nb[i] = pvat_old->euler_angles[i];
  Euler_to_CTM(true_eul_nb, true_C_b_n);

  printf("true_eul_nb: %lf, %lf, %lf\n", true_eul_nb[0], true_eul_nb[1], true_eul_nb[2]);

  /* Transposing "true_C_b_n" (in _n_b) to actual _b_n frame*/
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      C_Transp[j * 3 + i] = true_C_b_n[i * 3 + j];
    }
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      true_C_b_n[i * 3 + j] = C_Transp[i * 3 + j];
    }
  }

  printf("true_C_b_n from euler\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%.10lf ", true_C_b_n[i * 3 + j]);
    }
    printf("\n");
  }

  if (pvat_old->Nav_or_KF)
  {
    /* Previous sol was an KF Integrated */
    NED_to_ECEF(&true_L_b, &true_lambda_b, &true_h_b, true_v_eb_n, true_C_b_n,
                old_est_r_eb_e, old_est_v_eb_e, old_true_C_b_e);

    /* Attitude matrix from previous KF or Navigation */
    //for (i=0;i<9;i++) old_true_C_b_e[i]=pvat_old->C_b_e[i];
  }
  else
  {
    /* Previous sol. was Navigation only */
    /* Previous attitude from initialization */
    NED_to_ECEF(&true_L_b, &true_lambda_b, &true_h_b, true_v_eb_n, true_C_b_n,
                old_est_r_eb_e, old_est_v_eb_e, old_true_C_b_e);

    if (norm(pvat_old->C_b_e, 3) <= 0.0)
    {
      /* First solution, use what is coming  */
    }
    else
    {
      /* Attitude matrix from previous KF or Navigation */
      //  for (i=0;i<9;i++) old_true_C_b_e[i]=pvat_old->C_b_e[i];
    }
  }

  /* Initialize other matrices */
  for (j = 0; j < 2; j++)
    est_clock[j] = clock_offset_drift[j];

  /* Initialize estimated attitude solution */
  //Initialize_NED_attitude(true_C_b_n, initialization_errors, old_est_C_b_n);

  /* Initialize Kalman filter P matrix and IMU bias states */
  P_matrix = zeros(no_par, no_par); // in this case 17x17
  P_matrix_new = zeros(no_par, no_par);

  checkP = Sumdiag(pvat_old->P, no_par, 6, 9); /* For position only (as ignav)*/

  printf("CHECK FOR P: checkP:%lf > MAXVAR: %lf\n", checkP, MAXVAR);

  if (norm(out_errors, 17) > 0.02 && (checkP / 3) < MAXVAR)
  {
    printf("ERRROS ARE NOT ZERO \n");

    for (i = 0; i < 17; i++)
    {
      for (j = 0; j < 17; j++)
      {
        P_matrix[i * 17 + j] = pvat_old->P[i * 17 + j];
      }
    }
  }
  else
  {
    printf("INITIALIZING P_matrix TO DEFAULT\n");
    Initialize_TC_P_matrix(TC_KF_config, P_matrix); //IT WILL CHANGE FOR THE UNDIFF/UNCOMB MODEL
  }

  /* Generate IMU bias and clock output records
  out_IMU_bias_est(1,1) = old_time;
  out_IMU_bias_est(1,2:7) = est_IMU_bias';
  out_clock(1,1) = old_time;
  out_clock(1,2:3) = est_clock;  */

  /* Generate KF uncertainty record */

  /* Main loop */

  //for (i = 0; i < no_obs; i++) {
  /* Time interval */
  tor_i = time - old_time;

  if (old_time <= 0.000)
  {
    tor_i = 1.0; //1 second of time interval
    time_last_GNSS = time - tor_i;
  }

  /* Correct IMU errors */
  for (i = 0; i < 6; i++)
    est_IMU_bias[i] = out_IMU_bias_est[i];
  for (i = 0; i < 3; i++)
    meas_f_ib_b[i] = INS_measurements->f_ib_b[i];
  for (i = 0; i < 3; i++)
    meas_omega_ib_b[i] = INS_measurements->omega_ib_b[i];
  for (i = 0; i < 3; i++)
    meas_f_ib_b[i] = meas_f_ib_b[i] - est_IMU_bias[i];
  for (i = 0; i < 3; i++)
    meas_omega_ib_b[i] = meas_omega_ib_b[i] - est_IMU_bias[i + 3];

  /* Update estimated navigation solution */
  printf("\n *****************  NAV_EQUATIONS BEGINS ************************\n");
  Nav_equations_ECEF(tor_i,
                     old_est_r_eb_e, old_est_v_eb_e, old_true_C_b_e, meas_f_ib_b,
                     meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
  printf("\n *****************  NAV_EQUATIONS ENDS **************************\n");

  /* INS covariance filter propagation */
  propinsstateTC(P_matrix, est_r_eb_e, est_C_b_e, meas_f_ib_b, TC_KF_config, tor_i);

  printf("Norm of en_uncert.: t: %lf norm: %lf\n", GNSS_measurements->sec, norm(TC_KF_config->init_pos_unc_ned, 2));

  printf("GNSS.INS.Horizontal.velocities: %lf %lf \n", norm(pvagnss.v, 2), norm(est_v_eb_e, 2));
  /**/
  if (fabs(GNSS_measurements->sec - 243340.00) <= 0.0001 || fabs(GNSS_measurements->sec - 243341.00) <= 0.0001 ||
      fabs(GNSS_measurements->sec - 243339.00) <= 0.0001)
  {
    GNSS_measurements[0].gdop[0] = 3.0;
    printf("Epochs.taken.out: %lf\n", GNSS_measurements->sec);
  }

  /* Convert navigation solution to NED
  ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e,\
    &true_L_b,&est_lambda_b,&est_h_b,\
    est_v_eb_n,est_C_b_n);
    */

  if (fabs(norm(pvagnss.v, 2)) < 5)
  {
    printf("Epochs.GNSS_Vel_Jumps: %lf\n", GNSS_measurements->sec);
  }

  //&& fabs(norm(pvagnss.v,2)-norm(est_v_eb_e,2)) < 5

  /* Determine whether to run Kalman filter */
  if (fabs(GNSS_measurements->sec - INS_measurements->sec) < 0.0001 && GNSS_measurements[0].gdop[0] < 2.5 && no_GNSS_meas >= 4 &&
      norm(TC_KF_config->init_pos_unc_ned, 2) < 5.0)
  { //gdops: 2.4,2.8
    //no_GNSS_meas>=4 && norm(TC_KF_config->init_pos_unc_ned, 2)<4.0 &&

    /* KF time interval */
    tor_s = time - time_last_GNSS;

    /* Run Integration Kalman filter */
    printf("\n *******************  TC_KF_EPOCH BEGINS ************************\n");
    /**/ TC_KF_Epoch(GNSS_measurements, no_GNSS_meas, obs, nav, tor_s, est_C_b_e,
                     est_v_eb_e, est_r_eb_e, est_IMU_bias, est_clock, P_matrix,
                     meas_f_ib_b, est_L_b, TC_KF_config, est_C_b_e_new,
                     est_v_eb_e_new, est_r_eb_e_new, est_IMU_bias_new,
                     est_clock_new, P_matrix_new);
    printf("\n *******************  TC_KF_EPOCH ENDS **************************\n");

    printf("\n *****************  NAV_EQUATIONS CLOSED-LOOP BEGINS **************\n");
    /* Correct IMU errors
for (i=0;i<3;i++) meas_f_ib_b[i]=INS_measurements->f_ib_b[i];
for (i=0;i<3;i++) meas_omega_ib_b[i]=INS_measurements->omega_ib_b[i];
for (i=0;i<3;i++) meas_f_ib_b[i] = meas_f_ib_b[i] - est_IMU_bias_new[i];
for (i=0;i<3;i++) meas_omega_ib_b[i] = meas_omega_ib_b[i] - est_IMU_bias_new[i+3];

    Nav_equations_ECEF(tor_i,\
        old_est_r_eb_e, old_est_v_eb_e,old_true_C_b_e, meas_f_ib_b,\
        meas_omega_ib_b, est_r_eb_e_new, est_v_eb_e_new, est_C_b_e_new);
*/
    printf("\n *****************  NAV_EQUATIONS CLOSED-LOOP ENDS ****************\n");

    time_last_GNSS = GNSS_measurements->sec;

    pvat_new->Nav_or_KF = 1;

    /* Generate IMU bias and clock output records */
    fprintf(out_clock_file, "%lf %lf %lf %d\n", time, est_clock_new[0],
            est_clock_new[1], pvat_new->Nav_or_KF);

    for (i = 0; i < 6; i++)
      out_IMU_bias_est[i] = est_IMU_bias_new[i];

    /* Generate KF uncertainty output record */
    fprintf(out_KF_SD_file, "%lf\t", time);
    for (i = 0; i < no_par; i++)
    {
      for (j = 0; j < no_par; j++)
      {
        if (i == j)
        {
          out_errors[j] = P_matrix_new[i * no_par + j];
          fprintf(out_KF_SD_file, "%lf\t", sqrt(P_matrix_new[i * no_par + j]));
        }
      }
    }
    fprintf(out_KF_SD_file, "%d", pvat_new->Nav_or_KF);
    fprintf(out_KF_SD_file, "\n");

    /* Full weight matrix */
    for (i = 0; i < 17; i++)
    {
      for (j = 0; j < 17; j++)
      {
        pvat_new->P[i * 17 + j] = P_matrix_new[i * 17 + j];
      }
    }
  }
  else
  {
    pvat_new->Nav_or_KF = 0;
    if (fabs(GNSS_measurements->sec - INS_measurements->sec) < 0.0001)
    {
      printf("WAS.NOT.INTEGRATED: time: %lf, DOP:%lf<2.5, GNSSmeas: %d>=4, POsunc: %lf<3.0, VEL: %lf< 5m/s \n",
             GNSS_measurements->sec,
             GNSS_measurements[0].gdop[0], no_GNSS_meas, norm(TC_KF_config->init_pos_unc_ned, 2),
             (norm(pvagnss.v, 2) - norm(est_v_eb_e, 2)));
    }

    /* Generate KF uncertainty output record */
    fprintf(out_KF_SD_file, "%lf\t", time);
    for (i = 0; i < no_par; i++)
    {
      for (j = 0; j < no_par; j++)
      {
        if (i == j)
        {
          out_errors[j] = P_matrix[i * no_par + j];
          fprintf(out_KF_SD_file, "%lf\t", sqrt(P_matrix[i * no_par + j]));
        }
      }
    }
    fprintf(out_KF_SD_file, "%d", pvat_new->Nav_or_KF);
    fprintf(out_KF_SD_file, "\n");

    /* Full weight matrix */
    for (i = 0; i < 17; i++)
    {
      for (j = 0; j < 17; j++)
      {
        pvat_new->P[i * 17 + j] = P_matrix[i * 17 + j];
      }
    }
    for (i = 0; i < 6; i++)
      out_IMU_bias_est[i] = est_IMU_bias[i];
    for (i = 0; i < 3; i++)
      est_r_eb_e_new[i] = est_r_eb_e[i];
    for (i = 0; i < 3; i++)
      est_v_eb_e_new[i] = est_v_eb_e[i];
    for (i = 0; i < 9; i++)
      est_C_b_e_new[i] = est_C_b_e[i];
    for (i = 0; i < 2; i++)
      est_clock_new[i] = est_clock[i];

    /* Generate IMU bias and clock output records */
    fprintf(out_clock_file, "%lf %lf %lf %d\n", time, est_clock_new[0],
            est_clock_new[1], pvat_new->Nav_or_KF);
  }

  fprintf(out_IMU_bias_file, "%lf %lf %lf %lf %.10lf %.10lf %.10lf %d\n", time,
          out_IMU_bias_est[0], out_IMU_bias_est[1], out_IMU_bias_est[2],
          out_IMU_bias_est[3], out_IMU_bias_est[4], out_IMU_bias_est[5], pvat_new->Nav_or_KF);

  /* Convert navigation solution to NED  */
  ECEF_to_NED(est_r_eb_e_new, est_v_eb_e_new, est_C_b_e_new,
              &est_L_b, &est_lambda_b, &est_h_b,
              est_v_eb_n, est_C_b_n);

  printf("est_C_b_n\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%.13lf ", est_C_b_n[i * 3 + j]);
    }
    printf("\n");
  }
  printf("\n");

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      C_Transp[j * 3 + i] = est_C_b_n[i * 3 + j];
    }
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      est_C_b_n_T[i * 3 + j] = C_Transp[i * 3 + j];
    }
  }

  double q_b_n[4];
  DCM_to_quaternion(est_C_b_n, q_b_n);
  Quaternion_to_euler(q_b_n, pvat_new->euler_angles);

  printf("Euler from Quaternions: %lf, %lf, %lf\n", pvat_new->euler_angles[0],
         pvat_new->euler_angles[1], pvat_new->euler_angles[2]);

  pvat_new->latitude = est_L_b;
  pvat_new->longitude = est_lambda_b;
  pvat_new->height = est_h_b;
  for (i = 0; i < 3; i++)
    pvat_new->ned_velocity[i] = est_v_eb_n[i];
  for (i = 0; i < 9; i++)
    pvat_new->C_b_n[i] = est_C_b_n[i];
  CTM_to_Euler(pvat_new->euler_angles, est_C_b_n_T);
  pvat_new->time = time;

  printf("Euler from CTM: %lf, %lf, %lf\n", pvat_new->euler_angles[0],
         pvat_new->euler_angles[1], pvat_new->euler_angles[2]);

  /* Clock update */
  for (i = 0; i < 2; i++)
    out_clock[i] = est_clock_new[i];

  /* Passing Cbe from KF direclty to the next Navigation */
  for (i = 0; i < 9; i++)
    pvat_new->C_b_e[i] = est_C_b_e_new[i];
  for (i = 0; i < 3; i++)
    pvat_new->re[i] = est_r_eb_e_new[i];
  for (i = 0; i < 3; i++)
    pvat_new->ve[i] = est_v_eb_e_new[i];

  /* Free memory */
  free(P_matrix);
  free(P_matrix_new);
  printf("\n *****************  TC_INS/GNSS ENDS **************************\n");
}

/////////////////  MAIN FUNCTION ///////////////////////////////////////////////
/* Main function -------------------------------------------------------------*/
/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
//int main (void){
//int InsGnssCore (){ LATER, WHEN LINKED TO OTHER MAIN FUNCTION USE IT THIS WAY
extern void LC_INS_GNSS_core(double *GNSS_r_eb_e, double *gnss_xyz_ini_cov, double *gnss_ned_cov,
                             double *GNSS_v_eb_e, double *gnss_xyz_vel_cov,
                             double ini_pos_time, um7pack_t *imu_data, double imu_time_diff,
                             pva_t *PVA_old, imuraw_t *imuobsp, int Tact_or_Low_IMU, int tc_or_lc)
{
  int no_par = 15, no_meas = 6;
  double out_errors[15] = {0}, out_IMU_bias_est[6] = {0}, out_KF_SD[15];
  double meas_f_ib_b[3], meas_omega_ib_b[3];
  IMU_errors IMU_errors = {0};
  initialization_errors initialization_errors = {0};
  GNSS_config GNSS_config = {{0}};
  LC_KF_config LC_KF_config = {{0}};
  TC_KF_config TC_KF_config = {{0}};
  PVAT_solution pvat_old = {{0}};
  PVAT_solution pvat_new = {{0}};
  INS_measurements INS_meas = {0};
  double llh[3], enu_ini_vel[3], rr[3];
  double old_time = 0.0, last_GNSS_time = 0.0, curr_GNSS_time;
  int i, j, k;

  printf("\n *****************  LC INSGNSS CORE BEGINS ***********************\n");

  /* Load GNSS and INS errors and noises from configuration file -------------*/
  if (Tact_or_Low_IMU)
  {
    /* Tactical-grade imu */
    Initialize_INS_GNSS_LCKF_tactical_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,
                                                &LC_KF_config);
  }
  else
  {
    /*Low-grade imu */
    Initialize_INS_GNSS_LCKF_consumer_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,
                                                &LC_KF_config);
  }

  memset(&INS_meas, 0, sizeof(INS_measurements));
  /**/
  printf("IMU_data time: %f\n", imu_data->sec);
  printf("PVA_old.A: %lf %lf %lf\n", PVA_old->A[0], PVA_old->A[1], PVA_old->A[2]);

  /* Current INS time */
  INS_meas.sec = imu_data->sec;
  INS_meas.tt = imu_time_diff;

  /* Current GNSS time */
  curr_GNSS_time = ini_pos_time;

  /* Last GNSS or State time and Last imu time */
  if (PVA_old->sec <= 0.00)
  {
    old_time = pvat_old.time = PVA_old->sec = INS_meas.sec - 0.5; /* Initialize according to IMU rate*/
    INS_meas.tt = 0.5;
    old_time = (float)INS_meas.sec - INS_meas.tt;
    last_GNSS_time = ini_pos_time - 1.0; /* Initialize according to GNSS rate*/
  }
  else
  {
    old_time = pvat_old.time = PVA_old->sec;
    last_GNSS_time = PVA_old->t_s; /* last state time */
  }
  /*
    printf("Time: %lf  - timediff: %lf and old_time: %lf\n", imu_data->sec,imu_time_diff, old_time);
    printf("OLDTimes: %lf  %lf\n", old_time, pvat_old.time);*/

  /* Initial position and velocity from GNSS or previous PVA solution */
  if (imu_time_diff < 0.1)
  { /* If first solution, use GNSS */
    ecef2pos(GNSS_r_eb_e, llh);
    pvat_old.latitude = llh[0];
    pvat_old.longitude = llh[1];
    pvat_old.height = llh[2];
    ecef2enu(llh, GNSS_v_eb_e, enu_ini_vel); /* Here is ENU! */
                                             /* ENU to NED velocity from gnss solution */
    pvat_old.ned_velocity[0] = enu_ini_vel[1];
    pvat_old.ned_velocity[1] = enu_ini_vel[0];
    pvat_old.ned_velocity[2] = -enu_ini_vel[2];

    /* IMU bias */
    for (i = 0; i < 3; i++)
      out_IMU_bias_est[i] = IMU_errors.b_a[i];
    for (i = 0; i < 3; i++)
      out_IMU_bias_est[i + 3] = IMU_errors.b_g[i];
    for (i = 0; i < 3; i++)
      pvat_old.euler_angles[i] = PVA_old->A[i]; /* in _nb frame */
    /* Type of solution */
    PVA_old->Nav_or_KF = pvat_new.Nav_or_KF = 0;
    /* Attitude matrix */
    for (i = 0; i < 9; i++)
      pvat_old.C_b_e[i] = PVA_old->Cbe[i];
  }
  else
  { /* Otherwise, use previous PVA solution */
    pvat_old.latitude = PVA_old->r[0];
    pvat_old.longitude = PVA_old->r[1];
    pvat_old.height = PVA_old->r[2];
    for (i = 0; i < 3; i++)
      pvat_old.ned_velocity[i] = PVA_old->v[i]; /*already in NED */
    for (i = 0; i < 3; i++)
      pvat_old.euler_angles[i] = PVA_old->A[i]; /* in _nb frame */

    /* Attitude matrix */
    for (i = 0; i < 9; i++)
      pvat_old.C_b_e[i] = PVA_old->Cbe[i];

    /* IMU bias */
    for (i = 0; i < 6; i++)
      out_IMU_bias_est[i] = PVA_old->out_IMU_bias_est[i];
    for (i = 0; i < 17; i++)
      out_errors[i] = PVA_old->out_errors[i];
  }

  /* Type of preiou solution  */
  pvat_old.Nav_or_KF = PVA_old->Nav_or_KF;

  /* Current IMU measurement structure initialization */
  for (j = 0; j < 3; j++)
    INS_meas.omega_ib_b[j] = imu_data->g[j];
  for (j = 0; j < 3; j++)
    INS_meas.f_ib_b[j] = imu_data->a[j];

  /* IMU bias */
  for (i = 0; i < 6; i++)
    out_IMU_bias_est[i] = PVA_old->out_IMU_bias_est[i];
  for (i = 0; i < no_par; i++)
    out_errors[i] = PVA_old->out_errors[i];

  /* Full P matrix */
  for (i = 0; i < no_par; i++)
  {
    for (j = 0; j < no_par; j++)
    {
      pvat_old.P[i * no_par + j] = PVA_old->P[i * no_par + j];
    }
  }

  /* Position covariance from GNSS */
  for (j = 0; j < 3; j++)
    TC_KF_config.init_pos_unc[j] = sqrt(gnss_xyz_ini_cov[j]); /* For the filter */
  for (j = 0; j < 3; j++)
    LC_KF_config.init_pos_unc_ned[j] = sqrt(gnss_ned_cov[j]); /* For testing the integration condition */
  for (j = 0; j < 3; j++)
    LC_KF_config.init_vel_unc[j] = sqrt(gnss_xyz_vel_cov[j]);

  printf("LC_KF_config.init_pos_unc_ned[j]:%lf %lf %lf\n", LC_KF_config.init_pos_unc_ned[0],
         LC_KF_config.init_pos_unc_ned[1], LC_KF_config.init_pos_unc_ned[2]);

  /* Loosely coupled ECEF Inertial navigation and GNSS integrated navigation -*/
  Loosely_coupled_INS_GNSS(&INS_meas, no_meas, &pvat_old, no_par, old_time,
                           last_GNSS_time, curr_GNSS_time, GNSS_r_eb_e, GNSS_v_eb_e, &initialization_errors, &IMU_errors,
                           &GNSS_config, &LC_KF_config, &pvat_new, out_errors, out_IMU_bias_est,
                           out_KF_SD, meas_f_ib_b, meas_omega_ib_b);

  PVA_old->r[0] = pvat_new.latitude;
  PVA_old->r[1] = pvat_new.longitude;
  PVA_old->r[2] = pvat_new.height;
  //pos2ecef(llh, PVA_old->r);
  for (i = 0; i < 3; i++)
    PVA_old->v[i] = pvat_new.ned_velocity[i];
  for (i = 0; i < 3; i++)
    PVA_old->A[i] = pvat_new.euler_angles[i]; /*in ?????*/
  PVA_old->sec = pvat_new.time;

  for (i = 0; i < 9; i++)
    PVA_old->Cbe[i] = pvat_new.C_b_e[i];
  for (i = 0; i < 3; i++)
    PVA_old->re[i] = pvat_new.re[i];
  for (i = 0; i < 3; i++)
    PVA_old->ve[i] = pvat_new.ve[i];

  /* Full weight matrix  */
  for (i = 0; i < no_par; i++)
  {
    for (j = 0; j < no_par; j++)
    {
      PVA_old->P[i * no_par + j] = pvat_new.P[i * no_par + j];
    }
  }

  for (i = 0; i < 6; i++)
    PVA_old->out_IMU_bias_est[i] = out_IMU_bias_est[i];
  for (i = 0; i < no_par; i++)
    PVA_old->out_errors[i] = out_errors[i];

  /* Type of solution  */
  PVA_old->Nav_or_KF = pvat_new.Nav_or_KF;

  //PVA_old->t_s=curr_GNSS_time;

  /*
printf("Pllhnew: %lf, %lf, %lf\n",pvat_new.latitude*R2D,pvat_new.longitude*R2D,pvat_new.height );
printf("Vnew: %lf, %lf, %lf\n",pvat_new.ned_velocity[0],pvat_new.ned_velocity[1],pvat_new.ned_velocity[2] );
printf("Anew: %lf, %lf, %lf\n",pvat_new.euler_angles[0],pvat_new.euler_angles[1],pvat_new.euler_angles[2] );
printf("dtr, dtrs: %lf, %lf\n", PVA_old->clock_offset_drift[0],PVA_old->clock_offset_drift[1]);
*/
  /* Plots -------------------------------------------------------------------*/

  /* Write output profile and errors file ------------------------------------*/

  /* Free memory -------------------------------------------------------------*/

  printf("\n *****************  INSGNSS CORE ENDS *************************\n");
  return 0;
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
//int main (void){
//int InsGnssCore (){ LATER, WHEN LINKED TO OTHER MAIN FUNCTION USE IT THIS WAY
extern void TC_INS_GNSS_core(rtk_t *rtk, const obsd_t *obs, int n,
                             const nav_t *nav, um7pack_t *imu_data, double imu_time_diff, double *gnss_ned_cov,
                             double *gnss_vel_ned_cov, pva_t *PVA_old, int Tact_or_Low_IMU, int tc_or_lc)
{

  int no_par = 17, no_GNSS_meas;
  double out_errors[17] = {0.0}, out_IMU_bias_est[6] = {0.0}, out_clock[2], clock_offset[2], out_KF_SD[17];
  double meas_f_ib_b[3], meas_omega_ib_b[3];
  IMU_errors IMU_errors = {0};
  initialization_errors initialization_errors = {0};
  GNSS_config GNSS_config = {{0}};
  PVAT_solution pvat_old = {{0}};
  PVAT_solution pvat_new = {{0}};
  GNSS_measurements GNSS_measurements[n];
  INS_measurements INS_meas = {0};
  TC_KF_config TC_KF_config = {{0}};
  double *rs, *dts, *var, *azel, llh[3], enu_ini_vel[3];
  double r, rr[3], e[3], dion, vion, dtrp, vtrp, lam_L1, P, vs[3];
  double old_time = 0.0, last_GNSS_time = 0.0;
  int svh[MAXOBS], i, j, k, sys, flag[n], m;
  prcopt_t opt = rtk->opt;
  double x[20] = {0.0}, dtdx[3]; //This two are only used when Tropo gradientes are estimated in rtklib

  printf("\n *****************  INSGNSS CORE BEGINS ***********************\n");

  /* Checking input values
  printf("POS: %lf, %lf %lf\n",rtk->sol.rr[0],rtk->sol.rr[1], rtk->sol.rr[2] );
  printf("OBS: %lf, %lf %lf\n",obs[0].P[0],obs[8].L[0], obs[10].D[0] );
  printf("Number of obs: %d\n",n);
  printf("IMUDATA: ax:%lf, ay:%lf, az:%lf and gx:%lf, gy:%lf, gz:%lf\n",\
  imu_data->a[0], imu_data->a[1],\
  imu_data->a[2], imu_data->g[0],imu_data->g[1], imu_data->g[2] );*/

  /*  printf("P: %lf, %lf, %lf\n",PVA_old->re[0],PVA_old->re[1],PVA_old->re[2] );
  printf("V: %lf, %lf, %lf\n",PVA_old->v[0],PVA_old->v[1],PVA_old->v[2] );
  printf("A: %lf, %lf, %lf\n",PVA_old->A[0],PVA_old->A[1],PVA_old->A[2] );
  printf("dtr, dtrs: %lf, %lf\n", PVA_old->clock_offset_drift[0],PVA_old->clock_offset_drift[1]);*/

  /* Prepare GNSS and INS raw data into the proper structures --------------- */

  /* Load GNSS and INS errors and noises from configuration file -------------*/
  if (Tact_or_Low_IMU)
  {
    /* Tactical-grade imu */
    Initialize_INS_GNSS_TCKF_tactical_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,
                                                &TC_KF_config);
  }
  else
  {
    /*Low-grade imu */
    Initialize_INS_GNSS_TCKF_consumer_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,
                                                &TC_KF_config);
  }

  rs = mat(6, n);
  dts = mat(2, n);
  var = mat(1, n);
  azel = zeros(2, n);

  memset(&INS_meas, 0, sizeof(INS_measurements));
  memset(&GNSS_measurements, 0, sizeof(GNSS_measurements));
  /* current GNSS time */
  GNSS_measurements->sec = time2gpst(rtk->sol.time, NULL); /* time of week in (GPST) */
  GNSS_measurements->tt = rtk->tt;
  /* Current INS time */
  INS_meas.sec = imu_data->sec;
  INS_meas.tt = imu_time_diff;

  printf("IMU_data time: %f\n", imu_data->sec);
  printf("PVA_old.A: %lf %lf %lf\n", PVA_old->A[0], PVA_old->A[1], PVA_old->A[2]);

  /* Last GNSS or State time and Last imu time */
  if (PVA_old->sec <= 0.00)
  {
    old_time = pvat_old.time = PVA_old->sec = INS_meas.sec - 0.5; /* Initialize according to IMU rate*/
    //INS_meas.tt=0.5;
    old_time = (float)INS_meas.sec - INS_meas.tt;
    last_GNSS_time = GNSS_measurements->sec - 1.0; /* Initialize according to GNSS rate*/
  }
  else
  {
    old_time = pvat_old.time = PVA_old->sec;
    last_GNSS_time = PVA_old->t_s;
  }

  /* Initial position and velocity from GNSS or previous PVA solution */
  if (rtk->tt < 0.1)
  { /* If first solution, use GNSS */
    ecef2pos(rtk->sol.rr, llh);
    pvat_old.latitude = llh[0];
    pvat_old.longitude = llh[1];
    pvat_old.height = llh[2];
    ecef2enu(llh, rtk->sol.rr + 3, enu_ini_vel); /* Here is ENU! */
                                                 /* ENU to NED velocity from gnss solution */
    pvat_old.ned_velocity[0] = enu_ini_vel[1];
    pvat_old.ned_velocity[1] = enu_ini_vel[0];
    pvat_old.ned_velocity[2] = -enu_ini_vel[2];

    /* initialize clock offset and drift from gnss */
    clock_offset[0] = rtk->sol.dtr[0] * CLIGHT; //rtk->sol.dtr[0];//]x[3]; //or = rtk->sol.dtr[0];
    clock_offset[1] = rtk->sol.dtrr;

    /* IMU bias */
    for (i = 0; i < 3; i++)
      out_IMU_bias_est[i] = IMU_errors.b_a[i];
    for (i = 0; i < 3; i++)
      out_IMU_bias_est[i + 3] = IMU_errors.b_g[i];
    for (i = 0; i < 3; i++)
      pvat_old.euler_angles[i] = PVA_old->A[i]; /* in _nb frame */
    /* Type of solution */
    PVA_old->Nav_or_KF = pvat_new.Nav_or_KF = 0;
    /* Attitude matrix */
    for (i = 0; i < 9; i++)
      pvat_old.C_b_e[i] = PVA_old->Cbe[i];
  }
  else
  { /* Otherwise, use previous PVA solution */
    //ecef2pos(PVA_old->re,PVA_old->r);
    pvat_old.latitude = PVA_old->r[0];
    pvat_old.longitude = PVA_old->r[1];
    pvat_old.height = PVA_old->r[2];
    for (i = 0; i < 3; i++)
      pvat_old.ned_velocity[i] = PVA_old->v[i]; /*already in NED */
    for (i = 0; i < 3; i++)
      pvat_old.euler_angles[i] = PVA_old->A[i]; /* in _nb frame */

    /* initialize clock offset and drift from previous PVAT */
    clock_offset[0] = PVA_old->clock_offset_drift[0];
    clock_offset[1] = PVA_old->clock_offset_drift[1];

    /* Attitude matrix */
    for (i = 0; i < 9; i++)
      pvat_old.C_b_e[i] = PVA_old->Cbe[i];

    /* IMU bias */
    for (i = 0; i < 6; i++)
      out_IMU_bias_est[i] = PVA_old->out_IMU_bias_est[i];
    for (i = 0; i < 17; i++)
      out_errors[i] = PVA_old->out_errors[i];
  }

  /* Type of preiou solution  */
  pvat_old.Nav_or_KF = PVA_old->Nav_or_KF;

  /* Constraining height with GNSS ones since there are no jumps */
  ecef2pos(rtk->sol.rr, llh);
  //pvat_old.height = llh[2];

  /* initialize clock offset and drift from gnss  */
  clock_offset[0] = rtk->sol.dtr[0] * CLIGHT; //rtk->sol.dtr[0];//]x[3]; //or = rtk->sol.dtr[0];
  clock_offset[1] = rtk->sol.dtrr;

  /* From gyrocompassing and levelling from current IMU meas.
  for(j=0;j<3;j++) pvat_old.euler_angles[j]=imu_data->aea[j];*/

  /* Current IMU measurement structure initialization */
  for (j = 0; j < 3; j++)
    INS_meas.omega_ib_b[j] = imu_data->g[j];
  for (j = 0; j < 3; j++)
    INS_meas.f_ib_b[j] = imu_data->a[j];

  /* IMU bias */
  for (i = 0; i < 6; i++)
    out_IMU_bias_est[i] = PVA_old->out_IMU_bias_est[i];
  for (i = 0; i < 17; i++)
    out_errors[i] = PVA_old->out_errors[i];

  /* Full P matrix */
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      pvat_old.P[i * 17 + j] = PVA_old->P[i * 17 + j];
    }
  }

  /* Receiver position */
  for (i = 0; i < 3; i++)
    rr[i] = rtk->sol.rr[i];
  //dtr=rtk->sol.dtr[0]; //(sec).
  ecef2pos(rr, llh);

  /* satellite positions and clocks */
  /* rs[(0:2)+i*6] {x,y,z} is the satellite position in ECEF
     rs[(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
     dts[(0:1)+i*2] are the sat clock {bias,drift} (s|s/s)*/
  satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);
  m = 0;

  printf("Number of sat: %d and valid: %d\n", n, rtk->sol.ns);
  for (i = 0; i < n && i < MAXOBS; i++)
  {

    azel[i * 2] = azel[1 + i * 2] = 0.0;
    lam_L1 = nav->lam[obs[i].sat - 1][0];

    if (!(sys = satsys(obs[i].sat, NULL)))
    {
      flag[i] = 0;
      m++;
      printf("SYST. ERROR\n");
      continue;
    }

    /* reject duplicated observation data */
    if (i < n - 1 && i < MAXOBS - 1 && obs[i].sat == obs[i + 1].sat)
    {
      printf("duplicated observation data %s sat=%2d\n",
             time_str(obs[i].time, 3), obs[i].sat);
      i++;
      flag[i] = 0;
      m++;
      continue;
    }

    /* geometric distance/azimuth/elevation angle */
    //printf("geodist: %lf, %lf, %lf\n", geodist(rs+i*6,rr,e), satazel(llh,e,azel+i*2), opt.elmin );
    if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
        satazel(llh, e, azel + i * 2) < opt.elmin)
    {
      printf("ERROR GEOM.\n");
      flag[i] = 0;
      m++;
      continue;
    }

    /* psudorange with code bias correction   -> IT IS FOR P3
    if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;*/

    /* excluded satellite? */
    if (satexclude(obs[i].sat, svh[i], &opt))
    {
      flag[i] = 0;
      m++;
      printf("EXC. SAT. ERROR\n");
      continue;
    }

    /* ionospheric corrections */
    if (!ionocorr(obs[i].time, nav, obs[i].sat, llh, azel + i * 2,
                  IONOOPT_BRDC, &dion, &vion))
    {
      printf("ERROR IONO!\n");
      flag[i] = 0;
      m++;
      continue;
    } //opt.ionoopt

    /* GPS-L1 -> L1/B1 */
    if ((lam_L1 = nav->lam[obs[i].sat - 1][0]) > 0.0)
    {
      dion *= (lam_L1 / lam_carr[0]) * (lam_L1 / lam_carr[0]);
    }

    /* tropospheric corrections */
    if (!tropcorr(obs[i].time, nav, llh, azel + i * 2,
                  TROPOPT_SAAS, &dtrp, &vtrp))
    { //opt.tropopt

      printf("ERROR TROPO!\n");
      flag[i] = 0;
      m++;
      continue;
    }
    dtrp = prectrop(obs[i].time, llh, azel + i * 2, &opt, x, dtdx, &vtrp);
    printf("Tropo values: %lf\n", dtrp);

    if (obs[i].P[0] <= 0.0 || fabs(obs[i].D[0]) <= 0.0 || obs[i].L[0] <= 0.0)
    {
      printf("OBS BUG ERROR\n");
      flag[i] = 0;
      m++;
      continue;
    }

    /* pseudorange residual
      printf("Pres: %lf\n", obs[i].P[0]-(r+CLIGHT*rtk->sol.dtr[0]-CLIGHT*dts[i*2]+dion+dtrp));
      printf("Dres: %lf\n", -obs[i].D[j]*lam_L1+CLIGHT*dts[1+i*2]);
      //printf("RANGE: %lf, ION: %lf, Tropo: %lf, SAT clock: %lf, REC clk: %lf \
      REC.POS: %lf, %lf, %lf, Azel: %lf, PVA_oldpos: %lf,  %lf, %lf\n ", r, dion, dtrp, CLIGHT*dts[i*2], \
      CLIGHT*rtk->sol.dtr[0], llh[0], llh[1],llh[2],azel[i*2], PVA_old->r[0],PVA_old->r[1],PVA_old->r[2]);
      */
    //  printf("Pres: %lf\n", obs[i].P[j]-dion-dtrp+CLIGHT*dts[i*2]);
    //  printf("Dres: %lf\n", -obs[i].D[j]*lam_L1+CLIGHT*dts[1+i*2]);
    printf("RHO: %lf corected and uncorrected: %lf \n", r + CLIGHT * rtk->sol.dtr[0], r);
    printf("Pres: %lf\n", obs[i].P[0] - (r + CLIGHT * rtk->sol.dtr[0] - CLIGHT * dts[i * 2] + dion + dtrp));

    flag[i] = 1;
    GNSS_measurements[i].sat = obs[i].sat;
    GNSS_measurements[i].time = obs[i].time;
    for (j = 0; j < 2; j++)
      GNSS_measurements[i].P[j] = obs[i].P[j] - dion - dtrp + CLIGHT * dts[i * 2]; //-dion-dtrp//-clock_offset[0];
    for (j = 0; j < 2; j++)
      GNSS_measurements[i].L[j] = obs[i].L[j];
    for (j = 0; j < 2; j++)
      GNSS_measurements[i].D[j] = -obs[i].D[j] * lam_L1 + CLIGHT * dts[1 + i * 2]; //-clock_offset[1]; /* hz to m/s (radial velocity)*/
    for (j = 0; j < 3; j++)
      GNSS_measurements[i].Sat_r_eb_e[j] = rs[i * 6 + j];
    for (j = 0; j < 3; j++)
      GNSS_measurements[i].Sat_v_eb_e[j] = rs[i * 6 + (j + 3)];
  }

  /* Re-arranging GNSS_measurements structure to eliminate invalid  satellites */
  k = 0;
  for (i = 0; i < n; i++)
  {
    if (flag[i])
    {
      /* Valid satellite values */
      GNSS_measurements[i - k].sat = GNSS_measurements[i].sat;
      GNSS_measurements[i - k].time = GNSS_measurements[i].time;
      for (j = 0; j < 2; j++)
        GNSS_measurements[i - k].P[j] = GNSS_measurements[i].P[j];
      for (j = 0; j < 2; j++)
        GNSS_measurements[i - k].L[j] = GNSS_measurements[i].L[j];
      for (j = 0; j < 2; j++)
        GNSS_measurements[i - k].D[j] = GNSS_measurements[i].D[j]; /* hz to m/s (radial velocity)*/
      for (j = 0; j < 3; j++)
        GNSS_measurements[i - k].Sat_r_eb_e[j] = GNSS_measurements[i].Sat_r_eb_e[j];
      for (j = 0; j < 3; j++)
        GNSS_measurements[i - k].Sat_v_eb_e[j] = GNSS_measurements[i].Sat_v_eb_e[j];
    }
    else
    {
      /*Skip position - invalid satelite*/
      k++;
    }
  }

  /* Put zeros on the invalid positions?? */
  if (m > 0)
  {
    for (i = 1; i <= m; i++)
    {
      GNSS_measurements[n - i].sat = 0;
      for (j = 0; j < 2; j++)
        GNSS_measurements[n - i].P[j] = 0.0;
      for (j = 0; j < 2; j++)
        GNSS_measurements[n - i].L[j] = 0.0;
      for (j = 0; j < 2; j++)
        GNSS_measurements[n - i].D[j] = 0.0;
      for (j = 0; j < 3; j++)
        GNSS_measurements[n - i].Sat_r_eb_e[j] = 0.0;
      for (j = 0; j < 3; j++)
        GNSS_measurements[n - i].Sat_v_eb_e[j] = 0.0;
    }
  }

  /* Number of valid satellites */
  no_GNSS_meas = n - m;
  printf("Number of sat: after loop %d \n", no_GNSS_meas);

  /* Passing DOPS to the first satellite */
  for (i = 0; i < 4; i++)
    GNSS_measurements[0].gdop[i] = rtk->sol.gdop[i];
  /**/
  for (i = 0; i < no_GNSS_meas; i++)
  {
    printf("SAT: %2d, P: %lf, D: %lf, L: %lf\n", GNSS_measurements[i].sat,
           GNSS_measurements[i].P[0], GNSS_measurements[i].D[0], GNSS_measurements[i].L[0]);
  }

  /* Load GNSS and INS errors and noises from configuration file -------------*/

  /* Position and velocity covariance from GNSS */
  if (rtk->sol.qr[0] > 100.0 || rtk->sol.qr[1] > 100.0 || rtk->sol.qr[2] > 100.0)
  {
    /* FLAG TO CONTINUE INS NAVIGATION ONLY!! */
  }
  else
  {
    for (j = 0; j < 3; j++)
      TC_KF_config.init_pos_unc[j] = sqrt(rtk->sol.qr[j]);
  }

  if (rtk->sol.qrv[0] > 100.0 || rtk->sol.qrv[1] > 100.0 || rtk->sol.qrv[2] > 100.0 ||
      rtk->sol.qrv[0] <= 0.0 || rtk->sol.qrv[1] <= 0.0 || rtk->sol.qrv[2] <= 0.0)
  {
    /* FLAG TO CONTINUE INS NAVIGATION ONLY!! */
  }
  else
  {
    for (j = 0; j < 3; j++)
      TC_KF_config.init_vel_unc[j] = sqrt(rtk->sol.qrv[j]);
  }

  for (j = 0; j < 3; j++)
    TC_KF_config.init_pos_unc_ned[j] = sqrt(gnss_ned_cov[j]);

  /**/
  TC_KF_config.init_clock_offset_unc = sqrt(rtk->sol.qdtr); // It is in s should be in m ?
  TC_KF_config.init_clock_drift_unc = sqrt(rtk->sol.qdtrr); // In (s/s) it should be in m/s?
                                                            /**/

  /* Tightly coupled ECEF Inertial navigation and GNSS integrated navigation -*/
  Tightly_coupled_INS_GNSS(&INS_meas, &GNSS_measurements, no_GNSS_meas,
                           obs, nav,
                           &pvat_old, no_par, old_time, last_GNSS_time, clock_offset, &initialization_errors, &IMU_errors,
                           &GNSS_config, &TC_KF_config, &pvat_new, out_errors, out_IMU_bias_est,
                           out_clock, out_KF_SD, meas_f_ib_b, meas_omega_ib_b);

  PVA_old->r[0] = pvat_new.latitude;
  PVA_old->r[1] = pvat_new.longitude;
  PVA_old->r[2] = pvat_new.height;
  //pos2ecef(llh, PVA_old->r);
  for (i = 0; i < 3; i++)
    PVA_old->v[i] = pvat_new.ned_velocity[i];
  for (i = 0; i < 3; i++)
    PVA_old->A[i] = pvat_new.euler_angles[i]; /*in ?????*/
  PVA_old->sec = pvat_new.time;

  for (i = 0; i < 9; i++)
    PVA_old->Cbe[i] = pvat_new.C_b_e[i];
  for (i = 0; i < 3; i++)
    PVA_old->re[i] = pvat_new.re[i];
  for (i = 0; i < 3; i++)
    PVA_old->ve[i] = pvat_new.ve[i];

  /* Full weight matrix  */
  for (i = 0; i < 17; i++)
  {
    for (j = 0; j < 17; j++)
    {
      PVA_old->P[i * 17 + j] = pvat_new.P[i * 17 + j];
    }
  }

  for (i = 0; i < 2; i++)
    PVA_old->clock_offset_drift[i] = out_clock[i];
  for (i = 0; i < 6; i++)
    PVA_old->out_IMU_bias_est[i] = out_IMU_bias_est[i];
  for (i = 0; i < 17; i++)
    PVA_old->out_errors[i] = out_errors[i];

  PVA_old->t_s = last_GNSS_time;

  /* Type of solution */
  PVA_old->Nav_or_KF = pvat_new.Nav_or_KF;

  /* Plots -------------------------------------------------------------------*/

  /* Write output profile and errors file ------------------------------------*/

  /* Free memory -------------------------------------------------------------*/
  free(rs);
  free(dts);
  free(var);
  free(azel);
  printf("\n *****************  INSGNSS CORE ENDS *************************\n");
}

extern void imu_tactical_navigation(FILE *imu_file)
{
  char str[100];
  float t_prev = 0.0, t_curr, tor_i;
  double old_r_eb_e[3], old_v_eb_n[3], old_v_eb_e[3], old_C_b_e[9];
  double f_ib_b[3], omega_ib_b[3], f_ib_n[3], C_b_n[9], eul_nb_n[3];
  double r_eb_e[3], v_eb_e[3], C_b_e[9], eul_eb_e[3];
  double est_r_eb_e[3], est_v_eb_e[3], est_C_b_e[3];
  double C_Transp[9];
  double wiee[3], llh[3] = {0.0}, gan[3], q[4];
  um7pack_t imu = {0};
  double G = 9.80665;
  int i, j;
  imuraw_t imu_obs_prev = {0};
  um7pack_t imu_curr_meas = {0};
  pva_t PVA_prev_sol = {{0}};

  /* Update with previous solution */
  imu_obs_prev = imu_obs_global;
  PVA_prev_sol = pva_global;

  /* Initialize position, velocity and attitude */
  old_r_eb_e[0] = 1761300.949;
  old_r_eb_e[1] = -4078202.553;
  old_r_eb_e[2] = 4561403.792;

  /* Local or navigation-frame velocities and acceleration */
  old_v_eb_n[0] = 0.0; //-0.02;
  f_ib_n[0] = -0.144;  //E-W
  old_v_eb_n[1] = 0.0; //-0.01;
  f_ib_n[1] = -0.2455; //N-S
  old_v_eb_n[2] = 0.0; // 0.04;
  f_ib_n[2] = 0.605;   //U-D

  /* Attitude Initialization (Groves, 2013)  */

  /* Earth rotation vector in e-frame	*/
  wiee[0] = 0;
  wiee[1] = 0;
  wiee[2] = OMGE;

  ecef2pos(old_r_eb_e, llh);
  //llh[0]=45.950319743*D2R; llh[1]=-66.641372715*D2R;
  /* local apparent gravity vector */
  appgrav(llh, gan, wiee);

  /* Coarse alignment or use of previous solution? */
  /* Levelling and gyrocompassing, update imu->aea[] vector */
  for (i = 0; i < 3; i++)
    imu.a[i] = f_ib_n[i];
  for (i = 0; i < 3; i++)
    imu.v[i] = old_v_eb_n[i];

  coarseAlign(imu.a, imu.g, imu.v, gan, imu.aea);

  for (i = 0; i < 3; i++)
    eul_nb_n[i] = imu.aea[i];

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      C_Transp[j * 3 + i] = C_b_n[i * 3 + j];
    }
  }

  printf("POSb: %lf, %lf, %lf\n", old_r_eb_e[0], old_r_eb_e[1], old_r_eb_e[2]);
  printf("llh: %lf, %lf, %lf\n", llh[0], llh[1], llh[2]);

  NED_to_ECEF(&llh[0], &llh[1], &llh[2], old_v_eb_n, C_Transp,
              old_r_eb_e, old_v_eb_e, old_C_b_e);

  printf("POSa: %lf, %lf, %lf\n", old_r_eb_e[0], old_r_eb_e[1], old_r_eb_e[2]);

  printf("Here 1\n");
  /*
 printf("Pllhnew: %lf, %lf, %lf\n",llh[0]*R2D,llh[1]*R2D,llh );
 printf("eul_nb_n: %lf, %lf, %lf\n",eul_nb_n[0],eul_nb_n[1],eul_nb_n[2] );
 printf("C_b_n:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
  C_b_n[0],C_b_n[1],C_b_n[2],\
  C_b_n[3],C_b_n[4],C_b_n[5],\
  C_b_n[6],C_b_n[7],C_b_n[8]);

  printf("old_r_eb_e: %lf, %lf, %lf\n",old_r_eb_e[0],old_r_eb_e[1],old_r_eb_e[2] );
  printf("old_v_eb_e: %lf, %lf, %lf\n",old_v_eb_e[0],old_v_eb_e[1],old_v_eb_e[2] );
  printf("old_C_b_e:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
   old_C_b_e[0],old_C_b_e[1],old_C_b_e[2],\
   old_C_b_e[3],old_C_b_e[4],old_C_b_e[5],\
   old_C_b_e[6],old_C_b_e[7],old_C_b_e[8]);

*/
  fgets(str, 100, imu_file);
  sscanf(str, "%f %lf %lf %lf %lf %lf %lf", &t_prev, &f_ib_b[2],
         &f_ib_b[1], &f_ib_b[0], &omega_ib_b[2], &omega_ib_b[1], &omega_ib_b[0]);

  printf("time begninning: %f\n", t_prev);

  while (fgets(str, 100, imu_file) != NULL)
  {

    sscanf(str, "%f %lf %lf %lf %lf %lf %lf", &t_curr, &f_ib_b[2],
           &f_ib_b[1], &f_ib_b[0], &omega_ib_b[2], &omega_ib_b[1], &omega_ib_b[0]);

    tor_i = t_curr - t_prev;

    printf("Read time: %f and Diff: %f\n", t_curr, tor_i);

    /* Turn-on bias */
    /* Turn-on initial biases             //exp4*/
    f_ib_b[0] = (f_ib_b[0]) * G; //-0.0339421
    f_ib_b[1] = (f_ib_b[1]) * G; //+0.105076
    f_ib_b[2] = (f_ib_b[2]) * G; //+1.00335
    omega_ib_b[0] = (omega_ib_b[0]) * D2R;
    omega_ib_b[1] = (omega_ib_b[1]) * D2R;
    omega_ib_b[2] = (omega_ib_b[2]) * D2R;
    /*
  printf("old_r_eb_e: %lf, %lf, %lf\n",old_r_eb_e[0],old_r_eb_e[1],old_r_eb_e[2] );
  printf("old_v_eb_e: %lf, %lf, %lf\n",old_v_eb_e[0],old_v_eb_e[1],old_v_eb_e[2] );
  printf("old_C_b_e:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
   old_C_b_e[0],old_C_b_e[1],old_C_b_e[2],\
   old_C_b_e[3],old_C_b_e[4],old_C_b_e[5],\
   old_C_b_e[6],old_C_b_e[7],old_C_b_e[8]);
  printf("f_ib_b %lf, %lf, %lf\n",f_ib_b[0],f_ib_b[1],f_ib_b[2] );
  printf("omega_ib_b: %lf, %lf, %lf\n",omega_ib_b[0],omega_ib_b[1],omega_ib_b[2] );
  printf("tor_i: %f\n", tor_i);
*/

    printf("\n *****************  NAV_EQUATIONS BEGINS ************************\n");
    Nav_equations_ECEF(tor_i,
                       old_r_eb_e, old_v_eb_e, old_C_b_e, f_ib_b,
                       omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
    printf("\n *****************  NAV_EQUATIONS ENDS **************************\n");

    for (i = 0; i < 3; i++)
      old_r_eb_e[i] = est_r_eb_e[i];
    for (i = 0; i < 3; i++)
      old_v_eb_e[i] = est_v_eb_e[i];
    for (i = 0; i < 9; i++)
      old_C_b_e[i] = est_C_b_e[i];
    t_prev = t_curr;

    printf("POS: %lf, %lf, %lf\n", old_r_eb_e[0], old_r_eb_e[1], old_r_eb_e[2]);

    /* Convert navigation solution to NED  */
    ECEF_to_NED(est_r_eb_e, est_v_eb_e, old_C_b_e, &llh[0], &llh[1], &llh[2],
                old_v_eb_n, C_b_n);

    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        C_Transp[j * 3 + i] = C_b_n[i * 3 + j];
      }
    }
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        C_b_n[i * 3 + j] = C_Transp[i * 3 + j];
      }
    }

    /* Cbn to euler */
    CTM_to_Euler(eul_nb_n, C_b_n);

    ecef2pos(est_r_eb_e, llh);

    //Quaternion_to_euler(q, eul_nb_n);
    printf("eul_nb_n: %lf, %lf, %lf\n", eul_nb_n[0], eul_nb_n[1], eul_nb_n[2]);

    fprintf(out_PVA, "%f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            t_curr, llh[0] * R2D, llh[1] * R2D, llh[2],
            old_v_eb_e[0], old_v_eb_e[1], old_v_eb_e[2],
            eul_nb_n[0] * R2D, eul_nb_n[1] * R2D, eul_nb_n[2] * R2D);
  }
}

/* Name of function ------------------------------------------------------------
* Brief description
*Nav_equations_ECEF - Runs precision ECEF-frame inertial navigation
 %equations
 %
 % Software for use with "Principles of GNSS, Inertial, and Multisensor
 % Integrated Navigation Systems," Second Edition.
 %
 % This function created 1/4/2012 by Paul Groves
 %
 % Inputs:
 %   tor_i         time interval between epochs (s)
 %   old_r_eb_e    previous Cartesian position of body frame w.r.t. ECEF
 %                 frame, resolved along ECEF-frame axes (m)
 %   old_C_b_e     previous body-to-ECEF-frame coordinate transformation matrix
 %   old_v_eb_e    previous velocity of body frame w.r.t. ECEF frame, resolved
 %                 along ECEF-frame axes (m/s)
 %   f_ib_b        specific force of body frame w.r.t. ECEF frame, resolved
 %                 along body-frame axes, averaged over time interval (m/s^2)
 %   omega_ib_b    angular rate of body frame w.r.t. ECEF frame, resolved
 %                 about body-frame axes, averaged over time interval (rad/s)
 % Outputs:
 %   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
 %                 along ECEF-frame axes (m)
 %   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
 %                 ECEF-frame axes (m/s)
 %   C_b_e         body-to-ECEF-frame coordinate transformation matrix

 % Copyright 2012, Paul Groves
 % License: BSD; see license.txt for details
*-----------------------------------------------------------------------------*/
void Nav_equations_ECEF1(ins_states_t *ins)
{
  double tor_i;
  double omega_ie = 7.292115E-5; // Earth rotation rate (rad/s)
  double alpha_ie, mag_alpha;
  double C_Earth[9], alpha_ib_b[3], Alpha_ib_b[9];
  double Alpha_squared[9], second_term[9], first_term[9], C_new_old[9];
  double I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double C_aux[9], first_term_2[9], second_term_2[9], C_b_e_Cbb[9];
  double ave_C_b_e[9], Cbb[9], alpha_ie_vec[3], Alpha_ie[9], last_term[9];
  double f_ib_e[3], omega_ie_vec[3], Omega_ie[9], g[3], Omega_v_eb_e[3];
  double new_q_b_e[4], new_C_b_e[9];
  int i, j;

  tor_i = ins->dt;
  for (i = 0; i < 3; i++)
  {
    ins->data.fb[i] = ins->data.fb0[i] - ins->data.ba[i];
    ins->data.wibb[i] = ins->data.wibb0[i] - ins->data.bg[i];
  }
  /**/
  printf(" ********************* NAVIGATION EQUATIONS ******************\nINPUT: \n");
  printf("tor_i= %f;\n", tor_i);
  printf("ins->pre=[%lf; %lf; %lf];\n", ins->pre[0], ins->pre[1], ins->pre[2]);
  printf("llh: %lf, %lf, %lf\n", ins->prn[0]*R2D, ins->prn[1]*R2D, ins->prn[0]);
  printf("ins->pve=[%lf; %lf; %lf];\n", ins->pve[0], ins->pve[1], ins->pve[2]);
  printf("ins->data.fb=[%lf; %lf; %lf];\n", ins->data.fb[0], ins->data.fb[1], ins->data.fb[2]);
  printf("ins->data.wibb=[%lf, %lf, %lf];\n", ins->data.wibb[0], ins->data.wibb[1], ins->data.wibb[2]);
  printf("Ba=[%lf; %lf; %lf];\n", ins->data.ba[0], ins->data.ba[1], ins->data.ba[2]);
  printf("Bg=[%lf, %lf, %lf];\n", ins->data.bg[0], ins->data.bg[1], ins->data.bg[2]);
  printf("ins->pCbe=[");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf, ", ins->pCbe[i * 3 + j]); 
    }
    printf("];\n");
  }
   printf("ins->P\n");
  for (i = 0; i < ins->nx; i++)
  {
    for (j = 0; j < ins->nx; j++)
    {
      (i==j?printf("%lf ", ins->P[i * ins->nx + j]):0);
    }
  }
  printf("\n");

  /* Begins     */

  /* ATTITUDE UPDATE  */
  /* From (2.145) determine the Earth rotation over the update interval
  C_Earth = C_e_i' * old_C_e_i  */
  alpha_ie = omega_ie * tor_i;
  C_Earth[0] = cos(alpha_ie);
  C_Earth[3] = sin(alpha_ie);
  C_Earth[6] = 0.0;
  C_Earth[1] = -sin(alpha_ie);
  C_Earth[4] = cos(alpha_ie);
  C_Earth[7] = 0.0;
  C_Earth[2] = 0.0;
  C_Earth[5] = 0.0;
  C_Earth[8] = 1.0;

  omega_ie_vec[0] = 0;
  omega_ie_vec[1] = 0;
  omega_ie_vec[2] = omega_ie;

  /* Calculate attitude increment, magnitude, and skew-symmetric matrix  */
  for (i = 0; i < 3; i++)
    alpha_ib_b[i] = ins->data.wibb[i] * tor_i;
  mag_alpha = norm(alpha_ib_b, 3);
  skewsym3(alpha_ib_b, Alpha_ib_b);

  /* Attitude Update using Quaternion algebra 
  attitude_update(alpha_ib_b, omega_ie_vec, tor_i, ins->pCbe, new_q_b_e);
  Quaternion_to_DCM(new_q_b_e, new_C_b_e);
  */

  /* Obtain coordinate transformation matrix from the new attitude w.r.t. an
  inertial frame to the old using Rodrigues' formula, (5.73)  */
  matmul("NN", 3, 3, 3, 1.0, Alpha_ib_b, Alpha_ib_b, 0.0, Alpha_squared);

  for (i = 0; i < 9; i++)
    second_term[i] = (1 - cos(mag_alpha)) /
                     (mag_alpha * mag_alpha) * Alpha_squared[i];
  for (i = 0; i < 9; i++)
    first_term[i] = I[i] + sin(mag_alpha) /
                               mag_alpha * Alpha_ib_b[i];

  if (mag_alpha > 1.E-8)
  {
    for (i = 0; i < 9; i++)
      C_new_old[i] = first_term[i] + second_term[i];
  }
  else
  {
    for (i = 0; i < 9; i++)
      C_new_old[i] = I[i] + Alpha_ib_b[i];
  } // end if mag_alpha

  /* Update attitude using (5.75)  */
  matmul("NN", 3, 3, 3, 1.0, C_Earth, ins->pCbe, 0.0, C_aux);
  matmul("NN", 3, 3, 3, 1.0, C_aux, C_new_old, 0.0, ins->Cbe);

  /* SPECIFIC FORCE FRAME TRANSFORMATION
  % Calculate the average body-to-ECEF-frame coordinate transformation
  % matrix over the update interval using (5.84) and (5.85)  */
  for (i = 0; i < 9; i++)
    first_term_2[i] = I[i] + second_term[i];
  for (i = 0; i < 9; i++)
    second_term_2[i] = ((1 - (sin(mag_alpha) / mag_alpha)) /
                        (mag_alpha * mag_alpha)) *
                       Alpha_squared[i];
  for (i = 0; i < 9; i++)
    Cbb[i] = first_term_2[i] + second_term_2[i];
  alpha_ie_vec[0] = 0;
  alpha_ie_vec[1] = 0;
  alpha_ie_vec[2] = alpha_ie;
  skewsym3(alpha_ie_vec, Alpha_ie);
  matmul("NN", 3, 3, 3, 0.5, Alpha_ie, ins->pCbe, 0.0, last_term);

  if (mag_alpha > 1.E-8)
  {
    matmul("NN", 3, 3, 3, 1.0, ins->pCbe, Cbb, 0.0, C_b_e_Cbb);
    for (i = 0; i < 9; i++)
      ave_C_b_e[i] = C_b_e_Cbb[i] - last_term[i];
  }
  else
  {
    for (i = 0; i < 9; i++)
      ave_C_b_e[i] = ins->pCbe[i] - last_term[i];
  } //if mag_alpha

  /* Transform specific force to ECEF-frame resolving axes using (5.85) */
  matmul("NN", 3, 1, 3, 1.0, ave_C_b_e, ins->data.fb, 0.0, f_ib_e);

  /* UPDATE VELOCITY
  % From (5.36), */
  skewsym3(omega_ie_vec, Omega_ie);
  Gravity_ECEF(ins->pre, g);
  matmul("NN", 3, 1, 3, 1.0, Omega_ie, ins->pve, 0.0, Omega_v_eb_e);
  //matmul("NN", 1, 3, 3, 1.0, ins->pve, Omega_ie, 0.0, Omega_v_eb_e);
  for (i = 0; i < 3; i++){
    ins->data.fbe[i]=(f_ib_e[i] + g[i] - 2 * Omega_v_eb_e[i]);
    ins->ve[i] = ins->pve[i] + ins->data.fbe[i] * tor_i;
  }

  /* UPDATE CARTESIAN POSITION
   % From (5.38), */
  for (i = 0; i < 3; i++)
    ins->re[i] = ins->pre[i] + (ins->ve[i] + ins->pve[i]) * 0.5 * tor_i;
/**/
  printf("OUTPUT NAV: \n");
  double llh[3]={0.0};
  ecef2pos(ins->re,llh);
  printf("P: %lf, %lf, %lf\n", ins->re[0], ins->re[1], ins->re[2]);
  printf("llh: %lf, %lf, %lf\n", llh[0]*R2D, llh[1]*R2D, llh[2]);
  printf("V: %lf, %lf, %lf\n", ins->ve[0], ins->ve[1], ins->ve[2]);
  printf("Ba=[%lf; %lf; %lf];\n", ins->data.ba[0], ins->data.ba[1], ins->data.ba[2]);
  printf("Bg=[%lf, %lf, %lf];\n", ins->data.bg[0], ins->data.bg[1], ins->data.bg[2]);
  printf("ins->Cbe\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf ", ins->Cbe[i * 3 + j]);
    }
    printf("\n");
  }
  printf(" ********************* NAVIGATION EQUATIONS END ******************\n");
}

/* functions for initial error covariance matrix ----------------------------*/
extern void initP(int is, int ni, int nx, double unc, double unc0, double *P0)
{
  int i, j;
  for (i = is; i < is + ni; i++)
    for (j = 0; j < nx; j++)
    {
      if (j == i)
        P0[j + i * nx] = SQR(unc == 0.0 ? unc0 : unc);
      else
        P0[j + i * nx] = P0[i + j * nx] = 0.0;
    }
}

/* system noise covariance matrix--------------------------------------------*/
static void sysQ(int is, int n, int nx, double v, double dt, double *Q)
{
  int i;
  for (i = is; i < is + n; i++)
    Q[i + i * nx] = v * fabs(dt);
}

/* determine approximate system noise covariance matrix ---------------------*/
static void getQ(const insgnss_opt_t *opt, double dt, double *Q, int nx)
{
  trace(3, "getQ:\n");
  printf("getQ:\n"); 

  setzero(Q, nx, nx);

  sysQ(IA, NA, nx, opt->psd.gyro, dt, Q);
  sysQ(IV, NV, nx, opt->psd.accl, dt, Q);
  sysQ(iba, nba, nx, opt->psd.ba, dt, Q);
  sysQ(ibg, nbg, nx, opt->psd.bg, dt, Q);
  sysQ(irc, nrc, nx, opt->psd.clk, dt, Q);
  //sysQ(irr, nrr, nx, opt->psd.clkr, dt, Q);
  
}
/* process noise gain matrix-------------------------------------------------*/
static void getGn(const insgnss_opt_t *opt, const ins_states_t *ins, const double dt,
                  double *Gn, int nx)
{
  int nprn = NNPX;
  double *I = eye(3);

  setzero(Gn, nx, nprn);
  asi_blk_mat(Gn, nx, nprn, ins->Cbe, 3, 3, 0, 9);
  asi_blk_mat(Gn, nx, nprn, ins->Cbe, 3, 3, 3, 0);

  asi_blk_mat(Gn, nx, nprn, I, 3, 3, 9, 6);
  asi_blk_mat(Gn, nx, nprn, I, 3, 3, 12, 9);
  free(I);
}
/* process noise covariance matrix-------------------------------------------*/
static void getprn(const ins_states_t *ins, const insgnss_opt_t *opt, double dt, double *Q, int nx)
{
  int nprn = NNPX, i;
  double *Qn = zeros(nprn, nprn), *Gn = mat(nprn, nx);

  for (i = INAC; i < INGY + NNAC; i++)
    Qn[i + i * nprn] = opt->psd.gyro * fabs(dt);
  for (i = INGY; i < INGY + NNGY; i++)
    Qn[i + i * nprn] = opt->psd.accl * fabs(dt);
  for (i = INBA; i < INBA + NNBA; i++)
    Qn[i + i * nprn] = opt->psd.ba * fabs(dt);
  for (i = INBG; i < INBG + NNBG; i++)
    Qn[i + i * nprn] = opt->psd.bg * fabs(dt);

  getGn(opt, ins, fabs(dt), Gn, nx);
  matmul33("NNT", Gn, Qn, Gn, nx, nprn, nprn, nx, Q);
  free(Qn);
  free(Gn);
}
/* initial error covariance matrix-------------------------------------------*/
extern void getP0(const insgnss_opt_t *opt, double *P0, int nx)
{
  trace(3, "getP0:\n");
  printf("getP0 in propins \n");

  setzero(P0, nx, nx);

  initP(IA, NA, nx, opt->unc.att, UNC_ATT, P0);
  initP(IV, NV, nx, opt->unc.vel, UNC_VEL, P0);
  initP(IP, NP, nx, opt->unc.pos, UNC_POS, P0);
  initP(iba, nba, nx, opt->unc.ba, UNC_BA, P0);
  initP(ibg, nbg, nx, opt->unc.bg, UNC_BG, P0);
  initP(irc, nrc, nx, opt->unc.rc, UNC_CLK, P0);
  //initP(irr, nrr, nx, opt->unc.rr, UNC_CLKR, P0);
}
/* geocentric radius---------------------------------------------------------*/
extern double georadi(const double *pos)
{
  return RE_WGS84 / sqrt(1.0 - SQR(e_exc * sin(pos[0]))) *
         sqrt(SQR(cos(pos[0])) + SQR(1.0 - SQR(e_exc)) * SQR(sin(pos[0]))); //returns a vector
}
/* propagate matrix for stochastic parameters--------------------------------*/
static void stochasticPhi(int opt, int ix, int nix, int nx, double dt, double *phi)
{
  int i;
  if (nix <= 0)
    return;
  for (i = ix; i < ix + nix; i++)
  {
    if (opt == INS_RANDOM_CONS)
      phi[i + i * nx] = 1.0;
    if (opt == INS_RANDOM_WALK)
      phi[i + i * nx] = 1.0;
    if (opt == INS_GAUSS_MARKOV)
      phi[i + i * nx] = exp(-fabs(dt) / CORRETIME);
  }
}
/* system matrix for accl-bias,gyro-bias,sg,sa and so on --------------------*/
static void stochasticF(int opt, int ix, int nix, int nx, double *F)
{
  int i;
  if (nix <= 0)
    return;
  for (i = ix; i < ix + nix; i++)
  {
    if (opt == INS_RANDOM_CONS)
      F[i + i * nx] = 1E-10;
    if (opt == INS_RANDOM_WALK)
      F[i + i * nx] = 1E-10;
    if (opt == INS_GAUSS_MARKOV)
      F[i + i * nx] = -1.0 / CORRETIME;
  }
}
/* system matrix of ins states propagate in ecef-frame-----------------------*/
static void getF(const insgnss_opt_t *opt, const double *Cbe, const double *pos,
                 const double *omgb, const double *fib, double *F, int nx)
{
  int i, j;
  double F21[9], F23[9], I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}, omega[3], rn[3], ge[3], re;
  double W[18] = {0}, WC[18] = {0};

  printf("getF: %d\n", nx);

  setzero(F, nx, nx);

  matmul3("NN", Cbe, fib, omega);
  skewsym3(omega, F21);

  ecef2pos(pos, rn);
  Gravity_ECEF(pos, ge); //=pregrav(pos,ge);
  re = georadi(rn);
  matmul("NT", 3, 3, 1, -2.0 / (re * norm(pos, 3)), ge, pos, 0.0, F23);

  /* ins attitude system matrix */
  W[0] = omgb[1];
  W[3] = omgb[2];
  W[7] = omgb[0];
  W[10] = omgb[2];
  W[14] = omgb[0];
  W[17] = omgb[1];
  matmul("NN", 3, 6, 3, 1.0, Cbe, W, 0.0, WC);
  for (i = IA; i < IA + NA; i++)
  {
    for (j = IA; j < IA + NA; j++)
      F[i + j * nx] = -Omge[i - IA + (j - IA) * 3];
    for (j = ibg; j < ibg + nbg; j++)
      F[i + j * nx] = Cbe[i - IA + (j - ibg) * 3];
    //for (j=isg;j<isg+nsg;j++) F[i+j*nx]= Cbe [i-IA+(j-isg)*3]*omgb[j-isg];
    //for (j=irg;j<irg+nrg;j++) F[i+j*nx]= WC  [i-IA+(j-irg)*3];
  }
  /* ins velocity system matrix */
  W[0] = fib[1];
  W[3] = fib[2];
  W[7] = fib[0];
  W[10] = fib[2];
  W[14] = fib[0];
  W[17] = fib[1];
  matmul("NN", 3, 6, 3, 1.0, Cbe, W, 0.0, WC);

  for (i = IV; i < IV + NV; i++)
  {
    for (j = IA; j < IA + NA; j++)
      F[i + j * nx] = -F21[i - IV + (j - IA) * 3];
    for (j = IV; j < IV + NV; j++)
      F[i + j * nx] = -2.0 * Omge[i - IV + (j - IV) * 3];
    for (j = IP; j < IP + NP; j++)
      F[i + j * nx] = F23[i - IV + (j - IP) * 3];
    for (j = iba; j < iba + nba; j++)
      F[i + j * nx] = Cbe[i - IV + (j - iba) * 3];
    //for (j=isa;j<isa+nsa;j++) F[i+j*nx]= Cbe[i-IV+(j-isa)*3]*fib[j-isa];
    //for (j=ira;j<ira+nra;j++) F[i+j*nx]= WC [i-IV+(j-ira)*3];
  }

  /* ins position system matrix */
  for (i = IP; i < IP + NP; i++)
  {
    for (j = IV; j < IV + NV; j++)
      F[i + j * nx] = I[i - IP + (j - IV) * 3];
  }
  /* stochastic parameters system matrix */
  stochasticF(opt->baproopt, iba, nba, nx, F);
  stochasticF(opt->bgproopt, ibg, nbg, nx, F);
}
/* set matrix to eye-matrix--------------------------------------------------*/
extern void seteye(double *A, int n)
{
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i + j * n] = (i == j ? 1.0 : 0.0);
}
/* factorial operation-------------------------------------------------------*/
static double factorial(int n)
{
  if (n == 0)
    return 1;
  return n * factorial(n - 1);
}
/* exponential of a matrix----------------------------------------------------
 * args  :  double *A  I  a input matrix (nxn)
 *          int n      I  rows and cols of matrix A
 *          double *E  O  exponential of a matrix
 * return: none
 * --------------------------------------------------------------------------*/
extern void expmat(const double *A, int n, double *E)
{
  double s, *B, *C;
  int i, j, k;

  printf("expmat: %d\n", n);

  C = mat(n, n);
  B = mat(n, n);
  seteye(E, n);
  seteye(B, n);

  for (i = 0; i < ORDERS; i++)
  {

    s = 1.0 / factorial(i + 1);
    matmul("NN", n, n, n, s, B, A, 0.0, C);

    for (j = 0; j < n; j++)
    {
      for (k = 0; k < n; k++)
        E[j + k * n] += C[j + k * n];
    }
    matcpy(B, C, n, n);
  }
  free(B);
  free(C);
}
/* precise system propagate matrix-------------------------------------------*/
static void precPhi(const insgnss_opt_t *opt, double dt, const double *Cbe,
                    const double *pos, const double *omgb, const double *fib,
                    double *Phi, int nx)
{
  int i;
  double *F = zeros(nx, nx), *FF = zeros(nx, nx);
  double *FFF = zeros(nx, nx), *I = eye(nx);

  getF(opt, Cbe, pos, omgb, fib, F, nx);

  /* third-order approx */
  for (i = 0; i < nx * nx; i++)
    F[i] *= dt;
#if 0
    matmul("NN",nx,nx,nx,1.0,F,F,0.0,FF);
    matmul33("NNN",F,F,F,nx,nx,nx,nx,FFF);

    for (i=0;i<nx*nx;i++) {
        Phi[i]=I[i]+F[i]+0.5*FF[i]+1.0/6.0*FFF[i];
    }
#else
  expmat(F, nx, Phi);
#endif
  free(F);
  free(FF);
  free(FFF);
  free(I);
}
/* determine transition matrix(first-order approx: Phi=I+F*dt)---------------*/
static void getPhi1(const insgnss_opt_t *opt, double dt, const double *Cbe,
                    const double *pos, const double *omgb, const double *fib,
                    double *phi, int nx)
{
  int i, j;
  double omega[3] = {0}, T[9], ge[3], re, rn[3], W[18] = {0}, WC[18] = {0}, Cbv[9];

  trace(3, "getPhi1:\n");

  if (fabs(dt) > MAXDT)
  {
    trace(3, "large time difference between ins and gnss measurements\n");
  }
  /* initial phi to unit matrix */
  seteye(phi, nx);
  seteye(Cbv, 3);

  /* attitude transmit matrix */
  W[0] = omgb[1];
  W[3] = omgb[2];
  W[7] = omgb[0];
  W[10] = omgb[2];
  W[14] = omgb[0];
  W[17] = omgb[1];
  matmul("NN", 3, 6, 3, 1.0, Cbe, W, 0.0, WC);
  for (i = IA; i < IA + NA; i++)
  {
    for (j = IA; j < IA + NA; j++)
      phi[i + j * nx] -= Omge[i - IA + (j - IA) * 3] * dt;
    for (j = ibg; j < ibg + nbg; j++)
      phi[i + j * nx] = Cbe[i - IA + (j - ibg) * 3] * dt;
    //for (j=isg;j<isg+nsg;j++) phi[i+j*nx] =Cbe [i-IA+(j-isg)*3]*omgb[j-isg]*dt;
    //for (j=irg;j<irg+nrg;j++) phi[i+j*nx] =WC  [i-IA+(j-irg)*3]*dt;
  }
  /* velocity transmit matrix */
  matmul3("NN", Cbe, fib, omega);
  skewsym3(omega, T);
  ecef2pos(pos, rn);
  Gravity_ECEF(pos, ge);
  re = georadi(rn);

  W[0] = fib[1];
  W[3] = fib[2];
  W[7] = fib[0];
  W[10] = fib[2];
  W[14] = fib[0];
  W[17] = fib[1];
  matmul("NN", 3, 6, 3, 1.0, Cbe, W, 0.0, WC);

  for (i = IV; i < IV + NV; i++)
  {
    for (j = IV; j < IV + NV; j++)
      phi[i + j * nx] -= 2.0 * Omge[i - IV + (j - IV) * 3] * dt;
    for (j = IA; j < IA + NA; j++)
      phi[i + j * nx] = -T[i - IV + (j - IA) * 3] * dt;
    for (j = IP; j < IP + NP; j++)
      phi[i + j * nx] = -2.0 * dt / (re * norm(pos, 3)) * ge[i - IV] * pos[j - IP];
    for (j = iba; j < iba + nba; j++)
      phi[i + j * nx] = Cbe[i - IV + (j - iba) * 3] * dt;
    //for (j=isa;j<isa+nsa;j++) phi[i+j*nx] =Cbe[i-IV+(j-isa)*3]*fib[j-isa]*dt;
    //for (j=ira;j<ira+nra;j++) phi[i+j*nx] =WC [i-IV+(j-ira)*3]*dt;
  }
  /* position transmit matrix */
  for (i = IP; i < IP + NP; i++)
  {
    for (j = IV; j < IV + NV; j++)
      phi[i + j * nx] = (i - IP) == (j - IV) ? dt : 0.0;
  }
  /* propagate matrix for stochastic parameters */
  stochasticPhi(opt->baproopt, iba, nba, nx, dt, phi);
  stochasticPhi(opt->bgproopt, ibg, nbg, nx, dt, phi);
}
/* propagate state estimation error covariance-------------------------------*/
static void propP(const insgnss_opt_t *opt, const double *Q, const double *phi,
                  const double *P0, double *P, int nx)
{
  int i, j;
  double *PQ = mat(nx, nx), *Phi2 = mat(nx, nx);

  printf("PropP nx: %d\n", nx);

  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++) PQ[i + j * nx] = P0[i + j * nx] + 0.5 * Q[i + j * nx];
  }
  matmul33("NNT", phi, PQ, phi, nx, nx, nx, nx, Phi2);
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++)
      P[i + j * nx] = Phi2[i + j * nx] + 0.5 * Q[i + j * nx];
  }
  /* initialize every epoch for clock (white noise) */
  initP(irc, nrc, nx, opt->unc.rc, UNC_CLK, P);
  
  free(PQ);
  free(Phi2);
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr,
                      const double *dants, double phw, double *L, double *P,
                      double *Lc, double *Pc)
{
    const double *lam=nav->lam[obs->sat-1];
    double C1,C2;
    int i,sys;
    
    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        if (lam[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) continue;
        
        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*lam[i]-dants[i]-dantr[i]-phw*lam[i];
        P[i]=obs->P[i]       -dants[i]-dantr[i];
        
        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        if (obs->code[i]==CODE_L1C) {
            P[i]+=nav->cbias[obs->sat-1][1];
        }
        else if (obs->code[i]==CODE_L2C||obs->code[i]==CODE_L2X||
                 obs->code[i]==CODE_L2L||obs->code[i]==CODE_L2S) {
            P[i]+=nav->cbias[obs->sat-1][2];
#if 0
            L[i]-=0.25*lam[i]; /* 1/4 cycle-shift */
#endif
        }
    }
    /* iono-free LC */
    *Lc=*Pc=0.0;
    sys=satsys(obs->sat,NULL);
    i=(sys&(SYS_GAL|SYS_SBS))?2:1; /* L1/L2 or L1/L5 */
    if (lam[0]==0.0||lam[i]==0.0) return;
    
    C1= SQR(lam[i])/(SQR(lam[i])-SQR(lam[0]));
    C2=-SQR(lam[0])/(SQR(lam[i])-SQR(lam[0]));
    
#if 0
    /* P1-P2 dcb correction (P1->Pc,P2->Pc) */
    if (sys&(SYS_GPS|SYS_GLO|SYS_QZS)) {
        if (P[0]!=0.0) P[0]-=C2*nav->cbias[obs->sat-1][0];
        if (P[1]!=0.0) P[1]+=C1*nav->cbias[obs->sat-1][0];
    }
#endif
    if (L[0]!=0.0&&L[i]!=0.0) *Lc=C1*L[0]+C2*L[i];
    if (P[0]!=0.0&&P[i]!=0.0) *Pc=C1*P[0]+C2*P[i];
}

/* propagate state estimates noting that all states are zero due to closed-loop
 * correction----------------------------------------------------------------*/
static void propx(const insgnss_opt_t *opt, const double *x0, double *x)
{
  int i;
  for (i = 0; i < xnCl(); i++)
    x[i] = 1E-20;
}
/* updates phi,P,Q of ekf------------------------------------------------------- STTOPED HERE*********/
static void updstat(const insgnss_opt_t *opt, ins_states_t *ins, const double dt,
                    const double *x0, const double *P0, double *phi, double *P,
                    double *x, double *Q)
{
  int nx = ins->nx;

  struct timeval stop, start; 
  //printf("updstat\n");

  /* determine approximate system noise covariance matrix */
  //gettimeofday(&start, NULL); //do stuff
  opt->scalePN ? getprn(ins, opt, dt, Q, nx) : getQ(opt, dt, Q, nx);
  //gettimeofday(&stop, NULL); 
  //printf("took Q %lf\n", (stop.tv_usec - start.tv_usec)/1000.0);


  /* determine transition matrix
  * using last epoch ins states (first-order approx) */
  //gettimeofday(&start, NULL); //do stuff
  opt->exphi ? precPhi(opt, dt, ins->Cbe, ins->re, ins->data.wibb, ins->data.fb, phi, nx) : getPhi1(opt, dt, ins->Cbe, ins->re, ins->data.wibb, ins->data.fb, phi, nx);
  //gettimeofday(&stop, NULL); 
  //printf("took Phi %lf\n", (stop.tv_usec - start.tv_usec)/1000.0);

#if UPD_IN_EULER
  getphi_euler(opt, dt, ins->Cbe, ins->re, ins->omgb, ins->fb, phi);
#endif
  //gettimeofday(&start, NULL); //do stuff
  /* propagate state estimation error covariance */
  if (fabs(dt) >= MAXUPDTIMEINT){
    getP0(opt, P, nx);
  }
  else{
    propP(opt, Q, phi, P0, P, nx);
  }
 //gettimeofday(&stop, NULL); 
  //printf("took P %lf\n", (stop.tv_usec - start.tv_usec)/1000.0);
  /* propagate state estimates noting that
     * all states are zero due to close-loop correction */
  if (x)
    propx(opt, x0, x);

  /* predict info. */
  if (ins->P0)
    matcpy(ins->P0, P, nx, nx);
  if (ins->F)
    matcpy(ins->F, phi, nx, nx);
}

/* propagate ins states and its covariance matrix----------------------------
 * args  : insstate_t *ins  IO  ins states
 *         insopt_t *opt    I   ins options
 *         double dt        I   time difference between current and precious
 *         double *x        O   updates ins states
 *         double *P        O   upadtes ins states covariance matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void propinss(ins_states_t *ins, const insgnss_opt_t *opt, double dt,
                     double *x, double *P)
{
  int nx = ins->nx;
  double *phi, *Q;

  Q = mat(nx, nx);
  phi = mat(nx, nx);

  updstat(opt, ins, dt, ins->x, ins->P, phi, P, x, Q);

 int i, j;
//  printf("Phi:\n");
//   for (i = 0; i < ins->nx; i++)
//   {
//     for (j = 0; j < ins->nx; j++)
//     {
//       printf("%lf ", phi[i * ins->nx + j]);
//     }
//   }
//   printf("\n");

   printf("\nQ:\n");
  for (i = 0; i < 20; i++)
  {
    for (j = 0; j < 20; j++)
    {
      if(i==j) printf("%.15lf ", Q[i * ins->nx + j]); 
    }
  }
  printf("\n");
  
    
  free(Q);  
  free(phi);   
}

/* initialize ins/gnss parameter uncertainty with defaul values ------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
//extern void ig_paruncinit(insgnss_opt_t *insopt)
//{
  // int i;

  // trace(3, "ig_paruncernit :\n");

  // insopt->unc.att = UNC_ATT;
  // insopt->unc.vel = UNC_VEL;
  // insopt->unc.pos = UNC_POS;
  // insopt->unc.ba = UNC_BA;
  // insopt->unc.bg = UNC_BG;
  // insopt->unc.dt = UNC_DT;
  // insopt->unc.sg = UNC_SG;
  // insopt->unc.sa = UNC_SA;
  // insopt->unc.rg = UNC_RG;
  // insopt->unc.ra = UNC_RA;
  // insopt->unc.lever = UNC_LEVER;
  // insopt->unc.rc = UNC_CLK;
  // insopt->unc.rr = UNC_CLKR;
//}

/* Initialize Kalman filter uncertainties according to type of imu and processing */
void kf_par_unc_init(insgnss_opt_t *opt)
{
  if (opt->Tact_or_Low)
  {
    /* Tactical */
    opt->unc.att = D2R * 1;          /* Initial attitude uncertainty per axis (deg, converted to rad) */
    opt->unc.vel = 0.1;              /* Initial velocity uncertainty per axis (m/s) */
    opt->unc.pos = 10;               /* Initial position uncertainty per axis (m) */
    opt->unc.ba = 1000 * Mg2M;       /* Initial accelerometer bias uncertainty per instrument 
      (micro-g, converted /* to m/s^2) */
    opt->unc.bg = 10 * D2R / 3600.0; /* Initial gyro bias uncertainty per instrument (deg/hour, 
      converted to rad/sec)*/
    if (opt->mode)
    {
      /* Tightly */
      opt->unc.rc = 10.0;
      opt->unc.rr = 0.1;
    }
  }
  else
  {
    /* Consumer grade */
    opt->unc.att = D2R * 2;           /* Initial attitude uncertainty per axis (deg, converted to rad) */
    opt->unc.vel = 0.1;               /* Initial velocity uncertainty per axis (m/s) */
    opt->unc.pos = 10;                /* Initial position uncertainty per axis (m) */
    opt->unc.ba = 10000 * Mg2M;       /* Initial accelerometer bias uncertainty per instrument 
      (micro-g, converted /* to m/s^2) */
    opt->unc.bg = 200 * D2R / 3600.0; /* Initial gyro bias uncertainty per instrument (deg/hour, 
      converted to rad/sec)*/
    if (opt->mode)
    {
      /* Tightly */
      opt->unc.rc = 10.0;
      opt->unc.rr = 0.1;
    }
  }
}

/* Initialize Kalman filter noise information according to type of imu and processing */
void kf_noise_init(insgnss_opt_t *opt)
{
  if (opt->Tact_or_Low)
  {
    /* Tactical */
    opt->psd.accl = pow((200 * Mg2M), 2);      /* Accelerometer noise PSD (micro-g^2 per Hz, 
      converted to m^2 s^-3) */
    opt->psd.gyro = pow((0.02 * D2R / 60), 2); /* Gyro noise PSD (deg^2 per hour, converted 
      to rad^2/s) */
    opt->psd.ba = 1.0E-7;                      /* Accelerometer bias random walk PSD (m^2 s^-5) */
    opt->psd.bg = 2.0E-12;                     /* Gyro bias random walk PSD (rad^2 s^-3) */
    if (opt->mode)
    {
      /* Tightly */
      opt->psd.clk = 1;  /* Receiver clock frequency-drift PSD (m^2/s^3) */
      opt->psd.clkr = 1; /* Receiver clock phase-drift PSD (m^2/s) */
    }
  }
  else
  {
    /* Consumer grade */
    opt->psd.accl = 0.2 * 0.2;   /* Accelerometer noise PSD (micro-g^2 per Hz, 
      converted to m^2 s^-3) */
    opt->psd.gyro = 0.01 * 0.01; /* Gyro noise PSD (deg^2 per hour, converted 
      to rad^2/s) */
    /* % NOTE: A large noise PSD is modeled to account for the scale-factor and
       % cross-coupling errors that are not directly included in the Kalman filter
       % model*/
    opt->psd.ba = 1.0E-5;  /* Accelerometer bias random walk PSD (m^2 s^-5) */
    opt->psd.bg = 4.0E-11; /* Gyro bias random walk PSD (rad^2 s^-3) */
    if (opt->mode)
    {
      /* Tightly */
      opt->psd.clk = 1;  /* Receiver clock frequency-drift PSD (m^2/s^3) */
      opt->psd.clkr = 1; /* Receiver clock phase-drift PSD (m^2/s) */
    }
  }
}
/* check ins states covariance matrix----------------------------------------*/
static int chkpcov(int nx, const insgnss_opt_t *opt, double *P)
{
  int i;
  double var = 0.0;
  for (i = xiP(); i < xiP() + 3; i++){
    var += SQRT(P[i + i * nx]);
  }

  printf("checkpcov: var summation: %lf\n", var);

  if ((var / 3) > 100){
    if (P){
      printf("checkpcov: ok\n");
      getP0(opt, P, nx);
    }
  }
}

/* imu body position transform to gps antenna---------------------------------*/
extern void insp2antp(const ins_states_t *ins, insgnss_opt_t *insopt, double *rr)
{
    int i; double T[3];
    matmul3v("N",ins->Cbe,insopt->lever,T);
    for (i=0;i<3;i++) rr[i]=ins->re[i]+T[i];
}
/* initialize state and covariance -------------------------------------------*/
static void tcinitx(ins_states_t *ins,double xi, double var, int i)
{
    int j;
    ins->x[i]=xi;
    for (j=0;j<ins->nx;j++) {
        ins->P[i+j*ins->nx]=ins->P[j+i*ins->nx]=i==j?var:0.0;
    }
} 
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk, int nx, ins_states_t *ins)
{
    double *F,*P,*FP,*x,*xp,pos[3],Q[9]={0},Qv[9];
    int i,j,*ix;
    
    trace(3,"udpos_ppp:\n");
    printf("udpos_ppp:\n");
    
    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) tcinitx(ins,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(ins->x+xiP(),3)<=0.0) {
      printf("pos.ini.onetime:\n");
        for (i=0;i<3;i++) tcinitx(ins,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) tcinitx(ins,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) tcinitx(ins,1E-6,VAR_ACC,i);
        }
    }
    /* static ppp mode */
    if (rtk->opt.mode==PMODE_PPP_STATIC) {
        for (i=0;i<3;i++) {
            rtk->P[i*(1+nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt);
        }
        return;
    }
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
      printf("pos.ini.onetime: Without dynamics?\n");
        for (i=0;i<3;i++) {
            tcinitx(ins,rtk->sol.rr[i],VAR_POS,i);
        }
        return;
    }
    printf("Other options:\n");
    /* generate valid state index */
    ix=imat(nx,1);
    for (i=nx=0;i<nx;i++) {
        if (rtk->x[i]!=0.0&&rtk->P[i+i*nx]>0.0) ix[nx++]=i;
    }
    if (nx<9) {
        free(ix);
        return;
    }
    printf("Passed here?\n");
    /* state transition of position/velocity/acceleration */
    F=eye(nx); P=mat(nx,nx); FP=mat(nx,nx); x=mat(nx,1); xp=mat(nx,1);
    
    for (i=0;i<6;i++) {
        F[i+(i+3)*nx]=rtk->tt;
    }
    for (i=0;i<3;i++) {
        F[i+(i+6)*nx]=SQR(rtk->tt)/2.0;
    }
    for (i=0;i<nx;i++) {
        x[i]=rtk->x[ix[i]];
        for (j=0;j<nx;j++) {
            P[i+j*nx]=rtk->P[ix[i]+ix[j]*nx];
        }
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",nx,1,nx,1.0,F,x,0.0,xp);
    matmul("NN",nx,nx,nx,1.0,F,P,0.0,FP);
    matmul("NT",nx,nx,nx,1.0,FP,F,0.0,P);
    
    for (i=0;i<nx;i++) {
        rtk->x[ix[i]]=xp[i];
        for (j=0;j<nx;j++) {
            rtk->P[ix[i]+ix[j]*nx]=P[i+j*nx];
        }
    }
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(rtk->tt);
    Q[8]=SQR(rtk->opt.prn[4])*fabs(rtk->tt);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+6+(j+6)*nx]+=Qv[i+j*3];
    }
    free(ix); free(F); free(P); free(FP); free(x); free(xp);
}

/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk, ins_states_t *ins)
{
    double dtr;
    int i;
    
    trace(3,"udclk_ppp:\n");
    printf("udclk_ppp:\n");
    
    /* initialize every epoch for clock (white noise) */
    for (i=0;i<NSYS;i++) {
        if (rtk->opt.sateph==EPHOPT_PREC) {
            /* time of prec ephemeris is based gpst */
            /* negelect receiver inter-system bias  */
            dtr=rtk->sol.dtr[0];
        }
        else {
            dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];
        }
        printf("udclk here again?: dtr:%lf\n", CLIGHT*dtr);
        /* Only parameter update */
        //ins->x[xiRc()]=CLIGHT*dtr;
        tcinitx(ins,CLIGHT*dtr,VAR_CLK,xiRc());
    }

    /* initialize GLONASS */
    if(rtk->opt.navsys>5){
      dtr=rtk->sol.dtr[0]+rtk->sol.dtr[1];
      tcinitx(ins,CLIGHT*dtr,VAR_CLK,xiRc()+1);
    }

      // dtr=rtk->sol.dtr[0]+rtk->sol.dtr[1];
      // tcinitx(ins,CLIGHT*dtr,VAR_CLK,xiRc()+1);
    
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk, int nx, ins_states_t *ins)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=xiTr(&rtk->opt),j;
    
    trace(3,"udtrop_ppp:\n");
    printf("udtrop_ppp:\n");
    
    if (ins->x[i]==0.0) {
      
        ecef2pos(rtk->sol.rr,pos);
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);
        printf("If state is zero: %lf\n", ztd);
        tcinitx(ins,ztd,var,i);
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) tcinitx(ins,1E-6,VAR_GRA,j);
        }
    }
    else {
      printf("If state is not zero:\n");
        ins->P[i+i*nx]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt);
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                rtk->P[j+j*nx]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);
            }
        }
    }
}

/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, int nx, 
                        ins_states_t *ins)
{
    const double *lam;
    double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
    double ion,dantr[NFREQ]={0},dants[NFREQ]={0};
    int i,j,k,l,f,sat,slip[MAXOBS]={0},clk_jump=0;
    
    trace(3,"udbias  : n=%d\n",n);
    
    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
        rtk->ssat[i].slip[j]=0;
    }
    /* detect cycle slip by LLI */
    detslp_ll(rtk,obs,n);
    
    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk,obs,n,nav);
    
    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(rtk,obs,n,nav);
    
    ecef2pos(rtk->sol.rr,pos);
    
    for (f=0;f<1;f++) {
        
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout||
                rtk->opt.modear==ARMODE_INST||clk_jump) {
                tcinitx(ins,0.0,0.0,xiBs(&rtk->opt, i+1));
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=xiBs(&rtk->opt, sat);
            corr_meas(obs+i,nav,rtk->ssat[sat-1].azel,&rtk->opt,dantr,dants,
                      0.0,L,P,&Lc,&Pc);
            
            bias[i]=0.0;
            
            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                bias[i]=Lc-Pc;
                slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
            }
            else if (L[f]!=0.0&&P[f]!=0.0) {
                slip[i]=rtk->ssat[sat-1].slip[f];
                l=satsys(sat,NULL)==SYS_GAL?2:1;
                lam=nav->lam[sat-1];
                if (obs[i].P[0]==0.0||obs[i].P[l]==0.0||
                    lam[0]==0.0||lam[l]==0.0||lam[f]==0.0) continue;
                ion=(obs[i].P[0]-obs[i].P[l])/(1.0-SQR(lam[l]/lam[0]));
                bias[i]=L[f]-P[f]+2.0*ion*SQR(lam[f]/lam[0]);
            }
            if (ins->x[j]==0.0||slip[i]||bias[i]==0.0) continue;
            
            offset+=bias[i]-ins->x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherency */
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<MAXSAT;i++) {
                j=xiBs(&rtk->opt, i+1);
                if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k;
            }
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
        }
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=xiBs(&rtk->opt,sat);
            
            ins->P[j+j*nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);
            
            if (bias[i]==0.0||(ins->x[j]!=0.0&&!slip[i])) continue;
            
            /* reinitialize phase-bias if detecting cycle slip */
            tcinitx(ins,bias[i],VAR_BIAS,xiBs(&rtk->opt,sat));
            
            /* reset fix flags */
            for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0;
            
            trace(5,"udbias_ppp: sat=%2d bias=%.3f\n",sat,bias[i]);
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, int nx, 
                        ins_states_t *ins)
{
    trace(3,"udstate_ppp: n=%d\n",n);
    
    /* temporal update of position */
    //udpos_ppp(rtk,nx,ins);
    
    /* temporal update of clock */
    udclk_ppp(rtk, ins);
    //ins->x[i]=xi;dtr=rtk->sol.dtr[0];
    
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(rtk,nx,ins);
    }
    /* temporal update of phase-bias */
    udbias_ppp(rtk,obs,n,nav, nx, ins);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;
    
    trace(3,"testeclipse:\n");
    
    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);
    
    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;
        
        if ((r=norm(rs+i*6,3))<=0.0) continue;
        
        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;
        
        /* sun-earth-satellite angle */
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);
        
        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;
        
        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);
        
        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int type, const prcopt_t *opt)
{
    double a,b,a2,b2,fact=1.0;
    double sinel=sin(el);
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0);

    /* extended error model */
    if (type==1&&opt->exterr.ena[0]) { /* code */
        a=opt->exterr.cerr[i][0];
        b=opt->exterr.cerr[i][1];
        if (opt->ionoopt==IONOOPT_IFLC) {
            a2=opt->exterr.cerr[i][2];
            b2=opt->exterr.cerr[i][3];
            a=sqrt(SQR(2.55)*a*a+SQR(1.55)*a2*a2);
            b=sqrt(SQR(2.55)*b*b+SQR(1.55)*b2*b2);
        }
    }
    else if (type==0&&opt->exterr.ena[1]) { /* phase */
        a=opt->exterr.perr[i][0];
        b=opt->exterr.perr[i][1];
        if (opt->ionoopt==IONOOPT_IFLC) {
            a2=opt->exterr.perr[i][2];
            b2=opt->exterr.perr[i][3];
            a=sqrt(SQR(2.55)*a*a+SQR(1.55)*a2*a2);
            b=sqrt(SQR(2.55)*b*b+SQR(1.55)*b2*b2);
        }
    }
    else { /* normal error model */
        if (type==1) fact*=opt->eratio[0];
        fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
        if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
        a=fact*opt->err[1];
        b=fact*opt->err[2];
    }
    return a*a+b*b/sinel/sinel;
}


/* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, rtk_t *rtk, double *v, double *H, double *R,
                   double *azel,double *rpos, insgnss_opt_t *insopt,ins_states_t *insc)
{
  prcopt_t *opt=&rtk->opt;
    double r,rr[3],disp[3],pos[3],e[3],meas[2],dtdx[3],dantr[NFREQ]={0};
    double dants[NFREQ]={0},var[MAXOBS*2],dtrp=0.0,vart=0.0,varm[2]={0};
    int i,j,k,sat,sys,nv=0,nx=insc->nx,brk,tideopt;

    printf("res_ppp : n=%d nx=%d\n",n,nx);

    for (i=0;i<MAXSAT;i++) rtk->ssat[i].vsat[0]=0;

    if (insopt->Nav_or_KF){
      for (i=0;i<3;i++) rr[i]=insc->re[i];
    }else{
      for (i=0;i<3;i++) rr[i]=rtk->x[i];
      }
    for (i=0;i<3;i++) rr[i]=insc->re[i];
    //for (i=0;i<3;i++) rr[i]=rtk->sol.rr[i];

    /* earth tides correction */
    if (opt->tidecorr) {
        tideopt=1; /* 1:solid, 2:solid+otl+pole */
        tidedisp(gpst2utc(obs[0].time),rr,tideopt,&nav->erp,opt->odisp[0],
                 disp);
        for (i=0;i<3;i++) rr[i]+=disp[i];
    }
    ecef2pos(rr,pos);

    for (i=0;i<n&&i<MAXOBS;i++) {
      sat=obs[i].sat;

        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs) continue;

        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) continue;

        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* tropospheric delay correction */
        if (opt->tropopt==TROPOPT_SAAS) {
            dtrp=tropmodel(obs[i].time,pos,azel+i*2,REL_HUMI);
            vart=SQR(ERR_SAAS);
            printf("Tropo SAAST: %lf\n", dtrp);
        }
        else if (opt->tropopt==TROPOPT_SBAS) {
            dtrp=sbstropcorr(obs[i].time,pos,azel+i*2,&vart);
            printf("Tropo SBAS: %lf\n", dtrp);
        }
        else if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
            dtrp=prectrop(obs[i].time,pos,azel+i*2,opt,x+xiTr(opt),dtdx,&vart);
            printf("pos=%lf, %lf, %lf, azel=%6.1f x[%d]:%lf dtdx=%lf vart=%lf\n",
          pos[0]*R2D, pos[1]*R2D,pos[2],azel[i*2]*R2D,xiTr(opt),x[xiTr(opt)],dtdx,vart);
            printf("Tropo EST: %lf\n", dtrp);
        }
        else if (opt->tropopt==TROPOPT_COR||opt->tropopt==TROPOPT_CORG) {
            dtrp=prectrop(obs[i].time,pos,azel+i*2,opt,x,dtdx,&vart);
            printf("Tropo PRECTROP: %lf\n", dtrp);
        }
        /* satellite antenna model */
        if (opt->posopt[0]) {
            satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants);
        }
        /* receiver antenna model */
        antmodel(opt->pcvr,opt->antdel[0],azel+i*2,opt->posopt[1],dantr);

        /* phase windup correction */
        if (opt->posopt[2]) {
            windupcorr(rtk->sol.time,rs+i*6,rr,&rtk->ssat[sat-1].phw);
        }
        /* ionosphere and antenna phase corrected measurements */
        if (!corrmeas(obs+i,nav,pos,azel+i*2,&rtk->opt,dantr,dants,
                      rtk->ssat[sat-1].phw,meas,varm,&brk)) {
            continue;
        }

        /* satellite clock and tropospheric delay */
        r+=-CLIGHT*dts[i*2]+dtrp;

        printf("sat=%2d azel=%6.1f %5.1f dtrp=%.3f dantr=%6.3f %6.3f dants=%6.3f %6.3f phw=%6.3f\n",
              sat,azel[i*2]*R2D,azel[1+i*2]*R2D,dtrp,dantr[0],dantr[1],dants[0],
              dants[1],rtk->ssat[sat-1].phw);


     for (j=0;j<2;j++) { /* for phase and code */

         if (meas[j]==0.0) continue;

        printf("xip: %d, NP: %d, IP:%d\n", \
         xiP(),NP, IP);

         for (k=0;k<nx;k++) H[k+nx*nv]=0.0;

         v[nv]=meas[j]-r;
         printf("RES 0: %lf\n", v[nv]);

         for (k=xiP();k<NP+IP;k++) H[k+nx*nv]=-e[k-xiP()];

         if (sys!=SYS_GLO) {
             v[nv]-=x[xiRc()];
             H[xiRc()+nx*nv]=1.0;
         }
         else {
           
             v[nv]-=x[xiRc()+1];  //GLONASS
             H[xiRc()+1+nx*nv]=1.0;
             printf("RES GLO: %lf\n", v[nv]);
         }
         printf("RES 1: %lf\n", v[nv]);
         if (opt->tropopt>=TROPOPT_EST) {
             for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                 H[xiTr(opt)+k+nx*nv]=dtdx[k];
             }
         }
         if (j==0) {
             v[nv]-=x[xiBs(opt,obs[i].sat)];
             H[xiBs(opt,obs[i].sat)+nx*nv]=1.0;
         }
         printf("RES 2: %lf\n", v[nv]);

         var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],j,opt)+varm[j]+vare[i]+vart;

         printf("sat: %d, var[%d]: %lf, sys: %d, azel: %lf, varm: %lf, vare: %lf, vart: %lf\n", \
         obs[i].sat, nv, var[nv], sys, azel[1+i*2], varm[j], vare[i], vart);

         if (j==0) rtk->ssat[sat-1].resc[0]=v[nv];
         else      rtk->ssat[sat-1].resp[0]=v[nv];

         printf("Meas: %lf Modeled: %lf RES.: %lf\n", meas[j], r, v[nv]);

         /* test innovation */
#if 0
         if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
#else
         if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno&&sys!=SYS_GLO) {
#endif
             trace(2,"ppp outlier rejected %s sat=%2d type=%d v=%.3f\n",
                   time_str(obs[i].time,0),sat,j,v[nv]);
             rtk->ssat[sat-1].rejc[0]++;
             continue;
         }
         
         if (j==0) rtk->ssat[sat-1].vsat[0]=1;
         nv++;

     } // Phase and code loop (j)

     /* KF residuals output */
     fprintf(out_KF_residuals, "%lf %2d %lf %lf\n", time2gpst(obs[i].time, NULL),
            sat, v[nv-2], v[nv-1]);

   } //sat loop (i)


    for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
        R[i+j*nv]=i==j?var[i]:0.0;
    }
    trace(5,"x=\n"); tracemat(5,x, 1,nx,8,3);
    trace(5,"v=\n"); tracemat(5,v, 1,nv,8,3);
    trace(5,"H=\n"); tracemat(5,H,nx,nv,8,3);
    trace(5,"R=\n"); tracemat(5,R,nv,nv,8,5);
    return nv;
}

/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat, ins_states_t *ins)
{
    const prcopt_t *opt=&rtk->opt;
    double *xa,*Pa,*x,*P,rr[3];
    int i,j,tc,ip,na,nx,ivx,ivy,ivz,idtr, idtrr;

    nx=ins->nx;
    x=ins->x;
    P=ins->P;

    /* index of rover station position, velocity, and clock states */
    ip=xiP();
    ivx=xiV()+0;
    ivy=xiV()+1;
    ivz=xiV()+2;
    idtr=xiRc(&opt);
    idtrr=xiRr(&opt);
    
    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
            rtk->ssat[obs[i].sat-1].lock[j]++;
            rtk->ssat[obs[i].sat-1].outc[j]=0;
            if (j==0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;


      //   matcpy(rr,ins->re,1,3);
        
      //   for (i=0;i<3;i++) {
      //       rtk->sol.rr[i]=rr[i];
      //       rtk->sol.qr[i]=(float)P[i+ip+(i+ip)*nx];
      //   }
      //   rtk->sol.qr[3]=(float)P[1+ip];
      //   rtk->sol.qr[4]=(float)P[2+ip+ins->nx];
      //   rtk->sol.qr[5]=(float)P[2+ip];

      //   /* velocity and covariance */
      //   //if (rtk->opt.dynamics) {

      //       matcpy(rr,ins->ve,1,3);

      //       for (i=0;i<3;i++) rtk->sol.rr[3+i]=rr[i];

      //       rtk->sol.qrv[0]=(float)P[ivx+ivx*nx];
      //       rtk->sol.qrv[1]=(float)P[ivy+ivy*nx];
      //       rtk->sol.qrv[2]=(float)P[ivz+ivz*nx];

      //       rtk->sol.qrv[3]=(float)P[ivx+ivy*nx];
      //       rtk->sol.qrv[4]=(float)P[ivy+ivz*nx];
      //       rtk->sol.qrv[5]=(float)P[ivz+ivx*nx];
      //  // }
   
        /* clock offsets and drift */
        rtk->sol.dtr[0]=x[idtr]/CLIGHT; /* GPS in seconds */
        ins->dtr[0]=x[idtr]; /* in meters */
        rtk->sol.dtr[1]=(x[idtr+1]-x[idtr])/CLIGHT; /* GLO-GPS in seconds */
        ins->dtr[1]=x[idtr+1];  /* in meters */
       // ins->dtrr=x[idtrr];
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr[j]=obs[i].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
    }
}
/* precise point positioning -------------------------------------------------*/
extern int pppos1(rtk_t *rtk, const obsd_t *obs, ins_states_t *insp,
                    insgnss_opt_t *insopt, int n, const nav_t *nav)
{
    const prcopt_t *opt=&rtk->opt;
    double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,dr[3]={0},std[3];
    double *x,*P,rr[3];
    char str[32];
    int i,j,nv,info=0,svh[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE,tc;
    int nx=insp->nx;
    
    
    time2str(obs[0].time,str,2);
    printf("pppos   : time=%s nx=%d n=%d\n",str,insp->nx,n);
    printf("pppos inputs: ins.nx=%d, rtk->nx=%d, nsat:%d time: %lf\n",insp->nx,rtk->nx, n, insp->time);
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;

    /* temporal update of ekf states */
    udstate_ppp(rtk,obs,n,nav, insp->nx, insp);

    printf("GLonass clock: %lf\n", insp->dtr[1]);
  
    /* satellite positions and clocks */
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
    
    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        testeclipse(obs,n,nav,rs);
    }
    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rtk->x,opt->tidecorr==1?1:7,&nav->erp,
                 opt->odisp[0],dr);
    }
    nv=insp->nx*2;//rtk->opt.nf*2;
    xp=zeros(insp->nx,1); Pp=zeros(insp->nx,insp->nx);
    v=mat(nv,1); H=mat(insp->nx,nv); R=mat(nv,nv);

    printf("ppp observations: nv=%d, parameters=%d\n",nv, insp->nx);
    printf("Before PPP integration:\n");
  printf("x vector:\n");
  for (j = 0; j < nx; j++) printf("%lf ", insp->x[j]); 
/*
  printf("\nP\n");
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++)
    {
      (i==j?printf("%.15lf ", insp->P[i*nx + j]):0);
    }
  }
  printf("\n");
  */

    //x=insp->x;
    //P=insp->P;
    matcpy(xp,insp->x,nx,1);

    /* ins/gnss filter iteration */
    for (i=0;i<MAX_ITER;i++) {

        insp2antp(insp,insopt,rr);        
        //matcpy(xp,x,nx,1);
        matcpy(Pp,insp->P,nx,nx);
/*
          printf("\nPp before\n");
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++)
    {
      (i==j?printf("%.15lf ", Pp[i*nx + j]):0);
    }
  }
  printf("\n");
     */
        /* prefit residuals */
        if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel,rr,insopt,insp))) {
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);
            break;
        }
        printf("PPP nv: %d \n", nv);
        /* measurement update of ekf states */
        if ((info=filter(xp,Pp,H,v,R,nx,nv))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            info=0;
            break;
        }
     /* printf("H\n"); */
  for (i = 0; i < nv; i++)
  {
    for (j = 0; j < insp->nx; j++)
    {
      printf("%lf ", H[i * insp->nx + j]); //i * insp->nx + j , j * nv + i
    }
    printf("\n");
  }
/*
    printf("R\n");
  for (i = 0; i < nv; i++)
  {
    for (j = 0; j < nv; j++)
    {
      //printf("%lf ", R[j * insp->nx + i]);
      (i==j?printf("%lf ", R[j*nv + i]):0);
    }
    printf("\n");
  }
  */
  printf("v\n");
  for (i = 0; i < nv; i++) printf("%lf ", v[i]);
  printf("\n");
    printf("x vector after fisrt\n");
  for (j = 0; j < nx; j++) printf("%.15lf ", xp[j]); 
  printf("\nPp after\n");
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++)
    {
      (i==j?printf("%.15lf ", Pp[i*nx + j]):0);
    }
  }
  printf("\n");

        /* postfit residuals */
        if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel,rr,insopt,insp)) {
             printf("Postfit ok:\n");
            /* update state and covariance matrix */
            matcpy(insp->x,xp,nx,1);
            matcpy(insp->P,Pp,nx,nx);

            stat=SOLQ_PPP;
            insopt->Nav_or_KF=1;
      
            clp(insp,insopt,xp);
            for (j=0;j<xnCl();j++) xp[j]=0.0;
          
            break;
        }
    }
  printf("After PPP integration:\n");
 /* printf("H\n");
 
  for (i = 0; i < nv; i++)
  {
    for (j = 0; j < insp->nx; j++)
    {
      printf("%lf ", H[i * insp->nx + j]); //i * insp->nx + j , j * nv + i
    }
    printf("\n");
  }*/
/*
    printf("R\n");
  for (i = 0; i < nv; i++)
  {
    for (j = 0; j < nv; j++)
    {
      //printf("%lf ", R[j * insp->nx + i]);
      (i==j?printf("%lf ", R[j*nv + i]):0);
    }
    printf("\n");
  }*/
/*
  printf("v\n");
  for (i = 0; i < nv; i++) printf("%lf ", v[i]);
  printf("\n");
*/
  printf("x vector:\n");
  for (j = 0; j < nx; j++) printf("%.15lf ", xp[j]); 
/*
  printf("\nPp\n");
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++)
    {
      (i==j?printf("%.15lf ", Pp[i*nx + j]):0); 
    }
  }*/
  printf("\n");
  printf("\nP\n");
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < nx; j++)
    {
      (i==j?printf("%.15lf ", insp->P[i*nx + j]):0); 
    }
  }
  printf("\n");
  
    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",str,i);
        printf("%s ppp (%d) iteration overflows\n",str,i);
    }
    if (stat==SOLQ_PPP) {
      info=1; 
        printf("ppp solution update");
        insp2antp(insp,insopt,rr);

        /* update solution status */
        update_stat(rtk,obs,n,stat, insp);

    }
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);

    printf("ppp solution info: %d\n",info);
    return info;
}
/* converts a coordinate transformation matrix to the corresponding set of
 * Euler angles -------------------------------------------------------------
 * args    : double *Cnb      I matrix of ned-frame to body-frame
 *           double *rpy      O attitude {roll,pitch,yaw}
 * return  : none
 * ------------------------------------------------------------------------*/
extern void dcm2rpy(const double *Cnb,double *rpy)
{
    rpy[0]= atan2(Cnb[7],Cnb[8]);  /* roll */
    rpy[1]=-asin (Cnb[6]);         /* pitch */
    rpy[2]= atan2(Cnb[3],Cnb[0]);  /* yaw */
}
/* standard deviation of value series-----------------------------------------------
 * args    : double *val    I  value series
 *           int     n      I  number of value
 * return  : standard Deviation
 * ---------------------------------------------------------------------------*/
extern double stds(const double *val,int n)
{
    int i;
    double mv,std,*vv=mat(n,1);

    for (mv=0.0,i=0;i<n;i++) mv+=val[i]; mv/=n;
    for (i=0;i<n;i++) vv[i]=val[i]-mv;
    
    std=sqrt(dot(vv,vv,n)/n); free(vv); return std;
}
/* check vehicle whether is straight driving --------------------------------*/
extern int chksdri(const double *vel,int n)
{
    int i;
    double *head=mat(n,1),hstd=0.0;

    for (i=0;i<n;i++) {
        //head[i]=vel2head(vel+3*i);
        head[i]=atan2(vel[i*3+1],fabs(vel[i*3+0])<1E-4?1E-4:vel[i*3+0]);
      }
    hstd=stds(head,n);
    free(head); return hstd<8.0*D2R;
}
/* normalize angle-----------------------------------------------------------*/
extern double NORMANG(double ang)
{
    while (ang<0.0) ang+=360.0; return ang;
}
/* get attitude from ins states----------------------------------------------
 * args   :  insstate_t *ins  I  ins states
 *           double *rpy      O  ins attitude
 * return : none
 * --------------------------------------------------------------------------*/
extern void getatt(const ins_states_t *ins,double *rpy)
{
    double llh[3],C[9],Cnb[9];

    ecef2pos(ins->re,llh);
    ned2xyz(llh,C);
    matmul("TN",3,3,3,1.0,ins->Cbe,C,0.0,Cnb);
    dcm2rpy(Cnb,rpy);
}

/* re-check attitude---------------------------------------------------------
 * args   :  insstate_t *ins  IO  ins states
 *           imud_t    *imu   I   imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int rechkatt(ins_states_t *ins, const imuraw_t *imu)
{
    int NPOS=3;//insgnssopt.gnssw;
    int i,j;
    double dt,*vel,llh[3];
    double C[9],yaw,vn[3],rpy[3]={0};
    double vb[3],pvb[3],Cbe[9];

    vel=zeros(3,NPOS);

    trace(3,"rechkatt:\n");
    printf("rechkatt:\n");

    /* check gps solution status */
    for (i=0;i<NPOS;i++) {
        if (solw[i].stat==0) return 0;
    }
    /* recheck ins attitude if need */

    if (gnss_w_counter>=insgnssopt.gnssw) {

        /* velocity for trajectory*/
        if (norm(solw[i-1].rr+3, 3)){
          /* Velocity from solution */
          for (i=NPOS;i>0;i--) {
              for (j=0;j<3;j++) {
                vel[3*(NPOS-i)+j]=solw[insgnssopt.gnssw-i].rr[j+3];
               }
           }
         }else{
              /* Velocity from position */
              for (i=NPOS;i>=2;i--) {
                if ((dt=timediff(solw[i-1].time,solw[i-2].time))>3.0
                  ||fabs(dt)<=1E-5) {
                  continue;
                }
                for (j=0;j<3;j++) {
                 vel[3*(NPOS-i)+j]=(solw[i-1].rr[j]-solw[i-2].rr[j])/dt;
                }   
              }
           }

        
        for (i = 0; i < NPOS; i++)
        {
          printf("Solw.v: %lf %lf %lf \n", solw[i].rr[3],solw[i].rr[4],solw[i].rr[5] );
          printf("Vel: %lf %lf %lf \n", vel[3*i+0],vel[3*i+1],vel[3*i+2] );
        }
        
        /* check velocity whether is straight driving  */
        if (!chksdri(vel,NPOS-1)) {
            trace(2,"no straight driving\n");
            printf("no straight driving\n");
            return 0;
        }
        printf("straight driving\n");
        /* check velocity */
        if (norm(vel,3)>MAXVEL
            &&norm(imu->wibb,3)<MAXGYRO) {
              printf("Velocity and gyro turn check ok\n");

            /* velocity convert to attitude */
            ecef2pos(solw[NPOS-1].rr,llh);
            ned2xyz(llh,C);

            /* yaw */
            matmul("TN",3,1,3,1.0,C,vel,0.0,vn);
            yaw=NORMANG(vel2head(vn)*R2D);

            /* attitude for current ins states */
            getatt(ins,rpy);

            if (fabs(yaw-NORMANG(rpy[2]*R2D))<MAXANG) return 0;
#if 0
            rpy[2]=yaw*D2R; /* reset yaw */
#else
            rpy[2]=(NORMANG(rpy[2]*R2D)+yaw)/2.0*D2R;
#endif
            rpy2dcm(rpy,C);
            matt(C,3,3,ins->Cbn);

            ned2xyz(llh,C);
            matmul("NN",3,3,3,1.0,C,ins->Cbn,0.0,Cbe);
            matmul("TN",3,1,3,1.0,Cbe,ins->ve,0.0,vb);
            matmul("TN",3,1,3,1.0,ins->Cbe,ins->ve,0.0,pvb);

            /* check again */
            if (fabs(norm(vb,3)-norm(pvb,3))<MINVEL
                &&(fabs(vb[1])<fabs(pvb[1]))
                &&(fabs(vb[2])<fabs(pvb[2]))) {
                matcpy(ins->Cbe,Cbe,3,3);
                trace(3,"recheck attitude ok\n");
                return 1;
            }
        }
    }
    free(vel);
    trace(3,"no recheck attitude\n");
    return 0;
}
/* update ins states in n-frame----------------------------------------------*/
extern void update_ins_state_n(ins_states_t *ins)
{
    double Cne[9];

    /* position */
    ecef2pos(ins->re,ins->rn);

    /* attitude/velocity */
    ned2xyz(ins->rn,Cne);
    matmul("TN",3,3,3,1.0,Cne,ins->Cbe,0.0,ins->Cbn);
    /* get attitude in _nb from Cbe */
    getatt(ins,ins->an);
    matmul("TN",3,1,3,1.0,Cne,ins->ve ,0.0,ins->vn );

    /* acceleration */
    matmul("TN",3,1,3,1.0,Cne,ins->data.fbe,0.0,ins->data.fbn);
}
/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
//int main (void){
//int InsGnssCore (){ LATER, WHEN LINKED TO OTHER MAIN FUNCTION USE IT THIS WAY
extern int TC_INS_GNSS_core1(rtk_t *rtk, const obsd_t *obs, int n, nav_t *nav,
                             ins_states_t *insc, insgnss_opt_t *ig_opt,
                             int nav_or_int)
{
  int i,j, nx = insc->nx;
  double dt, gnss_time;

  printf("\n *****************  TC INSGNSS CORE BEGINS ***********************\n");

  prcopt_t *opt = &rtk->opt;

  //nx=xnRx(opt)+n;

  printf("Number of fixed par.: %d, Nsat:%d, Sum: %d parameters\n",xnRx(opt), n, nx );

  //insc->nx=nx;

  gnss_time = time2gpst(rtk->sol.time, NULL);

  /* Initialize KF parameters uncertainty */
  kf_par_unc_init(ig_opt);

  /* Initialize KF noise information */
  kf_noise_init(ig_opt);

  /* Initialize ins process noise indices */
  initPNindex(opt);

  /* Ins navigation */
  Nav_equations_ECEF1(insc);

       printf("ins->P before checkpcov\n");
  for (i = 0; i < insc->nx; i++)
  {
    for (j = 0; j < insc->nx; j++)
    {
      (i==j?printf("%lf ", insc->P[i * insc->nx + j]):0);
    }
  }
  printf("\n GNSS time: %lf", gnss_time);

  /* propagate ins states */
  printf("Prop ins\n");
  propinss(insc, ig_opt, insc->dt, insc->x, insc->P);
  printf("Prop ins end\n");


  /* Checking input values  */
  chkpcov(nx, ig_opt, insc->P);

  if (nav_or_int)
  {
    /* Integration */
    if (obs && insc && n){
      /* code */
      dt = fabs(gnss_time - insc->time);

      /* check synchronization */
      if (fabs(dt) > 3.0)
      {
        trace(2, "observation and imu sync error\n");
        ig_opt->Nav_or_KF = 0;
        return (0);
      }

      /* tightly coupled */
     if(ig_opt->Nav_or_KF=pppos1(rtk, obs, insc, ig_opt, n, nav)){
       printf("Tc ins/gnss integrated ok\n");
     }else printf("Tc ins/gnss integrated fail\n");

    }
  }else{
    /* Navigate only */
    ig_opt->Nav_or_KF = 0;

    /* ins state in n-frame */
      update_ins_state_n(insc);
  }

  /* Update integrated solution */
  if(ig_opt->Nav_or_KF){
    insc->ptctime=insc->time;

    /* recheck attitude */
      if(rechkatt(insc,&insc->data)){
        printf("rechecked attitude ok\n");
      } else printf("no rechecked attitude\n");

    /* ins state in n-frame */
      update_ins_state_n(insc);
  }

  printf("OUTPUT OF INTEGRATED OR NAVIGATED SOLUTION: %lf\n", insc->time);
  printf("P: %lf, %lf, %lf\n", insc->re[0], insc->re[1], insc->re[2]);
  printf("llh: %lf, %lf, %lf\n", insc->rn[0]*R2D, insc->rn[1]*R2D, insc->re[0]);
  printf("V: %lf, %lf, %lf\n", insc->ve[0], insc->ve[1], insc->ve[2]);
  printf("Ba=aft[%lf; %lf; %lf] t:%lf, %d\n", insc->data.ba[0], insc->data.ba[1], insc->data.ba[2],insc->time, ig_opt->Nav_or_KF);
  printf("Bg=aft[%lf, %lf, %lf]t:%lf, %d\n", insc->data.bg[0], insc->data.bg[1], insc->data.bg[2],insc->time, ig_opt->Nav_or_KF);
  printf("ins->Cbe\n");
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      printf("%lf ", insc->Cbe[i * 3 + j]);
    }
    printf("\n");
  }/*
   printf("ins->P\n");
  for (i = 0; i < insc->nx; i++)
  {
    for (j = 0; j < insc->nx; j++)
    {
      (i==j?printf("%lf ", insc->P[i * insc->nx + j]):0);
    }
  }
  printf("\n");*/

  /* Prepare GNSS and INS raw data into the proper structures --------------- */

  /* Plots -------------------------------------------------------------------*/

  /* Write output profile and errors file ------------------------------------*/

  /* Free memory -------------------------------------------------------------*/


  return 1;

  printf("\n *****************  TC INSGNSS CORE ENDS *************************\n");
}

/* Name of function ------------------------------------------------------------
* Brief description
* arguments  :
* datatype  name  I/O   description
* int        a     I   describe a (a unit)
* double    *b     O   describe b (b unit) {b components x,y,z}
*
* return : what does it return?
* notes  :
*-----------------------------------------------------------------------------*/
extern int LC_INS_GNSS_core1(rtk_t *rtk, const obsd_t *obs, int n, nav_t *nav,
                             ins_states_t *insc, insgnss_opt_t *ig_opt, int nav_or_int)
{

  /* Plots -------------------------------------------------------------------*/

  /* Write output profile and errors file ------------------------------------*/

  /* Free memory -------------------------------------------------------------*/

  printf("\n *****************  LCINSGNSS CORE ENDS *************************\n");
  return 0;
}
