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


/* struct declaration -------------------------------------------------------*/

/* Utility ----------------------------------------------------------*/
static inline int8_t sign(double val) {
  /* Corresponding sign() function in Matlab
  From: https://forum.arduino.cc/index.php?topic=37804.0 solution*/
 if (val < 0.0) return -1;
 if (val > 0.0) return 1;
 return 0;
}

void Initialize_INS_GNSS_LCKF_consumer_grade_IMU(initialization_errors *initialization_errors,\
  IMU_errors *IMU_errors, GNSS_config *GNSS_config, LC_KF_config *LC_KF_config){
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
double rad_to_deg = 1/deg_to_rad;
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
IMU_errors->M_a[0] = 50000* 1E-6; IMU_errors->M_a[1] = -15000* 1E-6;
IMU_errors->M_a[2] = 10000* 1E-6;
IMU_errors->M_a[3] =-7500* 1E-6; IMU_errors->M_a[4] = -60000* 1E-6;
IMU_errors->M_a[5] = 12500* 1E-6;
IMU_errors->M_a[6] =-12500* 1E-6; IMU_errors->M_a[7] = 5000* 1E-6;
IMU_errors->M_a[8] = 20000 * 1E-6;
/* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
IMU_errors->M_g[0] = 40000* 1E-6; IMU_errors->M_a[1] = -14000* 1E-6;
IMU_errors->M_g[2] = 12500* 1E-6;
IMU_errors->M_g[3] = 0* 1E-6; IMU_errors->M_a[4] = -30000* 1E-6;
IMU_errors->M_g[5] = -7500* 1E-6;
IMU_errors->M_g[6] = 0* 1E-6; IMU_errors->M_a[7] = 0* 1E-6;
IMU_errors->M_g[8] = -17500 * 1E-6;

/* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
IMU_errors->G_g[0] = 90 * deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[1] = -110* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[2] = -60* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[3] = -50* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[4] = 190* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[5] = -160* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[6] = 30* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[7] = 110* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[8] = -130* deg_to_rad / (3600 * 9.80665);

/* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) */
IMU_errors->accel_noise_root_PSD = 1000 * micro_g_to_meters_per_second_squared;
/* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) FROM um7 datasheet */
IMU_errors->accel_noise_root_PSD = 400 * micro_g_to_meters_per_second_squared;
/* Gyro noise root PSD (deg per root hour, converted to rad s^-0.5) */
IMU_errors->gyro_noise_root_PSD = 1 * deg_to_rad / 60.0;
/* Gyro noise root PSD (deg per second root Hz, converted to rad s^-0.5)  From um7 datasheet */
IMU_errors->gyro_noise_root_PSD = 0.005 * deg_to_rad ;
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
for (i=0;i<3;i++) LC_KF_config->init_att_unc[i] = D2R*2;
/* Initial velocity uncertainty per axis (m/s) */
for (i=0;i<3;i++) LC_KF_config->init_vel_unc[i] = 0.1;
/* Initial position uncertainty per axis (m) */
for (i=0;i<3;i++) LC_KF_config->init_pos_unc[i] = 10;
/* Initial accelerometer bias uncertainty per instrument (micro-g, converted
/* to m/s^2) */
LC_KF_config->init_b_a_unc = 10000 * micro_g_to_meters_per_second_squared;
/* Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)*/
LC_KF_config->init_b_g_unc = 200 * deg_to_rad / 3600;

/* Gyro noise PSD (deg^2 per hour, converted to rad^2/s) */
LC_KF_config->gyro_noise_PSD = 0.01*0.01;
/* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
LC_KF_config->accel_noise_PSD = 0.2*0.2;
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
void Initialize_INS_GNSS_TCKF_consumer_grade_IMU(initialization_errors *initialization_errors,\
  IMU_errors *IMU_errors, GNSS_config *GNSS_config,TC_KF_config *TC_KF_config){
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
double rad_to_deg = 1/deg_to_rad;
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
IMU_errors->M_a[0] = 50000* 1E-6; IMU_errors->M_a[1] = -15000* 1E-6;
IMU_errors->M_a[2] = 10000* 1E-6;
IMU_errors->M_a[3] =-7500* 1E-6; IMU_errors->M_a[4] = -60000* 1E-6;
IMU_errors->M_a[5] = 12500* 1E-6;
IMU_errors->M_a[6] =-12500* 1E-6; IMU_errors->M_a[7] = 5000* 1E-6;
IMU_errors->M_a[8] = 20000 * 1E-6;
/* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
IMU_errors->M_g[0] = 40000* 1E-6; IMU_errors->M_a[1] = -14000* 1E-6;
IMU_errors->M_g[2] = 12500* 1E-6;
IMU_errors->M_g[3] = 0* 1E-6; IMU_errors->M_a[4] = -30000* 1E-6;
IMU_errors->M_g[5] = -7500* 1E-6;
IMU_errors->M_g[6] = 0* 1E-6; IMU_errors->M_a[7] = 0* 1E-6;
IMU_errors->M_g[8] = -17500 * 1E-6;

/* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
IMU_errors->G_g[0] = 90 * deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[1] = -110* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[2] = -60* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[3] = -50* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[4] = 190* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[5] = -160* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[6] = 30* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[7] = 110* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[8] = -130* deg_to_rad / (3600 * 9.80665);

/* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) */
IMU_errors->accel_noise_root_PSD = 1000 * micro_g_to_meters_per_second_squared;
/* Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5) FROM um7 datasheet */
IMU_errors->accel_noise_root_PSD = 400 * micro_g_to_meters_per_second_squared;
/* Gyro noise root PSD (deg per root hour, converted to rad s^-0.5) */
IMU_errors->gyro_noise_root_PSD = 1 * deg_to_rad / 60.0;
/* Gyro noise root PSD (deg per second root Hz, converted to rad s^-0.5)  From um7 datasheet */
IMU_errors->gyro_noise_root_PSD = 0.005 * deg_to_rad ;
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
for (i=0;i<3;i++) TC_KF_config->init_att_unc[i] = D2R*2.0;
/* Initial velocity uncertainty per axis (m/s) */
for (i=0;i<3;i++) TC_KF_config->init_vel_unc[i] = 0.1;
/* Initial position uncertainty per axis (m) */
for (i=0;i<3;i++) TC_KF_config->init_pos_unc[i] = 10.0;
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
TC_KF_config->gyro_noise_PSD = 0.01*0.01;
/* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
TC_KF_config->accel_noise_PSD = 0.2*0.2;
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
void Initialize_INS_GNSS_LCKF_tactical_grade_IMU(initialization_errors *initialization_errors,\
  IMU_errors *IMU_errors, GNSS_config *GNSS_config,LC_KF_config *LC_KF_config){
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
double rad_to_deg = 1/deg_to_rad;
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
IMU_errors->M_a[0] = 500* 1E-6; IMU_errors->M_a[1] = -300* 1E-6;
IMU_errors->M_a[2] = 200* 1E-6;
IMU_errors->M_a[3] =-150* 1E-6; IMU_errors->M_a[4] = -600* 1E-6;
IMU_errors->M_a[5] = 250* 1E-6;
IMU_errors->M_a[6] =-250* 1E-6; IMU_errors->M_a[7] = 100* 1E-6;
IMU_errors->M_a[8] = 450 * 1E-6;
/* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
IMU_errors->M_g[0] = 400* 1E-6; IMU_errors->M_a[1] = -300* 1E-6;
IMU_errors->M_g[2] = 250* 1E-6;
IMU_errors->M_g[3] = 0* 1E-6; IMU_errors->M_a[4] = -300* 1E-6;
IMU_errors->M_g[5] = -150* 1E-6;
IMU_errors->M_g[6] = 0* 1E-6; IMU_errors->M_a[7] = 0* 1E-6;
IMU_errors->M_g[8] = -350 * 1E-6;

/* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
IMU_errors->G_g[0] = 0.9 * deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[1] = -1.1* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[2] = -0.6* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[3] = -0.5* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[4] = 1.9* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[5] = -1.6* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[6] = 0.3* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[7] = 1.1* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[8] = -1.3* deg_to_rad / (3600 * 9.80665);

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
for (i=0;i<3;i++) LC_KF_config->init_att_unc[i] = D2R*1;
/* Initial velocity uncertainty per axis (m/s) */
for (i=0;i<3;i++) LC_KF_config->init_vel_unc[i] = 0.1;
/* Initial position uncertainty per axis (m) */
for (i=0;i<3;i++) LC_KF_config->init_pos_unc[i] = 10;
/* Initial accelerometer bias uncertainty per instrument (micro-g, converted
/* to m/s^2) */
LC_KF_config->init_b_a_unc = 1000 * micro_g_to_meters_per_second_squared;
/* Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)*/
LC_KF_config->init_b_g_unc = 10 * deg_to_rad / 3600;

/* Gyro noise PSD (deg^2 per hour, converted to rad^2/s) */
LC_KF_config->gyro_noise_PSD = pow((0.02 * deg_to_rad / 60),2);
/* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
LC_KF_config->accel_noise_PSD = pow((200 * micro_g_to_meters_per_second_squared),2);
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
void Initialize_INS_GNSS_TCKF_tactical_grade_IMU(initialization_errors *initialization_errors,\
  IMU_errors *IMU_errors, GNSS_config *GNSS_config,TC_KF_config *TC_KF_config){
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
double rad_to_deg = 1/deg_to_rad;
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
IMU_errors->M_a[0] = 500* 1E-6; IMU_errors->M_a[1] = -300* 1E-6;
IMU_errors->M_a[2] = 200* 1E-6;
IMU_errors->M_a[3] =-150* 1E-6; IMU_errors->M_a[4] = -600* 1E-6;
IMU_errors->M_a[5] = 250* 1E-6;
IMU_errors->M_a[6] =-250* 1E-6; IMU_errors->M_a[7] = 100* 1E-6;
IMU_errors->M_a[8] = 450 * 1E-6;
/* Gyro scale factor and cross coupling errors (ppm, converted to unitless;
/* body axes) */
IMU_errors->M_g[0] = 400* 1E-6; IMU_errors->M_a[1] = -300* 1E-6;
IMU_errors->M_g[2] = 250* 1E-6;
IMU_errors->M_g[3] = 0* 1E-6; IMU_errors->M_a[4] = -300* 1E-6;
IMU_errors->M_g[5] = -150* 1E-6;
IMU_errors->M_g[6] = 0* 1E-6; IMU_errors->M_a[7] = 0* 1E-6;
IMU_errors->M_g[8] = -350 * 1E-6;

/* Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes) */
IMU_errors->G_g[0] = 0.9 * deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[1] = -1.1* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[2] = -0.6* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[3] = -0.5* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[4] = 1.9* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[5] = -1.6* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[6] = 0.3* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[7] = 1.1* deg_to_rad / (3600 * 9.80665);
IMU_errors->G_g[8] = -1.3* deg_to_rad / (3600 * 9.80665);

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
for (i=0;i<3;i++) TC_KF_config->init_att_unc[i] = D2R*1;
/* Initial velocity uncertainty per axis (m/s) */
for (i=0;i<3;i++) TC_KF_config->init_vel_unc[i] = 0.1;
/* Initial position uncertainty per axis (m) */
for (i=0;i<3;i++) TC_KF_config->init_pos_unc[i] = 10;
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
TC_KF_config->gyro_noise_PSD = pow((0.02 * deg_to_rad / 60),2);
/* Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3) */
TC_KF_config->accel_noise_PSD = pow((200 * micro_g_to_meters_per_second_squared),2);
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
void Skew_symmetric (double *vec, double *W){
  W[0]= 0.0;    W[1]=-vec[2]; W[2]= vec[1];
  W[3]= vec[2]; W[4]= 0.0;    W[5]=-vec[0];
  W[6]=-vec[1]; W[7]= vec[0]; W[8]= 0.0;
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
void Gravity_ECEF(double *r_eb_e, double *g){
  double mag_r, z_scale, gamma[3];

  /* Parameters  */
  double R_0 = RE_GRS80; /* WGS84 Equatorial radius in meters */
  double omega_ie = OMGE;  /* Earth rotation rate (rad/s)  */

  /* Begins  */

  /* Calculate distance from center of the Earth  */
  mag_r = norm(r_eb_e,3);

  /* If the input position is 0,0,0, produce a dummy output */
  if (mag_r==0){
      g[0] = 0; g[1] = 0; g[2] = 0;

  /* Calculate gravitational acceleration using (2.142)  */
  }else{
      z_scale = 5 * pow( (r_eb_e[2] / mag_r),2);
      gamma[0] = - (mu / pow(mag_r,3)) * (r_eb_e[0] + (1.5 * J_2 * pow((R_0 / mag_r),2)) *\
          ((1 - z_scale) * r_eb_e[0]));
      gamma[1] = - (mu / pow(mag_r,3)) * (r_eb_e[1] + (1.5 * J_2 * pow((R_0 / mag_r),2)) *\
          ((1 - z_scale) * r_eb_e[1]));
      gamma[2] = -(mu / pow(mag_r,3)) * (r_eb_e[2] + (1.5 * J_2 * pow((R_0 / mag_r),2)) *\
          ((3 - z_scale) * r_eb_e[2]));

      /* Add centripetal acceleration using (2.133)  */
      g[0] = gamma[0] + omega_ie*omega_ie * r_eb_e[0];
      g[1] = gamma[1] + omega_ie*omega_ie * r_eb_e[1];
      g[2] = gamma[2];

  }//end % if

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
void Euler_to_CTM(double *delta_eul_nb_n, double *C){\
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
  C[3] = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
  C[4] = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
  C[5] = sin_phi * cos_theta;
  C[6] = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
  C[7] = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
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

 /* From Groves Matlab code */

 /* Grove for attitude -> Cbn 2.24

 euler[0] = atan2(Cbn[7],Cbn[8]); /* Roll (x rotation) - phi */
//imu->aea[1]=-1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6])));
 euler[1]= -asin(Cbn[6]);  /* Pitch (y rotation) - theta */

 euler[2]= atan2(Cbn[3],Cbn[0]);  /* Yaw (z rotation) - psi */
/*
 printf("Euler 0: %lf\n", euler[0]);
 printf("Euler 1: %lf\n", euler[1]);
 printf("Euler 2: %lf\n",  euler[2]);  */

 /* Calculate Euler angles using 2.23 */

  euler[0] = atan2(Cbn[5],Cbn[8]); /* Roll (x rotation) - phi */
//imu->aea[1]=-1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6])));
  euler[1]= -asin(Cbn[2]);  /* Pitch (y rotation) - theta */

  euler[2]= atan2(Cbn[1],Cbn[0]);  /* Yaw (z rotation) - psi */


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
void ECEF_to_NED(double *r_eb_e, double *v_eb_e, double *C_b_e, \
 double *L_b, double *lambda_b, double *h_b, double *v_eb_n, \
 double *C_b_n){
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

   *lambda_b = atan2(r_eb_e[1],r_eb_e[0]);

   /* From (C.29) and (C.30)  */
   k1 = sqrt(1 - e_2) * fabs(r_eb_e[2]);

   k2 = e_2 * RE_GRS80;
   beta = sqrt( (r_eb_e[0]*r_eb_e[0]) + (r_eb_e[1]*r_eb_e[1]) );
   E = (k1 - k2) / beta;
   F = (k1 + k2) / beta;

   /* From (C.31)  */
   P = 4/3 * ((E*F) + 1);

   /* From (C.32)  */
   Q = 2 * ((E*E) - (F*F));

   /* From (C.33)  */
   D = (P*P*P) + (Q*Q);

   /* From (C.34)  */
   V = pow((sqrt(D) - Q),(1/3)) - pow((sqrt(D) + Q),(1/3));

   /* From (C.35) */
   G = 0.5 * (sqrt((E*E) + V) + E);

   /* From (C.36) */
   T = sqrt((G*G) + (F - (V * G)) / ((2 * G) - E)) - G;

   /* From (C.37)  */
   *L_b = sign(r_eb_e[2]) * atan((1 - (T*T)) / (2 * T * sqrt (1 - e_2)));

   /* From (C.38) */
   *h_b = (beta - (RE_GRS80 * T)) * cos(*L_b) +\
       (r_eb_e[2] - (sign(r_eb_e[2]) * RE_GRS80 * sqrt(1 - e_2))) * sin (*L_b);

   /* using RTKLIB because above computations are not consistent */
   ecef2pos(r_eb_e, pos);
   *L_b = pos[0];
   *lambda_b = pos[1];
   *h_b=pos[2];

   /* Calculate ECEF to NED coordinate transformation matrix using (2.150) */
   cos_lat = cos(*L_b);
   sin_lat = sin(*L_b);
   cos_long = cos(*lambda_b);
   sin_long = sin(*lambda_b);
   C_e_n[0] = -sin_lat * cos_long; C_e_n[1] = -sin_lat * sin_long; C_e_n[2] =  cos_lat;
   C_e_n[3] = -sin_long;           C_e_n[4] = cos_long;            C_e_n[5] = 0;
   C_e_n[6] = -cos_lat * cos_long; C_e_n[7] = -cos_lat * sin_long; C_e_n[8] =-sin_lat;

   /* Transform velocity using (2.73) */
   matmul("NN",3,1,3,1.0,C_e_n, v_eb_e, 0.0, v_eb_n);

   /* Transform attitude using (2.15) */
   matmul("NN",3,3,3,1.0,C_e_n, C_b_e, 0.0, C_b_n);

   /* Ends */


}

void NED_to_ECEF(double *L_b, double *lambda_b, double *h_b, double *v_eb_n,\
  double *C_b_n, double* r_eb_e, double *v_eb_e, double *C_b_e)
 {
/*NED_to_ECEF - Converts curvilinear to Cartesian position, velocity
%resolving axes from NED to ECEF and attitude from NED- to ECEF-referenced
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
     C_e_n[0] = -sin_lat * cos_long; C_e_n[1] = -sin_lat * sin_long; C_e_n[2] =  cos_lat;
     C_e_n[3] = -sin_long;           C_e_n[4] = cos_long;            C_e_n[5] = 0;
     C_e_n[6] = -cos_lat * cos_long; C_e_n[7] = -cos_lat * sin_long; C_e_n[8] =-sin_lat;

    /* Transform velocity using (2.73)  */
    matmul("TN",3,1,3,1.0, C_e_n, v_eb_n, 0.0, v_eb_e);

    /* Transform attitude using (2.15) */
    matmul("TN",3,3,3,1.0,C_e_n, C_b_n, 0.0, C_b_e);

    /* Ends */

}

/* function declaration ------------------------------------------------------*/

void pv_ECEF_to_NED(double *r_eb_e, double *v_eb_e, double *L_b,\
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
   *lambda_b = atan2(r_eb_e[1],r_eb_e[0]);

   /* From (C.29) and (C.30)  */
   k1 = sqrt(1 - e_2) * fabs (r_eb_e[2]);
   k2 = e_2 * RE_GRS80;
   beta = sqrt( (r_eb_e[0]*r_eb_e[0]) + (r_eb_e[1]*r_eb_e[1]) );
   E = (k1 - k2) / beta;
   F = (k1 + k2) / beta;

   /* From (C.31)  */
   P = 4/3 * (E*F + 1);

   /* From (C.32)  */
   Q = 2 * (E*E - F*F);

   /* From (C.33)  */
   D = P*P*P + Q*Q;

   /* From (C.34)  */
   V = pow((sqrt(D) - Q),(1/3)) - pow((sqrt(D) + Q),(1/3));

   /* From (C.35) */
   G = 0.5 * (sqrt(E*E + V) + E);

   /* From (C.36) */
   T = sqrt(G*G + (F - V * G) / (2 * G - E)) - G;

   /* From (C.37)  */
   *L_b = sign(r_eb_e[2]) * atan((1 - T*T) / (2 * T * sqrt (1 - e_2)));

   /* From (C.38) */
   *h_b = (beta - RE_GRS80 * T) * cos(*L_b) +\
       (r_eb_e[2] - sign(r_eb_e[2]) * RE_GRS80 * sqrt(1 - e_2)) * sin (*L_b);

    /* Using RTKLIB because computations are not consistent  */
    ecef2pos(r_eb_e, pos);
    *L_b = pos[0];
    *lambda_b = pos[1];
    *h_b=pos[2];

   /* Calculate ECEF to NED coordinate transformation matrix using (2.150) */
   cos_lat = cos(*L_b);
   sin_lat = sin(*L_b);
   cos_long = cos(*lambda_b);
   sin_long = sin(*lambda_b);
   C_e_n[0] = -sin_lat * cos_long; C_e_n[1] = -sin_lat * sin_long; C_e_n[2] =  cos_lat;
   C_e_n[3] = -sin_long;           C_e_n[4] = cos_long;            C_e_n[5] = 0;
   C_e_n[6] = -cos_lat * cos_long; C_e_n[7] = -cos_lat * sin_long; C_e_n[8] =-sin_lat;

   /* Transform velocity using (2.73) */
   matmul("NN",3,1,3,1.0,C_e_n, v_eb_e, 0.0, v_eb_n);

   /* Ends */

}

void pv_NED_to_ECEF(double *L_b, double *lambda_b, double *h_b, double *v_eb_n,\
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
     C_e_n[0] = -sin_lat * cos_long; C_e_n[1] = -sin_lat * sin_long; C_e_n[2] =  cos_lat;
     C_e_n[3] = -sin_long;           C_e_n[4] = cos_long;            C_e_n[5] = 0;
     C_e_n[6] = -cos_lat * cos_long; C_e_n[7] = -cos_lat * sin_long; C_e_n[8] =-sin_lat;

    /* Transform velocity using (2.73)  */
    matmul("TN",3,1,3,1.0, C_e_n, v_eb_n, 0.0, v_eb_e);

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
void Calculate_errors_NED(double *est_L_b, double *est_lambda_b, \
  double *est_h_b, double *est_v_eb_n, double *est_C_b_n, double *true_L_b,\
  double *true_lambda_b, double *true_h_b, double *true_v_eb_n, \
  double *true_C_b_n, double *delta_r_eb_n, double *delta_v_eb_n, \
  double *delta_eul_nb_n){

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
void Initialize_NED_attitude(double *C_b_n,\
  initialization_errors *initialization_errors, double *est_C_b_n){
    double delta_C_b_n[9];
    int i;

  /* Attitude initialization, using (5.109) and (5.111) */
  for(i=0;i<3;i++) initialization_errors->delta_eul_nb_n[i]=\
  -initialization_errors->delta_eul_nb_n[i];

  Euler_to_CTM(initialization_errors->delta_eul_nb_n, delta_C_b_n);

  /* est_C_b_n = delta_C_b_n * C_b_n; */
  matmul("NN",3,3,3,1.0,delta_C_b_n,C_b_n,0.0,est_C_b_n);

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
void Initialize_TC_P_matrix(TC_KF_config *TC_KF_config, double *P_matrix){
  int i,j, npar=17;
  double att_var[3], vel_var[3], pos_var[3], ba_var, bg_var, dt_off_var, dt_drift_var;


  /* Initialize error covariance matrix */
  // IT ASSUMES THE SAME VALUE FOR EACH AXES, HOWEVER IT MAY DIFFER!!!!!!!!!
  for (i=0;i<3;i++) {
    att_var[i] = TC_KF_config->init_att_unc[i] * TC_KF_config->init_att_unc[i];
    vel_var[i] = TC_KF_config->init_vel_unc[i] * TC_KF_config->init_vel_unc[i];
    pos_var[i] = TC_KF_config->init_pos_unc[i] * TC_KF_config->init_pos_unc[i];
  }

  ba_var = TC_KF_config->init_b_a_unc * TC_KF_config->init_b_a_unc;
  bg_var = TC_KF_config->init_b_g_unc * TC_KF_config->init_b_g_unc;
  dt_off_var = TC_KF_config->init_clock_offset_unc * TC_KF_config->init_clock_offset_unc;
  dt_drift_var = TC_KF_config->init_clock_drift_unc * TC_KF_config->init_clock_drift_unc;

  /* The 17-states {de,dv,dr,dba,dbg, dt, dtdot} system model order */

  for (i=0;i<3;i++) {
    printf("att: %lf, vel: %lf, pos: %lf\n", att_var[i], vel_var[i], pos_var[i]);
  }

  /* Attitude  */
   for (i = 0; i < 3; i++) {
     for (j = 0; j < 3; j++) {
       (i==j?P_matrix[i*npar+j]=att_var[i]:0.0);
     }
   }

  /* Velocity */
   for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      (i==j?P_matrix[i*npar+j]=vel_var[j-3]:0.0);
    }
  }

  /* Position */
  for (i = 6; i < 9; i++) {
    for (j = 6; j < 9; j++) {
      (i==j?P_matrix[i*npar+j]=pos_var[j-6]:0.0);
    }
  }

  /* Acc. bias  */
   for (i = 9; i < 12; i++) {
     for (j = 9; j < 12; j++) {
       (i==j?P_matrix[i*npar+j]=ba_var:0.0);
     }
   }

   /* Gyro. bias  */
    for (i = 12; i < 15; i++) {
      for (j = 12; j < 15; j++) {
        (i==j?P_matrix[i*npar+j]=bg_var:0.0);
      }
    }

   /* Clock offset and drift variances */
   P_matrix[15*npar+15]=dt_off_var;
   P_matrix[16*npar+16]=dt_drift_var;

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
void Initialize_LC_P_matrix(LC_KF_config *LC_KF_config, double *P_matrix){
  int i,j, npar=15;
  double att_var[3], vel_var[3], pos_var[3], ba_var, bg_var;


  /* Initialize error covariance matrix */
  // IT ASSUMES THE SAME VALUE FOR EACH AXES, HOWEVER IT MAY DIFFER!!!!!!!!!
  for (i=0;i<3;i++) {
    att_var[i] = LC_KF_config->init_att_unc[i] * LC_KF_config->init_att_unc[i];
    vel_var[i] = LC_KF_config->init_vel_unc[i] * LC_KF_config->init_vel_unc[i];
    pos_var[i] = LC_KF_config->init_pos_unc[i] * LC_KF_config->init_pos_unc[i];
  }

  ba_var = LC_KF_config->init_b_a_unc * LC_KF_config->init_b_a_unc;
  bg_var = LC_KF_config->init_b_g_unc * LC_KF_config->init_b_g_unc;


  /* The 15-states {de,dv,dr,dba,dbg} system model order */

  /* Attitude  */
   for (i = 0; i < 3; i++) {
     for (j = 0; j < 3; j++) {
       (i==j?P_matrix[i*npar+j]=att_var[j]:0.0);
     }
   }

  /* Velocity */
   for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      (i==j?P_matrix[i*npar+j]=vel_var[j-3]:0.0);
    }
  }

  /* Position */
  for (i = 6; i < 9; i++) {
    for (j = 6; j < 9; j++) {
      (i==j?P_matrix[i*npar+j]=pos_var[j-6]:0.0);
    }
  }

  /* Acc. bias  */
   for (i = 9; i < 12; i++) {
     for (j = 9; j < 12; j++) {
       (i==j?P_matrix[i*npar+j]=ba_var:0.0);
     }
   }

   /* Gyro. bias  */
    for (i = 12; i < 15; i++) {
      for (j = 12; j < 15; j++) {
        (i==j?P_matrix[i*npar+j]=bg_var:0.0);
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
quaternion_mult(double *p, double *q, double *pq){
  pq[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];
  pq[1] = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2];
  pq[2] = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1];
  pq[3] = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0];
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
void attitude_update(double *alpha_ib_b, double *omega_ie, float t, double* old_C_b_e, double *new_q_b_e)
{
 double old_q_b_e[4], q_less_plus[4], q_omega[4], q_aux[4], q_aux1[4];
 double mag_alpha,ac=0.0, as=0.0;
 int i;

 /* Form quaternions */

 /* Old Cbe quaternion */
 DCM_to_quaternion(old_C_b_e, old_q_b_e);

 /* Reverse quaternion to become old_q_e_b */
 old_q_b_e[0]=-old_q_b_e[0];

 /* Coeficients of q_less_plus quaternion */
 mag_alpha=norm(alpha_ib_b,3);
 ac=1-0.5*(mag_alpha*mag_alpha/4)+((1/24)*(pow(mag_alpha/2,4)));
 as=0.5-((1/12)*(pow(mag_alpha/2,2)));

 /* q_less_plus quaternion */
 q_less_plus[0]=ac;
 for (i=1;i<4;i++) q_less_plus[i]=as*alpha_ib_b[i-1];

 /*q_omega quaternion */
 q_omega[0]=0.0;
 for (i=1;i<4;i++) q_omega[i]= 0.5*omega_ie[i-1]*t;

 quaternion_mult(q_omega, old_q_b_e, q_aux);

 quaternion_mult(old_q_b_e, q_less_plus, q_aux1);

 /* The updated quternion */
 for (i=0;i<4;i++) new_q_b_e[i] = q_aux1[i] - q_aux[i];

 /* Reverse quaternion */
 new_q_b_e[0]=-new_q_b_e[0];

 }

/* Convert rotation matrix to quaternion ----------------- ---------------------
* description: Cbe or Cbeta_alpha to quaternion transformation
* args   : um7pack_t* imu       IO    imu structure
*
Reference: Groves (2013, pag.41)
*-----------------------------------------------------------------------------*/
void DCM_to_quaternion(double* C, double* q)
{
 /* Transformation between DCM Cbn and the quaternion is accomplished by: */
 q[0]=0.5*sqrt(1+C[0]+C[4]+C[8]);

 if (q[0]<0.001) {
   q[1]=0.5*sqrt(1+C[0]+C[4]+C[8]);
   q[0]=0.25*(C[5]-C[7])/(q[1]);
   q[2]=0.25*(C[3]-C[1])/(q[1]);
   q[3]=0.25*(C[6]-C[2])/(q[1]);
 }else{
   q[1]=0.25*(C[5]-C[7])/(q[0]);
   q[2]=0.25*(C[6]-C[2])/(q[0]);
   q[3]=0.25*(C[1]-C[3])/(q[0]);
 }

}

/* Quaternion to DCM ------------------------------------- ---------------------
* description: transform a quaternion to DCM Cbn amtrix
* args   : double* C[9]       IO    DCM Cbe matrix
* 	       double* q[4]	      I   quaternion
Reference: Groves (2013, pag.41)
*-----------------------------------------------------------------------------*/
void Quaternion_to_DCM(double* q, double* C)
{
  /* The transformation between the quaternion and the DCM Cbn is accomplished by */
  C[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]; C[1]=2*(q[1]*q[2]+q[3]*q[0]); C[2]=2*(q[1]*q[3]-q[2]*q[0]);
  C[3]=2*(q[1]*q[2]-q[3]*q[0]); C[4]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3]; C[5]=2*(q[2]*q[3]+q[1]*q[0]);
  C[6]=2*(q[1]*q[3]+q[2]*q[0]); C[7]=2*(q[2]*q[3]-q[1]*q[0]); C[8]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

/* Euler to quaternion ----------------------------------- ---------------------
* description: transform an euler vector to quaternion
* args   : double *eul_nb       I    Euler angles vector {3x1}
* 	       double* q_nb[4]	      IO   quaternion {4x1}
Reference: Groves (2013, pag.42)
*-----------------------------------------------------------------------------*/
Euler_to_quaternion(double *eul_nb, double *q_nb){
  double phi_half=eul_nb[0]/2, theta_half=eul_nb[1]/2, \
  yaw_half=eul_nb[2]/2;

  q_nb[0]=cos(phi_half)*cos(theta_half)*cos(yaw_half)+sin(phi_half)*sin(theta_half)*sin(yaw_half);
  q_nb[1]=sin(phi_half)*cos(theta_half)*cos(yaw_half)-cos(phi_half)*sin(theta_half)*sin(yaw_half);
  q_nb[2]=cos(phi_half)*sin(theta_half)*cos(yaw_half)+sin(phi_half)*cos(theta_half)*sin(yaw_half);
  q_nb[3]=cos(phi_half)*cos(theta_half)*sin(yaw_half)-sin(phi_half)*sin(theta_half)*cos(yaw_half);
}

/* Quaternion to Euler ---------------------------------- ---------------------
* description: transform an quaternion to euler vector
* args   : double *q_nb      I    quaternion {4x1}
* 	       double  *eul_nb      IO   Euler angles vector {3x1}
Reference: Groves (2013, pag.42)
*-----------------------------------------------------------------------------*/
Quaternion_to_euler(double *q, double *eul){

  eul[0]=atan2(2*(q[0]*q[1]+q[2]*q[3]),(1-2*q[1]-2*q[2]));
  eul[1]=asin(2*(q[0]*q[2]-q[1]*q[3]));
  eul[2]=atan2(2*(q[0]*q[3]+q[1]*q[2]),(1-2*q[2]-2*q[3]));

  printf("Quat.euler: %lf\n", 2*(q[0]*q[2]-q[1]*q[3]) );

}


/* Attitude update with quaternion --------------------- ---------------------
* description: quaternion approach for updating the attitude matrix Cbn
* args   : um7pack_t* imu       IO    imu structure
*
Reference: Shin (2001, pag.18)
*-----------------------------------------------------------------------------*/
void Cbn_from_rotations(um7pack_t* imu, double* C)
{
 double q[4], u, qdoti[16], qdot[4];

 /* Building quaternion from a rotation angle vector u{ux,uy,uz} */

 /* QUEST: Is it from the attitude Euler angles or  the rotation rate from gyroscopes measurements */
 u = sqrt(imu->aea[0]*imu->aea[0]+imu->aea[1]*imu->aea[1]+imu->aea[2]*imu->aea[2]);

 q[0]=(imu->g[0]/u)*sin(u/2);
 q[1]=(imu->g[1]/u)*sin(u/2);
 q[2]=(imu->g[2]/u)*sin(u/2);
 q[3]=cos(u/2);

/* Normality condition
 if( q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] > 1.000001 || q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] < 0.999999 )
 q[]=q[]/sqrt(qTq); //for each component
*/

 qdoti[0]=qdoti[5]=qdoti[10]=qdoti[15]=0;
 qdoti[1]=qdoti[11]=imu->g[2]; /*wz*/
 qdoti[4]=qdoti[14]=-imu->g[2]; /*-wz*/
 qdoti[7]=qdoti[8]=imu->g[1]; /*wy*/
 qdoti[13]=qdoti[2]=-imu->g[1]; /*-wy*/
 qdoti[3]=qdoti[6]=imu->g[0]; /*wx*/
 qdoti[9]=qdoti[12]=-imu->g[0]; /*-wx*/

 matmul("NN", 4, 1, 4, 0.5, qdoti, q, 0.0, qdot);

 /* The transformation between the quaternion and the DCM Cbn is accomplished by */
 C[0]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]-q[3]*q[3]; C[1]=2*(q[0]*q[1]-q[2]*q[3]); C[2]=2*(q[0]*q[2]-q[1]*q[3]);
 C[3]=2*(q[0]*q[1]-q[2]*q[3]); C[4]=q[1]*q[1]-q[0]*q[0]-q[2]*q[2]+q[3]*q[3]; C[5]=2*(q[1]*q[2]-q[0]*q[3]);
 C[6]=2*(q[0]*q[2]-q[1]*q[3]); C[7]=2*(q[1]*q[2]-q[0]*q[3]); C[8]=q[2]*q[2]-q[0]*q[0]-q[1]*q[1]+q[3]*q[3];

}

/* Cbn update from attitude error estimate --------------- ---------------------
* description: Corrects previous Cbe by the estimated attitude error to a new Cbe
* args   : double *est_delta_Psi     I    vector with attitude errors {3x1}
           double *est_C_b_e_old     I    old Cbe transformation matrix {3x3}
           double *est_C_b_e_new     O    new Cbe transformation matrix {3x3}
*
Reference: Grove (2013) Appendix E, section E.6.3 page E-14
*-----------------------------------------------------------------------------*/
void Quaternion_attitude_errror_correction(double *est_delta_Psi, \
  double *est_C_b_e_old, double *est_C_b_e_new){
  double q_new[4], q_old[4], M_delta_Psi[16];

  /* Obtain Cbe old quaternion */
  DCM_to_quaternion(est_C_b_e_old, q_old);

  printf("OLD_QUAT NORM %lf\n", norm(q_old,4) );


  M_delta_Psi[0]=M_delta_Psi[5]=M_delta_Psi[10]=M_delta_Psi[15]=1.0;
  M_delta_Psi[1]=M_delta_Psi[11]=-est_delta_Psi[0]; /*-wx*/
  M_delta_Psi[4]=M_delta_Psi[14]= est_delta_Psi[0]; /*wx*/
  M_delta_Psi[2]=M_delta_Psi[13]= -est_delta_Psi[1]; /*wy*/
  M_delta_Psi[7]=M_delta_Psi[8]= est_delta_Psi[1]; /*-wy*/
  M_delta_Psi[3]=M_delta_Psi[6]= -est_delta_Psi[2]; /*wz*/
  M_delta_Psi[9]=M_delta_Psi[12]= est_delta_Psi[2]; /*-wz*/

  matmul("NN", 4, 1, 4, -0.5, M_delta_Psi, q_old, 0.0, q_new);

  printf("NEW_QUAT NORM %lf\n", norm(q_new,4) );

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
 void TC_KF_Epoch(GNSS_measurements *GNSS_measurements, int no_meas,\
   float tor_s, double *est_C_b_e_old, double *est_v_eb_e_old, double *est_r_eb_e_old,\
   double *est_IMU_bias_old, double *est_clock_old, double *P_matrix_old, double *meas_f_ib_b,\
   double est_L_b_old, TC_KF_config *TC_KF_config, double *est_C_b_e_new,\
   double *est_v_eb_e_new, double *est_r_eb_e_new, double *est_IMU_bias_new,\
   double *est_clock_new, double *P_matrix_new){

    /* Constants (sone of these could be changed to inputs at a later date) */
    double c = 299792458; /* Speed of light in m/s */
    double omega_ie = OMGE;  /* Earth rotation rate in rad/s */
    double R_0 = RE_WGS84; /*WGS84 Equatorial radius in meters */
    double e = sqrt(e_2); /*WGS84 eccentricity                        */

    /* Begins */
    double omega_ie_vec[3], Omega_ie[9];
    double *Phi_matrix, *Q_prime_matrix, *x_est_propagated, *x_est_new;
    double *Q_, *Q_aux;
    double *P_matrix_propagated, *P_matrix, *P_aux, *u_as_e_T, *pred_meas;
    double *H_matrix, *R_matrix, *ones, *delta_r, *delta_z, approx_range, range;
    double range_rate, C_e_I[9], *I_meas, *K_matrix, *K_matrix_inv;
    double meas_f_ib_e[3], Skew_meas_f_ib_e[9], est_r_eb_e[3], delta_r_aux[3];
    double geocentric_radius, g[3], g_est_r_eb_e[9], *I, *I_par;
    double Omega_x_est_r_eb_old[3], Omega_x_GNSS_meas_r[3], v_plus_Omega_x_est_r[3];
    double v_plus_Omega_x_GNSS_meas_r[3], est_v_plus_Omega_x_est_r[3];
    double C_x_v[3],C_less_v_plus_Omega[3];
    double *F1, *Q1, *K_x_delta_z, *K_x_H, *I_less_k_x_H;
    double Skew_x_est_new[9], I_less_Skew[9];
    double rho, u[3], rate, vs[3], a[3], e_[3], cosel, azel[2*no_meas], pos[3];
    double E[9];
    int i,j, info;

    /* Initialize matrices and vectors */
    Phi_matrix=eye(17); Q_prime_matrix=zeros(17,17);
    x_est_propagated=zeros(17,1); x_est_new=zeros(17,1);
    P_matrix_propagated=mat(17,17); P_aux=mat(17,17); Q_aux=mat(17,17);
    P_matrix=mat(17,17); u_as_e_T=zeros(no_meas,3); pred_meas=zeros(no_meas,2);
    H_matrix=zeros((2 * no_meas),17); ones=mat(no_meas,1); I_meas=eye(no_meas);
    R_matrix=zeros((2 * no_meas),(2 * no_meas)); I_par=eye(17);
    K_matrix_inv=zeros((2 * no_meas),(2 * no_meas));
    K_matrix=mat(17,(2 * no_meas));
    delta_r=mat(no_meas,3); delta_z=mat(2 * no_meas,1); Q_=mat(17,17);
    F1=mat(17,2*no_meas); Q1=mat(2*no_meas,2*no_meas); K_x_delta_z=mat(17,1);
    K_x_H=mat(17,17); I_less_k_x_H=mat(17,17); I = eye(3);


    for (i=0;i<no_meas;i++) ones[i]=1.0;

    /* Skew symmetric matrix of Earth rate */
    omega_ie_vec[0]=0; omega_ie_vec[1]=0; omega_ie_vec[2]=OMGE;
    Skew_symmetric(omega_ie_vec, Omega_ie);

    /* SYSTEM PROPAGATION PHASE */

    /* 1. Determine transition matrix using (14.50) (first-order approx) */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        Phi_matrix[i*17+j] = Phi_matrix[i*17+j] - Omega_ie[i*3+j] * tor_s;
      }
    }

    for (i = 0; i < 3; i++) {
      for (j = 12; j < 15; j++) {
        Phi_matrix[i*17+j] = est_C_b_e_old[i*3+(j-12)] * tor_s;
      }
    }

    matmul("NN",3,1,3,1.0,est_C_b_e_old,meas_f_ib_b,0.0, meas_f_ib_e);
    Skew_symmetric(meas_f_ib_b,Skew_meas_f_ib_e);

    for (i = 3; i < 6; i++) {
      for (j = 0; j < 3; j++) {
        Phi_matrix[i*17+j] = -tor_s * Skew_meas_f_ib_e[(i-3)*3+j];
      }
    }

    for (i = 3; i < 6; i++) {
      for (j = 3; j < 6; j++) {
        Phi_matrix[i*17+j] = Phi_matrix[i*17+j] - 2 * Omega_ie[(i-3)*3+(j-3)] * tor_s;
      }
    }

    geocentric_radius = R_0 / sqrt(1 - pow((e * sin(est_L_b_old)),2)) *\
        sqrt(pow(cos(est_L_b_old),2) + pow((1 - e*e),2) * pow(sin(est_L_b_old),2)); /* from (2.137)*/
    Gravity_ECEF(est_r_eb_e_old, g);
    matmul("NT",3,3,1,1.0,g,est_r_eb_e,0.0, g_est_r_eb_e);
    for (i = 3; i < 6; i++) {
      for (j = 6; j < 9; j++) {
        Phi_matrix[i*17+j] = -tor_s * 2 /\
         geocentric_radius * g_est_r_eb_e[(i-3)*3+(j-6)] / (norm(est_r_eb_e_old,3));
      }
    }

    for (i = 3; i < 6; i++) {
      for (j = 9; j < 12; j++) {
        Phi_matrix[i*17+j] = est_C_b_e_old[(i-3)*3+(j-9)] * tor_s;
      }
    }

    for (i = 6; i < 9; i++) {
      for (j = 3; j < 6; j++) {
        Phi_matrix[i*17+j] = I[(i-6)*3+(j-3)] * tor_s;
      }
    }

    Phi_matrix[15*17+16] = tor_s;
/*
    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        printf("Phi_matrix:[%d] %lf\t", i*17+j,Phi_matrix[i*17+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
    /* 2. Determine approximate system noise covariance matrix using (14.82) */
    //Q_prime_matrix(1:3,1:3) = eye(3) * TC_KF_config.gyro_noise_PSD * tor_s;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        Q_prime_matrix[i*17+j] = (i==j?I[(i)*3+(j)] *\
         TC_KF_config->gyro_noise_PSD * tor_s:0.0);
      }
    }
    //Q_prime_matrix(4:6,4:6) = eye(3) * TC_KF_config.accel_noise_PSD * tor_s;
    for (i = 3; i < 6; i++) {
      for (j = 3; j < 6; j++) {
        Q_prime_matrix[i*17+j] = (i==j?I[(i-3)*3+(j-3)] *\
         TC_KF_config->accel_noise_PSD * tor_s:0.0);
      }
    }
   //Q_prime_matrix(10:12,10:12) = eye(3) * TC_KF_config.accel_bias_PSD * tor_s;
    for (i = 9; i < 12; i++) {
      for (j = 9; j < 12; j++) {
        Q_prime_matrix[i*17+j] = (i==j?I[(i-9)*3+(j-9)] *\
         TC_KF_config->accel_bias_PSD * tor_s:0.0);
      }
    }
   //Q_prime_matrix(13:15,13:15) = eye(3) * TC_KF_config.gyro_bias_PSD * tor_s;
    for (i = 12; i < 15; i++) {
      for (j = 12; j < 15; j++) {
        Q_prime_matrix[i*17+j] = (i==j?I[(i-12)*3+(j-12)] *\
         TC_KF_config->gyro_bias_PSD * tor_s:0.0);
      }
    }
    //Q_prime_matrix(16,16) = TC_KF_config.clock_phase_PSD * tor_s;
    Q_prime_matrix[15*17+15] = TC_KF_config->clock_phase_PSD * tor_s;

    //Q_prime_matrix(17,17) = TC_KF_config->clock_freq_PSD * tor_s;
    Q_prime_matrix[16*17+16] = TC_KF_config->clock_freq_PSD * tor_s;
/*
    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        printf("Q_prime_matrix:[%d] %lf\t", i*17+j,Q_prime_matrix[i*17+j]);
      }
      printf("\n");
    }
    printf("\n");
*/

    /* 3. Propagate state estimates using (3.14) noting that only the clock
    /* states are non-zero due to closed-loop correction.
    printf("est_clock: %lf, %lf\n", est_clock_old[0], est_clock_old[1]);
    printf("TOR_S: %lf\n", tor_s);*/

    //printf("REC. CLK OFF: %lf m AND DRIFT: %lf m/s \n",est_clock_old[0],est_clock_old[1]);

    /* Using values estimated in RTKLIB per epoch */
    x_est_propagated[15] = est_clock_old[0];// + est_clock_old[1] * tor_s;
    x_est_propagated[16] = est_clock_old[1];
    //x_est_propagated[15] = 0.0; //+ est_clock_old[1] * tor_s;
    //x_est_propagated[16] = 0.0;


    /* 4. Propagate state estimation error covariance matrix using (3.46) */
    //P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *\
    //    Phi_matrix' + 0.5 * Q_prime_matrix;
    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        P_aux[i*17+j]= P_matrix_old[i*17+j]+0.5*Q_prime_matrix[i*17+j];
      }
    }

    matmul("NT",17,17,17,1.0,P_aux,Phi_matrix,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
    matmul("NN",17,17,17,1.0,Phi_matrix,Q_aux,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        P_matrix_propagated[i*17+j]= Q_[i*17+j]+0.5*Q_prime_matrix[i*17+j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
      }
    }
/*

    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        printf("P_matrix_old[%d] %lf\t", i*17+j,P_matrix_old[i*17+j]);
      }
      printf("\n");
    }
    printf("\n");

    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        printf("Q_[%d] %lf\t", i*17+j,Q_[i*17+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
    /* MEASUREMENT UPDATE PHASE */

     ecef2pos(est_r_eb_e_old,pos);
     xyz2enu(pos,E);

    /* Loop measurements */
    for (j=0;j<no_meas;j++){  //no_meas represent the number of visible satellites
      azel[j*2]=azel[1+j*2]=0.0;

        /* Predict approx range */
        for(i=0;i<3;i++) delta_r[j*3+i] = GNSS_measurements[j].Sat_r_eb_e[i] - est_r_eb_e_old[i];
        //printf("\ndelta_r: %lf, %lf, %lf\n",delta_r[j*3], delta_r[j*3+1],delta_r[j*3+2]);

        rho=geodist(GNSS_measurements[j].Sat_r_eb_e, est_r_eb_e_old, u);
        satazel(pos,u,azel+j*2);

        approx_range = norm(delta_r+(j*3),3);
        //printf("approx_range: %lf\n",approx_range);

        /* Calculate frame rotation during signal transit time using (8.36) */
        C_e_I[0] = 1; C_e_I[1] = omega_ie * approx_range / c; C_e_I[2] = 0;
        C_e_I[3] = -omega_ie * approx_range / c; C_e_I[4] =  1; C_e_I[5] = 0;
        C_e_I[6] = 0; C_e_I[7] = 0; C_e_I[8] = 1;

        /* Predict pseudo-range using (9.165) */
        matmul("NT",3,1,3,1.0,C_e_I,GNSS_measurements[j].Sat_r_eb_e,0.0,delta_r_aux);
        for(i=0;i<3;i++) delta_r[j*3+i] = delta_r_aux[i] - est_r_eb_e_old[i];
        range = norm(delta_r+(j*3),3);


        printf("RHO RTKLIB:%lf AND RHO GROVES: %lf\n", rho, range);
        range=rho;

        pred_meas[j*2+0] = range + x_est_propagated[15];


        /* Predict line of sight */
        for(i=0;i<3;i++) u_as_e_T[j*3+i] = delta_r[j*3+i] / range;

        for(i=0;i<3;i++) u_as_e_T[j*3+i] = u[i];


        /* Predict pseudo-range rate using (9.165) */
        /*range_rate = u_as_e_T[j*3+i] * (C_e_I[] * (GNSS_measurements[j,6:8]' +...
            Omega_ie * GNSS_measurements[j,3:5]') - (est_v_eb_e_old +...
            Omega_ie * est_r_eb_e_old));*/

        /* RTKLIB range-rate */
        /* line-of-sight vector in ecef */
        cosel=cos(azel[1+j*2]);
        a[0]=sin(azel[j*2])*cosel;
        a[1]=cos(azel[j*2])*cosel;
        a[2]=sin(azel[1+j*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e_);

        /* satellite velocity relative to receiver in ecef */
        for (i=0;i<3;i++) vs[i]=GNSS_measurements[j].Sat_v_eb_e[i]-est_v_eb_e_old[i];


        /* range rate with earth rotation correction */
        rate=dot(vs,e_,3)+OMGE/CLIGHT*(GNSS_measurements[j].Sat_v_eb_e[1]*est_r_eb_e_old[0]+\
          GNSS_measurements[j].Sat_r_eb_e[1]*est_v_eb_e_old[0]-\
          GNSS_measurements[j].Sat_v_eb_e[0]*est_r_eb_e_old[1]-GNSS_measurements[j].Sat_r_eb_e[0]*est_v_eb_e_old[1]);

            /* doppler residual
            v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);*/

        matmul("NN",3,1,3,1.0,Omega_ie,est_r_eb_e_old,0.0, Omega_x_est_r_eb_old);
        matmul("NN",3,1,3,1.0,Omega_ie,GNSS_measurements[j].Sat_r_eb_e,0.0,Omega_x_GNSS_meas_r);

        for(i=0;i<3;i++) v_plus_Omega_x_GNSS_meas_r[i] = \
        GNSS_measurements[j].Sat_v_eb_e[i]+Omega_x_GNSS_meas_r[i];

        for(i=0;i<3;i++) est_v_plus_Omega_x_est_r[i] = \
        est_v_eb_e_old[i]+Omega_x_est_r_eb_old[i];

        matmul("NN",3,1,3,1.0,C_e_I,v_plus_Omega_x_GNSS_meas_r,0.0,C_x_v);

        for(i=0;i<3;i++) C_less_v_plus_Omega[i] = \
        C_x_v[i]-est_v_plus_Omega_x_est_r[i];

        matmul("NN",1,1,3,1.0,u_as_e_T+(j*3),C_less_v_plus_Omega,0.0,&range_rate);

        range_rate=rate;

        pred_meas[j*2+1] = range_rate + x_est_propagated[16];
        printf("range_rate[%d]:%lf AND RTKLIb: %lf\n",j*2+1,range_rate, rate );

    } //end for j


    /* 5. Set-up measurement matrix using (14.126) */

    //H_matrix(1:no_meas,7:9) = u_as_e_T(1:no_meas,1:3); - Position
    for (i = 0; i < no_meas; i++) {
      for (j = 6; j < 9; j++) {
        H_matrix[i*17+j] = u_as_e_T[i*3+(j-6)];
      }
    }

    //H_matrix(1:no_meas,16) = ones(no_meas,1);  - Clock offset
    for (i = 0; i < no_meas; i++) H_matrix[i*17+15] = 1.0;

    //H_matrix((no_meas + 1):(2 * no_meas),4:6) = u_as_e_T(1:no_meas,1:3); - Velocity
    for (i = no_meas; i < 2*no_meas; i++) {
      for (j = 3; j < 6; j++) {
        H_matrix[i*17+j] = u_as_e_T[(i-no_meas)*3+(j-3)];
      }
    }

    //H_matrix((no_meas + 1):(2 * no_meas),17) = ones(no_meas,1);  -  Clock drift
    for (i = no_meas; i < 2*no_meas; i++) {
      H_matrix[i*17+16] = 1.0;
    }
/*
    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < 17; j++) {
        printf("H[%d]:%lf\t",(int)(i*17+j),H_matrix[i*17+j]);
      }
      printf("\n");
    }
    printf("\n"); */

    /* 6. Set-up measurement noise covariance matrix assuming all measurements
    /* are independent and have equal variance for a given measurement type. */

    //R_matrix(1:no_meas,1:no_meas) = eye(no_meas) *TC_KF_config.pseudo_range_SD^2;
    for (i = 0; i < no_meas; i++) {
      for (j = 0; j < no_meas; j++) {
        R_matrix[i*(2*no_meas)+j] = I_meas[i*(no_meas)+j]*\
        pow(TC_KF_config->pseudo_range_SD,2);
      }
    }

    //R_matrix(1:no_meas,(no_meas + 1):(2 * no_meas)) = zeros(no_meas);
    for (i = 0; i < no_meas; i++) {
      for (j = no_meas; j < 2*no_meas; j++) {
        R_matrix[i*(2*no_meas)+j] = 0.0;
      }
    }

    //R_matrix((no_meas + 1):(2 * no_meas),1:no_meas) =  zeros(no_meas);
    for (i = no_meas; i < 2*no_meas; i++) {
      for (j = 0; j < no_meas; j++) {
        R_matrix[i*(2*no_meas)+j] = 0.0;
      }
    }

    //R_matrix((no_meas + 1):(2 * no_meas),(no_meas + 1):(2 * no_meas)) =...
    //  eye(no_meas) * TC_KF_config.range_rate_SD^2;
    for (i = no_meas; i < 2*no_meas; i++) {
      for (j = no_meas; j < 2*no_meas; j++) {
        R_matrix[i*(2*no_meas)+j] = I_meas[(i-no_meas)*(no_meas)+(j-no_meas)] *\
         pow(TC_KF_config->range_rate_SD,2);
      }
    }

    /* 7. Calculate Kalman gain using (3.21) */
   //K_matrix = P_matrix_propagated * H_matrix' * inv(H_matrix *...
    //    P_matrix_propagated * H_matrix' + R_matrix);
    matmul("NT",17,2*no_meas,17,1.0,P_matrix_propagated,H_matrix,0.0,F1);
    matmul("NN",2*no_meas,2*no_meas,17,1.0,H_matrix,F1,0.0,Q1);
    for (i = 0; i < 2*no_meas; i++) {
      for (j = 0; j < 2*no_meas; j++) {
        K_matrix_inv[i*(2*no_meas)+j] = Q1[i*(2*no_meas)+j]+ R_matrix[i*(2*no_meas)+j];
      }
    }
    /*
    for(i=0;i<17;i++) {
      for (j=0; j < 17; j++) {
        printf("P_matrix_propagated[%d]:%lf\t",(int)(i*(17)+j),P_matrix_propagated[i*(17)+j]);
      }
      printf("\n");
    }
    printf("\n");

    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < 17; j++) {
        printf("F1[%d]:%lf\t",(int)(i*(17)+j),F1[i*(17)+j]);
      }
      printf("\n");
    }
    printf("\n");

    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < 17; j++) {
        printf("H_matrix[%d]:%lf\t",(int)(i*(17)+j),H_matrix[i*(17)+j]);
      }
      printf("\n");
    }
    printf("\n");

    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < no_meas*2; j++) {
        printf("R_matrix[%d]:%lf\t",(int)(i*(no_meas*2)+j),R_matrix[i*(no_meas*2)+j]);
      }
      printf("\n");
    }
    printf("\n");

    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < no_meas*2; j++) {
        printf("Q1[%d]:%lf\t",(int)(i*(no_meas*2)+j),Q1[i*(no_meas*2)+j]);
      }
      printf("\n");
    }
    printf("\n");

    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < no_meas*2; j++) {
        printf("K_matrix_inv[%d]:%lf\t",(int)(i*(no_meas*2)+j),K_matrix_inv[i*(no_meas*2)+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
    if (!(info=matinv(K_matrix_inv,2*no_meas))) {
      printf("Invertion status, if 0 > info:error! info:%d\n", info);
        matmul("NN",17,2*no_meas,2*no_meas,1.0,F1,K_matrix_inv,0.0,K_matrix);
        printf("\n\n\n\n********************************** INVERTING MATRIX  **************************************\n\n\n\n\n");
    }else{
      printf("Invertion status, if 0 > info:error! info:%d\n", info);
      printf("\n\n\n\n********************************** ERROR INVERTING MATRIX**************************************\n\n\n\n\n");
    }

/*
    for(i=0;i<17;i++) {
      for (j=0; j < no_meas*2; j++) {
        printf("K_matrix[%d]:%lf\t",(int)(i*(no_meas*2)+j),K_matrix[i*(no_meas*2)+j]);
      }
      printf("\n");
    }
    printf("\n");
*/


    /* 8. Formulate measurement innovations using (14.119) */
    //delta_z(1:no_meas,1) = GNSS_measurements(1:no_meas,1) - pred_meas(1:no_meas,1);
    //delta_z((no_meas + 1):(2 * no_meas),1) = GNSS_measurements(1:no_meas,2) -...
        //pred_meas(1:no_meas,2);

    double delta_z_sum=0.0;
    for (i = 0; i < no_meas; i++) {
      delta_z[i] = GNSS_measurements[i].P[0] - pred_meas[i*2];
      delta_z_sum+=delta_z[i];
    }
    if (fabs(delta_z_sum/no_meas)>150000.0) {
      printf("******************************************************************CORR.: %lf >150000.0, CLK: %lf \n", fabs(delta_z_sum/no_meas), x_est_propagated[15] );
    }

    for (i = no_meas; i < 2*no_meas; i++) {
      delta_z[i] = GNSS_measurements[i-no_meas].D[0] - pred_meas[(i-no_meas)*2+1];
    }

    /* 9. Update state estimates using (3.24) */
    matmul("NN",17,1,2*no_meas,1.0,K_matrix,delta_z,0.0,K_x_delta_z);
    for (i=0;i<17;i++) x_est_new[i] = x_est_propagated[i] + K_x_delta_z[i];

    /* 10. Update state estimation error covariance matrix using (3.25) */
    //P_matrix_new = (eye(17) - K_matrix * H_matrix) * P_matrix_propagated;
    matmul("NN",17,17,2*no_meas,1.0,K_matrix,H_matrix,0.0,K_x_H);
    for (i=0;i<17;i++) {
      for (j=0;j<17;j++) {
        I_less_k_x_H[i*17+j]= I_par[i*17+j] - K_x_H[i*17+j];
      }
    }

    matmul("NN",17,17,17,1.0,I_less_k_x_H,P_matrix_propagated,0.0,P_matrix_new);

    /* CLOSED-LOOP CORRECTION */
/*
    for(i=0;i<no_meas;i++) {
      for (j = 0; j < 2; j++) {
       printf("pred_meas[%d]:%lf\t",i*2+j,K_matrix[i*2+j]);
      }
      printf("\n");
    }
    printf("\n");
*/
    for (i = 0; i < no_meas; i++) {
      printf("P1_in:,D1_in: %lf, %lf\n",GNSS_measurements[i].P[0],GNSS_measurements[i].D[0]);
    }

    for(i=0;i<no_meas;i++) {
      for (j = 0; j < 2; j++) {
       printf("delta_z[%d]:%lf\t",i*2+j,delta_z[i*2+j]);
      }
      printf("\n");
    }
    printf("\n");
/**/

/*
    for (i=0; i < 2*no_meas;i++) {
      printf("pred_meas[%d]:%lf, ",i,pred_meas[i]);
    }
    printf("\n");
*/

/**/
    for (i=0; i < 17;i++) {
      printf("x_est_propagated[%d]:%lf, ",i,x_est_propagated[i]);
    }
    printf("\n");

    for (i=0; i < 17;i++) {
      printf("k_x_delta_z[%d]:%lf, ",i,K_x_delta_z[i]);
    }
    printf("\n");


/*
    for(i=0;i<no_meas*2;i++) {
      for (j=0; j < no_meas*2; j++) {
        printf("Q1 inv[%d]:%lf\t",(int)(i*(no_meas*2)+j),Q1[i*(no_meas*2)+j]);
      }
      printf("\n");
    }
    printf("\n");
/*

    for(i=0;i<17;i++) {
      for (j=0; j < 17; j++) {
        printf("P_matrix_prop[%d]:%lf\t",(int)(i*(17)+j),P_matrix_propagated[i*(17)+j]);
      }
      printf("\n");
    }
    printf("\n");
*/


  for (i=0;i<17;i++) printf("x_est_new.[%d]:%lf\n",i, x_est_new[i]);


    /* Correct attitude, velocity, and position using (14.7-9) */

    /* Quaternion attitude correction */
    printf("est_C_b_e_old: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
    est_C_b_e_old[0],est_C_b_e_old[1],est_C_b_e_old[2],est_C_b_e_old[3],\
    est_C_b_e_old[4],est_C_b_e_old[5],est_C_b_e_old[6],est_C_b_e_old[7],est_C_b_e_old[8]);

    Quaternion_attitude_errror_correction(x_est_new, est_C_b_e_old, est_C_b_e_new);


    Skew_symmetric(x_est_new, Skew_x_est_new);

    for (i=0;i<9;i++) I_less_Skew[i]= I[i]-Skew_x_est_new[i];

    //matmul("NN",3,3,3,1.0,I_less_Skew,est_C_b_e_old,0.0,est_C_b_e_new);

    for (i=3;i<6;i++) est_v_eb_e_new[i-3] = est_v_eb_e_old[i-3] - x_est_new[i];
    for (i=6;i<9;i++) est_r_eb_e_new[i-6] = est_r_eb_e_old[i-6] - x_est_new[i];

    /* Update IMU bias and GNSS receiver clock estimates */
    for (i=0;i<6;i++) est_IMU_bias_new[i] = est_IMU_bias_old[i] + x_est_new[i+9];
    est_clock_new[0] = x_est_new[15];
    est_clock_new[1] = x_est_new[16];

    printf("est_C_b_e_new: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
    est_C_b_e_new[0],est_C_b_e_new[1],est_C_b_e_new[2],est_C_b_e_new[3],\
    est_C_b_e_new[4],est_C_b_e_new[5],est_C_b_e_new[6],est_C_b_e_new[7],est_C_b_e_new[8]);


    /* Ends */

    /* Freeing memory */
    free(Phi_matrix); free(Q_prime_matrix); free(I);
    free(x_est_propagated); free(x_est_new);
    free(P_matrix_propagated); free(P_aux); free(Q_aux);
    free(P_matrix); free(u_as_e_T); free(pred_meas);
    free(ones); free(I_par);
    free(K_matrix); free(delta_r); free(delta_z); free(Q_);
    free(F1); free(Q1); free(K_x_delta_z);
    free(I_meas); free(R_matrix); free (K_matrix_inv);
    free(K_x_H); free(I_less_k_x_H);
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
 void LC_KF_Epoch(double *GNSS_r_eb_e, double *GNSS_v_eb_e, int no_meas,\
   float tor_s, double *est_C_b_e_old, double *est_v_eb_e_old, double *est_r_eb_e_old,\
   double *est_IMU_bias_old, double *P_matrix_old, double *meas_f_ib_b,\
   double est_L_b_old, LC_KF_config *LC_KF_config, double *est_C_b_e_new,\
   double *est_v_eb_e_new, double *est_r_eb_e_new, double *est_IMU_bias_new,\
   double *P_matrix_new){

    /* Constants (sone of these could be changed to inputs at a later date) */
    double c = 299792458; /* Speed of light in m/s */
    double omega_ie = OMGE;  /* Earth rotation rate in rad/s */
    double R_0 = RE_WGS84; /*WGS84 Equatorial radius in meters */
    double e = sqrt(e_2); /*WGS84 eccentricity                        */

    /* Begins */
    double omega_ie_vec[3], Omega_ie[9];
    double *Phi_matrix, *Q_prime_matrix, *x_est_propagated, *x_est_new;
    double *Q_, *Q_aux;
    double *P_matrix_propagated, *P_matrix, *P_aux, *u_as_e_T, *pred_meas;
    double *H_matrix, *R_matrix, *ones, *delta_r, *delta_z, approx_range, range;
    double range_rate, C_e_I[9], *I_meas, *K_matrix, *K_matrix_inv;
    double meas_f_ib_e[3], Skew_meas_f_ib_e[9], est_r_eb_e[3], delta_r_aux[3];
    double geocentric_radius, g[3], g_est_r_eb_e[9], *I, *I_par;
    double Omega_x_est_r_eb_old[3], Omega_x_GNSS_meas_r[3], v_plus_Omega_x_est_r[3];
    double v_plus_Omega_x_GNSS_meas_r[3], est_v_plus_Omega_x_est_r[3];
    double C_x_v[3],C_less_v_plus_Omega[3];
    double *F1, *Q1, *K_x_delta_z, *K_x_H, *I_less_k_x_H;
    double Skew_x_est_new[9], I_less_Skew[9];
    double rho, u[3], rate, vs[3], a[3], e_[3], cosel, azel[2*no_meas], pos[3];
    double E[9];
    int i,j, info;

    /* Initialize matrices and vectors */
    Phi_matrix=eye(15); Q_prime_matrix=zeros(15,15);
    x_est_propagated=zeros(15,1); x_est_new=zeros(15,1);
    P_matrix_propagated=mat(15,15); P_aux=mat(15,15); Q_aux=mat(15,15);
    P_matrix=mat(15,15); u_as_e_T=zeros(no_meas,3); pred_meas=zeros(no_meas,2);
    H_matrix=zeros((2 * no_meas),15); ones=mat(no_meas,1); I_meas=eye(no_meas);
    R_matrix=zeros((2 * no_meas),(2 * no_meas)); I_par=eye(15);
    K_matrix_inv=zeros((2 * no_meas),(2 * no_meas));
    K_matrix=mat(15,(2 * no_meas));
    delta_r=mat(no_meas,3); delta_z=mat(2 * no_meas,1); Q_=mat(15,15);
    F1=mat(15,2*no_meas); Q1=mat(2*no_meas,2*no_meas); K_x_delta_z=mat(15,1);
    K_x_H=mat(15,15); I_less_k_x_H=mat(15,15); I = eye(3);


    for (i=0;i<no_meas;i++) ones[i]=1.0;

    /* Skew symmetric matrix of Earth rate */
    omega_ie_vec[0]=0; omega_ie_vec[1]=0; omega_ie_vec[2]=OMGE;
    Skew_symmetric(omega_ie_vec, Omega_ie);

    /* SYSTEM PROPAGATION PHASE */

    /* 1. Determine transition matrix using (14.50) (first-order approx) */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        Phi_matrix[i*15+j] = Phi_matrix[i*15+j] - Omega_ie[i*3+j] * tor_s;
      }
    }

    for (i = 0; i < 3; i++) {
      for (j = 12; j < 15; j++) {
        Phi_matrix[i*15+j] = est_C_b_e_old[i*3+(j-12)] * tor_s;
      }
    }

    matmul("NN",3,1,3,1.0,est_C_b_e_old,meas_f_ib_b,0.0, meas_f_ib_e);
    Skew_symmetric(meas_f_ib_b,Skew_meas_f_ib_e);

    for (i = 3; i < 6; i++) {
      for (j = 0; j < 3; j++) {
        Phi_matrix[i*15+j] = -tor_s * Skew_meas_f_ib_e[(i-3)*3+j];
      }
    }

    for (i = 3; i < 6; i++) {
      for (j = 3; j < 6; j++) {
        Phi_matrix[i*15+j] = Phi_matrix[i*15+j] - 2 * Omega_ie[(i-3)*3+(j-3)] * tor_s;
      }
    }

    geocentric_radius = R_0 / sqrt(1 - pow((e * sin(est_L_b_old)),2)) *\
        sqrt(pow(cos(est_L_b_old),2) + pow((1 - e*e),2) * pow(sin(est_L_b_old),2)); /* from (2.137)*/
    Gravity_ECEF(est_r_eb_e_old, g);
    matmul("NT",3,3,1,1.0,g,est_r_eb_e,0.0, g_est_r_eb_e);
    for (i = 3; i < 6; i++) {
      for (j = 6; j < 9; j++) {
        Phi_matrix[i*15+j] = -tor_s * 2 /\
         geocentric_radius * g_est_r_eb_e[(i-3)*3+(j-6)] / (norm(est_r_eb_e_old,3));
      }
    }

    for (i = 3; i < 6; i++) {
      for (j = 9; j < 12; j++) {
        Phi_matrix[i*15+j] = est_C_b_e_old[(i-3)*3+(j-9)] * tor_s;
      }
    }

    for (i = 6; i < 9; i++) {
      for (j = 3; j < 6; j++) {
        Phi_matrix[i*15+j] = I[(i-6)*3+(j-3)] * tor_s;
      }
    }


    /* 2. Determine approximate system noise covariance matrix using (14.82) */
    //Q_prime_matrix(1:3,1:3) = eye(3) * LC_KF_config.gyro_noise_PSD * tor_s;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        Q_prime_matrix[i*15+j] = (i==j?I[(i)*3+(j)] *\
         LC_KF_config->gyro_noise_PSD * tor_s:0.0);
      }
    }
    //Q_prime_matrix(4:6,4:6) = eye(3) * LC_KF_config.accel_noise_PSD * tor_s;
    for (i = 3; i < 6; i++) {
      for (j = 3; j < 6; j++) {
        Q_prime_matrix[i*15+j] = (i==j?I[(i-3)*3+(j-3)] *\
         LC_KF_config->accel_noise_PSD * tor_s:0.0);
      }
    }
   //Q_prime_matrix(10:12,10:12) = eye(3) * LC_KF_config.accel_bias_PSD * tor_s;
    for (i = 9; i < 12; i++) {
      for (j = 9; j < 12; j++) {
        Q_prime_matrix[i*15+j] = (i==j?I[(i-9)*3+(j-9)] *\
         LC_KF_config->accel_bias_PSD * tor_s:0.0);
      }
    }
   //Q_prime_matrix(13:15,13:15) = eye(3) * LC_KF_config.gyro_bias_PSD * tor_s;
    for (i = 12; i < 15; i++) {
      for (j = 12; j < 15; j++) {
        Q_prime_matrix[i*15+j] = (i==j?I[(i-12)*3+(j-12)] *\
         LC_KF_config->gyro_bias_PSD * tor_s:0.0);
      }
    }

    /* 3. Propagate state estimates using (3.14) noting that all states are zero
       due to closed-loop correction.  */
   for (i = 0; i < 15; i++) x_est_propagated[i] = 0.0;

    /* 4. Propagate state estimation error covariance matrix using (3.46) */
    //P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *\
    //    Phi_matrix' + 0.5 * Q_prime_matrix;
    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        P_aux[i*15+j]= P_matrix_old[i*15+j]+0.5*Q_prime_matrix[i*15+j];
      }
    }

    matmul("NT",15,15,15,1.0,P_aux,Phi_matrix,0.0,Q_aux); /* (Pp + 0.5*Q)*PHI^T */
    matmul("NN",15,15,15,1.0,Phi_matrix,Q_aux,0.0,Q_); /* PHI*(Pp + 0.5*Q)*PHI^T */

    for (i = 0; i < 15; i++) {
      for (j = 0; j < 15; j++) {
        P_matrix_propagated[i*15+j]= Q_[i*15+j]+0.5*Q_prime_matrix[i*15+j]; /*P_ = PHI*(Pp + 0.5*Q)*PHI^T + 0.5*Q*/
      }
    }

    /* MEASUREMENT UPDATE PHASE */

    /* 5. Set-up measurement matrix using (14.115)  */

   //H_matrix(1:3,7:9) = -eye(3);
   for (i = 0; i < 3; i++) {
     for (j = 6; j < 9; j++) {
       H_matrix[i*15+j] = (i==j?-1.0:0.0);
     }
   }

   //H_matrix(4:6,4:6) = -eye(3);
   for (i = 3; i < 6; i++) {
     for (j = 3; j < 6; j++) {
       H_matrix[i*15+j] = (i==j?-1.0:0.0);
     }
   }

   /* 6. Set-up measurement noise covariance matrix assuming all components of
    GNSS position and velocity are independent and have equal variance. */

    //R_matrix(1:3,1:3) = eye(3) * LC_KF_config.pos_meas_SD^2;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        R_matrix[i*(no_meas)+j] = I_meas[i*(no_meas)+j]*\
        pow(LC_KF_config->pos_meas_SD,2);
      }
    }

    //R_matrix(1:3,4:6) = zeros(3);
    for (i = 0; i < 3; i++) {
      for (j = 3; j < 6; j++) {
        R_matrix[i*(no_meas)+j] = 0.0;
      }
    }

    //R_matrix(4:6,1:3) = zeros(3);
    for (i = 3; i < 6; i++) {
      for (j = 0; j < 3; j++) {
        R_matrix[i*(no_meas)+j] = 0.0;
      }
    }

    //R_matrix(4:6,4:6) = eye(3) * LC_KF_config.vel_meas_SD^2;
    for (i = 3; i < 6; i++) {
      for (j = 3; j < 6; j++) {
        R_matrix[i*(no_meas)+j] = I_meas[(i-3)*(no_meas)+(j-3)] *\
         pow(LC_KF_config->vel_meas_SD,2);
      }
    }

    /* 7. Calculate Kalman gain using (3.21) */
   //K_matrix = P_matrix_propagated * H_matrix' * inv(H_matrix *...
    //    P_matrix_propagated * H_matrix' + R_matrix);
    matmul("NT",15,6,15,1.0,P_matrix_propagated,H_matrix,0.0,F1);
    matmul("NN",6,6,15,1.0,H_matrix,F1,0.0,Q1);
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {
        K_matrix_inv[i*(6)+j] = Q1[i*(6)+j]+ R_matrix[i*(6)+j];
      }
    }


    if (!(info=matinv(K_matrix_inv,6))) {
      printf("Invertion status, if 0 > info:error! info:%d\n", info);
        matmul("NN",15,6,6,1.0,F1,K_matrix_inv,0.0,K_matrix);
        printf("\n\n\n\n********************************** INVERTING MATRIX  **************************************\n\n\n\n\n");
    }else{
      printf("Invertion status, if 0 > info:error! info:%d\n", info);
      printf("\n\n\n\n********************************** ERROR INVERTING MATRIX**************************************\n\n\n\n\n");
    }


    /* 8. Formulate measurement innovations using (14.102), noting that zero
    % lever arm is assumed here  */

    for (i = 0; i < 3; i++) {
      delta_z[i] = GNSS_r_eb_e[i] -est_r_eb_e_old[i];
    }

    for (i = 3; i < 6; i++) {
      delta_z[i] = GNSS_v_eb_e[i-3] -est_v_eb_e_old[i-3];
    }


    /* 9. Update state estimates using (3.24) */
    matmul("NN",15,1,6,1.0,K_matrix,delta_z,0.0,K_x_delta_z);
    for (i=0;i<15;i++) x_est_new[i] = x_est_propagated[i] + K_x_delta_z[i];

    /* 10. Update state estimation error covariance matrix using (3.25) */
    //P_matrix_new = (eye(15) - K_matrix * H_matrix) * P_matrix_propagated;
    matmul("NN",15,15,6,1.0,K_matrix,H_matrix,0.0,K_x_H);
    for (i=0;i<15;i++) {
      for (j=0;j<15;j++) {
        I_less_k_x_H[i*15+j]= I_par[i*15+j] - K_x_H[i*15+j];
      }
    }

    matmul("NN",15,15,15,1.0,I_less_k_x_H,P_matrix_propagated,0.0,P_matrix_new);

    /* CLOSED-LOOP CORRECTION */


    for(i=0;i<no_meas;i++) {
       printf("delta_z[%d]:%lf\t",i,delta_z[i]);
    }
    printf("\n");

   for (i=0;i<15;i++) printf("x_est_new.[%d]:%lf\n",i, x_est_new[i]);


    /* Correct attitude, velocity, and position using (14.7-9) */

    /* Quaternion attitude correction */
    printf("est_C_b_e_old: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
    est_C_b_e_old[0],est_C_b_e_old[1],est_C_b_e_old[2],est_C_b_e_old[3],\
    est_C_b_e_old[4],est_C_b_e_old[5],est_C_b_e_old[6],est_C_b_e_old[7],est_C_b_e_old[8]);

    Quaternion_attitude_errror_correction(x_est_new, est_C_b_e_old, est_C_b_e_new);


    Skew_symmetric(x_est_new, Skew_x_est_new);

    for (i=0;i<9;i++) I_less_Skew[i]= I[i]-Skew_x_est_new[i];

    //matmul("NN",3,3,3,1.0,I_less_Skew,est_C_b_e_old,0.0,est_C_b_e_new);

    for (i=3;i<6;i++) est_v_eb_e_new[i-3] = est_v_eb_e_old[i-3] - x_est_new[i];
    for (i=6;i<9;i++) est_r_eb_e_new[i-6] = est_r_eb_e_old[i-6] - x_est_new[i];

    /* Update IMU bias and GNSS receiver clock estimates */
    for (i=0;i<6;i++) est_IMU_bias_new[i] = est_IMU_bias_old[i] + x_est_new[i+9];



    /* Ends */

    /* Freeing memory */
    free(Phi_matrix); free(Q_prime_matrix); free(I);
    free(x_est_propagated); free(x_est_new);
    free(P_matrix_propagated); free(P_aux); free(Q_aux);
    free(P_matrix); free(u_as_e_T); free(pred_meas);
    free(ones); free(I_par);
    free(K_matrix); free(delta_r); free(delta_z); free(Q_);
    free(F1); free(Q1); free(K_x_delta_z);
    free(I_meas); free(R_matrix); free (K_matrix_inv);
    free(K_x_H); free(I_less_k_x_H);
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
void Nav_equations_ECEF(float tor_i,\
    double *old_r_eb_e,double *old_v_eb_e, double *old_C_b_e, \
    double *f_ib_b, double *omega_ib_b, double *r_eb_e, \
    double *v_eb_e, double *C_b_e){

      /**/
      printf("NAVIGATION EQUATIONS\nOLD_P: %lf, %lf, %lf\n",old_r_eb_e[0],old_r_eb_e[1],old_r_eb_e[2] );
      printf("OLD_V: %lf, %lf, %lf\n",old_v_eb_e[0],old_v_eb_e[1],old_v_eb_e[2] );
      printf("ACC: %lf, %lf, %lf\n",f_ib_b[0],f_ib_b[1],f_ib_b[2] );
      printf("GYR: %lf, %lf, %lf\n",omega_ib_b[0],omega_ib_b[1],omega_ib_b[2] );


     /* parameters  */
     double omega_ie = 7.292115E-5;  // Earth rotation rate (rad/s)
     double alpha_ie, mag_alpha;
     double C_Earth[9], alpha_ib_b[3], Alpha_ib_b[9];
     double Alpha_squared[9], second_term[9], first_term[9], C_new_old[9];
     double I[9]={1,0,0,0,1,0,0,0,1};
     double C_aux[9], first_term_2[9], second_term_2[9], C_b_e_Cbb[9];
     double ave_C_b_e[9], Cbb[9], alpha_ie_vec[3],Alpha_ie[9], last_term[9];
     double f_ib_e[3], omega_ie_vec[3], Omega_ie[9], g[3], Omega_v_eb_e[3];
     double new_q_b_e[4], new_C_b_e[9];
     int i;

     /* Begins     */

     /* ATTITUDE UPDATE  */
     /* From (2.145) determine the Earth rotation over the update interval
      C_Earth = C_e_i' * old_C_e_i  */
     alpha_ie = omega_ie * tor_i;
     C_Earth[0]=cos(alpha_ie); C_Earth[1]=sin(alpha_ie); C_Earth[2]=0;
     C_Earth[3]=-sin(alpha_ie); C_Earth[4]=cos(alpha_ie); C_Earth[5]=0;
     C_Earth[6]=0; C_Earth[7]=0; C_Earth[8]=1;

     omega_ie_vec[0] = 0; omega_ie_vec[1] = 0; omega_ie_vec[2] = omega_ie;

     /* Calculate attitude increment, magnitude, and skew-symmetric matrix  */
     for(i=0;i<3;i++) alpha_ib_b[i] = omega_ib_b[i] * tor_i;
     mag_alpha = norm(alpha_ib_b,3);
     Skew_symmetric(alpha_ib_b, Alpha_ib_b);

     /* Attitude Update using Quaternion algebra */
     attitude_update(alpha_ib_b, omega_ie_vec, tor_i, old_C_b_e, new_q_b_e);

     Quaternion_to_DCM(new_q_b_e, new_C_b_e);

          printf("new_C_b_e: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
          new_C_b_e[0],new_C_b_e[1],new_C_b_e[2],new_C_b_e[3],\
          new_C_b_e[4],new_C_b_e[5],new_C_b_e[6],new_C_b_e[7],new_C_b_e[8]);

     /* Obtain coordinate transformation matrix from the new attitude w.r.t. an
      inertial frame to the old using Rodrigues' formula, (5.73)  */
      matmul("NT", 3, 3, 3, 1.0, Alpha_ib_b, Alpha_ib_b, 0.0, Alpha_squared);
      for(i=0;i<9;i++) second_term[i] = (1 - cos(mag_alpha)) /\
      (mag_alpha*mag_alpha) * Alpha_squared[i];
      for(i=0;i<9;i++) first_term[i] = I[i] + sin(mag_alpha) /\
      mag_alpha*Alpha_ib_b[i];

     if (mag_alpha>1.E-8){
       for(i=0;i<9;i++) C_new_old[i] = first_term[i]+second_term[i];
     }else{
       for (i=0;i<9;i++) C_new_old[i] = I[i]+Alpha_ib_b[i];
     }// end if mag_alpha

     /* Update attitude using (5.75)  */
     matmul("NN", 3, 3, 3, 1.0, C_Earth, old_C_b_e, 0.0, C_aux);
     matmul("NN", 3, 3, 3, 1.0, C_aux, C_new_old, 0.0, C_b_e);

     /* SPECIFIC FORCE FRAME TRANSFORMATION
     % Calculate the average body-to-ECEF-frame coordinate transformation
     % matrix over the update interval using (5.84) and (5.85)  */
     for(i=0;i<9;i++) first_term_2[i] = I[i] + second_term[i];
     for(i=0;i<9;i++) second_term_2[i] = ((1 - (sin(mag_alpha)/mag_alpha))/\
     (mag_alpha*mag_alpha))*Alpha_squared[i];
     for(i=0;i<9;i++) Cbb[i] = first_term_2[i]+second_term_2[i];
     alpha_ie_vec[0] = 0; alpha_ie_vec[1] = 0; alpha_ie_vec[2] = alpha_ie;
     Skew_symmetric(alpha_ie_vec,Alpha_ie);
     matmul("NN", 3, 3, 3, 0.5, Alpha_ie, old_C_b_e, 0.0, last_term);

     if (mag_alpha>1.E-8){
       matmul("NN", 3, 3, 3, 1.0, old_C_b_e, Cbb, 0.0, C_b_e_Cbb);
       for(i=0;i<9;i++) ave_C_b_e[i] = C_b_e_Cbb[i] - last_term[i];
     }else{
       for (i=0;i<9;i++) ave_C_b_e[i] = old_C_b_e[i] - last_term[i];
     } //if mag_alpha

     /* Transform specific force to ECEF-frame resolving axes using (5.85) */
     matmul("NN", 3, 1, 3, 1.0, ave_C_b_e, f_ib_b, 0.0, f_ib_e);

     /* UPDATE VELOCITY
     % From (5.36), */
     Skew_symmetric(omega_ie_vec,Omega_ie);
     Gravity_ECEF(old_r_eb_e, g);
     matmul("NN", 3, 1, 3, 1.0, Omega_ie, old_v_eb_e, 0.0, Omega_v_eb_e);
     for (i=0;i<3;i++) v_eb_e[i] = old_v_eb_e[i] + tor_i * (f_ib_e[i] + g[i] - \
         2 * Omega_v_eb_e[i]);

     /* UPDATE CARTESIAN POSITION
     % From (5.38), */
     for (i=0;i<3;i++) r_eb_e[i] = old_r_eb_e[i] + (v_eb_e[i] + old_v_eb_e[i]) * 0.5 * tor_i;
/**/
printf("tor_i: %f \n",tor_i);
    // printf("f_ib_e: %lf, %lf, %lf\n",f_ib_e[0],f_ib_e[1],f_ib_e[2] );
     printf("C_b_e: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
     C_b_e[0],C_b_e[1],C_b_e[2],C_b_e[3],\
     C_b_e[4],C_b_e[5],C_b_e[6],C_b_e[7],C_b_e[8]);
/*     printf("ave_C_b_e: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
     ave_C_b_e[0],ave_C_b_e[1],ave_C_b_e[2],ave_C_b_e[3],\
     ave_C_b_e[4],ave_C_b_e[5],ave_C_b_e[6],ave_C_b_e[7],ave_C_b_e[8]);
     printf("Omega_v_eb_e: %lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
     Omega_v_eb_e[0],Omega_v_eb_e[1],Omega_v_eb_e[2],Omega_v_eb_e[3],\
     Omega_v_eb_e[4],Omega_v_eb_e[5],Omega_v_eb_e[6],Omega_v_eb_e[7],Omega_v_eb_e[8]);
     printf("NAV_POS_NEW: %lf, %lf, %lf\n",r_eb_e[0],r_eb_e[1],r_eb_e[2] );
     printf("NAV_V_NEW: %lf, %lf, %lf\n",v_eb_e[0],v_eb_e[1],v_eb_e[2] );
*/
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
void Loosely_coupled_INS_GNSS(INS_measurements *INS_measurements, int no_meas,\
  PVAT_solution *pvat_old, int no_par, float old_time, float time_last_GNSS, \
  initialization_errors *initialization_errors, IMU_errors *IMU_errors,\
  GNSS_config *GNSS_config, LC_KF_config *LC_KF_config, PVAT_solution *pvat_new,\
  double *out_errors, double *out_IMU_bias_est, double *out_KF_SD, \
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
  double est_IMU_bias[6]={0};
  double quant_residuals[6]={0};
  double *P_matrix, *P_matrix_new;
  float GNSS_epoch, tor_i, time, tor_s;
  int i,j;
  double est_C_b_e_new[9], est_C_b_n_T[9], est_v_eb_e_new[3], est_r_eb_e_new[3];
  double est_IMU_bias_new[6], est_clock_new[2];
  double q_nb[4], llh[3];

  printf("\n *****************  LC_INS/GNSS BEGINS ************************\n");
/*
  printf("P: %lf, %lf, %lf\n",pvat_old->latitude*R2D,pvat_old->longitude*R2D,pvat_old->height );
  printf("V: %lf, %lf, %lf\n",pvat_old->ned_velocity[0],pvat_old->ned_velocity[1],pvat_old->ned_velocity[2] );
  printf("A: %lf, %lf, %lf\n",pvat_old->euler_angles[0],pvat_old->euler_angles[1],pvat_old->euler_angles[2] );*/

  /* Initialize true navigation solution */
  time = INS_measurements->sec;
  true_L_b = pvat_old->latitude;
  true_lambda_b = pvat_old->longitude;
  true_h_b = pvat_old->height;
  for (i=0;i<3;i++) true_v_eb_n[i] = pvat_old->ned_velocity[i];
  for (i=0;i<3;i++) true_eul_nb[i] = pvat_old->euler_angles[i];
  Euler_to_CTM(true_eul_nb, true_C_b_n);

  printf("true_C_b_n:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
   true_C_b_n[0],true_C_b_n[1],true_C_b_n[2],\
   true_C_b_n[3],true_C_b_n[4],true_C_b_n[5],\
   true_C_b_n[6],true_C_b_n[7],true_C_b_n[8]);

   printf("NORM OF CBN should be ~ 1: |CBN|: %lf\n",norm(true_C_b_n,9) );

   Euler_to_quaternion(true_eul_nb, q_nb);
   Quaternion_to_DCM(q_nb, true_C_b_n);

   printf("true_C_b_n:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
    true_C_b_n[0],true_C_b_n[1],true_C_b_n[2],\
    true_C_b_n[3],true_C_b_n[4],true_C_b_n[5],\
    true_C_b_n[6],true_C_b_n[7],true_C_b_n[8]);

   printf("NORM OF CBN should be ~ 1: |CBN|: %lf\n",norm(true_C_b_n,9) );

  printf("ATt or Euler_nb from pvat_old: %lf, %lf, %lf\n",true_eul_nb[0],true_eul_nb[1],true_eul_nb[2]);
  printf("V_eb_n from pvat_old: %lf, %lf, %lf\n",true_v_eb_n[0],true_v_eb_n[1],true_v_eb_n[2]);
  printf("LLH from pvat_old: %lf, %lf, %lf\n",true_L_b*R2D,true_lambda_b*R2D, true_h_b);

  NED_to_ECEF(&true_L_b,&true_lambda_b,&true_h_b,true_v_eb_n,true_C_b_n,\
    old_true_r_eb_e,old_true_v_eb_e,old_true_C_b_e);

    printf("old_true_C_b_e:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
     old_true_C_b_e[0],old_true_C_b_e[1],old_true_C_b_e[2],\
     old_true_C_b_e[3],old_true_C_b_e[4],old_true_C_b_e[5],\
     old_true_C_b_e[6],old_true_C_b_e[7],old_true_C_b_e[8]);

  /* Initialize other matrices */

  /* Determine Least-squares GNSS position solution */
  old_est_llh[0]=true_L_b; old_est_llh[1]=true_lambda_b; old_est_llh[2]=true_h_b;
  pos2ecef(old_est_llh, old_est_r_eb_e);
//  GNSS_LS_position_velocity(GNSS_measurements,no_GNSS_meas,\
  //  GNSS_config->init_est_r_ea_e,[0;0;0],old_est_r_eb_e,old_est_v_eb_e,est_clock);


  for (j=0;j<3;j++) old_est_v_eb_e[j]=old_true_v_eb_e[j];
  pv_ECEF_to_NED(old_est_r_eb_e,old_est_v_eb_e, &old_est_L_b,\
    &old_est_lambda_b,&old_est_h_b,old_est_v_eb_n);
  est_L_b = old_est_L_b;

  /* Initialize estimated attitude solution */
  //Initialize_NED_attitude(true_C_b_n, initialization_errors, old_est_C_b_n);
  for (j=0;j<9;j++) old_est_C_b_n[j]=true_C_b_n[j];

  NED_to_ECEF(&old_est_L_b,&old_est_lambda_b,&old_est_h_b,old_est_v_eb_n,\
    old_est_C_b_n,old_true_r_eb_e,old_true_v_eb_e,old_true_C_b_e);

  /* Initialize Kalman filter P matrix and IMU bias states */
  P_matrix=zeros(no_par,no_par); // in this case 15x15
  P_matrix_new=zeros(no_par,no_par);
  Initialize_LC_P_matrix(LC_KF_config, P_matrix);
  if (out_errors[0] > 0.0 && out_errors[1]>0.0 && out_errors[2]>0.0) {
    for (i = 0; i < no_par; i++) {
      for (j = 0; j < no_par; j++) {
        if (i==j) {
          P_matrix[i*no_par+j]=out_errors[j];
        }
      }
    }
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

  if (old_time<=0.000) {
    tor_i=1.0; //1 second of time interval
    time_last_GNSS=time-tor_i;
  }

  printf("TIMES: time, tor_i: %f, %f, \
  KF TIME DIFF: %f\n", time, tor_i, time - time_last_GNSS );

    /* Correct IMU errors */
    for (i=0;i<3;i++) meas_f_ib_b[i]=INS_measurements->f_ib_b[i];
    for (i=0;i<3;i++) meas_omega_ib_b[i]=INS_measurements->omega_ib_b[i];
    for (i=0;i<3;i++) meas_f_ib_b[i] = meas_f_ib_b[i] - est_IMU_bias[i];
    for (i=0;i<3;i++) meas_omega_ib_b[i] = meas_omega_ib_b[i] - est_IMU_bias[i+3];

  printf("BEFORE NAV est_r_eb_e:%lf, %lf, %lf\n", old_est_r_eb_e[0],\
  old_est_r_eb_e[1],old_est_r_eb_e[2]);
  /**/
      printf("est_v_eb_e:%lf, %lf, %lf\n", old_est_v_eb_e[0],old_est_v_eb_e[1],\
      old_est_v_eb_e[2]);
      //printf("est_v_eb_e:%lf, %lf, %lf\n", est_v_eb_e[0],est_v_eb_e[1],est_v_eb_e[2]);
      printf("old_true_C_b_e:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
       old_true_C_b_e[0],old_true_C_b_e[1],old_true_C_b_e[2],\
       old_true_C_b_e[3],old_true_C_b_e[4],old_true_C_b_e[5],\
       old_true_C_b_e[6],old_true_C_b_e[7],old_true_C_b_e[8]);

  printf("NORM OF old_true_CBE should be ~ 1: |CBE|: %lf\n",norm(old_true_C_b_e,9) );

    /* Update estimated navigation solution */
printf("\n *****************  NAV_EQUATIONS BEGINS ************************\n");
    Nav_equations_ECEF(tor_i,\
        old_est_r_eb_e, old_est_v_eb_e,old_true_C_b_e, meas_f_ib_b,\
        meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
printf("\n *****************  NAV_EQUATIONS ENDS **************************\n");

    /* Determine whether to update GNSS simulation and run Kalman filter */
/**/
   printf("BEFORE TCKF\nest_r_eb_e:%lf, %lf, %lf\n", est_r_eb_e[0],\
   est_r_eb_e[1],est_r_eb_e[2]);
/**/
    printf("est_v_eb_e:%lf, %lf, %lf\n", est_v_eb_e[0],est_v_eb_e[1],\
    est_v_eb_e[2]);
    //printf("est_v_eb_e:%lf, %lf, %lf\n", est_v_eb_e[0],est_v_eb_e[1],est_v_eb_e[2]);
    printf("est_C_b_e:%lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
     est_C_b_e[0],est_C_b_e[1],est_C_b_e[2],\
     est_C_b_e[3],est_C_b_e[4],est_C_b_e[5],\
     est_C_b_e[6],est_C_b_e[7],est_C_b_e[8]);

     printf("NORM OF est_CBE should be ~ 1: |CBE|: %lf\n",norm(est_C_b_e,9) );

     ecef2pos(est_r_eb_e, llh);

     printf("LLH: %lf, %lf, %lf\n", llh[0]*R2D, llh[1]*R2D, llh[2]);


 /* Only process if GNSS and INS timne differences are within 0.4 seconds */
//    if ( fabs(time - time_last_GNSS) < 0.4) {
printf("TIME DIFF OF GNSS: %f and IMU: %f = %f \n", time, INS_measurements->sec, fabs(time - INS_measurements->sec));
  if ( fabs(time - INS_measurements->sec) < 0.0001) {
      /* KF time interval */
      tor_s = time - time_last_GNSS;

      /* Run Integration Kalman filter */
printf("\n *******************  LC_KF_EPOCH BEGINS ************************\n");
    /**/  LC_KF_Epoch(old_true_r_eb_e,old_true_v_eb_e,no_meas, tor_s, est_C_b_e,\
            est_v_eb_e,est_r_eb_e,est_IMU_bias, P_matrix,\
            meas_f_ib_b,est_L_b,LC_KF_config, est_C_b_e_new,\
            est_v_eb_e_new, est_r_eb_e_new, est_IMU_bias_new,\
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

      if ( isnan(est_IMU_bias_new[0]) || isnan(est_IMU_bias_new[2]) ) {
         /* code */
       }else{
         fprintf(out_IMU_bias_file,"%f %lf %lf %lf %lf %lf %lf\n", time, \
         est_IMU_bias_new[0],est_IMU_bias_new[1],est_IMU_bias_new[2],\
         est_IMU_bias_new[3],est_IMU_bias_new[4], est_IMU_bias_new[5]);
       }


      for (i = 0; i < 6; i++) out_IMU_bias_est[i]=est_IMU_bias_new[i];

      /* Generate KF uncertainty output record */
      fprintf(out_KF_SD_file,"%f\t", time);
      for (i=0;i<no_par;i++) {
        for (j=0;j<no_par;j++) {
          if (i==j) {
            out_errors[j]=P_matrix_new[i*no_par+j];
            fprintf(out_KF_SD_file,"%lf\t", sqrt(P_matrix_new[i*no_par+j]) );
          }
        }
      }
      fprintf(out_KF_SD_file,"\n");


    }else{
      for (i=0;i<3;i++) est_r_eb_e_new[i]= est_r_eb_e[i];
      for (i=0;i<3;i++) est_v_eb_e_new[i]= est_v_eb_e[i];
      for (i=0;i<9;i++) est_C_b_e_new[i] = est_C_b_e[i];
      for (i=0;i<2;i++) est_clock_new[i] = est_clock[i];

    }

      /**/
      printf("After TCKF\nest_r_eb_e:%lf, %lf, %lf\n",\
      est_r_eb_e_new[0],est_r_eb_e_new[1],est_r_eb_e_new[2]);

      ecef2pos(est_r_eb_e_new, llh);
      printf("LLH: %lf, %lf, %lf\n", llh[0]*R2D, llh[1]*R2D, llh[2]);
      /**/
      printf("est_v_eb_e:%lf, %lf, %lf\n", est_v_eb_e_new[0],\
      est_v_eb_e_new[1],est_v_eb_e_new[2]);
          //printf("est_v_eb_e:%lf, %lf, %lf\n", est_v_eb_e[0],est_v_eb_e[1],est_v_eb_e[2]);
      printf("est_C_b_e_new:%lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
           est_C_b_e_new[0],est_C_b_e_new[1],est_C_b_e_new[2],\
           est_C_b_e_new[3],est_C_b_e_new[4],est_C_b_e_new[5],\
           est_C_b_e_new[6],est_C_b_e_new[7],est_C_b_e_new[8]);

   printf("NORM OF est_CBE_NEW should be ~ 1: |CBE|: %lf\n",norm(est_C_b_e_new,9) );


    /* Convert navigation solution to NED  */
    ECEF_to_NED(est_r_eb_e_new,est_v_eb_e_new,est_C_b_e_new,\
      &est_L_b,&est_lambda_b,&est_h_b,\
      est_v_eb_n,est_C_b_n);

      printf("est_C_b_n:%lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",\
       est_C_b_n[0],est_C_b_n[1],est_C_b_n[2],\
       est_C_b_n[3],est_C_b_n[4],est_C_b_n[5],\
       est_C_b_n[6],est_C_b_n[7],est_C_b_n[8]);

       printf("EST_CBN NORM |CBN|: %lf\n",norm(est_C_b_n,9) );

   for (i=0;i<3;i++){
     for (j=0;j<3;j++) {
       est_C_b_n_T[i*3+j]=est_C_b_n[j*3+i];
     }
   }

  double q_b_n[4];

   DCM_to_quaternion(est_C_b_n, q_b_n);
   Quaternion_to_euler(q_b_n, pvat_new->euler_angles);

   printf("Euler from Quaternions: %lf, %lf, %lf\n",pvat_new->euler_angles[0],\
 pvat_new->euler_angles[1],pvat_new->euler_angles[2] );

    pvat_new->latitude = est_L_b;
    pvat_new->longitude = est_lambda_b;
    pvat_new->height = est_h_b;
    for (i=0;i<3;i++) pvat_new->ned_velocity[i] = est_v_eb_n[i];
    CTM_to_Euler(pvat_new->euler_angles, est_C_b_n);
    pvat_new->time=time;

/*
    printf("Lat, lon, h: %lf, %lf, %lf\n",est_L_b,est_lambda_b,est_h_b );
    printf("LAT_lon_h_new: %lf, %lf, %lf\n",pvat_new->latitude,pvat_new->longitude,pvat_new->height );
    printf("est_r_eb_e_new: %lf, %lf, %lf\n",est_r_eb_e_new[0],est_r_eb_e_new[1],est_r_eb_e_new[2] );
    printf("est_v_eb_e_new: %lf, %lf, %lf\n",est_v_eb_e_new[0],est_v_eb_e_new[1],est_v_eb_e_new[2] );

    printf("V_new: %lf, %lf, %lf\n",pvat_new->ned_velocity[0],pvat_new->ned_velocity[1],pvat_new->ned_velocity[2] );
    printf("A_new: %lf, %lf, %lf\n",pvat_new->euler_angles[0],pvat_new->euler_angles[1],pvat_new->euler_angles[2] );
    printf("est_C_b_n: %lf, %lf, %lf\n",est_C_b_n[0],est_C_b_n[1],est_C_b_n[2] );
    printf("est_C_b_n: %lf, %lf, %lf\n",est_C_b_n[3],est_C_b_n[4],est_C_b_n[5] );
    printf("est_C_b_n: %lf, %lf, %lf\n",est_C_b_n[6],est_C_b_n[7],est_C_b_n[8] );/**/
/**/
    /* Generate output profile record
         out_profile(epoch,1) = time;
         out_profile(epoch,2) = est_L_b;
         out_profile(epoch,3) = est_lambda_b;
         out_profile(epoch,4) = est_h_b;
         out_profile(epoch,5:7) = est_v_eb_n'';
         out_profile(epoch,8:10) = CTM_to_Euler(est_C_b_n')';
    */

    /* Determine errors and generate output record
         [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
             est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
             true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
    /*
         out_errors(epoch,1) = time;
         out_errors(epoch,2:4) = delta_r_eb_n';
         out_errors(epoch,5:7) = delta_v_eb_n';
         out_errors(epoch,8:10) = delta_eul_nb_n'';
    */

    /* Reset old values
         old_time = time;
         for (i=0;i<3;i++) old_true_r_eb_e[i] = true_r_eb_e[i];
         for (i=0;i<3;i++) old_true_v_eb_e[i] = true_v_eb_e[i];
         for (i=0;i<9;i++) old_true_C_b_e[i] = true_C_b_e[i];
         for (i=0;i<3;i++) old_est_r_eb_e[i] = est_r_eb_e[i];
         for (i=0;i<3;i++) old_est_v_eb_e[i] = est_v_eb_e[i];
         for (i=0;i<9;i++) old_true_C_b_e[i] = est_C_b_e[i];
    */


  //} //end epoch loop

  /* Free memory */
  free(P_matrix); free(P_matrix_new);
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
void Tightly_coupled_INS_GNSS(INS_measurements *INS_measurements, \
  GNSS_measurements *GNSS_measurements,\
  int no_GNSS_meas, PVAT_solution *pvat_old, int no_par, float old_time,\
  float time_last_GNSS,  double *clock_offset_drift,\
  initialization_errors *initialization_errors, IMU_errors *IMU_errors,\
  GNSS_config *GNSS_config, TC_KF_config *TC_KF_config, PVAT_solution *pvat_new,\
  double *out_errors, double *out_IMU_bias_est, double *out_clock,\
  double *out_KF_SD, double *meas_f_ib_b, double *meas_omega_ib_b)
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
  double est_IMU_bias[6]={0};
  double quant_residuals[6]={0};
  double *P_matrix, *P_matrix_new;
  double GNSS_epoch, tor_i, time, tor_s;
  int i,j;
  double est_C_b_e_new[9], est_C_b_n_T[9], est_v_eb_e_new[3], est_r_eb_e_new[3];
  double est_IMU_bias_new[6], est_clock_new[2];
  double q_nb[4], llh[3];

  printf("\n *****************  TC_INS/GNSS BEGINS ************************\n");
/*
  printf("P: %lf, %lf, %lf\n",pvat_old->latitude*R2D,pvat_old->longitude*R2D,pvat_old->height );
  printf("V: %lf, %lf, %lf\n",pvat_old->ned_velocity[0],pvat_old->ned_velocity[1],pvat_old->ned_velocity[2] );
  printf("A: %lf, %lf, %lf\n",pvat_old->euler_angles[0],pvat_old->euler_angles[1],pvat_old->euler_angles[2] );*/

  /* Initialize true navigation solution */
  time = INS_measurements->sec;
  true_L_b = pvat_old->latitude;
  true_lambda_b = pvat_old->longitude;
  true_h_b = pvat_old->height;
  for (i=0;i<3;i++) true_v_eb_n[i] = pvat_old->ned_velocity[i];
  for (i=0;i<3;i++) true_eul_nb[i] = pvat_old->euler_angles[i];
  Euler_to_CTM(true_eul_nb, true_C_b_n);

//   Euler_to_quaternion(true_eul_nb, q_nb);
//   Quaternion_to_DCM(q_nb, true_C_b_n);

  NED_to_ECEF(&true_L_b,&true_lambda_b,&true_h_b,true_v_eb_n,true_C_b_n,\
    old_true_r_eb_e,old_true_v_eb_e,old_true_C_b_e);

  /* Initialize other matrices */

  /* Determine Least-squares GNSS position solution */
  old_est_llh[0]=true_L_b; old_est_llh[1]=true_lambda_b; old_est_llh[2]=true_h_b;
  pos2ecef(old_est_llh, old_est_r_eb_e);
//  GNSS_LS_position_velocity(GNSS_measurements,no_GNSS_meas,\
  //  GNSS_config->init_est_r_ea_e,[0;0;0],old_est_r_eb_e,old_est_v_eb_e,est_clock);

  for (j=0;j<2;j++) est_clock[j]=clock_offset_drift[j];

  for (j=0;j<3;j++) old_est_v_eb_e[j]=old_true_v_eb_e[j];
  pv_ECEF_to_NED(old_est_r_eb_e,old_est_v_eb_e, &old_est_L_b,\
    &old_est_lambda_b,&old_est_h_b,old_est_v_eb_n);
  est_L_b = old_est_L_b;

  /* Initialize estimated attitude solution */
  //Initialize_NED_attitude(true_C_b_n, initialization_errors, old_est_C_b_n);
  for (j=0;j<9;j++) old_est_C_b_n[j]=true_C_b_n[j];

  NED_to_ECEF(&old_est_L_b,&old_est_lambda_b,&old_est_h_b,old_est_v_eb_n,\
    old_est_C_b_n,old_true_r_eb_e,old_true_v_eb_e,old_true_C_b_e);

  /* Initialize Kalman filter P matrix and IMU bias states */
  P_matrix=zeros(no_par,no_par); // in this case 17x17
  P_matrix_new=zeros(no_par,no_par);
  Initialize_TC_P_matrix(TC_KF_config, P_matrix); //IT WILL CHANGE FOR THE UNDIFF/UNCOMB MODEL
/*
  for (i = 0; i < 17; i++) {
    for (j = 0; j < 17; j++) {
        printf("P_matrix_bef_KF[%d]:%lf\t",i*17+j, P_matrix[i*17+j]);
    }
    printf("\n");
  }
  printf("\n");
*/
  if (out_errors[0] > 0.0 && out_errors[1]>0.0 && out_errors[2]>0.0) {
    for (i = 0; i < 17; i++) {
      for (j = 0; j < 17; j++) {
        if (i==j) {
          P_matrix[i*17+j]=out_errors[j];
        }
      }
    }
  }
  /*
  for (i = 0; i < 17; i++) {
    for (j = 0; j < 17; j++) {
        printf("P_matrix_bef_KF[%d]:%lf\t",i*17+j, P_matrix[i*17+j]);
    }
    printf("\n");
  }
  printf("\n");*/


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

  if (old_time<=0.000) {
    tor_i=1.0; //1 second of time interval
    time_last_GNSS=time-tor_i;
  }

    /* Correct IMU errors */
    for (i=0;i<3;i++) meas_f_ib_b[i]=INS_measurements->f_ib_b[i];
    for (i=0;i<3;i++) meas_omega_ib_b[i]=INS_measurements->omega_ib_b[i];
    for (i=0;i<3;i++) meas_f_ib_b[i] = meas_f_ib_b[i] - est_IMU_bias[i];
    for (i=0;i<3;i++) meas_omega_ib_b[i] = meas_omega_ib_b[i] - est_IMU_bias[i+3];

    /* Update estimated navigation solution */
printf("\n *****************  NAV_EQUATIONS BEGINS ************************\n");
    Nav_equations_ECEF(tor_i,\
        old_est_r_eb_e, old_est_v_eb_e,old_true_C_b_e, meas_f_ib_b,\
        meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
printf("\n *****************  NAV_EQUATIONS ENDS **************************\n");

    /* Determine whether to update GNSS simulation and run Kalman filter */
    printf("GNSS-INS DIFF: %lf\n", GNSS_measurements->sec - INS_measurements->sec);
  if ( fabs(GNSS_measurements->sec - INS_measurements->sec) < 0.0001) {
      /* KF time interval */
      tor_s = time - time_last_GNSS;

      /* Run Integration Kalman filter */
printf("\n *******************  TC_KF_EPOCH BEGINS ************************\n");
    /**/  TC_KF_Epoch(GNSS_measurements,no_GNSS_meas,tor_s,est_C_b_e,\
            est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,\
            meas_f_ib_b,est_L_b,TC_KF_config, est_C_b_e_new,\
            est_v_eb_e_new, est_r_eb_e_new, est_IMU_bias_new,\
            est_clock_new, P_matrix_new);
printf("\n *******************  TC_KF_EPOCH ENDS **************************\n");
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
*/
    Nav_equations_ECEF(tor_i,\
        old_est_r_eb_e, old_est_v_eb_e,old_true_C_b_e, meas_f_ib_b,\
        meas_omega_ib_b, est_r_eb_e_new, est_v_eb_e_new, est_C_b_e_new);

printf("\n *****************  NAV_EQUATIONS CLOSED-LOOP ENDS ****************\n");

      time_last_GNSS = GNSS_measurements->sec;

      /* Generate IMU bias and clock output records */
      fprintf(out_clock_file,"%lf %lf %lf\n", time, est_clock_new[0],\
      est_clock_new[1]);
      fprintf(out_IMU_bias_file,"%lf %lf %lf %lf %lf %lf %lf\n", time, \
      est_IMU_bias_new[0],est_IMU_bias_new[1],est_IMU_bias_new[2],\
      est_IMU_bias_new[3],est_IMU_bias_new[4], est_IMU_bias_new[5] );

      for (i = 0; i < 6; i++) out_IMU_bias_est[i]=est_IMU_bias_new[i];

      /* */printf("CLOCK: %lf %lf %lf\n", time, est_clock_new[0],\
        est_clock_new[1]);
        printf("IMU_bias: %lf %lf %lf %lf %lf %lf %lf\n", time, \
        est_IMU_bias_new[0],est_IMU_bias_new[1],est_IMU_bias_new[2],\
        est_IMU_bias_new[3],est_IMU_bias_new[4], est_IMU_bias_new[5] );

      /* Generate KF uncertainty output record */
      fprintf(out_KF_SD_file,"%lf\t", time);
      for (i=0;i<no_par;i++) {
        for (j=0;j<no_par;j++) {
          if (i==j) {
            out_errors[j]=P_matrix_new[i*no_par+j];
            fprintf(out_KF_SD_file,"%lf\t", sqrt(P_matrix_new[i*no_par+j]) );
          //  printf("KF_SD[%d]: %lf\n",i,sqrt(P_matrix_new[i*no_par+j]));
          }
        }
      }
      fprintf(out_KF_SD_file,"\n");

    }else{
      for (i=0;i<3;i++) est_r_eb_e_new[i]= est_r_eb_e[i];
      for (i=0;i<3;i++) est_v_eb_e_new[i]= est_v_eb_e[i];
      for (i=0;i<9;i++) est_C_b_e_new[i] = est_C_b_e[i];
      for (i=0;i<2;i++) est_clock_new[i] = est_clock[i];

    }

    /* Convert navigation solution to NED  */
    ECEF_to_NED(est_r_eb_e_new,est_v_eb_e_new,est_C_b_e_new,\
      &est_L_b,&est_lambda_b,&est_h_b,\
      est_v_eb_n,est_C_b_n);

   for (i=0;i<3;i++){
     for (j=0;j<3;j++) {
       est_C_b_n_T[i*3+j]=est_C_b_n[j*3+i];
     }
   }


  double q_b_n[4];
   DCM_to_quaternion(est_C_b_n, q_b_n);

      for (i=0;i<3;i++){
        for (j=0;j<3;j++) {
          printf("[%d]:%lf\t",j*3+i, est_C_b_e_new[j*3+i] );
        }
        printf("\n");
      }

   Quaternion_to_euler(q_b_n, pvat_new->euler_angles);

   printf("Euler from Quaternions: %lf, %lf, %lf\n",pvat_new->euler_angles[0],\
    pvat_new->euler_angles[1],pvat_new->euler_angles[2] );

    pvat_new->latitude = est_L_b;
    pvat_new->longitude = est_lambda_b;
    pvat_new->height = est_h_b;
    for (i=0;i<3;i++) pvat_new->ned_velocity[i] = est_v_eb_n[i];
    CTM_to_Euler(pvat_new->euler_angles, est_C_b_n);
    pvat_new->time=time;

    printf("Euler from CTM: %lf, %lf, %lf\n",pvat_new->euler_angles[0],\
     pvat_new->euler_angles[1],pvat_new->euler_angles[2] );

    /* Clock update */
    for(i=0;i<2;i++) out_clock[i]=est_clock_new[i];


  /* Free memory */
  free(P_matrix); free(P_matrix_new);
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
extern void LC_INS_GNSS_core(double* GNSS_r_eb_e, double* gnss_xyz_ini_cov,
     double* GNSS_v_eb_e, double* gnss_xyz_vel_cov, \
     double ini_pos_time, um7pack_t *imu_data, double imu_time_diff,\
     pva_t *PVA_old, imuraw_t *imuobsp, int Tact_or_Low_IMU){
  int no_par=15, no_meas=6;
  double out_errors[15]={0}, out_IMU_bias_est[6]={0}, out_KF_SD[15];
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
  float old_time=0.0, last_GNSS_time=0.0;
  int i, j, k;

  printf("\n *****************  LC INSGNSS CORE BEGINS ***********************\n");


  memset(&INS_meas, 0, sizeof(INS_measurements));

  printf("IMU_data time: %f\n", imu_data->sec);

  /* Current INS time */
  INS_meas.sec=imu_data->sec;
  INS_meas.tt = imu_time_diff;

  /* Last GNSS or State time and Last imu time */
  if (PVA_old->sec<=0.00) {
    old_time=pvat_old.time=PVA_old->sec=INS_meas.sec-0.5; /* Initialize according to IMU rate*/
    INS_meas.tt=0.5;
    old_time = (float) INS_meas.sec - INS_meas.tt;
    last_GNSS_time = ini_pos_time -1.0; /* Initialize according to GNSS rate*/
  }else{
    old_time=pvat_old.time=PVA_old->sec;
    last_GNSS_time = PVA_old->t_s; /* last state time */
  }

  /* Initial position and velocity from GNSS or previous PVA solution */
  if (imu_time_diff<0.1) {        /* If first solution, use GNSS */
    ecef2pos(GNSS_r_eb_e,llh);
    pvat_old.latitude = llh[0];
    pvat_old.longitude = llh[1];
    pvat_old.height = llh[2];
    ecef2enu(llh,GNSS_v_eb_e,pvat_old.ned_velocity);
     /* From gyrocompassing and levelling from current IMU meas. */
     for(j=0;j<3;j++) pvat_old.euler_angles[j]=imu_data->aea[j];
     for(j=0;j<3;j++) pvat_old.euler_angles[j]=0.0;
  }else{                 /* Otherwise, use previous PVA solution */
    ecef2pos(PVA_old->r,llh);
    pvat_old.latitude = llh[0];
    pvat_old.longitude = llh[1];
    pvat_old.height = llh[2];
    for(i=0;i<3;i++) pvat_old.ned_velocity[i]=PVA_old->v[i]; /*already in NED */
    for(i=0;i<3;i++) pvat_old.euler_angles[i]=PVA_old->A[i];

  }

  /* Constraining height with GNSS ones since there are no jumps */
  ecef2pos(GNSS_r_eb_e,llh);
  pvat_old.height = llh[2];


  /* Current IMU measurement structure initialization */
  for(j=0;j<3;j++) INS_meas.omega_ib_b[j] = imu_data->g[j];
  for(j=0;j<3;j++) INS_meas.f_ib_b[j] = imu_data->a[j];

  /* IMU bias */
  for(i=0;i<6;i++) out_IMU_bias_est[i]=PVA_old->out_IMU_bias_est[i];
  for(i=0;i<17;i++) out_errors[i]=PVA_old->out_errors[i];


  /* Load GNSS and INS errors and noises from configuration file -------------*/
  if (Tact_or_Low_IMU) {
    /* Tactical-grade imu */
    Initialize_INS_GNSS_LCKF_tactical_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,\
    &LC_KF_config);
  }else{
    /*Low-grade imu */
    Initialize_INS_GNSS_LCKF_consumer_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,\
      &LC_KF_config);
  }

  /* Position covariance from GNSS */
  for (j = 0; j < 3; j++) LC_KF_config.init_pos_unc[j]=sqrt(gnss_xyz_ini_cov[j]);
  for (j = 0; j < 3; j++) LC_KF_config.init_vel_unc[j]=sqrt(gnss_xyz_vel_cov[j]);

  if (out_IMU_bias_est[0]>0.0 || out_IMU_bias_est[3]>0.0) {
    for(i=0;i<3;i++) IMU_errors.b_a[i]=out_IMU_bias_est[i];
    for(i=0;i<3;i++) IMU_errors.b_g[i]=out_IMU_bias_est[i+3];
  }
  if (out_errors[0]>0.0 || out_errors[1]>0.0) {
    for(i=0;i<3;i++) initialization_errors.delta_eul_nb_n[i]=out_errors[i];
  }


  /* Loosely coupled ECEF Inertial navigation and GNSS integrated navigation -*/
   Loosely_coupled_INS_GNSS(&INS_meas, no_meas, &pvat_old, no_par, old_time, \
     last_GNSS_time, &initialization_errors, &IMU_errors,\
      &GNSS_config, &LC_KF_config, &pvat_new, out_errors, out_IMU_bias_est,\
      out_KF_SD, meas_f_ib_b, meas_omega_ib_b);

      llh[0]=pvat_new.latitude;
      llh[1]=pvat_new.longitude;
      llh[2]=pvat_new.height;
      pos2ecef(llh, PVA_old->r);
      for(i=0;i<3;i++) PVA_old->v[i]=pvat_new.ned_velocity[i];
      for(i=0;i<3;i++) PVA_old->A[i]=pvat_new.euler_angles[i];
      PVA_old->sec=pvat_new.time;

      for(i=0;i<6;i++) PVA_old->out_IMU_bias_est[i]=out_IMU_bias_est[i];
      for(i=0;i<17;i++) PVA_old->out_errors[i]=out_errors[i];

      PVA_old->t_s=last_GNSS_time;


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
extern void TC_INS_GNSS_core(rtk_t *rtk, const obsd_t *obs, int n,\
  const nav_t *nav, um7pack_t *imu_data, double imu_time_diff, pva_t *PVA_old, \
  int Tact_or_Low_IMU){

  int no_par=17, no_GNSS_meas;
  double out_errors[17]={0}, out_IMU_bias_est[6]={0}, out_clock[2], clock_offset[2], out_KF_SD[17];
  double meas_f_ib_b[3], meas_omega_ib_b[3];
  IMU_errors IMU_errors = {0};
  initialization_errors initialization_errors = {0};
  GNSS_config GNSS_config = {{0}};
  TC_KF_config TC_KF_Epoch = {{0}};
  PVAT_solution pvat_old = {{0}};
  PVAT_solution pvat_new = {{0}};
  GNSS_measurements GNSS_measurements[n];
  INS_measurements INS_meas = {0};
  double *rs, *dts, *var, *azel, llh[3], enu_ini_vel[3];
  double r,rr[3],e[3], dion, vion, dtrp, vtrp, lam_L1, P, vs[3];
  double old_time=0.0, last_GNSS_time=0.0;
  int svh[MAXOBS], i, j, k, sys, flag[n], m;
  prcopt_t opt=rtk->opt;

  printf("\n *****************  INSGNSS CORE BEGINS ***********************\n");

  /* Checking input values
  printf("POS: %lf, %lf %lf\n",rtk->sol.rr[0],rtk->sol.rr[1], rtk->sol.rr[2] );
  printf("OBS: %lf, %lf %lf\n",obs[0].P[0],obs[8].L[0], obs[10].D[0] );
  printf("Number of obs: %d\n",n);
  printf("IMUDATA: ax:%lf, ay:%lf, az:%lf and gx:%lf, gy:%lf, gz:%lf\n",\
  imu_data->a[0], imu_data->a[1],\
  imu_data->a[2], imu_data->g[0],imu_data->g[1], imu_data->g[2] );*/

/*  printf("P: %lf, %lf, %lf\n",PVA_old->r[0],PVA_old->r[1],PVA_old->r[2] );
  printf("V: %lf, %lf, %lf\n",PVA_old->v[0],PVA_old->v[1],PVA_old->v[2] );
  printf("A: %lf, %lf, %lf\n",PVA_old->A[0],PVA_old->A[1],PVA_old->A[2] );
  printf("dtr, dtrs: %lf, %lf\n", PVA_old->clock_offset_drift[0],PVA_old->clock_offset_drift[1]);*/

  /* Prepare GNSS and INS raw data into the proper structures --------------- */

  rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);

  memset(&INS_meas, 0, sizeof(INS_measurements));
  memset(&GNSS_measurements, 0, sizeof(GNSS_measurements));
  /* current GNSS time */
  GNSS_measurements->sec= time2gpst(rtk->sol.time,NULL); /* time of week in (GPST) */
  GNSS_measurements->tt=rtk->tt;
  /* Current INS time */
  INS_meas.sec=imu_data->sec;
  INS_meas.tt=imu_time_diff;

  /* Last GNSS or State time and Last imu time */
  if (PVA_old->sec<=0.00) {
    old_time=pvat_old.time=PVA_old->sec=INS_meas.sec-0.5; /* Initialize according to IMU rate*/
    //INS_meas.tt=0.5;
    old_time = (float) INS_meas.sec - INS_meas.tt;
    last_GNSS_time = GNSS_measurements->sec -1.0; /* Initialize according to GNSS rate*/
  }else{
    old_time=pvat_old.time=PVA_old->sec;
    last_GNSS_time = PVA_old->t_s;
  }

  printf("OLD TIME: %f, %f, %f\n",old_time, pvat_old.time, PVA_old->sec);

  /* Initial position and velocity from GNSS or previous PVA solution */
  if (rtk->tt<0.1) {        /* If first solution, use GNSS */
    ecef2pos(rtk->sol.rr,llh);
    pvat_old.latitude = llh[0];
    pvat_old.longitude = llh[1];
    pvat_old.height = llh[2];
    ecef2enu(llh,rtk->sol.rr+3,enu_ini_vel); /* Here is ENU! */
    /* ENU to NED velocity from gnss solution */
     pvat_old.ned_velocity[0]= enu_ini_vel[1];
     pvat_old.ned_velocity[1]= enu_ini_vel[0];
     pvat_old.ned_velocity[2]=-enu_ini_vel[2];
     /* From gyrocompassing and levelling from current IMU meas. */
     for(j=0;j<3;j++) pvat_old.euler_angles[j]=imu_data->aea[j];
     for(j=0;j<3;j++) pvat_old.euler_angles[j]=0.0;
     /* initialize clock offset and drift from gnss */
     clock_offset[0] = rtk->sol.dtr[0]*CLIGHT;//rtk->sol.dtr[0];//]x[3]; //or = rtk->sol.dtr[0];
     clock_offset[1] = rtk->sol.dtrr;
  }else{                 /* Otherwise, use previous PVA solution */
    ecef2pos(PVA_old->r,llh);
    pvat_old.latitude = llh[0];
    pvat_old.longitude = llh[1];
    pvat_old.height = llh[2];
    for(i=0;i<3;i++) pvat_old.ned_velocity[i]=PVA_old->v[i]; /*already in NED */
    for(i=0;i<3;i++) pvat_old.euler_angles[i]=PVA_old->A[i];

    /* initialize clock offset and drift from previous PVAT */
    clock_offset[0] = PVA_old->clock_offset_drift[0];
    clock_offset[1] = PVA_old->clock_offset_drift[1];
  }

  /* Constraining height with GNSS ones since there are no jumps */
  ecef2pos(rtk->sol.rr,llh);
  pvat_old.height = llh[2];


  /* initialize clock offset and drift from gnss */
  clock_offset[0] = rtk->sol.dtr[0]*CLIGHT;//rtk->sol.dtr[0];//]x[3]; //or = rtk->sol.dtr[0];
  clock_offset[1] = rtk->sol.dtrr;

  printf("Clock OFF and DRIFT: %lf, %lf\n", clock_offset[0], clock_offset[1]);

  /* From gyrocompassing and levelling from current IMU meas.
  for(j=0;j<3;j++) pvat_old.euler_angles[j]=imu_data->aea[j];*/

  /* Current IMU measurement structure initialization */
  for(j=0;j<3;j++) INS_meas.omega_ib_b[j] = imu_data->g[j];
  for(j=0;j<3;j++) INS_meas.f_ib_b[j] = imu_data->a[j];

  /* IMU bias */
  for(i=0;i<6;i++) out_IMU_bias_est[i]=PVA_old->out_IMU_bias_est[i];
  for(i=0;i<17;i++) out_errors[i]=PVA_old->out_errors[i];

  /* Receiver position */
  for (i=0;i<3;i++) rr[i]=rtk->sol.rr[i];
  //dtr=rtk->sol.dtr[0]; //(sec).
  ecef2pos(rr,llh);

  /* satellite positions and clocks */
  /* rs[(0:2)+i*6] {x,y,z} is the satellite position in ECEF
     rs[(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
     dts[(0:1)+i*2] are the sat clock {bias,drift} (s|s/s)*/
  satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
  m=0;

  for(i=0;i<n&&i<MAXOBS;i++){

    azel[i*2]=azel[1+i*2]=0.0;
    lam_L1=nav->lam[obs[i].sat-1][0];

    if (!(sys=satsys(obs[i].sat,NULL))) {flag[i]=0;m++;continue;}

    /* reject duplicated observation data */
   if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
    printf("duplicated observation data %s sat=%2d\n",\
          time_str(obs[i].time,3),obs[i].sat);
       i++;
       flag[i]=0;
       m++;
       continue;
    }

    /* geometric distance/azimuth/elevation angle */
    //printf("geodist: %lf, %lf, %lf\n", geodist(rs+i*6,rr,e), satazel(llh,e,azel+i*2), opt.elmin );
    if ((r=geodist(rs+i*6,rr,e))<=0.0||\
            satazel(llh,e,azel+i*2)<opt.elmin){printf("ERROR GEOM.\n");
            flag[i]=0;m++;continue;}

    /* psudorange with code bias correction   -> IT IS FOR P3
    if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;*/

    /* excluded satellite? */
      if (satexclude(obs[i].sat,svh[i],&opt)) {flag[i]=0; m++; continue;}

    /* ionospheric corrections */
      if (!ionocorr(obs[i].time,nav,obs[i].sat,llh,azel+i*2,\
                    IONOOPT_BRDC,&dion,&vion)) {printf("ERROR IONO!\n");
                     flag[i]=0;m++;continue;} //opt.ionoopt

      /* GPS-L1 -> L1/B1 */
      if ((lam_L1=nav->lam[obs[i].sat-1][0])>0.0) {
          dion*=(lam_L1/lam_carr[0])*(lam_L1/lam_carr[0]);
      }

      /* tropospheric corrections */
      if (!tropcorr(obs[i].time,nav,llh,azel+i*2,\
                    TROPOPT_SAAS,&dtrp,&vtrp)) {//opt.tropopt
          printf("ERROR TROPO!\n");
          flag[i]=0;
          m++; continue;
      }

      if (obs[i].P[0] <= 0.0 || obs[i].D[0] <= 0.0 || obs[i].L[0] <= 0.0) {
          flag[i]=0;
          m++;
          continue;
      }

      /* pseudorange residual
      //printf("Pres: %lf\n", obs[i].P[0]-(r+CLIGHT*rtk->sol.dtr[0]-CLIGHT*dts[i*2]+dion+dtrp));
      //printf("Dres: %lf\n", -obs[i].D[j]*lam_L1+CLIGHT*dts[1+i*2]);
      printf("RANGE: %lf, ION: %lf, Tropo: %lf, SAT clock: %lf, REC clk: %lf \
      REC.POS: %lf, %lf, %lf, Azel: %lf, PVA_oldpos: %lf,  %lf, %lf\n ", r, dion, dtrp, CLIGHT*dts[i*2], \
      CLIGHT*rtk->sol.dtr[0], llh[0], llh[1],llh[2],azel[i*2], PVA_old->r[0],PVA_old->r[1],PVA_old->r[2]);
*/
    flag[i]=1;
    GNSS_measurements[i].sat = obs[i].sat;
    GNSS_measurements[i].time = obs[i].time;
    for(j=0;j<2;j++) GNSS_measurements[i].P[j] = obs[i].P[j]-dion-dtrp+CLIGHT*dts[i*2];//-clock_offset[0];
    for(j=0;j<2;j++) GNSS_measurements[i].L[j] = obs[i].L[j];
    for(j=0;j<2;j++) GNSS_measurements[i].D[j] = -obs[i].D[j]*lam_L1+CLIGHT*dts[1+i*2];//-clock_offset[1]; /* hz to m/s (radial velocity)*/
    for(j=0;j<3;j++) GNSS_measurements[i].Sat_r_eb_e[j]=rs[i*6+j];
    for(j=0;j<3;j++) GNSS_measurements[i].Sat_v_eb_e[j]=rs[i*6+(j+3)];
  }

  /* Number of valid satellites */
  no_GNSS_meas = n-m;

k=0;
for (i = 0; i < n; i++) {
  if (flag[i]) {
    /* Valid satellite values */
    GNSS_measurements[i-k].sat = GNSS_measurements[i].sat;
    GNSS_measurements[i-k].time = GNSS_measurements[i].time;
    for(j=0;j<2;j++) GNSS_measurements[i-k].P[j] = GNSS_measurements[i].P[j];
    for(j=0;j<2;j++) GNSS_measurements[i-k].L[j] = GNSS_measurements[i].L[j];
    for(j=0;j<2;j++) GNSS_measurements[i-k].D[j] = GNSS_measurements[i].D[j]; /* hz to m/s (radial velocity)*/
    for(j=0;j<3;j++) GNSS_measurements[i-k].Sat_r_eb_e[j]=GNSS_measurements[i].Sat_r_eb_e[j];
    for(j=0;j<3;j++) GNSS_measurements[i-k].Sat_v_eb_e[j]=GNSS_measurements[i].Sat_v_eb_e[j];
  }else{
    /*Skip position - invalid satelite*/
    k++;
  }
}

/* Put zeros on the invalid positions?? */
if (m>0) {
  for (i=1;i<=m;i++) {
    GNSS_measurements[n-i].sat = 0;
    for(j=0;j<2;j++) GNSS_measurements[n-i].P[j] = 0.0;
    for(j=0;j<2;j++) GNSS_measurements[n-i].L[j] = 0.0;
    for(j=0;j<2;j++) GNSS_measurements[n-i].D[j] = 0.0; /* hz to m/s (radial velocity)*/
    for(j=0;j<3;j++) GNSS_measurements[n-i].Sat_r_eb_e[j]=0.0;
    for(j=0;j<3;j++) GNSS_measurements[n-i].Sat_v_eb_e[j]=0.0;
  }
}

/* Load GNSS and INS errors and noises from configuration file -------------*/

if (Tact_or_Low_IMU) {
  /* Tactical-grade imu */
  Initialize_INS_GNSS_TCKF_tactical_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,\
  &TC_KF_Epoch);
}else{
  /*Low-grade imu */
  Initialize_INS_GNSS_TCKF_consumer_grade_IMU(&initialization_errors, &IMU_errors, &GNSS_config,\
  &TC_KF_Epoch);
}

  /* Position and velocity covariance from GNSS */
  /*
  if (rtk->sol.qrv[1]<0.0) {
    rtk->sol.qrv[1]=-1.0*rtk->sol.qrv[1];
  }

  for (j = 0; j < 3; j++) TC_KF_Epoch.init_pos_unc[j]=sqrt(rtk->sol.qr[j]);
  for (j = 0; j < 3; j++) TC_KF_Epoch.init_vel_unc[j]=sqrt(rtk->sol.qrv[j]);
  TC_KF_Epoch.init_clock_offset_unc=sqrt(rtk->sol.qdtr); // It is in s should be in m ?
  TC_KF_Epoch.init_clock_drift_unc=sqrt(rtk->sol.qdtrr); // In (s/s) it should be in m/s?
  //for (i = 0; i < 3; i++) =rtk->sol.qr[i];
*/
  if (out_IMU_bias_est[0]>0.0 || out_IMU_bias_est[3]>0.0) {
    for(i=0;i<3;i++) IMU_errors.b_a[i]=out_IMU_bias_est[i];
    for(i=0;i<3;i++) IMU_errors.b_g[i]=out_IMU_bias_est[i+3];
  }
  if (out_errors[0]>0.0 || out_errors[1]>0.0) {
    for(i=0;i<3;i++) initialization_errors.delta_eul_nb_n[i]=out_errors[i];
  }

  /* Tightly coupled ECEF Inertial navigation and GNSS integrated navigation -*/
   Tightly_coupled_INS_GNSS(&INS_meas, &GNSS_measurements, no_GNSS_meas,\
      &pvat_old, no_par, old_time, last_GNSS_time, clock_offset, &initialization_errors, &IMU_errors,\
      &GNSS_config, &TC_KF_Epoch, &pvat_new, out_errors, out_IMU_bias_est,\
      out_clock, out_KF_SD, meas_f_ib_b, meas_omega_ib_b);

      llh[0]=pvat_new.latitude;
      llh[1]=pvat_new.longitude;
      llh[2]=pvat_new.height;
      pos2ecef(llh, PVA_old->r);
      for(i=0;i<3;i++) PVA_old->v[i]=pvat_new.ned_velocity[i];
      for(i=0;i<3;i++) PVA_old->A[i]=pvat_new.euler_angles[i];
      PVA_old->sec=pvat_new.time;

      for(i=0;i<2;i++) PVA_old->clock_offset_drift[i]=out_clock[i];
      for(i=0;i<6;i++) PVA_old->out_IMU_bias_est[i]=out_IMU_bias_est[i];
      for(i=0;i<17;i++) PVA_old->out_errors[i]=out_errors[i];

      PVA_old->t_s=last_GNSS_time;

  /* Plots -------------------------------------------------------------------*/

  /* Write output profile and errors file ------------------------------------*/

  /* Free memory -------------------------------------------------------------*/
  free(rs); free(dts); free(var);
  free(azel);
  printf("\n *****************  INSGNSS CORE ENDS *************************\n");
}

extern void imu_tactical_navigation(FILE *imu_file) {
 char str[100];
 float t_prev=0.0,t_curr, tor_i;
 double old_r_eb_e[3], old_v_eb_n [3], old_v_eb_e [3], old_C_b_e[9];
 double f_ib_b[3], omega_ib_b[3], f_ib_n[3], C_b_n[9], eul_nb_n[3];
 double r_eb_e[3], v_eb_e[3], C_b_e[9], eul_eb_e[3];
 double est_r_eb_e[3], est_v_eb_e[3], est_C_b_e[3];
 double wiee[3], llh[3], gan[3], q[4];
 um7pack_t imu={0};
 double G = 9.80665;
 int i,j;
 imuraw_t imu_obs_prev={0};
 um7pack_t imu_curr_meas={0};
 pva_t PVA_prev_sol={{0}};

 /* Update with previous solution */
 imu_obs_prev=imu_obs_global;
 PVA_prev_sol=pva_global;

 /* Initialize position, velocity and attitude */
 old_r_eb_e[0]= 1761300.949;
 old_r_eb_e[1]= -4078202.553;
 old_r_eb_e[2]= 4561403.792;
 /* Local or navigation-frame velocities and acceleration */
 old_v_eb_n[0]=0.0;//-0.02;
 f_ib_n[0]= -0.144; //E-W
 old_v_eb_n[1]=0.0;//-0.01;
 f_ib_n[1]= -0.2455; //N-S
 old_v_eb_n[2]=0.0;// 0.04;
 f_ib_n[2]=  0.605; //U-D

 /* Attitude Initialization (Groves, 2013)  */

 /* Earth rotation vector in e-frame	*/
 wiee[0]=0;wiee[1]=0;wiee[2]=OMGE;

 ecef2pos(old_r_eb_e, llh);
 /* local apparent gravity vector */
 appgrav(llh, gan, wiee);

 /* Coarse alignment or use of previous solution? */
 /* Levelling and gyrocompassing, update imu->aea[] vector */
 for (i=0;i<3;i++) imu.a[i]=f_ib_n[i];
 for (i=0;i<3;i++) imu.v[i]=old_v_eb_n[i];
 coarseAlign(&imu, gan);

 for (i=0;i<3;i++) eul_nb_n[i]=imu.aea[i];

 Euler_to_quaternion(eul_nb_n, q);
 Quaternion_to_DCM(q, C_b_n); //Euler_to_CTM(eul_nb_n, C_b_n);

 NED_to_ECEF(&llh[0], &llh[1], &llh[2], old_v_eb_n, C_b_n, \
 old_r_eb_e, old_v_eb_e, old_C_b_e);
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

 while ( fgets(str, 100, imu_file)!= NULL ){

   sscanf(str, "%f %lf %lf %lf %lf %lf %lf", &t_curr, &f_ib_b[2],\
   &f_ib_b[0],&f_ib_b[1], &omega_ib_b[2],&omega_ib_b[0],&omega_ib_b[1]);

   tor_i=t_curr-t_prev;

   /* Turn-on bias */
   /* Turn-on initial biases             //exp4*/
  f_ib_b[0]=(f_ib_b[0])*G; //-0.0339421
  f_ib_b[1]=(f_ib_b[1])*G; //+0.105076
  f_ib_b[2]=(f_ib_b[2])*G; //+1.00335
  omega_ib_b[0]=(omega_ib_b[0])*D2R;
  omega_ib_b[1]=(omega_ib_b[1])*D2R;
  omega_ib_b[2]=(omega_ib_b[2])*D2R;
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


   /* Stationary-condition */
  if (norm(old_v_eb_n,3) > 0.5) {
    //  printf("IS MOVING\n");
  }else{//printf("IS STATIC\n");
    //for (j=0;j<3;j++) PVA_prev_sol.v[j]=0.0;
  }

   printf("\n *****************  NAV_EQUATIONS BEGINS ************************\n");
      Nav_equations_ECEF(tor_i,\
          old_r_eb_e, old_v_eb_e, old_C_b_e, f_ib_b,\
          omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
   printf("\n *****************  NAV_EQUATIONS ENDS **************************\n");

   for (i=0;i<3;i++) old_r_eb_e[i]=est_r_eb_e[i];
   for (i=0;i<3;i++) old_v_eb_e[i]=est_v_eb_e[i];
   for (i=0;i<9;i++) old_C_b_e[i]= est_C_b_e[i];
   t_prev=t_curr;

   ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, &llh[0], &llh[1], &llh[2],\
   old_v_eb_n, C_b_n);
/*
   printf("C_b_n:%lf, %lf, %lf, |%lf, %lf, %lf, |%lf, %lf, %lf\n",\
    C_b_n[0],C_b_n[1],C_b_n[2],\
    C_b_n[3],C_b_n[4],C_b_n[5],\
    C_b_n[6],C_b_n[7],C_b_n[8]);
*/
   ecef2pos(est_r_eb_e, llh);

   /* Cbn to euler */
   CTM_to_Euler(eul_nb_n, C_b_n);
   printf("eul_nb_n: %lf, %lf, %lf\n",eul_nb_n[0],eul_nb_n[1],eul_nb_n[2] );
   DCM_to_quaternion(C_b_n,q);

   //Quaternion_to_euler(q, eul_nb_n);
   printf("eul_nb_n: %lf, %lf, %lf\n",eul_nb_n[0],eul_nb_n[1],eul_nb_n[2] );


   fprintf(out_PVA,"%f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",\
   t_curr, llh[0]*R2D, llh[1]*R2D, llh[2],\
   old_v_eb_e[0], old_v_eb_e[1], old_v_eb_e[2],
   eul_nb_n[0]*R2D,eul_nb_n[1]*R2D,eul_nb_n[2]*R2D);

 }



}
