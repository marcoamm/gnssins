#include <stdio.h>
#include <stdint.h>
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

extern void Initialize_INS_GNSS_LCKF_consumer_grade_IMU(initialization_errors *initialization_errors,
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
extern void Initialize_INS_GNSS_TCKF_consumer_grade_IMU(initialization_errors *initialization_errors,
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
extern void Initialize_INS_GNSS_LCKF_tactical_grade_IMU(initialization_errors *initialization_errors,
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
extern void Initialize_INS_GNSS_TCKF_tactical_grade_IMU(initialization_errors *initialization_errors,
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