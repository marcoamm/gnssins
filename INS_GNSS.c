#include <stdio.h>
#include <stdlib.h>

#include "rtklib.h"
#include "INS_GNSS.h"



/* struct declaration -------------------------------------------------------*/



/* function declaration ------------------------------------------------------*/


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
 void Euler_to_CTM(double *delta_eul_nb_n, double *delta_C_b_n){


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
void Initialize_NED_attitude(double *true_C_b_n, \
  initialization_errors *initialization_errors, double *old_est_C_b_n){


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
void Nav_equations_ECEF(float *tor_i,\
    double *old_est_r_eb_e,double *old_est_v_eb_e, double *old_est_C_b_e, \
    double *meas_f_ib_b, double *meas_omega_ib_b, double *est_r_eb_e, \
    double *est_v_eb_e, double *est_C_b_e){



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
*-----------------------------------------------------------------------------
 void TC_KF_Epoch(GNSS_measurements, no_GNSS_meas, tor_s, est_C_b_e,
      est_v_eb_e, est_r_eb_e, est_IMU_bias, est_clock, P_matrix,
      meas_f_ib_b,est_L_b,TC_KF_config, est_C_b_e, est_v_eb_e, est_r_eb_e,
      est_IMU_bias, est_clock, P_matrix){

}*/

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
 void ECEF_to_NED(double *est_r_eb_e, double *est_v_eb_e, double *est_C_b_e, \
   double *est_L_b, double *est_lambda_b, double *est_h_b, double *est_v_eb_n, \
   double *est_C_b_n){



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
void Tightly_coupled_INS_GNSS(int no_epochs, \
  initialization_errors *initialization_errors, IMU_errors *IMU_errors,\
  GNSS_config *GNSS_config, TC_KF_config *TC_KF_config, double *out_profile,\
  double *out_errors, double *out_IMU_bias_est, double *out_clock,\
  double *out_KF_SD)
{
  

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
int main (void){
//int InsGnssCore (){ LATER, WHEN LINKED TO OTHER MAIN FUNCTION USE IT THIS WAY
  int no_epochs;
  double *out_profile, *out_errors, *out_IMU_bias_est, *out_clock, *out_KF_SD;
  IMU_errors IMU_errors = {0};
  initialization_errors initialization_errors = {0};
  GNSS_config GNSS_config = {0};
  TC_KF_config TC_KF_Epoch = {0};


  /* Prepare GNSS and INS raw data into the proper structures */

  /* Load GNSS and INS errors and noises from configuration file */

  /* Tightly coupled ECEF Inertial navigation and GNSS integrated navigation */
    Tightly_coupled_INS_GNSS(no_epochs, &initialization_errors, &IMU_errors, \
      &GNSS_config, &TC_KF_Epoch, out_profile, out_errors, out_IMU_bias_est, \
      out_clock, out_KF_SD);

  /* Plots */

  /* Write output profile and errors file */
  printf("FINISHED InsGNSScore !\n");
  return 0;
}
