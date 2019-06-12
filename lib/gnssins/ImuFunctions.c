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

/* Utility ----------------------------------------------------------*/
extern inline int8_t sign(double val)
{
  /* Corresponding sign() function in Matlab*/
  if (val < 0.0)
    return -1;
  if (val > 0.0)
    return 1;
  return 0;
}

/* Quaternion multiplication ---------------------------------------------------
* description: Multiplies p . q  quaternions
* args   : double *p      I    first quaternion {4x1}
*	         double *q		  I    second quaternion {4x1}
*	         double *pq     IO   resulting quaternion {4x1}
* Groves (2103) from Appendix E - at E.6.3 section, apge E-10
*-----------------------------------------------------------------------------*/
extern void quaternion_mult(double *p, double *q, double *pq)
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
extern void attitude_update(double *alpha_ib_b, double *omega_ie, float t, double *old_C_b_e, double *new_q_b_e)
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
extern void DCM_to_quaternion(double *C, double *q)
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
extern void Quaternion_to_DCM(double *q, double *C)
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
extern void Euler_to_quaternion(double *eul_nb, double *q_nb)
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
extern void Quaternion_to_euler(double *q, double *eul)
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
extern void Cbn_from_rotations(um7pack_t *imu, double *C)
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
extern void Quaternion_attitude_errror_correction(double *est_delta_Psi,
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