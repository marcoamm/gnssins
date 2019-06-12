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
extern void Initialize_NED_attitude(double *C_b_n,
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
extern void Initialize_TC_P_matrix(TC_KF_config *TC_KF_config, double *P_matrix)
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
extern void Initialize_LC_P_matrix(LC_KF_config *LC_KF_config, double *P_matrix)
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

/* initialize state and covariance -------------------------------------------*/
extern void tcinitx(ins_states_t *ins,double xi, double var, int i)
{
    int j;
    ins->x[i]=xi;
    for (j=0;j<ins->nx;j++) {
        ins->P[i+j*ins->nx]=ins->P[j+i*ins->nx]=i==j?var:0.0;
    }
} 
