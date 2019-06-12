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
 
  /* From Groves Matlab code */
  /* Calculate Euler angles using 2.23 */

  euler[0] = atan2(Cbn[5], Cbn[8]); /* Roll (x rotation) - phi */
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
