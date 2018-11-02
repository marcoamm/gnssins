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
void Tightly_coupled_INS_GNSS(in_profile,no_epochs,initialization_errors,
  IMU_errors,GNSS_config,TC_KF_config,
  //OUTPUT TOGETHER
out_profile,out_errors,out_IMU_bias_est,out_clock,out_KF_SD);

void TC_KF (){

}






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
int InsGnssCore (void)
{
  int no_epochs;
  /* Prepare GNSS and INS raw data into the proper structures */

  /* Load GNSS and INS errors and noises from configuration file */

  /* Tightly coupled ECEF Inertial navigation and GNSS integrated navigation */
[out_profile,out_errors,out_IMU_bias_est,out_clock,out_KF_SD] =...
    Tightly_coupled_INS_GNSS(in_profile,no_epochs,initialization_errors...
    ,IMU_errors,GNSS_config,TC_KF_config);

  /* Plots */

  /* Write output profile and errors file */



    return(0);
}
