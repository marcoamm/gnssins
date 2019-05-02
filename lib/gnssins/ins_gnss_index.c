/*-----------------------------------------------------------------------------
* ins-gnss-state.cc : ins-gnss coupled estimated states interface
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/02/06 1.0 new
*----------------------------------------------------------------------------*/
#include <rtklib.h>

/* get number of attitude states---------------------------------------------*/
extern int xnA() {return 3;}
/* get number of velocity states---------------------------------------------*/
extern int xnV() {return 3;}
/* get number of position states---------------------------------------------*/
extern int xnP() {return 3;}
/* get number of accl/gyro bias states---------------------------------------*/
extern int xnBa() {return 3;}
extern int xnBg() {return 3;}

/* get number of receiver clock bias (non-close-loop correction states)------*/
extern int xnRc(const prcopt_t* opt)
{
    return ((opt)->navsys<5?1:2);
}
/* get number of receiver clock drift (non-close-loop correction states) ----*/
extern int xnRr() {return 1;}

/* get number of tropo and phase bias are non-close-loop correction states
 * --------------------------------------------------------------------------*/
extern int xnT(const prcopt_t* opt)
{
    return ( (opt)->tropopt<TROPOPT_EST?0:1 );
}
extern int xnB(){return MAXSAT;}

extern int xnRx(const prcopt_t* opt)
{
    return (xnP ()+xnV ()+xnA ()+xnBa()+xnBg()+xnRc(opt)+xnRr()+xnT(opt));
}
/* get number of all states--------------------------------------------------*/
extern int xnX(const prcopt_t* opt) {return xnRx(opt)+xnB();}

extern int xiA () {return 0;}    /* attitude states index */
extern int xiV () {return 3;}    /* velocity states index */
extern int xiP () {return 6;}    /* position states index */
extern int xiBa() {return 9;}    /* accl bias state index */
extern int xiBg() {return 12;}    /* accl bias state index */

/* get index of receiver clock bias state------------------------------------*/
extern int xiRc()
{
    return xnA ()+xnV ()+xnP ()+xnBa()+xnBg();
}
/* get index of receiver clock drift-----------------------------------------*/
extern int xiRr(const prcopt_t* opt)
{
    return xiRc()+xnRc(opt);
}
/* get index of tropos (r:0=rov,1:ref)---------------------------------------*/
extern int xiTr(const prcopt_t* opt)
{
    return  xiRr(opt)+xnRr();
}
/* get index of phase bias (s:satno,f:freq)----------------------------------*/
extern int xiBs(const prcopt_t* opt,int s)
{
    return xiTr(opt)+s;
}