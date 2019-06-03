#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "rtklib.h"
#include "satinsmap.h"

/* global variables ----------------------------------------------------------*/

//FILE *fimu;          /* Imu datafile pointer  */


/* Local apparent gravity ----------------------------------------------------
* description: computes the gravity vector in the navigation frame (n-frame)
* args	:	double *pos[3]	I	{llh} vector [m] in e-frame
*	 	double *gan[3]	O	apparent gravity vector in n-frame [m/s*s]
* Reference: Chatfiled (1997, pag. 80)
*----------------------------------------------------------------------------*/
void appgrav (double* pos, double* gan, double* wiee){
 double ge[3], r[3];
 double auxvec[3], aux1vec[3], gae[3], Ren[9];
 int i;

 /* llh to xyz	*/
 pos2ecef(pos, r);

 /* apparent gravity vector in e-frame
 from normal gravity vector at the location and altitude	*/
 ge[0]=0;ge[1]=0;ge[2]=gn(pos[0], pos[2]);

 cross3(wiee, r, auxvec);
 cross3(wiee, auxvec, aux1vec);

 for (i=0; i<3;i++){gae[i]=ge[i]-aux1vec[i];}

 /* Rotation matrix from e to n-frame	*/
 Ce2n(pos[0],pos[1],Ren);

 /* Local apparent gravity vector in n-frame	*/
 matmul_row("NN", 3, 1, 3, 1.0, Ren, gae, 0.0, gan);

 for (i=0; i<3;i++) gan[i]=-gan[i];

 /* Forcing gan back to ge, as in Shin (2001) page 24 */
 for (i=0; i<3;i++) gan[i]=ge[i];


}

/* Project  the Earth vector to the n-frame ----------------------------------------------------
* description: computes the Earthrotation vector in n-frame
* args	:	double *pos[3]	I	{llh} vector [deg,m] in e-frame
*	 	double *wien[3]	O	rotation rate vector in n-frame
* Reference: Shin (2005, pag.13)
*----------------------------------------------------------------------------*/
void projvec(double* pos, double* wien){

 wien[0]=OMGE*cos(pos[0]);
 wien[1]=0;
 wien[2]=-OMGE*sin(pos[0]);

}

/* Change lat long rates from n wrt e-frame ----------------------------------------------------
* description: computes the Earth rotation vector in n-frame w.r.t e-frame
* args	:	double *pos[3]	I	{llh} position vector [deg,m] in e-frame
		double *v[3]	I	{vn,ve,vd} velocity vector [m] in n-frame
*	 	double *wenn[3]	O	rotation rate vector in n-frame wrt e-frame
* Reference: Shin (2001, pag.13)
*----------------------------------------------------------------------------*/
void raten(double* pos, double* v, double* wenn){

 wenn[0]= v[1]/(RN(pos[0])+pos[2]);
 wenn[1]=-v[0]/(RM(pos[0])+pos[2]);
 wenn[2]=-v[1]*tan(pos[0])/(RN(pos[0])+pos[2]);

}

/* winn ----------------------------------------------------
* description: computes winn
* args	:	double *r[3]	I	{lat,lon,h} geodesic position vector [rad,m] in e-frame
		double *v[3]	I	{ve,vn,vu} velocity vector [m] in n-frame
*	 	double *we[3]	O	earth rotation rate vector
*	 	double *winn[3]	O	rotation rate vector in n-frame wrt e-frame
* Reference: Shin (2005, pag. )
*----------------------------------------------------------------------------*/
void winncomp(double* r, double* v, double* winn){

 winn[0]= OMGE*cos(r[0]) +(v[0]/(RN(r[0])+r[2]));
 winn[1]=-v[1]/(RM(r[0])+r[2]);
 winn[2]=-OMGE*cos(r[0])-(v[0]*tan(r[0])/(RN(r[0])+r[2]));

}

/* Coarse alignment -----------------------------------------------------------
* description: Coarse alignment of roll, pitch and yaw
* args   : um7pack_t* imu       IO    imu structure
*	   double* gn	I	local gravity
           double acc[3]	I	INS accelerations from previous epoch {xyz}
*          double v[3]	I	GNSS/or INS/GNSS previous state velocities {NED}
           gan[3]  I  gravity vector
           euler_angles[3]  IO  attitude angles in _nb frame {roll, pitch, yaw}
Reference: Shin (2005, pag.42) and Groves (2013, pages 198 to 199)
*-----------------------------------------------------------------------------*/
//void coarseAlign(um7pack_t* imu, double* gan)
void coarseAlign(double *acc, double *gyr, double *vel, double *gan, double *euler_angles)
{
 double r1[9], r2[9], aux[3], sign=1.0;
 double sineyaw, coseyaw;

 /* Is it a static or kinematic alignment? Use GNSS velocity or any velocity in imu->v */
 if ( norm(vel,3) < 0.5){
   printf("LEVELLING.STATIC, %lf, %lf, %lf\n", vel[0], vel[1],vel[2]);
                                    //static alignment
  /* Then we can solve for the roll(x(phi)) and pitch(y(theta)) angles using accelerometers */
  //if(acc[2]<0.0){sign=-1.0;}
  //accea[0]= sign*(asin((acc[1])/norm(gan,3)));
  //accea[1]=-sign*(asin((acc[0])/norm(gan,3)));

  /* Levelling from Groves (2013) */
  /* roll - phi angle */
  euler_angles[0]= atan2(-acc[1],-acc[2]);

  /* pitch - theta angle */
  euler_angles[1]= atan (acc[0]/sqrt(acc[1]*acc[1]+acc[2]*acc[2]));

  /* Gyrocompasing by Groves(2013,page 199)  */
  /* yaw - psi angle */
  sineyaw=-gyr[1]*cos(euler_angles[0])+gyr[2]*sin(euler_angles[0]);
  coseyaw= gyr[0]*cos(euler_angles[1])+gyr[1]*sin(euler_angles[0])*sin(euler_angles[1])+\
  gyr[2]*cos(euler_angles[0])*sin(euler_angles[1]);
  euler_angles[2]=atan2(sineyaw,coseyaw);

}else{                            //kinematic alignment
 /* Roll can be initialized to zero with an uncertainty of +-5 degrees,
 in most cases, on the road */
 printf("LEVELLING.KINE, %lf, %lf, %lf\n", vel[0], vel[1],vel[2] );

 euler_angles[0]=0.0;

 /* Levelling from Groves (2013) */
 /* roll - phi angle */
 //euler_angles[0]= atan2(-acc[1],-acc[2]); /* CORRECT IT LATER*/

 /* The pitch(y(theta)) and heading(z(psi)) angles using velocities */
 euler_angles[1]= atan(vel[2]/sqrt(vel[0]*vel[0]+vel[1]*vel[1]));
 euler_angles[2]= atan(vel[2]/vel[1]);
 /* Can be improved using an EKF with LHU model, see Shin (2005) chapter 3	*/
 }

}

/* Determine imu errors --------------------------------------------------------
* description: determine the IMU calibration parameters (biases and scale factors)
* args   : um7pack_t* imu       IO    imu structure
*            double* l          I accel. observation vector {nx3} [x,y,z]
*            double* lg          I gyro. observation vector {nx3} [x,y,z]
*            int n              I number of observations
*	           double* gn	        I     normal gravity vector
*	           double* wiee	      I     earth rotation vector in e-frame
Reference: Shin (2001, pag. 55)
*-----------------------------------------------------------------------------*/
void imucalib(um7pack_t* imu, double* l, double* lg, int n, double* gn, double* wiee)
{
 int u=6,ug=3; /* Number of parameters for acc. (bias+scale factors) and
 gyro (g) (only bias) models*/
 int n0=(int)u*3,n0g=(int)ug*3; /* minimium # of observation for acc. and gyro*/
 int r=(int)n/3; /* Number of equations - 1 for each 3 set of IMU obs {ax,ay,az}
 or {gx,gy,gz} */
 double *A,*B,*w,*x,*v;  /*  acceleration quantities */
 double *Ag,*Bg,*wg,*xg,*vg; /*  gyroscopes quantities */
 double *M,*Mg,*N,*Ng,*u_,*ug_;
 double *v1,*v2,*v3,*vg1,*vg2,*vg3;
 double *la,*lga;
 double bw=0, bg=0; 	/* gyro and acc. scale biases	*/
 double sgx=0, sgy=0, sgz=0; 	/* acc. scale factors	*/
 double *xcpy,*witr;
 int i,j,k,p=0;

 xcpy=mat(u,1);
 la=mat(n,1);lga=mat(n,1);witr=mat(n,1);

 /* Correcting parameters: x=x0+xc */
for (i = 0; i < u; i++) xcpy[i]=imu->Ba[i];
for (i = 0; i < r; i++) {
  la[i*3]=l[i*3];
  la[i*3+1]=l[i*3+1];
  la[i*3+2]=l[i*3+2];
}

while (p != -1){
 A=mat(r,u);B=zeros(r,n);w=mat(n,1);x=mat(u,1);v=mat(n,1);M=zeros(r,r);N=zeros(u,u);
 u_=mat(u,1);
 Ag=mat(r,ug);Bg=zeros(r,n);wg=mat(n,1);xg=mat(ug,1);vg=mat(n,1);Mg=zeros(r,r);
 Ng=zeros(ug,ug);ug_=mat(ug,1);
 v1=mat(n,r);v2=mat(r,1);vg1=mat(n,r);vg2=mat(r,1);v3=mat(r,1);vg3=mat(r,1);

 /* Design matrices A and B for acceleremoters  */
 for (i = 0; i < r; i++) {
   /*dfg/d(biases {bgx,bgy,bgz})  */
   A[i*u]=-2*(l[i*3]-imu->Ba[0])/((1+imu->Ba[3])*(1+imu->Ba[3]));
   A[i*u+1]=-2*(l[i*3+1]-imu->Ba[1])/((1+imu->Ba[4])*(1+imu->Ba[4]));
   A[i*u+2]=-2*(l[i*3+2]-imu->Ba[2])/((1+imu->Ba[5])*(1+imu->Ba[5]));
   /*dfg/d(scale factors {sgx,sgy,sgz})  */
   A[i*u+3]=-2*(l[i*3]-imu->Ba[0])*(l[i*3]-imu->Ba[0])/((1+imu->Ba[3])*(1+imu->Ba[3])*(1+imu->Ba[3]));
   A[i*u+4]=-2*(l[i*3+1]-imu->Ba[1])*(l[i*3+1]-imu->Ba[1])/((1+imu->Ba[4])*(1+imu->Ba[4])*(1+imu->Ba[4]));
   A[i*u+5]=-2*(l[i*3+2]-imu->Ba[2])*(l[i*3+2]-imu->Ba[2])/((1+imu->Ba[5])*(1+imu->Ba[5])*(1+imu->Ba[5]));
 }

 for (i = 0; i < r; i++) {
   for (j = 0; j < n; j++) {
    // printf("%lf\t",A[i*n+j] );
   }
  // printf("\n");
 }

 for (i = 0; i < r; i++) {
   /*dfg/l(observables {lax,lay,laz})  */
     B[i*n+(i*3)]=2*(l[i*3]-imu->Ba[0])/((1+imu->Ba[3])*(1+imu->Ba[3]));
     B[i*n+(i*3)+1]=2*(l[i*3+1]-imu->Ba[1])/((1+imu->Ba[4])*(1+imu->Ba[4]));
     B[i*n+(i*3)+2]=2*(l[i*3+2]-imu->Ba[2])/((1+imu->Ba[5])*(1+imu->Ba[5]));
  }

  /* Misclosure vector for acceleration model  */
  if (p>0) {
    /* w = B(l-l0i)+F(x0i,l0i)code */
    for (i=0;i<r;i++){
      for (j=0;j<n;j++) {
        witr[i]+=B[i*n+j]*(l[i*3+j]-la[i*3+j]); /* Vector B*(l-l0i) */
      }
    }

    for (i = 0; i < r; i++) {
      w[i]=witr[i]+((la[i*3]-imu->Ba[0])/(1+imu->Ba[3]))*((la[i*3]-imu->Ba[0])/(1+imu->Ba[3]))+
      ((la[i*3+1]-imu->Ba[1])/(1+imu->Ba[4]))*((la[i*3+1]-imu->Ba[1])/(1+imu->Ba[4]))+
      ((la[i*3+2]-imu->Ba[2])/(1+imu->Ba[5]))*((la[i*3+2]-imu->Ba[2])/(1+imu->Ba[5]))-
      (sqrt(gn[0]*gn[0]+gn[1]*gn[1]+gn[2]*gn[2]))*(sqrt(gn[0]*gn[0]+gn[1]*gn[1]+gn[2]*gn[2]));
      //printf("w[%d]=%lf\n",i,w[i]);
    }
  }else{
    for (i = 0; i < r; i++) {
      w[i]=((l[i*3]-imu->Ba[0])/(1+imu->Ba[3]))*((l[i*3]-imu->Ba[0])/(1+imu->Ba[3]))+
      ((l[i*3+1]-imu->Ba[1])/(1+imu->Ba[4]))*((l[i*3+1]-imu->Ba[1])/(1+imu->Ba[4]))+
      ((l[i*3+2]-imu->Ba[2])/(1+imu->Ba[5]))*((l[i*3+2]-imu->Ba[2])/(1+imu->Ba[5]))-
      (sqrt(gn[0]*gn[0]+gn[1]*gn[1]+gn[2]*gn[2]))*(sqrt(gn[0]*gn[0]+gn[1]*gn[1]+gn[2]*gn[2]));
    //  printf("w[%d]=%lf\n",i,w[i]);
    }
 }

 /* M^-1 = (BB^t)^-1  */
 for (i = 0; i < r; i++) {
   for (j = 0; j < r; j++) {
     M[i*r+j]=(i==j?(1/(B[i*n+(i*3)]*B[i*n+(i*3)]+B[i*n+(i*3)+1]*B[i*n+(i*3)+1]+
     B[i*n+(i*3)+2]*B[i*n+(i*3)+2])):0.0);
   }
 }
 /* N */
 for (i = 0; i < u; i++) {
   for (j = 0; j < u; j++) {
     for (k = 0; k < r; k++) {
       N[i*u+j]+=A[k*u+i]*M[k*r+k]*A[k*u+j];
     }
   }
 }
 /* u */
   for (k = 0; k < u; k++) {
     for (i = 0; i < r; i++) {
       u_[k]+=A[i*u+k]*M[i*r+i]*w[i];
     }
   }

 /* Misclosure vector for gyroscopes model  */
 for (i = 0; i < r; i++) {
   wg[i]=((lg[i*3]-imu->Bg[0])*(lg[i*3]-imu->Bg[0])) +
        ((lg[i*3+1]-imu->Bg[1])*(lg[i*3+1]-imu->Bg[1]))+
        ((lg[i*3+2]-imu->Bg[2])*(lg[i*3+2]-imu->Bg[2]))-
   (sqrt(wiee[0]*wiee[0]+wiee[1]*wiee[1]+wiee[2]*wiee[2]))*(sqrt(wiee[0]*wiee[0]+wiee[1]*wiee[1]+wiee[2]*wiee[2]));
//printf("w[%d]=%lf\n",i,wg[i]);
 }

 /* Design matrices A and B for gyroscopes  */
 for (i = 0; i < r; i++) {
   /*dfw/d(biases {bwx,bwy,bwz})  */
   Ag[i*ug]=-2*(lg[i*3]-imu->Bg[0]);
   Ag[i*ug+1]=-2*(lg[i*3+1]-imu->Bg[0+1]);
   Ag[i*ug+2]=-2*(lg[i*3+2]-imu->Ba[0+2]);
   /*dfw/l(observables {lwx,lwy,lwz})  */
   Bg[i*n+(i*3)]=2*(lg[i*3]-imu->Bg[0]);
   Bg[i*n+(i*3)+1]=2*(lg[i*3+1]-imu->Bg[0+1]);
   Bg[i*n+(i*3)+2]=2*(lg[i*3+2]-imu->Bg[0+2]);
 }

 for (i = 0; i < r; i++) {
   for (j = 0; j < n; j++) {
     printf("%lf\t",Ag[i*n+j] );
   }
   printf("\n");
 }

 for (i = 0; i < r; i++) {
   for (j = 0; j < ug; j++) {
  //   printf("%lf\t",Ag[i*ug+j] );
   }
  // printf("\n");
 }

 /* Mb^-1 = (BB^t)^-1  */
 for (i = 0; i < r; i++) {
   for (j = 0; j < r; j++) {
     Mg[i*r+j]=(i==j?(1/(Bg[i*n+(i*3)]*Bg[i*n+(i*3)]+Bg[i*n+(i*3)+1]*Bg[i*n+(i*3)+1]
     +Bg[i*n+(i*3)+2]*Bg[i*n+(i*3)+2])):0.0);
   }
 }

 /* Ng */
 for (i = 0; i < ug; i++) {
   for (j = 0; j < ug; j++) {
     for (k = 0; k < r; k++) {
       Ng[i*ug+j]+=Ag[k*ug+i]*Mg[k*r+k]*Ag[k*ug+j];
     }
   }
 }
 /* ug_*/
   for (k = 0; k < ug; k++) {
     for (i = 0; i < r; i++) {
       ug_[k]+=Ag[i*ug+k]*Mg[i*r+i]*wg[i];
     }
   }

 /* least square estimation parameter correction xc ----*/
 matmul("NN",u,1,u,1.0,N,u_,0.0,x);
 matmul("NN",ug,1,ug,1.0,Ng,ug_,0.0,xg);

 /* Residuals vectors v=-P^-1.B^t.M^-1.(AX+W)*/
 matmul("TN",n,r,r,1.0,B,M,0.0,v1);
 matmul("TN",n,r,r,1.0,Bg,Mg,0.0,vg1);
 matmul("NN",r,1,u,1.0,A,x,0.0,v2);
 matmul("NN",r,1,ug,1.0,Ag,xg,0.0,vg2);
 for (i = 0; i < r; i++) v3[i]=v2[i]+w[i];
 for (i = 0; i < r; i++) vg3[i]=vg2[i]+wg[i];
 matmul("NN",n,1,r,1.0,v1,v3,0.0,v);
 matmul("NN",n,1,r,1.0,vg1,vg3,0.0,vg);

 for (i = 0; i < r; i++) {
   la[i*3]=l[i*3]-v[i*3];
   la[i*3+1]=l[i*3+1]-v[i*3+1];
   la[i*3+2]=l[i*3+2]-v[i*3+2];
   lga[i*3]=lg[i*3]-vg[i*3];
   lga[i*3+1]=lg[i*3+1]-vg[i*3+1];
   lga[i*3+2]=lg[i*3+2]-vg[i*3+2];
 }

 /* Correcting parameters: x=x0+xc */
for (i = 0; i < u; i++) xcpy[i]=imu->Ba[i]; /* For iteration control  */
for (i = 0; i < u; i++) imu->Ba[i]=imu->Ba[i]+x[i];
for (i = 0; i < ug; i++) imu->Bg[i]=imu->Bg[i]+xg[i];


free(A);free(B);free(w);free(x);free(v);free(M);free(N);free(u_);
free(Ag);free(Bg);free(wg);free(xg);free(vg);free(Mg);free(Ng);free(ug_);
free(v1);free(v2);free(v3);
free(vg1);free(vg2);free(vg3);
p++;

if (abs(xcpy[0]-imu->Ba[0]) && abs(xcpy[1]-imu->Ba[1]) &&
    abs(xcpy[2]-imu->Ba[2]) && abs(xcpy[3]-imu->Ba[3]) &&
    abs(xcpy[4]-imu->Ba[4]) && abs(xcpy[5]-imu->Ba[5]) <= 0.1 ) {
  p=-1;
 }
 //printf("%lf,%lf,%lf,%lf,%lf,%lf\n",abs(xcpy[0]-imu->Ba[0]),abs(xcpy[1]-imu->Ba[1]),
//abs(xcpy[2]-imu->Ba[2]),abs(xcpy[3]-imu->Ba[3]),abs(xcpy[4]-imu->Ba[4]),
//abs(xcpy[5]-imu->Ba[5]));

} /* while loop end  */

//printf("Acc par. \n");
 for (i = 0; i < u; i++) {
//     printf("%lf,",imu->Ba[i]);
 }
 //printf("Gyro par. \n");
  for (i = 0; i < u; i++) {
  //    printf("%lf,",imu->Bg[i]);
  }
//puts("\n");

free(la);free(lga);
free(xcpy);free(witr); 

}

/* Remove biases from IMU measurements -----------------------------------------
* description: account for gyro and accelerometers biases and scale factors
* args   : um7pack_t* imu       IO    imu structure
*	   double dt	        I     time increment of the time interval
Reference: Shin (2001, pag. 21)
*-----------------------------------------------------------------------------*/
void imuerrorcorr(um7pack_t* imu, double dt)
{
 int i;
        /* INSERT THE CORRECT VALUES! */
 double bw=0, bg=0; 	/* gyro and acc. scale biases	*/
 double sgx=0, sgy=0, sgz=0; 	/* acc. scale factors	*/
 double D[9], aux[3];

//printf("BEf.IMU acc. and gyro: %lf,%lf,%lf and %lf,%lf,%lf\n", imu->a[0],imu->a[1],imu->a[2],imu->g[0],imu->g[1],imu->g[2] );

 /* Corrected gyroscopes outputs	*/
 for (i=0; i<3;i++) imu->g[i]=imu->g[i]-imu->Bg[i]*dt;

 D[1]=D[2]=D[3]=D[5]=D[6]=D[7]=0;
 D[0]=1/(1+sgx);
 D[4]=1/(1+sgy);
 D[8]=1/(1+sgz);

 for (i=0; i<3;i++) aux[i]=imu->a[i]-imu->Ba[i]*dt;

 /* Corrected acceleremoters outputs	*/
 matmul("NN", 3, 1, 3, 1.0, D, aux, 0.0, imu->a);

 //printf("Aft.IMU acc. and gyro: %lf,%lf,%lf and %lf,%lf,%lf\n", imu->a[0],imu->a[1],imu->a[2],imu->g[0],imu->g[1],imu->g[2] );

}

/* Attitude integration -------------------------------------------------------
* description: Update the DCM Cbn from gyro corrected measurements
* args   : um7pack_t* imu       IO    current (tk) imu structure
*	   double dt		I     time difference from previous measu. (tk-1)
*	   double* wien[3], wenn[3]	I	rotation vectors
*	   double* C[9]		IO    receives an initial Cbn and outputs the updated one
Reference: Shin (2001, pag.22)
*-----------------------------------------------------------------------------*/
void attint(um7pack_t* imu, double dt, double* wien, double* wenn, double* C)
{
 double gnbb[3],gmag,waux[3],waux1[3],s,c,q[4],qp[4],Q[16], qaux[4];
 int i,j;

 /* Body angular increment wrt the n-frame  */
 for (i=0; i<3;i++) waux[i]=(wien[i]+wenn[i])*dt;
 matmul("NN",3, 1, 3, 1.0, C, waux, 0.0, waux1);
 for (i=0; i<3;i++) gnbb[i]=(imu->g[i]-waux1[i]);
 /* magnitude of the angular increment */
 gmag = sqrt(gnbb[0]*gnbb[0]+gnbb[1]*gnbb[1]+gnbb[2]*gnbb[2]);

 /* Updating gyro meas. from gibb to gnbb */
 for (i=0; i<3;i++) imu->g[i]=gnbb[i];

 /* Update the quaternion	*/
 s=1-(gmag*gmag/24)+(gmag*gmag*gmag*gmag/1920);
 c=-(gmag*gmag/4)+(gmag*gmag*gmag*gmag/192);

 Q[0]=Q[5]=Q[10]=Q[15]=c;
 Q[1]=Q[11]=s*gnbb[2]; /*wz*/
 Q[4]=Q[14]=-s*gnbb[2]; /*-wz*/
 Q[7]=Q[8]=s*gnbb[1]; /*wy*/
 Q[13]=Q[2]=-s*gnbb[1]; /*-wy*/
 Q[3]=Q[6]=s*gnbb[0]; /*wx*/
 Q[9]=Q[12]=-s*gnbb[0]; /*-wx*/

 //for (i=0; i<9;i++) C[i]=-C[i]; /* back to initial Cbn */
 dcm2quat(C, qp); /*  previous quaternion */

 matmul("NN", 4, 1, 4, 0.5, Q, qp, 0.0, qaux);

 for (i=0; i<4;i++) q[i]=qp[i]+qaux[i];

 /* The updated DCM Cbn is accomplished by */
 quat2dcm(C,q);
}

/* Attitude update from DCM to quaternion ---------------- ---------------------
* description: transform the DCM Cbn to quaternion
* args   : double* C[9]       I    DCM Cbn matrix
* 	   double* q[4]	      IO   quaternion
Reference: Shin (2001, pag.19)
*-----------------------------------------------------------------------------*/
void dcm2quat(double* C, double* q)
{
 /* Transformation between DCM Cbn and the quaternion is accomplished by: */
 q[0]=0.25*(C[7]-C[5])/(0.5*sqrt(1+C[0]+C[4]+C[8]));
 q[1]=0.25*(C[2]-C[6])/(0.5*sqrt(1+C[0]+C[4]+C[8]));
 q[2]=0.25*(C[3]-C[1])/(0.5*sqrt(1+C[0]+C[4]+C[8]));
 q[3]=0.5*sqrt(1+C[0]+C[4]+C[8]);
}

/* Quaternion to DCM ------------------------------------- ---------------------
* description: transform a quaternion to DCM Cbn amtrix
* args   : double* C[9]       IO    DCM Cbn matrix
* 	       double* q[4]	      I   quaternion
Reference: Shin (2001, pag.22)
*-----------------------------------------------------------------------------*/
void quat2dcm(double* C, double* q)
{
  /* The transformation between the quaternion and the DCM Cbn is accomplished by */
  C[0]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]-q[3]*q[3]; C[1]=2*(q[0]*q[1]-q[2]*q[3]); C[2]=2*(q[0]*q[2]-q[1]*q[3]);
  C[3]=2*(q[0]*q[1]-q[2]*q[3]); C[4]=q[1]*q[1]-q[0]*q[0]-q[2]*q[2]+q[3]*q[3]; C[5]=2*(q[1]*q[2]-q[0]*q[3]);
  C[6]=2*(q[0]*q[2]-q[1]*q[3]); C[7]=2*(q[1]*q[2]-q[0]*q[3]); C[8]=q[2]*q[2]-q[0]*q[0]-q[1]*q[1]+q[3]*q[3];
}

/* Convert quaternion to rotation matrix ----------------- ---------------------
* description: quaternion approach for updating the attitude
* args   : um7pack_t* imu       IO    imu structure
*
Reference: Infante (2016, pag.38)
*-----------------------------------------------------------------------------*/
//void quat2dcm(double* C, double* q){}

/* Attitude update with quaternion --------------------- ---------------------
* description: quaternion approach for updating the attitude matrix Cbn
* args   : um7pack_t* imu       IO    imu structure
*
Reference: Shin (2001, pag.18)
*-----------------------------------------------------------------------------*/
void quatconv(um7pack_t* imu, double* C)
{
 double q[4], u, qdoti[16], qdot[4];

 /* Building quaternion from a rotation angle vector u{ux,uy,uz} */

 /* QUESTION: Is it from the attitude Euler angles or  the rotation rate from gyroscopes measurements */
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

/* Acceleration (Specific force) transformation -------------------------------
* description: Transform acceleration (Specific force) to n-fram
* args   : um7pack_t* imu       IO    imu structure
*	   double* C[9]		I     Cbn rotation matrix
Reference: Shin (2001, pag.23)
*-----------------------------------------------------------------------------*/
void accb2n(um7pack_t* imu, double* C)
{
 double A[9], vecaux[3];
 int i,j;

 /* b-frame to n-frame applying first order sculling correction	*/
 A[0]=A[4]=A[8]=1;
 A[1]=0.5*imu->g[2];
 A[3]=-A[1];
 A[2]=0.5*imu->g[1];
 A[6]=-A[2];
 A[5]=0.5*imu->g[0];
 A[7]=-A[5];

 matmul("NN", 3, 1, 3, 1.0, A, imu->a, 0.0, vecaux);
 matmul("NN", 3, 1, 3, 1.0, C, vecaux, 0.0, imu->a);
}

/* Apply Coriolis and Gravity correction ----------------------------------------
* description: correct Coriolis and gravity on the velocity increments in the
* e-frame
* args   : um7pack_t* imu       IO    imu structure
*	   double* gvn[3]	I     normal gravity vector in n-frame
*	   double* wien[3]	I     Earth's rotation vector in n wrt i
*	   double* wenn[3]	I     Earth's rotation vector in n wrt e
*	   double dt		I     measurement time span (tk-tk-1)
Reference: Shin (2001, pag.24) and Groves(previous vers. of 2013 book, page 133)
*-----------------------------------------------------------------------------*/
void accgravcor(um7pack_t* imu, um7pack_t* imuprev, double* gvn, double* wien, double* wenn, double dt)
{
 int i;
 double vec[3], auxvec[3], v[3];

 /* Velocity increment is obtained by applying Coriolis and gravity correction:	*/
 for(i=0;i<3;i++) {vec[i]=2*wien[i]+wenn[i]; v[i]=imuprev->v[i]*dt;}
 cross3(vec, v, auxvec);

 /* With gravity correction */
 for(i=0;i<3;i++) imu->a[i]=imu->a[i]-auxvec[i]+gvn[i]*dt;

}

/* Velocity and Position Integration ------------------------------------------
* description: estimate velocity and position in n-frame
* args   : um7pack_t* imu       IO   current (tk) imu measurements in n-frame
*	   um7pack_t* imuprev	I    previous (tk-1) imu measurement in n-frame
*	   double dt	        I    time span between measurements
Reference: Shin (2001, pag.24)
*-----------------------------------------------------------------------------*/
void velposinteg(um7pack_t* imu, um7pack_t* imuprev, double dt){
 double D[9], vec[3], aux[3];
 int i;

 /* Assigning first measurements to the starting point */
 if(imuprev->sec == 0.00) imuprev->sec =imu->sec;

 printf("IMU-r: %lf, %lf, %lf\n",imu->r[0]*R2D,imu->r[1]*R2D,imu->r[2] );
 printf("IMU-v: %lf, %lf, %lf\n",imu->v[0],imu->v[1],imu->v[2]);
 printf("PREV-r: %lf, %lf, %lf\n",imuprev->r[0]*R2D,imuprev->r[1]*R2D,imuprev->r[2] );
 printf("PREV-v: %lf, %lf, %lf\n",imuprev->v[0],imuprev->v[1],imuprev->v[2]);

 /* The velocity integration can be performed: v(k) = v(k-1) + dv(k)	*/
 for (i=0;i<3;i++) imu->v[i] = imuprev->v[i] + imu->a[i];

 /* Positions are integrated using the second order Runge-Kutta method:	*/
 D[1]=D[2]=D[3]=D[5]=D[6]=D[7]=0.0;
 D[0]=1/(RM(imu->r[0])+imu->r[2]);
 D[4]=1/(RN(imu->r[0])+imu->r[2]);
 D[8]=-1;

 for (i=0;i<3;i++) vec[i]=(imuprev->v[i]+imu->v[i])*dt;
 matmul("NN", 3, 1, 3, 0.5, D, vec, 0.0, aux);
 printf("AUX: %lf, %lf, %lf\n",aux[0],aux[1],aux[2]);
/*
 imu->r[0]=(imuprev->r[0]+aux[0]);
 imu->r[1]=(imuprev->r[1]+aux[1]);
 imu->r[2]= imuprev->r[2]+aux[2];*/

/* From Groves (previous vers. of book) page 134 */
 imu->r[2]= imuprev->r[2]-dt/2*(imuprev->v[2]+imu->v[2]);
 imu->r[0]= imuprev->r[0]+dt/2*(imuprev->v[0]/(RN(imuprev->r[0])+imuprev->r[2])+
            imu->v[0]/(RN(imuprev->r[0])+imu->r[2]));
 imu->r[1]= imuprev->r[1]+dt/2*( imuprev->v[1]/((RM(imuprev->r[0])+imuprev->r[2])*cos(imuprev->r[0]))+
                       imu->v[1]/((RM(imu->r[0])+imu->r[2])*cos(imu->r[0])) );


}


/* IMU normal forces correction -------------------------------------------------
* description: treat the imu acceleration to compute velocity and position
* args   : float *ai       I   acceleration vector (xyz) in gravities and I-frame
* One epoch message type:
*-----------------------------------------------------------------------------*/
void imuacc(um7pack_t* imu){
 float g=9.8; // gravity velocity in m/s/s
 int i;

 /* Converting acceleration to m/s/s	*/
 for (i=0;i<3;i++) imu->a[i]=imu->a[i]*g;

 /* Removing normal forces (eq. 3 AN-1007)	*/
 imu->a[2]=imu->a[2]+g;

 return;
}

/* IMU velocity and position using accelerometer --------------------------------
* description: estimate velocity and positionn
* args   : um7pack_t* imu       I   imu structure in I-frame
* Document of reference: CH robotics AN-1007
*-----------------------------------------------------------------------------*/
void imuvelpos(um7pack_t* imu, um7pack_t* imuprev){
 float v[3], r[3];
 int i;

 if(imuprev->sec == 0) imuprev->sec =imu->sec;

 /* Estimating velocities: v(k) = v(k-1) + t*ai(k-1)	*/
 for (i=0;i<3;i++) imu->v[i] = imuprev->v[i]+ imuprev->a[i]*(imu->sec - imuprev->sec) + (imu->sec - imuprev->sec)*((imu->a[i] - imuprev->a[i]))/2.0;

 //for (i=0;i<3;i++) imu->v[i] = imuprev->v[i] + (imu->sec - imuprev->sec)*(imu->a[i]);

 /* Estimating positions from velocities: r(k) = r(k-1) + t*vi(k-1)	*/
 //for (i=0;i<3;i++) imu->r[i] = imuprev->r[i] + (imu->sec - imuprev->sec)*imu->v[i];

 for (i=0;i<3;i++) imu->r[i] = imuprev->r[i] + imuprev->v[i] + (imu->sec - imuprev->sec)*(imu->v[i]-imuprev->v[i]);
}

/* Accounting for initial biases ----------------------------------------------
* description: IMU
* args   : um7pack_t* imu       IO
*-----------------------------------------------------------------------------*/
void imucorr(um7pack_t* imu)
{
 float x_corr, y_corr, z_corr;
 float ax, bx, ay, by, az, bz;
 float corrx, corry, corrz;

 /* z-axis bias	*/
 /* Measuring the gravity biases in static periods	*/
 //z_corr = 1 - (0.883+0.884+0.885+0.882)/4; /* mean 3-axis resultant values */
 z_corr = -0.791173; // From matlab mean computation
 //(-0.799200-0.801300-0.795000-0.796600)/4; /* mean z acc on steady periods	*/
 /* should be zero when static	 */
 x_corr = 0.36628; // From matlab mean computation
 //(0.348500+0.346600+0.357600+0.36010)/4; /* mean x-acc on steady periods	*/
 y_corr = 0.125368; // From matlab mean computation
 //(0.127300+0.143600+0.122000+0.13140)/4; /* mean y-acc on steady periods	*/

 /* Mean value acceleration after detrending	*/
 x_corr = -0.000023;
 y_corr = 0.0000019585;
 z_corr = -1.000145;

 /* Detrending data with a fitted line equation (determined in Matlab app)	*/
 /* x,y,z acceleration meas. straight line coeficients	without bias correction */
 ax=0.00000062; bx=0.06218; ay=0.000004871; by=-2.264; az=0.000005335; bz=-3.408;
 imu->a[2]= imu->a[2]-(1); //gravity account

 /* x,y,z acceleration meas. straight line coeficients after bias correction
 ax=0.000006076; bx=-2.98; ay=0.00004774; by=-23.42; az=0.00005228; bz=-25.65;*/
 /* Correction value (f(x)), or value on the line	*/
 corrx=(ax*imu->sec)+bx; corry=(ay*imu->sec)+by; corrz=(az*imu->sec)+bz;
 //printf("Corrections x: %f, y: %f, z: %f\n", corrx, corry, corrz);

 /* Correcting measurements	*/
 imu->a[0]= imu->a[0]-corrx;
 imu->a[1]= imu->a[1]-corry;
 imu->a[2]= imu->a[2]-corrz;

 /* Correcting biases in the measurements 	*/
 imu->a[0]= imu->a[0]-x_corr;
 imu->a[1]= imu->a[1]-y_corr;
 imu->a[2]= imu->a[2]-(z_corr+1);

 return;
}

/* Compute Euler from Cbn ------------------------------------------------------
* description: Compute euler angles from the DCM Cbn
* args   : um7pack_t* imu       IO
* ref.: (Shin, 2001, page 12) and Groves(prev version of 2013, page 28)
*-----------------------------------------------------------------------------*/
void Cbn2euler(um7pack_t* imu, double *Cbn)
{
  /* Roll (x rotation) - phi */
  imu->aea[0]=atan2(Cbn[7],Cbn[8]);
  /* Pitch (y rotation) - theta */
//  imu->aea[1]=-1/(tan(Cbn[6]/sqrt(1-Cbn[6]*Cbn[6])));
  imu->aea[1]=-asin(Cbn[6]);
  /* Yaw (z rotation) - psi */
  imu->aea[2]=atan2(Cbn[3],Cbn[0]);
}

/* Parse and organize the MEMs-IMU received packets ----------------------------
* description: IMU
* args   : char *strline       I   char with one line containing one sensor IMU readings
* One epoch message type:
489629.200,1946,COM3,36,$PCHRS,0,101.633,3.32,0.33,-0.49,*7E
489629.200,1946,COM3,42,$PCHRS,1,101.633,0.3717,0.1030,-0.7725,*77
489629.200,1946,COM3,42,$PCHRS,2,101.633,0.4922,0.2708,-1.6283,*7D, where
secofweek sync, GPS week, port, ?, header filed, count, sensor time, x,y,z,checksum
    O  imu.status -1: meassage not complete, 0: error reading, 1: ok
*-----------------------------------------------------------------------------*/
void parseimudata(char* strline, um7pack_t* imu)
{
   char *token;
   int i=1, sscanstat;
   float sensor_x,sensor_y,sensor_z;
   double imucurrtime;
   double G = 9.80665;

  sscanstat=sscanf(strline, "%lf,%d,%*4s,%*2d,%*6s,%d,%f,%*f,%*f,%*f,%*s",
  &imu->sec, &imu->gpsw, &imu->count, &imu->internal_time);

  /* Checking if received message is complete and synchronized	*/
  if (sscanstat < 4) {
    imu->status=-1;
    return;
  }

  sscanstat=sscanf(strline, "%lf,%d,%*4s,%*2d,%*6s,%d,%*f,%f,%f,%f,%3c",
  &imu->sec, &imu->gpsw, &imu->count, &sensor_x, &sensor_y, &sensor_z, &imu->checksum);
  /*printf("%6.1f, %d, %d, %f, %f, %f, %c \n",
  imu->sec, imu->gpsw, imu->count, sensor_x, sensor_y, sensor_z, imu->checksum);*/
  imu->time=gpst2time(imu->gpsw,imu->sec);

 /* Check if message was read successfully	*/
 if (sscanstat < 7) {
   imu->status=0;
   return;
 }

 if (imu->count == 0){
  // printf("Gyro data: %d \n", imu->count);
   imu->g[0]=sensor_x*D2R; imu->g[1]=sensor_y*D2R;imu->g[2]=sensor_z*D2R; /* rad/s */
   imu->tgyr=imu->internal_time;
 }else if (imu->count == 1){
       	//  printf("Accelerometer data: %d \n", imu->count);
	  imu->a[0]=sensor_x*G;imu->a[1]=sensor_y*G;imu->a[2]=sensor_z*G; /* in gravities to m/s2 */
    imu->tacc=imu->internal_time;
       }else {
	 // printf("Magnetometer data: %d \n", imu->count);
	  imu->m[0]=sensor_x;imu->m[1]=sensor_y;imu->m[2]=sensor_z; /* uT unitless */
    imu->tmag=imu->internal_time;
       }

   /* Checking if data is complete for processing	*/
   /* Using time of each data measurement */
   if (imu->tgyr==0.0 || imu->tacc==0.0 ) {
     /* Not complete yet */
     imu->status=0;
     return;
   }

   if (imu->tgyr==imu->tacc) {
     /* Data complete and synchronized */
     imu->status=1;
   }else{
     /* Data not synchronized */
     imu->status=0;
   }
   //if ( imu->count==2 ) imu->status = 1;

}

/* INS-15 error-model system (transition) matrix, Groves (2013, page 587) */
extern void inssysmatrix(double *PHI, double *G, int nx, pva_t *pva_global,
  imuraw_t *imuobsc, double dt){
  double *F,*I;
  int i,j;
  double M,N,R,lat,h,vn,ve,vd,g;
  double tbax,tbay,tbaz,tbgx,tbgy,tbgz; /*acc. and gyro correlation times */
  double fb[3],fn[3],Cbn[9],winn[3];

  /* Values extracted from Angrisano (2010,page 93) for low-cost IMU Crista */
  tbax=tbay=tbaz=270; // [s]
  tbgx=tbgy=tbgz=350; // [s]

  lat=pva_global->r[0]; h=pva_global->r[2];
  ve=pva_global->v[0];vn=pva_global->v[1];vd=pva_global->v[2];
  M=RM(lat);
  N=RN(lat); /* Groves uses RN for M and RE for N*/
  R=sqrt(M*N);
  g=gn(lat,h);

  for (i=0;i<3;i++) fb[i]=imuobsc->fb[i];
//  Cb2n(pva_global.A[0],pva_global.A[1],pva_global.A[2],Cbn);
  for (i=0;i<9;i++) Cbn[i]=pva_global->Cbn[i];
  matmul("NN",3,1,3,1.0,Cbn,fb,0.0,fn);

  winn[0]= OMGE*cos(lat) +(ve/(RN(lat)+h));
  winn[1]=-vn/(RM(lat)+h);
  winn[2]=-OMGE*cos(lat)-(ve*tan(lat)/(RN(lat)+h));

  F=zeros(nx,nx);
  I=eye(nx);

  /* Building matrix F in navigation-frame */

    /* Attitude error dynamics Fee=(winn X),Fev, Fer, 0(3), Cbn */
    /* 1st line */
    F[nx*0+1]=winn[2];
    F[nx*0+2]=-winn[1];
    F[nx*0+4]=1/(N+h);
    F[nx*0+6]=-OMGE*sin(lat);
    F[nx*0+8]=-ve/((N+h)*(N+h));
    F[nx*0+12]=-Cbn[0];
    F[nx*0+13]=-Cbn[1];
    F[nx*0+14]=-Cbn[2];

    /* 2nd line */
    F[nx*1+0]=-winn[2];
    F[nx*1+2]=winn[0];
    F[nx*1+3]=-1/(M+h);
    F[nx*1+8]=vn/((M+h)*(M+h));
    F[nx*1+12]=-Cbn[3];
    F[nx*1+13]=-Cbn[4];
    F[nx*1+14]=-Cbn[5];

    /* 3rd line */
    F[nx*2+0]=winn[1];
    F[nx*2+1]=-winn[0];
    F[nx*2+4]=-tan(lat)/(N+h);
    F[nx*2+6]=-OMGE*cos(lat)-(ve/((N+h)*cos(lat)*cos(lat)));
    F[nx*2+8]=ve*tan(lat)/((N+h)*(N+h));
    F[nx*2+12]=-Cbn[6];
    F[nx*2+13]=-Cbn[7];
    F[nx*2+14]=-Cbn[8];

    /* Velocity error dynamics Fvr,Fvv, Fve=(f X), Cbn, 0(3)*/
    /* 1st line */
    F[nx*3+1]=-fn[2];
    F[nx*3+2]= fn[1];
    F[nx*3+3]=vd/(M+h);
    F[nx*3+4]=-2*OMGE*sin(lat)-2*(ve*tan(lat)/(N+h));
    F[nx*3+5]=vn/(M+h);
    F[nx*3+6]=-2*ve*OMGE*cos(lat)-((ve*ve)/((N+h)*cos(lat)*cos(lat)));
    F[nx*3+8]=-vn*vd/((M+h)*(M+h))+(ve*ve*tan(lat)/((N+h)*(N+h)));
    F[nx*3+9]=Cbn[0];
    F[nx*3+10]=Cbn[1];
    F[nx*3+11]=Cbn[2];

    /* 2nd line */
    F[nx*4+0]=fn[2];
    F[nx*4+2]=-fn[0];
    F[nx*4+3]=2*OMGE*sin(lat)+(ve*tan(lat)/(N+h));
    F[nx*4+4]=vd+vn*tan(lat)/(N+h);
    F[nx*4+5]=2*OMGE*cos(lat)+(ve/(N+h));
    F[nx*4+6]=2*OMGE*(vn*cos(lat)-vd*sin(lat))+(ve*vn/((N+h)*cos(lat)*cos(lat)));
    F[nx*4+8]=-ve*vd/((N+h)*(N+h))-(vn*ve*tan(lat)/((N+h)*(N+h)));
    F[nx*4+9]=Cbn[3];
    F[nx*4+10]=Cbn[4];
    F[nx*4+11]=Cbn[5];

    /* 3rd line */
    F[nx*5+0]=-fn[1];
    F[nx*5+1]=fn[0];
    F[nx*5+3]=-2*(vn/(M+h));
    F[nx*5+4]=-2*OMGE*cos(lat)-2*(ve/(N+h));
    F[nx*5+6]=2*ve*OMGE*sin(lat);
    F[nx*5+8]=ve*ve/((N+h)*(N+h))+(vn*vn/((M+h)*(M+h)))-(2*g/(R+h));
    F[nx*5+9]=Cbn[6];
    F[nx*5+10]=Cbn[7];
    F[nx*5+11]=Cbn[8];

    /* Position Error dynamics 0(3),Frv,Frr,0(3),0(3) */
    /* 1st line */
    F[nx*6+3]=1/(M+h);
    F[nx*6+8]=-vn/((M+h)*(M+h));

    /* 2nd line */
    F[nx*7+4]=1/((N+h)*cos(lat));
    F[nx*7+6]=ve*sin(lat)/((N+h)*cos(lat)*cos(lat));
    F[nx*7+8]=-ve/((N+h)*(N+h)*cos(lat));

    /* 3rd line */
    F[nx*8+5]=-1;

    /* Bias  */
    F[nx*9+9]=-1/tbax;
    F[nx*10+10]=-1/tbay;
    F[nx*11+11]=-1/tbaz;
    F[nx*12+12]=-1/tbgx;
    F[nx*13+13]=-1/tbgy;
    F[nx*14+14]=-1/tbgz;

  /* The transition matrix must be found using the inverse Laplace parameter
  Shin (2001, page 37) and Groves(2013) equation 14.72                     */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < nx; j++) {
      PHI[i*nx+j]=I[i*nx+j]+F[i*nx+j]*dt;
    }
  }

  /* G matrix from (Angrisano 2001,page 90)  */
  for (i = 3; i < 6; i++) {
    for (j = 0; j < 3; j++) {
      G[i*(nx-3)+j]=Cbn[(i-3)*3+j];
    }
  }
  for (i = 6; i < 9; i++) {
    for (j = 3; j < 6; j++) {
      G[i*(nx-3)+j]=Cbn[(i-6)*3+(j-3)];
    }
  }
  G[114]=G[127]=G[140]=G[153]=G[166]=G[179]=1;


  free(F); free(I);
}

/* INS-15 error-model system (transition) matrix, Groves (2013, page 587) */
extern void inssysmatrix1(double *PHI, int nx, pva_t *pva_global, imuraw_t *imuobsc, double dt){
  double M,N,R,lat,h,vn,ve,vd,g;
  double tbax,tbay,tbaz,tbgx,tbgy,tbgz; /*acc. and gyro correlation times */
  double fb[3],fn[3],skyfn[9],Cbn[9],winn[3], skywinn[9];
  double *F11,*F12,*F13,*F21,*F22,*F23,*F32,*F33;
  double *F,*I;
  int i,j;


  F=zeros(nx,nx);
  I=eye(nx);
  F11=zeros(3,3);F12=zeros(3,3);F13=zeros(3,3);
  F21=zeros(3,3);F22=zeros(3,3);F23=zeros(3,3);
  F32=zeros(3,3);F33=zeros(3,3);

  /* Values extracted from Angrisano (2010,page 93) for low-cost IMU Crista */
  tbax=tbay=tbaz=270; // [s]
  tbgx=tbgy=tbgz=350; // [s]

  lat=pva_global->r[0]; h=pva_global->r[2];
  ve=pva_global->v[0];vn=pva_global->v[1];vd=pva_global->v[2];
  M=RM(lat);
  N=RN(lat);
  R=sqrt(M*N);
  g=gn(lat,h);

  for (i=0;i<3;i++) fb[i]=imuobsc->fb[i];
//  Cb2n(pva_global.A[0],pva_global.A[1],pva_global.A[2],Cbn);
  for (i=0;i<9;i++) Cbn[i]=pva_global->Cbn[i];
  matmul("NN",3,1,3,1.0,Cbn,fb,0.0,fn);

  winn[0]= OMGE*cos(lat) +(ve/(RN(lat)+h));
  winn[1]=-vn/(RM(lat)+h);
  winn[2]=-OMGE*cos(lat)-(ve*tan(lat)/(RN(lat)+h));

  vec2skew(winn, skywinn);


  /* Building F matrix based on Groves(2013) page 587 */

  /* F11  (14.64) */
  for (i=0;i<9;i++) F11[i]=-skywinn[i];

  /* F12  (14.65) */
  F12[1]=-1/(N+h);
  F12[3]= 1/(M+h);
  F12[7]= tan(lat)/(N+h);

  /* F13  (14.66) */
  F13[0]= OMGE*sin(lat);
  F13[2]= ve/((N+h)*(N+h));
  F13[5]=-vn/((M+h)*(M+h));
  F13[6]= OMGE*cos(lat)+(ve/((N+h)*cos(lat)*cos(lat)));
  F13[8]= -ve*tan(lat)/((N+h)*(N+h));

  /* F21  (14.67) */
  vec2skew(fn,skyfn);
  for (i=0;i<9;i++) F21[i]=-skyfn[i];

  /* F22 (14.68) */
  F22[0]= vd/(M+h);
  F22[1]=-2*OMGE*sin(lat)-2*(ve*tan(lat)/(N+h));
  F22[2]= vn/(M+h);
  F22[3]= 2*OMGE*sin(lat)+(ve*tan(lat)/(N+h));;
  F22[4]= (vd+vn*tan(lat))/(N+h);
  F22[5]= 2*OMGE*cos(lat)+(ve/(N+h));
  F22[6]=-2*(vn/(M+h));
  F22[7]=-2*OMGE*cos(lat)-2*(ve/(N+h));

  /* F23 (14.69) */
  F23[0]=-2*ve*OMGE*cos(lat)-((ve*ve)/((N+h)*cos(lat)*cos(lat)));
  F23[2]=-vn*vd/((M+h)*(M+h))+(ve*ve*tan(lat)/((N+h)*(N+h)));
  F23[3]= 2*OMGE*(vn*cos(lat)-vd*sin(lat))+(ve*vn/((N+h)*cos(lat)*cos(lat)));
  F23[5]=-ve*vd/((N+h)*(N+h))-(vn*ve*tan(lat)/((N+h)*(N+h)));
  F23[6]= 2*ve*OMGE*sin(lat);
  F23[8]= ve*ve/((N+h)*(N+h))+(vn*vn/((M+h)*(M+h)))-(2*g/(R+h));

  /* F32 (14.70) */
  F32[0]= 1/(M+h);
  F32[4]= 1/((N+h)*cos(lat));
  F32[8]=-1;

  /* F33 (14.71) */
  F33[2]=-vn/((M+h)*(M+h));
  F33[3]= ve*sin(lat)/((N+h)*cos(lat)*cos(lat));
  F33[5]=-ve/((N+h)*(N+h)*cos(lat));


  /* Building matrix F in navigation-frame for consumer-grade snesors (14.73) */
  for (i = 0; i < 3; i++) {
    for (j = 12; j < nx; j++) {
      F[i*nx+j]=Cbn[(i)*3+(j-12)];
    }
  }

  for (i = 3; i < 6; i++) {
    for (j = 0; j < 3; j++) {
      F[i*nx+j]=F21[(i-3)*3+(j)];
    }
  }

  for (i=0;i<9;i++) F23[i]=0.0;
  F23[8]=-2*g/reeS(lat);

  for (i = 3; i < 6; i++) {
    for (j = 6; j < 9; j++) {
      F[i*nx+j]=F23[(i-3)*3+(j-6)];
    }
  }

  for (i = 3; i < 6; i++) {
    for (j = 9; j < 12; j++) {
      F[i*nx+j]=Cbn[(i-3)*3+(j-9)];
    }
  }

  for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      F[i*nx+j]=F32[(i-3)*3+(j-3)];
    }
  }

  /* The transition matrix must be found using the inverse Laplace parameter
  Shin (2001, page 37) and Groves(2013) equation 14.72                     */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < nx; j++) {
      PHI[i*nx+j]=I[i*nx+j]+F[i*nx+j]*dt;
    }
  }

  free(F); free(I);
  free(F11);free(F12);free(F13);
  free(F21);free(F22);free(F23);
  free(F32);free(F33);
}

/* System noise covariance matrix Q, Shin (2001, page 38) and
Angrisano (2010,page 93), and Groves (2013, page 591)       */
extern void sysnoise1(double *Q, double dt, int nx){
  double gyro_noise_PSD, accel_noise_PSD, accel_bias_PSD, gyro_bias_PSD;
  double pos_meas_SD, vel_meas_SD;
  double micro_g_to_meters_per_second_squared=9.80665E-6;
  int i,j;

  // Gyro noise PSD (deg^2 per hour, converted to rad^2/s)
  gyro_noise_PSD = (0.02 * D2R / 60)*(0.02 * D2R / 60);
  // Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)
  accel_noise_PSD = (200 *micro_g_to_meters_per_second_squared)*
                   (200 *micro_g_to_meters_per_second_squared);
  // Accelerometer bias random walk PSD (m^2 s^-5)
  accel_bias_PSD = 1.0E-7;
  // Gyro bias random walk PSD (rad^2 s^-3)
  gyro_bias_PSD = 2.0E-12;

  // Position measurement noise SD per axis (m)
  pos_meas_SD = 2.5;
  // Velocity measurement noise SD per axis (m/s)
  vel_meas_SD = 0.1;

  // Determine approximate system noise covariance matrix using (14.82)
  /* The noise matrix Q will be */

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      (i==j?Q[i*(nx)+j]=gyro_noise_PSD * dt:0.0);
    }
  }

  for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      (i==j?Q[i*(nx)+j]=accel_noise_PSD * dt:0.0);
    }
  }

  for (i = 6; i < 9; i++) {
    for (j = 6; j < 9; j++) {
      (i==j?Q[i*(nx)+j]=0.0:0.0);
    }
  }


  for (i = 9; i < 12; i++) {
    for (j = 9; j < 12; j++) {
      (i==j?Q[i*(nx)+j]=accel_bias_PSD * dt:0.0);
    }
  }

  for (i = 12; i < 15; i++) {
    for (j = 12; j < 15; j++) {
      (i==j?Q[i*(nx)+j]=gyro_bias_PSD * dt:0.0);
    }
  }

}

/* System noise covariance matrix Q, Groves (2010,page 592) */
extern void sysnoise(double *Q, double dt, int nx){
  double qa[3],qg[3],qba[3],qbg[3]; /* acc, gyro and repec. biases spectral densities */
  double stda[3],stdg[3],stdba[3],stdbg[3]; /* standard deviations */
  double tba,tbg; /*acc. and gyro correlation times */
  int i,j;

  /* Values extracted from Angrisano (2010,page 93) for low-cost IMU Crista */
  for (i=0;i<3;i++) qa[i]=(300E-6)*(300E-6); // [g/sqrt(1Hz)]
  for (i=0;i<3;i++) qg[i]=(220)*(220); // [deg/h/sqrt(1Hz)]

  for (i=0;i<3;i++) stdba[i]=0.0077; // [m/s^2]
  tba=270; // [s]
  for (i=0;i<3;i++) stdbg[i]=192*D2R*3600; // [deg/h] must it be transformed to rad/s
  tbg=350; // [s]

  /* The biases spectral densities of the Gass-Markov process noises are */
  for (i=0;i<3;i++) qba[i]=2*stdba[i]*stdba[i]/tba; // IN Groves there is no: 2*..()!
  for (i=0;i<3;i++) qbg[i]=2*stdbg[i]*stdbg[i]/tbg;

//  printf("Spectral dessities*********:\n %lf, %lf, %lf, %lf\n",sqrt(qa[0]),sqrt(qg[0]),sqrt(qba[0]),sqrt(qbg[0]) );

  /* The noise matrix Q will be */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      (i==j?Q[i*(nx)+j]=dt*qg[j]:0.0);
    }
  }

  for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      (i==j?Q[i*(nx)+j]=qa[j-3]*dt:0.0);
    }
  }

  for (i = 9; i < 12; i++) {
    for (j = 9; j < 12; j++) {
      (i==j?Q[i*(nx)+j]=qbg[j-9]:0.0);
    }
  }

  for (i = 12; i < 15; i++) {
    for (j = 12; j < 15; j++) {
      (i==j?Q[i*(nx)+j]=qba[j-12]:0.0);
    }
  }


}

/* Initialize attitude and acc. bias covariances due to the coarse alignment -
 Groves (2013, page 198 and 596)
 * modifies Pp matrix                                                             */
extern void attbiasP0(imuraw_t *imu, pva_t *pvap, double *Pp, int nx){
  double T[9], aux[9], aux1[9];
  double *I,*Zro,*zro;
  double *A,*B,*C,*D;
  double *a,*b,*c;
  double fx,fy,fz,fyz,fxyz2;
  double stdyaw,stdba[3];
  int i,j;

  /* Heading and accelerometer bias uncertainties  */
  for (i = 0; i < 3; i++) stdba[i]=imu->stdba[i];
  stdyaw=pvap->P[2];

  fx=imu->fb[0];
  fy=imu->fb[1];
  fz=imu->fb[2];
  fyz=sqrt(fy*fy+fz*fz);
  fxyz2=fx*fx+fy*fy+fz*fz;

  I=eye(3); Zro=zeros(3,3); zro=zeros(1,3);
  A=zeros(6,4);B=eye(6);C=zeros(4,6);D=eye(6);
  a=mat(4,6);b=mat(6,6);c=mat(6,6);

  T[0]=1;T[1]=T[3]=T[6]=0.0;
  T[2]=-sin(pvap->A[0]);
  T[4]=cos(pvap->A[1]);
  T[5]=sin(pvap->A[1])*cos(pvap->A[0]);
  T[7]=-sin(pvap->A[1]);
  T[8]=cos(pvap->A[1])*cos(pvap->A[0]);

  matmul("NN",3,3,3,1.0,pvap->Cbn,T,0.0,aux);
  matmul("TT",3,3,3,1.0,T,pvap->Cbn,0.0,aux1);

  A[1]=fz/(fy*fy+fz*fz);
  A[2]=-fy/(fy*fy+fz*fz);
  A[4]=fyz/fxyz2;
  A[5]=-fx*fy/(fxyz2*fyz);
  A[6]=-fx*fz/(fxyz2*fyz);
  A[22]=A[17]=A[12]=A[11]=1;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      B[i*3+j]=aux[i*3+j];
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      (i==j?C[i*4+j]=stdba[j]*stdba[j]:0.0);
    }
  }
  C[4*6-1]=stdyaw*stdyaw;


  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      D[i*3+j]=aux1[i*3+j];
    }
  }

    matmul("TN",4,6,6,1.0,A,D,0.0,a);
    matmul("NN",6,6,4,1.0,C,a,0.0,b);
    matmul("NN",6,6,6,1.0,B,b,0.0,c);

 /* Update attitude and acc. bias uncertainties, where c is a 6x6 matrix with
 attitude and acc. bias variances and covariances */

 /* The 15-states {de,dv,dr,dba,dbg} system model order */

/* Attitude  */
 for (i = 0; i < 3; i++) {
   for (j = 0; j < 3; j++) {
     (i==j?Pp[i*nx+j]=c[i*6+j]:0.0);
   }
 }

 /* Acc. bias  */
 for (i = 9; i < 12; i++) {
   for (j = 9; j < 12; j++) {
     (i==j?Pp[i*nx+j]=c[(i-6)*6+(j-6)]:0.0);
   }
 }

 /* Should it account for the error correlations??? */
    free(I); free(Zro); free(zro);
    free(A);free(B);free(C);free(D);
    free(a);free(b);free(c);
}

/* Initialize velocity and position state covariance  Shin (2001, page 41)*/
extern void posvelP0 (double *P, int nx, rtk_t *rtk){
  int i,j;
  /* ins 15-states {de,dv,dr,dba,dbg} system model */

  /* Velocity and position uncertainties from GNSS */
  for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      (i==j?P[i*nx+j]=0.0:0.0);
    }
  }

  for (i = 6; i < 9; i++) {
    for (j = 6; j < 9; j++) {
      (i==j?P[i*nx+j]=rtk->sol.qr[i-6]:0.0);
    }
  }

}

/* Generic initialization of P matrix for LC GNSS/INS ------------------------
   From INS-GNSS_Demo_3.m and Initialize_LC_P_matrix.m function
   (Groves, 2013 - MATLAB code)                                               */
extern void initialize_P0(double *P_matrix, int n){
  double init_att_unc,init_vel_unc,init_pos_unc, init_b_a_unc, init_b_g_unc;
  double micro_g_to_meters_per_second_squared=9.80665E-6;
  int i,j;

  // Initial attitude uncertainty per axis (deg, converted to rad)
  init_att_unc = 1*D2R;
  // Initial velocity uncertainty per axis (m/s)
  init_vel_unc = 0.1;
  // Initial position uncertainty per axis (m)
  init_pos_unc = 10;
  // Initial accelerometer bias uncertainty per instrument (micro-g, converted
  // to m/s^2)
  init_b_a_unc = 1000 * micro_g_to_meters_per_second_squared;
  // Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
  init_b_g_unc = 10 * D2R / 3600;

  // Initialize error covariance matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        (i==j?P_matrix[i*n+j]=init_att_unc*init_att_unc:0.0);
    }
  }

  for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
        (i==j?P_matrix[i*n+j]=init_vel_unc*init_vel_unc:0.0);
    }
  }

  for (i = 6; i < 9; i++) {
    for (j = 6; j < 9; j++) {
        (i==j?P_matrix[i*n+j]=init_pos_unc*init_pos_unc:0.0);
    }
  }
  /* For local-navigation 2*pos_unc factor for the vertical components
  Groves(2013, page 595)                                                      */
  P_matrix[8*n+8]=2*(init_pos_unc*init_pos_unc);

  for (i = 9; i < 12; i++) {
    for (j = 9; j < 12; j++) {
        (i==j?P_matrix[i*n+j]=init_b_a_unc*init_b_a_unc:0.0);
    }
  }

  for (i = 12; i < 15; i++) {
    for (j = 12; j < 15; j++) {
        (i==j?P_matrix[i*n+j]=init_b_g_unc*init_b_g_unc:0.0);
    }
  }
}

/* INS/GNSS measurement vector from Shin (2001, pages 40-44) */
extern void measvec(double *z, double* xyz_ini_pos, pva_t *pvap, imuraw_t *imu, double *l){
  double D[9],rgnss[6];
  double wibb[3], wien[3], wenn[3];
  double Wibb[9], Wien[9], Wenn[9];
  double a[3],b[3];
  double c[3],d[3],e[3],f[9],g[3],h[3];
  int i;

  /* Current GNSS position and velpocity */
   /*ecef to geodetic postion	*/
   ecef2pos(xyz_ini_pos, rgnss);
   /* gnss velocity */
  for (i=0;i<3;i++) rgnss[i+3]=pvagnss.v[i];

  D[1]=D[2]=D[3]=D[5]=D[6]=D[7]=0.0;
  D[0]=1/(RM(pvap->r[0])+pvap->r[2]);
  D[4]=1/(RN(pvap->r[0])+pvap->r[2])*cos(pvap->r[0]);
  D[8]=-1;

  for ( i = 0; i < 3; i++) wibb[i]=imu->fb[i];
  vec2skew(wibb,Wibb);

  /* Project wiee (Earth rotation vector) to n-frame	*/
  projvec(pvap->r, wien);
  vec2skew(wien,Wien);

  /* Transport rate in terms of lat long in n-frame	*/
  raten(pvap->r, pvap->v, wenn);
  vec2skew(wenn,Wenn);

  matmul("NN",3,1,3,1.0,pvap->Cbn,l,0.0,b);
  matmul("NN",3,1,3,1.0,D,b,0.0,a);

  /* Measurement vector - Position differences */
  z[0]=(RM(pvap->r[0])+pvap->r[2])*(pvap->r[0] - rgnss[0]- a[0]);
  z[1]=(RN(pvap->r[0])+pvap->r[2])*cos(pvap->r[0])*(pvap->r[1] - rgnss[1]- a[1]);
  z[2]=(pvap->r[2] - rgnss[2]- a[2]);

  matmul("NN",3,1,3,1.0,Wibb,l,0.0,c);
  matmul("NN",3,1,3,1.0,pvap->Cbn,c,0.0,d);
  matmul("NN",3,1,3,1.0,pvap->Cbn,l,0.0,e);
  for (i=0;i<9;i++) f[i]=Wien[i]+Wenn[i];
  matmul("NN",3,1,3,1.0,f,e,0.0,g);
  for (i=0;i<3;i++) h[i]=g[i]-d[i];
/*
  printf("MEASVEC pvagnss: %lf, %lf, %lf\n",pvagnss.v[0],pvagnss.v[1],pvagnss.v[2]);
  printf("MEASVEC RGNSS: %lf, %lf, %lf\n",rgnss[3],rgnss[4],rgnss[5]);
  printf("MEASVEC V: %lf, %lf, %lf\n",pvap->v[0],pvap->v[1],pvap->v[2]);   */

  /* Measurement vector - Velocity differences */
  z[3]=(pvap->v[0] - rgnss[3] + h[0]);
  z[4]=(pvap->v[1] - rgnss[4] + h[1]);
  z[5]=(pvap->v[2] - rgnss[5] + h[2]);

  /* According to Groves(2013, page 600) if the lever arm is not very large(?)
  it may be neglected, becoming
  z[0]=(pvap->r[0] - rgnss[0]);
  z[1]=(pvap->r[1] - rgnss[1]);
  z[2]=(pvap->r[2] - rgnss[2]);
  z[3]=(pvap->v[0] - rgnss[3]);
  z[4]=(pvap->v[1] - rgnss[4]);
  z[5]=(pvap->v[2] - rgnss[5]);  */
}

/* Initialize INS/GNSS design matrix H  Shin (2001, page 41) adapted to
Groves (2013, page 599) parameter order  */
extern void measmatrixH (double *H, int n, int nx, pva_t *pva_global){
  double r[3];
  int i,j;

  /* ins 15-states {de,dv,dr,dba,dbg} system model */

  /* geodetic postion	*/
  for (i = 0; i < 3; i++) r[i]=pva_global->r[i];

  /* Velocity elements in H */
  for (i = 3; i < 6; i++) {
    for (j = 3; j < 6; j++) {
      (i==j?H[i*nx+j]=1:0.0);
    }
  }


 /* Using H matrix from Groves(2013) page 601 (14.115)+Shin terms for making it
  stable */
 for (i = 0; i < nx*n; i++) H[i]=0.0;

 /* Position elements in H */
 for (i = 0; i < 3; i++) {
   for (j = 6; j < 9; j++) {
     (i==(j-6)?H[i*nx+j]=-1.0:0.0);
   }
 }

 /* Position elements in H */
 H[0*nx+6]=-(RM(r[0])+r[2]);
 H[1*nx+7]=-(RN(r[0])+r[2])*cos(r[0]);
 //H[2*nx+8]=1;

/* Velocity elements in H */
 for (i = 3; i < n; i++) {
   for (j = 3; j < 6; j++) {
     (i==j?H[i*nx+j]=-1.0:0.0);
   }
 }

}

/* Shin(2001, page 41)  -------------------------------------------------*/
extern void measnoiseR (double *R, int nv, double* gnss_xyz_ini_cov){
  double pos_meas_SD, vel_meas_SD;
  int i,j;

 /* Gnss position covariance */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      (i==j?R[i*nv+j]=gnss_xyz_ini_cov[j]:0.0); /* CORRECT HERE!- it is cov for xyz */
    }
  }

  /* Gnss velocity covariance ???? */
   for (i = 3; i < 6; i++) {
     for (j = 3; j < 6; j++) {
       (i==j?R[i*nv+j]=gnss_xyz_ini_cov[j]:0.0); /* CORRECT HERE!- it is the position cov */
     }
   }

   /* Using a fixed value based on Groves(2013) LC_KF_Epoch.m function and
   INS_GNSS_Demo_3.m*/

   // Position measurement noise SD per axis (m)
   pos_meas_SD = 2.5;
   // Velocity measurement noise SD per axis (m/s)
   vel_meas_SD = 0.1;


  // 6. Set-up measurement noise covariance matrix assuming all components of
  // GNSS position and velocity are independent and have equal variance.
  for (i=0;i<nv*nv;i++) R[i]=0.0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      (i==j?R[i*nv+j]=pos_meas_SD*pos_meas_SD:0.0);
    }
  }

  for (i = 3; i < nv; i++) {
    for (j = 3; j < nv; j++) {
      (i==j?R[i*nv+j]=vel_meas_SD*vel_meas_SD:0.0);
    }
  }


}

/* Return imu file pointer one complete set of observation  */
extern void imufileback(){
  int j;
  char c;

  for (j = 0; j <= 3 ; j++) {  // Returns 2 sets of obs
    while (c != '\n'){
      fseek(fimu, -1L, SEEK_CUR);
      c = fgetc(fimu);
      fseek(fimu, -1L, SEEK_CUR);
    }
    c=0;
  }
  /* to match the begining of the line after the \n */
  fseek(fimu, 1, SEEK_CUR);
}

/* Velocity from the GNSS coordinates */
void velfromgnss(um7pack_t *imu, float gnsscurrent_time){
  double imuxyz[3],vxyz[3], llh[3], venu[3];
  int i;

  pos2ecef(imu->r, imuxyz);

  /* cartesian velocity */
  if ( norm(pvagnss.r,3)==0.0) {
    /* use previous gnss velocity */
    imu->v[0]=pvagnss.v[0]=0.00;
    imu->v[1]=pvagnss.v[1]=0.00;
    imu->v[2]=pvagnss.v[2]=0.00;
  }else{
    for (i=0;i<3;i++) vxyz[i]=(imuxyz[i]-pvagnss.r[i])/(gnsscurrent_time-pvagnss.sec);
    ecef2enu(imu->r,vxyz,venu);
    venu[2]=-venu[2]; /* Down */
    /* E,N,D velocity */
    imu->v[0]=venu[1];
    imu->v[1]=venu[0];
    imu->v[2]=venu[2];
    pvagnss.v[0]=venu[1];
    pvagnss.v[1]=venu[0];
    pvagnss.v[2]=venu[2];
  }

}

/* INS navigation mechanization -------------------------------------------------
* description: INS navigaiton mehcanization process. Read imu raw data and
* compute a position, velocity and attitude solution
* args   :
      rtk_t rtk I   gnss solution
      pva_t *pvap IO  position, velocity, and attitude structure
      imuraw_t *imuobsp IO  imu measurement and bias structure
* return: status 1: full solution computed, 0: some error
*-----------------------------------------------------------------------------*/
extern int insnav (rtk_t *rtk, um7pack_t *imu, pva_t *pvap, imuraw_t *imuobsp)
{
 um7pack_t imuprev={0};
 int i=0,j=0,count=0;
 FILE* facc;
 FILE* fvp;
 //facc = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/Inertial_acceleration1.txt","a");
 //fvp = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/Inertial_velpos1.txt","a");
 double pos[3],gan[3], wiee[3], wien[3], wenn[3], winn[3], Cen[9], Cbn[9], dt;
 float t0=404616.600;
 float tstatini=1.19349,tstatend=7.42128; /* Static period initial and end times*/
 float gm=0.0,gmcomput=8.227047;
 double lacc[12*3],lgyr[12*3],*x,*A,*v,*w;


  if(imu->status==1){ /* imu package complete check	*/

    /* Initial position (l,l,h) and velocity {ve,vn,vd} from GNSS */
      /* ecef to geodetic postion	*/
      ecef2pos(rtk->sol.rr, pos);
      for (i=0;i<3;i++) imu->r[i]=pos[i];
      //for (i=0;i<3;i++) imu->v[i]=rtk->sol.rr[i+3];
      velfromgnss(imu,time2gpst(rtk->sol.time,NULL)); /* In case there is no doppler observable */

      for (i=0;i<3;i++) pvagnss.r[i]=rtk->sol.rr[i];
      pvagnss.sec=time2gpst(rtk->sol.time,NULL);

      /* Previous PV solution or from GNSS in case is the first   */
      if ( norm(pvap->r,3) == 0.00 && norm(pvap->v,3) == 0.0) {
        for (i=0;i<3;i++) imuprev.r[i]=imu->r[i];
        for (i=0;i<3;i++) imuprev.v[i]=imu->v[i];
      }else{
        for (i=0;i<3;i++) imuprev.r[i]=pvap->r[i];
        for (i=0;i<3;i++) imuprev.v[i]=pvap->v[i];
      }

      /* Using the GNSS velocity as previous solution */
      for (i=0;i<3;i++) imuprev.v[i]=imu->v[i];

      /* Earth rotation vector in e-frame	*/
      wiee[0]=0;wiee[1]=0;wiee[2]=OMGE;

      /* Time span from previous measurement */
      if (!imuobsp->sec) {
        imuprev.sec= imu->sec;
        dt = imu->sec-imuprev.sec; /* In fraction of seconds */
      }else{
      imuprev.sec= imuobsp->sec;
      dt = imu->sec-imuprev.sec; /* In fraction of seconds */
      //printf("Time span: %lf \n",dt);
    }
  //  printf("IMU time prev: %lf and current: %lf, diff: %lf\n",imuprev.sec,imu->sec,dt );

      /* local apparent gravity vector */
      appgrav(imu->r, gan, wiee);

      /* Initial biases
      imu->a[0]=imu->a[0]-0.462904;
      imu->a[1]=imu->a[1]-0.972083;
      imu->a[2]=imu->a[2]-1.64328;    */

    //  fprintf(facc,"%lf %lf %lf %lf \n",imu->sec, imu->a[0], imu->a[1], imu->a[2]); /* accelerations	*/

      /* Measuring average gravity during static period
      if ( tstatini<(imu->sec-t0)/60.0 || (imu->sec-t0)/60.0<tstatend ) {

      for (i=0;i<3;i++) lacc[count*3+i]=imu->a[i];
      for (i=0;i<3;i++) lgyr[count*3+i]=imu->g[i];

      gm+=sqrt(imu->a[0]*imu->a[0]+imu->a[1]*imu->a[1]+imu->a[2]*imu->a[2]);
      count++;
      if (count>=12) { // Determining acc biases with 6 dof
          // perform sequential adjustment
          imucalib(&imu, lacc, lgyr, count*3, gan, wiee);
          count=0;
          memset(lacc, 0.0, sizeof(lacc));
          memset(lgyr, 0.0, sizeof(lgyr));
        }
      }                                                   */

    /*    IMU error compensation    */
  /*  printf("BEFORE CORRECTION\n");
      printf("acc.: %lf, %lf, %lf\n",imu->a[0],imu->a[1],imu->a[2] );
      printf("Gyro.: %lf, %lf, %lf\n",imu->g[0],imu->g[1],imu->g[2] );
      printf("Mag.: %lf, %lf, %lf\n",imu->m[0],imu->m[1],imu->m[2] );
      printf("P: %lf, %lf, %lf\n",imu->r[0],imu->r[1],imu->r[2] );
      printf("V: %lf, %lf, %lf\n",imu->v[0],imu->v[1],imu->v[2] );
      printf("A: %lf, %lf, %lf\n",imu->aea[0],imu->aea[1],imu->aea[2] );  */

      /* Acc. and Gyro biases=(bias+drift) are fedback from latest ING/GNSS solution */
  //  printf("%lf,%lf,%lf,Timediff:%lf\n", imuobsp->ba[0],imuobsp->ba[1], imuobsp->ba[2],dt );
      for(i=0;i<6;i++) imu->Ba[i]=imuobsp->ba[i];
      for(i=0;i<6;i++) imu->Bg[i]=imuobsp->bg[i];

      imuerrorcorr(imu, dt);

/*   printf("AFTER\n");
      printf("acc.: %lf, %lf, %lf\n",imu->a[0],imu->a[1],imu->a[2] );
      printf("Gyro.: %lf, %lf, %lf\n",imu->g[0],imu->g[1],imu->g[2] );
      printf("Mag.: %lf, %lf, %lf\n",imu->m[0],imu->m[1],imu->m[2] );
      printf("P: %lf, %lf, %lf\n",imu->r[0]*R2D,imu->r[1]*R2D,imu->r[2] );
      printf("V: %lf, %lf, %lf\n",imu->v[0],imu->v[1],imu->v[2] );
      printf("A: %lf, %lf, %lf\n",imu->aea[0],imu->aea[1],imu->aea[2] ); */

      /* Update imuobs in b-frame for KF use later */
      for(i=0;i<3;i++) imuobsp->fb[i]=imu->a[i];
      for(i=0;i<3;i++) imuobsp->wibb[i]=imu->g[i];

  /* Attitude Initialization (Groves, 2013)  */

     /* Coarse alignment or use of previous solution */
      if (norm(pvap->A,3)==0.0) {
      //  printf("HERE1\n");
        /* Levelling and gyrocompassing, update imu->aea[] vector */
        coarseAlign(imu->a, imu->g, imu->v, gan, imu->aea);
        //coarseAlign(imu, gan);
      }else{
      //  printf("HERE2\n");
        for (i=0;i<3;i++) imuprev.aea[i]=pvap->A[i];
      }

  /* Update attitude  */
      /* Project wiee (Earth rotation vector) to n-frame	*/
      projvec(imuprev.r, wien);

      /* Transport rate in terms of lat long in n-frame	*/
      raten(imuprev.r, imuprev.v, wenn);

      /* winn	*/
      winncomp(imu->r, imu->v, winn);

      /* Initial rotation matrix from e to n-frame */
      Cb2n(imu->aea[0],imu->aea[1],imu->aea[2],Cbn);

      /* Rotation matrix from e to n-frame	?*/
      Ce2n(pos[0],pos[1],Cen);

      /* Attitude integration, updating Cbn matrix	*/
      attint(imu, dt, wien, wenn, Cbn);

  /* Update position and velocity  */
      /* Acceleration transformation, updates imu->a[] */
      accb2n(imu, Cbn);

      /* Acceleration gravitation and coriolis correction, updates imu->a[] */
      accgravcor(imu, &imuprev, gan, wien, wenn, dt);

      /* Estimate velocity and position in n-frame	*/
      velposinteg(imu, &imuprev, dt);

  /* Update navigation solution - PVA  */

      Cbn2euler(imu,Cbn);
      for (i=0;i<3;i++) pvap->A[i]=imu->aea[i];
      for (i=0;i<9;i++) pvap->Cbn[i]=Cbn[i];
      for (i=0;i<3;i++) pvap->v[i]=imu->v[i];
      for (i=0;i<3;i++) pvap->r[i]=imu->r[i];
      pvap->time = imu->time;


  /*    printf("acc.: %lf, %lf, %lf\n",imuobsp->fb[0],imuobsp->fb[1],imuobsp->fb[2] );
      printf("Gyro.: %lf, %lf, %lf\n",imuobsp->wibb[0],imuobsp->wibb[1],imuobsp->wibb[2] );
      printf("P: %lf, %lf, %lf\n",pvap->r[0]*R2D,pvap->r[1]*R2D,pvap->r[2] );
      printf("V: %lf, %lf, %lf\n",pvap->v[0],pvap->v[1],pvap->v[2] );
      printf("A: %lf, %lf, %lf\n",pvap->A[0],pvap->A[1],pvap->A[2] );        /*     */

    // printf("HERE 5 Acc. %f %lf %lf %lf\n",imu->sec, imu->a[0], imu->a[1], imu->a[2] );
    // printf("HERE 5 Gyro. %f %lf %lf %lf\n",imu->sec, imu->g[0], imu->g[1], imu->g[2] );
    // printf("HERE 5 Att. %f %lf %lf %lf\n",imu->sec, imu->aea[0], imu->aea[1], imu->aea[2] );

    /* Correct measurements	*/
    // imucorr(&imu);
    /* Body-to-Inertial frame convertion */
    //imubdy2in(&imu, &imuprev);
    /* Acceleration normal correction and convertion to m/s/s	*/
    //imuacc(&imu);
    //fprintf(facc,"%f %f %f %f \n",imu->sec, imu->a[0], imu->a[1], imu->a[2]); /* accelerations	*/
    /* Estimate position and velocity	*/
    //imuvelpos(&imu, &imuprev);

    /* Printing imu data into a file 	*/
  //  fprintf(facc,"%f %lf %lf %lf \n",imu->sec, imu->a[0], imu->a[1], imu->a[2]); /* accelerations	*/
    //fprintf(facc,"%lf %lf %lf %lf \n",imu->sec, imu->g[0], imu->g[1], imu->g[2]); /* rate-rotations	*/
    //fprintf(facc,"%lf %lf %lf %lf %lf %lf %lf \n",imu->sec, imu->g[0], imu->g[1], imu->g[2], imu->aea[0], imu->aea[1], imu->aea[2]); /* rate and Euler-rotations	*/
    //fprintf(facc,"%f %f %f %f \n",imu->sec, imu->m[0], imu->m[1], imu->m[2]); /* magnetometers	*/

    /* Printing velocity and/or positions into a file 	*/
  //  fprintf(fvp,"%f %lf %lf %lf \n",imu->sec, imu->v[0], imu->v[1], imu->v[2]); /* vx,vy,vz */
    //fprintf(fvp,"%lf %lf\n", imu->sec, imu->r[2]);// imu->r[1]); // imu->r[2]); /* up */
    //fprintf(fvp,"%lf %lf\n", imu->r[0], imu->r[1]);// imu->r[1]); // imu->r[2]); /* EN */

    /* Update measurement info	*/
        imuobsp->sec=  imu->sec;
        imuobsp->time= imu->time;
        memset(&imu, 0, sizeof(um7pack_t));
   }else {
     imuobsp->sec=  imu->sec;
     imuobsp->time= imu->time;
    // fclose(facc);
    // fclose(fvp);
     return -1;
   }
//printf("Average gravity while static: %lf\n", gm/count );

 //fclose(facc);
 //fclose(fvp);
 /*Printing imu data*/
// imuaccplot();
 //imugyroplot();
 //imumagplot();
 //imueulerplot();
 /* Printing imu pos and vel	*/
 //imuvelplot();
 //imuposplot();
 //imuUpplot();
 return 1;
}

/* INS navigation mechanization -------------------------------------------------
* description: INS navigaiton mehcanization process. Read imu raw data and
* compute a position, velocity and attitude solution
* args   :
      rtk_t rtk I   gnss solution
      pva_t *pvap IO  position, velocity, and attitude structure
      imuraw_t *imuobsp IO  imu measurement and bias structure
* return: status 1: full solution computed, 0: some error
* Based on Groves (2013 and his Nav_equations_NED.m code)
*-----------------------------------------------------------------------------*/
extern int insnav1 (double* gnss_xyz_ini_pos, double* gnss_enu_vel,
             float ini_pos_time, um7pack_t *imu_curr_meas, pva_t *pvap, imuraw_t *imuobsp)
{
 um7pack_t imuprev={0};
 int i=0,j=0,count=0;
 double pos[3],gan[3], wiee[3], wien[3], wenn[3], winn[3], Cen[9], Cbn[9], dt;

 double alpha_ib_b[3], mag_alpha, Alpha_ib_b[9];
 double old_wenn[3], *I=eye(3), omega_sum_vec[3], omega_sum_matrix[9];
 double old_C_b_n[9], ave_C_b_n[9];
 double cte, cte1, first_term[9], second_term[9], Alpha_squared[9];
 double sumterms[9], first_final_term[9], second_final_term[9];
 double f_ib_b[3],f_ib_n[3];
 double v_eb_n[3], old_v_eb_n[3],omega_sum2_vec[3],omega_sum2_matrix[9];
 double v_second_vector[3], v_inter_vector[3], v_inter_minus_second[3];
 double old_h_b, old_L_b, old_lambda_b, h_b, L_b, lambda_b;
 double llh[3], C_new_old[9], C_b_n[9], C_aux[9], omegas_sum_vec[3];
 double omegas_sum_matrix[9], att_first_term[9];
 double G = 9.80665;

 //fprintf(facc,"%f %lf %lf %lf\n",imu_curr_meas->sec, imu_curr_meas->a[0], imu_curr_meas->a[1], imu_curr_meas->a[2]); // accelerations
 //fprintf(fgyr,"%f %lf %lf %lf\n",imu_curr_meas->sec, imu_curr_meas->g[0]*R2D, imu_curr_meas->g[1]*R2D, imu_curr_meas->g[2]*R2D); // rate-rotations


  if(imu_curr_meas->status==1){ /* imu package complete check	*/

   /* External position and velocity for initialization */
    if ( norm(pvap->r,3) == 0.0 ) {
  //    printf("\n Current GNSS p\n");
      ecef2pos(gnss_xyz_ini_pos, pos);
      for (i=0;i<3;i++) imuprev.r[i]=pos[i];
      for (i=0;i<3;i++) imuprev.v[i]=0.0;
      }else{
    //    printf("\n Previous GNSS p\n");
        ecef2pos(pvagnss.r, pos);
        for (i=0;i<3;i++) imuprev.r[i]=pos[i];
        for (i=0;i<3;i++) imuprev.v[i]=pvagnss.v[i];
      } /* previous gnss sol */

      ecef2pos(gnss_xyz_ini_pos, pos); /* Forcing gnss as initial positioning */
      for (i=0;i<3;i++) imuprev.r[i]=pos[i];

      for (i=0;i<3;i++) imuprev.v[i]=pvagnss.v[i];
      // When imumaps only
      //imuprev.v[0]=gnss_enu_vel[1]; imuprev.v[1]=gnss_enu_vel[0];imuprev.v[2]=-gnss_enu_vel[2];


    /* Initial position (l,l,h) and velocity {vn,ve,vd} from GNSS */
      /* ecef to geodetic postion	*/
      ecef2pos(gnss_xyz_ini_pos, pos);
      for (i=0;i<3;i++) imu_curr_meas->r[i]=pos[i];

      /* NED velocities */
      // When imumap only
      //imu_curr_meas->v[0]=gnss_enu_vel[1]; imu_curr_meas->v[1]=gnss_enu_vel[0];imu_curr_meas->v[2]=-gnss_enu_vel[2];

      /* Velocity - in case there is no doppler observable */
      //velfromgnss(imu,ini_pos_time);

      /* Earth rotation vector in e-frame	*/
      wiee[0]=0;wiee[1]=0;wiee[2]=OMGE;

      /* Time span from previous measurement */
      if (!imuobsp->sec) {
        imuprev.sec= imu_curr_meas->sec;
        dt = imu_curr_meas->sec-imuprev.sec; /* In fraction of seconds */
      }else{
      imuprev.sec= imuobsp->sec;
      dt = imu_curr_meas->sec-imuprev.sec; /* In fraction of seconds */
    }

    /* Previous measurements  */
    for (i=0;i<3;i++){
    imuprev.a[i]=imuobsp->fb[i];
    imuprev.g[i]=imuobsp->wibb[i];
    }

    /* Turn-on initial biases             //exp4
    imu_curr_meas->a[0]=imu_curr_meas->a[0]-0.4453528;//-0.455972;//*G;
    imu_curr_meas->a[1]=imu_curr_meas->a[1]-1.162248;//-1.02027;//-0.972083*G;
    imu_curr_meas->a[2]=imu_curr_meas->a[2];//+8.0764127;//+8.1495;//-1.64328*G; */
    // exp1
    imu_curr_meas->a[0]=imu_curr_meas->a[0]-0.455972;//*G;
    imu_curr_meas->a[1]=imu_curr_meas->a[1]-1.02027;//*G;
    imu_curr_meas->a[2]=imu_curr_meas->a[2]+8.1495;//-1.64328*G;*/
  /* 2nd experiment
    imu_curr_meas->a[0]=imu_curr_meas->a[0]-0.455972;
    imu_curr_meas->a[1]=imu_curr_meas->a[1]-1.41527;
    imu_curr_meas->a[2]=imu_curr_meas->a[2]+8.20592;//+8.18169;//-1.64328*G;*/
    /* Down aixs gravity compensation */
    //imu_curr_meas->a[2]=G+imu_curr_meas->a[2];

      /* local apparent gravity vector */
      appgrav(imu_curr_meas->r, gan, wiee);

      /* Acc. and Gyro biases=(bias+drift) are fedback from latest ING/GNSS solution */
  //  printf("%lf,%lf,%lf,Timediff:%lf\n", imuobsp->ba[0],imuobsp->ba[1], imuobsp->ba[2],dt );
      for(i=0;i<6;i++) imu_curr_meas->Ba[i]=imuobsp->ba[i];
      for(i=0;i<6;i++) imu_curr_meas->Bg[i]=imuobsp->bg[i];

      imuerrorcorr(imu_curr_meas, dt);

    /* Update imuobs measurements in b-frame for KF use later */
      for(i=0;i<3;i++) imuobsp->fb[i]=imu_curr_meas->a[i];
      for(i=0;i<3;i++) imuobsp->wibb[i]=imu_curr_meas->g[i];

  /* Attitude Initialization (Groves, 2013)  */

     /* Coarse alignment or use of previous solution */
      if (norm(pvap->A,3) == 0.0 || norm(pvap->r,3) == 0.0) {
        /* Levelling and gyrocompassing, update imu_curr_meas->aea[] vector */
        //coarseAlign(imu_curr_meas, gan);
        coarseAlign(imu_curr_meas->a, imu_curr_meas->g, imu_curr_meas->v, gan, imu_curr_meas->aea);
        /* Initial rotation matrix from e to n-frame */
        Cb2n(imu_curr_meas->aea[0],imu_curr_meas->aea[1],imu_curr_meas->aea[2],old_C_b_n);
      }else{
        /* Previous attitude  */
        //for (i=0;i<9;i++) old_C_b_n[i]=pvap->Cbn[i];
        //for (i=0;i<3;i++) imuprev.aea[i]=pvap->A[i];

        /* Heading use only */
        //imuprev.aea[2]=pvap->A[2];
        /* Using previous meas to initialize  */
        coarseAlign(&imuprev.a, &imuprev.g, &imuprev.v, gan, &imuprev.aea);
        //coarseAlign(&imuprev, gan);
        /* Initial rotation matrix from e to n-frame */
        Cb2n(imuprev.aea[0],imuprev.aea[1],imuprev.aea[2],old_C_b_n);
      }

      /* Small angle approximation Groves (previous version of 2013) page 29
      for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
        (i==j?old_C_b_n[i*3+j]=1-old_C_b_n[i*3+j]:(old_C_b_n[i*3+j]=old_C_b_n[i*3+j]));
      }
    }*/

  /* Preliminaries */
  /* Calculate attitude increment, magnitude, and skew-symmetric matrix  */
  for (i=0;i<3;i++) alpha_ib_b[i] = imu_curr_meas->g[i]*dt;
  mag_alpha = norm(alpha_ib_b,3);
  vec2skew(alpha_ib_b, Alpha_ib_b);

  /* From (2.123) , determine the angular rate of the ECEF frame
  w.r.t the ECI frame, resolved about NED  */
  projvec(imuprev.r, wien);

  /* From (5.44), determine the angular rate of the NED frame
  w.r.t the ECEF frame, resolved about NED  */
  raten(imuprev.r, imuprev.v, old_wenn);

  for (i=0;i<3;i++) omega_sum_vec[i] = wien[i]+old_wenn[i];
    vec2skew(omega_sum_vec, omega_sum_matrix);

  /* Specific-force frame transformation */
  /* Calculate the average body-to-ECEF-frame coordinate transformation
  matrix over the update interval using (5.84) and (5.86)    */
if (mag_alpha>1.E-8){
  cte = (1 - cos(mag_alpha)) / (mag_alpha*mag_alpha);
  cte1 = (1 - (sin(mag_alpha) / mag_alpha)) / (mag_alpha*mag_alpha);
  for (i=0;i<9;i++) first_term[i]=I[i]+cte*Alpha_ib_b[i];
  matmul("NT", 3, 3, 3, 1.0, Alpha_ib_b, Alpha_ib_b, 0.0, Alpha_squared);
  for (i=0;i<9;i++) second_term[i]=cte1*Alpha_squared[i];
  for (i=0;i<9;i++) sumterms[i]= first_term[i] + second_term[i];
  matmul("NN", 3, 3, 3, 1.0, old_C_b_n, sumterms, 0.0, first_final_term);
  matmul("NN", 3, 3, 3, 0.5, omega_sum_matrix, old_C_b_n, 0.0, second_final_term);
  for (i=0;i<9;i++) ave_C_b_n[i] = first_final_term[i] - second_final_term[i]*dt;
}else{
  matmul("NN", 3, 3, 3, 0.5, omega_sum_matrix, old_C_b_n, 0.0, second_final_term);
  for (i=0;i<9;i++) ave_C_b_n[i] = old_C_b_n[i] - second_final_term[i]*dt;
}

 /* Transform specific force to ECEF-frame resolving axes using (5.86) */
 for (i=0;i<3;i++) f_ib_b[i]=imu_curr_meas->a[i];
 matmul("NN", 3, 1, 3, 1.0, ave_C_b_n, f_ib_b, 0.0, f_ib_n);

 /* UPDATE VELOCITY    */
  // From (5.54)
  for (i=0;i<3;i++) omega_sum2_vec[i] = old_wenn[i] + 2*wien[i];
  vec2skew(omega_sum2_vec, omega_sum2_matrix);
  for (i=0;i<3;i++) old_v_eb_n[i] = imuprev.v[i];
  matmul("NN", 3, 1, 3, 1.0, omega_sum2_matrix, old_v_eb_n, 0.0, v_second_vector);
  for (i=0;i<3;i++) v_inter_vector[i] = f_ib_n[i] + gan[i];
  for (i=0;i<3;i++) v_inter_minus_second[i] = dt*(v_inter_vector[i] - v_second_vector[i]);

  for (i=0;i<3;i++) v_eb_n[i] = old_v_eb_n[i] + v_inter_minus_second[i];

 /* UPDATE CURVILINEAR POSITION   */
 old_L_b=imuprev.r[0];
 old_lambda_b=imuprev.r[1];
 old_h_b=imuprev.r[2];

  // Update height using (5.56)
  h_b = old_h_b - 0.5 * dt * (old_v_eb_n[2] + v_eb_n[2]);

  // Update latitude using (5.56)
  L_b = old_L_b + 0.5 * dt * (old_v_eb_n[0] / (RM(old_L_b) + old_h_b) +
     v_eb_n[0] / (RM(old_L_b) + h_b));

 // Update longitude using (5.56)
 lambda_b = old_lambda_b + 0.5 * dt * (old_v_eb_n[1] / ((RN(old_L_b) +
    old_h_b) * cos(old_L_b)) + v_eb_n[1] / ((RN(L_b) + h_b) * cos(L_b)));


 /* ATTITUDE UPDATE               */
  /* From (5.44), determine the angular rate of the NED frame
  w.r.t the ECEF frame, resolved about NED  */
  llh[0]=L_b; llh[1]=lambda_b; llh[2]=h_b;
  raten(llh, v_eb_n, wenn);

  /* Obtain coordinate transformation matrix from the new attitude w.r.t. an
  inertial frame to the old using Rodrigues' formula, (5.73)   */
  if (mag_alpha>1.E-8){
    cte = sin(mag_alpha) / mag_alpha;
    cte1 = (1 - cos(mag_alpha)) / (mag_alpha*mag_alpha);
    for (i=0;i<9;i++) first_term[i]=I[i]+(cte*Alpha_ib_b[i]);
    matmul("NT", 3, 3, 3, 1.0, Alpha_ib_b, Alpha_ib_b, 0.0, Alpha_squared);
    for (i=0;i<9;i++) second_term[i]=cte1*Alpha_squared[i];
    for (i=0;i<9;i++) C_new_old[i] =  first_term[i] + second_term[i];
  }else{
    for (i=0;i<9;i++) C_new_old[i] = I[i]+Alpha_ib_b[i];
  }

  /* Update attitude using (5.77)         */
  matmul("NN", 3, 3, 3, 1.0, old_C_b_n, C_new_old, 0.0, C_aux);
  for (i=0;i<3;i++) omegas_sum_vec[i] = dt*(wien[i]+0.5*old_wenn[i]+0.5*wenn[i]);
  vec2skew(omegas_sum_vec, omegas_sum_matrix);
  for (i=0;i<9;i++) att_first_term[i] = I[i]-omegas_sum_matrix[i];
  matmul("NN", 3, 3, 3, 1.0, att_first_term, C_aux, 0.0, C_b_n);


  /* Update navigation solution - PVA  */
      Cbn2euler(imu_curr_meas,C_b_n);
  //printf("A - endinsnav: %lf, %lf, %lf\n",imu_curr_meas->aea[0],imu_curr_meas->aea[1],imu_curr_meas->aea[2] );

      for (i=0;i<3;i++) pvap->A[i]=imu_curr_meas->aea[i];
      for (i=0;i<9;i++) pvap->Cbn[i]=C_b_n[i];
      for (i=0;i<3;i++) pvap->v[i]=v_eb_n[i];
      for (i=0;i<3;i++) pvap->r[i]=llh[i];
      pvap->time = imu_curr_meas->time;

    /* Update measurement info	*/
        imuobsp->sec=  imu_curr_meas->sec;
        imuobsp->time= imu_curr_meas->time;
        //memset(&imu_curr_meas, 0, sizeof(um7pack_t));
   }else {
     printf("INS not processed!\n");
     *pvap=pva_global;
     pvap->time = imu_curr_meas->time;
     imuobsp->sec=  imu_curr_meas->sec;
     imuobsp->time= imu_curr_meas->time;

     return -1;
   }



 return 1;
}
