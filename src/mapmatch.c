/**/
#include "rtklib.h"
#include "satinsmap.h"

/* ppp.c functions  */
#define SQR(x)      ((x)*(x))
#define NP(opt)     ((opt)->dynamics?9:3) /* number of pos solution */
#define IC(s,opt)   (NP(opt)+(s))      /* state index of clocks (s=0:gps,1:glo) */
#define IT(opt)     (IC(0,opt)+NSYS)   /* state index of tropos */
#define NR(opt)     (IT(opt)+((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3)))
                                       /* number of solutions */
#define IB(s,opt)   (NR(opt)+(s)-1)    /* state index of phase bias */
#define NX(opt)     (IB(MAXSAT,opt)+1) /* number of estimated states */

/* measurement error variance ------------------------------------------------*/
static double varerr1(int sat, int sys, double el, int type, const prcopt_t *opt)
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
/* Map-Matching functions ----------------------------------------------------*/

/* Initial reference lane buffer ---------------------------------------------*/
/* Fetch 15m of positions surrounding the the current lane
 file pointer (last closest point) and fill lane structure.
 args:	I	none
        O	modify the lane.buff structure
 obs.: In case is at the begining or end of lane file,
 the buffer is filled accordingly -------------------------------------------*/
void extern inibuff (){
 char str[150];
 double r[3];
 int i,j,k=0;

 /* Clean lane.buff vector in case fill size change	*/
 memset(lane.buffer, 0, sizeof(lane.buffer));

 //printf("Lane position: %d \n",ftell(fp_lane) );

 /* Reading the reference lane points 	*/
 /* Checking if file pointer is in the beginning or end of file	*/
 if ( ftell(fp_lane) <= BUFFSIZE*lane.sizeline || ftell(fp_lane) > lane.size-BUFFSIZE*lane.sizeline ){
  if (ftell(fp_lane) <= BUFFSIZE*lane.sizeline){ 		/* Beginning of file */
   rewind(fp_lane); /* returning to the begninning of reference lanes file */
   for (i=0; i< BUFFSIZE*2 ; i++){
	  fgets(str, 150, fp_lane);
	  sscanf(str, "%lf %lf %lf", &r[0], &r[1], &r[2]);
	  for(j=0;j<3;j++) {
	    lane.buffer[i*3+j] = r[j];
	  }
    k++;
   }
 }else{
   fseek(fp_lane, -(BUFFSIZE+1)*lane.sizeline, SEEK_CUR); /* End of file */
   i=0;
   while (fgets(str, 150, fp_lane) != NULL){
    sscanf(str, "%lf %lf %lf", &r[0], &r[1], &r[2]);
    for(j=0;j<3;j++) lane.buffer[i*3+j] = r[j];
	  i++;
    k++;
   }
  }
}else{						/* Anywhere in the file */
  /* returning BUFFSIZE positions from the closest point position */
   fseek(fp_lane, -(BUFFSIZE+1)*lane.sizeline, SEEK_CUR);
   for (i=0; i< BUFFSIZE*2 ; i++){
	  fgets(str, 150, fp_lane);
	  sscanf(str, "%lf %lf %lf", &r[0], &r[1], &r[2]);
	  for(j=0;j<3;j++) {
		 lane.buffer[i*3+j] = r[j];
	  }
    k++;
  }
 }
 lane.sizebuffer=k;
}

/* Closest point (3D) ----------------------------------------------------------
* Compute the closest point between the initial buffer space positions and the
* rover SPP position
* args:	I double rov_p[]	SPP cartesian coordinates (3x1) [m]
*        O int lane.currsearch 	closest p. buffer position
*-----------------------------------------------------------------------------*/
void clstpt(double *rov_p){
 double src_area, r[3], pos[3];
 double d[BUFFSIZE*2], mindist=0.0;

 int i=0, j=0;
 char str[150];

 /* Determine the closest point position 	*/
 for (i=0; i< lane.sizebuffer ; i++){
	for (j=0;j<3;j++) r[j]=lane.buffer[i*3+j];
	d[i] = dist(rov_p, r);
 }
 /* Closest point position in the buffer size	*/
 lane.currsearch = min(d, BUFFSIZE*2, &mindist);

 /* Repositioning the pointer in the lane file	*/
 fseek(fp_lane, (-(BUFFSIZE*2-lane.currsearch) )*lane.sizeline, SEEK_CUR);
}

/* Azimuth ---------------------------------------------------------------------
* Compute the azimuth between two points from the lane buffer
* args:	I space	enu Lane buffer (lane.sizebufferx3) [m]
*       O double az[2]  vector with azimuth and its quadrant [rad,1,2,3,4]
*-----------------------------------------------------------------------------*/
void azmt(double *space, double *az){
  double enu[3],enun[3];
  double DE,DN,delta;
  double azm,q; //azimuth and its quadrant
  int j;

  /* Closest point positions and next neighbour point */
	for (j=0;j<3;j++) enu[j]=space[lane.currsearch*3+j];

  if(lane.currsearch==lane.sizebuffer) {
    for (j=0;j<3;j++) enun[j]=space[(lane.currsearch-1)*3+j];
    DE=enu[0]-enun[0]; DN=enu[1]-enun[1];
  }
  else {
    for (j=0;j<3;j++) enun[j]=space[(lane.currsearch-1)*3+j];
    DE=enun[0]-enu[0]; DN=enun[1]-enu[1];
  }


  delta = ( DE / DN );
  //                                angle in radians
  if (delta<0){ // 2nd or 4th Q
      if ( DN < 0){
          // 2nd quadrant
          azm = 180*(PI/180) + atan(delta); //adding because the value is negative
          q=2;
        }else{
          // 4th quadrant
          azm = 360*(PI/180) + atan(delta);//adding because the value is negative
          q=4;
        }
   }else if ( DE < 0){   // 1st or 3rd Q
          // 3rd quadrant
          azm = 180*(PI/180) + atan(delta);
          q=3;
        }
      else{ // 1st quadrant
          azm = atan(delta);
          q=1;
      }

  az[0]=azm;
  az[1]=q;
  return az;
}

/* 2D confidence ellipse parameters --------------------------------------------
* Compute the parameters of a 2D error ellipse with 95% confidence
* args:	I	double* cov SPP en covariance matrix (4x1)
*       O double* maj_axis ellipse semi-major axis
*       O double* min_axis ellipse semi-minor axis
*       O double* alpha the ellipse orientation from x(e) axis
*-----------------------------------------------------------------------------*/
void ell2D95(double* cov, double* maj_axis, double* min_axis, double* alpha){
 double eval[2],evec[4];
 double M,smax,smin,theta;
 int i,j;

  /* semi-major and minor axes from eigenvalues and vector
  if (0 != rs(2, cov, eval, 1, evec) ){
    printf("Error computing eigenvalue and vector\n");
  }*/
      theta = 0.5*atan( 2*cov[1]/(cov[0] - cov[3]));
      M = sqrt(  4*(cov[1]*cov[1]) + (cov[0] - cov[3])*(cov[0] - cov[3]) );
      smax = 0.5*( cov[0] + cov[3] + M);
      smin = 0.5*( cov[0] + cov[3] - M);

//printf("Ell par: %lf,%lf,%lf\n",cov[0],cov[1],cov[3]);

  /* The 95% confidence ellipse semi-major and minor lengths are */
  *maj_axis = (sqrt( smax ) )*2.447;
  *min_axis = (sqrt( smin ) )*2.447;
  *alpha = theta;
  //printf("Ell par: %lf,%lf,%lf \n",*maj_axis,*min_axis,*alpha);
  /*
    *maj_axis = 2*sqrt(5.991*eval[1]);
    *min_axis = 2*sqrt(5.991*eval[0]);
    *alpha = atan(evec[4]/evec[2]); //atan( Vmax(y or n) / Vmax(x or e) )*/
  //printf("Eig. min,max: %lf,%lf \n",eval[0],eval[1]);
  //printf("Ellipse par.: %lf,%lf,%lf\n",*maj_axis,*min_axis,*alpha);
  //printf("Ellipse pa1.: %lf,%lf,%lf\n",(sqrt( smax ) )*2.447,(sqrt( smin ) )*2.447,theta);
}

/* Curve fitting ---------------------------------------------------------------
* Fits a curve over some points using least squares to determine its parameters
* args:	I	double* pts set of enu points (BUFFSIZE*2 x 3)
*       O double* a,b,c parameters of a quadratic function,such as:
* f(x)=a*x^2+b*x+c, in this case: n = a*e^2+b*e+c
*-----------------------------------------------------------------------------*/
void curvfit(double* pts, double *a, double *b, double *c){
  int i,j,n,smpsize = 40; // In this case 40*0.2 = 20 meters from te closest point
  double* y,*x,*A,*At;
  double *Q, X[3];

  /* Getting 10m radius sample from the closest point */
  if (lane.currsearch - (smpsize/2) <= 0 ) { //clstp at the beginning
    n=lane.currsearch + (smpsize/2);
    y=mat(n,1);
    x=mat(n,1);
    A=mat(n,3);
    At=mat(n,3);
    j=0;
    for (i = 0; i < lane.currsearch+(smpsize/2); i++) {
        x[j]=pts[i*3];
        y[j]=pts[i*3+1];
        j++;
      //  printf(" %lf, %lf\n", y[i],x[i]);
    }
  } else {
    if (lane.currsearch+(smpsize/2) >= lane.sizebuffer) { //clst at the end
      n=lane.sizebuffer-(lane.currsearch-(smpsize/2));
      y=mat(n,1);
      x=mat(n,1);
      A=mat(n,3);
      At=mat(n,3);
      j=0;
      for (i = lane.currsearch-(smpsize/2); i < lane.sizebuffer; i++) {
          x[j]=pts[i*3];
          y[j]=pts[i*3+1];
        //  printf(" %lf, %lf\n", y[j],x[j]);
          j++;
      }
    }else{ //clst anywhere in the middle
      n=(lane.currsearch+(smpsize/2))-(lane.currsearch-(smpsize/2));
      y=mat(n,1);
      x=mat(n,1);
      A=mat(n,3);
      At=mat(n,3);
      j=0;
      for (i = lane.currsearch-(smpsize/2); i < lane.currsearch+(smpsize/2); i++) {
          x[j]=pts[i*3];
          y[j]=pts[i*3+1];
        //  printf(" %lf, %lf\n", y[j],x[j]);
          j++;
      }
    }
  }

  /* Design matrix  */
  for (i = 0; i < n; i++) {
    A[i*3]  = x[i]*x[i];
    A[i*3+1]= x[i];
    A[i*3+2]= 1;
  }
  /* Transpose A to lsq() function */
  for (i = 0; i < n; i++) {
    for (j = 0; j < 3; j++) {
      At[j*n+i]=A[i*3+j];
    }
  }

 Q=mat(3,3);

 /* LSQ  */
  if (0 != lsq(At, y, 3, n, X, Q) ){
    printf("Error estimating curve parameters\n");
  }

  *a=X[0];
  *b=X[1];
  *c=X[2];
  //printf("Curve a,b,c par.: %8.10lf,%8.10lf,%8.10lf \n",X[0],X[1],X[2]);
  //printf("Curve a,b,c par.: %8.10lf,%8.10lf,%8.10lf \n",*a,*b,*c);
  free(x);free(y);free(A);free(At);free(Q);
  memset(X, 0, sizeof(X));
}

/* Ellipse curve instersection -------------------------------------------------
* Determine if lanes and ellipse intersect and returns the lane positions that
* are within it.
* args:	I	double e0,n0 ellipse origin points in east,north local coordinates
*       I double a,b,alpha ellipse major and minor axis and orientation par.
*       I double* space Lane buffer in enu coordinates (lane.sizebuffer x 3)
*       O int* intersecflags flag vector with the position of the lane points
*              that intersects(=1) and the ones that do not intersect(=0)
*-----------------------------------------------------------------------------*/
void ellcurvIntersec(double e0, double n0, double a, double b, double alpha,
  double* space, int* intersecflags){
    double ep,np,de,dn,beta,t,eell,nell;
    int i,j;

    FILE* fp1;
    //fp1=fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/exp1_curvecorrep.txt","a");

    /* Determining angle between search points and ellipse origin */
     for (i=0; i< lane.sizebuffer ; i++){
       ep=space[i*3+0];//e
       np=space[i*3+1];//n
       de=ep-e0;
       dn=np-n0;

       /* Lane-ellipse center angle - beta */
       if (dn>0) {
         if (de>0) {
           beta = atan(dn/de);
         }else beta = 90*(D2R)+atan(de/dn);
       }else if (de<0) {
         beta = 180*(D2R)+atan(dn/de);
       }else {beta = 270*(D2R)+atan(de/dn);}

       /* Transform the Lane-ellipse center angle to the ellipse orientation */
      if ( abs(beta-alpha) < 0.00001 ) {
        t=alpha;
      }else if (beta-alpha>0.0) {
        t=alpha+(beta-alpha);
      }else {t=360*(D2R)-(alpha-beta);}

      /* Ellipse correspondent point  */
      eell=e0+a*cos(t);
      nell=n0+b*sin(t);
     //eell=e0+a*cos(alpha)*cos(t)-b*sin(alpha)*sin(t);
      //nell=n0+a*cos(alpha)*sin(t)+b*cos(t)*sin(alpha);
      if(sqrt((e0-ep)*(e0-ep)+(n0-np)*(n0-np))<
      sqrt((e0-eell)*(e0-eell)+(n0-nell)*(n0-nell))) {
        /* Intersects or is within the Ellipse  */
        intersecflags[i]=1;
        //printf("FLAG 1!!!!!!!!!!!!!!!!!!!!\n");
      }else{
        /* Does not intersect or is within the Ellipse  */
        intersecflags[i]=0;
        //printf("FLAG 0!!!!!!!!!!!!!!!!!!!!\n");
      }
    }

    for ( i = 0; i < lane.sizebuffer; i++) {
      if(intersecflags[i]==1){//printf("Cand.:%d\n",intersecflags[i]);
      }
    }

double k;
    /* For printing purposes  */
    for (k = 0.0; k < 2*3.141592; k=k+0.25) {
      //eell=e0+a*cos(alpha+k);
      //nell=n0+b*sin(alpha+k);
      eell=e0+a*cos(k)*cos(alpha)-b*sin(k)*sin(alpha);
      nell=n0+a*cos(k)*sin(alpha)+b*cos(alpha)*sin(k);
      //fprintf(fp1, "%lf %lf\n",eell,nell);
    }

    //fprintf(fp1, "\n");
    //fclose(fp1);
}
/* Enhanced search space function ---------------------------------------------
* Compute the search space based on the SPP rover position and its uncertaity
* args:	I	double r SPP enu position (3x1)
*       I	double  Q SPP enu full covariance matrix (9x1)
*       I double* space Lane buffer in enu coordinates (lane.sizebuffer x 3)
*       O searchcandidates[] a position flag (1) vector of the candidates [lane.sizebufferx1]
*-----------------------------------------------------------------------------*/
void ensearchspace (double *r, double *Q, double* space, int* searchcandidates){
 double en[2],enQ[4];
 double flag[BUFFSIZE*2];
 double maj_axis,min_axis,alpha;
 double a,b,c;
 int i,j,count=0;

 FILE *fp1;
 //fp1=fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/exp1_enh_space_qntts.txt","a");

 en[0]=r[0];
 en[1]=r[1];

 for (i = 0; i < 2; i++) {
   for (j = 0; j < 2; j++) {
     enQ[i*2+j]=(double)Q[i*3+j];
   }
 }

 /* Compute confidence ellipse parameters */
 ell2D95(enQ,&maj_axis,&min_axis,&alpha);

 /* Determine the curve parameters of surspacerounding points to the closest point */
 //curvfit(space,&a,&b,&c);

 /* Intercept ellipsoid area with the curve and determine the enhanced area   */
 ellcurvIntersec(en[0],en[1],maj_axis,min_axis,alpha,space,searchcandidates);

//fprintf(fp1,"%lf %lf %lf %lf %lf\n",en[0], en[1], maj_axis, min_axis, alpha*R2D);
//fclose(fp1);
}

/* solution to covariance ----------------------------------------------------*/
static void soltocov(const sol_t sol, double *P)
{
    P[0]     =sol.qr[0]; /* xx or ee */
    P[4]     =sol.qr[1]; /* yy or nn */
    P[8]     =sol.qr[2]; /* zz or uu */
    P[1]=P[3]=sol.qr[3]; /* xy or en */
    P[5]=P[7]=sol.qr[4]; /* yz or nu */
    P[2]=P[6]=sol.qr[5]; /* zx or ue */
}
/* Locking phase ---------------------------------------------------------------
* Determine a better rover position within the ehanced search space with phase
* observables per satellite
* args:	I rtk_t *rtk rtk structure
*       I rr[3] Closest point position in ecef [m]
*       I enspc[nx3] Enhanced search space pointer, n is the number of
*       candidates, in ecef [m]
*       I int nc size of Enhanced search space
*       I	rtk_t *sol solution structure with the SPP position and covariance
*       I const obsd_t *obs structure containing the observables for one epoch
*       I const obsd_t *nav structure containing the satellite positions
*       I int n number of observation data
*       I const *rs satellite position and velocity vector
*       I const double *dts satellite clock bias
*
* obs.:dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
* returns an improved rover position in ecef coordinates [3x1]
*-----------------------------------------------------------------------------*/
double lckphase(rtk_t *rtk, double* rr, double* enspc, int nc, const obsd_t *obs,
  const nav_t *nav, int n, const double *rs, const double *dts, double trop,
  int satn, double *meas, double *vare, double *var, double vart, double *varm,
  double *dtdx, double *v, double *H, double *R, double *azel, double *x, int sat,
double *mmcand){

    double *r1,*N, *Res, minres=0;
    double L1,L2,r=0.0;
    float dtsat;
    int svh[MAXOBS], posmin;
    const double *lam=nav->lam[obs->sat-1];

    prcopt_t *opt=&rtk->opt;
    double disp[3],pos[3],e[3],dantr[NFREQ]={0};
    double dants[NFREQ]={0},dtrp=0.0;
    int i,j,k,sys,nv=0,nx=rtk->nx,brk,tideopt;

    double *v1,*H1,*R1,*azel1,*x1,*var1;

    nv=n*rtk->opt.nf*2;
    v1=mat(nc,2);
    //matcpy(v1,v,nv,1);

  /*  azel1=zeros(2,n); var1=mat(1,n);
    x1=mat(rtk->nx,1);
    matcpy(x1,x,rtk->nx,1);
    nv=n*rtk->opt.nf*2; v1=mat(nv,1); H1=mat(rtk->nx,nv); R1=mat(nv,nv);

    matcpy(v1,v,nv,1);
    matcpy(H1,H,rtk->nx,nv);
    matcpy(var1,var,MAXOBS,2); matcpy(azel1,azel,2,n);
    matcpy(R1,R,nv,nv);*/

    nv=0;

    N=mat(2,1); /* For 1 satellite - L1 and L2  */
    r1=mat(1,nc);
    Res=mat(1,nc);

    /* Compute the approximate geometric ranges between candidates and satellite
    positions */
       //printf("GOT here in phaselock?, Number of candidates: %d\n", nc);

    if(nc>0){
        for (i = 0; i < nc; i++) {
          /* geometric distance/azimuth/elevation angle */
          if ((r1[i]=geodist(rs,enspc+i*3,e))<=0.0) continue;
          //  printf("Dist: %lf, Geom: %lf\n",dist(rs,enspc+i*3), geodist(rs,enspc+i*3,e) );
          //printf("LCK ENSPCE: %lf %lf %lf\n",*(enspc+i*3),*(enspc+i*3+1),*(enspc+i*3+2));
        //  printf("LCKSPC-XYZ: %lf %lf %lf\n",enspc[i*3+0], enspc[i*3+1], enspc[i*3+2]);

        //  printf("r-P: %lf, r-L: %lf\n", meas[1]-(r[i]-CLIGHT*(*dts)+(*trop)),meas[0]-(r[i]-CLIGHT*(*dts)+(*trop)));
        r=(double)r1[i];
      //  printf("Geom. cand.: %lf, %lf \n",r1[i],r);

         for (j = 0; j < 2; j++) {

           if (meas[j]==0.0) continue;

           v1[nv]=meas[j]-r;

           if (sys!=SYS_GLO) {
               v1[nv]-=x[IC(0,opt)];
               //H1[IC(0,opt)+nx*nv]=1.0;
           }
           else {
               v1[nv]-=x[IC(1,opt)];
              // H1[IC(1,opt)+nx*nv]=1.0;
           }

           if (j==0) {
               v1[nv]-=x[IB(obs[satn].sat,opt)];
               //H1[IB(obs[satn].sat,opt)+nx*nv]=1.0;
           }

          // if (j==0) printf("Phase res.: %lf\n", v1[nv]);
        //   else      //printf("Code res.: %lf\n", v1[nv]);

          nv++;
         }


  /*       CPOIED FROM PPP_RES FUNCTION, USE IT CAREFULLY TO JUDGE THE Residuals  */
    /*    for (j=0;j<2;j++) { // for phase and code

            if (meas[j]==0.0) continue;

            for (k=0;k<nx;k++) H1[k+nx*nv]=0.0;

            v1[nv]=meas[j]-r;

            for (k=0;k<3;k++) H1[k+nx*nv]=-e[k];

            if (sys!=SYS_GLO) {
                v1[nv]-=x[IC(0,opt)];
                H1[IC(0,opt)+nx*nv]=1.0;
            }
            else {
                v1[nv]-=x[IC(1,opt)];
                H1[IC(1,opt)+nx*nv]=1.0;
            }
            if (opt->tropopt>=TROPOPT_EST) {
                for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                    H1[IT(opt)+k+nx*nv]=dtdx[k];
                }
            }
            if (j==0) {
                v1[nv]-=x[IB(obs[satn].sat,opt)];
                H1[IB(obs[satn].sat,opt)+nx*nv]=1.0;
            }
            var1[nv]=varerr(obs[satn].sat,sys,azel1[1+satn*2],j,opt)+varm[j]+vare[satn]+vart;

            if (j==0) rtk->ssat[sat-1].resc[0]=v[nv];
            else      rtk->ssat[sat-1].resp[0]=v[nv];
//printf("GOT here in phaselock 15\n");
            // test innovation
#if 0
            if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
#else
            if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno&&sys!=SYS_GLO) {
#endif
                trace(2,"ppp outlier rejected %s sat=%2d type=%d v=%.3f\n",
                      time_str(obs[satn].time,0),sat,j,v[nv]);
                printf("PPP OUTLIER DETECTED: sat=%2d\n",sat);
                //rtk->ssat[sat-1].rejc[0]++;
                continue;
            }
            if (j==0) //rtk->ssat[sat-1].vsat[0]=1;

            printf("Meas: %lf, r: %lf and v[nv]: %lf\n", meas[j], r, v1[nv]);
            nv++;
         }*/

      }
        /* Judge here which geometric distance candidate on r[] offers the smallest
        residuals? per candidate, using the phase measurement  */
        /* Compute "ambiguity" terms: Nij = frac_phase-rho(geom. dist candidates)*/
      //  printf("Sat: %d: Nof cand.%d\n", obs->sat+1, nc);

    //  for (i = 0; i < nc; i++) {
    //    printf("Phase: %lf and Code: %lf measurements res.\n",v1[i*2],v1[i*2+1]);
    //  }
        posmin = min(v1,nc,&minres);
        r=(double)r1[posmin];
      // printf("LCK candidateXYZ: %lf, %lf, %lf\n",enspc[posmin*3], enspc[posmin*3+1], enspc[posmin*3+2] );

        for (k = 0; k < 3; k++) mmcand[k]=enspc[posmin*3+k];
      //  printf("Smaller residuals at lanepos: %d, Resvec: %d - %lf, respective Geom.: %lf\n",
      //  posmin, min(v1,nc,&minres), minres, r);

    //   L1=obs[i].L[0]*lam[0];
    //   L2=obs[i].L[1]*lam[1];
    //   N[i]=L1-r[i];
      // printf("Sat: %d,\nL1 and L2: %lf, %lf\nP1 and P2: %lf, %lf",
    //   obs[i].sat,L1,L2,obs[i].P[0],obs[i].P[1]);
      // printf("Computed rho and P: %lf, %lf\n",r[i],obs[i].P[0]);
  //  }
    }
    //  printf("End of phaselock, r:%lf\n",r);

      free(r1); free(N); free(Res);
      free(v1);
      //free(v1); free(H1); free(R1);free(azel1);free(x1);free(var1);
      return r;
}

/* Map-Matching ----------------------------------------------------------------
* Find a position, within the enhanced search space, that
* best represents the receiver's true position
* args:	I	rtk_t *rtk solution structure with the SPP position and covariance
*       I const obsd_t *obs structure containing the observables for one epoch
*       I int n number of observed satellites
*       I const obsd_t *nav structure containing the satellite positions
* returns an improved rover position in ecef coordinates [3x1]
*-----------------------------------------------------------------------------*/
extern int match (rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
  double *enhspce){//input rtk_t
 double mm_r[3], r_spp[3], pos[3], pos1[3], pos2[3], enuspp[3],enu[3],r[3];
 double pos0[3],r0[3],enu0[3],dpos[3],enuQ[9],P[9],dr[3];
 double *enhaux;
 double lat0=45.94312323*D2R, lon0=-66.64106835*D2R, h0=52.7309; // BMO corner 1
 double space[BUFFSIZE*2*3];
 int ensrcspcflag[BUFFSIZE*2*3], count;
 float r_sppq[6];
 int i,j,k;
 FILE *fp1;
 double rs1[3],dts1[2],trop1;
 double rechgt[3]; /* Receiver height {e,n,u} */
 double rho=0.0, tow;
 double epoch[6],e[3],azq[2]; /* vector of day/time {year,month,day,hour,min,sec} */
 double recbody[3]; /* Platform gnss antena offsets in body-frame
 {x(nav.direc),y(right),z(down)} in meters */
 recbody[0]=0;recbody[1]=0.113;recbody[2]=1.629;

 //fp1 = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/exp1Closest_point_report.txt","a");
 //fp1 = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/spp_pos_test.txt","a");
 //fp1 = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/exp1_enhancedspace.txt","a");
 //fp1 = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/exp1_mmpoint.txt","a");
// fp1 = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/clst_lane_heights_rtklib4.txt","a");
 //fp1 = fopen("/home/emerson/Desktop/Connected_folders/SatInsMap/out/3dclosestpointdist_rtklib.txt","a");

 tow=time2gpst(rtk->sol.time,NULL);
 time2epoch(rtk->sol.time, epoch);

 for(i=0;i<3;i++) r_spp[i]=rtk->sol.rr[i];
 for(i=0;i<6;i++) r_sppq[i]=rtk->sol.qr[i];

 /* Fetching data to the first buffer structure (lane.buffer)	*/
 inibuff();

 /* Closest point computation  			*/
 clstpt(&r_spp);

 /* Lane to local and lane offset with the antenna height amount  */
 /* llh to enu from lat0, h0 */
 pos0[0]=lat0;pos0[1]=lon0;pos0[2]=h0; /* Manually defining origin but it can be done
 this way: Enu baseline form - example from solution.c Rtklib function
 outenu(): */
 pos2ecef(pos0, r0);

 /* Vehicle height */
 //for (j=0;j<3;j++) rechgt[j]=rtk->opt.antdel[0][j];

 /* Lane buffer conversion to local {enu}  */
  for (i=0; i< lane.sizebuffer; i++){
    for (j=0;j<3;j++) r[j]=lane.buffer[i*3+j];
    //pos1[2]=pos1[2]+rechgt[2]; /* Account the vehicle height  into the Lanes */
    //pos2ecef(pos1, r);
    //for (j=0;j<3;j++) lane.buffer[i*3+j]=r[j]; /* updating lane buffer  */
    for (j = 0; j < 3; j++) dr[j]=r[j]-r0[j];
    ecef2enu(pos0, dr, enu);
    //d2lgs(lat0, h0, pos1, enu);
    //printf("ENU-d2l: %lf,%lf,%lf\n",enu[0],enu[1],enu[2])
    //for (j = 0; j < 3; j++) dr[j]=r[j]-r0[j];
    //ecef2enu(pos, dr, enu);
    //printf("ENU-rtk: %lf,%lf,%lf\n",enu[0],enu[1],enu[2]);
    for (j=0; j < 3; j++) space[i*3+j] = enu[j];
 }

 /* SPP position and covariance to local (enu) */
 for (j = 0; j < 3; j++) dr[j]=r_spp[j]-r0[j];
 ecef2enu(pos0, dr, enuspp);
 ecef2pos(r_spp, pos2);
 //d2lgs(lat0, h0, pos2, enuspp);
 soltocov(rtk->sol,P);
 covenu(pos2, P, enuQ);


  /* Closest point positions to local */
	for (j=0;j<3;j++) r[j]=lane.buffer[lane.currsearch*3+j];
  for (j = 0; j < 3; j++) dr[j]=r[j]-r0[j];
  ecef2enu(pos0, dr, enu);
  ecef2pos(r, pos1);
  //d2lgs(lat0, h0, pos1, enu);

  /* Closest point enu check */
  //fprintf(fp1, "%lf %lf\n",enu[0],enu[1]);
//  fprintf(fp1, "%lf %lf\n", (epoch[3]+epoch[4]/60+epoch[5]/3600), pos1[2]);
  //fprintf(fp1, "%lf %lf %lf \n", tow, space[2],enuspp[2]);

 /* Enhanced search space computation       */
 ensearchspace(enuspp, enuQ, space, &ensrcspcflag);

count=0;
 for (i = 0; i < lane.sizebuffer; i++) {
   if (ensrcspcflag[i]==1) {
     count++;
   }
 }

 //enhspce=mat(count,3);

 k=0;
 for (i = 0; i < lane.sizebuffer; i++) {
   if (ensrcspcflag[i]==1) {
       for (j=0;j<3;j++) r[j]=lane.buffer[i*3+j];
       for (j=0;j<3;j++) {
       enhspce[k] = r[j];
       k++;
       }
     //fprintf(fp1, "%lf %lf\n",space[i*3],space[i*3+1]);
   }
}

/* Lane buffer conversion to local {enu} */
 for (i=0; i< lane.sizebuffer; i++){
   for (j=0;j<3;j++) r[j]=lane.buffer[i*3+j];
   for (j = 0; j < 3; j++) dr[j]=r[j]-r0[j];
   ecef2enu(pos0, dr, enu);
   for (j=0; j < 3; j++) space[i*3+j] = enu[j];
}

 /* Starting point of search is the closest point	*/
 for(i=0;i<3;i++) pos[i]=lane.buffer[lane.currsearch*3+i];//[_+1],[_+2]

 /* Distance between SPP-Closest point check */
 //fprintf(fp1, "%lf %lf\n", (epoch[3]+epoch[4]/60+epoch[5]/3600), dist(r_spp, pos) );

 //printf("GEOM. DIST. to SPP: %lf\n", geodist(rs,r_spp,e) );

 /* Locking phase */
 //printf("CLSPXYZ: %lf, %lf, %lf\n",pos[0],pos[1],pos[2]);

// rho=lckphase(rtk, pos, enhspce, count, obs, nav, n, rs, dts, trop, satn, meas,
// vare,var,vart,varm,dtdx,v,H,R,azel,x,sat,mmcand);

// free(enhspce);
 //fclose(fp1);
 //printf("Geom at Match: %lf\n",rho);
 return count;
 //printf("Finished match \n");
}
