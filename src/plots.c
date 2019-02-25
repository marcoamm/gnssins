/*  */
#include "rtklib.h"
#include "satinsmap.h"

char *outpath[] = {"../out/"};

/* struct declaration -------------------------------------------------------*/
  typedef struct {	/* Satellite residuos states structure	*/
    float *P;	/* Pseudorange residuos (m) */
    float *L;	/* Carrier phase residuos (m) */
    float *ep;	/* Epoch time of solution in TOW (sec)	*/
 }satres;

 typedef struct {	/* Clock bias states structure	*/
    float *clk1;	/* GPS clock bias (ns)	*/
    float *clk2;	/* GLONASS clock bias (ns)	*/
    float *ep;	/* Epoch time of solution in TOW (sec)	*/
 }clkres;

 typedef struct {	/* Troposphere parameter states structure	*/
    float *ztd;	/* zenith total delay (m) float	*/
    float *ztdf;	/* zenith total delay (m) fixed	*/
    float *ep;	/* Epoch time of solution in TOW (sec)	*/
 }tropres;

/* MEMs IMU functions -------------------------------------------------------*/

/* Data IMU plotting ----------------------------------------------------------
* description: read data from imu functions and plot with GNuplot
* args   :
*-----------------------------------------------------------------------------*/
extern void imuaccplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(min)' \n"
              "set ylabel 'ACCELERATIONS (m/s/s)' \n"
	      //"set xrange [0:40] \n"
 	      //"set term postscript eps enhanced color\n"
	      //"set output '%s accelerations3.ps'\n"
	      "set autoscale \n"
        //"set yrange [-2.5:2.5] \n"
	      "set grid \n", outpath[0]);

 fprintf(gp,"plot '%s' u ($1):($4) w l title \"acc. z\" ,\
       '%s' u ($1):($2) w l title \"acc. x\" ,\
       '%s' u ($1):($3) w l title \"acc. y\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imugyroplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(min)' \n"
              "set ylabel 'RATE GYROSCOPES (deg/s)' \n"
	      "set xrange [0:40] \n"
 	      //"set term postscript eps enhanced color\n"
	      //"set output '%s gyroscopes3.ps'\n"
	      "set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($7) w l title \"gyro. z\" ,\
          '%s' u ($1):($5) w l title \"gyro. x\" ,\
          '%s' u ($1):($6) w l title \"gyro. y\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imumagplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(min)' \n"
 	      "set ylabel 'UNITLESS ()' \n"
	      "set xrange [0:40] \n"
 	      //"set term postscript eps enhanced color\n"
	      //"set output '%s magnetometers3.ps'\n"
	      "set autoscale \n"
	      "set grid \n",outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($2) w l title \"mag. x\" ,\
         '%s' u ($1):($3) w l title \"mag. y\" ,\
         '%s' u ($1):($4) w l title \"mag. z\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imueulerplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(min)' \n"
              "set ylabel 'ATTITUDE ANGLE (deg)' \n"
	      "set yrange [-200:200] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s euler_angles.ps'\n"
	      //"set autoscale \n"
        "set xrange [242141.328125:244778.140625] \n"
	      "set grid \n", outpath[0]);

  /* printing Euler angles */
  fprintf(gp,"plot '%s' u ($1):($10) w l title \"yaw(z)\" ,\
         '%s' u ($1):($8) w l title \"roll(x)\" ,\
         '%s' u ($1):($9) w l title \"pitch(y)\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imuposplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'y(lat)(m)' \n"
	      "set ylabel 'x(long)(m)' \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s positions_higway.ps'\n"
        /* General view  */
        "set yrange [45.942:45.9432] \n"  //BMO field general view
        "set xrange [-66.6418:-66.6402] \n"
        /* #point 1
        "set yrange [45.943:45.9432] \n"
        "set xrange [-66.6412:-66.6410] \n"*/
        //"set autoscale \n"
        "set yrange [45.93:45.98] \n"  // Kinematic positioning course dataset general view
        "set xrange [-66.675:-66.63] \n"
        //"set yrange [45.958:45.963] \n"  // Downtown Interruptions
        //"set xrange [-66.64:-66.635] \n"
        //"set yrange [45.932:45.937] \n"  // Highway
        //"set xrange [-66.650:-66.655] \n"
	      "set grid \n", outpath[0]);

        // Obstructions: 396546.04 < $1 < , Off-road: 396838.72 - end
       // With points: w points pointtype 1.4
/**/
  fprintf(gp,"plot '../out/PPP_car_back(copy).pos' u ($4):($3) with points pointsize 0.4 lt 3 title \" SPP \", \
   '%s' u ($3):($2) w lp lt 2 pt 1 ps 0.2 title \" ins/gnss position \"", filename);
/*
  fprintf(gp,"plot '../data/Lanes_5_llh.csv' u ($2):($1) with points pointtype 1.4 ps 0.3 lc rgb \"red\" title \" Ground-truth \" ,\
  '%s' u ($3): ( $1>396546.04 && $1<396838.72 ? ($2) : 1/0 ) w lp lt 2 pt 1 ps 0.2 title \" ins/gnss position \" ,\
  '../out/out_ppp_tllh.txt' u ($3): ( $1>396546.04 && $1<396838.72 ? ($2) : 1/0 ) with points pointsize 0.4 lt 3 title \" SPP \" ", filename);
*/

/*
  fprintf(gp,"plot '%s' u ($3):($2) w l ls .5 title \"ins/map position\" ,\
         '%s' u ($6):($5) w points pointtype 1.4 ps 0.2 title \"lane position\" \n", filename, filename);
*/
 fflush(gp);
 fclose(gp);
 return;
}

extern void imuUpplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(min)' \n"
              "set ylabel 'UP(m)' \n"
 	      //"set term postscript eps enhanced color\n"
	      //"set output '%s up3.ps'\n"
	      //"set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u (($1-39477.0)/60):($4) w l ls .5 title \"ins/map up or height\" ,\
         '%s' u (($1-39477.0)/60):($7) w points pointtype 1.4 ps 0.1 title \"lane up or height\" \n", filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imuvelplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'VELOCITIES (m/s)' \n"
	      "set yrange [-20:20] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s velocities.ps'\n"
        "set xrange [242141.328125:244778.140625] \n"
	      //"set autoscale \n"
        //"set yrange [-10:10] \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($5) w l title \"vel. N\" ,\
         '%s' u ($1):($6) w l title \"vel. E\" ,\
         '%s' u ($1):($7) w l title \"vel. D\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imuaccbiasplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'IMU ACC. BIAS' \n"
	      "set yrange [-0.05:0.05] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s IMU_acc_bias.ps'\n"
	      //"set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($2) w l title \"acc. bias x\" ,\
         '%s' u ($1):($3) w l title \"acc. bias y\" ,\
         '%s' u ($1):($4) w l title \"acc. bias z\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void imugyrobiasplot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'IMU GYROS. BIAS' \n"
	      "set yrange [-0.05:0.05] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s IMU_gyro_bias.ps'\n"
	      //"set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($5) w l title \"gyr. bias x\" ,\
         '%s' u ($1):($6) w l title \"gyr. bias y\" ,\
         '%s' u ($1):($7) w l title \"gyr. bias z\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void KF_att_stds_plot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'ATTITUDE STDs' \n"
	      "set yrange [0:2] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s KF_att_std.ps'\n"
	      //"set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($2) w l title \"roll(x) std\" ,\
         '%s' u ($1):($3) w l title \"phi(y) std\" ,\
         '%s' u ($1):($4) w l title \"yaw(z) std\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void KF_vel_stds_plot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'VELOCITY STDs' \n"
	      "set yrange [0:2] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s KF_vel_std.ps.ps'\n"
	      //"set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($5) w l title \"x vel. std\" ,\
         '%s' u ($1):($6) w l title \"y vel. std\" ,\
         '%s' u ($1):($7) w l title \"z vel. std\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void KF_pos_stds_plot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'POSITION STDs' \n"
	      "set yrange [0:2] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s KF_pos_std.ps'\n"
	      //"set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($8) w l title \"lat. std\" ,\
         '%s' u ($1):($9) w l title \"long. std\" ,\
         '%s' u ($1):($10) w l title \"height std\" \n", filename, filename, filename);

 fflush(gp);
 fclose(gp);
 return;
}

extern void KF_clock_plot(char* filename){
 FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - PLotting INS data  */
  fprintf(gp, "set xlabel 'TIME(s)' \n"
	      "set ylabel 'CLOCK' \n"
	      "set yrange [0:2] \n"
 	      "set term postscript eps enhanced color\n"
	      "set output '%s KF_clk.ps'\n"
	      "set autoscale \n"
	      "set grid \n", outpath[0]);

  fprintf(gp,"plot '%s' u ($1):($2) w l title \"clock\" \n", filename);

 fflush(gp);
 fclose(gp);
 return;
}


/* map-matching plotting functions ------------------------------------*/

/* Receive map-matching quantities for plot  --------------------------*/
extern void mapmatchplot (){//double t, obsd_t *obs, int n
	int i=0, j=0, k=0, sizef;
  FILE * gp = popen ("gnuplot -persistent", "w");
	FILE *f1, *f2, *f3,*f4;
	f4 = fopen("../data/exp1_Lanes_enu.txt", "r");
  f1 = fopen("../out/3dcand.txt", "r");
	f2 = fopen("../out/spp.txt", "r");
	f3 = fopen("../out/finalpos.txt", "r");

	double lat, lon,h, x, y, z, pos[3], pos1[3], r[3], pos0[3], enu[3],t;
	char str[150];
  double lat0=45.943131456*D2R, lon0=-66.641066493*D2R, h0=57.5157;
  pos0[0]=lat0;pos0[1]=lon0;pos0[2]=h0;

	/* Gnuplot plot commands   - Receiver trajectory  */
	fprintf(gp, "set xlabel 'E' \n"
                "set ylabel 'N' \n"
		//"set term postscript eps enhanced color\n"
		//"set output '%s 1st_run_analysis.ps' \n"
		//"set title '1st run around BMO field' \n"
    //"set yrange [-60:70]\n"
    //"set xrange [-21500:-21370]\n"
    //"set yrange [-7:-10]\n"
  //  "set xrange [-21475:-21455]\n"
    "set xrange [-66.6978:-66.5717]\n"
    "set yrange [45.875:45.95]\n"
  //  "set autoscale \n"
		"set grid \n",outpath[0]);

  fprintf(gp,"plot '../out/gnssins_pos.txt' u ($3):($2) with points pointtype 1 ps 0.3 lc rgb \"red\" title \" GNSS/INS \" \n");

/*
	fprintf(gp,"plot '-' with points pt 2 ps 0.1 lc rgb \"red \" title \"Ref. Trajectory \", \
 		'-' w points ps 0.5 title \" Closest P. \", \
		 '-' w lp ls 2 title \" PPP sol. \" \n"
   ); //,'-' w lp ls 6 title \" Closest Pt. \" \n ");
   fprintf(gp,"plot '-' with points pt 2 ps 0.05 lc rgb \"black \" title \"Ground-truth \", \
   '-' with points pt 2 ps 0.05 lc rgb \"red \" title \"Final PPP-kine Sol \", \
       '-' w points pt 2 ps 0.25 title \" SPP solution \", \
      '-' w points pt 1 ps 0.25 title \" Best lane candidate \" \n");*/
/*
      fprintf(gp,"plot '-' with points pt 2 ps 0.05 lc rgb \"red \" title \"Final PPP-kine Sol \", \
          '-' w points pt 2 ps 0.25 title \" SPP solution \", \
         '-' w points pt 1 ps 0.25 title \" Best lane candidate \" \n");*/


  while ( fgets(str, 150, f4)!= NULL ){
  		sscanf(str, "%lf %lf %lf", &enu[0], &enu[1], &enu[2]);
      //fprintf(gp, "%lf %lf\n", enu[1],enu[0]);
  }
  //fflush(gp);
  //fprintf(gp, "e\n");

	while ( fgets(str, 150, f3)!= NULL ){
		sscanf(str, "%lf %lf %lf %lf", &t,&r[0], &r[1], &r[2]);
    if (16.55<t&&t<16.65) { //16.55<t&&t<16.65,16.65<t&&t<16.75,16.75<t&&t<16.883
      ecef2pos(r,pos);
      d2lgs(lat0, h0, pos, enu);
      //ecef2enu(pos0, pos, enu);
    //  fprintf(gp, "%lf %lf\n", enu[1],enu[0]);
    }
	 }
	//fflush(gp);
	//fprintf(gp, "e\n");

	while ( fgets(str, 150, f2)!= NULL ){
    sscanf(str, "%lf %lf %lf %lf", &t,&r[0], &r[1], &r[2]);
    if (16.55<t&&t<16.65) {
      ecef2pos(r,pos);
      d2lgs(lat0, h0, pos, enu);
      //ecef2enu(pos0, pos, enu);
    //  fprintf(gp, "%lf %lf\n", enu[1],enu[0]);
    }
	 }
	//fflush(gp);
	//fprintf(gp, "e\n");

	while ( fgets(str, 150, f1)!= NULL ){
    sscanf(str, "%lf %lf %lf %lf", &t,&r[0], &r[1], &r[2]);
  if (16.55<t&&t<16.65) {
    ecef2pos(r,pos);
    d2lgs(lat0, h0, pos, enu);
    //ecef2enu(pos0, pos, enu);
  //  fprintf(gp, "%lf %lf\n", enu[1],enu[0]);
    }
	}
	fflush(gp);
	fclose(gp);
	fclose(f1);	fclose(f2);	fclose(f3); fclose(f4);
}

/* Receive map-matching quantities for 2D trajectory plot  --------------------------*/
extern void mapmatchplot1 (){//double t, obsd_t *obs, int n
  FILE * gp = popen ("gnuplot -persistent", "w");

	/* Gnuplot plot commands   - Receiver trajectory  */
	fprintf(gp, "set xlabel 'Local east [m]' offset 0,-1,0 font 'Arial,19' \n"
                "set ylabel 'Local north [m]' offset -2,0,0 font 'Arial,19' \n"
                "set tics font 'Arial,19' \n" //adjust tics size
                "set key font 'Arial,19' \n" //ajusdt legend size
                "set label font 'Arial,19' \n"
                "set bmargin at screen 0.15 \n"
		//"set term postscript eps enhanced color\n"
    //"set term postscript eps font ',15'\n"
    //"set output '%sENconvergences.ps' \n"
		//"set output '%sEN_1stReconvergence.ps' \n"
    //"set output '%sEN_2ndReconvergence.ps' \n"
    //"set output '%sEN_3rdReconvergence_font.eps' \n"
    "set key top right \n"
    /* #1
    "set xrange [-2:2] \n"
    "set yrange [-2:2] \n" */
    /* #12
    "set xrange [25:30] \n"
    "set yrange [-35:-30] \n"*/
    /*#2
    "set xrange [60:65] \n"
    "set yrange [-80:-75] \n" */
    /*#21
    "set xrange [20:30] \n"
    "set yrange [-110:-100] \n"  */
    /*#3
    "set xrange [10:13] \n"
    "set yrange [-122:-119] \n"*/
    /*#3 4 */
    "set xrange [-30:-5] \n"
    "set yrange [-100:-75] \n"
    /*#4
    "set xrange [-55:-50] \n"
    "set yrange [-44:-40] \n" */
    /*#4 1
    "set xrange [-30:-5] \n"
    "set yrange [-30:-5] \n"*/
		//"set autoscale \n"
    /*"set xrange [60:65] \n"
    "set yrange [-80:-75] \n"*/
		"set grid \n",outpath[0]);

    /* EN convergence plots   ($2):( $1<16.4833? ($3) : 1/0 )
    static  $1<16.57
    1s conv. 16.783<$1 &&  $1<16.8 / 2 conv $1>16.80 / 3rd conv 16.850<$1<16.866  */
     fprintf(gp,"plot 1/0 with points pointtype 1.4 ps 1.1 lc rgb \"red\" title \" Ground-truth \" ,\
    '../out/Lanes_enu1.txt' u ($1):($2) with points pointtype 1.4 ps 0.3 lc rgb \"red\" notitle ,\
    '../out/PPK_enu.txt' u ($2):( $1>16.8  ? ($3) : 1/0 ) with points pointsize 1.1 lt 2 title \" RTK \" ,\
    '../out/PPP_enu.txt' u ($2):( $1>16.8 ? ($3) : 1/0 ) with points pointsize 1.1 lt 3 title \" PPP-only \" ,\
    '../out/PPP_cand_enu.txt' u ($2):( $1>16.8 ? ($3) : 1/0 ) with points pointsize 1.1 lt 4 title \" PPP-mod \" \n");

	//fprintf(gp,"plot '/home/emerson/Desktop/rtk_simul/Lanes1_XYZ_2016_ITRF08.txt' u 3:4 with points pt 2 ps 0.2 lc rgb \"red\" title \" Ref. Trajectory \" \n ");

/*  fprintf(gp,"plot '../data/kin_Rtklib_enu.txt' using 1:2 w lp ls 8 ps 0.5 lc \"blue \" title \" PPK position \" ,\
                  '../data/Lanes_enu_rtklib1.txt' u 1:2 with lines lt rgb \"black \" title \" Lanes \",\
                   '../data/Lanes_enu_rtklib3.txt' u 1:2 with lines lt rgb \"red \" title \" Lanes corrected \" \n");  */

	//fprintf(gp,"plot '-' with points pt 2 ps 0.1 lc rgb \"red \" title \"Ref. Trajectory \", \n");
/*  fprintf(gp,"plot '../out/exp1_enh_space_qntts.txt' using 1:2 w lp ls 8 ps 0.5 lc \"blue \" title \" SPP position \" , \
                  '../data/exp1_Lanes_enu.txt' u 1:2 with lines lt rgb \"red \" title \" Reference trajectory \",\
                  '../out/exp1_enhancedspace.txt' u 1:2 with lines lt rgb \"green \" title \" Search space \",\
                  '../out/exp1_mmpoint.txt' using 1:2 w p ps 0.15 lc \"black \" title \" Closest point \"  \n");
                //  '../out/exp1_curvecorrep.txt' u 1:2 with lines lt rgb \"blue\" title \" Curve respective points \" \n");
                  //'../out/mmpoint (copy).txt' using 1:2 w p ps 0.15 lc \"black \" title \" Closest point \"  \n");*/
 //'../out/enh_space_qntts3.txt' u 1:2:3:4:5 with ellipses lc \"blue \" title \" 95 conf. ellipse \",  with points pointsize 0.1
  fflush(gp);
 	fclose(gp);
  }

/* Plot height series                                                         */
  extern void heightplot (){//double t, obsd_t *obs, int n
    FILE * gp = popen ("gnuplot -persistent", "w");
    double h1,h2;

  	/* Gnuplot plot commands   - Receiver trajectory  */
  	fprintf(gp, "set xlabel 'Time [UTC]' offset 0,-1,0 font 'Arial,19' \n"
                "set ylabel 'North off-track distance [m]' offset -2,0,0 font 'Arial,19' \n"
                "set tics font 'Arial,19' \n"
                "set key font 'Arial,19' \n"
                "set bmargin at screen 0.15 \n"
                //"set title 'UP COMPONENT COMPARISON' \n"
                //"set title 'RTK/GROUND TRUTH ENU SEPARATION (RAPID-STATIC)' \n"
              //  "set datafile separator ',' \n"
  		//"set term postscript eps enhanced color\n"
  		//"set output '%s ENU_sep_rapidstatic.eps' \n"
      "set style line 2 lc rgb 'blue' lt 6 lw 0.5 ps 0.1 pt 8\n"
      //"set xrange [300:650] \n"
      "set xrange [13.9903:14.1411] \n"
      "set xrange [16.38:16.915] \n"  //static end:16:38:16.516 - obs 16.79:16.866
      //"set yrange [53.5:55.5] \n"
      //"set yrange [40:80] \n"
      "set yrange [-0.05:0.7] \n"
      //"set yrange [52.0:60.5] \n"
      "set key top center\n"
  		//"set autoscale \n"
      "set grid ytics mytics \n"  //draw lines for each ytics and mytics
      //"set mytics 3 \n"     // set the spacing for the mytics
  		"set grid \n",outpath[0]);

      /* Height convergence plots             // Time management COnv: $1<16.516 ($1):( $1<16.516? ($4): 1/0 )
       fprintf(gp,"plot 1/0 with points pointsize 1.1 lt 2 title \" RTK \" ,\
                        '../out/PPK_enu.txt' u ($1):( $1<16.79? ($4): 1/0 ) with points pointsize 0.5 lt 2 notitle,\
                        1/0 with points pointsize 1.1 lt 3 title \" PPP-only \" ,\
                        '../out/PPP_enu.txt' u ($1):( $1<16.79? ($4): 1/0 ) with points pointsize 0.5 lt 3 notitle ,\
                        1/0 with points pointsize 1.1 lt 4 title \" PPP-mod \" ,\
                        '../out/PPP_cand_enu.txt' u ($1):( $1<16.79? ($4): 1/0 ) with points pointsize 0.5 lt 4 notitle \n");
*/

      /* 2D distances to the Ground-truth
      fprintf(gp,"plot 1/0 with points pointsize 1.1 lt 4 title \" PPP-mod \" ,\
                 '../out/PPPmod2Ddist.txt' u 1:3 with points pointsize 0.7 lt 4 notitle,\
                 1/0 with points pointsize 1.1 lt 2 title \" RTK \" ,\
                 '../out/PPK2Ddist.txt' u 1:3 with points pointsize 0.7 lt 2 notitle,\
                 1/0 with points pointsize 1.1 lt 3 title \" PPP-only \" ,\
                 '../out/PPPonly2Ddist.txt' u 1:3 with points pointsize 0.7 lt 3 notitle \n");
   */

  // fprintf(gp,"plot '../data/PPK_CLST_enu_diff2.txt' u 1:(sqrt( ($2)*($2)+($3)*($3)+($4)*($4) ) ) with points pointsize 0.5 lt 2 title \" 3D dist \" \n");


  /* fprintf(gp,"plot '../data/PPK_CLST_enu_diff5.txt' u 1:2 with points pointsize 0.5 lt 2 title \" E \" ,\
                    '../data/PPK_CLST_enu_diff5.txt' u 1:3 with points pointsize 0.5 lt 3 title \" N \" ,\
                    '../data/PPK_CLST_enu_diff5.txt' u 1:4 with points pointsize 0.5 lt 4 title \" UP \" \n");*/

/*   fprintf(gp,"plot '../out/clst_lane_heights_rtklib1.txt' u 1:($2) with lines ls 0.3 lt rgb 'blue' title \"Ground truth (at rover's height)\" ,\
    '../out/clst_lane_heights_rtklib4.txt' u 1:($2) with lines ls 0.3 lt rgb 'red' title \"Ground truth (on the ground) \" ,\
                     '../data/Rtklib_kin_heightcorr.txt' u 1:($2) with points pointsize 0.6 lt 4 title \"RTK - Rtklib \",\
                     '../data/Rtklib_PPP_estTrop_kin_height.txt' with points pointsize 0.6 lt 2 title \"PPP (Trop est.) - Rtklib \" ,\
                     '../data/Rtklib_PPP_kin_height.txt' u 1:($2) with points pointsize 0.6 lt 3 title \"PPP (Trop. model Saast.)- Rtklib \" ,\
                     '../data/NrCan_kin_height1.txt' u 1:($2) with points pointsize 0.6 lt 5 title \"PPP - NRCan \" ,\
                     '../data/Rtklib_PPP_kin_height.txt' u 1:($2) with points pointsize 0.6 lt 1 title \"PPP (Trop. est.) - GAPS \" ,\
                     '../data/MgT_kin_height.txt' u 1:($2) with points pointsize 0.6 lt 6 title \"RTK (on the ground) - Magnet Tools \" ,\
                     '../data/exp1_PPP_back_height.pos' u 1:($2) with points pointsize 0.6 lt 7 title \"PPP (backward) - Rtklib \" \n"); */


    fflush(gp);
   	fclose(gp);
}

/* Plot corrdinate conmparisons series                                                         */
  extern void coordcompplot (){//double t, obsd_t *obs, int n
    FILE * gp = popen ("gnuplot -persistent", "w");

  	/* Gnuplot plot commands   - Receiver trajectory  */
  	fprintf(gp,
      "set multiplot layout 3,1 title 'X,Y, and Z coordinate comparison' \n"
      "set xlabel 'TIME(UTC)' \n"
  		"set autoscale \n"
  		"set grid \n",outpath[0]);

      fprintf(gp, "set xrange [16.55:16.883] \n");
      fprintf(gp, "set yrange [1761520:1761650] \n");
    fprintf(gp,
      "plot '../out/3dcand.txt' u 1:2 with lines title \" map-match cand.\" ,\
      '../out/spp.txt' u 1:2 with lines title \" SPP \",\
      '../out/finalpos.txt' u 1:2 with lines title \" PPP \" \n");

    //  fflush(gp);
fprintf(gp, "set yrange [-4078850:-4078730] \n");
      fprintf(gp,
        "plot '../out/3dcand.txt' u 1:3 with lines title \" map-match cand.\" ,\
        '../out/spp.txt' u 1:3 with lines title \" SPP \",\
        '../out/finalpos.txt' u 1:3 with lines title \" PPP \" \n");

      //  fflush(gp);
fprintf(gp, "set yrange [4560800:4560900] \n");
        fprintf(gp,
          "plot '../out/3dcand.txt' u 1:4 with lines title \" map-match cand.\" ,\
          '../out/spp.txt' u 1:4 with lines title \" SPP \",\
          '../out/finalpos.txt' u 1:4 with lines title \" PPP \" \n");

          fprintf(gp, "unset multiplot\n");
          fprintf(gp, "set term postscript eps enhanced color \n"
            		"set output '%s cand_spp_final_pos_analysis.ps' \n");

    fflush(gp);
   	fclose(gp);
}

/* Receive map-matching quantities for plot  --------------------------*/
extern void coordcompcalcplot (){//double t, obsd_t *obs, int n
	int i=0, j=0, k=0, sizef;
  FILE * gp = popen ("gnuplot -persistent", "w");
	FILE *f1, *f2, *f3;
  f1 = fopen("../out/3dcand.txt", "r");
	f2 = fopen("../out/spp.txt", "r");
	f3 = fopen("../out/ppk_llh.txt", "r");

	double lat, lon,h, x, y, z, pos[3], pos1[3], pos2[3], r[3], pos0[3], enu[3];
  double t,t1,t2,r1[3],r2[3],enu1[3],enu2[3];
	char str[150];
  double lat0=45.943134902*D2R, lon0=-66.641067859*D2R, h0=52.5909;
  pos0[0]=lat0;pos0[1]=lon0;pos0[2]=h0;

	/* Gnuplot plot commands   - Receiver trajectory  */
	fprintf(gp,
  "set multiplot layout 3,1 title 'PPK/SPP/MAP local e,n,u time series comparison' \n"
  "set xlabel 'TIME(UTC)' \n"
  "set autoscale \n"
  "set grid \n");

  fprintf(gp, "set yrange [-100:100] \n");
  fprintf(gp, "set autoscale \n");
  fprintf(gp, "set xrange [16.550:16.883] \n");

  fprintf(gp,"plot '-' with lines title \" PPK east \" ,\
     '-' with lines title \" MAP east \",\
     '-' with lines title \" SPP east \" \n");

  while ( fgets(str, 150, f3)!= NULL ){
   sscanf(str, "%lf  %lf %lf %lf", &t, &pos[0], &pos[1], &pos[2]);
   pos2ecef(pos,r);
   ecef2enu(pos0, r, enu);
   fprintf(gp, "%lf %lf\n", t, enu[0]);
  }
  fflush(gp);
  fprintf(gp, "e\n");
  while ( fgets(str, 150, f1)!= NULL ){
     sscanf(str, "%lf %lf %lf %lf", &t1, &r1[0], &r1[1], &r1[2]);
     ecef2enu(pos0, r1, enu1);
     fprintf(gp, "%lf %lf\n", t1, enu1[0]);
   }
   fflush(gp);
   fprintf(gp, "e\n");
   while ( fgets(str, 150, f2)!= NULL ){
      sscanf(str, "%lf %lf %lf %lf", &t2, &r2[0], &r2[1], &r2[2]);
      ecef2enu(pos0, r2, enu2);
      fprintf(gp, "%lf %lf\n", t2, enu2[0]);
    }
   fflush(gp);
   fprintf(gp, "e\n");
  rewind(f1);rewind(f2);rewind(f3);

     /* Second plot*/
     fprintf(gp, "set yrange [-21500:-21350] \n");
     fprintf(gp, "set autoscale \n");
     fprintf(gp, "set xrange [16.550:16.883] \n");

     fprintf(gp,"plot '-' with lines title \" PPK north \" ,\
        '-' with lines title \" MAP north \",\
        '-' with lines title \" SPP north \" \n");

    while ( fgets(str, 150, f3)!= NULL ){
      sscanf(str, "%lf  %lf %lf %lf", &t, &pos[0], &pos[1], &pos[2]);
      pos2ecef(pos,r);
      ecef2enu(pos0, r, enu);
      fprintf(gp, "%lf %lf\n", t, enu[1]);
     }
     fflush(gp);
     fprintf(gp, "e\n");
     while ( fgets(str, 150, f1)!= NULL ){
        sscanf(str, "%lf %lf %lf %lf", &t1, &r1[0], &r1[1], &r1[2]);
        ecef2enu(pos0, r1, enu1);
        fprintf(gp, "%lf %lf\n", t1, enu1[1]);
      }
      fflush(gp);
      fprintf(gp, "e\n");
      while ( fgets(str, 150, f2)!= NULL ){
         sscanf(str, "%lf %lf %lf %lf", &t2, &r2[0], &r2[1], &r2[2]);
         ecef2enu(pos0, r2, enu2);
         fprintf(gp, "%lf %lf\n", t2, enu2[1]);
       }
      fflush(gp);
      fprintf(gp, "e\n");
  rewind(f1);rewind(f2);rewind(f3);

  /* Third plot */
  fprintf(gp, "set yrange [6367150:6367160] \n");
  fprintf(gp, "set autoscale\n");
  fprintf(gp, "set xrange [16.550:16.883] \n");

  fprintf(gp,"plot '-' with lines title \" PPK up \" ,\
     '-' with lines title \" MAP up \",\
     '-' with lines title \" SPP up \" \n");

 while ( fgets(str, 150, f3)!= NULL ){
   sscanf(str, "%lf  %lf %lf %lf", &t, &pos[0], &pos[1], &pos[2]);
   pos2ecef(pos,r);
   ecef2enu(pos0, r, enu);
   fprintf(gp, "%lf %lf\n", t, enu[2]);
  }
  fflush(gp);
  fprintf(gp, "e\n");
  while ( fgets(str, 150, f1)!= NULL ){
     sscanf(str, "%lf %lf %lf %lf", &t1, &r1[0], &r1[1], &r1[2]);
     ecef2enu(pos0, r1, enu1);
     fprintf(gp, "%lf %lf\n", t1, enu1[2]);
   }
   fflush(gp);
   fprintf(gp, "e\n");
   while ( fgets(str, 150, f2)!= NULL ){
      sscanf(str, "%lf %lf %lf %lf", &t2, &r2[0], &r2[1], &r2[2]);
      ecef2enu(pos0, r2, enu2);
      fprintf(gp, "%lf %lf\n", t2, enu2[2]);
    }
	//fprintf(gp, "%g %lf \n", t, obs[1].L[1]);
  fprintf(gp, "unset multiplot \n");
	fflush(gp);
	fclose(gp);
	fclose(f1);	fclose(f2);	fclose(f3);
}


/* Residuals visualization from RTKLIB kalman filter --------------------------------*/

/* Receive quantities for plot -----------------------------------------------*/
void satresplot (satres* sat, int n, int nprn, int *PRN, float ti, float tf){
 int i=0, j=0;
 char prn[nprn];
 FILE * gp = popen ("gnuplot -persistent", "w");
 FILE * gp1 = popen ("gnuplot -persistent", "w");

 /* Gnuplot plot commands   - Pseudorange  */
 fprintf(gp,  "set xlabel 'TOW (sec)' \n"
              "set ylabel 'PSEUDORANGE RESIDUAL (m)' \n"
	            "set autoscale \n"
              "set key horiz \n"
	            "set key bot center \n"
              "set xrange [%g:%g] \n"
              //"set term postscript eps enhanced color\n"
	            //"set output 'Pseudorange_res.ps'\n"
	            "set grid \n", ti, tf);
  /* Gnuplot plot commands   - Carrier phase  */
 fprintf(gp1, "set xlabel 'TOW (sec)' \n"
              "set ylabel 'CARRIER PHASE RESIDUAL (m)' \n"
	            "set autoscale \n"
              "set key horiz \n"
	            "set key bot center \n"
              "set xrange [%g:%g] \n"
              //"set term postscript eps enhanced color\n"
	            //"set output 'Carrier_res.ps'\n"
	            "set grid \n", ti, tf);
                                              //  FIX PLOT TITTLES FOR PRNs****
 fprintf(gp,"plot for[i=1:%d] '-' with points ps 0.3 title \" PRN \".i  \n", nprn);
 fprintf(gp1,"plot for[i=1:%d] '-' with points ps 0.3 title \" PRN \".i  \n", nprn);

 for (i=0;i<nprn;i++){
  for(j=0; j<n; j++)
        {
             fprintf(gp,"%f %f\n", sat[PRN[i]-1].ep[j], sat[PRN[i]-1].P[j]);
	     fprintf(gp1,"%f %f\n", sat[PRN[i]-1].ep[j], sat[PRN[i]-1].L[j]);
  }
  fflush(gp);
  fprintf(gp, "e\n");
  fflush(gp1);
  fprintf(gp1, "e\n");
 }
 fflush(gp);
 fclose(gp);
 fflush(gp1);
 fclose(gp1);
}

void clkplot(clkres clkr, int n, float ti, float tf){
 int i;
 FILE * gp = popen ("gnuplot -persistent", "w");

 /* Gnuplot plot commands   - Clock  */
 fprintf(gp, "set xlabel 'TOW (sec)' \n"
             "set ylabel 'CLOCK (ns)' \n"
	     "set autoscale \n"
             "set key horiz \n"
	     "set key bot center \n"
             "set xrange [%g:%g] \n"
             "set term postscript eps enhanced color\n"
	     "set output 'Recclpck_states.ps'\n"
	     "set grid \n", ti, tf);


 fprintf(gp,"plot '-' with lines ls 2 title \" Receiver clock \" \n");

 for(i=0; i<n; i++){
   fprintf(gp,"%f %f\n", clkr.ep[i], clkr.clk1[i]);
 }
 fflush(gp);
 fclose(gp);
}

void tropplot(tropres tropo, int n, float ti, float tf){
int i;
 FILE * gp = popen ("gnuplot -persistent", "w");

 /* Gnuplot plot commands   - Troposphere  */
 fprintf(gp, "set xlabel 'TOW (sec)' \n"
             "set ylabel 'ZTD (m)' \n"
	     "set autoscale \n"
             "set key horiz \n"
	     "set key bot center \n"
             "set xrange [%g:%g] \n"
             "set term postscript eps enhanced color\n"
	     "set output 'Tropo_states.ps'\n"
	     "set grid \n", ti, tf);


 fprintf(gp,"plot '-' w l title \" Zenith tropospheric delay (float) \", \
                  '-' w l title \" Zenith tropospheric delay (fix) \" \n");

 for(i=0; i<n; i++){
   fprintf(gp,"%f %f\n", tropo.ep[i], tropo.ztd[i]);
 }
 fflush(gp);
 fprintf(gp, "e\n");
 for(i=0; i<n; i++){
   fprintf(gp,"%f %f\n", tropo.ep[i], tropo.ztdf[i]);
 }

 fflush(gp);
 fclose(gp);
}

/* Read stat file and print results --------------------------------------------
* read stat file from rtklib
* args   : FILE res	I	file pointer to rtklib stat file
* return : none
*-----------------------------------------------------------------------------*/
extern void resprint(FILE* res){
 float ti,tf, t, P, L, clk, ztd, ztdf;
 int i,j,n=0, prn, tcur, PRN[MAXSAT], nprn;
 char buff[150], marksat[]="$SAT", markclk[]="$CLK";
 char marktrop[]="$TROP", str;
 satres satres[MAXSAT];
 clkres clkr={};
 tropres tropo={};

/* Get initial and final times	*/
 fgets(buff, 150, res);
 sscanf(buff, "%*4s,%*4d,%f", &ti);
 while ( fgets(buff, 150, res)!= NULL ){ /*puts(buff);*/ }
 sscanf(buff, "%*4s,%*4d,%f", &tf);
 //rewind(fp_lane);

 n = (int)(tf-ti)/1; /* Number of epochs tf-ti/t, t is the obs. rate	*/

 /* Declaring satellite residuals structure with proper sizes	*/
 for (i=0;i<MAXSAT;i++){
 satres[i].ep = (float *)calloc(n,  sizeof(float));
 satres[i].P = (float *)calloc(n,  sizeof(float));
 satres[i].L = (float *)calloc(n,  sizeof(float));
 }
 /* Allocating memory for clk and trop structures	*/
 clkr.clk1 = (float *)calloc(n,  sizeof(float));
 clkr.ep = (float *)calloc(n,  sizeof(float));
 tropo.ztd = (float *)calloc(n,  sizeof(float));
 tropo.ztdf = (float *)calloc(n,  sizeof(float));
 tropo.ep = (float *)calloc(n,  sizeof(float));

 /* Run through the file looking for the marks 	*/
 rewind(res);
 while ( fgets(buff, 150, res)!= NULL ){
 //	puts(buff);
 	if(strstr(buff,marksat)){ /* Sat res. */
         sscanf(buff, "%*4s,%*4d,%f,%*1s%2d,%*d,%*f,%*f,%f,%f", &t, &prn,&P,&L);
         tcur = (int)(t-ti)/30;
	 satres[prn-1].ep[tcur]=t;
	 satres[prn-1].P[tcur]=P;
	 satres[prn-1].L[tcur]=L;
	}else if( strstr(buff,markclk) ){ /* Clock */
		sscanf(buff, "%*4s,%*4d,%f,%*d,%*d,%f", &t, &clk);
         	tcur = (int)(t-ti)/30;
	 	clkr.clk1[tcur]= clk;
	 	clkr.ep[tcur]= t;
	      } else if( strstr(buff,marktrop) ){ /* Tropo */
			sscanf(buff, "%*5s,%*4d,%f,%*d,%*d,%f,%f", &t, &ztd,&ztdf);
         		tcur = (int)(t-ti)/30;
	 		tropo.ztd[tcur] = ztd;
 			tropo.ztdf[tcur] = ztdf;
 			tropo.ep[tcur] = t;
		 }
 }

 /* Counting PRNs and assigning to an array for better plot control	*/
 nprn=0;
 j=0;
 for (i=0;i<MAXSAT;i++){
  for (j=0;j<n;j++){
   if (satres[i].P[j] != 0){
    PRN[nprn]=i+1;
    nprn++;
    j=n;
   }
  }
 }

 /* Calling plot functions	*/
 //satresplot(satres, n, nprn, PRN, ti, tf);
 //clkplot(clkr, n, ti, tf);
 //tropplot(tropo, n, ti, tf);

 /* Free memory	*/
 for (i=0;i<MAXSAT;i++){
  free(satres[i].ep);  free(satres[i].P);  free(satres[i].L);
 }
 free(clkr.clk1); free(clkr.ep);
 free(tropo.ztd); free(tropo.ztdf); free(tropo.ep);

 return;
}


/* Read stat file and print results --------------------------------------------
* read stat file from rtklib
* args   : FILE res	I	file pointer to rtklib stat file
* return : none
*-----------------------------------------------------------------------------*/
extern void resprint1(){
  FILE * gp = popen ("gnuplot -persistent", "w");
  //FILE * gp1 = popen ("gnuplot -persistent", "w");

              /* Gnuplot plot commands   - Pseudorange  */
              fprintf(gp,  "set datafile separator ',' \n"
                           "set xlabel 'TOW (sec)' \n"
                           "set ylabel 'PSEUDORANGE RESIDUAL (m)' \n"
                          // "set term postscript eps enhanced color\n"
                          // "set output '%s mapmatch_analysis_test3.ps' \n"
             	            "set autoscale \n"
                           "set key horiz \n"
             	            "set key bot center \n"
                           //"set term postscript eps enhanced color\n"
             	            //"set output 'Pseudorange_res.ps'\n"
             	            "set grid \n",outpath[0]);
               /* Gnuplot plot commands   - Carrier phase
              fprintf(gp1, "set xlabel 'TOW (sec)' \n"
                           "set ylabel 'CARRIER PHASE RESIDUAL (m)' \n"
             	            "set autoscale \n"
                           "set key horiz \n"
             	            "set key bot center \n"
                           "set xrange [%g:%g] \n"
                           //"set term postscript eps enhanced color\n"
             	            //"set output 'Carrier_res.ps'\n"
             	            "set grid \n", ti, tf);*/
                                                           //  FIX PLOT TITTLES FOR PRNs****
/* fprintf(gp,"plot for[i=1:%d] '-' with points ps 0.3 title \" PRN \".i  \n", nprn);
 //  fprintf(gp1,"plot for[i=1:%d] '-' with points ps 0.3 title \" PRN \".i  \n", nprn);

              for (i=0;i<nprn;i++){
               for(j=0; j<n; j++)
                     {
                          fprintf(gp,"%f %f\n", sat[PRN[i]-1].ep[j], sat[PRN[i]-1].P[j]);
             	     fprintf(gp1,"%f %f\n", sat[PRN[i]-1].ep[j], sat[PRN[i]-1].L[j]);
               }
               fflush(gp);
               fprintf(gp, "e\n");
               //fflush(gp1);
            //   fprintf(gp1, "e\n");
              }
              fflush(gp);
              fclose(gp);
              fflush(gp1);
              fclose(gp1);*/

  fprintf(gp,"plot '../out/PPP_kine.pos.stat' using 3:( strcol(1) eq '$SAT' &&strcol(4) eq 'G10' ? $8 : 1/0 ) with lines title strcol(4) \n");
//( 'stringcolumn(1)' eq '$SAT' ? 8 : 1/0 )
  fflush(gp);
   fclose(gp);

}

extern void resprint2 (){
 FILE *fp;
 fp = fopen("../out/Rescomp_nrcan.txt","r");
 char str[150];
 double ti,tf,t,taux=0.0,vl1,vl2,vp1,vp2;
 int i,j,sat=0,aux=0,count=0,eq,dif;
 int prn[MAXSAT]={0};

/* Getting first sat to count the rest */
 fgets(str, 150, fp);
 sscanf(str, "%lf %d", &t, &sat);
 prn[0]=sat;
 count=1;//one value in the vector position count-1
 ti=t;
/* Counting the number of satellite from file and storing at prn[] vector */
 while ( fgets(str, 150, fp)!= NULL ){
  sscanf(str, "%lf %d", &t, &sat);

    eq=0;dif=0;
    for (i = 0; i < count; i++) {
      if (sat == prn[i]) {
        eq++;
        break;
      }else dif++;
   }
   if (i==count) { //Add a new one
     prn[count]=sat;
     count++;
   }
 }
 /* Getting final time */
 tf=t;
 rewind(fp);
ti=16.39;tf=16.883;

for (i = 0; i < count; i++) {
  FILE * gp = popen ("gnuplot -persistent", "w");
  FILE * gp1 = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - Pseudorange */
  fprintf(gp,  "set xlabel 'TOW (sec)' \n"
               "set ylabel 'PSEUDORANGE RESIDUAL (m)' \n"
               "set title 'PSEUDORANGE RESIDUAL PRN%d' \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key bot center \n"
               "set xrange [%g:%g] \n"
               //"set term postscript eps enhanced color\n"
 	            //"set output 'Pseudorange_res.ps'\n"
 	            "set grid \n", prn[i],ti, tf);
   /* Gnuplot plot commands   - Carrier phase */
  fprintf(gp1, "set xlabel 'TOW (sec)' \n"
               "set ylabel 'CARRIER PHASE RESIDUAL (m)' \n"
               "set title 'CARRIER PHASE RESIDUAL PRN%d' \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key bot center \n"
               "set xrange [%g:%g] \n"
               //"set term postscript eps enhanced color\n"
 	            //"set output 'Carrier_res.ps'\n"
 	            "set grid \n", prn[i],ti, tf);

fprintf(gp,"plot '../out/Rescomp_nrcan.txt'using 1:( $2==(%d)? $4 : 1/0 ) with points ps 0.3 title \" SPP miscl. \" , \
'../out/Rescomp_nrcan.txt' using 1:( $2==(%d)? $6 : 1/0 ) with points ps 0.3 title \" Candidate miscl.\" \n",prn[i],prn[i]);

//fprintf(gp,"plot '../out/Rescomp.txt' using 1:( $2==(%d)? $6 : 1/0 ) with points ps 0.3 title \" Candidate miscl.\" \n",prn[i]);

fprintf(gp1,"plot '../out/Rescomp_nrcan.txt' using 1:( $2==(%d)? $3 : 1/0 ) with points ps 0.3 title \" SPP miscl. \" , \
 '../out/Rescomp_nrcan.txt' using 1:( $2==(%d)? $5 : 1/0 ) with points ps 0.3 title \" Candidate miscl.\" \n",prn[i],prn[i]);

  fflush(gp);
  fclose(gp);
  fflush(gp1);
  fclose(gp1);
}

 fclose(fp);
}

extern void resprint3 (){
 FILE *fp;
 fp = fopen("../out/Resppos_2iter.txt","r");
 char str[150];
 double ti,tf,t,taux=0.0,vl1,vl2,vp1,vp2;
 int i,j,sat=0,aux=0,count=0,countcand=0,eq,dif,cdd;
 int prn[MAXSAT]={0},cand[BUFFSIZE*2];

/* Getting first sat to count the rest */
 fgets(str, 150, fp);
 sscanf(str, "%lf %*d %2d", &t, &sat);
 prn[0]=sat;
 count=1;//one value in the vector position count-1
 ti=t;
/* Counting the number of satellite from file and storing at prn[] vector */
 while ( fgets(str, 150, fp)!= NULL ){
  sscanf(str, "%*lf %*d %2d", &sat);

    eq=0;dif=0;
    for (i = 0; i < count; i++) {
      if (sat == prn[i]) {
        eq++;
        break;
      }else dif++;
   }
   if (i==count) { //Add a new one
     prn[count]=sat;
     count++;
   }
 }
 /* Getting final time */
 tf=t;
 rewind(fp);
//ti=16.39;tf=16.883;
//printf("%lf,%lf\n", ti,tf);

/* Counting the number of candidates  */
rewind(fp);
cand[0]=0;
t=ti;
while (t==ti){
  fgets(str, 150, fp);
   sscanf(str, "%lf %d %*1s%2d", &t, &cdd);
     eq=0;dif=0;
     for (i = 0; i < countcand; i++) {
       if (cdd == cand[i]) {
         eq++;
         break;
       }else dif++;
    }
    if (i==countcand) { //Add a new one
      cand[countcand]=cdd;
      countcand++;
    }
}
printf(" # of Sat: %d \n", count);
for (i = 0; i < count; i++) {
  printf("%d,", prn[i]);
}
printf(" \n# of Cand: %d \n", countcand);
for (i = 0; i < countcand; i++) {
  printf("%d,", cand[i]);
}


t==ti;
t=16.803500;
  for (i = 0; i < count; i++) { //PRN loop
    FILE * gp = popen ("gnuplot -persistent", "w");
    FILE * gp1 = popen ("gnuplot -persistent", "w");

    /* Gnuplot plot commands   - Pseudorange */
    fprintf(gp,  "set xlabel 'CANDIDATES'' \n"
                 "set ylabel '(m)' \n"
                 "set title 'G%d PSEUDORANGE RESIDUALS AT t=%g' \n"
   	            "set autoscale \n"
                 "set key horiz \n"
   	            "set key top center \n"
                // "set xrange [%d:%d] \n"
                 //"set yrange [-1:1] \n"
                 //"set term postscript eps enhanced color\n"
   	            //"set output 'Pseudorange_res.ps'\n"
   	            "set grid \n", prn[i],t,cand[0]-2, cand[countcand-1]+8);
     /* Gnuplot plot commands   - Carrier phase */
    fprintf(gp1, "set xlabel 'CANDIDATES' \n"
                 "set ylabel '(m)' \n"
                 "set title 'G%d CARRIER PHASE RESIDUALS AT t=%g' \n"
   	            "set autoscale \n"
                 "set key horiz \n"
   	            "set key top center \n"
                // "set xrange [%d:%d] \n"
                 //"set yrange [-0.2:0.2] \n"
                 //"set term postscript eps enhanced color\n"
   	            //"set output 'Carrier_res.ps'\n"
   	            "set grid \n", prn[i],t,cand[0]-2, cand[countcand-1]+8);

    //while (t==ti){

      fprintf(gp,"plot '../out/Resppos_2iter.txt' \
      using ($2):( $1==(%g)&&$3==(%d) ? $4 : 1/0 ) with points pointsize 1.3 lt rgb 'blue' \
      title \" Candidates residuals \" \n",t,prn[i]);

      fprintf(gp1,"plot '../out/Resppos_2iter.txt' \
      using ($2):( $1==(%g)&&$3==(%d)? $5 : 1/0 ) with points pointsize 1.3 lt rgb 'red' \
      title \" Candidates residuals \"  \n",t,prn[i]);

     fflush(gp);
     fclose(gp);
     fflush(gp1);
     fclose(gp1);
//    }

  }


 fclose(fp);
}
/* Plot PPP and Best Candidates vtpv value for each epoch */
extern void resprint4 (){
 double ti,tf;
ti=16.39;tf=16.866;

  FILE * gp = popen ("gnuplot -persistent", "w");
  FILE * gp1 = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - Pseudorange */
  fprintf(gp,  "set xlabel 'Hours' \n"
               "set terminal postscript enhanced \n"
               "set ylabel '{/Symbol c}^2_{k}' \n"
               "set title 'PSEUDORANGE NORMALIZED RESIDUALS {/Symbol c}^2_{P3}'  \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
               "set xrange [%g:%g] \n"
            //   "set yrange [-0.005:0.07] \n"
               "set yrange [-0.001:2] \n"
               "set term postscript eps enhanced color\n"
 	             "set output '../out/Pseudorange4_normres_ppp_cand_comp.ps'\n"
 	            "set grid \n", ti, tf);
   /* Gnuplot plot commands   - Carrier phase */
  fprintf(gp1, "set xlabel 'Hours' \n"
               "set terminal postscript enhanced \n"
               "set ylabel '{/Symbol c}^2_{k}' \n"
               "set title 'CARRIER-PHASE NORMALIZED RESIDUALS {/Symbol c}^2_{L3}' \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
              "set xrange [%g:%g] \n"
            //  "set yrange [-0.001:0.03] \n"
               "set yrange [-0.001:0.2] \n"
               "set term postscript eps enhanced color\n"
 	            "set output '../out/Carrier_normres_ppp_cand_comp.ps'\n"
 	            "set grid \n", ti, tf);


fprintf(gp,"plot '../out/Resppos3_vtpv_ppp.txt' using 1:3 with lines ls 0.3 lt 5 title \"PPP-only \" ,\
                 '../out/Resppos3_vtpv.txt' using 1:3 with points pointsize 0.3 lt 2 title \"PPP-mod\" \n");

fprintf(gp1,"plot '../out/Resppos3_vtpv_ppp.txt' using 1:2 with lines ls 0.3 lt 5 title \"PPP-only \" ,\
                  '../out/Resppos3_vtpv.txt' using 1:2 with points pointsize 0.3 lt 2 title \"PPP-mod\" \n");

/*
      fprintf(gp,"plot '../out/Resppos3_vtpv_ppp.txt' using 1:3 with lines pointsize 0.3 lt 4 title \"PPP \"  \n");

      fprintf(gp1,"plot '../out/Resppos3_vtpv_ppp.txt' using 1:2 with lines pointsize 0.3 lt 2 title \"PPP \"  \n");
*/
  fflush(gp);
  fclose(gp);
  fflush(gp1);
  fclose(gp1);
}

/* Read RTKLIB position solution file and plot # of satellite                 */
extern void nsat (){//double t, obsd_t *obs, int n
    FILE * gp = popen ("gnuplot -persistent", "w");

  	/* Gnuplot plot commands   - Receiver trajectory  */
  	fprintf(gp, "set xlabel 'Time (hours)' font 'Arial,19' \n"
                //"set ylabel 'Number of satellites' \n"
                //"set title 'SPP STD AND NUMBER OF SATELLITES' \n"
  	//	"set term postscript eps enhanced color\n"
  	//	"set output '%s sppstd_nsat.ps' \n"
    "set tics font 'Arial,19' \n"
    "set key font 'Arial,19' \n"
     "set size ratio 0.6 1\n"
      "set xrange [16.38:16.889] \n"
      "set yrange [0:10] \n"
      "set key top center \n"
  		//"set autoscale \n"
  		"set grid \n");


      fprintf(gp,"plot '../out/PPP_cand_exp1 (copy).pos' u 3:8 with linespoints ls 0.8 lt -1 title \"Number of satellites \" \n");

    fflush(gp);
   	fclose(gp);
}

/* Plot series                                                         */
  extern void sppstdnsat (){//double t, obsd_t *obs, int n
    FILE * gp = popen ("gnuplot -persistent", "w");
    double h1,h2;

  	/* Gnuplot plot commands   - Receiver trajectory  */
  	fprintf(gp, "set xlabel 'Time (hours)' \n"
                "set ylabel '[m]' \n"
                "set title 'SPP STD AND NUMBER OF SATELLITES' \n"
  		"set term postscript eps enhanced color\n"
  		"set output '%s sppstd_nsat.ps' \n"
    //  "set xrange [16.39:16.91] \n"
    //  "set x2range [16.39:16.91] \n"
      "set xrange [16.7:16.91] \n"
      "set x2range [16.7:16.91] \n"
      "set ylabel 'std' \n"
      "set yrange [0:15] \n"
      "set y2range [0:20] \n"
      "set y2label '# of satellites' \n"
      "set ytics nomirror \n"
      "set y2tics\n"
      "set tics out\n"
      "set key top center \n"
      //"set autoscale y2\n"
  		//"set autoscale \n"
      //"set grid ytics mytics \n"  //draw lines for each ytics and mytics
      //"set mytics 3 \n"     // set the spacing for the mytics
  		"set grid \n",outpath[0]);


      fprintf(gp,"plot '../out/sppstd_nsat.txt' u 1:2 with lines ls 0.5 lt 1 title \" X-std \", \
                        '../out/sppstd_nsat.txt' u 1:3 with lines ls 0.5 lt 2 title \" Y-std \", \
                        '../out/sppstd_nsat.txt' u 1:4 with lines ls 0.5 lt 3 title \" Z-std \", \
                        '../out/sppstd_nsat.txt' u 1:5 with lines ls 0.5 lt 4 title \" # SAT \" axes x2y2 \n");

    fflush(gp);
   	fclose(gp);
}

/* Plot PPP and Best Candidates vtpv value for each epoch */
extern void resinnovbounds (){
 double ti,tf;
ti=16.39;tf=16.866;

  FILE * gp = popen ("gnuplot -persistent", "w");
  FILE * gp1 = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - Pseudorange */
  fprintf(gp,  "set xlabel 'Time[UTC]' offset 0,-1,0 font 'Arial,19' \n"
               //"set terminal postscript enhanced \n"
               "set ylabel '[m]' font 'Arial,19' \n"
               "set tics font 'Arial,19' \n"
               "set key font 'Arial,19' \n"
               "set bmargin at screen 0.15 \n"
               //"set title 'PSEUDORANGE RESIDUALS v_{k} '  \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
               "set xrange [%g:%g] \n"
              // "set yrange [-0.005:0.07] \n"
               "set yrange [-20:20] \n"
               "set xtics 16.45,0.1,16.85 \n"
               //"set term postscript eps enhanced color\n"
 	             //"set output '../out/Pseudorange_v.ps'\n"
 	            "set grid \n", ti, tf);
   /* Gnuplot plot commands   - Carrier phase */
  fprintf(gp1, "set xlabel 'Time[UTC]' offset 0,-1,0 font 'Arial,19' \n"
               //"set terminal postscript enhanced \n"
               "set ylabel '[m]' font 'Arial,19'\n"
               "set tics font 'Arial,19' \n"
               "set key font 'Arial,19' \n"
               "set bmargin at screen 0.15 \n"
               "set xtics 16.45,0.1,16.85 \n"
               //"set title 'CARRIER-PHASE RESIDUALS v_{k} ' \n"
 	            "set autoscale \n"
               "set key horiz \n"
              // "set xdata time\n"
              // "set timefmt '%H.hh' \n"
               "set xrange [%g:%g] \n"
              // "set format x '%H:%M' \n"
              // "set timefmt '%H:%M'\n"
 	            "set key top center \n"
              //"set xrange [%g:%g] \n"
             // "set yrange [-0.001:0.03] \n"
              "set yrange [-0.1:0.1] \n"
              //"set term postscript eps enhanced color\n"
 	            //"set output '../out/Carrier_vs.ps'\n"
 	            "set grid \n", ti, tf);
/*
fprintf(gp,"plot '../out/Res_pseud_inov_ppp.txt' using 1:2 with lines ls 0.1 lt 3 title \" v_{k} PPP \" ,\
                 '../out/Res_pseud_inov.txt' using 1:2 with points pointsize 0.1 lt 4 title \" v_{k} candidates \" ,\
                 '../out/Res_pseud_inov.txt' using 1:($3/2) with lines ls 0.001 dt 3 title \" +/- 1 {/Symbol s} \", \
                 '../out/Res_pseud_inov.txt' using 1:($4/2) with lines ls 0.001 dt 3 notitle \n");

 fprintf(gp1,"plot '../out/Res_phase_inov_ppp.txt' using 1:2 with lines ls 0.1 lt 5 title \" v_{k} PPP \" ,\
                   '../out/Res_phase_inov.txt' using 1:2 with points pointsize 0.1 lt 2 title \" v_{k} candidates \" ,\
                   '../out/Res_phase_inov.txt' using 1:($3/2) with lines ls 0.001 dt 3 title \" +/- 1 {/Symbol s} \", \
                   '../out/Res_phase_inov.txt' using 1:($4/2) with lines ls 0.001 dt 3 notitle \n");           */

                   fprintf(gp,"plot 1/0 with lines ls 0.5 lt 5 title \"PPP-only \" ,\
                  '../out/Res_pseud_inov_ppp.txt' using 1:2 with lines ls 0.1 lt 5 notitle ,\
                  1/0 with points pointsize 1.1 lt 2 title \"PPP-mod \" ,\
                  '../out/Res_pseud_inov.txt' using 1:2 with points pointsize 0.5 lt 2 notitle \n");

                  fprintf(gp1,"plot 1/0 with lines ls 0.5 lt 5 title \"PPP-only \" ,\
                 '../out/Res_phase_inov_ppp.txt' using 1:2 with lines ls 0.1 lt 5 notitle ,\
                 1/0 with points pointsize 1.1 lt 2 title \"PPP-mod \" ,\
                 '../out/Res_phase_inov.txt' using 1:2 with points pointsize 0.5 lt 2 notitle \n");

  fflush(gp);
  fclose(gp);
  fflush(gp1);
  fclose(gp1);
}

/* Separate residuals */
extern void resinnovqchisq (){
 double ti,tf;
ti=16.39;tf=16.883;

  FILE * gp = popen ("gnuplot -persistent", "w");
  FILE * gp1 = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - Pseudorange */
  fprintf(gp,  "set xlabel 'Epoch_{k} [hours]' \n"
               "set terminal postscript enhanced \n"
               "set ylabel '{/Symbol c}^2' \n"
               "set title 'PSEUDORANGE CHI-SQUARE RESIDUALS TEST'  \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
               "set xrange [%g:%g] \n"
              // "set yrange [-0.005:0.07] \n"
               "set yrange [-50:50] \n"
               "set term postscript eps enhanced color\n"
 	             "set output '../out/Pseudorange_inov_chisqr.ps'\n"
 	            "set grid \n", ti, tf);
   /* Gnuplot plot commands   - Carrier phase */
  fprintf(gp1, "set xlabel 'Epoch_{k} [hours]' \n"
               "set terminal postscript enhanced \n"
               "set ylabel '{/Symbol c}^2' \n"
               "set title 'CARRIER PHASE CHI-SQUARE RESIDUALS TEST' \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
              "set xrange [%g:%g] \n"
            //  "set yrange [-0.001:0.03] \n"
               "set yrange [-3:10] \n"
              "set term postscript eps enhanced color\n"
 	            "set output '../out/Carrier_inov_chisqr.ps'\n"
 	            "set grid \n", ti, tf);

fprintf(gp,"plot '../out/Resppos_vtpv_chisq.txt' using 1:3 with lines ls 0.1 lt 3 title \" v^{T}.R^{-1}.v / (n-u)\" ,\
                 '../out/Resppos_vtpv_chisq.txt' using 1:($4) with lines ls 0.001 dt 3 title \" {/Symbol c}^2_{(n-u,0.05)} \" \n");

 fprintf(gp1,"plot '../out/Resppos_vtpv_chisq.txt' using 1:2 with lines ls 0.1 lt 5 title \" v^{T}.R^{-1}.v / (n-u)  \" ,\
                   '../out/Resppos_vtpv_chisq.txt' using 1:($4) with lines ls 0.001 dt 3 title \" {/Symbol c}^2_{(n-u,0.05)} \" \n");

  fflush(gp);
  fclose(gp);
  fflush(gp1);
  fclose(gp1);
}

/* Separate residuals */
extern void resinnovqchisqtger (){
 double ti,tf;
ti=16.39;tf=16.883;

  FILE * gp = popen ("gnuplot -persistent", "w");


  /* Gnuplot plot commands   - Pseudorange */
  fprintf(gp,  "set xlabel 'Epoch_{k} [hours]' \n"
               "set terminal postscript enhanced \n"
               "set ylabel '{/Symbol c}^2' \n"
               "set title 'CHI-SQUARE RESIDUALS TEST'  \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
               "set xrange [%g:%g] \n"
              // "set yrange [-0.005:0.07] \n"
               "set yrange [-50:50] \n"
               "set term postscript eps enhanced color\n"
 	             "set output '../out/P3L3_inov_chisqr.ps'\n"
 	            "set grid \n", ti, tf);


fprintf(gp,"plot '../out/Resppos_vtpv_chisq_tger.txt' using 1:2 with lines ls 0.1 lt 3 title \" v^{T}.R^{-1}.v / (n-u)\" ,\
                 '../out/Resppos_vtpv_chisq_tger.txt' using 1:3 with lines ls 0.001 dt 3 title \" {/Symbol c}^2_{(n-u,0.05)} \" \n");

  fflush(gp);
  fclose(gp);

}

extern void amb_and_var_plot (){
  FILE * gp = popen ("gnuplot -persistent", "w");

  /* Gnuplot plot commands   - Pseudorange */
  fprintf(gp,  "set xlabel 'Epoch_{k} [hours]' \n"
               "set terminal postscript enhanced \n"
               //"set ylabel 'N_{L3} [m]' \n"
              // "set title 'IONOSPHERE-FREE AMBIGUITY PRN 30'  \n"
               "set ylabel '{/Symbol s}_{N_{L3}} [m]' \n"
               "set title 'IONOSPHERE-FREE AMBIGUITY STD PRN 30'  \n"
 	            "set autoscale \n"
               "set key horiz \n"
 	            "set key top center \n"
              // "set xrange [] \n"
              // "set yrange [-0.005:0.030] \n"
               "set yrange [0:5] \n"
               "set term postscript eps enhanced color\n"
 	             "set output '../out/prn30_ambstd_comp.ps'\n"
 	            "set grid \n");

/*
fprintf(gp,"plot '../out/amb_sol_ppp_only.txt' using 1:( ($2)==30 ? ($3):1/0) with points ps 0.3 title \" PPP-only \" ,\
                 '../out/amb_sol.txt' using 1:( ($2)==30 ? ($3):1/0) with points ps 0.3 title \" PPP-mod \" \n");
*/

fprintf(gp,"plot '../out/amb_sol_ppp_only.txt' using 1:( ($2)==30 ? ($4):1/0) with points pointsize 1 lt 3 title \" PPP-only \" ,\
                   '../out/amb_sol.txt' using 1:( ($2)==30 ? ($4):1/0) with points pointsize 0.7 lt 1 lc rgb \"goldenrod\" title \" PPP-mod \" \n");
  fflush(gp);
  fclose(gp);

}
