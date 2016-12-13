#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double  delta  = 3e-2; // delta for off-momentum aperture


bool    bda = false;
char    mfile_name[max_str];


void err_and_corr(const char *param_file)
{
  bool      cod, cav, rad, aper;
  char      str1[max_str], str2[max_str];
  int       j, k, n, N;
  long int  lastn, lastpos;
  double    m_dbeta[2], s_dbeta[2], m_dnu[2], s_dnu[2], B_10L;

  // save state for dynamics
  cav = globval.Cavity_on; rad = globval.radiation; aper = globval.Aperture_on;

  // set state for correction
  globval.Cavity_on   = false; globval.radiation = false;
  globval.Aperture_on = false;

  get_param(param_file);

  Ring_GetTwiss(true, 0.0); printglob();

  if (strcmp(ae_file, "") != 0) {
    gcmat(1); gcmat(2); gtcmat(1); gtcmat(2);
  }

  // coupling correction and dispersion wave initialization
  if (n_lin > 0) ini_skew_cor(disp_wave_y);

  if (strcmp(ae_file, "") != 0) {
    cod = CorrectCOD_N(ae_file, n_orbit, n_scale, 1);
    printf("\n");
    if (cod)
      printf("orbit correction completed\n");
    else
      chk_cod(cod, "error_and_correction");
  }
      
  GetEmittance(ElemIndex("cav"), true);

  if (n_lin > 0) corr_eps_y();

  if (strcmp(ap_file, "") != 0) LoadApers(ap_file, 1.0, 1.0);

  if (strcmp(fe_file, "") != 0) LoadFieldErr(fe_file, false, 1.0, true);

  Ring_GetTwiss(true, 0.0); printglob();

  GetEmittance(ElemIndex("cav"), true);

  prt_beamsizes();

  prt_lat("linlat_err.out", globval.bpm, true);
  prtmfile("flat_file_err.dat");

  globval.Cavity_on = cav; globval.radiation = rad; globval.Aperture_on = aper;
}


int main(int argc, char *argv[])
{
  int   j;
  FILE  *fp;

  // const double Qb = 1.3e-9, eps_x = 0.85e-9, eps_y = 0.008e-9;
  // const double sigma_s = 4.4e-3, sigma_delta = 0.84e-3;

  // 3 DWs, V_RF = 3.30 MV.
  // const double Qb = 1.3e-9, eps_x = 0.90e-9, eps_y = 0.008e-9;
  // const double sigma_s = 3.79e-3, sigma_delta = 0.76e-3;

  // Commissioning.
  const double Qb = 0.06e-9, eps_x = 2.0e-9, eps_y = 0.01e-9;
  const double sigma_s = 3.4e-3, sigma_delta = 0.52e-3;


  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.Aperture_on = false;

  if (strcmp("n", argv[1]) == 0) {
    n_track = 512;
    printf("\n");
    printf("track for one synchrotron oscillation\n");
  } else if (strcmp("y", argv[1]) == 0) {
    globval.radiation = true;
    n_track = 10000;
    printf("\n");
    printf("track for one damping time\n");
  } else {
    printf("bad command line: %s\n", argv[1]);
    exit(1);
  }

  err_and_corr("param.dat");

  globval.delta_RF = 3.0e-2;

  Touschek(Qb, globval.delta_RF, eps_x, eps_y, sigma_delta, sigma_s);
      
  double  sum_delta[globval.Cell_nLoc+1][2];
  double  sum2_delta[globval.Cell_nLoc+1][2];
//  double  mean_delta_s[globval.Cell_nLoc+1][2];
//  double  sigma_delta_s[globval.Cell_nLoc+1][2];

  // initialize momentum aperture arrays
  for(j = 0; j <= globval.Cell_nLoc; j++){
    sum_delta[j][X_] = 0.0; sum_delta[j][Y_] = 0.0;
    sum2_delta[j][X_] = 0.0; sum2_delta[j][Y_] = 0.0;
  }
 
  Touschek(Qb, globval.delta_RF, false,
	   eps_x, eps_y, sigma_delta, sigma_s,
	   n_track, true, sum_delta, sum2_delta);

  fp = file_write("mom_aper.out"); 
  for(j = 0; j <= globval.Cell_nLoc; j++)
    fprintf(fp, "%4d %7.2f %5.3f %6.3f\n",
	    j, Cell[j].S, 1e2*sum_delta[j][X_], 1e2*sum_delta[j][Y_]);
  fclose(fp);
}
