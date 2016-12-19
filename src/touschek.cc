#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void err_and_corr(const string &param_file)
{
  int             j;
  param_data_type params;
  orb_corr_type   orb_corr[2];
  FILE            *fp;

  const double Qb = 5e-9, eps_x = 16e-12, eps_y = 16e-12;
  const double sigma_s = 1e-2, sigma_delta = 1e-3;

  const string file_name = "mom_aper.out";

  params.err_and_corr_init(param_file, orb_corr);

  globval.delta_RF = 10e-2;

  globval.Cavity_on = true;

  Touschek(Qb, globval.delta_RF, eps_x, eps_y, sigma_delta, sigma_s);
      
  double  sum_delta[globval.Cell_nLoc+1][2];
  double  sum2_delta[globval.Cell_nLoc+1][2];
//  double  mean_delta_s[globval.Cell_nLoc+1][2];
//  double  sigma_delta_s[globval.Cell_nLoc+1][2];

  for(j = 0; j <= globval.Cell_nLoc; j++){
    sum_delta[j][X_] = 0.0; sum_delta[j][Y_] = 0.0;
    sum2_delta[j][X_] = 0.0; sum2_delta[j][Y_] = 0.0;
  }
 
  Touschek(Qb, globval.delta_RF, false,
	   eps_x, eps_y, sigma_delta, sigma_s,
	   params.n_track_DA, false, sum_delta, sum2_delta);

  fp = file_write((file_name).c_str()); 
  for(j = 0; j <= globval.Cell_nLoc; j++)
    fprintf(fp, "%4d %7.2f %5.3f %6.3f\n",
	    j, Cell[j].S, 1e2*sum_delta[j][X_], 1e2*sum_delta[j][Y_]);
  fclose(fp);

  params.err_and_corr_exit(orb_corr);
}


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.Aperture_on = false;

  if (argc < 1) {
    printf("*** bad command line\n");
    exit(1);
  }

  err_and_corr(argv[1]);
}
