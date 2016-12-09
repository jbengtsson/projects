#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


// instantiate templates
// template class vector<std::string>;
// template class vector<long int>;


void err_and_corr(const char *param_file,
		  const double delta, const int n_delta)
{
  bool            cod = false, cav, rad, aper;
  char            str[max_str];
  int             j, k;
  double          DA, x_min[2], x_max[2], d[n_delta+1], x_hat;
  double          DA_m[n_delta+1], DA_s[n_delta+1];
  double          x_hat_m[n_delta+1][2], x_hat_s[n_delta+1][2];
  param_data_type params;
  orb_corr_type   orb_corr[2];
  FILE            *DA_real = NULL, *fp[n_delta+1];

  const int    n_cell = 20;
  const double scl    = 0.5;
  
  // Save state.
  cav = globval.Cavity_on; rad = globval.radiation; aper = globval.Aperture_on;

  // Set state.
  globval.Cavity_on   = false; globval.radiation = false;
  globval.Aperture_on = false;

  params.get_param(param_file);

  n_track = 512;

  Ring_GetTwiss(true, 0.0); printglob();

  // Store optics function values at sextupoles.
  params.get_bare();

  if (params.DA_bare) get_DA_bare(delta, n_delta);

  globval.Cavity_on = false;

  if (params.ae_file != "") {
    if (params.bba) {
      // Beam based alignment
      params.Align_BPMs(Quad);
    }

    cod_ini(params.bpm_Fam_names, params.corr_Fam_names, orb_corr);
  }
 
  // Linear coupling correction and dispersion wave initialization.
  if (params.n_lin > 0) params.ini_skew_cor(params.disp_wave_y);

  for (j = 0; j <= n_delta; j++) {
    d[j] = (n_delta > 0)? (double)j/(double)n_delta*delta : 0.0;
    sprintf(str, "DA_real_%4.2f.out", 1e2*d[j]); fp[j] = file_write(str);

    DA_m[j] = 0.0; DA_s[j] = 0.0;

    for (k = 0; k <= 1; k++) {
      x_hat_m[j][k] = 0.0; x_hat_s[j][k] = 0.0;
    }
  }

  DA_real = file_write("DA_real.out");
  fprintf(DA_real, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  Cell[globval.Cell_nLoc].Beta[X_],
	  Cell[globval.Cell_nLoc].Beta[Y_]);
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# Real lattice\n");
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# deltaP        DA               Ax              Ay"
	  "            x^              y^\n");
  fprintf(DA_real, "#   %%          mm^2              mm              mm"
	    "            mm              mm\n");

  for (j = 1; j <= params.n_stat; j++) {
    globval.Cavity_on = false;

    if (params.fe_file != "")
      params.LoadFieldErr(false, 1e0, true);

    if (params.ae_file != "") {
      cod = cod_correct(n_cell, params.n_orbit, scl, j, orb_corr, params);
    }
      
    // cod = CorrectCOD(n_orbit, 0.3);

    printf("\n");
    if (cod) {
      printf("err_and_corr: orbit correction completed\n");
      if (j == 1) prt_cod("cod.out", globval.bpm, true);
 
      Ring_GetTwiss(true, 0.0); printglob();

      GetEmittance(ElemIndex("cav"), true);

      if (params.n_lin > 0) params.corr_eps_y();

      if (params.ap_file != "") params.LoadApers(1.0, 1.0);

      Ring_GetTwiss(true, 0.0); printglob();

      GetEmittance(ElemIndex("cav"), true);

      prt_beamsizes();

      globval.Cavity_on = true;
      for (k = 0; k <= n_delta; k++) {
	DA = get_dynap(fp[k], 10e-3, d[k], n_track, 0.1e-3, n_aper,
		       x_min, x_max);
	DA_m[k] += DA; DA_s[k] += sqr(DA);
	x_hat = (x_max[X_]-x_min[X_])/2.0;
	x_hat_m[k][X_] += x_hat; x_hat_s[k][X_] += sqr(x_hat);
	x_hat = x_max[Y_];
	x_hat_m[k][Y_] += x_hat; x_hat_s[k][Y_] += sqr(x_hat);
      }

     if (params.n_lin > 0)
	// reset skew quads
	// set_bn_design_fam(globval.qt, Quad, 0.0, 0.0);

      if (params.N_calls > 0) params.reset_quads();  
    } else
      chk_cod(cod, "error_and_correction");
  }

  for (j = 0; j <= n_delta; j++) {
    fclose(fp[j]);

    get_mean_sigma(params.n_stat, DA_m[j], DA_s[j]);
    for (k = 0; k <= 1; k++)
      get_mean_sigma(params.n_stat, x_hat_m[j][k], x_hat_s[j][k]);

    fprintf(DA_real, "  %5.2f  %6.1f \xB1 %6.1f"
	    "   %5.1f \xB1 %5.2f    %5.1f \xB1 %5.2f"
	    "   %4.1f \xB1 %5.2f     %4.1f \xB1 %5.2f\n", 
	    d[j]*1e2,
	    1e6*DA_m[j], 1e6*DA_s[j],
	    1e6*sqr(x_hat_m[j][X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat_s[j][X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat_m[j][Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e6*sqr(x_hat_s[j][Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e3*x_hat_m[j][X_], 1e3*x_hat_s[j][X_],
	    1e3*x_hat_m[j][Y_], 1e3*x_hat_s[j][Y_]);
  }

  fclose(DA_real);


  if (params.ae_file != "") {
    for (j = 0; j < 2; j++)
      orb_corr[j].dealloc();
  }

  // Reset state.
  globval.Cavity_on = cav; globval.radiation = rad; globval.Aperture_on = aper;
}


int main(int argc, char *argv[])
{
  int    n_delta;
  double delta;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.Aperture_on = false;

  if (argc < 3) {
    printf("*** bad command line\n");
    exit(1);
  }

  sscanf(argv[2], "%d", &n_delta);
  sscanf(argv[3], "%le", &delta);
  printf("\n");
  printf("n = %d, delta = %3.1f%%\n", n_delta, 1e2*delta);

  err_and_corr(argv[1], delta, n_delta);
}
