#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double  delta  = 2.5e-2; // delta for off-momentum aperture

const char    home_dir[] = "/home/bengtsson";


void err_and_corr(const char *param_file, const double delta, const int mode)
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

  prt_lat("linlat.out", globval.bpm, true);

  get_bare();

  if (strcmp(ae_file, "") != 0) {
    gcmat(1); gcmat(2); gtcmat(1); gtcmat(2);
  }

  // Coupling correction and dispersion wave initialization.
  if (n_lin > 0) ini_skew_cor(disp_wave_y);

  if (strcmp(fe_file, "") != 0) LoadFieldErr(fe_file, false, 1.0, true);

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

  Ring_GetTwiss(true, 0.0); printglob();

  GetEmittance(ElemIndex("cav"), true);

  prt_beamsizes();

  prt_lat("linlat_err.out", globval.bpm, true);
  prtmfile("flat_file_err.dat");

  if (mode == 1) {
    cout << "fmap" << endl;
    fmap(n_x, n_y, n_tr, x_max_FMA, y_max_FMA, 0.0, true, false);
  } else if (mode == 2) {
    cout << "fmapdp" << endl;
    fmapdp(n_x, n_dp, n_tr, x_max_FMA, -delta_FMA, 0.1e-3, true, false);
  } else {
    cout << "bad param.: mode = " << mode << endl;
    exit(1);
  }

  globval.Cavity_on = cav; globval.radiation = rad; globval.Aperture_on = aper;
}


int main(int argc, char *argv[])
{

  globval.H_exact     = false; globval.quad_fringe = false;
  globval.Cavity_on   = false; globval.radiation   = false;
  globval.emittance   = false; globval.IBS         = false;
  globval.pathlength  = false; globval.Aperture_on = false;

  err_and_corr("param.dat", delta, atoi(argv[1]));
}
