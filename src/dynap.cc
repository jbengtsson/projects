#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void err_and_corr_ini(const string &param_file, param_data_type &params,
		      orb_corr_type orb_corr[])
{
  // Set state.
  globval.Cavity_on   = false; globval.radiation = false;
  globval.Aperture_on = false;

  params.get_param(param_file);

  Ring_GetTwiss(true, 0.0); printglob();

  // Store optics function values at sextupoles.
  params.get_bare();

  globval.Cavity_on = false;

  if (params.ae_file != "") {
    if (params.bba) {
      // Beam based alignment
      params.Align_BPMs(Quad);
    }

    cod_ini(params.bpm_Fam_names, params.corr_Fam_names, orb_corr);
  }

  if (params.n_lin > 0) params.ini_skew_cor(params.disp_wave_y);
}


void err_and_corr(const string &param_file)
{
  bool            cav, rad, aper;
  int             j;
  param_data_type params;
  orb_corr_type   orb_corr[2];
  DA_data_type    DA;

  cav = globval.Cavity_on; rad = globval.radiation; aper = globval.Aperture_on;

  err_and_corr_ini(param_file, params, orb_corr);

  if (params.DA_bare) DA.get_DA_bare(params);

  DA.get_DA_real(params, orb_corr);

  if (params.ae_file != "") {
    for (j = 0; j < 2; j++)
      orb_corr[j].dealloc();
  }

  // Reset state.
  globval.Cavity_on = cav; globval.radiation = rad; globval.Aperture_on = aper;
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
