#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void err_and_corr(const string &param_file, const int mode)
{
  param_data_type params;
  orb_corr_type   orb_corr[2];

  params.err_and_corr_init(param_file, orb_corr);

  if (mode == 1) {
    cout << "fmap" << endl;
    fmap(params.n_x, params.n_y, params.n_tr, params.x_max_FMA,
	 params.y_max_FMA, 0e0, true, false);
  } else if (mode == 2) {
    cout << "fmapdp" << endl;
    fmapdp(params.n_x, params.n_dp, params.n_tr, params.x_max_FMA,
	   -params.delta_FMA, 0.1e-3, true, false);
  } else {
    cout << "bad param.: mode = " << mode << endl;
    exit(1);
  }

  params.err_and_corr_exit(orb_corr);
}


int main(int argc, char *argv[])
{
  globval.H_exact     = false; globval.quad_fringe = false;
  globval.Cavity_on   = false; globval.radiation   = false;
  globval.emittance   = false; globval.IBS         = false;
  globval.pathlength  = false; globval.Aperture_on = false;

  if (argc < 2) {
    printf("*** bad command line\n");
    exit(1);
  }

  err_and_corr(argv[1], atoi(argv[2]));
}
