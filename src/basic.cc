#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


int main(int argc, char *argv[])
{
  int           k;
  double        alpha[2], beta[2], eta[2], etap[2];
  ostringstream str;

  const double delta = 3e-2;

  const int   n_bpm_Fam = 1, n_hcorr_Fam = 1, n_vcorr_Fam = 1;

  const string  bpm_names[n_bpm_Fam]     = { "BPM" };
  const string  hcorr_names[n_hcorr_Fam] = { "CHV" };
  const string  vcorr_names[n_vcorr_Fam] = { "CHV" };

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

//   no_sxt();

  // globval.Cavity_on = true;

  Ring_GetTwiss(true, 0.0); printglob();

  if (true) get_alphac2();

// GetEmittance(ElemIndex("cav"), true);

  if (false) {
    for (k = 0; k < 2; k++) {
      alpha[k] = Cell[0].Alpha[k]; beta[k] = Cell[0].Beta[k];
      eta[k] = 0e0; etap[k] = 0e0;
    }

    ttwiss(alpha, beta, eta, etap, 0.0);
  }

  if (false) {
//     str << home_dir << "/Thor-2.0/thor/wrk/fit_isoch.dat";
//     get_bn2(str.str(), "get_b2.dat", 0, true);

    str.str(""); str << "/home/johan/Thor-2.0/thor/wrk/fit_alpha.dat";
    get_bn2(str.str(), "get_b3.dat", 0, true);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_lat("chromlat.out", globval.bpm, true, 10);

  prtmfile("flat_file.dat");

  if (false) {
    globval.Cavity_on = true; n_track = 512;
    n_aper = 25;
    get_dynap(delta, true);
  }
}
