#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


int main(int argc, char *argv[])
{
  int           k;
  ostringstream str;

  const double delta = 3e-2;

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

  if (false) get_alphac2();

// GetEmittance(ElemIndex("cav"), true);

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
