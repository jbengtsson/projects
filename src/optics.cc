#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


int main(int argc, char *argv[])
{
  int           qf, qd, sf, sd;
  double        b2[2], a2, b3[2], a3;
  ostringstream str;

  const double delta = 3e-2;
  const double nu[] = {101.2/20.0, 27.32/20.0};

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

  GetEmittance(ElemIndex("cav"), true);

  if (false) {
//     str << home_dir << "/Thor-2.0/thor/wrk/fit_isoch.dat";
//     get_bn2(str.str(), "get_b2.dat", 0, true);

    str.str(""); str << "/home/johan/Thor-2.0/thor/wrk/fit_alpha.dat";
    get_bn2(str.str(), "get_b3.dat", 0, true);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (false) {
    qf = ElemIndex("qf"); qd = ElemIndex("bh");
    FitTune(qf, qd, nu[X_], nu[Y_]);
    get_bn_design_elem(qf, 1, Quad, b2[0], a2);
    get_bn_design_elem(qd, 1, Quad, b2[1], a2);

    printf("\nnu_x = %8.5f nu_y = %8.5f\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("  qf = %8.5f  bh = %8.5f\n", b2[0], b2[1]);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (false) {
    sf = ElemIndex("sf"); sd = ElemIndex("sd");
    FitChrom(sf, sd, 0e0, 0e0);
    get_bn_design_elem(sf, 1, Sext, b3[0], a3);
    get_bn_design_elem(sd, 1, Sext, b3[1], a3);

    printf("\nsf = %10.5f, sd = %10.5f", b3[0], b3[1]);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_lat("chromlat.out", globval.bpm, true, 10);

  prtmfile("flat_file.dat");

  if (false) {
    globval.Cavity_on = true;
    get_dynap(delta, 25, 512, true);
  }
}
