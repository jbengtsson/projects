#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.EPU = false;

  dnu_dA(4e-3, 4e-3, 0e0, 25); get_ksi2(3e-2);
}
