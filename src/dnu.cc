#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


int main(int argc, char *argv[])
{

  const long  seed = 1121;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
//  idprset(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile("flat_file.dat");

  globval.EPU = false;

  dnu_dA(4e-3, 4e-3, 0.0, 25); get_ksi2(3.0e-2);

}
