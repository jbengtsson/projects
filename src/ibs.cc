#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void prt_ZAP(const int n)
{
  long int  k;
  FILE      *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n",
	  globval.Cell_nLoc+1, n*Cell[globval.Cell_nLoc].S);
  fprintf(outf, "One super period\n");

  for (k = 0; k <= globval.Cell_nLoc; k++)
    fprintf(outf, "%10.5f %8.5f %9.6f %8.5f %7.3f %8.5f %8.5f %7.5f\n",
	    Cell[k].S,
	    Cell[k].Beta[X_], Cell[k].Alpha[X_],
	    Cell[k].Beta[Y_], Cell[k].Alpha[Y_],
	    Cell[k].Eta[X_], Cell[k].Etap[X_], Cell[k].maxampl[X_][1]);

  fprintf(outf, "0\n");

  fclose(outf);
}


void chk_IBS(void)
{
  int       k;
  double    chi;
  std::ofstream  outf;

  const double  chi_min = 1e-9, chi_max = 1e-5, fact = 2.0;


  file_wr(outf, "chk_IBS.out");

  k = 0; chi_m = chi_min;
  while (chi_m <= chi_max) {
    k++;

    chi = f_IBS(chi_m);

    outf << std::scientific << std::setprecision(5)
	 << std::setw(4) << k << std::setw(12) << chi_m << std::setw(12) << chi;

    chi = get_int_IBS();

    outf << std::scientific << std::setprecision(5)
	 << std::setw(12) << chi << std::endl;

    chi_m *= fact;
  }

  outf.close();
}


void get_Touschek(void)
{
  long int      k;
  int           j;
  double        gamma_z, eps[3];
  double        sum_delta[globval.Cell_nLoc+1][2];
  double        sum2_delta[globval.Cell_nLoc+1][2];
  FILE          *fp;

  const double  Qb = 1.3e-9;

  globval.eps[X_] = 2.040e-9;
  globval.eps[Y_] = 8e-12;
  globval.eps[Z_] = 1.516e-6;

  globval.alpha_z = 2.502e-02; globval.beta_z = 5.733;

  if (true) {
    printf("\n");
    printf("alpha: %11.3e %11.3e %11.3e\n",
	   globval.alpha_rad[X_], globval.alpha_rad[Y_],
	   globval.alpha_rad[Z_]);

    for (j = 0; j < 3; j++)
      eps[j] = globval.eps[j];

    for (j = 1; j <= 0; j++)
      IBS(Qb, globval.eps, eps, true, true);

    for (j = 0; j < 3; j++)
      eps[j] = globval.eps[j];

    for (j = 1; j <= 5; j++)
      IBS_BM(Qb, globval.eps, eps, true, true);
  }

  if (true) {
    globval.delta_RF = 3.0e-2;
//     Touschek(Qb, globval.delta_RF, 0.85e-9, 8e-12, 0.84e-3, 4.4e-3);
    gamma_z = (1.0+sqr(globval.alpha_z))/globval.beta_z;
//     Touschek(Qb, globval.delta_RF, eps[X_], eps[Y_],
// 	     sqrt(gamma_z*eps[Z_]), sqrt(globval.beta_z*eps[Z_]));
    Touschek(Qb, 3.03e-2, 2.446e-9, 9.595e-12, 0.580e-3, 3.33e-3);

    if (false) {
      // initialize momentum aperture arrays
      for(k = 0; k <= globval.Cell_nLoc; k++){
	sum_delta[k][0] = 0.0; sum_delta[k][1] = 0.0;
	sum2_delta[k][0] = 0.0; sum2_delta[k][1] = 0.0;
      }

      Touschek(1.3e-9, globval.delta_RF, false,
	       0.85e-9, 0.008e-9, 0.84e-3, 4.4e-3,
	       512, true, sum_delta, sum2_delta);

      fp = file_write("mom_aper.out");
      for(k = 0; k <= globval.Cell_nLoc; k++)
	fprintf(fp, "%4ld %7.2f %5.3f %6.3f\n",
		k, Cell[k].S, 1e2*sum_delta[k][0], 1e2*sum_delta[k][1]);
    }
  }
}


int main(int argc, char *argv[])
{
  const long  seed = 1121;

  iniranf(seed); setrancut(1.0);

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
//  idprset(-1);

  if (true)
    Read_Lattice(argv[1]);
  else {
    globval.Energy = 3e0;
    rdmfile(argv[1]);
  }

  Ring_GetTwiss(true, 0.0); printglob();

//   prt_ZAP(15);

  Ring_GetTwiss(true, 0.0); printglob();

  prtmfile("flat_file.dat");

  GetEmittance(ElemIndex("cav"), true);

//  chk_IBS();

  get_Touschek();
}
