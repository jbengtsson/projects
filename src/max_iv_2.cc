#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


double get_eps_x1(void)
{
  // Evaluate emittance and damping partition from the synchrotron integrals.
  bool         cav, emit;
  long int     lastpos;
  double       eps_x;
  ss_vect<tps> A;

  const bool prt = false;

  cav = globval.Cavity_on; emit = globval.emittance;

  globval.Cavity_on = false; globval.emittance = false;

  Ring_GetTwiss(false, 0e0);

  putlinmat(6, globval.Ascr, A);

  // prt_lin_map(3, A);

  globval.emittance = true;

  Cell_Pass(0, globval.Cell_nLoc, A, lastpos);

  eps_x = 1470e0*pow(globval.Energy, 2)*I5/(I2-I4);

  if (prt)
    printf("eps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	   1e3*eps_x, 1e0-I4/I2, 2e0+I4/I2);

  globval.Cavity_on = cav; globval.emittance = emit;

  return eps_x;
}


void quad_scan(const int n,
	       const char *qf, const double dk_qf,
	       const char *qd, const double dk_qd)
{
  int    i, j, qf_num, qd_num;
  double k_qf, k_qd, k_qf_step, k_qd_step, a2, eps_x;
  FILE   *outf;

  const char file_name[] = "quad_scan.out";

  outf = fopen(file_name, "w");

  qf_num = ElemIndex(qf); qd_num = ElemIndex(qd);
  get_bn_design_elem(qf_num, 1, Quad, k_qf, a2);
  get_bn_design_elem(qd_num, 1, Quad, k_qd, a2);
  k_qf_step = dk_qf/n; k_qf -= dk_qf;
  k_qd_step = dk_qd/n; k_qd += dk_qd;
  for (i = -n; i <= n; i++) {
    set_bn_design_fam(qf_num, Quad, k_qf, a2);
    k_qd -= (2*n+1)*k_qd_step;
    for (j = -n; j <= n; j++) {
      set_bn_design_fam(qd_num, Quad, k_qd, a2);
      Ring_GetTwiss(true, 0e0);
      if (globval.stable)
	eps_x = get_eps_x1();
      else {
	globval.TotalTune[X_] = NAN; globval.TotalTune[Y_] = NAN;
	globval.Chrom[X_] = NAN; globval.Chrom[Y_] = NAN; eps_x = NAN;
      }
      fprintf(outf, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	      k_qf, k_qd,
	      globval.TotalTune[X_], globval.TotalTune[Y_],
	      globval.Chrom[X_], globval.Chrom[Y_], 1e3*eps_x);
      k_qd += k_qd_step;
    }
    fprintf(outf, "\n");
    k_qf += k_qf_step;
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  Read_Lattice(argv[1]);

  no_sxt();

  Ring_GetTwiss(true, 0e0); printglob();

  get_eps_x1();
  GetEmittance(ElemIndex("cav"), true);

  quad_scan(10, "qf", 3e0, "bh", 2e0);
}
