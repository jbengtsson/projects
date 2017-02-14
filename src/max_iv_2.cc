#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


struct param_type {
private:

public:
  int                 n_prm;
  double              bn_tol, svd_cut, step;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn, double *bn_lim);
  void set_prm(double *bn) const;
};


param_type b2_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(ElemIndex(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_prm = Fnum.size();
}


void param_type::ini_prm(double *bn, double *bn_lim)
{
  int    i;
  double an;

  n_prm = Fnum.size();
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      // Multipole.
      get_bn_design_elem(Fnum[i-1], 1, n[i-1], bn[i], an);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Location.
      // bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1]);
      ;
  }
}


void param_type::set_prm(double *bn) const
{
  int i, j;

  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      for (j = 1; j <= GetnKid(Fnum[i-1]); j++)
	set_bn_design_elem(Fnum[i-1], j, n[i-1], bn[i], 0e0);
    else if (n[i-1] == -1)
      set_L(Fnum[i-1], bn[i]);
    else if (n[i-1] == -2)
      // set_bn_s(-Fnum[i-1], n[i-1], bn[i]);
      ;
  }
}


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    S += Cell[j].Elem.PL; Cell[j].S = S;
  }
}


double get_eps_x1(void)
{
  // Evaluate emittance and damping partition from the synchrotron integrals.
  bool         cav, emit;
  long int     lastpos;
  double       eps_x;
  ss_vect<tps> A;

  const bool prt = true;

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
  // Parametric scan of gradients for unit cell.
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


void prt_match(const param_type &b2_prms, const double *b2)
{
  double l5h, l6, l7h, l8;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "bm:  bending, l = 0.166667, t = 0.5, k = %8.5f, t1 = 0.0"
	  ", t2 = 0.0,\n     gap = 0.00, N = Nbend, Method = Meth;\n", b2[1]);
  fprintf(outf, "qfe: quadrupole, l = 0.15, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[2]);
  fprintf(outf, "qde: quadrupole, l = 0.1,  k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[3]);
  fprintf(outf, "qm:  quadrupole, l = 0.15, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[4]);

  l5h = get_L(ElemIndex("l5h"), 1); l6 = get_L(ElemIndex("l6"), 1);
  l7h = get_L(ElemIndex("l7h"), 1); l8 = get_L(ElemIndex("l8"), 1);
  
  fprintf(outf, "\nl5h: drift, l = %7.5f;\n", l5h);
  fprintf(outf, "l6:  drift, l = %7.5f;\n", l6);
  fprintf(outf, "l7h: drift, l = %7.5f;\n", l7h);
  fprintf(outf, "l8:  drift, l = %7.5f;\n", l8);

  fclose(outf);
}


double f_match(double *b2)
{
  // Minimize linear linear chromaticity for super period.
  // Lattice: super period.

  static double chi2_ref = 1e30;

  int          i, loc1, loc2, loc3;
  double       tr[2], chi2;

  const double scl_eta = 1e5,
	       beta0[] = {3.0,     3.0},
               beta1[] = {1.35633, 1.91479};

  b2_prms.set_prm(b2);

  // Zero dispersion after 2nd BM.
  loc1 = Elem_GetPos(ElemIndex("bm"), 2);
  // Optics at center of Unit Cell.
  loc2 = Elem_GetPos(ElemIndex("sfh"), 1);
  // Optics at end of super period.
  loc3 = globval.Cell_nLoc;

  // Turn on cavity; or else dispersion will be zeroed in Ring_Twiss.
  globval.Cavity_on = true;
  Ring_GetTwiss(false, 0e0);

  tr[X_] = globval.OneTurnMat[x_][x_] + globval.OneTurnMat[px_][px_];
  tr[Y_] = globval.OneTurnMat[y_][y_] + globval.OneTurnMat[py_][py_];
  // printf("trace: %6.3f %6.3f\n", tr[X_], tr[Y_]);

  chi2 = 0e0;
  chi2 += sqr(scl_eta*Cell[loc1].Eta[X_]);
  chi2 += sqr(scl_eta*Cell[loc1].Etap[X_]);
  chi2 += sqr(Cell[loc2].Beta[X_]-beta1[X_]);
  chi2 += sqr(Cell[loc2].Beta[Y_]-beta1[Y_]);
  chi2 += sqr(Cell[loc3].Beta[X_]-beta0[X_]);
  chi2 += sqr(Cell[loc3].Beta[Y_]-beta0[Y_]);
  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0)) chi2 += 1e10;
  for (i = 1; i <= b2_prms.n_prm; i++) {
    if (b2[i] < b2_prms.bn_min[i-1]) chi2 += 1e10;
    if (fabs(b2[i]) > b2_prms.bn_max[i-1]) chi2 += 1e10;
  }

  if (chi2 < chi2_ref) {
    printf("\nchi2: %12.5e, %12.5e\n", chi2, chi2_ref);
    printf("b:    %10.3e %10.3e %8.3f %8.3f %8.3f %8.3f\n",
	   Cell[loc1].Eta[X_], Cell[loc1].Etap[X_],
	   Cell[loc2].Beta[X_], Cell[loc2].Beta[Y_],
	   Cell[loc3].Beta[X_], Cell[loc3].Beta[Y_]);
    printf("b2s: ");
    for (i = 1; i <= b2_prms.n_prm; i++)
      printf("%9.5f", b2[i]);
    printf("\n");

    prtmfile("flat_file.fit");
    prt_match(b2_prms, b2);

    get_S();
    prt_lat("linlat.out", globval.bpm, true);
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void opt_match(param_type &b2_prms)
{
  // Minimize linear linear chromaticity for super period.
  // Lattice: super period.

  int    n_b2, i, j, iter;
  double *b2, *b2_lim, **xi, fret;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_match);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


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

  no_sxt();

  Ring_GetTwiss(true, 0e0); printglob();

  if (false) {
    get_eps_x1();
    GetEmittance(ElemIndex("cav"), true);
  }

  if (false) quad_scan(10, "qf", 3e0, "bh", 2e0);

  if (true) {
    b2_prms.add_prm("bm",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qfe",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qde",  2, 0.0, 30.0, 1.0);
    b2_prms.add_prm("qm",   2, 0.0, 25.0, 1.0);

    b2_prms.add_prm("l5h", -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l6",  -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l7h", -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l8",  -1, 0.05, 1.0,  0.01);

    // b2_prms.add_prm("bm",  -2,  0.05, 0.5, 1.0);
    // b2_prms.add_prm("qde", -2,  0.05, 0.5, 1.0);
    // b2_prms.add_prm("qm",  -2,  0.05, 0.5, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.svd_cut = 1e-8; b2_prms.step = 0.001;

    no_sxt();
    opt_match(b2_prms);
  }
}
