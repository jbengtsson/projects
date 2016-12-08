#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

extern bool freq_map;
extern int  N_Fam, Q_Fam[];

// dynamic aperture
const double delta  = 3e-2; // delta for off-momentum aperture


double get_eps_x1(void)
{
  bool             cav, emit;
  long int         lastpos;
  double           eps_x;
  ss_vect<tps>     A;

  cav = globval.Cavity_on; emit = globval.emittance;

  globval.Cavity_on = false; globval.emittance = false;

  Ring_GetTwiss(false, 0e0);

  putlinmat(6, globval.Ascr, A);

  prt_lin_map(3, A);

  globval.emittance = true;

  Cell_Pass(0, globval.Cell_nLoc, A, lastpos);

  eps_x = 1470e0*pow(globval.Energy, 2)*I5/(I2-I4);

  printf("\neps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	 1e3*eps_x, 1e0-I4/I2, 2e0+I4/I2);

  globval.Cavity_on = cav; globval.emittance = emit;

  return eps_x;
}


void get_b3s(char *file_name)
{
  char   line[max_str], str[max_str], str1[max_str], *token, *name, *p;
  int    n_prm, Fnum, Knum, order;
  double bnL, bn, C, L;
  FILE   *inf, *fp_lat;

  const bool prt = true;

  inf = file_read(file_name); fp_lat = file_write("get_b3s.lat");

  no_sxt();

  if (prt) printf("\nreading values: %s\n\n", file_name);

  fprintf(fp_lat, "\n");
  n_prm = 0;
  while (fgets(line, max_str, inf) != NULL) {
    if (strcmp(line, "\n") == 0) break;
    n_prm++;
    name = strtok_r(line, "(", &p);
    rm_space(name);
    strcpy(str, name); Fnum = ElemIndex(str);
    strcpy(str1, name); upr_case(str1);
    token = strtok_r(NULL, ")", &p); sscanf(token, "%d", &Knum);
    strtok_r(NULL, "=", &p); token = strtok_r(NULL, "\n", &p);
    sscanf(token, "%lf %d", &bnL, &order);
    if (prt) printf("%6s(%2d) = %10.6f %d\n", name, Knum, bnL, order);
   if (Fnum != 0) {
      if (order == 0)
        SetL(Fnum, bnL);
      else
        SetbnL(Fnum, order, bnL);

      L = GetL(Fnum, 1);
      if (Knum == 1) {
	if (order == 0)
	  fprintf(fp_lat, "%s: Drift, L = %8.6f;\n", str1, bnL);
	else {
	  bn = (L != 0.0)? bnL/L : bnL;
	  if (order == Quad)
	    fprintf(fp_lat, "%s: Quadrupole, L = %8.6f, K = %10.6f, N = Nquad"
		    ", Method = Meth;\n", str1, L, bn);
	  else if (order == Sext)
	    fprintf(fp_lat, "%s: Sextupole, L = %8.6f, K = %11.6f"
		    ", N = Nsext, Method = Meth;\n", str1, L, bn);
	  else {
	    fprintf(fp_lat, "%s: Multipole, L = %8.6f"
		    ", N = 1, Method = Meth,\n", str1, L);
	    fprintf(fp_lat, "     HOM = (%d, %13.6f, %3.1f);\n",
		    order, bn, 0.0);
	  }
	}
      }
    } else {
      printf("element %s not found\n", name);
      exit_(1);
    }
  }
  if (prt) printf("\n");

  C = Cell[globval.Cell_nLoc].S; recalc_S();
  if (prt)
    printf("New Cell Length: %5.3f (%5.3f)\n", Cell[globval.Cell_nLoc].S, C);

  fclose(inf); fclose(fp_lat);
}


void prt_dlat(const char *fname, const bool all)
{
  long int      i = 0;
  double        code = 0.0;
  FILE          *outf;

  outf = file_write(fname);

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if (all ||
	((Cell[i].Elem.Pkind == Mpole) &&
	 (Cell[i].Elem.M->n_design == Sext))) {
      switch (Cell[i].Elem.Pkind) {
      case drift:
	code = 0.0;
	break;
      case Mpole:
	if (Cell[i].Elem.M->Pirho != 0.0)
	  code = 0.5;
	else if (Cell[i].Elem.M->PBpar[Quad+HOMmax] != 0)
	  code = sgn(Cell[i].Elem.M->PBpar[Quad+HOMmax]);
	else if (Cell[i].Elem.M->PBpar[Sext+HOMmax] != 0)
	  code = 1.5*sgn(Cell[i].Elem.M->PBpar[Sext+HOMmax]);
	else if (Cell[i].Fnum == globval.bpm)
	  code = 2.0;
	else
	  code = 0.0;
	break;
      default:
	code = 0.0;
	break;
      }

      fprintf(outf, "%4ld %15s %6.2f %4.1f %6.2f %7.4f %6.2f %7.4f\n",
	      i, Cell[i].Elem.PName, Cell[i].S, code,
	      1e2*(Cell[i].Beta[X_]-betas0_[i][X_])/betas0_[i][X_],
	      Cell[i].Nu[X_]-nus0_[i][X_],
	      1e2*(Cell[i].Beta[Y_]-betas0_[i][Y_])/betas0_[i][Y_],
	      Cell[i].Nu[Y_]-nus0_[i][Y_]);
    }
  }

  fclose(outf);
}


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

  if (false) {
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


void get_IBS(const int n, const double ds, const double Qb, const double eps[])
{
  int    i, j;
  double eps0[3], eps1[3], sigma0_s, sigma0_delta, sigma_s, sigma_delta;
  FILE   *outf;

  const char file_name[] = "ibs.out";
  const int n_iter = 20;

  outf = file_write(file_name);

  for (j = 0; j < 3; j++)
    eps0[j] = eps[j];

  // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
  sigma0_s = sqrt(globval.beta_z*eps0[Z_]); sigma0_delta = eps0[Z_]/sigma0_s;

  fprintf(outf, "# sigma0_s eps_x    ratio_x  eps_y    ratio_y"
	  "  sigma_s  sigma_delta\n");
  fprintf(outf, "#   [cm]  [pm.rad]          [pm.rad]   [cm]\n");
  for (i = 1; i <= n; i++) {
    for (j = 0; j < 3; j++)
      eps1[j] = eps0[j];

    printf("\n");
    for (j = 1; j <= n_iter; j++) {
      if (j == 1) {
	printf("\nIBS %d:\n", j);
	IBS_BM(Qb, eps, eps1, true, true);
      } else if ((j == n_iter-1) || (j == n_iter)) {
	printf("\nIBS %d:\n", j);
	IBS_BM(Qb, eps0, eps1, false, true);
      } else
	IBS_BM(Qb, eps0, eps1, false, false);
    }

    sigma_s = sqrt(globval.beta_z*eps1[Z_]);
    sigma_delta = eps1[Z_]/sigma_s;
    fprintf(outf, "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %12.3e\n",
	    1e2*sigma0_s, 1e12*eps1[X_], eps1[X_]/eps0[X_], 1e12*eps1[Y_],
	    eps1[Y_]/eps0[Y_], 1e2*sigma_s, sigma_delta);
    fflush(outf);

    // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
    sigma0_s += ds; eps0[Z_] = sigma0_s*sigma0_delta;
    globval.beta_z = sqr(sigma0_s)/eps0[Z_];
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{
  int           qf, qd, sf, sd, j;
  double        b2[2], a2, b3[2], a3;
  ss_vect<tps>  map;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  const double Qb   = 5e-9, sigma_s = 1e-2, sigma_delta = 1e-3;
  const double nu[] = {101.9/20.0, 27.6/20.0};

  // disable from TPSALib- and LieLib log messages
//  idprset(-1);

  if (true)
    Read_Lattice(argv[1]);
  else {
    globval.Energy = 3e0;
    rdmfile(argv[1]);
  }

  if (false) no_sxt();

  if (false) get_b3s(argv[2]);

  Ring_GetTwiss(true, 0.0); printglob();

  if (false) {
    qf = ElemIndex("qf"); qd = ElemIndex("bh");
    FitTune(qf, qd, nu[X_], nu[Y_]);
    get_bn_design_elem(qf, 1, Quad, b2[0], a2);
    get_bn_design_elem(qd, 1, Quad, b2[1], a2);

    printf("\nnu_x = %7.5f nu_y = %7.5f\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("  qf = %8.5f   qd = %8.5f\n", b2[0], b2[1]);

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
  prt_chrom_lat();

  prtmfile("flat_file.dat");

  if (false) prt_ZAP(20);

  get_eps_x1();
  GetEmittance(ElemIndex("cav"), true);

  if (false) {
    printf("\nalpha_z = %11.3e, beta_z = %10.3e\n",
	   globval.alpha_z,  globval.beta_z);

    globval.eps[Y_] = globval.eps[X_];

    // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
    globval.eps[Z_] = sigma_s*sigma_delta;
    globval.beta_z  = sqr(sigma_s)/globval.eps[Z_];

    printf("\nbeta_z = %5.3f, eps_z = %10.3e\n",
	   globval.beta_z, globval.eps[Z_]);

    get_IBS(21, 1e-2, Qb, globval.eps);
  }

//   get_Touschek();

  if (false) {
    globval.Cavity_on = true;
    n_aper = 25;
    get_dynap(delta, true);
  }
}
