#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void get_bn2(const string file_name1, const string file_name2, int n,
	     const bool prt)
{
  char   line[max_str], str[max_str], str1[max_str], *token, *name, *p;
  int    n_prm, Fnum, Knum, order;
  double bnL, bn, C, L;
  FILE   *inf, *fp_lat;

  inf = file_read(file_name1.c_str()); fp_lat = file_write(file_name2.c_str());

  // if n = 0: go to last data set
  if (n == 0) {
    while (fgets(line, max_str, inf) != NULL )
      if (strstr(line, "n = ") != NULL)	sscanf(line, "n = %d", &n);

    fclose(inf); inf = file_read(file_name1.c_str());
  }

  if (prt) {
    printf("\n");
    printf("reading values (n=%d): %s\n", n, file_name1.c_str());
    printf("\n");
  }

  sprintf(str, "n = %d", n);
  do
    fgets(line, max_str, inf);
  while (strstr(line, str) == NULL);

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
	    fprintf(fp_lat, "%s: Quadrupole, L = %8.6f, K = %19.16f,"
		    " N = Nquad, Method = Meth;\n", str1, L, bn);
	  else if (order == Sext)
	    fprintf(fp_lat, "%s: Sextupole, L = %8.6f, K = %19.16f"
		    ", N = Nsext, Method = Meth;\n", str1, L, bn);
	  else {
	    fprintf(fp_lat, "%s: Multipole, L = %8.6f"
		    ", N = 1, Method = Meth,\n", str1, L);
	    fprintf(fp_lat, "     HOM = (%d, %19.16f, %3.1f);\n",
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


int main(int argc, char *argv[])
{
  int           k;
  double        alpha[2], beta[2], eta[2], etap[2];
  ostringstream str;

  // Work around.
  std::vector<std::string> dummy[2];
  dummy[X_].push_back("1"); dummy[Y_].push_back("2");

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
}
