#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void prt_name(FILE * outf, const char *name)
{
    int j, k, len;

    len = strlen(name);

    j = 0;
    do {
	fprintf(outf, "%c", name[j]);
	j++;
    } while ((j < len) && (name[j] != ' '));
    fprintf(outf, ",");
    for (k = j; k < len; k++)
	fprintf(outf, "%c", name[k]);
}


void prt_lat_maxlab(const char *fname, const int Fnum, const bool all)
{
  // Generate CSV file for the linear optics.
  long int i = 0;
  FILE *outf;

  outf = fopen(fname, "w");
  fprintf(outf, "#        name               s      type"
	  "    alphax    betax      nux       etax     etapx");
  fprintf(outf, "      alphay    betay      nuy      etay      etapy\n");
  fprintf(outf, "#                          [m]"
	  "                        [m]                 [m]");
  fprintf(outf, "                            [m]                [m]\n");
  fprintf(outf, "#\n");

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if (all || (Cell[i].Fnum == Fnum)) {
      fprintf(outf, "%4ld, ", i);
      prt_name(outf, Cell[i].Elem.PName);
      fprintf(outf, " %9.5f, %4.1f,"
           " %9.5f, %8.5f, %8.5f, %8.5f, %8.5f,"
           " %9.5f, %8.5f, %8.5f, %8.5f, %8.5f\n",
           Cell[i].S, get_code(Cell[i]),
           Cell[i].Alpha[X_], Cell[i].Beta[X_], Cell[i].Nu[X_],
           Cell[i].Eta[X_], Cell[i].Etap[X_],
           Cell[i].Alpha[Y_], Cell[i].Beta[Y_], Cell[i].Nu[Y_],
           Cell[i].Eta[Y_], Cell[i].Etap[Y_]);
    }
  }

  fclose(outf);
}


void get_cod_rms(const double dx, const double dy,
		 const int n_seed, const bool all)
{
  bool                cod;
  int                 i, j, k, n, n_cod;
  std::vector<double> x1[6], x2[6], x_mean[6], x_sigma[6];
  FILE                *fp;

  const int n_cod_corr = 5;

  globval.Cavity_on = false;

  for (j = 0; j <= globval.Cell_nLoc; j++)
    for (k = 0; k < 6; k++) {
      x1[k].push_back(0e0); x2[k].push_back(0e0);
    }
  
  fp = file_write("cod_rms.out");
  
  n_cod = 0;
  for (i = 0; i < n_seed; i++) {
    printf("\norb_corr: seed no %d\n", i+1);

    misalign_rms_type(Dip,  dx, dy, 0e0, true);
    misalign_rms_type(Quad, dx, dy, 0e0, true);
    
    cod = orb_corr(n_cod_corr);

    if (cod) {
      n_cod++;

      n = 0;
      for (j = 0; j <= globval.Cell_nLoc; j++)
	if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		    (Cell[j].Elem.M->n_design == Sext))) {
	  n++;
	  for (k = 0; k < 6; k++) {
	    x1[k][n-1] += Cell[j].BeamPos[k];
	    x2[k][n-1] += sqr(Cell[j].BeamPos[k]);
	  }
	}
    } else
      printf("orb_corr: failed\n");

    // Reset orbit trims.
    set_bn_design_fam(globval.hcorr, Dip, 0e0, 0e0);
    set_bn_design_fam(globval.vcorr, Dip, 0e0, 0e0);
  }

  printf("\nget_cod_rms: no of seeds %d, no of cods %d\n", n_seed, n_cod);

  n = 0;
  for (j = 0; j <= globval.Cell_nLoc; j++)
    if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		(Cell[j].Elem.M->n_design == Sext))) {
      n++;
      for (k = 0; k < 6; k++) {
	x_mean[k].push_back(x1[k][n-1]/n_cod);
	x_sigma[k].push_back(sqrt((n_cod*x2[k][n-1]-sqr(x1[k][n-1]))
				  /(n_cod*(n_cod-1.0))));
      }
      fprintf(fp, "%8.3f %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      Cell[j].S, get_code(Cell[j]),
	      1e3*x_mean[x_][n-1], 1e3*x_sigma[x_][n-1],
	      1e3*x_mean[y_][n-1], 1e3*x_sigma[y_][n-1]);
    } else
      fprintf(fp, "%8.3f %6.2f\n", Cell[j].S, get_code(Cell[j]));
  
  fclose(fp);
}


int main(int argc, char *argv[])
{
  int           qf, qd, sf, sd;
  double        b2[2], a2, b3[2], b3L[2], a3, a3L;
  ostringstream str;

  const long   seed = 1121;
  const double delta = 5e-2;
  const double nu[] = {102.22/20.0, 68.18/20.0};

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (false) no_sxt();

  Ring_GetTwiss(true, 0e0); printglob();

  if (false) get_alphac2();

  prtmfile("flat_file.dat");

  // prt_lat_maxlab("m4-20121107-430-bare.out", globval.bpm, true);
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_lat("chromlat.out", globval.bpm, true, 10);

  if (false) {
    iniranf(seed); setrancut(1e0);
    globval.bpm = ElemIndex("mon");
    globval.hcorr = ElemIndex("ch"); globval.vcorr = ElemIndex("cv");

    gcmat(globval.bpm, globval.hcorr, 1);
    gcmat(globval.bpm, globval.vcorr, 2);

    get_cod_rms(50e-6, 50e-6, 100, true);

    exit(0);
  }

  GetEmittance(ElemIndex("cav"), true);

  if (false) {
    qf = ElemIndex("qfe"); qd = ElemIndex("qde");
    FitTune(qf, qd, nu[X_], nu[Y_]);
    get_bn_design_elem(qf, 1, Quad, b2[0], a2);
    get_bn_design_elem(qd, 1, Quad, b2[1], a2);

    printf("\nnu_x = %8.5f nu_y = %8.5f\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("  qfe = %8.5f  qde = %8.5f\n", b2[0], b2[1]);

    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    sf = ElemIndex("sfh"); sd = ElemIndex("sd");
    FitChrom(sf, sd, 0e0, 0e0);
    get_bn_design_elem(sf, 1, Sext, b3[0], a3);
    get_bn_design_elem(sd, 1, Sext, b3[1], a3);
    get_bnL_design_elem(sf, 1, Sext, b3L[0], a3L);
    get_bnL_design_elem(sd, 1, Sext, b3L[1], a3L);

    printf("\nsfh = %10.5f (%10.5f), sd = %10.5f (%10.5f)\n",
	   b3[0], b3L[0], b3[1], b3L[1]);

    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    globval.Cavity_on = false;
    get_dynap(delta, 25, 512, true);
  }
}
