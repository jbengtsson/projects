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


int main(int argc, char *argv[])
{
  int           qf, qd, sf, sd;
  double        b2[2], a2, b3[2], b3L[2], a3, a3L;
  ostringstream str;

  const double delta = 3e-2;
  const double nu[] = {102.65/20.0, 68.4/20.0};

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (true) no_sxt();

  Ring_GetTwiss(true, 0e0); printglob();

  // prt_lat_maxlab("m4-20121107-430-bare.out", globval.bpm, true);
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_lat("chromlat.out", globval.bpm, true, 10);

  prtmfile("flat_file.dat");

  if (false) get_alphac2();

  GetEmittance(ElemIndex("cav"), true);

  if (false) {
    qf = ElemIndex("qfe"); qd = ElemIndex("qde");
    FitTune(qf, qd, nu[X_], nu[Y_]);
    get_bn_design_elem(qf, 1, Quad, b2[0], a2);
    get_bn_design_elem(qd, 1, Quad, b2[1], a2);

    printf("\nnu_x = %8.5f nu_y = %8.5f\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("  qfe = %8.5f  qde = %8.5f\n", b2[0], b2[1]);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (false) {
    sf = ElemIndex("sfh"); sd = ElemIndex("sd");
    FitChrom(sf, sd, 0e0, 0e0);
    get_bn_design_elem(sf, 1, Sext, b3[0], a3);
    get_bn_design_elem(sd, 1, Sext, b3[1], a3);
    get_bnL_design_elem(sf, 1, Sext, b3L[0], a3L);
    get_bnL_design_elem(sd, 1, Sext, b3L[1], a3L);

    printf("\nsf = %10.5f (%10.5f), sd = %10.5f (%10.5f)",
	   b3[0], b3L[0], b3[1], b3L[1]);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (false) {
    globval.Cavity_on = true;
    get_dynap(delta, 25, 512, true);
  }
}
