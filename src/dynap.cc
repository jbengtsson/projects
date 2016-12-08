#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


class orb_corr_type {
private:
  double                **A, *w, **U, **V, *b, *x, eps;

public:
  bool                  hor, periodic;
  int                   m, n;
  std::vector<long int> bpms, corrs;

  void alloc(const long int i0, const long int i1, const long int i2,
	     const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> &corr_Fam_names,
	     const bool hor, const bool periodic, const double eps);
  void alloc(const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> &corr_Fam_names,
	     const bool hor, const bool periodic, const double eps);
  void dealloc(void);

  void get_trm_mat(void);
  void get_orm_mat(void);
  void svd_decomp(void);
  void solve(const double scl) const;
  void clr_trims(void);
};


std::vector<long int> get_elem(const long int i0, const long int i1,
			       const std::vector<std::string> &names)
{
  long int              loc;
  int                   j, k, Fnum;
  std::vector<long int> elems;

  for (j = 0; j < (int)names.size(); j++) {
    Fnum = ElemIndex(names[j]);
    for (k = 1; k <= GetnKid(Fnum); k++) {
      loc = Elem_GetPos(Fnum, k);
      if ((i0 <= loc) && (loc <= i1)) elems.push_back(loc);
    }
  }

  return elems;
}


void orb_corr_type::alloc(const long int i0, const long int i1,
			  const long int i2,
			  const std::vector<string> &bpm_Fam_names,
			  const std::vector<string> &corr_Fam_names,
			  const bool hor, const bool periodic,
			  const double eps)
{
  this->hor = hor; this->periodic = periodic; this->eps = eps;

  bpms  = get_elem(i0, i2, bpm_Fam_names);
  corrs = get_elem(i0, i1, corr_Fam_names);

  m = bpms.size(); n = corrs.size();

  A = dmatrix(1, m, 1, n); U = dmatrix(1, m, 1, n);
  w = dvector(1, n);       V = dmatrix(1, n, 1, n);
  b = dvector(1, m);       x = dvector(1, n);

  if (!periodic)
    get_trm_mat();
  else
    get_orm_mat();

  svd_decomp();

  printf("\nalloc: n_bpm = %d, n_corr = %d\n", m, n);
}


void orb_corr_type::alloc(const std::vector<string> &bpm_Fam_names,
			  const std::vector<string> &corr_Fam_names,
			  const bool hor, const bool periodic,
			  const double eps)
{
  alloc(0, globval.Cell_nLoc, globval.Cell_nLoc, bpm_Fam_names, corr_Fam_names,
	hor, periodic, eps);
}


void orb_corr_type::dealloc(void)
{
  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dvector(w, 1, n);       free_dmatrix(V, 1, n, 1, n);
  free_dvector(b, 1, m);       free_dvector(x, 1, n);
}


void orb_corr_type::get_trm_mat(void)
{
  int      plane, i, j;
  long int loc_bpm, loc_corr;
  double   betai, betaj, nui, nuj;

  plane = (hor)? 0 : 1;

  for (i = 0; i < m; i++) {
    loc_bpm = bpms[i];
    betai = Cell[loc_bpm].Beta[plane]; nui = Cell[loc_bpm].Nu[plane];
    for (j = 0; j < n; j++) {
      loc_corr = corrs[j];
      betaj = Cell[loc_corr].Beta[plane]; nuj = Cell[loc_corr].Nu[plane];
      A[i+1][j+1] = (loc_bpm > loc_corr)?
	sqrt(betai*betaj)*sin(2.0*M_PI*(nui-nuj)) : 0e0;
    }
  }
}


void orb_corr_type::get_orm_mat(void)
{
  int      plane, i, j;
  long int loc_bpm, loc_corr;
  double   nu, betai, betaj, nui, nuj, spiq;

  plane = (hor)? 0 : 1;

  nu = globval.TotalTune[plane]; spiq = sin(M_PI*nu);

  for (i = 0; i < m; i++) {
    loc_bpm = bpms[i];
    betai = Cell[loc_bpm].Beta[plane]; nui = Cell[loc_bpm].Nu[plane];
    for (j = 0; j < n; j++) {
      loc_corr = corrs[j];
      betaj = Cell[loc_corr].Beta[plane]; nuj = Cell[loc_corr].Nu[plane];
      A[i+1][j+1] = 
	sqrt(betai*betaj)/(2.0*spiq)*cos(nu*M_PI-fabs(2.0*M_PI*(nui-nuj)));
    }
  }
}


void orb_corr_type::svd_decomp(void)
{
  int i, j;

  const int n_prt = 5;

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
      U[i][j] = A[i][j];

  dsvdcmp(U, m, n, w, V);

  printf("\ngtcmat singular values:\n");
  for (j = 1; j <= n; j++) {
    printf("%11.3e", w[j]);
    if (w[j] < eps) {
      w[j] = 0e0;
      printf(" (zeroed)");
    }
    if (j % n_prt == 0) printf("\n");
  }
  if (n % n_prt != 0) printf("\n");
}


void orb_corr_type::solve(const double scl) const
{
  int      plane, j;
  long int loc;

  plane = (hor)? 0 : 1;

  for (j = 0; j < m; j++) {
    loc = bpms[j];
    b[j+1] = -Cell[loc].BeamPos[2*plane] + Cell[loc].dS[plane];
  }
      
  dsvbksb(U, w, V, m, n, b, x);

  for (j = 0; j < n; j++) {
    loc = corrs[j];
    if (plane == 0)
      set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip,
			   -scl*x[j+1], 0e0);
    else
      set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip,
			   0e0, scl*x[j+1]);
  }
}


void orb_corr_type::clr_trims(void)
{
  long int loc;
  int      j;

  for (j = 0; j < n; j++) {
    loc = corrs[j];
    set_bnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, 0e0, 0e0);
  }
}


void codstat(double mean[], double sigma[], double xmax[], const long lastpos,
	     const bool all, const std::vector<long int> &bpms)
{
  long    i, n, loc;
  int     j;
  Vector2 sum, sum2;

  for (j = 0; j < 2; j++) {
    sum[j] = 0e0; sum2[j] = 0e0; xmax[j] = 0e0;
  }

  n = 0;
  if (all) {
    for (i = 0; i < lastpos; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	sum[j]  += Cell[i].BeamPos[j*2];
	sum2[j] += sqr(Cell[i].BeamPos[j*2]);
	xmax[j] =  max(xmax[j], fabs(Cell[i].BeamPos[j*2]));
      }
    }
  } else {
    for (i = 0; i < (int)bpms.size(); i++) {
      n++;
      for (j = 0; j < 2; j++) {
	loc = bpms[i];
	sum[j]  += Cell[loc].BeamPos[j*2];
	sum2[j] += sqr(Cell[loc].BeamPos[j*2]);
	xmax[j] =  max(xmax[j], fabs(Cell[loc].BeamPos[j*2]));
      }
    }
  }

  for (j = 0; j < 2; j++) {
    if (n != 0)
      mean[j] = sum[j] / n;
    else
      mean[j] = -1e0;
    if (n != 0 && n != 1) {
      sigma[j] = (n*sum2[j]-sqr(sum[j]))/(n*(n-1e0));
    } else
      sigma[j] = 0e0;
    if (sigma[j] >= 0e0)
      sigma[j] = sqrt(sigma[j]);
    else
      sigma[j] = -1e0;
  }
}


void thread_beam(const int n_cell, const string &Fam_name,
		 const std::vector<string> &bpm_Fam_names,
		 const std::vector<string> corr_Fam_names[],
		 const int n_orbit, const double scl)
{
  // Thread beam one super period at the time.

  long int              lastpos, i0, i1, i2, j1, j2 = 0;
  int                   i, j, Fnum;
  Vector2               mean, sigma, max;
  ss_vect<double>       ps;
  orb_corr_type         orb_corr[2];
  std::vector<long int> bpms;

  const double eps = 1e-4;

  Fnum = ElemIndex(Fam_name);
  i0 = Elem_GetPos(Fnum, 1); i1 = Elem_GetPos(Fnum, 2);
  i2 = Elem_GetPos(Fnum, 4);
  for (j = 0; j < 2; j++) {
    orb_corr[j].alloc(i0, i1, i2, bpm_Fam_names, corr_Fam_names[j],
		      j == 0, false, eps);
  }

  ps.zero(); Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  codstat(mean, sigma, max, lastpos, true, orb_corr[X_].bpms);
  printf("\nInitial rms trajectory (all):   x = %7.1e mm, y = %7.1e mm\n",
	 1e3*sigma[X_], 1e3*sigma[Y_]);

   for (i = 0; i < n_cell; i++) {
    i0 = Elem_GetPos(Fnum, 2*i+1); i1 = Elem_GetPos(Fnum, 2*i+2);
    i2 = Elem_GetPos(Fnum, 2*i+4);
    for (j = 0; j < 2; j++) {
      orb_corr[j].corrs = get_elem(i0, i1, corr_Fam_names[j]);
      if (i != n_cell-1)
	orb_corr[j].bpms = get_elem(i0, i2, bpm_Fam_names);
      else {
	orb_corr[j].bpms = get_elem(i0, i1, bpm_Fam_names);
	j1 = Elem_GetPos(Fnum, 1); j2 = Elem_GetPos(Fnum, 2);
	bpms = get_elem(j1, j2,	bpm_Fam_names);
	orb_corr[j].bpms.insert(orb_corr[j].bpms.end(),
				bpms.begin(), bpms.end());
      }
    }

    for (j = 1; j <= n_orbit; j++) {
      ps.zero();
      Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
      if (i != n_cell-1) Cell_Pass(0, j2, ps, lastpos);
      orb_corr[0].solve(scl); orb_corr[1].solve(scl);
    }
  }

   ps.zero(); Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
   codstat(mean, sigma, max, lastpos, true, orb_corr[X_].bpms);
   printf("Corrected rms trajectory (all): x = %7.1e mm, y = %7.1e mm\n",
	  1e3*sigma[X_], 1e3*sigma[Y_]);

  for (j = 0; j < 2; j++)
    orb_corr[j].dealloc();
}


void cod_ini(const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> corr_Fam_names[],
	     orb_corr_type orb_corr[])
{
  int j;

  const double eps = 1e-4;

  for (j = 0; j < 2; j++)
    orb_corr[j].alloc(bpm_Fam_names, corr_Fam_names[j], j == 0, true, eps);
}


bool cod_correct(const int n_orbit, const double scl, orb_corr_type orb_corr[])
{
  bool     cod = false;
  long int lastpos;
  int      j;
  Vector2  mean, sigma, max;

  for (j = 1; j <= n_orbit; j++) {
    cod = getcod(0e0, lastpos);
    if (cod) {
      if (j == 1) {
	codstat(mean, sigma, max, globval.Cell_nLoc, true, orb_corr[X_].bpms);
	printf("\nInitial rms cod (all):          x = %7.1e mm, y = %7.1e mm\n",
	       1e3*sigma[X_], 1e3*sigma[Y_]);
	codstat(mean, sigma, max, globval.Cell_nLoc, false, orb_corr[X_].bpms);
	printf("Initial rms cod (bpms):         x = %7.1e mm, y = %7.1e mm\n",
	       1e3*sigma[X_], 1e3*sigma[Y_]);
      }

      orb_corr[0].solve(scl); orb_corr[1].solve(scl);

      codstat(mean, sigma, max, globval.Cell_nLoc, false, orb_corr[X_].bpms);
      printf("Corrected rms orbit (bpms):     x = %7.1e mm, y = %7.1e mm\n",
	     1e3*sigma[X_], 1e3*sigma[Y_]);
    } else
      printf("\ncod_correct failed");
  }

  codstat(mean, sigma, max, globval.Cell_nLoc, true, orb_corr[X_].bpms);
  printf("Corrected rms orbit (all):      x = %7.1e mm, y = %7.1e mm\n",
	 1e3*sigma[X_], 1e3*sigma[Y_]);

  return cod;
}


bool cod_correct(const int n_cell, const int n_orbit, const double scl,
		 const int k, orb_corr_type orb_corr[],
		 const param_data_type &params)
{
  bool            cod = false;
  long int        lastpos;
  ss_vect<double> ps;

  // Load misalignments; set seed, no scaling of rms errors.
  params.LoadAlignTol(false, 1e0, true, k);

  if (params.bba) {
    // Beam based alignment
    params.Align_BPMs(Quad);
  }

  orb_corr[X_].clr_trims(); orb_corr[Y_].clr_trims();
  
  cod = getcod(0e0, lastpos);

  if (!cod) {
    orb_corr[X_].clr_trims(); orb_corr[Y_].clr_trims();
    thread_beam(n_cell, "ls", params.bpm_Fam_names, params.corr_Fam_names,
		n_orbit, scl);

    // prt_cod("cod.out", globval.bpm, true);    
  }

  cod = cod_correct(n_orbit, scl, orb_corr);

  prt_cod("cod.out", globval.bpm, true);    

  return cod;
}


bool track(const double x, const double px, const double y, const double py,
	   const double delta, const long int n, const double f_rf,
	   const bool prt)
{
  long int         i, lastpos;
  ss_vect<double>  ps;
  std::ofstream         os;

  ps[x_] = x; ps[px_] = px; ps[y_] = y; ps[py_] = py;
  ps[delta_] = delta; ps[ct_] = 0.0;

  if (prt) {
    os.open("track.out", std::ios::out);
    os << "# Tracking with Thor" << std::endl;
    os << "#" << std::endl;
    os << "#  n       x           p_x          y            p_y  "
       << "       delta         cdt" << std::endl;
    if (f_rf == 0.0) {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [mm]" << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	 << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	 << std::setw(24) << 1e2*ps[delta_] 
	 << std::setw(24) << 1e3*ps[ct_] << std::endl;
    } else {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [deg]" << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	 << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	 << std::setw(24) << 1e2*ps[delta_] 
	 << std::setw(24) << 2.0*f_rf*180.0*ps[ct_]/c0 << std::endl;
    }
    os << "#" << std::endl;
  }

  for (i = 1; i <= n; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    if (lastpos == globval.Cell_nLoc) {
      if (prt) {
	if (f_rf == 0.0)
	  os << std::scientific << std::setprecision(16)
	     << std::setw(4) << i
	     << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	     << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	     << std::setw(24) << 1e2*ps[delta_] 
	     << std::setw(24) << 1e3*ps[ct_] << std::endl;
	else
	  os <<std:: scientific << std::setprecision(16)
	     << std::setw(4) << i
	     << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	     << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	     << std::setw(24) << 1e2*ps[delta_] 
	     << std::setw(24) << 2.0*f_rf*180.0*ps[ct_]/c0 << std::endl;
      }
    } else
      return false;
  }
  if (prt) os.close();

  return true;
}


void get_r_stable(double &r, const double phi, const double delta,
		  const long int n, const double eps)
{
  /* Binary search for dynamic aperture. */
  bool    lost = false;
  double  r_min = 0.0, r_max = r;

  while (!lost ) {
    lost = ! track(r_max*cos(phi), 0.0, r_max*sin(phi), 0.0, delta,
		   n, 0, false);
    r_max *= 2.0;
  }
  while (r_max-r_min >= eps) {
    r = r_min + (r_max-r_min)/2.0;
    lost = !track(r*cos(phi), 0.0, r*sin(phi), 0.0, delta, n, 0, false);
    if (!lost)
      r_min = r;
    else
      r_max = r;
  }
  r = r_min + (r_max-r_min)/2.0;
}


double get_dynap(FILE *fp, const double r, const double delta, const int n,
		 const double eps, const int n_pts,
		 double x_min[], double x_max[])
{
  /* Determine the dynamic aperture by tracking.
     Assumes mid-plane symmetry.                                              */

  int           i, j;
  double        r1, phi, x0[2] = {0e0, 0e0}, x1[2], x2[2], DA;

  fprintf(fp, "\n");
  fprintf(fp, "# Dynamic Aperture:\n");
  fprintf(fp, "#    x      y\n");
  fprintf(fp, "#   [mm]   [mm]\n");
  fprintf(fp, "#\n");

  for (i = 0; i < 2; i++) {
    x_min[i] = 0.0; x_max[i] = 0.0;
  }

  DA = 0.0; r1 = r;
  for (i = 0; i < n_pts; i++) {
    phi = i*pi/(n_pts-1);
    if (i == 0)
      phi = 1e-3;
    else if (i == n_pts-1)
      phi -= 1e-3;
    get_r_stable(r1, phi, delta, n, eps);
    x2[X_] = r1*cos(phi); x2[Y_] = r1*sin(phi);
    for (j = 0; j <= 1; j++) {
      x_min[j] = min(x2[j], x_min[j]); x_max[j] = max(x2[j], x_max[j]);
    }
    if (i == 0) {
      x0[X_] = x2[X_]; x0[Y_] = x2[Y_];
    } else
      DA += x1[X_]*x2[Y_] - x2[X_]*x1[Y_];

    fprintf(fp, "  %6.2f %6.2f\n", 1e3*x2[X_], 1e3*x2[Y_]);

    x1[X_] = x2[X_]; x1[Y_] = x2[Y_];
  }
  DA += x2[X_]*x0[Y_] - x0[X_]*x2[Y_];
  // x2 from mid-plane symmetry
  DA = fabs(DA)/sqrt(Cell[globval.Cell_nLoc].Beta[X_]
       *Cell[globval.Cell_nLoc].Beta[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "# DA^ = %6.1f mm^2"
	  ", x^ = %6.2f - %5.2f mm, y^ = %6.2f - %5.2f mm\n",
	  1e6*DA, 1e3*x_min[X_], 1e3*x_max[X_], 1e3*x_min[Y_], 1e3*x_max[Y_]);

  fflush(fp);

  return DA;
} 

void get_DA_bare(const double delta, const int n_delta)
{
  char    str[max_str];
  int     j;
  double  DA, x_min[2], x_max[2], x_hat[2], d;
  FILE    *DA_bare, *fp;

  globval.Cavity_on = true;

  DA_bare = file_write("DA_bare.out");

  fprintf(DA_bare, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  Cell[globval.Cell_nLoc].Beta[X_], Cell[globval.Cell_nLoc].Beta[Y_]);
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# Ideal lattice\n");
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# delta   DA      Ax        Ay      x^    y^\n");
  fprintf(DA_bare, "#  [%%]  [mm^2] [mm.mrad] [mm.mrad] [mm]  [mm]\n");

  for (j = 0; j <= n_delta; j++) {
    d = (n_delta > 0)? (double)j/(double)n_delta*delta : 0.0;

    sprintf(str, "DA_bare_%4.2f.out", 1e2*d); fp = file_write(str);

    DA = get_dynap(fp, 10e-3, d, n_track, 0.1e-3, n_aper, x_min, x_max); 

    fclose(fp);

    x_hat[X_] = (x_max[X_]-x_min[X_])/2.0; x_hat[Y_] = x_max[Y_];

    fprintf(DA_bare, "  %5.2f %6.1f   %4.1f      %4.1f   %4.1f  %4.1f\n", 
	    1e2*d, 1e6*DA,
	    1e6*sqr(x_hat[X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat[Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e3*x_hat[X_], 1e3*x_hat[Y_]);
  
    fflush(DA_bare);
  }

  fclose(DA_bare);
}


void get_mean_sigma(const int n, double &m, double &s)
{

  m = (n > 0)? m/n : 0e0; 
  s = (n > 1)? sqrt((s-n*sqr(m))/(n-1)) : 0e0;
}


void err_and_corr(const char *param_file,
		  const double delta, const int n_delta)
{
  bool            cod = false, cav, rad, aper;
  char            str[max_str];
  int             j, k;
  double          DA, x_min[2], x_max[2], d[n_delta+1], x_hat;
  double          DA_m[n_delta+1], DA_s[n_delta+1];
  double          x_hat_m[n_delta+1][2], x_hat_s[n_delta+1][2];
  param_data_type params;
  orb_corr_type   orb_corr[2];
  FILE            *DA_real = NULL, *fp[n_delta+1];

  const int    n_cell = 20;
  const double scl    = 0.5;
  
  // Save state.
  cav = globval.Cavity_on; rad = globval.radiation; aper = globval.Aperture_on;

  // Set state.
  globval.Cavity_on   = false; globval.radiation = false;
  globval.Aperture_on = false;

  params.get_param(param_file);

  n_track = 512;

  Ring_GetTwiss(true, 0.0); printglob();

  // Store optics function values at sextupoles.
  params.get_bare();

  if (params.DA_bare) get_DA_bare(delta, n_delta);

  globval.Cavity_on = false;

  if (params.ae_file != "") {
    if (params.bba) {
      // Beam based alignment
      params.Align_BPMs(Quad);
    }

    cod_ini(params.bpm_Fam_names, params.corr_Fam_names, orb_corr);
  }
 
  // Linear coupling correction and dispersion wave initialization.
  if (params.n_lin > 0) params.ini_skew_cor(params.disp_wave_y);

  for (j = 0; j <= n_delta; j++) {
    d[j] = (n_delta > 0)? (double)j/(double)n_delta*delta : 0.0;
    sprintf(str, "DA_real_%4.2f.out", 1e2*d[j]); fp[j] = file_write(str);

    DA_m[j] = 0.0; DA_s[j] = 0.0;

    for (k = 0; k <= 1; k++) {
      x_hat_m[j][k] = 0.0; x_hat_s[j][k] = 0.0;
    }
  }

  DA_real = file_write("DA_real.out");
  fprintf(DA_real, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  Cell[globval.Cell_nLoc].Beta[X_],
	  Cell[globval.Cell_nLoc].Beta[Y_]);
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# Real lattice\n");
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# deltaP        DA               Ax              Ay"
	  "            x^              y^\n");
  fprintf(DA_real, "#   %%          mm^2              mm              mm"
	    "            mm              mm\n");

  for (j = 1; j <= params.n_stat; j++) {
    globval.Cavity_on = false;

    if (params.fe_file != "")
      params.LoadFieldErr(false, 1e0, true);

    if (params.ae_file != "") {
      cod = cod_correct(n_cell, params.n_orbit, scl, j, orb_corr, params);
    }
      
    // cod = CorrectCOD(n_orbit, 0.3);

    printf("\n");
    if (cod) {
      printf("err_and_corr: orbit correction completed\n");
      if (j == 1) prt_cod("cod.out", globval.bpm, true);
 
      Ring_GetTwiss(true, 0.0); printglob();

      GetEmittance(ElemIndex("cav"), true);

      if (params.n_lin > 0) params.corr_eps_y();

      if (params.ap_file != "") params.LoadApers(1.0, 1.0);

      Ring_GetTwiss(true, 0.0); printglob();

      GetEmittance(ElemIndex("cav"), true);

      prt_beamsizes();

      globval.Cavity_on = true;
      for (k = 0; k <= n_delta; k++) {
	DA = get_dynap(fp[k], 10e-3, d[k], n_track, 0.1e-3, n_aper,
		       x_min, x_max);
	DA_m[k] += DA; DA_s[k] += sqr(DA);
	x_hat = (x_max[X_]-x_min[X_])/2.0;
	x_hat_m[k][X_] += x_hat; x_hat_s[k][X_] += sqr(x_hat);
	x_hat = x_max[Y_];
	x_hat_m[k][Y_] += x_hat; x_hat_s[k][Y_] += sqr(x_hat);
      }

      // reset orbit trims
      zero_trims();

      if (params.n_lin > 0)
	// reset skew quads
	set_bn_design_fam(globval.qt, Quad, 0.0, 0.0);

      if (params.N_calls > 0) params.reset_quads();  
    } else
      chk_cod(cod, "error_and_correction");
  }

  for (j = 0; j <= n_delta; j++) {
    fclose(fp[j]);

    get_mean_sigma(params.n_stat, DA_m[j], DA_s[j]);
    for (k = 0; k <= 1; k++)
      get_mean_sigma(params.n_stat, x_hat_m[j][k], x_hat_s[j][k]);

    fprintf(DA_real, "  %5.2f  %6.1f \xB1 %6.1f"
	    "   %5.1f \xB1 %5.2f    %5.1f \xB1 %5.2f"
	    "   %4.1f \xB1 %5.2f     %4.1f \xB1 %5.2f\n", 
	    d[j]*1e2,
	    1e6*DA_m[j], 1e6*DA_s[j],
	    1e6*sqr(x_hat_m[j][X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat_s[j][X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat_m[j][Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e6*sqr(x_hat_s[j][Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e3*x_hat_m[j][X_], 1e3*x_hat_s[j][X_],
	    1e3*x_hat_m[j][Y_], 1e3*x_hat_s[j][Y_]);
  }

  fclose(DA_real);


  if (params.ae_file != "") {
    for (j = 0; j < 2; j++)
      orb_corr[j].dealloc();
  }

  // Reset state.
  globval.Cavity_on = cav; globval.radiation = rad; globval.Aperture_on = aper;
}


int main(int argc, char *argv[])
{
  int    n_delta;
  double delta;

  // Problem with STL vector class.
  std::vector<std::string> dummy[2];
  dummy[X_].push_back("1"); dummy[Y_].push_back("2");

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.Aperture_on = false;

  if (argc < 3) {
    printf("*** bad command line\n");
    exit(1);
  }

  sscanf(argv[2], "%d", &n_delta);
  sscanf(argv[3], "%le", &delta);
  printf("\n");
  printf("n = %d, delta = %3.1f%%\n", n_delta, 1e2*delta);

  err_and_corr(argv[1], delta, n_delta);
}
