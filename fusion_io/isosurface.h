#ifndef FIO_ISOSURFACE_H
#define FIO_ISOSURFACE_H

int fio_find_val(fio_field* f, const double val, double* x,
		 const double tol, const double max_step, 
		 const int dim, const double* axis, fio_hint h=0);
int fio_isosurface(fio_field* f, const double val, const double* guess,
		   const double* enclose, const double dl, const double tol,
		   const double max_step, const double dphi, 
		   int* n, double*** path, fio_hint h=0);
int fio_isosurface_2d(fio_field* f, const double val, const double* guess,
		      const double* axis, const double* norm, 
		      const double dl, const double tol,
		      const double max_step,
		      int* n, double*** path, fio_hint h=0);
int fio_gridify_loop(const int m0, double** path0, const double* axis, 
		     const int m, double** path, double* theta, 
		     const int param=0);
int fio_geom_isosurface(fio_field* f, const double val, const double* guess,
			   const double* axis, const double tol,
			   const double max_step,
			   const double nphi, const double ntheta, 
			   double*** path, fio_hint h=0);
int fio_gridify_surface(const int m0, double** path0, const double* axis, 
			const int nphi, const int ntheta, 
			double** path, double* phi, double* theta);
int fio_gridded_isosurface(fio_field* f, const double val, const double* guess,
			   const double* axis, 
			   const double dl_tor, const double dl_pol, 
			   const double tol, const double max_step,
			   const int nphi, const int ntheta, 
			   double* phi, double* theta,
			   double** path, const char* label, fio_hint h=0);

int fio_q_at_surface(fio_field* f, const int n, double** x, double* q,
		     double* bpol, fio_hint h=0);
int fio_surface_average(fio_field* f, const int n, double** x, double* a,
			double* bpol, fio_hint h=0);


#endif
