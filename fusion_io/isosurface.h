#ifndef FIO_ISOSURFACE_H
#define FIO_ISOSURFACE_H

int fio_find_val_2d(fio_field* f, const double val, double* x,
		    const double tol, const double max_step, fio_hint h=0);
int fio_isosurface(fio_field* f, const double val, const double* guess,
		   const double* enclose, const double dl, const double tol,
		   const double max_step, const double dphi, 
		   int* n, double*** path, fio_hint h=0);
int fio_gridify_surface(const int m0, double** path0, const double* axis, 
			const int nphi, const int ntheta, 
			double** path, double* phi, double* theta);

#endif
