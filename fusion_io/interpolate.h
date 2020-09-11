#ifndef INTERPOLATE_H
#define INTERPOLATE_H

bool bicubic_interpolation_coeffs(const double** x, const int m, const int n,
				  const int i, const int j, double a[4][4]);
bool cubic_interpolation_coeffs(const double* x, const int m, const int i,
				double* a);
bool cubic_interpolation_coeffs_cyc(const double* x, const int m, const int i,
				double* a);
bool cubic_interpolation(const int m, const double* p, const double p0, 
			 const double* f, double* f0, bool cyc=false);
bool shift_array(const int m, double* x, const int s);
bool reverse_array(const int m, double* x);


#endif
