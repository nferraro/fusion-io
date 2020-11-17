#ifndef TRACE_INTEGRATOR_S_H
#define TRACE_INTEGRATOR_S_H

#include <fusion_io_source.h>

#include <fstream>
#include <deque>

//#include "trace_source.h"
//#include "interpolation_source.h"

struct trace_source {
  void* hint;
  fio_source* source;
  fio_field* field;
  fio_field* rst;
  fio_field* zst;
  fio_field* psi_norm;
  double magaxis[2];
  trace_source() {
    source = 0;
    field = 0;
    rst = 0;
    zst = 0;
    psi_norm = 0;
    hint = 0;
    magaxis[0] = 0.;
    magaxis[1] = 0.;
  }
  void free() {
    if(hint) {
      if(source) source->deallocate_search_hint(&hint);
      hint = 0;
    }
    if(psi_norm) {
      delete(psi_norm);
      psi_norm = 0;
    }
    if(field) {
      delete(field);
      field = 0;
    }
    if(rst) {
      delete(rst);
      rst = 0;
    }
    if(zst) {
      delete(zst);
      zst = 0;
    }
    if(source) {
      delete(source);
      source = 0;
    }
  }
};

typedef std::deque<trace_source> trace_source_list;

class trace_integrator_s {
  double R, Phi, Z;

  double dr[4], dz[4];
  bool reverse;

  std::fstream file;

 public:

  bool toroidal;
  double period;
  int nplanes;
  int tpts;

  struct integrator_data {
    int toroidal_transits;
    int poloidal_transits;
    double q;
    double distance;
  };

  double plane;
  trace_source_list sources;

  trace_integrator_s();
  virtual ~trace_integrator_s();

  bool load();
  bool eval(const double r, const double phi, const double z, 
	    double* br, double* bphi, double* bz);
  bool eval_psi(const double r, const double phi, const double z, 
	    double* psi);
  bool rhs(const double r, const double phi, const double z, 
	    double* dx, double* dy);
  bool center(double* r0, double* z0) const;
  bool get_surface(const double r0, const double phi0, const double z0,
		   const double ds, double** r, double** z, int* n);

  void set_reverse(const bool r);

  bool set_pos(const double r, const double phi, const double z);
  void get_pos(double *r, double* phi, double* z) const
  { *r = R; *phi = Phi; *z = Z; }
  bool step_euler(double dphi);
  bool step_rk3(double dphi);
  bool step_rk4(double dphi);
  bool step_predcorr(double dphi);
  bool integrate(int transits, int steps_per_transit, 
		 integrator_data* data=0);

  void plot_point();

  bool open_file(const char*);
  bool close_file();
};

#endif
