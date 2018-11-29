#ifndef MARS_SOURCE_H
#define MARS_SOURCE_H

#include "fusion_io_source.h"

class mars_source : public fio_source {
 public:
  int nrp1, nchi;
  double aspect, r0exp, b0exp, bfac, vfac;
  double* psiiso;

  struct mars_eqdata {
  private:
    int nr, nc;

  public:
    double* cse;     // radial coordinate s
    double* peq;
    double* t;
    double* ttp;
    double* ppeq;
    double* dpsids;
    double** req;
    double** zeq;
    double** rja;
    double** g11l;
    double** g22l;
    double** g33l;
    double** g12l;
    double** rdcdz;
    double** rdsdz;
    double** rbz;
    double* tp;

    mars_eqdata(const int, const int);
    ~mars_eqdata();
  } *eqdata, *eqdatam;

  struct mars_eigdata {
    int maxm, nrp;
    bool extra;
    
    double* dpsids[3];
    double* t[3];
    double* mnum_r[3];
    double* mnum_i[3];
    double** val_r[3];
    double** val_i[3];

    mars_eigdata(const int, const int, const bool);
    ~mars_eigdata();
  } *xplasma, *vplasma, *bplasma;

  int load_eigdata(const char*, mars_eigdata**, const bool ex);

 public:
  mars_source();
  virtual ~mars_source();

  int open(const char*);
  int close();

  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);
  int get_series(const series_type, fio_series**);
  int get_available_fields(fio_field_list*) const;
};


#endif 
