#ifndef NIMROD_FIELD_H
#define NIMROD_FIELD_H

#ifdef FUSIONIO_ENABLE_NIMROD

#include "fusion_io_field.h"
#include <string>

class nimrod_field : public fio_field {
  std::string name;
  int nqty;
  int eq;

 public:
  nimrod_field(const std::string& n, int q, int e)
    : name(n), nqty(q), eq(e) { }
  virtual ~nimrod_field() { }

  fio_field* clone() const { return new nimrod_field(name, nqty, eq); }
  int dimension() const { return nqty; }
  int eval(const double*, double*, fio_hint = 0);
  int eval_deriv(const double*, double*, fio_hint = 0);
};

#endif // FUSIONIO_ENABLE_NIMROD
#endif // NIMROD_FIELD_H
