#ifndef FUSION_IO_FIELD_H
#define FUSION_IO_FIELD_H

#include "fusion_io_defs.h"
#include <vector>

typedef int field_attribute;
typedef int field_parameter;
typedef void* fio_hint;

class fio_field {
 public:
  virtual ~fio_field()
    { }

  virtual fio_field* clone() const = 0;
  virtual int dimension() const = 0;
  virtual int eval(const double*, double*, fio_hint =0) = 0;

  // partial derivatives (first index = partial deriv, second index = coord)
  virtual int eval_deriv(const double*, double*, fio_hint =0)
  {  return FIO_UNSUPPORTED; }

  virtual int get_real_parameter(const field_parameter, double*)
  {  return FIO_UNSUPPORTED; }

  fio_field& operator+(const fio_field&);
  fio_field& operator*(const fio_field&);

  int find_val_on_line(const double, const double*, const double*, 
		       double*, fio_hint =0, const double=1e-4);
};

#endif
