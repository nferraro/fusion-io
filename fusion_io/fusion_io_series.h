#ifndef FIO_SERIES
#define FIO_SERIES


class fio_series {
 public:
  virtual ~fio_series()
    { }

  virtual int eval(const double, double*) = 0;
  virtual int bounds(double*, double*) const = 0;
};


class fio_scalar_series : public fio_series {
  double data;

 public:
  fio_scalar_series(const double d) { data = d; }
  virtual ~fio_scalar_series()
    { }

  virtual int eval(const double, double*);
  virtual int bounds(double*, double*) const;
};


class fio_array_series : public fio_series {
  double* data;
  int n;

 public:
  fio_array_series(const int, const double*);
  virtual ~fio_array_series();

  virtual int eval(const double, double*);
  virtual int bounds(double*, double*) const;
};


#endif
