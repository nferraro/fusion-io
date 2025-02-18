#ifndef FIO_C_INTERFACE_H
#define FIO_C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif
int fio_add_field(const int, const int, const int, const double);
int fio_allocate_search_hint(const int, void**);
int fio_close_field(const int);
int fio_close_series(const int);
int fio_close_source(const int);
int fio_create_compound_field(int*);
int fio_deallocate_search_hint(const int, void**);
int fio_eval_field(const int, const double*, double*, void*);
int fio_eval_field_deriv(const int, const double*, double*, void*);
int fio_eval_series(const int, const double, double*);
int fio_get_options(const int);
int fio_get_available_fields(const int, int*, int**);
int fio_get_field(const int, const int, int*);
int fio_get_series(const int, const int, int*);
int fio_get_series_bounds(const int, double*, double*);
int fio_open_source(const int, const char*, int*);
int fio_set_int_option(const int, const int);
int fio_set_str_option(const int, const char*);
int fio_set_real_option(const int, const double);
int fio_get_int_parameter(const int, const int,  int*);
int fio_get_real_parameter(const int, const int, double*);
int fio_get_real_field_parameter(const int, const int, double*);
int fio_set_quiet_option(const bool);
int fio_sizeof_search_hint(const int);
#ifdef __cplusplus
}
#endif

#endif
