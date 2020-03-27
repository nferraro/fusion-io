module fusion_io
  use iso_c_binding
  implicit none

#include "fusion_io_defs.h"

  type fio_search_hint
    type(c_ptr) :: ptr
  end type fio_search_hint

  interface
    integer(c_int) function fio_add_field(icfield, ifield, op, fac) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: icfield, ifield, op
      real(c_double), intent(in), value :: fac
    end function fio_add_field

    integer(c_int) function fio_allocate_search_hint(isrc, s) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc
      type(c_ptr) :: s
    end function fio_allocate_search_hint

    integer(c_int) function fio_close_field(ifield) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: ifield
    end function fio_close_field

    integer(c_int) function fio_close_series(iseries) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: iseries
    end function fio_close_series

    integer(c_int) function fio_close_source(ifield) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: ifield
    end function fio_close_source

    integer(c_int) function fio_create_compound_field(ifield) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: ifield
    end function fio_create_compound_field

    integer(c_int) function fio_deallocate_search_hint(isrc, s) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc
      type(c_ptr) :: s
    end function fio_deallocate_search_hint

    integer(c_int) function fio_eval_field(ifield, x, v, s) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: ifield
      real(c_double), intent(in), dimension(*) :: x
      real(c_double), dimension(*) :: v
      type(c_ptr), value :: s
    end function fio_eval_field

    integer(c_int) function fio_eval_field_deriv(ifield, x, v, s) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: ifield
      real(c_double), intent(in), dimension(*) :: x
      real(c_double), dimension(*) :: v
      type(c_ptr), value :: s
    end function fio_eval_field_deriv

    integer(c_int) function fio_eval_series(iseries, x, v) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: iseries
      real(c_double), intent(in), value :: x
      real(c_double) :: v
    end function fio_eval_series

    integer(c_int) function fio_get_available_fields(isrc, n, f) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc
      integer(c_int) :: n
      type(c_ptr) :: f
    end function fio_get_available_fields

    integer(c_int) function fio_get_field(isrc, type, handle) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc, type
      integer(c_int) :: handle
    end function fio_get_field

    integer(c_int) function fio_get_int_parameter(isrc, t, p) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc, t
      integer(c_int) :: p
    end function fio_get_int_parameter

    integer(c_int) function fio_get_options(isrc) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc
    end function fio_get_options

    integer(c_int) function fio_get_real_parameter(isrc, t, p) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc, t
      real(c_double) :: p
    end function fio_get_real_parameter

    integer(c_int) function fio_get_real_field_parameter(ifield, t, p) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: ifield, t
      real(c_double) :: p
    end function fio_get_real_field_parameter

    integer(c_int) function fio_get_series(isrc, type, handle) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: isrc, type
      integer(c_int) :: handle
    end function fio_get_series

    integer(c_int) function fio_get_series_bounds(iseries, tmin, tmax) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: iseries
      real(c_double) :: tmin, tmax
    end function fio_get_series_bounds

    integer(c_int) function fio_open_source(itype, filename, handle) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: itype
      character(len=1,kind=c_char), dimension(*), intent(in) :: filename
      integer(c_int) :: handle
    end function fio_open_source

    integer(c_int) function fio_set_int_option(iopt, v) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: iopt, v
    end function fio_set_int_option

    integer(c_int) function fio_set_real_option(iopt, v) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: iopt
      real(c_double), intent(in), value :: v
    end function fio_set_real_option

    integer(c_int) function fio_set_str_option(iopt, v) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: iopt
      character(len=1,kind=c_char), dimension(*), intent(in) :: v
    end function fio_set_str_option
  end interface

contains

  subroutine fio_add_field_f(icfield, ifield, iop, fac, ierr)
    integer, intent(in) :: icfield
    integer, intent(in) :: ifield
    integer, intent(in) :: iop
    real, intent(in) :: fac
    integer, intent(out) :: ierr
    ierr = fio_add_field(icfield, ifield, iop, fac)
  end subroutine fio_add_field_f

  subroutine fio_allocate_search_hint_f(isrc, hint, ierr)
    integer, intent(in) :: isrc
    type(fio_search_hint), intent(out) :: hint
    integer, intent(out) :: ierr
    ierr = fio_allocate_search_hint(isrc, hint%ptr)
  end subroutine fio_allocate_search_hint_f

  subroutine fio_close_field_f(ifield, ierr)
    integer, intent(in) :: ifield
    integer, intent(out) :: ierr
    ierr = fio_close_field(ifield)
  end subroutine fio_close_field_f

  subroutine fio_close_series_f(iseries, ierr)
    integer, intent(in) :: iseries
    integer, intent(out) :: ierr
    ierr = fio_close_series(iseries)
  end subroutine fio_close_series_f

  subroutine fio_close_source_f(isrc, ierr)
    integer, intent(in) :: isrc
    integer, intent(out) :: ierr
    ierr = fio_close_source(isrc)
  end subroutine fio_close_source_f

  subroutine fio_create_compound_field_f(ifield, ierr)
    integer, intent(out) :: ifield
    integer, intent(out) :: ierr
    ierr = fio_create_compound_field(ifield)
  end subroutine fio_create_compound_field_f

  subroutine fio_deallocate_search_hint_f(isrc, hint, ierr)
    integer, intent(in) :: isrc
    type(fio_search_hint), intent(inout) :: hint
    integer, intent(out) :: ierr
    ierr = fio_deallocate_search_hint(isrc, hint%ptr)
  end subroutine fio_deallocate_search_hint_f

  subroutine fio_eval_field_f(ifield, x, v, ierr, hint)
    integer, intent(in) :: ifield
    real, intent(in), dimension(*) :: x
    real, intent(out), dimension(*) :: v
    integer, intent(out) :: ierr
    type(fio_search_hint), intent(inout), optional :: hint
    if(present(hint)) then
      ierr = fio_eval_field(ifield, x, v, hint%ptr)
    else
      ierr = fio_eval_field(ifield, x, v, c_null_ptr)
    endif
  end subroutine fio_eval_field_f

  subroutine fio_eval_field_deriv_f(ifield, x, v, ierr, hint)
    integer, intent(in) :: ifield
    real, intent(in), dimension(*) :: x
    real, intent(out), dimension(*) :: v
    integer, intent(out) :: ierr
    type(fio_search_hint), intent(inout), optional :: hint
    if(present(hint)) then
      ierr = fio_eval_field_deriv(ifield, x, v, hint%ptr)
    else
      ierr = fio_eval_field_deriv(ifield, x, v, c_null_ptr)
    endif
  end subroutine fio_eval_field_deriv_f

  subroutine fio_eval_series_f(iseries, x, v, ierr)
    integer, intent(in) :: iseries
    real, intent(in) :: x
    real, intent(out) :: v
    integer, intent(out) :: ierr
    ierr = fio_eval_series(iseries, x, v)
  end subroutine fio_eval_series_f

  subroutine fio_get_field_f(isrc, itype, ifield, ierr)
    integer, intent(in) :: isrc, itype
    integer, intent(out) :: ifield, ierr
    ierr = fio_get_field(isrc, itype, ifield)
  end subroutine fio_get_field_f

  subroutine fio_get_options_f(isrc, ierr)
    integer, intent(in) :: isrc
    integer, intent(out) :: ierr
    ierr = fio_get_options(isrc)
  end subroutine fio_get_options_f

  subroutine fio_get_int_parameter_f(isrc, t, p, ierr)
    integer, intent(in) :: isrc, t
    integer, intent(out) :: p
    integer, intent(out) :: ierr
    ierr = fio_get_int_parameter(isrc, t, p)
  end subroutine fio_get_int_parameter_f

  subroutine fio_get_real_parameter_f(isrc, t, p, ierr)
    integer, intent(in) :: isrc, t
    real, intent(out) :: p
    integer, intent(out) :: ierr
    ierr = fio_get_real_parameter(isrc, t, p)
  end subroutine fio_get_real_parameter_f

  subroutine fio_get_real_field_parameter_f(ifield, t, p, ierr)
    integer, intent(in) :: ifield, t
    real, intent(out) :: p
    integer, intent(out) :: ierr
    ierr = fio_get_real_field_parameter(ifield, t, p)
  end subroutine fio_get_real_field_parameter_f

  subroutine fio_get_series_f(isrc, itype, iseries, ierr)
    integer, intent(in) :: isrc, itype
    integer, intent(out) :: iseries, ierr
    ierr = fio_get_series(isrc, itype, iseries)
  end subroutine fio_get_series_f

  subroutine fio_get_series_bounds_f(iseries, tmin, tmax, ierr)
    integer, intent(in) :: iseries
    real, intent(out) :: tmin, tmax
    integer, intent(out) :: ierr
    ierr = fio_get_series_bounds(iseries, tmin, tmax)
  end subroutine fio_get_series_bounds_f

  subroutine fio_open_source_f(type, filename, isrc, ierr)
    integer, intent(in) :: type
    character(len=*), intent(in) :: filename
    integer, intent(out) :: isrc, ierr
    ierr = fio_open_source(type, filename//c_null_char, isrc)
  end subroutine fio_open_source_f

  subroutine fio_set_int_option_f(iopt, val, ierr)
    integer, intent(in) :: iopt
    integer, intent(in) :: val
    integer, intent(out) :: ierr
    ierr = fio_set_int_option(iopt, val)
  end subroutine fio_set_int_option_f

  subroutine fio_set_real_option_f(iopt, val, ierr)
    integer, intent(in) :: iopt
    real, intent(in) :: val
    integer, intent(out) :: ierr
    ierr = fio_set_real_option(iopt, val)
  end subroutine fio_set_real_option_f

  subroutine fio_set_str_option_f(iopt, val, ierr)
    integer, intent(in) :: iopt
    character(len=*), intent(in) :: val
    integer, intent(out) :: ierr
    ierr = fio_set_str_option(iopt, val//c_null_char)
  end subroutine fio_set_str_option_f
end module fusion_io
