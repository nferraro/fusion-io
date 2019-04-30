!==========================================================================
! fio_example
! ~~~~~~~~~~~
! In this example, several files from linear calculations are loaded.
! The fields are evaluated at several points, and the results are summed.
!==========================================================================
program fio_example
  use fusion_io

  implicit none

  ! factor by which to multiply perturbed (linear) part of solution
  real, parameter :: factor = 1.

  ! time slice to read
  integer, parameter :: timeslice = 1

  character(len=256) :: filename_m3dc1, filename_efit, filename_gato, filename_mars
  integer :: isrc, ipres, ine, ini, imag, ij
  integer :: isrc_efit, imag_efit, ipres_efit, ij_efit
  integer :: isrc_gato, imag_gato, ipres_gato
  integer :: isrc_mars, imag_mars
  integer :: ipsi_axis, ipsi_lcfs
  real :: psi_axis, psi_lcfs

  real :: p(1), ne(1), ni(1), b(3), x(3), curr(3), curr2(3), db(9)

  integer, parameter :: ifile = 22
  integer, parameter :: npts = 40
  real :: R0, R1, Z0, Z1, phi0, phi1, period, tmin, tmax, time
  integer :: i, j, cs, ntime, ierr
  type(fio_search_hint) :: hint

  filename_m3dc1 = 'C1.h5'

  ! read files and fields
  print *, 'Reading ', trim(filename_m3dc1)
  call fio_open_source_f(FIO_M3DC1_SOURCE, trim(filename_m3dc1), isrc, ierr)
  if(ierr.ne.0) goto 100

  ! Print information about coordinate system
  call fio_get_int_parameter_f(isrc, FIO_GEOMETRY, cs, ierr)
  call fio_get_real_parameter_f(isrc, FIO_PERIOD, period, ierr)
  if(cs.eq.FIO_CARTESIAN) then
     print *, 'Using CARTESIAN coordinate system'
     print *, 'Toroidal period (in m) = ', period
  else if(cs.eq.FIO_CYLINDRICAL) then
     print *, 'Using CYLINDRICAL coordinate system'
     print *, 'Toroidal period (in rad) = ', period
  else 
     print *, 'ERROR: Unrecognized coordinate system'
  end if

  call fio_get_int_parameter_f(isrc, FIO_NUM_TIMESLICES, ntime, ierr)
  print *, 'Number of timeslices: ', ntime

  ! Set options appropriate to this source
  call fio_get_options_f(isrc, ierr)
  call fio_set_int_option_f(FIO_TIMESLICE, timeslice, ierr)
  call fio_set_real_option_f(FIO_LINEAR_SCALE, factor, ierr)

  ! By default, full fields (equilibrium + perturbed) will be read
  ! to read only the perturbed fields, use
  ! call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)

  ! read fields
  ! magnetic field and total pressure are species-independent
  call fio_get_field_f(isrc, FIO_TOTAL_PRESSURE, ipres, ierr);
  call fio_get_field_f(isrc, FIO_MAGNETIC_FIELD, imag, ierr);
  call fio_get_field_f(isrc, FIO_CURRENT_DENSITY, ij, ierr);

  ! density is species dependent; specify electrons
  call fio_set_int_option_f(FIO_SPECIES, FIO_ELECTRON, ierr)
  call fio_get_field_f(isrc, FIO_DENSITY, ine,ierr);
  
  call fio_set_int_option_f(FIO_SPECIES, FIO_MAIN_ION, ierr)
  call fio_get_field_f(isrc, FIO_DENSITY, ini,ierr);

  ! Determine time of timeslice
  call fio_get_real_field_parameter_f(ini, FIO_TIME, time, ierr);
  print *, 'Density is being evaluated at t = ', time


  ! Evaluate the flux at the magnetic axis and lcfs
  call fio_get_series_f(isrc, FIO_MAGAXIS_PSI, ipsi_axis, ierr)
  call fio_get_series_f(isrc, FIO_LCFS_PSI, ipsi_lcfs, ierr)

  call fio_eval_series_f(ipsi_axis, 0., psi_axis, ierr)
  call fio_eval_series_f(ipsi_lcfs, 0., psi_lcfs, ierr)
  print *, 'Psi at magnetic axis: ', psi_axis
  print *, 'Psi at lcfs: ', psi_lcfs

  call fio_get_series_bounds_f(ipsi_axis, tmin, tmax, ierr)
  print *, 'Time domain bounds for magnetic axis data: ', tmin, tmax

  call fio_close_series_f(ipsi_axis, ierr)
  call fio_close_series_f(ipsi_lcfs, ierr)

  
!  ! open efit file
!  filename_efit = &
!       '/p/tsc/nferraro/data/DIII-D/126006/3600_efit06/orlov/mesh21a_kap6_amu6_n=3/geqdsk'
!  call fio_open_source_f(FIO_GEQDSK_SOURCE,trim(filename_efit),isrc_efit,ierr)
!  if(ierr.ne.0) goto 100

!  call fio_get_options_f(isrc_efit, ierr)
!  call fio_get_field_f(isrc_efit, FIO_MAGNETIC_FIELD, imag_efit, ierr);
!  call fio_get_field_f(isrc_efit, FIO_TOTAL_PRESSURE, ipres_efit, ierr);
!  call fio_get_field_f(isrc_efit, FIO_CURRENT_DENSITY, ij_efit, ierr);

!!$  ! open gato file
!!$  filename_gato = '/Users/ferraro/data/GATO/diagnostics.dat'
!!$  call fio_open_source_f(FIO_GATO_SOURCE,trim(filename_gato),isrc_gato,ierr)
!!$  if(ierr.ne.0) goto 100
!!$  call fio_get_options_f(isrc_gato, ierr)
!!$  call fio_get_field_f(isrc_gato, FIO_MAGNETIC_FIELD, imag_gato, ierr);
!!$  call fio_get_field_f(isrc_gato, FIO_TOTAL_PRESSURE, ipres_gato, ierr);

!!$  ! open mars file
!!$  filename_mars = ''
!!$  call fio_open_source_f(FIO_MARS_SOURCE,trim(filename_mars),isrc_mars,ierr)
!!$  if(ierr.ne.FIO_SUCCESS) then
!!$     print *, 'Error reading MARS source'
!!$     goto 100
!!$  end if
!!$
!!$  call fio_get_options_f(isrc_mars, ierr)
!!$  call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)
!!$
!!$  call fio_get_field_f(isrc_mars, FIO_MAGNETIC_FIELD, imag_mars, ierr);
!!$  if(ierr.ne.FIO_SUCCESS) then
!!$     print *, 'Error reading MARS magnetic field'
!!$     goto 100
!!$  end if 

!!$  x(1) = 1.02
!!$  x(2) = 0.
!!$  x(3) = 0.05
!!$  print *, 'Evaluating field'
!!$  call fio_eval_field_f(imag_mars, x, b, ierr)
!!$  print *, 'Done. ierr = ', ierr

  ! initialize hint
  call fio_allocate_search_hint_f(isrc, hint, ierr)

  R0 = 0.4;
  R1 = 2.4;
  Z0 = -1.0;
  Z1 =  1.0;
  phi0 = 0.;
  phi1 = 0.;

!!$  open(ifile, file='b.fio', action='write')

  do i=1, npts
     do j=1, npts
        x(1) = R0 + (R1-R0)*(i-1)/(npts-1);
        x(2) = 0.
        x(3) = Z0 + (Z1-Z0)*(j-1)/(npts-1);

!        write(*, '("(",3F12.4,"):")') x

!!$        call fio_eval_field_f(imag_mars, x, b, ierr)
!        write(*, '("        efit b = ", 1p3E12.4)') b
!!$     
!!$        write(ifile, '(6e14.5)') x(1), x(3), x(2), b(1), b(3), b(2)
!!$
     call fio_eval_field_f(ipres, x, p, ierr, hint=hint)
     call fio_eval_field_f(ine, x, ne, ierr, hint=hint)
     call fio_eval_field_f(ini, x, ni, ierr, hint=hint)
     call fio_eval_field_f(imag, x, b, ierr, hint=hint)
     call fio_eval_field_f(ij, x, curr, ierr, hint=hint)
     call fio_eval_field_deriv_f(imag, x, db, ierr, hint=hint)

     write(*, '("        pressure = ",1pE12.4)') p
     write(*, '("        electron density = ",1pE12.4)') ne
     write(*, '("        ion density = ",1pE12.4)') ni
     write(*, '("        total B = ",1p3E12.4)') b
     write(*, '("        total J = ",1p3E12.4)') curr
     write(*, '("        dB/dR = ",1p3E12.4)') &
          db(FIO_DR_R), db(FIO_DR_PHI), db(FIO_DR_Z)
     write(*, '("        dB/dPhi = ",1p3E12.4)') &
          db(FIO_DPHI_R), db(FIO_DPHI_PHI), db(FIO_DPHI_Z)
     write(*, '("        dB/dZ = ",1p3E12.4)') &
          db(FIO_DZ_R), db(FIO_DZ_PHI), db(FIO_DZ_Z)

!     b = 0.
!     p = 0.
!     curr2 = 0.
!     call fio_eval_field_f(ipres_efit, x, p, ierr)
!     call fio_eval_field_f(imag_efit, x, b, ierr)
!     call fio_eval_field_f(ij_efit, x, curr2, ierr)
!     write(*, '("        efit press = ", 1pE12.4)') p
!     write(*, '("        efit b = ", 1p3E12.4)') b
!     write(*, '("        efit j = ", 1p3E12.4)') curr2

!!$     call fio_eval_field_f(ipres_gato, x, p, ierr)
!!$     call fio_eval_field_f(imag_gato, x, b, ierr)
!!$     write(*, '("        gato press = ", 1pE12.4)') p
!!$     write(*, '("        gato b = ", 1p3E12.4)') b
     end do

!!$     write(ifile, *)
  end do

!!$  close(ifile)

100 continue 

  call fio_deallocate_search_hint_f(isrc, hint, ierr)

  call fio_close_field_f(ine, ierr)
  call fio_close_field_f(ini, ierr)
  call fio_close_field_f(ipres, ierr)
  call fio_close_field_f(imag, ierr)
  call fio_close_field_f(ij, ierr)
  call fio_close_source_f(isrc, ierr)

!!$  call fio_close_field_f(imag_mars, ierr)
!!$  call fio_close_source_f(isrc_mars, ierr)


!  call fio_close_field_f(imag_efit, ierr)
!  call fio_close_field_f(ipres_efit, ierr)
!  call fio_close_field_f(ij_efit, ierr)
!  call fio_close_source_f(isrc_efit, ierr)

!!$  call fio_close_field_f(imag_gato, ierr)
!!$  call fio_close_field_f(ipres_gato, ierr)
!!$  call fio_close_source_f(isrc_gato, ierr)
end program fio_example
