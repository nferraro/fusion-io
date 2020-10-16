!==========================================================================
! c1_lp
! ~~~~~~~~~~~
!
!
!==========================================================================
program c1_lp
  use fusion_io
  use hdf5

  implicit none

  ! Command lind arguments
  character(len=256) :: filename ! C1.h5 file
  integer            :: slice    ! time slice to read
  integer            :: ipellet  ! pellet number to read

  ! Optional namelist variables
  real :: R0, phi0, Z0  ! position to analyze from (default to pellet location)
  real :: L             ! length to traverse along field line before eval
  integer :: Nstep      ! number of steps along field line
  real :: rpol          ! radius for averaging in poloidal plane
  integer :: Npol       ! number of points to average in poloidal plane
  character(len=*), parameter :: input_nml = 'C1LP.nml'
  namelist /PARAMS/ R0, phi0, Z0, L, Nstep, rpol, Npol

  character(len=256) :: arg
  integer :: isrc, ine, ite, imag

  real :: ne_avg(1), Te_avg(1), B_avg(3)
  real :: ne(1), Te(1), B(3)

  !RK4 advance
  real :: dl
  real :: xp(3), xm(3), xt(3), k1(3), k2(3), k3(3), k4(3)

  ! parameters used in sunflower model for poloidal avg
  integer :: Bpts
  real :: r, th

  integer :: i, j, ierr, fu
  type(fio_search_hint) :: hint

  integer(HID_T) :: file_id, root_id, group_id, attr_id
  integer(HID_T) :: dset_id, space_id, mem_id
  character(len=8) :: time_group_name
  integer(HSIZE_T) :: dim1(1), dim2(2)
  integer :: itime
  integer(HSIZE_T) :: coord(2,1)
  real :: buf(1)

  ! Interpret arguments if provided
  if(iargc() .ge. 1) then
     call getarg(1, filename)
  else
     filename = 'C1.h5'
  end if
  if(iargc() .ge. 2) then
     call getarg(2, arg)
     read(arg,*,iostat=ierr) slice
  else
     slice = 0
  end if
  if(iargc() .ge. 3) then
     call getarg(3, arg)
     read(arg,*,iostat=ierr) ipellet
  else
     ipellet = 1
  end if

  ! Default values for namelist variables
  Npol = 100
  L = 0.15
  Nstep = 100

  !----------------------------------
  ! Read pellet location of C1.h5
  !----------------------------------
  call h5open_f(ierr)
  call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr)
  call h5gopen_f(file_id, "/", root_id, ierr)

  ! get time step for this slice
  write(time_group_name, '("time_",I3.3)') slice
  call h5gopen_f(root_id, time_group_name, group_id, ierr)
  call h5aopen_name_f(group_id, "ntimestep", attr_id, ierr)
  dim1 = 1
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, itime, dim1, ierr)
  call h5aclose_f(attr_id, ierr)
  call h5gclose_f(group_id, ierr)


  ! read the pellet parameters at this time
  call h5gopen_f(root_id, "pellet", group_id, ierr)
  dim2 = 1
  call h5screate_simple_f(2, dim2, mem_id, ierr)
  coord(1,1) = ipellet
  coord(2,1) = itime+1 ! itime starts at 0, but indexing from 1

  ! R0 = pellet_r
  call h5dopen_f(group_id, "pellet_r", dset_id, ierr)
  call h5dget_space_f(dset_id, space_id, ierr)
  call h5sselect_elements_f(space_id, H5S_SELECT_SET_F, 2, 1, coord, ierr)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dim2, ierr, &
       file_space_id=space_id, mem_space_id=mem_id)
  call h5sclose_f(space_id, ierr)
  call h5dclose_f(dset_id, ierr)
  R0 = buf(1)

  ! phi0 = pellet_phi
  call h5dopen_f(group_id, "pellet_phi", dset_id, ierr)
  call h5dget_space_f(dset_id, space_id, ierr)
  call h5sselect_elements_f(space_id, H5S_SELECT_SET_F, 2, 1, coord, ierr)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dim2, ierr, &
       file_space_id=space_id, mem_space_id=mem_id)
  call h5sclose_f(space_id, ierr)
  call h5dclose_f(dset_id, ierr)
  phi0 = buf(1)

  ! Z0 = pellet_z
  call h5dopen_f(group_id, "pellet_z", dset_id, ierr)
  call h5dget_space_f(dset_id, space_id, ierr)
  call h5sselect_elements_f(space_id, H5S_SELECT_SET_F, 2, 1, coord, ierr)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dim2, ierr, &
       file_space_id=space_id, mem_space_id=mem_id)
  call h5sclose_f(space_id, ierr)
  call h5dclose_f(dset_id, ierr)
  Z0 = buf(1)

  ! rpol = r_p*cloud_pel
  call h5dopen_f(group_id, "r_p", dset_id, ierr)
  call h5dget_space_f(dset_id, space_id, ierr)
  call h5sselect_elements_f(space_id, H5S_SELECT_SET_F, 2, 1, coord, ierr)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dim2, ierr, &
       file_space_id=space_id, mem_space_id=mem_id)
  call h5sclose_f(space_id, ierr)
  call h5dclose_f(dset_id, ierr)
  rpol = buf(1)
  call h5dopen_f(group_id, "cloud_pel", dset_id, ierr)
  call h5dget_space_f(dset_id, space_id, ierr)
  call h5sselect_elements_f(space_id, H5S_SELECT_SET_F, 2, 1, coord, ierr)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dim2, ierr, &
       file_space_id=space_id, mem_space_id=mem_id)
  call h5sclose_f(space_id, ierr)
  call h5dclose_f(dset_id, ierr)
  rpol = rpol*buf(1)

  call h5sclose_f(mem_id, ierr)
  call h5gclose_f(group_id, ierr)
  call h5gclose_f(root_id, ierr)
  call h5fclose_f(file_id, ierr)
  call h5close_f(ierr)

  print *, R0, phi0, Z0, rpol

  ! Open C1.h5 file
  call fio_open_source_f(FIO_M3DC1_SOURCE, trim(filename), isrc, ierr)
  if(ierr.ne.0) goto 100

  ! Set options appropriate to this source
  call fio_get_options_f(isrc, ierr)
  call fio_set_int_option_f(FIO_TIMESLICE, slice, ierr)

  ! Overwrite with namelist
  open(action='read', file=input_nml, iostat=ierr, newunit=fu)
  if(ierr.eq.0) then
     read(nml=PARAMS, iostat=ierr, unit=fu)
     close(fu)
  end if

  ! read fields
  call fio_get_field_f(isrc, FIO_MAGNETIC_FIELD, imag, ierr)
  call fio_set_int_option_f(FIO_SPECIES, FIO_ELECTRON, ierr)
  call fio_get_field_f(isrc, FIO_DENSITY, ine, ierr)
  call fio_get_field_f(isrc, FIO_TEMPERATURE, ite, ierr)

  ! initialize hint
  call fio_allocate_search_hint_f(isrc, hint, ierr)

  dl = L/Nstep
  Bpts = nint(sqrt(real(Npol)))
  ne_avg = 0.
  Te_avg = 0.
  B_avg = 0.
  do i=1, Npol

     ! sunflower arrangement to distribute starting points poloidally
     ! stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle
     if(i>(Npol-Bpts)) then
        ! sqrt(Npol) points on the boundary
        r = rpol
     else
        ! other radii distributed as such
        r = rpol*sqrt(i-0.5)/sqrt(Npol-(Bpts+1.)/2.)
     end if
     ! angle is 2*pi*i/(golden ratio)^2
     th = 8.*(4.*atan(1.))*i/(sqrt(5.)+1.)**2

     ! this position in the direction of B
     xp(1) = R0 + r*cos(th)
     xp(2) = phi0
     xp(3) = Z0 + r*sin(th)

     ! this one in the opposite direction
     xm = xp

     ! step along this field line using RK4
     ! Cylindrical coordinates so
     ! dR/ds = Br/|B|
     ! dphi/ds = Bphi/(R*|B|)
     !
     do j=1, Nstep

        ! with B
        xt = xp
        call fio_eval_field_f(imag, xt, k1, ierr, hint=hint)
        k1 = k1/norm2(k1)
        k1(2) = k1(2)/xt(1)

        xt = xp + 0.5*dl*k1
        call fio_eval_field_f(imag, xt, k2, ierr, hint=hint)
        k2 = k2/norm2(k2)
        k2(2) = k2(2)/xt(1)

        xt = xp + 0.5*dl*k2
        call fio_eval_field_f(imag, xt, k3, ierr, hint=hint)
        k3 = k3/norm2(k3)
        k3(2) = k3(2)/xt(1)

        xt = xp + dl*k3
        call fio_eval_field_f(imag, xt, k4, ierr, hint=hint)
        k4 = k4/norm2(k4)
        k4(2) = k4(2)/xt(1)

        xp = xp + dl*(k1 + 2.*k2 + 2.*k3 + k4)/6.

     end do

     do j=1, Nstep
        ! opposite B
        xt = xm
        call fio_eval_field_f(imag, xt, k1, ierr, hint=hint)
        k1 = k1/norm2(k1)
        k1(2) = k1(2)/xt(1)

        xt = xm - 0.5*dl*k1
        call fio_eval_field_f(imag, xt, k2, ierr, hint=hint)
        k2 = k2/norm2(k2)
        k2(2) = k2(2)/xt(1)

        xt = xm - 0.5*dl*k2
        call fio_eval_field_f(imag, xt, k3, ierr, hint=hint)
        k3 = k3/norm2(k3)
        k3(2) = k3(2)/xt(1)

        xt = xm - dl*k3
        call fio_eval_field_f(imag, xt, k4, ierr, hint=hint)
        k4 = k4/norm2(k4)
        k4(2) = k4(2)/xt(1)

        xm = xm - dl*(k1 + 2.*k2 + 2.*k3 + k4)/6.

     end do

     call fio_eval_field_f(ine, xp, ne, ierr, hint=hint)
     ne_avg = ne_avg + ne
     call fio_eval_field_f(ine, xm, ne, ierr, hint=hint)
     ne_avg = ne_avg + ne
     call fio_eval_field_f(ite, xp, Te, ierr, hint=hint)
     Te_avg = Te_avg + Te
     call fio_eval_field_f(ite, xm, Te, ierr, hint=hint)
     Te_avg = Te_avg + Te
     call fio_eval_field_f(imag, xp, B, ierr, hint=hint)
     B_avg = B_avg + B
     call fio_eval_field_f(imag, xm, B, ierr, hint=hint)
     B_avg = B_avg + B

  end do
  ne_avg = ne_avg/(2.*Npol)
  Te_avg = Te_avg/(2.*Npol)
  B_avg  = B_avg/(2.*Npol)
  write(*, '("        electron density = ",1pE12.4)') ne_avg
  write(*, '("        electron temperature = ",1pE12.4)') Te_avg
  write(*, '("        total B = ",1pE12.4)') norm2(B_avg)

100 continue

  call fio_deallocate_search_hint_f(isrc, hint, ierr)

  call fio_close_field_f(ine, ierr)
  call fio_close_field_f(ite, ierr)
  call fio_close_field_f(imag, ierr)
  call fio_close_source_f(isrc, ierr)

end program c1_lp
