import fio_py

filename = "C1.h5"

isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)

fio_py.get_options(isrc)
ia = fio_py.get_field(isrc, fio_py.FIO_VECTOR_POTENTIAL)
imag = fio_py.get_field(isrc, fio_py.FIO_MAGNETIC_FIELD)
ipres = fio_py.get_field(isrc, fio_py.FIO_TOTAL_PRESSURE)
ivel = fio_py.get_field(isrc, fio_py.FIO_FLUID_VELOCITY)
ij = fio_py.get_field(isrc, fio_py.FIO_CURRENT_DENSITY)

fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_ELECTRON)
idens = fio_py.get_field(isrc, fio_py.FIO_DENSITY)

fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_MAIN_ION)
ipi = fio_py.get_field(isrc, fio_py.FIO_PRESSURE)

list = fio_py.get_available_fields(isrc)

print('Available fields:')
for field in list:
    print(' ', fio_py.get_field_name(field))

cs = fio_py.get_int_parameter(isrc, fio_py.FIO_GEOMETRY)
if cs == fio_py.FIO_CYLINDRICAL:
    print('Using CYLINDRICAL coordinate system')
else :
    print('Using CARTESIAN coordinate system')

ntime = fio_py.get_int_parameter(isrc, fio_py.FIO_NUM_TIMESLICES)
print('Number of time slices: ', ntime)

period = fio_py.get_real_parameter(isrc, fio_py.FIO_PERIOD)
print('Toroidal period = ', period)

ipsi_axis = fio_py.get_series(isrc, fio_py.FIO_MAGAXIS_PSI)
ipsi_lcfs = fio_py.get_series(isrc, fio_py.FIO_LCFS_PSI)

psi_axis = fio_py.eval_series(ipsi_axis, 0.)
psi_lcfs = fio_py.eval_series(ipsi_lcfs, 0.)
print('Psi at magnetic axis: ', psi_axis)
print('Psi at lcfs: ', psi_lcfs)

fio_py.close_series(ipsi_axis)
fio_py.close_series(ipsi_lcfs)

hint = fio_py.allocate_hint(isrc);

x = (1.6, 0., 0.)
(ar, aphi, az) = fio_py.eval_vector_field(ia, x, hint)
(br, bphi, bz) = fio_py.eval_vector_field(imag, x, hint)
(jr, jphi, jz) = fio_py.eval_vector_field(ij, x, hint)
p = fio_py.eval_scalar_field(ipres, x, hint)
pi = fio_py.eval_scalar_field(ipi, x, hint)
ne = fio_py.eval_scalar_field(idens, x, hint)
(vr, vphi, vz) = fio_py.eval_vector_field(ivel, x, hint)
psi = aphi*x[0]
psi_norm = (psi - psi_axis) / (psi_lcfs - psi_axis)

fio_py.free_hint(hint)

print('At x = ', x)
print(' poloidal flux = ', psi, ' Wb')
print(' normalized poloidal flux = ', psi_norm)
print(' vector potential = ', (ar, aphi, az), ' Wb/m')
print(' magnetic field = ', (br, bphi, bz), ' T')
print(' current density = ', (jr, jphi, jz), ' A/m^2')
print(' fluid velocity = ', (vr, vphi, vz), ' m/s')
print(' toroidal rotation = ', vphi/x[0], ' rad / s' )
print(' pressure = ', p, ' Pa')
print(' ion pressure = ', pi, ' Pa')
print(' density = ', ne, ' ')

fio_py.close_field(imag)
fio_py.close_field(ipres)
fio_py.close_source(isrc)

