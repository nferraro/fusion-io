import fio_py
import sys

import faulthandler
faulthandler.enable()

# Example script for using NIMROD interface in fusion-io
# Usage: python fusionio_example_python_nimrod.py <dumpfile>

if len(sys.argv) < 2:
    print("Usage: python fusionio_example_python_nimrod.py <dumpfile>")
    sys.exit(1)

filename = sys.argv[1]

# Open NIMROD source
isrc = fio_py.open_source(fio_py.FIO_NIMROD_SOURCE, filename)
fio_py.get_options(isrc)

# Get available fields
fields = fio_py.get_available_fields(isrc)
print('Available fields:')
for f in fields:
    print(' ', fio_py.get_field_name(f))

# Get fields
imag = fio_py.get_field(isrc, fio_py.FIO_MAGNETIC_FIELD)
ivel = fio_py.get_field(isrc, fio_py.FIO_FLUID_VELOCITY)
ipres = fio_py.get_field(isrc, fio_py.FIO_PRESSURE)
idens = fio_py.get_field(isrc, fio_py.FIO_DENSITY)

# Evaluate at a point (R, Phi, Z)
x = (1.5, 0.0, 0.0)
hint = fio_py.allocate_hint(isrc)

print(f'\nEvaluating at {x}:')

try:
    br, bphi, bz = fio_py.eval_vector_field(imag, x, hint)
    print(f' Magnetic Field (R, Phi, Z): ({br}, {bphi}, {bz})')
except Exception as e:
    print(f' Error evaluating magnetic field: {e}')

try:
    vr, vphi, vz = fio_py.eval_vector_field(ivel, x, hint)
    print(f' Fluid Velocity (R, Phi, Z): ({vr}, {vphi}, {vz})')
except Exception as e:
    print(f' Error evaluating fluid velocity: {e}')

try:
    p = fio_py.eval_scalar_field(ipres, x, hint)
    print(f' Pressure: {p}')
except Exception as e:
    print(f' Error evaluating pressure: {e}')

try:
    n = fio_py.eval_scalar_field(idens, x, hint)
    print(f' Density: {n}')
except Exception as e:
    print(f' Error evaluating density: {e}')

# Evaluate derivatives
try:
    db = fio_py.eval_vector_field_deriv(imag, x, hint)
    print('\n Magnetic Field Derivatives:')
    print(f'  dB/dR: {db[0:3]}')
    print(f'  dB/dPhi: {db[3:6]}')
    print(f'  dB/dZ: {db[6:9]}')
except Exception as e:
    print(f' Error evaluating derivatives: {e}')

fio_py.free_hint(hint)
fio_py.close_field(imag)
fio_py.close_field(ivel)
fio_py.close_field(ipres)
fio_py.close_field(idens)
fio_py.close_source(isrc)
