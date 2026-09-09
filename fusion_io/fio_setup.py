from distutils.core import setup, Extension, os

root = os.environ['FIO_ROOT']
arch = os.environ['FIO_ARCH']
srcdir= os.environ['SRCDIR']
copt = os.environ['CFLAGS']
archflags = os.environ['ARCHFLAGS']
ldflags = os.environ.get('LDFLAGS', '')

# Split ARCHFLAGS into library search paths (-L...) and library flags (-l...).
# We need to wrap the -l flags with --whole-archive so that the linker pulls
# in ALL object files from libfusionio.a and libm3dc1.a, not just those that
# resolve undefined symbols.  This is required because the RTTI typeinfo for
# m3dc1_numpy_source (and its base classes) lives inside those static archives
# but is never referenced by name from python_interface.cpp, so the linker
# would silently drop it without --whole-archive.
#
# The resulting linker command will look like:
#   -L.../m3dc1_lib/_arch -L.../fusion_io/_arch
#   -Wl,--whole-archive -lfusionio -lm3dc1 -Wl,--no-whole-archive
#   -lhdf5 -lz ... (other non-archive flags)
#
# We also append LDFLAGS (e.g. -Wl,--export-dynamic from make.inc.stellar).

arch_parts = archflags.split()
ld_parts    = ldflags.split()

lib_search = []   # -L... flags
lib_flags  = []   # -l... flags that refer to our static archives
other_flags = []  # everything else (HDF5, -lstdc++, rpath, etc.)

# Names of the static archives we want to force-include entirely.
whole_archive_libs = {'fusionio', 'm3dc1'}

for part in arch_parts:
    if part.startswith('-L'):
        lib_search.append(part)
    elif part.startswith('-l'):
        libname = part[2:]
        if libname in whole_archive_libs:
            lib_flags.append(part)
        else:
            other_flags.append(part)
    else:
        other_flags.append(part)

# Build the extra_link_args list:
#   1. Library search paths
#   2. --whole-archive wrapper around our static archives
#   3. --no-whole-archive to stop wrapping
#   4. Everything else (HDF5, rpath, etc.)
#   5. LDFLAGS (e.g. -Wl,--export-dynamic)
extra_link = (
    lib_search
    + (['-Wl,--whole-archive'] + lib_flags + ['-Wl,--no-whole-archive']
       if lib_flags else [])
    + other_flags
    + ld_parts
)

fio_module = Extension('fio_py',
                       sources = [srcdir+'/python_interface.cpp'],
                       include_dirs = [root+'/m3dc1_lib'],
                       libraries = [],
                       library_dirs = [root+'/m3dc1_lib/_'+arch,
                                       root+'/fusion_io/_'+arch],
                       extra_link_args = extra_link)

setup (name = 'fio_py',
       author = 'N.M. Ferraro', 
       version = '1.1',
       description = 'Fusion IO Package',
       ext_modules = [fio_module],
       py_modules = ['fpy'])
