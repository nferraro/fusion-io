from distutils.core import setup, Extension, os

root = os.environ['FIO_ROOT']
arch = os.environ['FIO_ARCH']
srcdir= os.environ['SRCDIR']
copt = os.environ['CFLAGS']
archflags = os.environ['ARCHFLAGS']

fio_module = Extension('fio_py',
                       sources = [srcdir+'/python_interface.cpp'],
                       include_dirs = [root+'/m3dc1_lib'],
                       libraries = ['fusionio', 'm3dc1'],
                       library_dirs = [root+'/m3dc1_lib/_'+arch,
                                       root+'/fusion_io/_'+arch],
                       extra_link_args = [archflags])

setup (name = 'fio_py',
       author = 'N.M. Ferraro', 
       version = '1.0',
       description = 'Fusion IO Package',
       ext_modules = [fio_module],
       py_modules = ['fpy'])
