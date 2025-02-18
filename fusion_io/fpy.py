# fpy.py: wrapper functions for the sim_data library to make accessing
# simulation output more pythonic.
#
#
# coded by Christopher Berg Smiet on 18 January 2019
# Edited by Ralf Mackenbach on the 14th of August 2019
# Mesh added by Andreas Kleiner on 20 August 2019
# Edited by Brendan C. Lyons on 9 December 2020
# csmiet@pppl.gov
# rmackenb@pppl.gov
# akleiner@pppl.gov
# lyonsbc@fusion.gat.com

import fio_py
import numpy as np
import h5py
import os
try:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except:
    from scipy.integrate import cumtrapz

class sim_data:
    """
    Class that accesses the fusion-io functions in a more pythonic way.
    Invoke an element of the class as follows:
    sim = sim_data(filename='../C1.h5')
    This creates a simulation object. This object has four attributes:
    - Fields
    - Time traces
    - Mesh
    - Diagnostics
    - Constants

    --Fields--
    Fields objects can be created by, for example,
    magnetic_field = sim.get_field('B',time=10)
    This field can in turn be evaluated using
    magnetic_field.evaluate((R,phi,Z)).
    !! Fields are given in SI units !!
    **Technical info**
    A list of all strings allowed to call on fields can be
    found in sim.available_fields. This dictionary is made so that
    one string contains all information on the type of field
    (scalar/vector), and the species. The fields make use of
    the sim_data library to interpolate fields in a consistent
    way.

    -- Time traces--
    Time trace objects can be created by, for example,
    dt = sim.get_time_trace('dt').
    The object dt then has two attributes:
    dt.time, this is an array of all the times where dt is evaluated
    dt.values, this is an array of all values of dt.
    To plot some arbitrary time trace plt.plot(dt.time,dt.values) is enough.
    A list of all callable traces can be found by invoking
    sim.available_traces.
    !! Time series are given in M3D-C1 units !!
    **Technical info**
    The time traces are read directly from the C1.h5 file. This is different
    from how the fields are evaluated (which use the sim_data library).

    --Mesh--
    Mesh object can be created by, for example
    mesh = sim.get_mesh(time=1).
    This object has a few attributes:
    mesh.nplanes shows how many toroidal planes were used
    mesh.elements contains info on the actual locations. The specifics
    of this are explained in S. Jardin's 2004 paper. But useful info is:
    mesh.elements[:,4] are the R locations of the mesh-points
    mesh.elements[:,5] are the Z locations of the mesh-points
    **Technical info**
    The mesh is read directly from the time_xxx.h5 files.

    --Diagnostics--
    The diagnostic object can be created by, for example
    diagnostic = sim.get_diagnostics(diagnostic)
    The diagnostic object is meant for non-physical parameters such as
    number of iterations, timings for each solve, and time_slice - time correspondensies.
    The diagnostic object has two attributes:
    !! Diagnostic are given in M3D-C1 units !!
    diagnostic.diagnostic - contains the value of the diagnostic. Can be a simple array (i.e. for
                            'time slices', or several objects i.e. to see a breakdown for timings)
    diagnostic.x_axis     - contains the values of the x-axis of the diagnostic. This is alfven times
                            for 'time slices' and iteration number for 'timings'
    To see all available diagnostics, check sim.available_diagnostics
    **Technical info**
    The diagnostic object is a bit of a hodge-podge of different methods. They all call on the
    C1.h5 file, and some make use of string manipulations.

    --Constants--
    On subclass containing the constants can be created by,
    for example,
    constants = sim.get_constants()
    In this subclass, stuff like gamma, R0, and the version
    number of the M3DC1 simulation are stored.
    """

    def __init__(self, filename='C1.h5', filetype='m3dc1', verbose=False, time=0, fast=False):
        """
        Initializes the fusion-io bindings to a file.

        Keyworded arguments:

        **filename**
            the name of the file which is to be read. Can include path.

        **filetype**
            the type of file.

        **verbose**
            if true, gives extra info on file read. the fio source does some printing anyways.

        **time**
            sets the time variable during initialization. Not necessarily needed as is also set
            when reading a field

        **fast**
            Workaround to keep mesh connectivity in memory. If True, magnetic field will be loaded
            upon initialization and the mesh connectivity will not be recalculated until the object
            is destroyed.
        """
        if filetype == 'm3dc1':
            ifiletype = fio_py.FIO_M3DC1_SOURCE
        else:
            print('Sorry, cannot do that Dave. Only supports M3DC1 for now. \n'
                'Please feel free to add to my functionality and add to the sim_data library!')

        # Dictionary contains the field abbreviation, the type of field, and the species
        self.typedict = {'j' : ('current density',  'vector',  None , 'simple'),
                         'ni': ('density',          'scalar', 'main ion' , 'simple'),
                         'ne': ('density',          'scalar', 'electron' , 'simple'),
                         'v' : ('fluid velocity',   'vector',  None , 'simple'),
                         'B' : ('magnetic field',   'vector',  None , 'simple'),
                         'p' : ('total pressure',   'scalar',  None , 'simple'),
                         'pi': ('pressure',         'scalar', 'main ion' , 'simple'),
                         'pe': ('pressure',         'scalar', 'electron' , 'simple'),
                         'alpha': ('alpha',         'scalar', None , 'simple'),
                         'ti': ('temperature',      'scalar', 'main ion' , 'simple'),
                         'te': ('temperature',      'scalar', 'electron' , 'simple'),
                         'A' : ('vector potential', 'vector',  None , 'simple'),
                         'gradA' : ('grad vector potential', 'tensor',  None , 'simple'),
                         'E' : ('electric field',   'vector',  None , 'simple'),
                         'psi' : ('psi',   'scalar',  None , 'composite'),
                         'kprad_rad' : ('total radiation',   'scalar',  None , 'simple')
                         }
        self.available_fields = self.typedict
        self.filename = os.path.abspath(filename)
        self._all_attrs       = h5py.File(self.filename, 'r')
        self._all_attrs_list  = list(self._all_attrs.keys())
        self._all_traces = self._all_attrs['scalars']
        self.available_traces = list(self._all_traces.keys())
        self.available_traces.extend(['bharmonics','keharmonics'])
        self.available_diagnostics = ['slice times','timings','iterations']
        self.fields = []
        self.seriess = []
        self._isrc  = fio_py.open_source(ifiletype, self.filename)
        self.hint   = fio_py.allocate_hint(self._isrc)
        fio_py.get_options(self._isrc)
        self.ntime = self._all_attrs.attrs["ntime"]
        self.set_timeslice(time)
        if fast:
            self._imag = fio_py.get_field(self._isrc, fio_py.FIO_MAGNETIC_FIELD)
            self.fields.append(self._imag)
        self._cs = fio_py.get_int_parameter(self._isrc, fio_py.FIO_GEOMETRY)
        self._iavailable_fields = fio_py.get_available_fields(self._isrc)
        #available fields is a dictionary from names to assigned integers
        self._available_fields = dict(zip([fio_py.get_field_name(nr) for nr in self._iavailable_fields], self._iavailable_fields))
        self.ntor = fio_py.get_int_parameter(self._isrc, fio_py.FIO_TOROIDAL_MODE)
        self.period = fio_py.get_real_parameter(self._isrc, fio_py.FIO_PERIOD)
        self.fc = None #used to store flux coodinate object upon calculation of flux coordinates (see class definition below and flux_coordinates.py)

        if verbose:
            print('Available fields:')
            if self._cs == fio_py.FIO_CYLINDRICAL:
                print('Using CYLINDRICAL coordinate system')
            else :
                print('Using CARTESIAN coordinate system')
            print('Number of time slices: ', self.ntime)
            print('Toroidal period = ', self.period)

    def __del__(self):
        print('deleting simulation object and closing {} fields'.format(
            len(self.fields)+len(self.seriess)))
        for ifield in self.fields:
            fio_py.close_field(ifield)
        for iseries in self.seriess:
            fio_py.close_series(iseries)
        fio_py.close_source(self._isrc)

    def set_timeslice(self,time):
        if isinstance(time,str):
            if time == 'last':
                time = self.ntime - 1
                print('last time slice = '+str(time))
        time = int(time)
        self.timeslice = time
        fio_py.set_int_option(fio_py.FIO_TIMESLICE, time)


    def get_time_trace(self,scalar,ipellet=None):
        """
        Makes an object containing a time-array, and an
        array with the correspond values of the physical
        quantity. Call on them using
        trace.time
        trace.values
        """
        return self.time_trace(scalar,sim_data=self,ipellet=ipellet)

    def get_diagnostic(self, diagnostic):
        """
        returns a diagnostics object.

        Available diagnostics:
        *slice times*
            relation between slice numbers and time in M3DC1 Alfven units
        *timings*
            Clock time needed per iteration of the sim
        *iterations*
            Number of iterations needed per time advance
        """
        return self.diagnostic(self,diagnostic)

    def get_field(self, field, time=None):
        """
        Returns a field object.

        Keyworded arguments:

        **field**
            Contains field name. Allowed fieldnames mysim.typedict

        **time**
            Timeslice for field to be read out
            If none, it will read for self.timeslice
        """

        if field in self.typedict.keys():
            return self.field(self, field=self.typedict[field][0], time=time, species=self.typedict[field][2], ftype=self.typedict[field][1])
        else:
            return self.field(self, field=field, time=time, species=None, ftype='scalar')

    def get_mesh(self, time=None, quiet=False):
        """
        Returns a mesh object.

        Keyworded arguments:

        **time**
            Timeslice for field to be read out
            If none, it will read for self.timeslice
        """
        return self.mesh(self, time, quiet=quiet)

    def get_signal(self,signame):
        """
        Return a signal object
        """
        return self.signal(self,signame)

    def get_constants(self):
        return self.constants(self)

    class field:
        """
        Field class: the init sets up the data access, and its bound methods
        return to you what you want. ex:
        myfield =mysim.get_field('...')
        Allowed strings are:
        'j'   - current density
        'ni'  - ion density
        'ne'  - electron density
        'v'   - fluid velocity
        'B'   - magnetic field
        'p'   - total pressure
        'pi'  - ion pressure
        'pe'  - electron pressure
        'ti'  - ion temperature
        'te'  - electron temperature
        'A'   - vector potential
        'E'   - electric field
        'kprad_rad' - total radiation
        Fields are evaluated using:
        myfield.evaluate((r,phi,theta))
        """
        def __init__(self, sim_data, field, time=None, species=None, ftype='scalar'):
            self.sim_data = sim_data
            self.field = field
            if time is not None:
                self.sim_data.set_timeslice(time)
            if species == 'electron':
                fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_ELECTRON)
            elif species == 'main ion':
                fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_MAIN_ION)
            if field in sim_data._available_fields.keys():
                itype  = sim_data._available_fields[field]
            else:
                itype = fio_py.FIO_SCALAR_FIELD
                fio_py.set_str_option(fio_py.FIO_FIELD_NAME, field)
            self._ifield = fio_py.get_field(sim_data._isrc, itype)
            self.sim_data.fields.append(self._ifield)
            self.ftype = ftype

        def __del__(self):
            self.sim_data.fields.remove(self._ifield)
            fio_py.close_field(self._ifield)

        def evaluate(self, x):
            """
            Evaluates the required field, returns a one-tuple (scalar quantity required)
            or a three-tuple (vector) or a nine-tuple (tensor) evaluated at the location x

            Arguments:

            **x**
                Three-tuple with the R, phi, z coordinate where the field is to be evaluated.
            """
            if self.ftype == 'vector':
                try:
                    return fio_py.eval_vector_field(self._ifield, x, self.sim_data.hint)
                except:
                    return (None,None,None)
            if self.ftype == 'tensor':
                try:
                    return fio_py.eval_tensor_field(self._ifield, x, self.sim_data.hint)
                except:
                    return (None,None,None,None,None,None,None,None,None)
            elif self.ftype == 'scalar':
                try:
                    return (fio_py.eval_scalar_field(self._ifield, x, self.sim_data.hint),)
                except:
                    return (None,)
            else:
                print('ftype not recognized!')
        
        def evaluate_deriv(self, x):
            """
            Evaluates the derivative of the required field, returns a three-tuple (scalar quantity required)
            or a nine-tuple (vector) evaluated at the location x

            Arguments:

            **x**
                Three-tuple with the R, phi, z coordinate where the field is to be evaluated.
            """
            if self.ftype == 'vector':
                try:
                    return fio_py.eval_vector_field_deriv(self._ifield, x, self.sim_data.hint)
                except:
                    return (None,None,None,None,None,None,None,None,None)
            elif self.ftype == 'scalar':
                try:
                    return (fio_py.eval_scalar_field_deriv(self._ifield, x, self.sim_data.hint),)
                except:
                    return (None,None,None)
            else:
                print('ftype not recognized!')

    class mesh:
        """
        Mesh class: init routine calls read_mesh() that reads the mesh from time slice file
        and stores the elements in an array. It also reads the output version and nplanes.
        The latter variable is needed for plotting 3D meshes.
        """
        def __init__(self, sim_data, time=None, quiet=False):
            self.sim_data = sim_data
            self._quiet = quiet
            if time is not None:
                self.sim_data.set_timeslice(time)
            self.elements, self.version, self.nplanes = self.read_mesh()

        def __eq__(self, other):
            return (np.array_equal(self.elements,other.elements) and
                    (self.version==other.version) and
                    (self.nplanes==other.nplanes))

        def read_mesh(self):
            if self.sim_data.timeslice == -1:
                group = "equilibrium"
            else:
                group = "time_%03d"%self.sim_data.timeslice

            f = self.sim_data._all_attrs
            mesh  = np.asarray(f[group+'/mesh/elements'])

            if "version" in f[group].attrs:
                version = f[group].attrs["version"]
            else:
                version = self.sim_data.get_constants().version

            nplanes= f[group+'/mesh'].attrs["nplanes"]

            if not self._quiet:
                print('Output version: '+str(version))
                print('Mesh shape: '+str(mesh.shape))

            return mesh, version, nplanes

    class signal:
        """
        Signal class: for diagnostic signals from magnetic probes and flux loops
        """
        def __init__(self, sim_data, signame):
            self.sim_data = sim_data
            self._all_attrs  = self.sim_data._all_attrs
            self._all_attrs_list = self.sim_data._all_attrs_list
            self.sigvalues = np.asarray(self.sim_data._all_attrs[signame+'/value'][:])
            

    class time_trace:
        """
        time_trace class: init routine calls read_scalars() that reads all possible scalars
        from the C1.h5 file.
        """
        def __init__(self, scalar, time=None, sim_data=None, ipellet=None):
            if (time is None) and (sim_data is None):
                raise ValueError('time or sim_data must not be None')

            if time is None:
                time = np.asarray(sim_data._all_traces['time'])

            if isinstance(scalar,str):
                if scalar in ['bharmonics','keharmonics']:
                    values = np.asarray(sim_data._all_attrs['%s/%s'%(scalar,scalar)])
                elif 'pellet/%s'%scalar in sim_data._all_attrs:
                    values = np.asarray(sim_data._all_attrs['pellet/%s'%scalar])
                    if ipellet is not None:
                        values = values[:,ipellet]
                else:
                    values = np.asarray(sim_data._all_traces[scalar])
            else:
                values = np.asarray(scalar)

            if len(time) != len(values):
                if abs(len(time) - len(values))==1:
                    print('WARNING: len(time) and len(values) differ by 1. Simulation might currently be running.')
                    if len(time) > len(values):
                        self.time = time[:len(values)]
                        self.values = values
                else:
                    #raise ValueError('time_trace time and values have different lengths. len(time) = '+str(len(time)) + ', len(values) = ' + str(len(values)))
                    print('WARNING: time_trace time and values have different lengths. len(time) = '+str(len(time)) + ', len(values) = ' + str(len(values)))
                    if len(time) > len(values):
                        self.time = time[:len(values)]
                        self.values = values
            else:
                if values.ndim == 1:
                    not_nan = ~np.isnan(values)
                    self.time = time[not_nan]
                    self.values = values[not_nan]
                else:
                    self.time = time
                    self.values = values

        # Addition
        def __add__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                values = self.values + other.values
            else:
                values = self.values + other
            return sim_data.time_trace(values,time=self.time)

        __radd__ = __add__

        def __iadd__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                self.values = self.values + other.values
            else:
                self.values = self.values + other
            return self

        # Subtraction
        def __sub__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                values = self.values - other.values
            else:
                values = self.values - other
            return sim_data.time_trace(values,time=self.time)

        def __rsub__(self,other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                values = other.values - self.values
            else:
                values = other - self.values
            return sim_data.time_trace(values,time=self.time)

        def __isub__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                self.values = self.values - other.values
            else:
                self.values = self.values - other
            return self

        # Multiplication
        def __mul__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                values = self.values*other.values
            else:
                values = self.values*other
            return sim_data.time_trace(values,time=self.time)

        __rmul__ = __mul__

        def __imul__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                self.values = self.values*other.values
            else:
                self.values = self.values*other
            return self

        # Division
        def __truediv__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                values = self.values/other.values
            else:
                values = self.values/other
            return sim_data.time_trace(values,time=self.time)

        def __rtruediv__(self,other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                values = other.values/self.values
            else:
                values = other/self.values
            return sim_data.time_trace(values,time=self.time)

        def __itruediv__(self, other):
            if isinstance(other,type(self)):
                if not np.array_equal(self.time,other.time):
                    raise ValueError('operands share different .time attributes')
                self.values = self.values/other.values
            else:
                self.values = self.values/other
            return self

        # Exponetiation
        def __pow__(self, other):
            return sim_data.time_trace(self.values**other,time=self.time)

        #__rpow__ is intentionally NotImplemented

        def __ipow__(self, other):
            self.values = self.values**other
            return self

        def __abs__(self):
            return sim_data.time_trace(abs(self.values),time=self.time)
        def __pos__(self):
            return sim_data.time_trace(+self.values,time=self.time)
        def __neg__(self):
            return sim_data.time_trace(-self.values,time=self.time)
        
        
        def cum_int(self,nts):
            values = cumtrapz(self.values[:nts],self.time[:nts],initial=0.0)
            return sim_data.time_trace(values,time=self.time)

    class diagnostic:
        def __init__(self, sim_data, diagnostic):
            self.sim_data = sim_data
            # Multiple cases depending on the diagnostic

            # Case corresponding to knowing the alfven times of the time slices
            if diagnostic == 'slice times':
                self.diagnostic = np.arange(sim_data.ntime)
                self.x_axis = np.zeros(sim_data.ntime)
                for t in self.diagnostic:
                    # Retrieve the time attribute
                    self.x_axis[t] = sim_data._all_attrs['time_%03d'%t].attrs['time']

            # Case corresponding to the timings of each iteration
            # self.diagnostic holds several objects, so you can see a breakdown of the
            # timings as well.
            # x_axis contains the total time step length
            if diagnostic == 'timings':
                self.diagnostic = {}
                for timing in sim_data._all_attrs['timings']:
                    self.diagnostic[timing] = np.array(sim_data._all_attrs['timings'][timing])
                self.x_axis = np.arange(self.diagnostic['t_onestep'].shape[0])+1

            # Case corresponding to the iterations of each step
            # self.diagnostic is a set of arrays for different iterations
            # x_axis contains iteration number
            if diagnostic =='iterations':
                self.diagnostic = np.asarray(sim_data._all_attrs['kspits/kspits'])
                self.x_axis     = np.arange(np.shape(self.diagnostic[:,0])[0])+1


    class constants:
        def __init__(self, sim_data):
            self.sim_data = sim_data
            self.R0        = sim_data._all_attrs.attrs['rzero']
            self.B0        = sim_data._all_attrs.attrs['bzero']
            self.gamma     = sim_data._all_attrs.attrs['gam']
            self.amupar    = sim_data._all_attrs.attrs['amupar']
            self.version   = sim_data._all_attrs.attrs['version']
            self.numvar    = sim_data._all_attrs.attrs['numvar']
            self.itor      = sim_data._all_attrs.attrs['itor']
            self.is3D      = sim_data._all_attrs.attrs['3d']


class flux_coordinates:
    """
    Class that represents a flux surface coordinate system, e.g. PEST, Boozer or Hamada coordinates.
    """
    def __init__(self, m,n,rpath,zpath,axis,omega,psi,psin,period,theta,jac,q,area,polarea,dV,fcoords,V,phi,itor,r0,current,dpsi_dpsin,points):
        """
        Initializes the flux_coordinates.
        """
        self.m = m
        self.n = n
        self.rpath = rpath
        self.zpath = zpath
        self.rma = axis[0]
        self.zma = axis[1]
        self.omega = omega
        self.psi = psi
        self.psi_norm = np.asarray(psin)
        self.flux_pol = -period*(psi-psi[0])
        self.theta = theta
        self.j = jac
        self.q = q
        self.area = area
        self.polarea = polarea
        self.dV_dchi = dV
        self.fcoords = fcoords
        #self.pest = 
        #self.boozer = 
        #self.hamada = 
        self.V = V
        self.flux_tor = phi
        self.phi_norm = phi/phi[n-1]
        self.rho = np.sqrt(phi/phi[n-1])
        self.period = period
        self.itor = itor
        self.r0 = r0
        self.current = current
        self.dpsi_dchi = dpsi_dpsin
        self.points = points
