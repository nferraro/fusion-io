# fpy.py: wrapper functions for the sim_data library to make accessing
# simulation output more pythonic.
#
#
# coded by Christopher Berg Smiet on 18 January 2019
# Edited by Ralf Mackenbach on the 14th of August 2019
# Mesh added by Andreas Kleiner on 20 August 2019
# csmiet@pppl.gov
# rmackenb@pppl.gov
# akleiner@pppl.gov

import fio_py
import numpy as np
import h5py

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
    dt = sim.get_time_traces('dt').
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

    def __init__(self, filename='C1.h5', filetype='m3dc1', verbose=False, time = 0):
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
        """
        if filetype == 'm3dc1':
            ifiletype = fio_py.FIO_M3DC1_SOURCE
        else:
            print('Sorry, cannot do that Dave. Only supports M3DC1 for now. \n'
                'Please feel free to add to my functionality and add to the sim_data library!')

        # Dictionary contains the field abbreviation, the type of field, and the species
        self.typedict = {'j' : ('current density',  'vector',  None ),
                         'ni': ('density',          'scalar', 'main ion'),
                         'ne': ('density',          'scalar', 'electron'),
                         'v' : ('fluid velocity',   'vector',  None ),
                         'B' : ('magnetic field',   'vector',  None ),
                         'p' : ('total pressure',   'scalar',  None ),
                         'pi': ('pressure',         'scalar', 'main ion'),
                         'pe': ('pressure',         'scalar', 'electron'),
                         'alpha': ('alpha',         'scalar', None),
                         'ti': ('temperature',      'scalar', 'main ion'),
                         'te': ('temperature',      'scalar', 'electron'),
                         'A' : ('vector potential', 'vector',  None ),
                         'gradA' : ('grad vector potential', 'tensor',  None ),
                         'E' : ('electric field',   'vector',  None )}
        self.available_fields = self.typedict
        self._all_attrs       = h5py.File(filename, 'r')
        self._all_attrs_list  = list(h5py.File(filename, 'r').keys())
        self._all_traces = h5py.File(filename, 'r')['scalars']
        self.available_traces = list(self._all_traces.keys())
        self.available_traces.extend(['bharmonics','keharmonics'])
        self.available_diagnostics = ['slice times','timings','iterations']
        self.fields = []
        self.seriess = []
        self._isrc  = fio_py.open_source(ifiletype, filename)
        self.hint   = fio_py.allocate_hint(self._isrc)
        fio_py.get_options(self._isrc)
        #self.ntime = fio_py.get_int_parameter(self._isrc, fio_py.FIO_NUM_TIMESLICES)
        self.ntime = self._all_attrs.attrs["ntime"]
        if time == 'last':
            time = self.ntime - 1
            print('last time slice = '+str(time))
        fio_py.set_int_option(fio_py.FIO_TIMESLICE, int(time))
        self._imag = fio_py.get_field(self._isrc, fio_py.FIO_MAGNETIC_FIELD)
        self.fields.append(self._imag)
        self._cs = fio_py.get_int_parameter(self._isrc, fio_py.FIO_GEOMETRY)
        self._iavailable_fields = fio_py.get_available_fields(self._isrc)
        #available fields is a dictionary from names to assigned integers
        self._available_fields = dict(zip([fio_py.get_field_name(nr) for nr in self._iavailable_fields], self._iavailable_fields))
        self.ntor = fio_py.get_int_parameter(self._isrc, fio_py.FIO_TOROIDAL_MODE)
        self.period = fio_py.get_real_parameter(self._isrc, fio_py.FIO_PERIOD)
        self.fc = None #used to store flux coodinate object upon calculation of flux coordinates (see class definition below and flux_coordinates.py)
        #if time == -1:
        #    self.timeslice = int(0)
        #else:
        self.timeslice = int(time)
        self.time = fio_py.get_real_field_parameter(self._imag, fio_py.FIO_TIME)
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

    def get_time_traces(self,scalar):
        """
        Makes an object containing a time-array, and an
        array with the correspond values of the physical
        quantity. Call on them using
        trace.time
        trace.values
        """
        return self.time_traces(self,scalar)

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

    def get_field(self, field, time):
        """
        Returns a field object.

        Keyworded arguments:

        **field**
            Contains field name. Allowed fieldnames mysim.typedict

        **time**
            Timeslice for field to be read out
        """

        return self.field(self, field=self.typedict[field][0], time=time, species=self.typedict[field][2], ftype=self.typedict[field][1])

    def get_mesh(self, time=0):
        """
        Returns a mesh object.

        Keyworded arguments:

        **time**
            Timeslice for field to be read out
        """
        return self.mesh(self, time)

    def get_signal(self,filename,signame):
        """
        Return a signal object
        """
        return self.signal(self,filename,signame)

    def get_constants(self):
        return self.constants(self)

    class field:
        """
        Field class: the init sets up the data access, and its bound methods
        return to you what you want. ex:
        myfield =mysim.get_field('...')
        Allowed strings are:
        'j'  - current density
        'ni' - ion density
        'ne' - electron density
        'v'  - fluid velocity
        'B'  - magnetic field
        'p'  - total pressure
        'pi' - ion pressure
        'pe' - electron pressure
        'ti' - ion temperature
        'te' - electron temperature
        'A'  - vector potential
        'E'  - electric field
        Fields are evaluated using:
        myfield.evaluate((r,phi,theta))
        """
        def __init__(self, sim_data, field, time, species, ftype):
            self.sim_data = sim_data
            self.field = field
            fio_py.set_int_option(fio_py.FIO_TIMESLICE, time)
            if   species == 'electron':
                fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_ELECTRON)
            elif species == 'main ion':
                fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_MAIN_ION)
            itype  = sim_data._available_fields[field]
            self._ifield = fio_py.get_field(sim_data._isrc, itype)
            self.sim_data.fields.append(self._ifield)
            self.time  = fio_py.get_real_field_parameter(self._ifield, fio_py.FIO_TIME)
            self.ftype = ftype

        def __del__(self):
            fio_py.close_field(self._ifield)

        def evaluate(self, x):
            """
            Evaluates the required field, returns a one-tuple (scalar quantity required)
            or a three-tuple (vector) evaluated at the location x

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

    class mesh:
        """
        Mesh class: init routine calls read_mesh() that reads the mesh from time slice file
        and stores the elements in an array. It also reads the output version and nplanes.
        The latter variable is needed for plotting 3D meshes.
        """
        def __init__(self, sim_data, time=0):
            self.sim_data = sim_data
            fio_py.set_int_option(fio_py.FIO_TIMESLICE, int(time))
            itype  = sim_data._available_fields['magnetic field']
            #'magnetic field' is just used as a dummy field to read time
            self.time = fio_py.get_real_field_parameter(fio_py.get_field(sim_data._isrc, itype), fio_py.FIO_TIME)
            self.elements, self.version, self.nplanes = self.read_mesh(sim_data, time)

        def read_mesh(self, sim_data, time):
            self.sim_data = sim_data
            timestr = str(time)
            if timestr == '-1':
                fname = 'equilibrium.h5'
            else:
                fname = "time_"+timestr.zfill(3)+'.h5'
            f     = h5py.File(fname, 'r')
            mesh  = np.asarray(f['mesh/elements'])

            version = f.attrs["version"]
            print('Output version: '+str(version))

            meshshape = mesh.shape
            print('Mesh shape: '+str(meshshape))

            dset = f["mesh"]
            nplanes = dset.attrs["nplanes"]

            return mesh, version, nplanes

    class signal:
        """
        Signal class: for diagnostic signals from magnetic probes and flux loops
        """
        def __init__(self, sim_data, filename, signame):
            self.sim_data = sim_data
            self._all_attrs       = h5py.File(filename, 'r')
            self._all_attrs_list  = list(h5py.File(filename, 'r').keys())
            self.sigvalues = h5py.File(filename, 'r')[signame+'/value']
            

    class time_traces:
        """
        time_trace class: init routine calls read_scalars() that reads all possible scalars
        from the C1.h5 file.
        """
        def __init__(self, sim_data, scalar):
            self.sim_data = sim_data
            self.time      = self.sim_data._all_traces['time'][()]
            if scalar is not 'bharmonics' and scalar is not 'keharmonics':
                self.values    = self.sim_data._all_traces[scalar][()]
            if scalar is 'bharmonics':
                self.values    = self.sim_data._all_attrs['bharmonics/bharmonics'][()]
            if scalar is 'keharmonics':
                self.values    = self.sim_data._all_attrs['keharmonics/keharmonics'][()]

    class diagnostic:
        def __init__(self, sim_data, diagnostic):
            self.sim_data = sim_data
            # Multiple cases depending on the diagnostic

            # Case corresponding to knowing the alfven times of the time slices
            if diagnostic == 'slice times':
                times = []
                for string in sim_data._all_attrs_list:
                    if string.startswith('time_') == True:
                        times.append(string)
                        times.sort()
                alfven_times = []
                for time in times:
                    # Retrieve the time attribute
                    alfven_times.append(sim_data._all_attrs[time].attrs['time'])

                for (idx,time) in enumerate(times):
                    times[idx] = int(time.replace('time_',''))

                self.diagnostic = times
                self.x_axis     = alfven_times

            # Case corresponding to the timings of each iteration
            # self.diagnostic holds several objects, so you can see a breakdown of the
            # timings as well.
            # x_axis contains the iteration number
            if diagnostic == 'timings':
                self.diagnostic = sim_data._all_attrs['timings']
                self.x_axis     = np.arange(self.diagnostic['t_onestep'].shape[0])+1

            # Case corresponding to the iterations of each step
            # self.diagnostic is a set of arrays for different iterations
            # x_axis contains iteration number
            if diagnostic =='iterations':
                self.diagnostic = sim_data._all_attrs['kspits/kspits'][()]
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
    def __init__(self, m,n,rpath,zpath,axis,omega,psi,psin,period,theta,jac,q,area,dV,fcoords,V,phi,itor,r0,current,dpsi_dpsin):
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
        self.psi_norm = psin
        self.flux_pol = -period*(psi-psi[0])
        self.theta = theta
        self.j = jac
        self.q = q
        self.area = area
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
