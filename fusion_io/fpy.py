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
import os
import io
import re
import hashlib
import shlex
import sqlite3
import subprocess
import sys
import zlib

try:
    import h5py
    _have_h5py = True
except ImportError:
    _have_h5py = False

try:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except ImportError:
    from scipy.integrate import cumtrapz


# ---------------------------------------------------------------------------
# Remote scripts executed on the remote host via SSH
# ---------------------------------------------------------------------------

# Reads HDF5 root-level or group-level attributes and returns them as JSON.
_ATTRS_SCRIPT = r"""
import json, sys
import h5py
path, objpath = sys.argv[1], sys.argv[2]
def plain(v):
    try:
        return v.item() if getattr(v, "shape", None) == () else v.tolist()
    except AttributeError:
        return v.decode() if isinstance(v, bytes) else v
with h5py.File(path, "r") as f:
    obj = f[objpath] if objpath != "/" else f
    attrs = {k: plain(v) for k, v in obj.attrs.items()}
sys.stdout.write("###JSON###" + json.dumps(attrs))
"""

# Reads one HDF5 variable (optionally sliced) and writes it as a .npy stream.
_REMOTE_READ_SCRIPT = r"""
import io, sys
import h5py, numpy as np
path, var, slc = sys.argv[1], sys.argv[2], sys.argv[3]
with h5py.File(path, "r") as f:
    d = f[var]
    if slc != "-":
        start, _, stop = slc.partition(":")
        data = d[int(start or 0):(int(stop) if stop else None)]
    else:
        data = d[()]
buf = io.BytesIO()
np.save(buf, np.asarray(data), allow_pickle=False)
sys.stdout.buffer.write(buf.getvalue())
sys.stdout.buffer.flush()
"""

_NPY_MAGIC = b"\x93NUMPY"
_TIME_RE = re.compile(r"time[_.-]?(\d+)\.h5$", re.IGNORECASE)


class _AcaBackend:
    """
    Backend for reading M3D-C1 data from a remote cluster via an HPC campaign
    archive (.aca file), without using ADIOS2.

    A .aca file is a SQLite database that indexes HDF5 datasets living on one
    or more remote hosts. For every HDF5 dataset it also embeds a
    zlib-compressed metadata skeleton: the full HDF5 group/dataset tree, all
    attributes, and the actual values of small datasets (<= 128 elements), but
    no large array data.

    Strategy
    --------
    - Metadata and small datasets: read from the embedded HDF5 skeleton
      (local, no network access).
    - File-root attributes: fetched from the remote host via SSH + h5py on
      first access, then cached in memory.
    - Large arrays (time traces, mesh data, field coefficients): streamed
      from the remote host via SSH + h5py and cached locally as .npy files.
    - Field evaluation: field coefficient arrays and mesh arrays are streamed
      per-field per-timeslice and injected into a fio_py numpy source via
      fio_py.set_numpy_field() and fio_py.set_numpy_mesh(). The C++ fusion-io
      code then evaluates the finite element basis functions exactly as it
      would for a local file, with no full HDF5 file transfer required.
    """

    def __init__(self, aca_filename, campaign_simulation='', login='',
                 cache_dir='.cache/aca', remote_python='python3',
                 remote_setup='module load python >/dev/null 2>&1 || true',
                 verbose=False):
        self._login = login
        self._cache_dir = cache_dir
        self._remote_python = remote_python
        self._remote_setup = remote_setup
        self._verbose = verbose

        if not os.path.isfile(aca_filename):
            raise FileNotFoundError(
                'Campaign archive not found: {}'.format(aca_filename))

        self._con = sqlite3.connect(
            'file:{}?mode=ro'.format(aca_filename), uri=True)
        cur = self._con.cursor()
        row = cur.execute(
            'SELECT id, name, version FROM info').fetchone()
        self._archive_id, self._archive_name, self._version = row

        tables = {r[0] for r in cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        self._has_repfiles = 'repfiles' in tables
        if 'replica' not in tables:
            raise RuntimeError(
                'Unsupported archive schema (version {}). '
                'Run: hpc_campaign manager <aca> upgrade'.format(
                    self._version))

        self._simulation = campaign_simulation.strip('/')
        self._datasets = self._load_datasets()

        prefix = self._simulation + '/' if self._simulation else ''
        self._files = {}   # basename -> dataset dict
        for ds in self._datasets:
            if ds['fileformat'] == 'HDF5' and ds['name'].startswith(prefix):
                self._files[os.path.basename(ds['name'])] = ds
        if not self._files:
            available = sorted(
                {os.path.dirname(d['name']) for d in self._datasets})
            raise KeyError(
                "No HDF5 files under '{}' in {}. "
                "Registered simulations: {}".format(
                    self._simulation, aca_filename, available))

        self._timeslices = sorted(
            (int(_TIME_RE.search(n).group(1)), n)
            for n in self._files if _TIME_RE.search(n))

        if 'C1.h5' in self._files:
            self._default_file = 'C1.h5'
        elif self._timeslices:
            self._default_file = self._timeslices[-1][1]
        else:
            self._default_file = next(iter(self._files))

        self._skeletons = {}          # ds['name'] -> local skeleton path
        self._root_attrs_cache = {}   # file basename -> dict of attrs
        self._mesh_injected = False   # whether mesh has been injected yet

    # ------------------------------------------------------------------
    # SQLite helpers
    # ------------------------------------------------------------------

    def _load_datasets(self):
        cur = self._con.cursor()
        hosts = {
            r[0]: (r[1], r[2])
            for r in cur.execute(
                'SELECT rowid, hostname, longhostname FROM host')}
        dirs = {
            r[0]: (r[1], r[2])
            for r in cur.execute(
                'SELECT rowid, hostid, name FROM directory')}

        datasets = {}
        for did, uuid, name, fmt, deltime in cur.execute(
                'SELECT rowid, uuid, name, fileformat, deltime FROM dataset'):
            if deltime:
                continue
            datasets[did] = {
                'id': did, 'uuid': uuid, 'name': name,
                'fileformat': fmt, 'replicas': [], 'embedded_files': []}

        for rid, dsid, hostid, dirid, name, deltime, size in cur.execute(
                'SELECT rowid, datasetid, hostid, dirid, name, '
                'deltime, size FROM replica'):
            ds = datasets.get(dsid)
            if ds is None or deltime:
                continue
            hostname, longhost = hosts.get(hostid, ('?', '?'))
            _, dirname = dirs.get(dirid, (None, '?'))
            ds['replicas'].append({
                'id': rid, 'hostname': hostname, 'longhost': longhost,
                'directory': dirname, 'name': name, 'size': size or 0,
                'remote_path': os.path.join(dirname, name)})

        if self._has_repfiles:
            file_q = (
                'SELECT r.datasetid, f.fileid, f.name, f.lenorig, '
                'f.lencompressed FROM file f '
                'JOIN repfiles rf ON f.fileid = rf.fileid '
                'JOIN replica r ON rf.replicaid = r.rowid')
        else:
            file_q = (
                'SELECT r.datasetid, f.rowid, f.name, f.lenorig, '
                'f.lencompressed FROM file f '
                'JOIN replica r ON f.replicaid = r.rowid')
        for dsid, fid, name, lenorig, lencomp in cur.execute(file_q):
            ds = datasets.get(dsid)
            if ds is None:
                continue
            if all(e['id'] != fid for e in ds['embedded_files']):
                ds['embedded_files'].append({
                    'id': fid, 'name': name,
                    'len_orig': lenorig or 0,
                    'len_compressed': lencomp or 0})

        return list(datasets.values())

    def _get_dataset(self, basename):
        if basename not in self._files:
            raise KeyError(
                "'{}' not in this simulation; files: {}".format(
                    basename, list(self._files)))
        return self._files[basename]

    def _live_replica(self, ds):
        reps = [r for r in ds['replicas']]
        if not reps:
            raise RuntimeError(
                'No live replica for {}'.format(ds['name']))
        return reps[0]

    # ------------------------------------------------------------------
    # Embedded skeleton extraction
    # ------------------------------------------------------------------

    def _skeleton(self, basename=None):
        basename = basename or self._default_file
        ds = self._get_dataset(basename)
        if ds['name'] not in self._skeletons:
            cur = self._con.cursor()
            written = []
            for ef in ds['embedded_files']:
                fid = ef['id']
                if self._has_repfiles:
                    row = cur.execute(
                        'SELECT data, compression FROM file '
                        'WHERE fileid = ?', (fid,)).fetchone()
                else:
                    row = cur.execute(
                        'SELECT data, compression FROM file '
                        'WHERE rowid = ?', (fid,)).fetchone()
                if row is None or row[0] is None:
                    continue
                blob, compression = row
                data = zlib.decompress(blob) if compression else bytes(blob)
                dest = os.path.join(
                    self._cache_dir, 'skeletons', ds['name'],
                    os.path.basename(ef['name']))
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with open(dest, 'wb') as f:
                    f.write(data)
                written.append(dest)
            if not written:
                raise RuntimeError(
                    'No embedded metadata in .aca for {}'.format(ds['name']))
            self._skeletons[ds['name']] = written[0]
        return self._skeletons[ds['name']]

    # ------------------------------------------------------------------
    # SSH helpers
    # ------------------------------------------------------------------

    def _remote_run(self, script, argv, what):
        inner = (
            '{setup}; exec {python} - {args}'.format(
                setup=self._remote_setup,
                python=self._remote_python,
                args=' '.join(shlex.quote(a) for a in argv)))
        if self._login:
            cmd = ['ssh', self._login, 'bash', '-lc', shlex.quote(inner)]
        else:
            cmd = ['bash', '-c', inner]
        proc = subprocess.run(
            cmd, input=script.encode(), capture_output=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                'Remote {} failed:\n{}'.format(
                    what, proc.stderr.decode(errors='replace')[-2000:]))
        return proc.stdout

    def _var_cache_path(self, ds, variable, slc):
        rep_mt = max((r['id'] for r in ds['replicas']), default=0)
        key = hashlib.sha1(
            '{}:{}:{}:{}'.format(
                ds['uuid'], rep_mt, variable, slc or '').encode()
        ).hexdigest()[:16]
        return os.path.join(self._cache_dir, 'vars', '{}.npy'.format(key))

    def _log(self, msg):
        if self._verbose:
            print('[fpy._AcaBackend] {}'.format(msg), file=sys.stderr)

    # ------------------------------------------------------------------
    # Array reading
    # ------------------------------------------------------------------

    def read(self, variable, basename=None, slc=None):
        """
        Return one HDF5 variable as a numpy array.

        Resolution order:
        1. Embedded value in the skeleton (small datasets, full reads only).
        2. Per-variable disk cache from a previous streamed read.
        3. Stream from the remote host via SSH + h5py.
        """
        ds = self._get_dataset(basename or self._default_file)
        variable = variable.lstrip('/')

        # 1. Embedded skeleton (small datasets only)
        if slc is None and _have_h5py:
            try:
                with h5py.File(self._skeleton(basename), 'r') as f:
                    if variable in f and isinstance(f[variable], h5py.Dataset):
                        obj = f[variable]
                        if obj.size <= 128:
                            self._log(
                                '{}: from embedded skeleton'.format(variable))
                            return obj[()]
            except Exception:
                pass

        # 2. Per-variable disk cache
        cache_file = self._var_cache_path(ds, variable, slc)
        if os.path.isfile(cache_file):
            self._log(
                '{}: from variable cache {}'.format(variable, cache_file))
            return np.load(cache_file, allow_pickle=False)

        # 3. Stream from remote
        rep = self._live_replica(ds)
        self._log('{}: streaming from {}:{}'.format(
            variable, rep['hostname'], rep['remote_path']))
        out = self._remote_run(
            _REMOTE_READ_SCRIPT,
            [rep['remote_path'], variable, slc or '-'],
            "read of '{}' from {}".format(variable, rep['remote_path']))
        pos = out.find(_NPY_MAGIC)
        if pos < 0:
            raise RuntimeError(
                "No array data returned for '{}'".format(variable))
        data = np.load(io.BytesIO(out[pos:]), allow_pickle=False)
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        np.save(cache_file, data, allow_pickle=False)
        return data

    # ------------------------------------------------------------------
    # Attribute reading
    # ------------------------------------------------------------------

    def _fetch_root_attrs(self, basename=None):
        import json
        key = basename or self._default_file
        if key not in self._root_attrs_cache:
            ds = self._get_dataset(key)
            rep = self._live_replica(ds)
            out = self._remote_run(
                _ATTRS_SCRIPT,
                [rep['remote_path'], '/'],
                'root attrs of {}'.format(ds['name']))
            pos = out.find(b'###JSON###')
            if pos < 0:
                self._root_attrs_cache[key] = {}
            else:
                self._root_attrs_cache[key] = json.loads(
                    out[pos + len(b'###JSON###'):].decode())
        return self._root_attrs_cache[key]

    def get_root_attr(self, name, basename=None):
        """Read a scalar attribute from the file root."""
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(basename), 'r') as f:
                    if name in f.attrs:
                        v = f.attrs[name]
                        return v.item() if hasattr(v, 'item') else v
            except Exception:
                pass
        return self._fetch_root_attrs(basename).get(name)

    def get_group_attr(self, group, name, basename=None):
        """Read a scalar attribute from a group inside the file."""
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(basename), 'r') as f:
                    if group in f and name in f[group].attrs:
                        v = f[group].attrs[name]
                        return v.item() if hasattr(v, 'item') else v
            except Exception:
                pass
        import json
        ds = self._get_dataset(basename or self._default_file)
        rep = self._live_replica(ds)
        out = self._remote_run(
            _ATTRS_SCRIPT,
            [rep['remote_path'], group],
            'attrs of {} in {}'.format(group, ds['name']))
        pos = out.find(b'###JSON###')
        if pos < 0:
            return None
        attrs = json.loads(out[pos + len(b'###JSON###'):].decode())
        return attrs.get(name)

    def list_group(self, group, basename=None):
        """Return the names of direct children of a group in the skeleton."""
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(basename), 'r') as f:
                    if group in f:
                        return list(f[group].keys())
            except Exception:
                pass
        return []

    def group_exists(self, group, basename=None):
        """Return True if the group exists in the skeleton."""
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(basename), 'r') as f:
                    return group in f
            except Exception:
                pass
        return False

    # ------------------------------------------------------------------
    # Injection helpers
    # ------------------------------------------------------------------

    def _ts_file_and_group(self, timeslice):
        """Return (basename, hdf5_group_name) for a given integer timeslice."""
        if timeslice < 0:
            if self._timeslices:
                timeslice = self._timeslices[timeslice][0]
            else:
                return self._default_file, 'equilibrium'
        group = 'time_{:03d}'.format(timeslice)
        if self._timeslices:
            for num, fname in self._timeslices:
                if num == timeslice:
                    return fname, group
            idx = timeslice % len(self._timeslices)
            return self._timeslices[idx][1], group
        return self._default_file, group

    def inject_parameters(self, isrc):
        """
        Read physical normalization parameters from the remote C1.h5 root
        attributes and inject them into the fio_py numpy source.
        """
        def _f(name, default=0.0):
            v = self.get_root_attr(name)
            return float(v) if v is not None else default

        def _i(name, default=0):
            v = self.get_root_attr(name)
            return int(v) if v is not None else default

        fio_py.set_numpy_parameters(
            isrc,
            _f('bzero'), _f('rzero'),
            _i('ntime'), _i('ntor'),
            _f('ion_mass', 1.0), _f('z_ion', 1.0),
            _i('itor', 1), _i('numvar'), _i('version'),
            _f('n0_norm'), _f('l0_norm'), _f('b0_norm'))

    def inject_mesh(self, isrc, timeslice):
        """
        Stream the mesh arrays from the appropriate time-slice HDF5 file on
        the remote host and inject them into the fio_py numpy source.

        The mesh is the same for all time slices in M3D-C1, so this only
        needs to be done once. Subsequent calls are no-ops.
        """
        if self._mesh_injected:
            return

        ts_file, group = self._ts_file_and_group(timeslice)
        self._log('Injecting mesh from {} / {}'.format(ts_file, group))

        elements = self.read(group + '/mesh/elements', basename=ts_file)

        version = int(self.get_root_attr('version') or 0)

        nelms = elements.shape[0]

        a  = np.ascontiguousarray(elements[:, 0], dtype=np.float64)
        b  = np.ascontiguousarray(elements[:, 1], dtype=np.float64)
        c  = np.ascontiguousarray(elements[:, 2], dtype=np.float64)
        # Column 3 is the angle; C++ stores cos/sin separately.
        co = np.ascontiguousarray(np.cos(elements[:, 3]), dtype=np.float64)
        sn = np.ascontiguousarray(np.sin(elements[:, 3]), dtype=np.float64)
        x  = np.ascontiguousarray(elements[:, 4], dtype=np.float64)
        z  = np.ascontiguousarray(elements[:, 5], dtype=np.float64)
        bound = np.ascontiguousarray(
            elements[:, 6].astype(np.int32), dtype=np.int32)
        if version >= 15:
            region = np.ascontiguousarray(
                elements[:, 7].astype(np.int32), dtype=np.int32)
        else:
            region = np.zeros(nelms, dtype=np.int32)

        nplanes_attr = self.get_group_attr(
            group + '/mesh', 'nplanes', basename=ts_file)
        nplanes = int(nplanes_attr) if nplanes_attr is not None else 1

        itor = int(self.get_root_attr('itor') or 1)
        if itor == 1:
            period = 2.0 * np.pi
        else:
            rzero  = float(self.get_root_attr('rzero') or 1.0)
            L0_cgs = float(self.get_root_attr('l0_norm') or 1.0)
            L0_mks = L0_cgs / 100.0
            period = 2.0 * np.pi * rzero * L0_mks

        fio_py.set_numpy_mesh(
            isrc, nelms,
            a.tobytes(), b.tobytes(), c.tobytes(),
            co.tobytes(), sn.tobytes(),
            x.tobytes(), z.tobytes(),
            bound.tobytes(), region.tobytes(),
            nplanes, period, itor)

        self._mesh_injected = True
        self._log('Mesh injected: {} elements, {} planes'.format(
            nelms, nplanes))

    def inject_field(self, isrc, timeslice, hdf5_name):
        """
        Stream the field coefficient array for one field at one timeslice
        from the remote host and inject it into the fio_py numpy source.

        hdf5_name: the name of the field dataset inside the HDF5 fields group
                   (e.g. 'psi', 'den', 'te').
        """
        ts_file, group = self._ts_file_and_group(timeslice)
        real_path = '{}/fields/{}'.format(group, hdf5_name)
        self._log('Injecting field {} timeslice {} from {}'.format(
            hdf5_name, timeslice, ts_file))

        real_data = self.read(real_path, basename=ts_file)
        real_flat = np.ascontiguousarray(real_data.ravel(), dtype=np.float64)

        is_complex = bool(int(self.get_root_attr('icomplex') or 0))
        imag_flat = None
        if is_complex:
            imag_path = '{}/fields/{}_i'.format(group, hdf5_name)
            try:
                imag_data = self.read(imag_path, basename=ts_file)
                imag_flat = np.ascontiguousarray(
                    imag_data.ravel(), dtype=np.float64)
            except Exception:
                imag_flat = None

        ntor = int(self.get_root_attr('ntor') or 0)

        fio_py.set_numpy_field(
            isrc, timeslice, hdf5_name,
            real_flat.tobytes(),
            imag_flat.tobytes() if imag_flat is not None else None,
            ntor if imag_flat is not None else 0)

    def inject_scalars(self, isrc):
        """
        Stream all scalar time series from the remote C1.h5 and inject them
        into the fio_py numpy source. Also injects 'slice_times'.
        """
        scalar_names = self.list_group('scalars')
        for name in scalar_names:
            try:
                arr = self.read('scalars/{}'.format(name))
                fio_py.set_numpy_scalar(
                    isrc, name,
                    np.ascontiguousarray(arr, dtype=np.float64).tobytes())
            except Exception as e:
                self._log(
                    'Could not inject scalar {}: {}'.format(name, e))

        # Inject slice_times for series evaluation.
        ntime = int(self.get_root_attr('ntime') or 0)
        slice_times = []
        for t in range(ntime):
            tv = self.get_group_attr('time_{:03d}'.format(t), 'time')
            slice_times.append(float(tv) if tv is not None else 0.0)
        if slice_times:
            arr = np.array(slice_times, dtype=np.float64)
            fio_py.set_numpy_scalar(
                isrc, 'slice_times',
                np.ascontiguousarray(arr).tobytes())

    def close(self):
        self._con.close()

    def __del__(self):
        try:
            self._con.close()
        except Exception:
            pass


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
    The time traces are read directly from the C1.h5 file for local mode.
    For campaign archive mode they are streamed from the remote host.

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
    number of iterations, timings for each solve, and time_slice - time
    correspondensies.
    The diagnostic object has two attributes:
    !! Diagnostic are given in M3D-C1 units !!
    diagnostic.diagnostic - contains the value of the diagnostic. Can be a
                            simple array (i.e. for 'time slices', or several
                            objects i.e. to see a breakdown for timings)
    diagnostic.x_axis     - contains the values of the x-axis of the
                            diagnostic. This is alfven times for 'time slices'
                            and iteration number for 'timings'
    To see all available diagnostics, check sim.available_diagnostics

    --Constants--
    On subclass containing the constants can be created by,
    for example,
    constants = sim.get_constants()
    In this subclass, stuff like gamma, R0, and the version
    number of the M3DC1 simulation are stored.

    --Campaign Archive Mode (.aca files)--
    To stream M3D-C1 data from a remote cluster via an HPC campaign
    archive (.aca file), use filetype='adios2':

        sim = sim_data(
            filename='mastu_45272_vped477.aca',
            filetype='adios2',
            campaign_simulation='99/1f_eqrotnc-C_eta_x1/n40',
            time=-1,
            login='user@perlmutter.nersc.gov',
        )

    Field coefficient arrays and mesh data are streamed on demand from the
    remote host via SSH + h5py and injected directly into the fusion-io C++
    finite element evaluation machinery. No ADIOS2 installation is required.
    Only the data actually needed for the requested field and time slice is
    transferred.
    """

    def __init__(self, filename='C1.h5', filetype='m3dc1', verbose=False,
                 time=0, fast=False, campaign_simulation='',
                 login='', cache_dir='.cache/aca',
                 remote_python='python3',
                 remote_setup='module load python >/dev/null 2>&1 || true'):
        """
        Initializes the fusion-io bindings to a file.

        Keyworded arguments:

        **filename**
            the name of the file which is to be read. Can include path.
            For filetype='adios2', this should be the path to the .aca
            campaign archive file.

        **filetype**
            the type of file. One of:
            - 'm3dc1'  : local M3D-C1 HDF5 file (default)
            - 'nimrod' : local NIMROD file
            - 'adios2' : M3D-C1 data accessed via an HPC campaign archive
                         (.aca file). Field coefficients and mesh arrays are
                         streamed via SSH + h5py and injected into the C++
                         fusion-io finite element evaluator. No ADIOS2
                         installation is required despite the name.

        **verbose**
            if true, gives extra info on file read.

        **time**
            sets the time variable during initialization.

        **fast**
            Workaround to keep mesh connectivity in memory.

        **campaign_simulation**
            Only used when filetype='adios2'. Identifies which simulation
            within the campaign archive to stream.

        **login**
            Only used when filetype='adios2'. SSH login string for the remote
            host, e.g. 'user@perlmutter.nersc.gov'.

        **cache_dir**
            Directory for extracted skeleton files and cached .npy arrays.

        **remote_python**
            Python interpreter to invoke on the remote host.

        **remote_setup**
            Shell command run before the remote Python script.
        """
        self._isrc = None
        self.fields = []
        self.seriess = []
        self._aca = None
        self._filetype = filetype

        if filetype == 'm3dc1':
            ifiletype = fio_py.FIO_M3DC1_SOURCE
        elif filetype == 'nimrod':
            if getattr(fio_py, 'FIO_HAS_NIMROD', 0):
                ifiletype = fio_py.FIO_NIMROD_SOURCE
            else:
                raise ValueError(
                    'fusion-io was not compiled with NIMROD support.')
        elif filetype == 'adios2':
            ifiletype = None   # numpy source, opened below
        else:
            raise ValueError(
                'Unsupported filetype: {}. '
                'Supported: m3dc1, nimrod, adios2.'.format(filetype))

        # Dictionary contains the field abbreviation, the type of field,
        # and the species
        self.typedict = {
            'j'        : ('current density',       'vector', None,       'simple'),
            'ni'       : ('density',               'scalar', 'main ion', 'simple'),
            'ne'       : ('density',               'scalar', 'electron', 'simple'),
            'v'        : ('fluid velocity',        'vector', None,       'simple'),
            'B'        : ('magnetic field',        'vector', None,       'simple'),
            'p'        : ('total pressure',        'scalar', None,       'simple'),
            'pi'       : ('pressure',              'scalar', 'main ion', 'simple'),
            'pe'       : ('pressure',              'scalar', 'electron', 'simple'),
            'alpha'    : ('alpha',                 'scalar', None,       'simple'),
            'ti'       : ('temperature',           'scalar', 'main ion', 'simple'),
            'te'       : ('temperature',           'scalar', 'electron', 'simple'),
            'A'        : ('vector potential',      'vector', None,       'simple'),
            'gradA'    : ('grad vector potential', 'tensor', None,       'simple'),
            'E'        : ('electric field',        'vector', None,       'simple'),
            'psi'      : ('psi',                   'scalar', None,       'composite'),
            'kprad_rad': ('total radiation',       'scalar', None,       'simple'),
        }
        self.available_fields = self.typedict
        self.filename = os.path.abspath(filename)

        if filetype == 'adios2':
            if not _have_h5py:
                raise ImportError(
                    'h5py is required for campaign archive mode. '
                    'Install with: pip install h5py')
            self._aca = _AcaBackend(
                aca_filename=self.filename,
                campaign_simulation=campaign_simulation,
                login=login,
                cache_dir=cache_dir,
                remote_python=remote_python,
                remote_setup=remote_setup,
                verbose=verbose)

            if verbose:
                print('Campaign archive: streaming from simulation '
                      '"{}" in {}'.format(campaign_simulation, self.filename))

            # Open a numpy source and inject parameters + scalars.
            self._isrc = fio_py.open_numpy_source()
            self._aca.inject_parameters(self._isrc)
            self._aca.inject_scalars(self._isrc)

            ntime_attr = self._aca.get_root_attr('ntime')
            self.ntime = int(ntime_attr) if ntime_attr is not None else 0
            self.available_traces = self._aca.list_group('scalars')
            self.available_traces.extend(['bharmonics', 'keharmonics'])
            self._all_attrs = None
            self._all_traces = None

        else:
            if not _have_h5py:
                raise ImportError(
                    'h5py is required to read local M3D-C1 files. '
                    'Install with: pip install h5py')
            self._all_attrs      = h5py.File(self.filename, 'r')
            self._all_attrs_list = list(self._all_attrs.keys())
            self._all_traces     = self._all_attrs['scalars']
            self.available_traces = list(self._all_traces.keys())
            self.available_traces.extend(['bharmonics', 'keharmonics'])
            self.ntime = self._all_attrs.attrs['ntime']
            self._isrc = fio_py.open_source(ifiletype, self.filename)

        self.available_diagnostics = ['slice times', 'timings', 'iterations']
        self.hint = fio_py.allocate_hint(self._isrc)
        fio_py.get_options(self._isrc)

        # For campaign mode, inject mesh before querying geometry parameters.
        time_int = int(time) if not isinstance(time, str) else -1
        if filetype == 'adios2':
            self._aca.inject_mesh(self._isrc, time_int)

        self.set_timeslice(time)
        if fast:
            self._imag = fio_py.get_field(self._isrc, fio_py.FIO_MAGNETIC_FIELD)
            self.fields.append(self._imag)
        self._cs = fio_py.get_int_parameter(self._isrc, fio_py.FIO_GEOMETRY)
        self._iavailable_fields = fio_py.get_available_fields(self._isrc)
        # available fields is a dictionary from names to assigned integers
        self._available_fields = dict(zip(
            [fio_py.get_field_name(nr) for nr in self._iavailable_fields],
            self._iavailable_fields))
        self.ntor = fio_py.get_int_parameter(self._isrc, fio_py.FIO_TOROIDAL_MODE)
        self.period = fio_py.get_real_parameter(self._isrc, fio_py.FIO_PERIOD)
        self.fc = None

        if verbose:
            print('Available fields:')
            if self._cs == fio_py.FIO_CYLINDRICAL:
                print('Using CYLINDRICAL coordinate system')
            else:
                print('Using CARTESIAN coordinate system')
            print('Number of time slices: ', self.ntime)
            print('Toroidal period = ', self.period)

    def __del__(self):
        print('deleting simulation object and closing {} fields'.format(
            len(self.fields) + len(self.seriess)))
        for ifield in self.fields:
            fio_py.close_field(ifield)
        for iseries in self.seriess:
            fio_py.close_series(iseries)
        if self._isrc is not None:
            fio_py.close_source(self._isrc)
        if self._aca is not None:
            self._aca.close()

    def set_timeslice(self, time):
        if isinstance(time, str):
            if time == 'last':
                time = self.ntime - 1
                print('last time slice = ' + str(time))
        time = int(time)
        self.timeslice = time
        fio_py.set_int_option(fio_py.FIO_TIMESLICE, time)

    def get_time_trace(self, scalar, ipellet=None):
        return self.time_trace(scalar, sim_data=self, ipellet=ipellet)

    def get_diagnostic(self, diagnostic):
        return self.diagnostic(self, diagnostic)

    def get_field(self, field, time=None):
        if field in self.typedict.keys():
            return self.field(self, field=self.typedict[field][0], time=time,
                              species=self.typedict[field][2],
                              ftype=self.typedict[field][1])
        else:
            return self.field(self, field=field, time=time, species=None,
                              ftype='scalar')

    def get_mesh(self, time=None, quiet=False):
        return self.mesh(self, time, quiet=quiet)

    def get_signal(self, signame):
        return self.signal(self, signame)

    def get_constants(self):
        return self.constants(self)

    def _get_attr(self, group, name):
        """Read a scalar attribute, dispatching to h5py or _AcaBackend."""
        if self._filetype == 'adios2':
            if group:
                return self._aca.get_group_attr(group, name)
            else:
                return self._aca.get_root_attr(name)
        else:
            if group:
                return self._all_attrs[group].attrs[name]
            else:
                return self._all_attrs.attrs[name]

    def _get_dataset(self, path):
        """Read a dataset as a numpy array."""
        if self._filetype == 'adios2':
            return self._aca.read(path)
        else:
            return np.asarray(self._all_attrs[path])

    def _list_group(self, group):
        """List the keys of a group."""
        if self._filetype == 'adios2':
            return self._aca.list_group(group)
        else:
            return list(self._all_attrs[group].keys())

    def _group_exists(self, group):
        """Check whether a group exists."""
        if self._filetype == 'adios2':
            return self._aca.group_exists(group)
        else:
            return group in self._all_attrs

    def _ensure_field_injected(self, hdf5_name, timeslice):
        """
        Ensure that the named field coefficients for the given timeslice have
        been injected into the numpy source. Called by the field class before
        requesting field evaluation. For local files this is a no-op.
        """
        if self._filetype != 'adios2':
            return
        self._aca.inject_field(self._isrc, timeslice, hdf5_name)

    class field:
        """
        Field class: the init sets up the data access, and its bound methods
        return to you what you want. ex:
        myfield = mysim.get_field('...')
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
        def __init__(self, sim_data, field, time=None, species=None,
                     ftype='scalar'):
            self.sim_data = sim_data
            self.field = field
            if time is not None:
                self.sim_data.set_timeslice(time)
            if species == 'electron':
                fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_ELECTRON)
            elif species == 'main ion':
                fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_MAIN_ION)

            # For campaign mode, stream and inject the field coefficients
            # before calling get_field so the C++ evaluator has the data.
            if sim_data._filetype == 'adios2':
                hdf5_map = {
                    'psi'            : 'psi',
                    'magnetic field' : 'psi',
                    'current density': 'J',
                    'fluid velocity' : 'v',
                    'density'        : ('ne' if species == 'electron'
                                        else 'den'),
                    'pressure'       : ('Pe' if species == 'electron'
                                        else 'pi'),
                    'total pressure' : 'P',
                    'temperature'    : ('te' if species == 'electron'
                                        else 'ti'),
                    'vector potential': 'A',
                    'electric field' : 'E',
                    'total radiation': 'kprad_rad',
                    'alpha'          : 'alpha',
                }
                hdf5_name = hdf5_map.get(field, field)
                sim_data._ensure_field_injected(
                    hdf5_name, sim_data.timeslice)

            if field in sim_data._available_fields.keys():
                itype = sim_data._available_fields[field]
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
            Evaluates the required field, returns a one-tuple (scalar quantity
            required) or a three-tuple (vector) or a nine-tuple (tensor)
            evaluated at the location x.

            Arguments:

            **x**
                Three-tuple with the R, phi, z coordinate where the field is
                to be evaluated.
            """
            if self.ftype == 'vector':
                try:
                    return fio_py.eval_vector_field(
                        self._ifield, x, self.sim_data.hint)
                except Exception:
                    return (None, None, None)
            if self.ftype == 'tensor':
                try:
                    return fio_py.eval_tensor_field(
                        self._ifield, x, self.sim_data.hint)
                except Exception:
                    return (None,) * 9
            elif self.ftype == 'scalar':
                try:
                    return (fio_py.eval_scalar_field(
                        self._ifield, x, self.sim_data.hint),)
                except Exception:
                    return (None,)
            else:
                print('ftype not recognized!')

        def evaluate_deriv(self, x):
            """
            Evaluates the derivative of the required field, returns a
            three-tuple (scalar) or a nine-tuple (vector) evaluated at x.

            Arguments:

            **x**
                Three-tuple with the R, phi, z coordinate where the field is
                to be evaluated.
            """
            if self.ftype == 'vector':
                try:
                    return fio_py.eval_vector_field_deriv(
                        self._ifield, x, self.sim_data.hint)
                except Exception:
                    return (None,) * 9
            elif self.ftype == 'scalar':
                try:
                    return (fio_py.eval_scalar_field_deriv(
                        self._ifield, x, self.sim_data.hint),)
                except Exception:
                    return (None, None, None)
            else:
                print('ftype not recognized!')

    class mesh:
        """
        Mesh class: init routine calls read_mesh() that reads the mesh from
        the time slice file and stores the elements in an array. It also reads
        the output version and nplanes.
        """
        def __init__(self, sim_data, time=None, quiet=False):
            self.sim_data = sim_data
            self._quiet = quiet
            if time is not None:
                self.sim_data.set_timeslice(time)
            self.elements, self.version, self.nplanes = self.read_mesh()

        def __eq__(self, other):
            return (np.array_equal(self.elements, other.elements) and
                    (self.version == other.version) and
                    (self.nplanes == other.nplanes))

        def read_mesh(self):
            if self.sim_data.timeslice == -1:
                group = 'equilibrium'
            else:
                group = 'time_{:03d}'.format(self.sim_data.timeslice)

            if self.sim_data._filetype == 'adios2':
                mesh = self.sim_data._get_dataset(group + '/mesh/elements')
                version = self.sim_data._get_attr(group, 'version')
                if version is None:
                    version = self.sim_data.get_constants().version
                nplanes = self.sim_data._get_attr(group + '/mesh', 'nplanes')
            else:
                f = self.sim_data._all_attrs
                mesh = np.asarray(f[group + '/mesh/elements'])
                if 'version' in f[group].attrs:
                    version = f[group].attrs['version']
                else:
                    version = self.sim_data.get_constants().version
                nplanes = f[group + '/mesh'].attrs['nplanes']

            if not self._quiet:
                print('Output version: ' + str(version))
                print('Mesh shape: ' + str(mesh.shape))

            return mesh, version, nplanes

    class signal:
        """
        Signal class: for diagnostic signals from magnetic probes and flux
        loops. Not available in campaign archive mode.
        """
        def __init__(self, sim_data, signame):
            self.sim_data = sim_data
            if sim_data._filetype == 'adios2':
                raise RuntimeError(
                    'get_signal() is not available in campaign archive mode.')
            self._all_attrs = self.sim_data._all_attrs
            self._all_attrs_list = self.sim_data._all_attrs_list
            self.sigvalues = np.asarray(
                self.sim_data._all_attrs[signame + '/value'][:])

    class time_trace:
        """
        time_trace class: reads scalar time traces from C1.h5 (local) or
        from the remote host via SSH + h5py (campaign archive mode).
        """
        def __init__(self, scalar, time=None, sim_data=None, ipellet=None):
            if (time is None) and (sim_data is None):
                raise ValueError('time or sim_data must not be None')

            if time is None:
                if sim_data._filetype == 'adios2':
                    time = sim_data._aca.read('scalars/time')
                else:
                    time = np.asarray(sim_data._all_traces['time'])

            if isinstance(scalar, str):
                if scalar in ['bharmonics', 'keharmonics']:
                    values = sim_data._get_dataset(
                        '{}/{}'.format(scalar, scalar))
                elif sim_data._filetype == 'adios2':
                    if (sim_data._aca.group_exists('pellet') and
                            sim_data._aca.group_exists(
                                'pellet/{}'.format(scalar))):
                        values = sim_data._aca.read(
                            'pellet/{}'.format(scalar))
                        if ipellet is not None:
                            values = values[:, ipellet]
                    else:
                        values = sim_data._get_dataset(
                            'scalars/{}'.format(scalar))
                else:
                    if 'pellet/{}'.format(scalar) in sim_data._all_attrs:
                        values = np.asarray(
                            sim_data._all_attrs[
                                'pellet/{}'.format(scalar)])
                        if ipellet is not None:
                            values = values[:, ipellet]
                    else:
                        values = np.asarray(sim_data._all_traces[scalar])
            else:
                values = np.asarray(scalar)

            if len(time) != len(values):
                if abs(len(time) - len(values)) == 1:
                    print('WARNING: len(time) and len(values) differ by 1. '
                          'Simulation might currently be running.')
                    if len(time) > len(values):
                        self.time = time[:len(values)]
                        self.values = values
                else:
                    print('WARNING: time_trace time and values have different '
                          'lengths. len(time) = ' + str(len(time)) +
                          ', len(values) = ' + str(len(values)))
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
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                values = self.values + other.values
            else:
                values = self.values + other
            return sim_data.time_trace(values, time=self.time)

        __radd__ = __add__

        def __iadd__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                self.values = self.values + other.values
            else:
                self.values = self.values + other
            return self

        # Subtraction
        def __sub__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                values = self.values - other.values
            else:
                values = self.values - other
            return sim_data.time_trace(values, time=self.time)

        def __rsub__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                values = other.values - self.values
            else:
                values = other - self.values
            return sim_data.time_trace(values, time=self.time)

        def __isub__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                self.values = self.values - other.values
            else:
                self.values = self.values - other
            return self

        # Multiplication
        def __mul__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                values = self.values * other.values
            else:
                values = self.values * other
            return sim_data.time_trace(values, time=self.time)

        __rmul__ = __mul__

        def __imul__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                self.values = self.values * other.values
            else:
                self.values = self.values * other
            return self

        # Division
        def __truediv__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                values = self.values / other.values
            else:
                values = self.values / other
            return sim_data.time_trace(values, time=self.time)

        def __rtruediv__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                values = other.values / self.values
            else:
                values = other / self.values
            return sim_data.time_trace(values, time=self.time)

        def __itruediv__(self, other):
            if isinstance(other, type(self)):
                if not np.array_equal(self.time, other.time):
                    raise ValueError(
                        'operands share different .time attributes')
                self.values = self.values / other.values
            else:
                self.values = self.values / other
            return self

        # Exponentiation
        def __pow__(self, other):
            return sim_data.time_trace(self.values ** other, time=self.time)

        # __rpow__ is intentionally NotImplemented

        def __ipow__(self, other):
            self.values = self.values ** other
            return self

        def __abs__(self):
            return sim_data.time_trace(abs(self.values), time=self.time)

        def __pos__(self):
            return sim_data.time_trace(+self.values, time=self.time)

        def __neg__(self):
            return sim_data.time_trace(-self.values, time=self.time)

        def cum_int(self, nts):
            values = cumtrapz(self.values[:nts], self.time[:nts], initial=0.0)
            return sim_data.time_trace(values, time=self.time)

    class diagnostic:
        def __init__(self, sim_data, diagnostic):
            self.sim_data = sim_data

            # Case corresponding to knowing the alfven times of the time slices
            if diagnostic == 'slice times':
                self.diagnostic = np.arange(sim_data.ntime)
                self.x_axis = np.zeros(sim_data.ntime)
                for t in self.diagnostic:
                    self.x_axis[t] = sim_data._get_attr(
                        'time_{:03d}'.format(t), 'time')

            # Case corresponding to the timings of each iteration
            if diagnostic == 'timings':
                self.diagnostic = {}
                for timing in sim_data._list_group('timings'):
                    self.diagnostic[timing] = sim_data._get_dataset(
                        'timings/' + timing)
                self.x_axis = np.arange(
                    self.diagnostic['t_onestep'].shape[0]) + 1

            # Case corresponding to the iterations of each step
            if diagnostic == 'iterations':
                self.diagnostic = sim_data._get_dataset('kspits/kspits')
                self.x_axis = np.arange(
                    np.shape(self.diagnostic[:, 0])[0]) + 1

    class constants:
        def __init__(self, sim_data):
            self.sim_data = sim_data
            self.R0      = sim_data._get_attr('', 'rzero')
            self.B0      = sim_data._get_attr('', 'bzero')
            self.gamma   = sim_data._get_attr('', 'gam')
            self.amupar  = sim_data._get_attr('', 'amupar')
            self.version = sim_data._get_attr('', 'version')
            self.numvar  = sim_data._get_attr('', 'numvar')
            self.itor    = sim_data._get_attr('', 'itor')
            self.is3D    = sim_data._get_attr('', '3d')


class flux_coordinates:
    """
    Class that represents a flux surface coordinate system, e.g. PEST, Boozer
    or Hamada coordinates.
    """
    def __init__(self, m, n, rpath, zpath, axis, omega, psi, psin, period,
                 theta, jac, q, area, polarea, dV, fcoords, V, phi, itor,
                 r0, current, dpsi_dpsin, points):
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
        self.flux_pol = -period * (psi - psi[0])
        self.theta = theta
        self.j = jac
        self.q = q
        self.area = area
        self.polarea = polarea
        self.dV_dchi = dV
        self.fcoords = fcoords
        self.V = V
        self.flux_tor = phi
        self.phi_norm = phi / phi[n - 1]
        self.rho = np.sqrt(phi / phi[n - 1])
        self.period = period
        self.itor = itor
        self.r0 = r0
        self.current = current
        self.dpsi_dchi = dpsi_dpsin
        self.points = points
