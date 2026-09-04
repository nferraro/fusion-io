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
import tempfile
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

# Reads ALL datasets under a given HDF5 group and writes them as a single
# .npz archive (one array per dataset name) to stdout.
_REMOTE_READ_GROUP_SCRIPT = r"""
import io, sys
import h5py, numpy as np
path, group = sys.argv[1], sys.argv[2]
arrays = {}
with h5py.File(path, "r") as f:
    if group in f:
        grp = f[group]
        for name in grp.keys():
            try:
                arrays[name] = np.asarray(grp[name])
            except Exception:
                pass
buf = io.BytesIO()
np.savez(buf, **arrays)
sys.stdout.buffer.write(buf.getvalue())
sys.stdout.buffer.flush()
"""

# Lists all dataset names under a given HDF5 group and writes them as
# newline-separated text to stdout.  Used to discover available field names.
_REMOTE_LIST_GROUP_SCRIPT = r"""
import sys
import h5py
path, group = sys.argv[1], sys.argv[2]
with h5py.File(path, "r") as f:
    if group in f:
        for name in f[group].keys():
            sys.stdout.write(name + "\n")
"""

_NPY_MAGIC = b"\x93NUMPY"
_NPZ_MAGIC = b"PK"          # zip/npz files start with PK
_TIME_RE   = re.compile(r"time[_.-]?(\d+)\.h5$", re.IGNORECASE)

# Matches 'equilibrium' or 'time_NNN' as the first path component.
_TS_PREFIX_RE = re.compile(
    r'^(equilibrium|time_(\d+))(/|$)')

# Patterns in SSH/Lmod stderr that are noise and should not cause a failure.
_STDERR_NOISE_RES = [
    re.compile(r'Lmod has detected the following error'),
    re.compile(r'module\s+--ignore_cache'),
    re.compile(r'module spider'),
    re.compile(r'The following module\(s\) are unknown'),
    re.compile(r'Swap failed'),
    re.compile(r'automatically replacing'),
    re.compile(r'have been reloaded with a version change'),
    re.compile(r'PrgEnv-'),
    re.compile(r'cray-mpich'),
    re.compile(r'^\s*$'),           # blank lines
    re.compile(r'"gsl"'),           # gsl module not found
    re.compile(r'Please check the spelling or version number'),
    re.compile(r'Also try'),
    re.compile(r'it may help to try'),
    re.compile(r'\$ module'),
    re.compile(r'Also make sure'),
    re.compile(r'modulefiles written in TCL'),
    re.compile(r'#%Module'),
    re.compile(r'gpu/'),
    re.compile(r'cpu/'),
]

# ---------------------------------------------------------------------------
# HDF5 field dependency map
# ---------------------------------------------------------------------------
_FIELD_HDF5_DEPS = {
    'magnetic field'       : ('psi', 'I', 'f', 'fp'),
    'vector potential'     : ('psi', 'f'),
    'grad vector potential': ('psi', 'f'),
    'current density'      : ('psi', 'I', 'f', 'fp'),
    'fluid velocity'       : ('phi', 'V', 'chi'),
    'electric field'       : ('E_R', 'E_PHI', 'E_Z'),
    'pressure'             : ('Pe', 'P'),
}

_SIMPLE_FIELD_HDF5_NAME = {
    'psi'            : 'psi',
    'total pressure' : 'P',
    'total radiation': 'kprad_rad',
    'alpha'          : 'alpha',
}

_HDF5_FIELD_FALLBACKS = {
    'psi'   : ('Psi', 'flux', 'psi_n'),
    'I'     : ('i', 'Itor'),
    'f'     : ('F',),
    'fp'    : ('Fp', 'f_prime'),
    'phi'   : ('Phi', 'PHI'),
    'V'     : ('v_tor',),
    'chi'   : ('Chi',),
    'E_R'   : ('e_r',),
    'E_PHI' : ('e_phi', 'E_phi'),
    'E_Z'   : ('e_z',),
    'Pe'    : ('pe', 'P_e'),
    'P'     : ('p', 'ptot'),
    'ne'    : ('Ne', 'n_e'),
    'den'   : ('ni', 'n_i', 'density'),
    'te'    : ('Te', 'T_e'),
    'ti'    : ('Ti', 'T_i'),
}


def _required_attrs_for(basename):
    """
    Return the tuple of root attribute names that must be present in the
    cached attrs for *basename*.

    Timeslice files (equilibrium.h5, time_NNN.h5) must have 'ntimestep'.
    C1.h5 must have 'ntime'.
    If any are missing from the skeleton, they are fetched from the remote.
    """
    if basename == 'C1.h5':
        return ('ntime', 'bzero', 'version')
    else:
        return ('ntimestep', 'version')


def _get_fio_field_type(field_name):
    """
    Return the fio_py integer field-type constant for *field_name*, or None
    if the name is not a known composite field type.
    """
    _map = {
        'magnetic field'       : 'FIO_MAGNETIC_FIELD',
        'vector potential'     : 'FIO_VECTOR_POTENTIAL',
        'grad vector potential': 'FIO_GRAD_VECTOR_POTENTIAL',
        'current density'      : 'FIO_CURRENT_DENSITY',
        'fluid velocity'       : 'FIO_FLUID_VELOCITY',
        'electric field'       : 'FIO_ELECTRIC_FIELD',
        'pressure'             : 'FIO_PRESSURE',
        'total pressure'       : 'FIO_TOTAL_PRESSURE',
        'density'              : 'FIO_DENSITY',
        'temperature'          : 'FIO_TEMPERATURE',
        'psi'                  : 'FIO_POLOIDAL_FLUX',
        'total radiation'      : 'FIO_TOTAL_RADIATION',
        'alpha'                : 'FIO_ALPHA',
    }
    attr = _map.get(field_name)
    if attr is None:
        return None
    return getattr(fio_py, attr, None)


def _is_real_error(stderr_text):
    for line in stderr_text.splitlines():
        if not any(pat.search(line) for pat in _STDERR_NOISE_RES):
            return True
    return False


def _safe_get_field_name(nr):
    try:
        return fio_py.get_field_name(nr)
    except Exception:
        return None


def _split_ts_prefix(path):
    m = _TS_PREFIX_RE.match(path)
    if m is None:
        return None, path
    prefix = m.group(1)
    remainder = path[m.end():]
    ts_basename = prefix + '.h5'
    return ts_basename, remainder


# ---------------------------------------------------------------------------
# _LocalH5Proxy  –  h5py-compatible proxy for local M3D-C1 HDF5 files
# ---------------------------------------------------------------------------

class _LocalH5AttrsProxy:
    """
    Proxy for h5py AttributeManager that searches multiple HDF5 files.

    For root-level access, searches the timeslice file (e.g. equilibrium.h5
    or time_002.h5) first, then C1.h5.  This mirrors the behaviour of
    _AcaAttrsProxy so that attributes like 'ntimestep' (which live in the
    timeslice file, not in C1.h5) are found correctly.
    """

    def __init__(self, h5_files, group=''):
        """
        Parameters
        ----------
        h5_files : list of h5py.File
            Open h5py file handles to search, in priority order (first wins).
        group : str
            HDF5 group path whose attributes to expose.  Empty string or '/'
            means root attributes.
        """
        self._h5_files = h5_files
        self._group = group.strip('/')

    def _is_root(self):
        return self._group == ''

    def _get_attrs_obj(self, h5f):
        """Return the h5py AttributeManager for the relevant group/root."""
        if self._is_root():
            return h5f.attrs
        if self._group in h5f:
            return h5f[self._group].attrs
        return None

    def _merged(self):
        merged = {}
        for h5f in reversed(self._h5_files):
            try:
                ao = self._get_attrs_obj(h5f)
                if ao is not None:
                    merged.update({k: v for k, v in ao.items()})
            except Exception:
                pass
        return merged

    def __getitem__(self, key):
        for h5f in self._h5_files:
            try:
                ao = self._get_attrs_obj(h5f)
                if ao is not None and key in ao:
                    return ao[key]
            except Exception:
                pass
        raise KeyError(
            "Unable to synchronously open attribute "
            "(can't locate attribute: '{}')".format(key))

    def __contains__(self, key):
        for h5f in self._h5_files:
            try:
                ao = self._get_attrs_obj(h5f)
                if ao is not None and key in ao:
                    return True
            except Exception:
                pass
        return False

    def keys(self):
        return self._merged().keys()

    def items(self):
        return self._merged().items()


class _LocalH5Proxy:
    """
    Top-level proxy that mimics an h5py File handle for local M3D-C1 data.

    Wraps both the main C1.h5 file and the current timeslice file
    (e.g. equilibrium.h5 or time_002.h5).  Root attributes (.attrs) are
    searched in the timeslice file first, then C1.h5, so that attributes
    like 'ntimestep' (which live in the timeslice file) are found correctly.

    Dataset/group access delegates to whichever open file contains the path,
    trying the timeslice file first, then C1.h5.

    This mirrors the behaviour of _AcaH5Proxy for remote/campaign data so
    that analysis code works identically for both local and remote files.
    """

    def __init__(self, c1_file, ts_file=None):
        """
        Parameters
        ----------
        c1_file : h5py.File
            Open handle to C1.h5.
        ts_file : h5py.File or None
            Open handle to the timeslice file (equilibrium.h5 / time_NNN.h5).
            If None, only C1.h5 is used.
        """
        self._c1 = c1_file
        self._ts = ts_file
        # Priority order: timeslice first, then C1.
        self._files = [f for f in [ts_file, c1_file] if f is not None]

    @property
    def filename(self):
        return self._c1.filename

    @property
    def attrs(self):
        return _LocalH5AttrsProxy(self._files, group='')

    def _resolve(self, path):
        """Return the first open file that contains *path*."""
        path = path.strip('/')
        for h5f in self._files:
            try:
                if path in h5f:
                    return h5f, path
            except Exception:
                pass
        # Fall back to C1.h5 even if path not found (will raise naturally).
        return self._c1, path

    def __getitem__(self, path):
        h5f, inner = self._resolve(path)
        return h5f[inner]

    def __contains__(self, path):
        path = path.strip('/')
        for h5f in self._files:
            try:
                if path in h5f:
                    return True
            except Exception:
                pass
        return False

    def keys(self):
        seen = set()
        result = []
        for h5f in self._files:
            try:
                for k in h5f.keys():
                    if k not in seen:
                        seen.add(k)
                        result.append(k)
            except Exception:
                pass
        return result

    def visititems(self, func):
        """Mimic h5py File.visititems — visits items in the timeslice file."""
        primary = self._ts if self._ts is not None else self._c1
        try:
            primary.visititems(func)
        except Exception:
            pass

    def update_timeslice(self, ts_file):
        """
        Replace the timeslice file handle.

        Called by sim_data.set_timeslice() when the user switches time steps,
        so that _all_attrs always reflects the current timeslice.
        """
        self._ts = ts_file
        self._files = [f for f in [ts_file, self._c1] if f is not None]


# ---------------------------------------------------------------------------
# _AcaH5Proxy  –  h5py-compatible proxy for campaign archive data
# ---------------------------------------------------------------------------

class _AcaAttrsProxy:
    """
    Proxy for h5py AttributeManager.

    Supports: attrs[key], key in attrs, attrs.keys(), attrs.items()

    Routes attribute reads through _AcaBackend, trying each of the provided
    basenames in order until the attribute is found.

    For root-level access (group == '' or '/'), delegates to
    _aca._fetch_root_attrs(basename) so that the remote fallback logic is
    triggered when required attributes (e.g. 'ntimestep') are missing from
    the embedded skeleton.
    """
    def __init__(self, aca, basenames, group):
        self._aca = aca
        if isinstance(basenames, str):
            basenames = [basenames]
        self._basenames = basenames
        self._group = group

    def _is_root(self):
        return self._group in ('', '/')

    def _fetch_from(self, basename):
        """Fetch attrs dict from one file."""
        if self._is_root():
            # Route through _fetch_root_attrs so the remote fallback for
            # missing required attributes (e.g. 'ntimestep') is triggered.
            return self._aca._fetch_root_attrs(basename)
        else:
            skeleton = self._aca._skeleton(basename)
            if _have_h5py:
                try:
                    with h5py.File(skeleton, 'r') as f:
                        if self._group in f:
                            def plain(v):
                                try:
                                    return (v.item()
                                            if getattr(v, 'shape', None) == ()
                                            else v.tolist())
                                except AttributeError:
                                    return (v.decode()
                                            if isinstance(v, bytes) else v)
                            return {k: plain(v)
                                    for k, v in f[self._group].attrs.items()}
                except Exception:
                    pass
            return {}

    def _fetch(self):
        merged = {}
        for basename in reversed(self._basenames):
            try:
                merged.update(self._fetch_from(basename))
            except Exception:
                pass
        return merged

    def __getitem__(self, key):
        if self._is_root():
            # Use _fetch_root_attrs (with remote fallback) for each basename.
            for basename in self._basenames:
                try:
                    d = self._aca._fetch_root_attrs(basename)
                    if key in d:
                        return d[key]
                except Exception:
                    pass
        else:
            for basename in self._basenames:
                try:
                    d = self._fetch_from(basename)
                    if key in d:
                        return d[key]
                except Exception:
                    pass
        raise KeyError(
            "Unable to synchronously open attribute "
            "(can't locate attribute: '{}')".format(key))

    def __contains__(self, key):
        if self._is_root():
            for basename in self._basenames:
                try:
                    if key in self._aca._fetch_root_attrs(basename):
                        return True
                except Exception:
                    pass
            return False
        else:
            for basename in self._basenames:
                try:
                    if key in self._fetch_from(basename):
                        return True
                except Exception:
                    pass
            return False

    def keys(self):
        return self._fetch().keys()

    def items(self):
        return self._fetch().items()


class _AcaGroupProxy:
    """
    Proxy for an h5py Group or Dataset inside a campaign archive.

    Supports:
      - group[key]          -> _AcaGroupProxy or numpy array
      - key in group        -> bool
      - group.keys()        -> list of child names
      - group.attrs         -> _AcaAttrsProxy
      - group.visititems()  -> like h5py Group.visititems
      - group.filename      -> str  (for error messages)
    """
    def __init__(self, aca, basename, path):
        self._aca = aca
        self._basename = basename
        self._path = path.strip('/')

    @property
    def filename(self):
        return self._basename

    @property
    def attrs(self):
        return _AcaAttrsProxy(self._aca, self._basename, self._path)

    def _child_path(self, key):
        if self._path:
            return self._path + '/' + key
        return key

    def __getitem__(self, key):
        child_path = self._child_path(key)
        try:
            data = self._aca.read(child_path, basename=self._basename)
            return data
        except Exception:
            pass
        return _AcaGroupProxy(self._aca, self._basename, child_path)

    def __contains__(self, key):
        child_path = self._child_path(key)
        skeleton = self._aca._skeleton(self._basename)
        if _have_h5py:
            try:
                with h5py.File(skeleton, 'r') as f:
                    return child_path in f
            except Exception:
                pass
        return False

    def keys(self):
        skeleton = self._aca._skeleton(self._basename)
        if _have_h5py:
            try:
                with h5py.File(skeleton, 'r') as f:
                    if self._path in f:
                        return list(f[self._path].keys())
                    elif self._path == '':
                        return list(f.keys())
            except Exception:
                pass
        return []

    def visititems(self, func):
        """Mimic h5py Group.visititems — calls func(name, obj) for all items."""
        skeleton = self._aca._skeleton(self._basename)
        if _have_h5py:
            try:
                with h5py.File(skeleton, 'r') as f:
                    root = f[self._path] if self._path else f
                    root.visititems(func)
            except Exception:
                pass


class _AcaH5Proxy:
    """
    Top-level proxy that mimics an h5py File handle for a campaign archive.

    This allows existing analysis code that does:
        h5file = sim._all_attrs
        h5file.attrs['version']        # root attribute from C1.h5
        h5file.attrs['ntimestep']      # root attribute from equilibrium.h5
        h5file['scalars/time']         # dataset from C1.h5
        'scalars' in h5file
        h5file.visititems(print)

    to work transparently whether the data is local or remote, without any
    changes to the analysis scripts.

    Root attributes (.attrs) searches both the timeslice file
    (e.g. equilibrium.h5) AND C1.h5, so that both 'ntimestep' (from
    equilibrium.h5) and 'bzero', 'ntime', etc. (from C1.h5) are accessible
    via sim._all_attrs.attrs[key] without any changes to calling code.

    If the embedded skeleton is missing required attributes (e.g. 'ntimestep'),
    _fetch_root_attrs() automatically fetches them from the remote host.
    """

    def __init__(self, aca, timeslice_basename='equilibrium.h5',
                 scalar_basename='C1.h5'):
        self._aca = aca
        self._ts_basename = timeslice_basename
        self._sc_basename = scalar_basename

    def _resolve(self, path):
        path = path.strip('/')
        if path.startswith('scalars') or path.startswith('mag_probes'):
            return self._sc_basename, path
        basename, remainder = self._aca._route(path)
        return basename, remainder

    @property
    def filename(self):
        return self._ts_basename

    @property
    def attrs(self):
        # Search timeslice file first (for 'ntimestep', 'time', etc.),
        # then C1.h5 (for 'bzero', 'ntime', 'version', etc.).
        # _AcaAttrsProxy will route root-level access through
        # _fetch_root_attrs (with remote fallback) for each basename.
        return _AcaAttrsProxy(
            self._aca,
            [self._ts_basename, self._sc_basename],
            '')

    def __getitem__(self, path):
        basename, inner_path = self._resolve(path)
        try:
            data = self._aca.read(inner_path, basename=basename)
            return data
        except Exception:
            pass
        return _AcaGroupProxy(self._aca, basename, inner_path)

    def __contains__(self, path):
        basename, inner_path = self._resolve(path)
        skeleton = self._aca._skeleton(basename)
        if _have_h5py:
            try:
                with h5py.File(skeleton, 'r') as f:
                    return inner_path in f
            except Exception:
                pass
        return False

    def keys(self):
        skeleton = self._aca._skeleton(self._ts_basename)
        if _have_h5py:
            try:
                with h5py.File(skeleton, 'r') as f:
                    return list(f.keys())
            except Exception:
                pass
        return []

    def visititems(self, func):
        """Mimic h5py File.visititems."""
        skeleton = self._aca._skeleton(self._ts_basename)
        if _have_h5py:
            try:
                with h5py.File(skeleton, 'r') as f:
                    f.visititems(func)
            except Exception:
                pass


class _SshMux:
    """
    Manages a single persistent SSH ControlMaster connection so that all
    subsequent SSH calls reuse the same authenticated session.

    If *login* is empty the mux is a no-op.
    """

    def __init__(self, login):
        self._login = login
        self._socket = None
        self._proc   = None
        self._tmpdir = None

    def __enter__(self):
        if not self._login:
            return self
        self._tmpdir = tempfile.mkdtemp(prefix='fpy_ssh_mux_')
        self._socket = os.path.join(self._tmpdir, 'ctl')
        cmd = [
            'ssh',
            '-o', 'ControlMaster=yes',
            '-o', 'ControlPath={}'.format(self._socket),
            '-o', 'ControlPersist=yes',
            '-o', 'BatchMode=no',
            '-N',
            self._login,
        ]
        self._proc = subprocess.Popen(
            cmd,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        import time as _time
        for _ in range(60):
            if os.path.exists(self._socket):
                break
            _time.sleep(0.5)
        return self

    def __exit__(self, *_):
        if self._socket and os.path.exists(self._socket):
            subprocess.run(
                ['ssh',
                 '-o', 'ControlPath={}'.format(self._socket),
                 '-O', 'exit',
                 self._login],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        if self._proc is not None:
            try:
                self._proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._proc.kill()
        if self._tmpdir:
            import shutil
            shutil.rmtree(self._tmpdir, ignore_errors=True)

    def ssh_opts(self):
        if self._socket:
            return [
                '-o', 'ControlMaster=no',
                '-o', 'ControlPath={}'.format(self._socket),
            ]
        return []


class _AcaBackend:
    """
    Backend for reading M3D-C1 data from a remote cluster via an HPC campaign
    archive (.aca file).

    Skeleton selection
    ------------------
    Each dataset in the .aca archive may have multiple embedded files.  We
    select the embedded file whose basename matches the dataset basename
    (e.g. 'equilibrium.h5').  If no such file is found we fall back to any
    embedded .h5 file, then to the first embedded file.

    Root attribute completeness
    ---------------------------
    The embedded skeleton may be missing root attributes that are present in
    the actual remote file (e.g. 'ntimestep' in equilibrium.h5).  After
    reading the skeleton, _fetch_root_attrs() checks whether all required
    attributes are present.  If any are missing, it fetches the complete set
    of root attributes from the remote host via SSH and merges them into the
    cached dict.  This ensures that readParameter('ntimestep', h5file=sim._all_attrs)
    always works correctly without any changes to calling code.
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
        if row is None:
            raise RuntimeError(
                'Campaign archive info table is empty: {}'.format(aca_filename))
        self._archive_id, self._archive_name, self._version = row

        self._simulation = campaign_simulation.strip('/')
        self._datasets = self._load_datasets()

        prefix = self._simulation + '/' if self._simulation else ''
        self._files = {}
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

        self._skeletons = {}
        self._root_attrs_cache = {}
        self._fields_cache = {}
        self._mesh_injected = False
        self._injected_fields = set()

        self._mux = _SshMux(login)
        self._mux.__enter__()

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
        for row in cur.execute(
                'SELECT rowid, name, uuid, fileformat, deltime FROM dataset'):
            did, name, uuid, fmt, deltime = row
            if deltime:
                continue
            datasets[did] = {
                'id': did, 'uuid': uuid or name, 'name': name,
                'fileformat': fmt or '', 'replicas': [],
                'embedded_files': []}

        replica_rowids = {}
        for row in cur.execute(
                'SELECT rowid, datasetid, hostid, dirid, name, '
                'deltime, size FROM replica'):
            rid, dsid, hostid, dirid, rname, deltime, size = row
            ds = datasets.get(dsid)
            if ds is None or deltime:
                continue
            hostname, longhost = hosts.get(hostid, ('?', '?'))
            _, dirname = dirs.get(dirid, (None, '?'))
            rep_entry = {
                'id': rid,
                'hostname': hostname,
                'longhost': longhost,
                'directory': dirname,
                'name': rname,
                'size': size or 0,
                'remote_path': os.path.join(dirname or '', rname or '')}
            ds['replicas'].append(rep_entry)
            replica_rowids[rid] = ds

        for row in cur.execute(
                'SELECT rf.replicaid, f.fileid, f.name, f.lenorig, '
                'f.lencompressed '
                'FROM repfiles rf JOIN file f ON rf.fileid = f.fileid'):
            replicaid, fid, fname, lenorig, lencomp = row
            ds = replica_rowids.get(replicaid)
            if ds is None:
                continue
            if all(e['id'] != fid for e in ds['embedded_files']):
                ds['embedded_files'].append({
                    'id': fid, 'name': fname,
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

    def _pick_skeleton_embedded_file(self, ds, basename):
        ef_list = ds['embedded_files']
        if not ef_list:
            return None
        target = basename.lower()
        for ef in ef_list:
            if os.path.basename(ef['name']).lower() == target:
                return ef
        for ef in ef_list:
            if ef['name'].lower().endswith('.h5'):
                return ef
        return ef_list[0]

    def _skeleton(self, basename=None):
        basename = basename or self._default_file
        ds = self._get_dataset(basename)
        if ds['name'] not in self._skeletons:
            cur = self._con.cursor()
            ef = self._pick_skeleton_embedded_file(ds, basename)
            if ef is None:
                raise RuntimeError(
                    'No embedded metadata in .aca for {}'.format(ds['name']))
            fid = ef['id']
            row = cur.execute(
                'SELECT data, compression FROM file '
                'WHERE fileid = ?', (fid,)).fetchone()
            if row is None or row[0] is None:
                raise RuntimeError(
                    'No embedded data for fileid {} in {}'.format(
                        fid, ds['name']))
            blob, compression = row
            data = zlib.decompress(blob) if compression else bytes(blob)
            dest = os.path.join(
                self._cache_dir, 'skeletons', ds['name'], basename)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, 'wb') as fh:
                fh.write(data)
            self._skeletons[ds['name']] = dest

        skeleton_path = self._skeletons[ds['name']]

        # Also make the skeleton available in the CWD under its bare basename
        # so that code like read_h5.py that does h5py.File('equilibrium.h5')
        # can find it without knowing about the campaign archive.
        cwd_path = os.path.join(os.getcwd(), basename)
        if not os.path.exists(cwd_path):
            try:
                os.symlink(os.path.abspath(skeleton_path), cwd_path)
                self._log(
                    'skeleton symlink: {} -> {}'.format(
                        cwd_path, skeleton_path))
            except (OSError, NotImplementedError):
                import shutil
                shutil.copy2(skeleton_path, cwd_path)
                self._log(
                    'skeleton copy: {} -> {}'.format(
                        skeleton_path, cwd_path))

        return skeleton_path

    # ------------------------------------------------------------------
    # SSH helpers
    # ------------------------------------------------------------------

    def _remote_run(self, script, argv, what):
        parts = []
        if self._remote_setup:
            parts.append(self._remote_setup)
        parts.append(
            'exec {python} - {args}'.format(
                python=self._remote_python,
                args=' '.join(shlex.quote(a) for a in argv)))
        inner = '; '.join(parts)

        if self._login:
            self._log('Connecting to {} for: {}'.format(self._login, what))
            cmd = (['ssh']
                   + self._mux.ssh_opts()
                   + [self._login, 'bash', '-lc', shlex.quote(inner)])
        else:
            cmd = ['bash', '-c', inner]

        proc = subprocess.run(
            cmd, input=script.encode(), capture_output=True, check=False)

        stderr_text = proc.stderr.decode(errors='replace')

        if proc.returncode != 0 and _is_real_error(stderr_text):
            raise RuntimeError(
                'Remote {} failed:\n{}'.format(what, stderr_text[-2000:]))

        if proc.returncode != 0:
            self._log(
                'WARNING: remote {} exited with code {} '
                '(stderr was only module-system noise, continuing)'.format(
                    what, proc.returncode))

        return proc.stdout

    def _var_cache_path(self, ds, variable, slc):
        key = hashlib.sha1(
            '{}:{}:{}'.format(
                ds['uuid'], variable, slc or '').encode()
        ).hexdigest()[:16]
        return os.path.join(self._cache_dir, 'vars', '{}.npy'.format(key))

    def _group_cache_path(self, ds, group):
        key = hashlib.sha1(
            '{}:group:{}'.format(ds['uuid'], group).encode()
        ).hexdigest()[:16]
        return os.path.join(self._cache_dir, 'vars', '{}.npz'.format(key))

    def _log(self, msg):
        if self._verbose:
            print('[fpy._AcaBackend] {}'.format(msg), file=sys.stderr)

    # ------------------------------------------------------------------
    # Path routing
    # ------------------------------------------------------------------

    def _route(self, path, basename=None):
        ts_basename, remainder = _split_ts_prefix(path.lstrip('/'))
        if ts_basename is not None:
            if ts_basename in self._files:
                return ts_basename, remainder or '.'
            return basename or self._default_file, path.lstrip('/')
        return basename or self._default_file, path.lstrip('/')

    # ------------------------------------------------------------------
    # Field name discovery
    # ------------------------------------------------------------------

    def list_fields(self, basename):
        if basename in self._fields_cache:
            return self._fields_cache[basename]

        names = set()

        if _have_h5py:
            try:
                with h5py.File(self._skeleton(basename), 'r') as f:
                    if 'fields' in f:
                        names = set(f['fields'].keys())
                        self._fields_cache[basename] = names
                        self._log(
                            'fields in {}: {} (from skeleton)'.format(
                                basename, sorted(names)))
                        return names
            except Exception:
                pass

        ds = self._get_dataset(basename)
        rep = self._live_replica(ds)
        out = self._remote_run(
            _REMOTE_LIST_GROUP_SCRIPT,
            [rep['remote_path'], 'fields'],
            'list fields in {}'.format(basename))
        names = {
            line.strip() for line in out.decode(errors='replace').splitlines()
            if line.strip()}
        self._fields_cache[basename] = names
        self._log('fields in {}: {}'.format(basename, sorted(names)))
        return names

    # ------------------------------------------------------------------
    # Array reading
    # ------------------------------------------------------------------

    def read(self, variable, basename=None, slc=None):
        eff_basename, variable = self._route(variable, basename)
        ds = self._get_dataset(eff_basename)

        if slc is None and _have_h5py:
            try:
                with h5py.File(self._skeleton(eff_basename), 'r') as f:
                    if variable in f and isinstance(f[variable], h5py.Dataset):
                        obj = f[variable]
                        if obj.size <= 128:
                            self._log(
                                '{}: from embedded skeleton'.format(variable))
                            return obj[()]
            except Exception:
                pass

        cache_file = self._var_cache_path(ds, variable, slc)
        if os.path.isfile(cache_file):
            self._log(
                '{}: from variable cache {}'.format(variable, cache_file))
            return np.load(cache_file, allow_pickle=False)

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

    def read_group(self, group, basename=None):
        eff_basename, group = self._route(group, basename)
        ds = self._get_dataset(eff_basename)

        cache_file = self._group_cache_path(ds, group)
        if os.path.isfile(cache_file):
            self._log('group {}: from cache {}'.format(group, cache_file))
            npz = np.load(cache_file, allow_pickle=False)
            return dict(npz)

        if _have_h5py:
            try:
                with h5py.File(self._skeleton(eff_basename), 'r') as f:
                    if group in f:
                        grp = f[group]
                        result = {}
                        all_small = True
                        for name in grp.keys():
                            obj = grp[name]
                            if isinstance(obj, h5py.Dataset) and obj.size <= 128:
                                result[name] = obj[()]
                            else:
                                all_small = False
                                break
                        if all_small and result:
                            self._log(
                                'group {}: all from skeleton'.format(group))
                            return result
            except Exception:
                pass

        rep = self._live_replica(ds)
        self._log('group {}: streaming from {}:{}'.format(
            group, rep['hostname'], rep['remote_path']))
        out = self._remote_run(
            _REMOTE_READ_GROUP_SCRIPT,
            [rep['remote_path'], group],
            "group read of '{}' from {}".format(group, rep['remote_path']))

        pos = out.find(_NPZ_MAGIC)
        if pos < 0:
            raise RuntimeError(
                "No npz data returned for group '{}'".format(group))
        npz = np.load(io.BytesIO(out[pos:]), allow_pickle=False)
        result = dict(npz)
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        np.savez(cache_file, **result)
        return result

    # ------------------------------------------------------------------
    # Attribute reading
    # ------------------------------------------------------------------

    def _fetch_root_attrs_remote(self, basename):
        """
        Fetch ALL root attributes directly from the remote host for *basename*.

        This bypasses the skeleton cache and always contacts the remote host.
        Used when the skeleton is missing required attributes (e.g. 'ntimestep'
        is present in the actual equilibrium.h5 but absent from the embedded
        skeleton in the .aca archive).
        """
        import json
        ds = self._get_dataset(basename)
        rep = self._live_replica(ds)
        self._log(
            'Fetching root attrs remotely for {} from {}'.format(
                basename, rep['remote_path']))
        out = self._remote_run(
            _ATTRS_SCRIPT,
            [rep['remote_path'], '/'],
            'root attrs of {}'.format(ds['name']))
        pos = out.find(b'###JSON###')
        if pos < 0:
            return {}
        return json.loads(out[pos + len(b'###JSON###'):].decode())

    def _fetch_root_attrs(self, basename=None):
        """
        Fetch root attributes for *basename*, with automatic remote fallback
        if the skeleton is missing required attributes.

        Strategy:
        1. Return from in-memory cache if already fetched.
        2. Try the embedded skeleton (fast, no SSH).
        3. Check whether all required attributes are present using
           _required_attrs_for(basename).
        4. If any are missing, fetch the complete attribute set from the
           remote host and merge into the result.  This handles the common
           case where the .aca skeleton was created without 'ntimestep'.
        5. Cache and return the merged result.
        """
        import json
        key = basename or self._default_file
        if key in self._root_attrs_cache:
            return self._root_attrs_cache[key]

        attrs = {}

        # Step 1: try skeleton.
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(key), 'r') as f:
                    def plain(v):
                        try:
                            return (v.item()
                                    if getattr(v, 'shape', None) == ()
                                    else v.tolist())
                        except AttributeError:
                            return (v.decode()
                                    if isinstance(v, bytes) else v)
                    attrs = {k: plain(v) for k, v in f.attrs.items()}
                    self._log(
                        'root attrs for {}: from skeleton'.format(key))
            except Exception:
                pass

        # Step 2: check completeness.
        required = _required_attrs_for(key)
        missing = [a for a in required if a not in attrs]

        if missing:
            # Step 3: fetch from remote and merge.
            self._log(
                'root attrs for {}: skeleton missing {}; '
                'fetching from remote'.format(key, missing))
            try:
                remote_attrs = self._fetch_root_attrs_remote(key)
                attrs.update(remote_attrs)
            except Exception as e:
                self._log(
                    'WARNING: could not fetch remote attrs for {}: {}'.format(
                        key, e))

        self._root_attrs_cache[key] = attrs
        return attrs

    def get_root_attr(self, name, basename=None):
        return self._fetch_root_attrs(basename).get(name)

    def get_group_attr(self, group, name, basename=None):
        eff_basename, eff_group = self._route(group, basename)
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(eff_basename), 'r') as f:
                    if eff_group in f and name in f[eff_group].attrs:
                        v = f[eff_group].attrs[name]
                        return v.item() if hasattr(v, 'item') else v
            except Exception:
                pass
        import json
        ds = self._get_dataset(eff_basename)
        rep = self._live_replica(ds)
        out = self._remote_run(
            _ATTRS_SCRIPT,
            [rep['remote_path'], eff_group],
            'attrs of {} in {}'.format(eff_group, ds['name']))
        pos = out.find(b'###JSON###')
        if pos < 0:
            return None
        attrs = json.loads(out[pos + len(b'###JSON###'):].decode())
        return attrs.get(name)

    def list_group(self, group, basename=None):
        eff_basename, eff_group = self._route(group, basename)
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(eff_basename), 'r') as f:
                    if eff_group in f:
                        return list(f[eff_group].keys())
            except Exception:
                pass
        return []

    def group_exists(self, group, basename=None):
        eff_basename, eff_group = self._route(group, basename)
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(eff_basename), 'r') as f:
                    return eff_group in f
            except Exception:
                pass
        return False

    def get_skeleton_path(self, basename):
        return self._skeleton(basename)

    def ensure_all_skeletons_in_cwd(self):
        for basename in self._files:
            if basename == 'C1.h5':
                continue
            try:
                self._skeleton(basename)
            except Exception as e:
                self._log(
                    'WARNING: could not extract skeleton for {}: {}'.format(
                        basename, e))

    # ------------------------------------------------------------------
    # Injection helpers
    # ------------------------------------------------------------------

    def _ts_file_and_group(self, timeslice):
        if timeslice == -1:
            ts_file = 'equilibrium.h5'
        else:
            ts_file = 'time_{:03d}.h5'.format(timeslice)

        if ts_file not in self._files:
            if self._timeslices:
                ts_file = self._timeslices[-1][1]
            else:
                ts_file = self._default_file

        return ts_file, ''

    def inject_parameters(self, isrc):
        attrs = self._fetch_root_attrs()

        def _f(name, default=0.0):
            v = attrs.get(name)
            return float(v) if v is not None else default

        def _i(name, default=0):
            v = attrs.get(name)
            return int(v) if v is not None else default

        fio_py.set_numpy_parameters(
            isrc,
            _f('bzero'), _f('rzero'),
            _i('ntime'), _i('ntor'),
            _f('ion_mass', 1.0), _f('z_ion', 1.0),
            _i('itor', 1), _i('numvar'), _i('version'),
            _f('n0_norm'), _f('l0_norm'), _f('b0_norm'))

    def inject_mesh(self, isrc, timeslice):
        if self._mesh_injected:
            return

        ts_file, _ = self._ts_file_and_group(timeslice)
        self._log('Injecting mesh from {}'.format(ts_file))

        elements = self.read('mesh/elements', basename=ts_file)
        version = int(self.get_root_attr('version') or 0)
        nelms = elements.shape[0]

        a  = np.ascontiguousarray(elements[:, 0], dtype=np.float64)
        b  = np.ascontiguousarray(elements[:, 1], dtype=np.float64)
        c  = np.ascontiguousarray(elements[:, 2], dtype=np.float64)
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

        nplanes_attr = self.get_group_attr('mesh', 'nplanes', basename=ts_file)
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
        ts_file, _ = self._ts_file_and_group(timeslice)

        injection_key = (ts_file, hdf5_name)
        if injection_key in self._injected_fields:
            self._log(
                "field '{}' already injected from {}; skipping".format(
                    hdf5_name, ts_file))
            return

        fallbacks = _HDF5_FIELD_FALLBACKS.get(hdf5_name, ())

        available = self.list_fields(ts_file)
        actual_name = None
        for candidate in (hdf5_name,) + tuple(fallbacks):
            if candidate in available:
                actual_name = candidate
                if candidate != hdf5_name:
                    self._log(
                        "field '{}' not found in {}; using '{}' instead".format(
                            hdf5_name, ts_file, candidate))
                break

        if actual_name is None:
            self._log(
                "WARNING: field '{}' (and fallbacks {}) not found in '{}'; "
                "skipping injection".format(hdf5_name, fallbacks, ts_file))
            return

        real_path = 'fields/{}'.format(actual_name)
        self._log("Injecting HDF5 field '{}' (as '{}') timeslice {} "
                  "from {}".format(actual_name, hdf5_name, timeslice, ts_file))

        real_data = self.read(real_path, basename=ts_file)
        real_flat = np.ascontiguousarray(real_data.ravel(), dtype=np.float64)

        is_complex = bool(int(self.get_root_attr('icomplex') or 0))
        imag_flat = None
        if is_complex:
            imag_path = 'fields/{}_i'.format(actual_name)
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

        self._injected_fields.add(injection_key)

    def inject_scalars(self, isrc):
        self._log('Fetching all scalars in one SSH call')
        try:
            scalar_arrays = self.read_group('scalars')
        except Exception as e:
            self._log('Could not batch-read scalars: {}'.format(e))
            scalar_arrays = {}

        for name, arr in scalar_arrays.items():
            try:
                fio_py.set_numpy_scalar(
                    isrc, name,
                    np.ascontiguousarray(arr, dtype=np.float64).tobytes())
            except Exception as e:
                self._log(
                    'Could not inject scalar {}: {}'.format(name, e))

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
        self._mux.__exit__(None, None, None)
        self._con.close()

    def __del__(self):
        try:
            self._mux.__exit__(None, None, None)
        except Exception:
            pass
        try:
            self._con.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Helper: open a local timeslice HDF5 file given a timeslice index
# ---------------------------------------------------------------------------

def _local_ts_filename(c1_dir, timeslice):
    """
    Return the path to the local timeslice HDF5 file for *timeslice*.

    timeslice == -1  ->  equilibrium.h5
    timeslice >= 0   ->  time_NNN.h5
    """
    if timeslice == -1:
        basename = 'equilibrium.h5'
    else:
        basename = 'time_{:03d}.h5'.format(timeslice)
    return os.path.join(c1_dir, basename)


class sim_data:
    """
    Class that accesses the fusion-io functions in a more pythonic way.

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

    In adios2 mode, sim._all_attrs is an _AcaH5Proxy object that mimics
    an h5py File handle.  Code that accesses sim._all_attrs['scalars/time'],
    sim._all_attrs.attrs['version'], sim._all_attrs.attrs['ntimestep'],
    'scalars' in sim._all_attrs, etc. works transparently for both local
    and remote data without any changes to calling code.

    The proxy searches both the timeslice file (equilibrium.h5) and C1.h5
    for root attributes, so that both 'ntimestep' (from equilibrium.h5) and
    'bzero', 'ntime', etc. (from C1.h5) are accessible via
    sim._all_attrs.attrs[key].  If the embedded skeleton is missing required
    attributes, they are fetched from the remote host automatically.

    --Local Mode (.h5 files)--
    For local M3D-C1 data, sim._all_attrs is a _LocalH5Proxy object that
    wraps both C1.h5 and the current timeslice file (equilibrium.h5 or
    time_NNN.h5).  Root attributes are searched in the timeslice file first,
    then C1.h5, so that attributes like 'ntimestep' are found correctly.
    When set_timeslice() is called, the proxy is updated to reflect the new
    timeslice file.

    Timeslice convention
    --------------------
    time=-1     -> equilibrium.h5
    time=0      -> time_000.h5
    time=N      -> time_NNN.h5   (N > 0)
    time='last' -> highest available timeslice
    """

    def __init__(self, filename='C1.h5', filetype='m3dc1', verbose=False,
                 time=0, fast=False, campaign_simulation='',
                 login='', cache_dir='.cache/aca',
                 remote_python='python3',
                 remote_setup='module load python >/dev/null 2>&1 || true'):
        self._isrc = None
        self.fields = []
        self.seriess = []
        self._aca = None
        self._filetype = filetype
        self._c1_h5 = None       # open h5py.File handle for C1.h5 (local mode)
        self._ts_h5 = None       # open h5py.File handle for timeslice (local mode)

        if filetype == 'm3dc1':
            ifiletype = fio_py.FIO_M3DC1_SOURCE
        elif filetype == 'nimrod':
            if getattr(fio_py, 'FIO_HAS_NIMROD', 0):
                ifiletype = fio_py.FIO_NIMROD_SOURCE
            else:
                raise ValueError(
                    'fusion-io was not compiled with NIMROD support.')
        elif filetype == 'adios2':
            ifiletype = None
        else:
            raise ValueError(
                'Unsupported filetype: {}. '
                'Supported: m3dc1, nimrod, adios2.'.format(filetype))

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

            self._isrc = fio_py.open_numpy_source()
            self._aca.inject_parameters(self._isrc)
            self._aca.inject_scalars(self._isrc)

            ntime_attr = self._aca.get_root_attr('ntime')
            self.ntime = int(ntime_attr) if ntime_attr is not None else 0
            self.available_traces = self._aca.list_group('scalars')
            self.available_traces.extend(['bharmonics', 'keharmonics'])

            # Extract skeleton files for all timeslice HDF5 files and make
            # them available in the CWD so that third-party m3dc1 code that
            # opens files by bare name (e.g. readParameter(fname='equilibrium.h5'))
            # can find them and read attributes like 'ntimestep', 'time', etc.
            self._aca.ensure_all_skeletons_in_cwd()

            # Determine the timeslice file for the requested time.
            if isinstance(time, str) and time == 'last':
                time_int = self.ntime - 1
            else:
                time_int = int(time)
            ts_file, _ = self._aca._ts_file_and_group(time_int)

            # Set _all_attrs to an _AcaH5Proxy so that all existing analysis
            # code that accesses sim._all_attrs like an h5py File handle works
            # transparently for both local and remote data.
            self._all_attrs = _AcaH5Proxy(
                self._aca,
                timeslice_basename=ts_file,
                scalar_basename='C1.h5')
            self._all_attrs_list = self._all_attrs.keys()
            self._all_traces = self._all_attrs['scalars']

        else:
            if not _have_h5py:
                raise ImportError(
                    'h5py is required to read local M3D-C1 files. '
                    'Install with: pip install h5py')

            # Open C1.h5 and keep the handle alive for the lifetime of sim_data.
            self._c1_h5 = h5py.File(self.filename, 'r')

            # Determine the timeslice integer so we can open the timeslice file.
            if isinstance(time, str) and time == 'last':
                _ntime_tmp = self._c1_h5.attrs.get('ntime', 0)
                time_int = int(_ntime_tmp) - 1
            else:
                time_int = int(time)

            # Open the timeslice file (equilibrium.h5 / time_NNN.h5) if it
            # exists alongside C1.h5.
            c1_dir = os.path.dirname(self.filename)
            ts_path = _local_ts_filename(c1_dir, time_int)
            if os.path.isfile(ts_path):
                self._ts_h5 = h5py.File(ts_path, 'r')
            else:
                self._ts_h5 = None

            # Wrap both files in a _LocalH5Proxy so that attribute lookups
            # search the timeslice file first (for 'ntimestep', etc.) and
            # then C1.h5 (for 'bzero', 'ntime', etc.).
            self._all_attrs = _LocalH5Proxy(self._c1_h5, self._ts_h5)
            self._all_attrs_list = self._all_attrs.keys()
            self._all_traces     = self._c1_h5['scalars']
            self.available_traces = list(self._all_traces.keys())
            self.available_traces.extend(['bharmonics', 'keharmonics'])
            self.ntime = self._c1_h5.attrs['ntime']
            self._isrc = fio_py.open_source(ifiletype, self.filename)

        self.available_diagnostics = ['slice times', 'timings', 'iterations']
        self.hint = fio_py.allocate_hint(self._isrc)
        fio_py.get_options(self._isrc)

        if filetype == 'adios2':
            self._aca.inject_mesh(self._isrc, time_int)

        self.set_timeslice(time)
        if fast:
            self._imag = fio_py.get_field(self._isrc, fio_py.FIO_MAGNETIC_FIELD)
            self.fields.append(self._imag)
        self._cs = fio_py.get_int_parameter(self._isrc, fio_py.FIO_GEOMETRY)
        self._iavailable_fields = fio_py.get_available_fields(self._isrc)
        self._available_fields = {}
        for nr in self._iavailable_fields:
            name = _safe_get_field_name(nr)
            if name is not None:
                self._available_fields[name] = nr
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
        # Close local HDF5 file handles.
        if self._ts_h5 is not None:
            try:
                self._ts_h5.close()
            except Exception:
                pass
        if self._c1_h5 is not None:
            try:
                self._c1_h5.close()
            except Exception:
                pass

    def set_timeslice(self, time):
        if isinstance(time, str):
            if time == 'last':
                time = self.ntime - 1
                print('last time slice = ' + str(time))
        time = int(time)
        self.timeslice = time
        fio_py.set_int_option(fio_py.FIO_TIMESLICE, time)

        # For local mode, update the timeslice file handle in the proxy so
        # that _all_attrs.attrs always reflects the current timeslice file.
        if self._filetype != 'adios2' and isinstance(
                self._all_attrs, _LocalH5Proxy):
            c1_dir = os.path.dirname(self.filename)
            ts_path = _local_ts_filename(c1_dir, time)
            # Close the old timeslice handle if it is different from the new one.
            new_ts_h5 = None
            if os.path.isfile(ts_path):
                # Avoid reopening the same file if it hasn't changed.
                current_path = (
                    self._ts_h5.filename if self._ts_h5 is not None else None)
                if current_path != ts_path:
                    if self._ts_h5 is not None:
                        try:
                            self._ts_h5.close()
                        except Exception:
                            pass
                    new_ts_h5 = h5py.File(ts_path, 'r')
                    self._ts_h5 = new_ts_h5
                else:
                    new_ts_h5 = self._ts_h5
            else:
                if self._ts_h5 is not None:
                    try:
                        self._ts_h5.close()
                    except Exception:
                        pass
                    self._ts_h5 = None
            self._all_attrs.update_timeslice(new_ts_h5)

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
        if self._filetype == 'adios2':
            return self._aca.read(path)
        else:
            return np.asarray(self._all_attrs[path])

    def _list_group(self, group):
        if self._filetype == 'adios2':
            return self._aca.list_group(group)
        else:
            return list(self._all_attrs[group].keys())

    def _group_exists(self, group):
        if self._filetype == 'adios2':
            return self._aca.group_exists(group)
        else:
            return group in self._all_attrs

    def _ensure_fields_injected(self, field_name, species, timeslice):
        if self._filetype != 'adios2':
            return

        if field_name in _FIELD_HDF5_DEPS:
            for dep in _FIELD_HDF5_DEPS[field_name]:
                self._aca.inject_field(self._isrc, timeslice, dep)
        else:
            if field_name == 'density':
                hdf5_name = 'ne' if species == 'electron' else 'den'
            elif field_name == 'temperature':
                hdf5_name = 'te' if species == 'electron' else 'ti'
            elif field_name == 'psi':
                hdf5_name = 'psi'
            else:
                hdf5_name = _SIMPLE_FIELD_HDF5_NAME.get(field_name, field_name)
            self._aca.inject_field(self._isrc, timeslice, hdf5_name)

    class field:
        """
        Field class: the init sets up the data access, and its bound methods
        return to you what you want. ex:
        myfield = mysim.get_field('...')
        Allowed strings are:
        'j'   - current density  (C++ loads: psi, I, f/fp)
        'ni'  - ion density      (C++ loads: den)
        'ne'  - electron density (C++ loads: ne)
        'v'   - fluid velocity   (C++ loads: phi, V, chi)
        'B'   - magnetic field   (C++ loads: psi, I, f/fp)
        'p'   - total pressure   (C++ loads: P)
        'pi'  - ion pressure     (C++ loads: P, Pe; computes P-Pe)
        'pe'  - electron pressure(C++ loads: Pe)
        'ti'  - ion temperature  (C++ loads: ti)
        'te'  - electron temperature (C++ loads: te)
        'A'   - vector potential (C++ loads: psi, f)
        'E'   - electric field   (C++ loads: E_R, E_PHI, E_Z)
        'psi' - poloidal flux    (C++ loads: psi)
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

            if sim_data._filetype == 'adios2':
                sim_data._ensure_fields_injected(
                    field, species, sim_data.timeslice)

            itype = None
            if sim_data._filetype == 'adios2':
                itype = _get_fio_field_type(field)

            if itype is None:
                if field in sim_data._available_fields:
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
            Evaluates the required field, returns a one-tuple (scalar),
            three-tuple (vector) or nine-tuple (tensor) at location x.

            **x**
                Three-tuple (R, phi, Z).
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
            Evaluates the derivative of the field at location x.

            **x**
                Three-tuple (R, phi, Z).
            """
            if self.ftype == 'vector':
                try:
                    return fio_py.eval_vector_field_deriv(
                        self._ifield, x, self.sim_data.hint)
                except Exception as e:
                    print('WARNING: eval_vector_field_deriv failed at {}: {}'.format(x, e))
                    return (None,) * 9
            elif self.ftype == 'scalar':
                try:
                    return (fio_py.eval_scalar_field_deriv(
                        self._ifield, x, self.sim_data.hint),)
                except Exception as e:
                    print('WARNING: eval_scalar_field_deriv failed at {}: {}'.format(x, e))
                    return (None, None, None)
            else:
                print('ftype not recognized!')

    class mesh:
        """
        Mesh class: reads the mesh from the timeslice file.
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
        Signal class: for diagnostic signals. Not available in campaign mode.
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
        time_trace class: reads scalar time traces.
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

        def __pow__(self, other):
            return sim_data.time_trace(self.values ** other, time=self.time)

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

            if diagnostic == 'slice times':
                self.diagnostic = np.arange(sim_data.ntime)
                self.x_axis = np.zeros(sim_data.ntime)
                for t in self.diagnostic:
                    self.x_axis[t] = sim_data._get_attr(
                        'time_{:03d}'.format(t), 'time')

            if diagnostic == 'timings':
                self.diagnostic = {}
                for timing in sim_data._list_group('timings'):
                    self.diagnostic[timing] = sim_data._get_dataset(
                        'timings/' + timing)
                self.x_axis = np.arange(
                    self.diagnostic['t_onestep'].shape[0]) + 1

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
