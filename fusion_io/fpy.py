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
]


def _is_real_error(stderr_text):
    """
    Return True only if *stderr_text* contains lines that are not just
    Lmod/module-system noise.
    """
    for line in stderr_text.splitlines():
        if not any(pat.search(line) for pat in _STDERR_NOISE_RES):
            return True
    return False


def _safe_get_field_name(nr):
    """
    Wrapper around fio_py.get_field_name that returns None instead of
    raising/crashing when the field integer ID is not recognised by the
    current source (e.g. a numpy source that only knows a subset of IDs).
    """
    try:
        return fio_py.get_field_name(nr)
    except Exception:
        return None


def _split_ts_prefix(path):
    """
    If *path* starts with 'equilibrium/' or 'time_NNN/', return
    (ts_basename, remainder) where ts_basename is the corresponding HDF5
    filename ('equilibrium.h5' or 'time_NNN.h5') and remainder is the rest
    of the path with the prefix stripped.

    If *path* has no such prefix, return (None, path).

    Examples
    --------
    'equilibrium/mesh/elements' -> ('equilibrium.h5', 'mesh/elements')
    'time_002/fields/psi'       -> ('time_002.h5',    'fields/psi')
    'scalars/time'              -> (None,              'scalars/time')
    """
    m = _TS_PREFIX_RE.match(path)
    if m is None:
        return None, path
    prefix = m.group(1)          # 'equilibrium' or 'time_002'
    remainder = path[m.end():]   # everything after the prefix (and slash)
    ts_basename = prefix + '.h5'
    return ts_basename, remainder


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
        """Return extra SSH options to reuse the master connection."""
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

    Actual .aca SQLite schema
    -------------------------
    info      (id TEXT, name TEXT, version TEXT, modtime INT)
    host      (hostname TEXT PK, longhostname TEXT, modtime INT, deltime INT,
               default_protocol TEXT)
               -- replica.hostid references host.rowid  (implicit integer rowid)
    directory (hostid INT, name TEXT, modtime INT, deltime INT,
               PK (hostid, name))
               -- replica.dirid references directory.rowid
    dataset   (name TEXT PK, uuid TEXT, modtime INT, deltime INT,
               fileformat TEXT, tsid INT, tsorder INT)
               -- replica.datasetid references dataset.rowid
    replica   (datasetid INT, hostid INT, dirid INT, archiveid INT,
               name TEXT, modtime INT, deltime INT, keyid INT, size INT,
               PK (datasetid, hostid, dirid, archiveid, name))
               -- repfiles.replicaid references replica.rowid
    file      (fileid INTEGER PK, name TEXT, compression INT, lenorig INT,
               lencompressed INT, modtime INT, checksum TEXT, data BLOB)
    repfiles  (replicaid INT, fileid INT, PK (replicaid, fileid))

    Timeslice / filename convention
    --------------------------------
      timeslice == -1  ->  equilibrium.h5   (data at file root)
      timeslice ==  0  ->  time_000.h5      (data at file root)
      timeslice ==  N  ->  time_NNN.h5      (data at file root)
      'last'           ->  resolved upstream to ntime-1

    Path routing
    ------------
    Paths that start with 'equilibrium/' or 'time_NNN/' are automatically
    routed to the corresponding timeslice file with the prefix stripped, so
    callers (including sim_data.mesh.read_mesh) can use the same
    'group/subpath' convention regardless of whether data lives in C1.h5 or
    in per-timeslice files.
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

        # 'id' is TEXT in this schema.
        row = cur.execute(
            'SELECT id, name, version FROM info').fetchone()
        if row is None:
            raise RuntimeError(
                'Campaign archive info table is empty: {}'.format(aca_filename))
        self._archive_id, self._archive_name, self._version = row

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
        self._fields_cache = {}       # basename -> set of field names in 'fields/'
        self._mesh_injected = False

        self._mux = _SshMux(login)
        self._mux.__enter__()

    # ------------------------------------------------------------------
    # SQLite helpers
    # ------------------------------------------------------------------

    def _load_datasets(self):
        """
        Load all dataset/replica/file metadata from the .aca SQLite database.

        Key points about the schema:
        - host, directory, dataset, and replica all use SQLite's implicit
          integer rowid as the foreign-key target, even though their declared
          PRIMARY KEY columns are text/composite.  We must SELECT rowid
          explicitly to get those integer ids.
        - repfiles.replicaid references replica.rowid.
        - replica.datasetid references dataset.rowid.
        - replica.hostid    references host.rowid.
        - replica.dirid     references directory.rowid.
        """
        cur = self._con.cursor()

        # host: rowid -> (hostname, longhostname)
        hosts = {
            r[0]: (r[1], r[2])
            for r in cur.execute(
                'SELECT rowid, hostname, longhostname FROM host')}

        # directory: rowid -> (hostid, name)
        dirs = {
            r[0]: (r[1], r[2])
            for r in cur.execute(
                'SELECT rowid, hostid, name FROM directory')}

        # dataset: rowid -> dict  (name is the TEXT primary key)
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

        # replica: rowid -> dataset entry
        # replica.datasetid -> dataset.rowid
        # replica.hostid    -> host.rowid
        # replica.dirid     -> directory.rowid
        replica_rowids = {}   # replica.rowid -> dataset dict (for repfiles)
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

        # repfiles: replicaid -> fileid  (replicaid = replica.rowid)
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

    def _skeleton(self, basename=None):
        basename = basename or self._default_file
        ds = self._get_dataset(basename)
        if ds['name'] not in self._skeletons:
            cur = self._con.cursor()
            written = []
            for ef in ds['embedded_files']:
                fid = ef['id']
                row = cur.execute(
                    'SELECT data, compression FROM file '
                    'WHERE fileid = ?', (fid,)).fetchone()
                if row is None or row[0] is None:
                    continue
                blob, compression = row
                data = zlib.decompress(blob) if compression else bytes(blob)
                dest = os.path.join(
                    self._cache_dir, 'skeletons', ds['name'],
                    os.path.basename(ef['name']))
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with open(dest, 'wb') as fh:
                    fh.write(data)
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
        """
        Run *script* on the remote host via SSH, passing *argv* as arguments.

        Prints a connection message when verbose=True.
        Lmod/module-system noise in stderr is ignored; only genuine errors
        (Python tracebacks, h5py errors, etc.) cause a RuntimeError.
        """
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
        """
        Resolve a logical HDF5 path to (effective_basename, hdf5_path).

        If *path* starts with 'equilibrium/' or 'time_NNN/', the prefix is
        stripped and the corresponding timeslice file is used as the basename,
        regardless of what *basename* was passed in.

        If the resolved timeslice file is not registered in the archive, fall
        back to the default file and keep the full path (for C1.h5 layouts
        where timeslice data is stored under group names).
        """
        ts_basename, remainder = _split_ts_prefix(path.lstrip('/'))
        if ts_basename is not None:
            # Check whether the timeslice file is actually registered.
            if ts_basename in self._files:
                return ts_basename, remainder or '.'
            # Fall back: use whatever basename was requested (or default).
            return basename or self._default_file, path.lstrip('/')
        return basename or self._default_file, path.lstrip('/')

    # ------------------------------------------------------------------
    # Field name discovery
    # ------------------------------------------------------------------

    def list_fields(self, basename):
        """
        Return the set of dataset names present under the 'fields/' group in
        the given timeslice file.  Results are cached in memory.

        This is used by inject_field to discover the actual name of a field
        when the canonical name (e.g. 'A') is not found in the file.
        """
        if basename in self._fields_cache:
            return self._fields_cache[basename]

        names = set()

        # Try the embedded skeleton first (no SSH needed).
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

        # Fall back to a remote listing.
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

    def resolve_field_name(self, hdf5_name, basename, fallbacks=()):
        """
        Return the actual dataset name for a field inside the 'fields/' group
        of *basename*, trying *hdf5_name* first and then each name in
        *fallbacks*.  Raises KeyError if none are found.

        This handles cases where the canonical name used in fpy (e.g. 'A')
        differs from the actual name in the HDF5 file (e.g. 'Aphi').
        """
        available = self.list_fields(basename)
        for candidate in (hdf5_name,) + tuple(fallbacks):
            if candidate in available:
                if candidate != hdf5_name:
                    self._log(
                        "field '{}' not found in {}; using '{}' instead".format(
                            hdf5_name, basename, candidate))
                return candidate
        raise KeyError(
            "Field '{}' (and fallbacks {}) not found in 'fields/' of {}. "
            "Available: {}".format(
                hdf5_name, fallbacks, basename, sorted(available)))

    # ------------------------------------------------------------------
    # Array reading
    # ------------------------------------------------------------------

    def read(self, variable, basename=None, slc=None):
        """
        Return one HDF5 variable as a numpy array.

        If *variable* starts with 'equilibrium/' or 'time_NNN/', the path is
        automatically routed to the correct timeslice file with the prefix
        stripped.

        Resolution order:
        1. Embedded value in the skeleton (small datasets, full reads only).
        2. Per-variable disk cache from a previous streamed read.
        3. Stream from the remote host via SSH + h5py.
        """
        eff_basename, variable = self._route(variable, basename)
        ds = self._get_dataset(eff_basename)

        # 1. Embedded skeleton (small datasets only)
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

    def read_group(self, group, basename=None):
        """
        Fetch all datasets under *group* in a single SSH call and return them
        as a dict {name: ndarray}.  Results are cached as a .npz file.
        """
        eff_basename, group = self._route(group, basename)
        ds = self._get_dataset(eff_basename)

        cache_file = self._group_cache_path(ds, group)
        if os.path.isfile(cache_file):
            self._log('group {}: from cache {}'.format(group, cache_file))
            npz = np.load(cache_file, allow_pickle=False)
            return dict(npz)

        # Try embedded skeleton for small groups.
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

        # Stream the whole group in one SSH call.
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

    def _fetch_root_attrs(self, basename=None):
        """Fetch ALL root attributes in a single call and cache them."""
        import json
        key = basename or self._default_file
        if key not in self._root_attrs_cache:
            if _have_h5py:
                try:
                    with h5py.File(self._skeleton(basename), 'r') as f:
                        def plain(v):
                            try:
                                return (v.item()
                                        if getattr(v, 'shape', None) == ()
                                        else v.tolist())
                            except AttributeError:
                                return (v.decode()
                                        if isinstance(v, bytes) else v)
                        self._root_attrs_cache[key] = {
                            k: plain(v) for k, v in f.attrs.items()}
                        self._log(
                            'root attrs for {}: from skeleton'.format(key))
                        return self._root_attrs_cache[key]
                except Exception:
                    pass

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
        return self._fetch_root_attrs(basename).get(name)

    def get_group_attr(self, group, name, basename=None):
        """
        Read a scalar attribute from a group inside the file.

        If *group* starts with 'equilibrium/' or 'time_NNN/', the path is
        routed to the correct timeslice file automatically.
        """
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
        """Return the names of direct children of a group in the skeleton."""
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
        """Return True if the group exists in the skeleton."""
        eff_basename, eff_group = self._route(group, basename)
        if _have_h5py:
            try:
                with h5py.File(self._skeleton(eff_basename), 'r') as f:
                    return eff_group in f
            except Exception:
                pass
        return False

    # ------------------------------------------------------------------
    # Injection helpers
    # ------------------------------------------------------------------

    def _ts_file_and_group(self, timeslice):
        """
        Return (basename, '') for a given integer timeslice.

          timeslice == -1  ->  equilibrium.h5
          timeslice ==  0  ->  time_000.h5
          timeslice ==  N  ->  time_NNN.h5

        Data always lives at the file root inside timeslice files.
        Falls back to the last known timeslice file if the expected file is
        not registered in the archive.
        """
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
        """
        Read physical normalization parameters from the C1.h5 root attributes
        and inject them into the fio_py numpy source.
        All attributes are fetched in a single call to _fetch_root_attrs().
        """
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
        """
        Stream the mesh arrays from the appropriate timeslice HDF5 file and
        inject them into the fio_py numpy source.
        Only done once; subsequent calls are no-ops.
        """
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
        """
        Stream the field coefficient array for one field at one timeslice
        and inject it into the fio_py numpy source.

        The actual dataset name in the HDF5 file is resolved via
        resolve_field_name(), which checks the 'fields/' group of the
        timeslice file and tries known fallback names when the canonical
        name is not present.  This handles cases like 'A' vs 'Aphi'.
        """
        ts_file, _ = self._ts_file_and_group(timeslice)

        # Known fallback names for fields whose canonical fpy name differs
        # from the actual HDF5 dataset name in some M3D-C1 output versions.
        _fallbacks = {
            'A'  : ('Aphi', 'Avec', 'vecpot'),
            'J'  : ('Jtor', 'jtor', 'j'),
            'den': ('ni', 'n_i', 'density'),
            'Pe' : ('pe', 'P_e'),
            'P'  : ('p', 'ptot'),
        }
        fallbacks = _fallbacks.get(hdf5_name, ())

        try:
            actual_name = self.resolve_field_name(hdf5_name, ts_file, fallbacks)
        except KeyError as e:
            raise RuntimeError(str(e))

        real_path = 'fields/{}'.format(actual_name)
        self._log('Injecting field {} (as {}) timeslice {} from {}'.format(
            hdf5_name, actual_name, timeslice, ts_file))

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

    def inject_scalars(self, isrc):
        """
        Stream all scalar time series from C1.h5 in a single SSH call and
        inject them into the fio_py numpy source.
        """
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

    Timeslice convention
    --------------------
    time=-1     -> equilibrium.h5
    time=0      -> time_000.h5
    time=N      -> time_NNN.h5   (N > 0)
    time='last' -> highest available timeslice

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
            if true, gives extra info on file read. In campaign archive mode,
            also prints a message each time a connection to the remote host
            is established.

        **time**
            sets the time variable during initialization.
            -1        -> equilibrium.h5
             0        -> time_000.h5
             N        -> time_NNN.h5
            'last'    -> highest available timeslice

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
            Shell command run before the remote Python script. Module-load
            errors in this command are automatically ignored so that HPC
            systems that print Lmod warnings do not cause failures.
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

        if filetype == 'adios2':
            if isinstance(time, str) and time == 'last':
                time_int = self.ntime - 1
            else:
                time_int = int(time)
            self._aca.inject_mesh(self._isrc, time_int)

        self.set_timeslice(time)
        if fast:
            self._imag = fio_py.get_field(self._isrc, fio_py.FIO_MAGNETIC_FIELD)
            self.fields.append(self._imag)
        self._cs = fio_py.get_int_parameter(self._isrc, fio_py.FIO_GEOMETRY)
        self._iavailable_fields = fio_py.get_available_fields(self._isrc)
        # Build _available_fields dict, skipping any field IDs whose name
        # cannot be resolved (e.g. IDs not supported by the numpy source).
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
        """
        Read a scalar attribute, dispatching to h5py or _AcaBackend.

        In campaign archive mode, paths starting with 'equilibrium' or
        'time_NNN' are automatically routed to the correct timeslice file
        via _AcaBackend.get_group_attr (which calls _route internally).
        """
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
        """
        Read a dataset as a numpy array.

        In campaign archive mode, paths starting with 'equilibrium/' or
        'time_NNN/' are automatically routed to the correct timeslice file
        via _AcaBackend.read (which calls _route internally).
        """
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
        Ensure field coefficients for the given timeslice have been injected.
        For local files this is a no-op.
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
                # _get_dataset and _get_attr route 'equilibrium/...' and
                # 'time_NNN/...' paths to the correct timeslice file.
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
