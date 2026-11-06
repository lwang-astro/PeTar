# read snapshot and obtain multiple systems
import collections
from scipy import spatial as sp
from sdar.base import *
from sdar.functions import *
from sdar.ar import SDARInterruptBinary
from sdar.ar import SDARData as sdar_SDARData
from sdar.particle import SimpleParticle, ParticleGroup
from sdar.particle import Binary as sdar_Binary
from sdar.hermite import HermiteData as hermite_HermiteData
import sdar.group as hermite_group
from .bse import *
from .dsm import *
from .astro_units import G_ASTRO
from sdar.base import enablePathInit, AmbiguousSchemaError, detectAsciiFile

# alias kept for existing tools; single source: tools/analysis/astro_units.py
# (generated together with src/astro_units.hpp by tools/generate_astro_units.py)
G_MSUN_PC_MYR = G_ASTRO # pc^3/(Msun*Myr^2)
G_HENON=1 # Henon unit
HEADER_OFFSET=24 # header offset in bytes for snapshots with the BINARY format
HEADER_OFFSET_F128=32 # header offset in bytes for snapshots with the BINARY format
HEADER_OFFSET_WITH_CM=72 # header offset with center-of-the-mass data in bytes for snapshots with the BINARY format
HEADER_OFFSET_WITH_CM_F128=128 # header offset with center-of-the-mass data in bytes for snapshots with the BINARY format

class PeTarDataHeader():
    """ Petar snapshot data header
    members:
        file_id: int 
           file id
        n: int 
           number of particles
        time: float 
           time of snapshot
        *pos_offset: list of float (length of 3)
           position offset of particle system
        *vel_offset: list of float (length of 3)
           velocity offset of particle system
  
        pos_offset and vel_offset only exist when keyword argument 'external_mode' is not none
    """

    # binary header variants: (header bytes, c.m.-offset entries present, float type)
    BINARY_VARIANTS = ((24, False, np.float64), (72, True, np.float64),
                       (32, False, np.float128), (128, True, np.float128))
    # ASCII header line item counts: (file_id N time) or (+ 6 c.m. offsets)
    ASCII_HEADER_ITEMS = (3, 9)

    @classmethod
    def asciiHeaderInfo(cls, path):
        """ (external_present, first_data_row_ncols) of an ASCII snapshot, or
        None when the first line is not a snapshot header (file_id N time
        [cm offsets]) - e.g. headerless event/stream tables. """
        with open(str(path)) as fp:
            header = fp.readline().split()
            first_row = fp.readline().split()
        if len(header) not in cls.ASCII_HEADER_ITEMS:
            return None
        try:
            int(header[0]); int(header[1]); float(header[2])
        except (ValueError, IndexError):
            return None
        if not first_row:
            raise ValueError('%s: snapshot has no particle rows; pass explicit reader kwargs' % path)
        return (len(header) == 9, len(first_row))

    @classmethod
    def binaryHeaderVariants(cls, path):
        """ (n_glb, [(hsize, cm_present, float_type), ...]) of a binary
        snapshot: reads the leading (file_id, N) ints (fingerprint sanity)
        and returns the header variants that remain size-plausible. None
        when the file carries no snapshot header (derived/stream files and
        non-snapshot data start with data, whose leading bytes do not pass
        the fingerprint). """
        import os
        path = str(path)
        size = os.path.getsize(path)
        if size < 16:
            return None
        head = np.fromfile(path, dtype=np.int64, count=2)
        file_id, n_glb = int(head[0]), int(head[1])
        if not (0 <= file_id <= 10**8) or (n_glb <= 0):
            return None
        variants = [(h, cm, ft) for (h, cm, ft) in cls.BINARY_VARIANTS if size > h]
        return (n_glb, variants) if variants else None

    def __init__(self, _filename=None, **kwargs):
        """ Initial data header
        
        Parameters:
        -----------
        _filename: string
            PeTar snapshot file name to read the header, if not provide, all members are initialized to zero (None)
        kwargs: dict
            Keyword arguments:
            float_type: type (np.float64)
                floating point data type
            snapshot_format: string (binary)
                Data format of snapshot files: binary or ascii
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
                If not none, this option indicates the pos_offset and vel_offset exists 
        """
        self.file_id = int(0)
        self.n = int(0)
        self.time = 0.0
        self.pos_offset=[0.0,0.0,0.0]
        self.vel_offset=[0.0,0.0,0.0]
        self.offset_flag=False
        
        if (_filename!=None): self.read(_filename,**kwargs)

    def read(self, _filename, **kwargs):
        """ Read snapshot file to obtain the header information

        Parameters:
        -----------
        _filename: string
            PeTar snapshot file name to read the header
        kwargs: dict
            Keyword arguments:
            float_type: type (np.float64)
                floating point data type
            snapshot_format: string (binary)
                Data format of snapshot files: binary or ascii
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
                If not none, this option indicates the pos_offset and vel_offset exists 
        """
        snapshot_format='binary'
        if ('snapshot_format' in kwargs.keys()): snapshot_format=kwargs['snapshot_format']
        if ('external_mode' in kwargs.keys()):
            if (kwargs['external_mode']!='none'): self.offset_flag=True
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        if (snapshot_format=='ascii'):
            fp = open(_filename, 'r')
            header=fp.readline()
            header_items=header.split()
            if (self.offset_flag):
                if (len(header_items)!=9):
                    raise ValueError('Snapshot header item number mismatch! Need 9 (file_id, N, time, xcm, ycm, zcm, vxcm, vycm, vzcm), got %d. Make sure the external_mode keyword set correctly.' % len(header_items))

                file_id, n_glb, t, x, y, z, vx, vy, vz = header_items
                fp.close()

                self.file_id = int(file_id)
                self.n = int(n_glb)
                self.time = float_type(t)
                self.pos_offset = [float_type(x),float_type(y),float_type(z)]
                self.vel_offset = [float_type(vx),float_type(vy),float_type(vz)]
            else:
                if (len(header_items)!=3):
                    raise ValueError('Snapshot header item number mismatch! Need 3 (file_id, N, time), got %d. Make sure the external_mode keyword set correctly.' % len(header_items))

                file_id, n_glb, t = header_items
                fp.close()

                self.file_id = int(file_id)
                self.n = int(n_glb)
                self.time = float_type(t)

        else:
            if (self.offset_flag):
                fp = np.fromfile(_filename, dtype=np.dtype([('file_id',np.int64),('n_glb',np.int64),('time',float_type),('x',float_type),('y',float_type),('z',float_type),('vx',float_type),('vy',float_type),('vz',float_type)]),count=1)
                self.file_id = fp['file_id'][0]
                self.n = fp['n_glb'][0]
                self.time = fp['time'][0]
                self.pos_offset = np.array([fp['x'][0], fp['y'][0], fp['z'][0]])
                self.vel_offset = np.array([fp['vx'][0], fp['vy'][0], fp['vz'][0]])
            else:
                fp = np.fromfile(_filename, dtype=np.dtype([('file_id',np.int64),('n_glb',np.int64),('time',float_type)]),count=1)
                self.file_id = fp['file_id'][0]
                self.n = fp['n_glb'][0]
                self.time = fp['time'][0]

    def savetxt(self, fname, **kwargs):
        """ Save class member data to a file
        Use the getherDataToArray and then numpy.savetxt

        Parameters
        ----------
        fname: string of filename or file handler
        kwargs: dict
            keyword arguments for numpy.savetxt
        """
        offset_flag=False
        close_flag = False
        if (type(fname)==str):
            f = open(fname, 'w')
            close_flag = True
        elif (hasattr(fname, 'write')):
            f = fname
        if (self.offset_flag):
            f.write("%d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n" % 
                    (self.file_id, self.n, self.time, 
                     *self.pos_offset, *self.vel_offset))
        else:
            f.write("%d %d %.20g\n" % (self.file_id, self.n, self.time))

        if close_flag:
            f.close()

    def tofile(self, fname):
        """ Write class member data to a file using numpy.save
        Use numpy.save to write data, the dtype is defined by keys (members)

        Parameters
        ----------
        fname: string of filename or file header
        kwargs: dict
            keyword arguments for numpy.save, notice dtype is already defined, do not provide that
        """

        close_flag = False
        if (type(fname)==str):
            f = open(fname, 'wb')
            close_flag = True
        elif (hasattr(fname, 'write')):
            f = fname

        import struct
        if (self.offset_flag):
            header_buffer = struct.pack('qqddddddd', self.file_id, self.n, self.time, *self.pos_offset, *self.vel_offset)
        else:
            header_buffer = struct.pack('qqd', self.file_id, self.n, self.time)
        f.write(header_buffer)

        if close_flag:
            f.close()
        

    def toSkyCoord(self, **kwargs):
        """ generate astropy.coordinates.SkyCoord data in galactocentric frame

        Parameters
        -----------------
        kwargs: dict()
            pos_unit: astropy.units (units.pc)
                 position unit of the particle data
            vel_unit: astropy.units (units.pc/units.Myr)
                 velocity unit of the particle data
            galcen_distance: floating with length units (8.0*units.kpc [Galpy])
                 galactic central distance of the Sun
            z_sun: floating with length units (15.0*units.pc [Galpy])
                 z direction distance of the Sun
            galcen_v_sun: astropy.coordinates.CartesianDifferential ([10.0, 235.0, 7.0]*units.km/units.s [Galpy])
                 velocity of the Sun

        Return
        ----------------
        core_g: astropy.coordinates.SkyCoord
            core c.m. data using SkyCoord
        """
        import astropy 
        from astropy.coordinates import SkyCoord  # High-level coordinates
        from astropy.coordinates import ICRS, Galactic, Galactocentric, FK4, FK5  # Low-level frames
        from astropy.coordinates import Angle, Latitude, Longitude  # Angles
        from astropy.coordinates import CartesianDifferential
        import astropy.units as u

        pos_unit = u.pc
        if ('pos_unit' in kwargs.keys()): pos_unit = kwargs['pos_unit']
        vel_unit = u.pc/u.Myr
        if ('vel_unit' in kwargs.keys()): vel_unit = kwargs['vel_unit']

        parameters={'galcen_distance':8.0*u.kpc, 'z_sun':15.*u.pc, 'galcen_v_sun':CartesianDifferential([10.0,235.,7.]*u.km/u.s)}
        for key in parameters.keys():
            if key in kwargs.keys():
                parameters[key] = kwargs[key]

        sky = SkyCoord(x=self.pos_offset[0]*pos_unit, 
                       y=self.pos_offset[1]*pos_unit, 
                       z=self.pos_offset[2]*pos_unit, 
                       v_x=self.vel_offset[0]*vel_unit,
                       v_y=self.vel_offset[1]*vel_unit,
                       v_z=self.vel_offset[2]*vel_unit,
                       frame='galactocentric', representation_type='cartesian', **parameters)
        return sky
    
#: derived per-snapshot products sharing the particle record layout
_DERIVED_SUFFIXES = ('.single', '.binary', '.triple', '.quadruple')

#: interrupt_mode values sharing one reader schema; 'bse' represents them all
BSE_FAMILY = ('bse', 'mobse', 'bseEmp')


def _isSnapshotFamilyName(name):
    """ True for petar snapshot-family file names (snapshots, derived
    products, run-level streams) - used only for the strict no-silent-
    misparse policy, not for sidecar lookup (that uses the generic walk). """
    import re
    if any(name.endswith(s) for s in _DERIVED_SUFFIXES) or name.endswith('.npy'):
        return True
    if re.search(r'\.\d+$', name):
        return True
    return bool(re.search(r'\.(?:esc_single|esc_binary|lagr|core|status|group(?:\.\d+)?\.n\d+'
                          r'|sse(?:\.(?:type_change|sn_kick))?|bse_status|interrupt'
                          r'|bse(?:\.(?:type_change|sn_kick|gw_kick|dynamic_merge|binary_merge'
                          r'|hyperbolic_tde|binary_tde|tide|gw_tide_merge))?)$', name))


def _refineLabelByRunFiles(path, value):
    """ Values detected structurally can be family representatives (e.g. the
    external-potential variant): probe the solver-written parameter sidecars
    of each run-prefix stem (generic walk) to name the actual variant, else
    keep the detected value """
    import os
    directory = os.path.dirname(str(path)) or '.'
    name = os.path.basename(str(path))
    if name.endswith('.npy'):
        name = name[:-len('.npy')]
    parts = name.split('.')
    for n in range(len(parts), 0, -1):
        stem = '.'.join(parts[:n])
        for variant in ('galpy', 'agama'):
            if os.path.isfile(os.path.join(directory, stem + '.par.' + variant)):
                return variant
    return value


class ParticleRecordLayoutMixin:
    """ Layout machinery for petar readers whose FILE RECORDS embed the
    particle-record layout (the snapshot family: snapshots, .single/.binary/
    .triple/.quadruple, escapers, group files). Their records can only be
    parsed after the particle layout combination is pinned, using the
    snapshot-family evidence chain (sidecar -> sibling snapshot header ->
    divisibility). Readers with their own column layouts (Lagrangian, Core,
    SSE/BSE event tables, Status, SDAR/Hermite data) do NOT use this mixin -
    the generic sdar.base machinery suffices for them.

    The mixin itself contains NO parameter names: everything layout-specific
    is read from class-level declarations -

      LAYOUT_PARAMS          (sdar.base) - enumerable layout parameters,
                              first value = default/most plausible
      HEADER_PINNED_PARAMS   - physical snapshot-header facts that pin a
                              parameter: fact name -> parameter name
                              ('cm': 'external_mode', 'float128': 'float_type')
      LABEL_REFINE_PARAMS    - parameters whose detected value is only a
                              family representative and should be refined by
                              probing the run directory (galpy/agama)
      PINNING_PARAMS         - parameters that pin the main layout; given
                              explicitly they exempt a detection failure from
                              the strict raise policy

    Mixed into every petar reader whose records carry the particle layout
    (particle chain, Binary, escapers, GroupInfo); each class detects with
    ITS OWN record layout and declarations. """

    LAYOUT_PARAMS = {}
    #: header-fact -> layout-parameter bindings (overridden by declaring classes)
    HEADER_PINNED_PARAMS = {}
    #: parameters refined by run-directory probing after structure detection
    LABEL_REFINE_PARAMS = ()
    #: explicitly-given parameters that pin the layout (failure-policy exemption)
    PINNING_PARAMS = ()

    @classmethod
    def _resolveTierMatches(cls, path, matches, facts):
        """ Plausibility-tier resolution over raw structure matches (fewest
        non-default layout parameters wins; same-tier collision or no match
        raises AmbiguousSchemaError; suppressed alternatives warn). """
        if not matches:
            raise AmbiguousSchemaError(
                'No known reader schema matches %s (%s): the output is likely from a PeTar '
                'version or build variant whose layout this reader does not know. Report it '
                'as a version gap instead of working around it.' % (path, facts))
        tier = cls.layoutTier
        min_tier = min(tier(m) for m in matches)
        best = [m for m in matches if tier(m) == min_tier]
        if len(best) > 1:
            cand = ', '.join(str({k: v for k, v in m.items()}) for m in best)
            keys = ', '.join(cls.LAYOUT_PARAMS.keys())
            raise AmbiguousSchemaError(
                'Cannot uniquely determine the reader schema of %s (%s): candidates [%s]. '
                'Pass explicit layout kwargs (%s) or provide/regenerate the .meta sidecar '
                'next to the data.' % (path, facts, cand, keys))
        if min_tier < max(tier(m) for m in matches):
            suppressed = [str({k: v for k, v in m.items()})
                          for m in matches if tier(m) > min_tier]
            import warnings as _w
            _w.warn('Schema detection of %s (%s) resolved to the default-build layout %s; the '
                    'file could also match rarer builds [%s]. If the producing solver used '
                    'those flags, pass explicit kwargs or a .meta sidecar.'
                    % (path, facts, best[0], '; '.join(suppressed)))
        return dict(best[0])

    @classmethod
    def _factCheck(cls, kw, fact, value):
        """ generic header-fact filter: the parameter bound to `fact` (see
        HEADER_PINNED_PARAMS) must be consistent with the measured fact """
        param = cls.HEADER_PINNED_PARAMS.get(fact)
        if param is None:
            return True
        kw_value = kw.get(param)
        if fact == 'cm':
            return ((kw_value is not None) and (kw_value != cls.LAYOUT_PARAMS[param][0])) == value
        if fact == 'float128':
            return (kw_value == np.float128) == value
        return True

    @classmethod
    def _snapshotMatches(cls, path):
        """ Raw (pre-resolution) structure evidence for a snapshot-family
        file: (matches, pinned, facts) where `pinned` maps fact-name ->
        measured value for every header fact the file pins (ASCII header
        item count -> cm; binary header variant -> cm/float128). Derived
        products fall back to their raw sibling; isolated headerless
        products give weak divisibility evidence (unique-match enforced by
        the caller). """
        import os
        path = str(path)
        if detectAsciiFile(path):
            info = PeTarDataHeader.asciiHeaderInfo(path)
            if info is None:
                raise ValueError('%s: first line is not a PeTar snapshot header' % path)
            ext_pinned, ncol = info
            matches = [kw for kw in cls._layoutMatches(ncols=ncol)
                       if cls._factCheck(kw, 'cm', ext_pinned)]
            return matches, {'cm': ext_pinned}, 'ASCII snapshot, %d columns' % ncol
        head = PeTarDataHeader.binaryHeaderVariants(path)
        if head is None:
            name = os.path.basename(path)
            if name.endswith('.npy') and os.path.isfile(path[:-len('.npy')]):
                return cls._snapshotMatches(path[:-len('.npy')])
            derived = next((s for s in _DERIVED_SUFFIXES if name.endswith(s)), None)
            if (derived is not None) and os.path.isfile(path[:-len(derived)]):
                return cls._snapshotMatches(path[:-len(derived)])
            size = os.path.getsize(path)
            if size == 0:
                raise ValueError('%s is empty' % path)
            divisors = set()
            for d in range(1, int(size ** 0.5) + 1):
                if size % d == 0:
                    divisors.add(d)
                    divisors.add(size // d)
            return cls._layoutMatches(record_size=divisors), {}, \
                'headerless file, %d bytes' % size
        n_glb, variants = head
        size = os.path.getsize(path)
        matches = []
        cm_vals, f128_vals = set(), set()
        for hsize, cm_present, ftype in variants:
            body = size - hsize
            if (body <= 0) or (body % n_glb != 0):
                continue
            for kw in cls._layoutMatches(record_size=body // n_glb):
                if cls._factCheck(kw, 'cm', cm_present) and cls._factCheck(kw, 'float128', ftype == np.float128):
                    matches.append(kw)
                    cm_vals.add(cm_present)
                    f128_vals.add(ftype == np.float128)
        pinned = {}
        if len(cm_vals) == 1:
            pinned['cm'] = cm_vals.pop()
        if len(f128_vals) == 1:
            pinned['float128'] = f128_vals.pop()
        return matches, pinned, 'binary snapshot, %d bytes, N = %d' % (size, n_glb)

    @classmethod
    def _detectSnapshotLayout(cls, path):
        """ Resolved layout kwargs for a snapshot-family file: raw structure
        evidence from _snapshotMatches, plausibility-tier resolution, then
        label refinement (e.g. galpy/agama) by run-directory probing. """
        matches, _pinned, facts = cls._snapshotMatches(str(path))
        if facts.startswith('headerless'):
            # weak divisibility evidence: several layouts can divide one
            # file size, so a plausibility pick could misparse silently;
            # the match must be unique
            if len(matches) != 1:
                cand = ', '.join(str({k: v for k, v in m.items()}) for m in matches) or 'none'
                raise AmbiguousSchemaError(
                    'Cannot uniquely determine the reader schema of %s (%s): candidates [%s]. '
                    'Pass explicit kwargs or provide/regenerate the .meta sidecar next to the data.'
                    % (path, facts, cand))
            kw = dict(matches[0])
        else:
            kw = cls._resolveTierMatches(str(path), matches, facts)
        for param in cls.LABEL_REFINE_PARAMS:
            default = cls.LAYOUT_PARAMS.get(param, (None,))[0]
            if (param in kw) and (kw[param] != default):
                kw[param] = _refineLabelByRunFiles(str(path), kw[param])
        return kw

    @classmethod
    def _verifySidecarAgainstStructure(cls, path, kwargs, sidecar_layout=None):
        """ Cross-validate sidecar-provided layout kwargs with the file
        structure, with distinct diagnoses: (a) the SIDECAR-declared layout
        is not among the structure candidates -> the sidecar was likely
        copied from another run (ValueError); (b) the sidecar agrees with
        the structure but explicitly given kwargs contradict both -> the
        error names the offending explicit parameter(s) (e.g. --spin-1d on
        3D-spin data) instead of blaming the sidecar. Generic shape check:
        the layout instantiated with THIS class must match one of the
        measured candidate shapes (ALL raw candidates - aliasing must not
        refute a sidecar that is among them); header-pinned facts are
        checked through HEADER_PINNED_PARAMS. Weak (headerless) evidence
        cannot refute. Skipped when the structure cannot refute. """
        import os
        try:
            matches, pinned, facts = cls._snapshotMatches(str(path))
        except (AmbiguousSchemaError, ValueError):
            return
        if (not matches) or facts.startswith('headerless'):
            return
        ascii_file = detectAsciiFile(str(path))

        def shape_of(layout_kwargs):
            try:
                obj = cls(**dict(layout_kwargs))
            except Exception:
                return None
            if ascii_file:
                return obj.ncols
            return np.dtype(obj.collectDtype()).itemsize

        ref_shapes = set()
        for m in matches:
            s = shape_of({k: v for k, v in m.items() if k in cls.LAYOUT_PARAMS})
            if s is not None:
                ref_shapes.add(s)
        if not ref_shapes:
            return

        sidecar_kwargs = {k: v for k, v in (sidecar_layout or {}).items() if k in cls.LAYOUT_PARAMS}
        merged_kwargs = {k: v for k, v in kwargs.items() if k in cls.LAYOUT_PARAMS}
        if not sidecar_kwargs:
            return
        # kwargs at this point = explicit values with sidecar values merged
        # via setdefault: anything differing from (or absent in) the sidecar
        # came from the caller (CLI flags)
        explicit_conflicts = {
            k: v for k, v in merged_kwargs.items()
            if k not in sidecar_kwargs or sidecar_kwargs[k] != v
        }

        def reject(key, meta_value, struct_value):
            raise ValueError(
                'data.meta disagrees with the file structure of %s (%s: meta=%s, structure=%s). '
                'The sidecar was likely produced by a different run; remove or regenerate it.'
                % (path, key, meta_value, struct_value))

        def reject_explicit(key, explicit_value, struct_value):
            raise ValueError(
                'explicitly given layout parameter(s) disagree with the file structure of %s '
                '(%s: explicit=%s, structure=%s). Check flags like --spin-1d / -i (--interrupt-mode) '
                '/ -t (--external-mode) / --use-mpfrc against the data.'
                % (path, key, explicit_value, struct_value))

        side_shape = shape_of(sidecar_kwargs)
        merged_shape = shape_of(merged_kwargs)
        shape_key = 'layout columns' if ascii_file else 'layout record size'
        if side_shape not in ref_shapes:
            reject(shape_key, side_shape, sorted(ref_shapes))
        if (merged_shape not in ref_shapes) and explicit_conflicts:
            reject_explicit(shape_key + ' ' + ', '.join('%s=%r' % kv for kv in explicit_conflicts.items()), merged_shape, sorted(ref_shapes))

        for fact, value in pinned.items():
            param = cls.HEADER_PINNED_PARAMS.get(fact)
            if param is None or param not in sidecar_kwargs or param not in merged_kwargs:
                continue
            if fact == 'cm':
                side_on = sidecar_kwargs[param] != cls.LAYOUT_PARAMS[param][0]
                merged_on = merged_kwargs[param] != cls.LAYOUT_PARAMS[param][0]
                if side_on != value:
                    reject(param, sidecar_kwargs[param], 'external' if value else cls.LAYOUT_PARAMS[param][0])
                if merged_on != value and param in explicit_conflicts:
                    reject_explicit(param, merged_kwargs[param], 'external' if value else cls.LAYOUT_PARAMS[param][0])
            elif fact == 'float128':
                if (sidecar_kwargs[param] == np.float128) != value:
                    reject(param, sidecar_kwargs[param], np.float128 if value else np.float64)
                if (merged_kwargs[param] == np.float128) != value and param in explicit_conflicts:
                    reject_explicit(param, merged_kwargs[param], np.float128 if value else np.float64)

    @classmethod
    def resolvePathInput(cls, _dat, kwargs):
        """ petar concise-constructor entry: sidecar first (generic walk,
        legacy bridges, values merged explicitly), then snapshot structure
        detection with the strict failure policy; explicit kwargs and the
        ``meta=`` injection always win. """
        import os
        if not isinstance(_dat, (str, bytes, os.PathLike)):
            return _dat, kwargs, None
        path = str(_dat)
        if not os.path.isfile(path):
            return _dat, kwargs, None
        meta_arg = kwargs.pop('meta', 'auto')
        if isinstance(meta_arg, dict):
            for key, value in meta_arg.items():
                kwargs.setdefault(key, value)
            return None, kwargs, path
        if (meta_arg is not False) and (meta_arg != 'none'):
            if (meta_arg is not None) and (meta_arg != 'auto'):
                meta_file = str(meta_arg)
                if not os.path.isfile(meta_file):
                    raise ValueError('meta sidecar not found: %s' % meta_file)
            else:
                meta_file = cls.findLayoutMeta(path)
            if meta_file is not None:
                try:
                    layout = cls.readLayoutMeta(meta_file)
                except ValueError:
                    layout = None
                if not layout:
                    # transitional bridges: sidecars written before the
                    # layout.* key convention carry unprefixed keys
                    layout = _legacyPetarLayoutMeta(meta_file)
                if layout:
                    # petar readers test kwarg truthiness, not presence, so
                    # merging values equal to the defaults is safe and keeps
                    # the modes explicit for CLI auto-fill
                    for name, value in layout.items():
                        kwargs.setdefault(name, value)
                    cls._verifySidecarAgainstStructure(path, kwargs, sidecar_layout=layout)
                    return None, kwargs, path
        try:
            detected = cls._detectSnapshotLayout(path)
        except (AmbiguousSchemaError, ValueError) as error:
            cls._onLayoutDetectionFailure(path, kwargs, error)
            return None, kwargs, path
        for key, value in detected.items():
            kwargs.setdefault(key, value)
        return None, kwargs, path

    @classmethod
    def _onLayoutDetectionFailure(cls, path, kwargs, error):
        """ petar policy: for snapshot-family files a default-kwargs
        fallback can misparse SILENTLY when the byte size happens to align,
        so detection failure raises unless the caller pinned the main
        layout (all PINNING_PARAMS given explicitly); other file types
        (e.g. hermite .log dumps) warn. """
        import os
        pinning_given = all(p in kwargs for p in cls.PINNING_PARAMS)
        if (not pinning_given) and _isSnapshotFamilyName(os.path.basename(path)):
            keys = ', '.join(list(cls.LAYOUT_PARAMS.keys()) + ['...'])
            raise AmbiguousSchemaError(
                'Cannot auto-detect the reader schema for the snapshot-family file %s (%s). '
                'Pass the layout kwargs (%s) explicitly, or read it next to its run outputs '
                '(the <prefix>.meta sidecar or the raw snapshots make the schema recoverable); '
                'a fallback to default kwargs could misparse the file silently.'
                % (path, error, keys)) from error
        import warnings as _w
        _w.warn('Cannot auto-detect the reader schema for %s (%s); constructing with the '
                'given/default keyword arguments. Pass layout kwargs explicitly if the read fails.'
                % (path, error))
        return kwargs

    def readAutoFile(self, path):
        """ petar format dispatch: .npy -> load(); ASCII -> loadtxt with the
        snapshot header line skipped when present; binary -> fromfile with
        the snapshot header offset when the leading (file_id, N) fingerprint
        and the record size identify a headered snapshot, else offset 0. """
        import os
        path = str(path)
        if path.endswith('.npy'):
            self.load(path)
            return self
        if detectAsciiFile(path):
            skip = 1 if PeTarDataHeader.asciiHeaderInfo(path) is not None else 0
            self.loadtxt(path, skiprows=skip)
            return self
        head = PeTarDataHeader.binaryHeaderVariants(path)
        if head is not None:
            n_glb, variants = head
            size = os.path.getsize(path)
            rec = np.dtype(self.collectDtype()).itemsize
            for hsize, _cm, _ft in variants:
                body = size - hsize
                if (body > 0) and (body % n_glb == 0) and (body // n_glb == rec):
                    self.fromfile(path, offset=hsize)
                    return self
        self.fromfile(path)
        return self


def _legacyPetarLayoutMeta(meta_file):
    """ Transitional bridge: layout kwargs from a pre-release petar-meta/1
    sidecar (unprefixed keys, '0'/'1' booleans). Removed once all run
    directories carry output-meta/1 sidecars with layout.* keys. """
    meta = readMetaData(meta_file)
    layout = {}
    if meta.get('interrupt_mode'):
        layout['interrupt_mode'] = meta['interrupt_mode']
    if meta.get('external_mode'):
        layout['external_mode'] = meta['external_mode']
    for key in ('use_mpfrc', 'collect_sp_acc', 'spin_3d'):
        if key in meta:
            layout[key] = meta[key] in ('1', 'true', 'True')
    if meta.get('float_type') == 'float128':
        layout['float_type'] = np.float128
    return layout


class BaseParticle(ParticleRecordLayoutMixin, SimpleParticle):
    """ Base particle type of PeTar
        The members include simple particle information, binary status and stellar evolution data

    keys: (class members)
        Members inherited from SimpleParticle: mass (1D), pos (2D,3), *pos_high (2D,3) vel (2D,3) 
            see help(petar.SimpleParticle)
        binary_state (1D): binary interruption state
        if (keyword argument 'interrupt_mode' == 'merger', 'bse', 'bseEmp', 'mobse', 'dsm', 'sevn'):
            radius:        (1D): radius for merger checker
            dm:            (1D): mass loss
            time_record    (1D): last time of interruption check
            time_interrupt (1D): next interruption time
        if (keyword argument 'interrupt_mode' == 'bse', 'bseEmp', 'mobse'):
            star  (SSEStarParameter): BSE based stellar evolution parameters
        if (keyword argument 'interrupt_mode' == 'dsm'):
            star  (DSMStarParameter): Disk star merger parameters
    """

    # Layout-parameter declaration for the generic schema detection in
    # DictNpArrayMix (sdar.base): the constructor keyword arguments that change
    # the column layout, with all allowed values (first = default/most
    # plausible). Combinations producing identical columns are collapsed
    # automatically, so parameters that only matter for some modes (e.g.
    # spin_3d for the bse family) can be declared unconditionally.
    LAYOUT_PARAMS = {
        'interrupt_mode': ('none', 'merger', 'bse', 'dsm', 'sevn'),
        'use_mpfrc': (False, True),
        'spin_3d': (True, False),
        'float_type': (np.float64, np.float128),
    }

    # snapshot-header facts that pin a layout parameter (see ParticleRecordLayoutMixin)
    HEADER_PINNED_PARAMS = {'cm': 'external_mode', 'float128': 'float_type'}
    #: structurally-detected family representatives refined by run-dir probing
    LABEL_REFINE_PARAMS = ('external_mode',)
    #: explicitly-given parameters that pin the layout (failure exemption)
    PINNING_PARAMS = ('interrupt_mode', 'external_mode')

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
               This option indicates whether columns of stellar evolution exist
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
            spin_3d: bool (True)
                when interrupt_mode is bse/bseEmp/mobse: width of the star.spin column, 3 (current) or 1;
                use False for outputs produced before Dec 2024 (1D spin)
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys_bstat = [['binary_state',np.int64]]
        keys_se  = [['radius',float_type],['dm',float_type],['time_record',float_type],['time_interrupt',float_type]]    
        
        keys = keys_bstat
        if ('interrupt_mode' in kwargs.keys()):
            if (kwargs['interrupt_mode'] in ('base', 'merger')):
                keys = keys_bstat+keys_se
            elif ('bse' in kwargs['interrupt_mode']):
                keys = keys_bstat+keys_se+[['star',SSEStarParameter]]
            elif (kwargs['interrupt_mode']=='dsm'):
                keys = keys_bstat+keys_se+[['star',DSMStarParameter]]
            elif ('sevn' in kwargs['interrupt_mode']):
                keys = keys_bstat+keys_se+[['star', SEVNStarParameter]]

        SimpleParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class HardParticle(BaseParticle):
    """ Hard particle type of PeTar
        The member include BaseParticle and searching radius, id, artificial particle data and changeover radii
        
    keys: (class members)
        Members inherited from BaseParticle: see help(petar.BaseParticle)
             Please set the keyword argument 'interrupt_mode' to determine the members of stellar evolution data
        r_search (1D): searching radius
        id       (1D): identification
        mass_bk  (1D): artificial particle parameter 1 
        status   (1D): artificial particle parameter 2
        r_in     (1D): changeover function inner boundary
        r_out    (1D): changeover function outer boundary
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
               This option indicates whether columns of stellar evolution exist
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
            spin_3d: bool (True)
                when interrupt_mode is bse/bseEmp/mobse: width of the star.spin column; use False for pre-Dec-2024 outputs
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['r_search',float_type], ['id',np.int64], ['mass_bk',np.float64], ['status',np.float64], ['r_in',float_type], ['r_out',float_type]]

        BaseParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class HermiteParticle(HardParticle):
    """ Hermite particle type of PeTar
        The member include HardParticle and dt, time, acc, jerk and pot
        
    keys: (class members)
        Members inherited from HardParticle: see help(petar.HardParticle)
             Please set the keyword argument 'interrupt_mode' to determine the members of stellar evolution data
        dt    (1D): time step
        time  (1D): current time
        acc   (2D,3): acceleration x, y, z
        jerk  (2D,3): acceleration derivative x, y, z
        pot   (1D): potential
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
               This option indicates whether columns of stellar evolution exist
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
            spin_3d: bool (True)
                when interrupt_mode is bse/bseEmp/mobse: width of the star.spin column; use False for pre-Dec-2024 outputs
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['dt',float_type],['time',float_type],['acc',(float_type,3)],['jerk',(float_type,3)],['pot',float_type]]

        HardParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

@enablePathInit
class Particle(HardParticle):
    """ (Soft) Particle type of PeTar, also used in snapshot
        The particle data of PeTar. Depending on the compile configuration of PeTar, 
        The data structures (columns) of the particle snapshots are different.
        Using the correct keyword arguments in the initialization to control the member definition (Keys)

    keys: (class members)
        Members inherited from HardParticle: see help(petar.HardParticle)
             Please set the keyword argument 'interrupt_mode' to determine the members of stellar evolution data
             Please set the keyword argument 'use_mpfrc' to determine whether high-precision parts of particle position are included
        acc_soft (2D,3): long-range interaction acceleration (particle-tree) x, y, z
        if (keyword argument 'collect_sp_acc' == True'):
            acc_sp (2D,3): superparticle acceleration (only used when superparticle is enabled) x, y, z
        pot      (1D): total potential
        pot_soft (1D): long-range interaction potential
        if (keyword argument 'external_mode' != 'none'):
             pot_ext  (1D): external potential (only exist when keyword argument 'external_mode' is not 'none')
        n_nb:    (1D): number of neighbors (short-interaction)
    """

    # snapshot-level layout parameters on top of the particle-chain ones;
    # see BaseParticle.LAYOUT_PARAMS for the mechanism
    LAYOUT_PARAMS = dict(HardParticle.LAYOUT_PARAMS,
                         external_mode=('none', 'galpy', 'agama'),
                         collect_sp_acc=(False, True))

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
               This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
               PeTar external mode (set in configure): galpy, agama, none 
               This option indicates whether the column of externa potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            collect_sp_acc: bool (False)
               If true, the superparticle acceleration is collected and the column acc_sp exists            spin_3d: bool (True)
                when interrupt_mode is bse/bseEmp/mobse: width of the star.spin column; use False for pre-Dec-2024 outputs            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['acc_soft',(float_type,3)]]
        if ('collect_sp_acc' in kwargs.keys()):
            if (kwargs['collect_sp_acc']):
                keys += [['acc_sp',(float_type,3)]]
        keys += [['pot',float_type], ['pot_soft',float_type]]
        if ('external_mode' in kwargs.keys()):
            if (kwargs['external_mode']!='none'):
                keys += [['pot_ext',float_type]]
        keys += [['n_nb',np.int64]]

        HardParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

    def calcEtot(self):
        """ Calculate total energy and add it as the member, etot
        """
        etot = self.ekin + self.mass*self.pot
        self.addNewMember('etot',etot)

@enablePathInit
class Binary(ParticleRecordLayoutMixin, sdar_Binary):
    # members are Particles: the snapshot layout parameters and header-fact
    # bindings apply
    LAYOUT_PARAMS = dict(Particle.LAYOUT_PARAMS)
    HEADER_PINNED_PARAMS = dict(Particle.HEADER_PINNED_PARAMS)
    LABEL_REFINE_PARAMS = Particle.LABEL_REFINE_PARAMS
    PINNING_PARAMS = Particle.PINNING_PARAMS
    """ Binary class
        The binary (tree) data. Depending on the definition of two members 
        (keyword argument member_particle_type(|_one|_two), 
        The binary can refer to any type of multiple system.

    Keys:
        The final keys depends on kwargs of initial function
  
        kwargs['simple_mode'] (bool)
            True: (default)
                mass (1D): total mass of two components
                pos  (2D,3): c.m. position x, y, z
                vel  (2D,3): c.m. velocity vx, vy, vz
                rrel (1D): relative distance
                semi (1D): semi-major axis
                ecc  (1D): eccentricity
                p1   (member_particle_type_one) component one
                p2   (member_particle_type_two) component two
            False:
                mass (1D): total mass of two components
                pos  (2D,3): c.m. position x, y, z
                vel  (2D,3): c.m. velocity vx, vy, vz
                m1   (1D): component 1 mass
                m2   (1D): component 2 mass
                rrel (1D): relative distance
                semi (1D): semi-major axis
                am   (2D,3): specific angular momemtum x, y, z
                L    (2D,3): angular momemtum x, y, z
                eccvec  (2D,3): eccentric vector
                incline (1D): inclination
                rot_horizon (1D): frame rotational angle in x-y plane (longitude of ascending node)
                ecc  (1D): eccentricity
                rot_self (1D): frame rotational angle in orbital plane (argument of periapsis)
                ecca (1D): eccentric anomaly
                period (1D): period
                t_peri (1D): time to peri-center
                p1 (member_particle_type_one) component one
                p2 (member_particle_type_two) component two

        The member_particle_type(|_one|_two) is given by keyword arguments:
           'member_particle_type' (for both members),'member_particle_type_one','member_particle_type_two'.
        In default, it is petar.Particle.
        If a type (e.g., petar.Particle) is given, the member is a single star.
        If a list with two members (e.g., [petar.Particle, petar.Particle]) is given, 
        the member is a binary with two single stars.
        A hierarchical list can be provided, e.g., [petar.Particle, [petar.Particle, petar.Particle]]
        to indicate a triple system.
    """

    def __init__ (self, _p1=None, _p2=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _p1: particle data | 2D numpy.ndarray | petar.Binary | None
            If _p1 is a particle type data or a petar.Binary type data, it is treated as the first component of binary
            If _p1 is a petar.Binary type data and _p2 is None, the class instance is initialized by copy the data of _p1.
            If _p1 is None, initialize class with empty data
        _p2: particle data | None
            If _p2 is a particle type data or a petar.Binary type data, it is treated as the second component of binary
            If _p2 is None, _p1 should be petar.Binary data or None
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members

        keyword arguments:
            simple_mode: bool (True)
                If True, only calculate semi and ecc, save computing time significantly
            G: float (1.0)
                Gravitational constant
            interrupt_mode: string (none)
                PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
                This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
                This option indicates whether the column of external potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            member_particle_type: type or list (SimpleParticle)
                Type of component particle (both)
            member_particle_type_one: type or list (SimpleParticle)
                Type of 1st component
            member_particle_type_two: type or list (SimpleParticle)
                Type of 2nd component 
            float_type: type (np.float64)
                floating point data type
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_particle_type'] = Particle
        super().__init__(_p1, _p2, _offset, _append, **kwargs_local)



@enablePathInit
class InterruptBinary(SDARInterruptBinary):
    # records embed HardParticle columns: the particle-chain layout
    # parameters apply (sidecar + ASCII column-count detection)
    LAYOUT_PARAMS = dict(BaseParticle.LAYOUT_PARAMS)

    """ Data of stellar evolution interrupted binary in base mode
        Inherit from sdar.ar.SDARInterruptBinary
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ Initial InterruptBinary class
        Parameters
        ----------
        ----------
        keyword arguments:
            particle_type: type (HardParticle)
                particle data type
        """

        if (not 'particle_type' in kwargs.keys()):
            kwargs['particle_type'] = HardParticle
        particle_type = kwargs['particle_type']

        SDARInterruptBinary.__init__(self, _dat, _offset, _append, **kwargs)

def calculateParticleCMDict(pcm, _p1, _p2):
    """ Calculate the center-of-the-mass of two particle sets
    
    Parameters
    ----------
    _p1: inherited SimpleParticle
        particle set 1
    _p2: inherited SimpleParticle 
        particle set 2, should have the same size as _p1
    pcm: dict 
        particle center-of-the-mass, should include keys: 'mass','pos','vel'.
    """
    if (issubclass(type(_p1), SimpleParticle)) & (issubclass(type(_p2),SimpleParticle)):
        pcm['mass'] = _p1.mass + _p2.mass
        pcm['pos']  = (_p1.mass[:,None]*_p1.pos + _p2.mass[:,None]*_p2.pos)/pcm['mass'][:,None]
        pcm['vel']  = (_p1.mass[:,None]*_p1.vel + _p2.mass[:,None]*_p2.vel)/pcm['mass'][:,None]
    elif (isinstance(_p1, collections.OrderedDict)) & (isinstance(_p2,collections.OrderedDict)) | (isinstance(_p1, dict)) & (isinstance(_p2, dict)):
        pcm['mass'] = _p1['mass'] + _p2['mass']
        pcm['pos']  = (_p1['mass'][:,None]*_p1['pos'] + _p2['mass'][:,None]*_p2['pos'])/pcm['mass'][:,None]
        pcm['vel']  = (_p1['mass'][:,None]*_p1['vel'] + _p2['mass'][:,None]*_p2['vel'])/pcm['mass'][:,None]
    else:
        raise ValueError('Initial fail, date type should be Particle or collections.OrderDict, given',type(_p1))

@enablePathInit
class GroupInfo(ParticleRecordLayoutMixin, hermite_group.GroupInfo):
    # members are HardParticles: the particle-chain layout parameters and
    # header-fact bindings apply
    LAYOUT_PARAMS = dict(HardParticle.LAYOUT_PARAMS)
    HEADER_PINNED_PARAMS = dict(HardParticle.HEADER_PINNED_PARAMS)
    LABEL_REFINE_PARAMS = HardParticle.LABEL_REFINE_PARAMS
    PINNING_PARAMS = HardParticle.PINNING_PARAMS
    """ Group information output from PeTar
    Keys: (class members)
        type (1D): group type, 0: new group; 1: end group
        n    (1D): number of members in group (should be consistent with keyword argument N
        time (1D): current time
        pos  (2D,3): position of the group c.m. in the framework of the global system (without shift of global system c.m. if external_mode is on)
        vel  (2D,3): velocity of the group c.m. in the framework of the global system (without shift of global system c.m.)
        bin[X] (BinaryTreeSDAR): members of the group in a hierarchical binary tree
               Here X indicates the order. 0 represents the root (outer most) binary; 1,2,3 ... are inner binaries
               For a triple, bin0 is outer binary, bin1 is inner binary.
               p2 of bin0 is the c.m. of bin1, the id of p2 is the minimum id from the two components in bin1.
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ GroupInfo class for petar, inherit from hermite_group.GroupInfo

            Parameters
            ----------
            keyword arguments:
                member_particle_type: type (HardParticle)
                    Type of component particle, do not change this!
                interrupt_mode: string (none)
                    PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
                    This option indicates whether columns of stellar evolution exist
                external_mode: string (none)
                    PeTar external mode (set in configure): galpy, agama, none 
                    This option indicates whether the column of external potential exist
                use_mpfrc: bool (False)
                    If true, add three columns of pos_high indicating the high-precision parts of position
                spin_3d: bool (True)
                    when interrupt_mode is bse/bseEmp/mobse: width of the star.spin column; use False for pre-Dec-2024 outputs
                float_type: type (np.float64)
                    floating point data type
                N: int (2)
                    Number of members of one group
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_particle_type'] = HardParticle
        super().__init__(_dat, _offset, _append, **kwargs_local)

    
@enablePathInit
class SDARData(sdar_SDARData):
    """ SDARData class for petar, inherit from sdar.SDARData

        Keys: (class members)
            time (1D): current evolved time (counting from zero)
            de (1D): physical energy error
            etot_ref (1D): initial total energy
            ekin (1D): kinetic energy
            epot (1D): potential energy
            gt_drift (1D): time tranformation for drift step
            H (1D): extened phase space Hamiltonian
            if (keyword argument 'include_H_approx' == True):
                H_approx (1D): approximated phase space Hamiltonian
            de_interrupt (1D): energy change due to interruption
            dH_interrupt (1D): H change due to interruption
            perturber (perturber_type): perturber data, depending on the keyword argument 'perturber_type'
                                        if not given, this key is not included
            info (SDARInfo): SDAR information shown as follows:
                ds (1D): integration step
                time_offset (1D): time offset to obtain the actual time (time_offset + time)
                r_break_crit (1D): distance criterion to break group (used in Hermite)
            if (keyword argument 'g_func' == True):
                g_func (1D): g-function active state: 0=standard LogH, 1=the g-function method of this build (per-binary macro; old data files used codes 1-4)
            profile (SDARProfile): SDAR profile
            if (keyword argument 'slowdown' == True):
                de_sd (1D): slowdown energy error
                etot_sd (1D): slowdown energy
                ekin_sd (1D): slowdown kinetic energy
                epot_sd (1D): slowdown potential energy
                de_sd_change (1D): slowdown energy change
                dH_sd_change (1D): slowdown H change
                de_sd_interrupt (1D): slowdown energy change due to interruption
                dH_sd_interrupt (1D): slowdown H change due to interruption
                sd (SlowDownGroup): slowdown data
            particles (ParticleGroup): particle group, depending on the keyword argument 'member_type' and 'cm_type'

    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ SDARData class for petar, inherit from sdar.SDARData

            Parameters
            ----------
            keyword arguments:
                member_type: member particle type (HardParticle)
                cm_type: c.m. particle type (HermiteParticle)
                perturber_type: perturber particle type (None)
                N_particle: int (0)
                    Number of particles, determined from file if not provided
                slowdown: bool (False)
                    if True, add slowdown keys
                N_sd: int (0)
                    Number of slowdown pair, used when slowdown='on'
                time_measure: bool (False)
                    if True, add time measure keys in profile
                include_H_approx: bool (False)
                    if True, add H_approx key
                interrupt_mode: string (none)
                    PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
                    This option indicates whether columns of stellar evolution exist
                external_mode: string (none)
                    PeTar external mode (set in configure): galpy, agama, none 
                    This option indicates whether the column of external potential exist
                use_mpfrc: bool (False)
                    If true, add three columns of pos_high indicating the high-precision parts of position
                float_type: type (np.float64)
                    floating point data type
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_type'] = HardParticle
        kwargs_local['cm_type'] = HermiteParticle

        super().__init__(_dat, _offset, _append, **kwargs_local)

@enablePathInit
class HermiteData(hermite_HermiteData):
    """ HermiteData class for petar, inherit from sdar.HermiteData
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        time_offset (1D): time offset to calculate the global time (time+time_offset)
        energy_phy (HermiteEnergy): physical energy data
        energy_sd (HermiteEnergy): slowdown energy data
        sd (SlowDownGroup): slowdown data
        profile (HermiteProfile): hermite profile
        particles (ParticleGroup): particle group, depending on the keyword argument 'member_type' and 'cm_type'
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ HermiteData class for petar, inherit from sdar.HermiteData

            Parameters
            ----------
            keyword arguments:
                member_type: member particle type (HermiteParticle)
                cm_type: c.m. particle type (HermiteParticle)
                N_particle: int (0)
                    Number of members of one group
                N_sd: int (0)
                    Number of slowdown pairs
                time_measure: bool (False)
                    if True, add time measure keys in profile
                interrupt_mode: string (none)
                    PeTar interrupt mode (set in configure): merger, bse, mobse, bseEmp, sevn, none, dsm
                    This option indicates whether columns of stellar evolution exist
                external_mode: string (none)
                    PeTar external mode (set in configure): galpy, agama, none 
                    This option indicates whether the column of external potential exist
                use_mpfrc: bool (False)
                    If true, add three columns of pos_high indicating the high-precision parts of position
                float_type: type (np.float64)
                    floating point data type
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_type'] = HermiteParticle
        kwargs_local['cm_type'] = HermiteParticle

        super().__init__(_dat, _offset, _append, **kwargs_local)

# ---------------------------------------------------- snapshot + conveniences

class Snapshot:
    """ A PeTar snapshot: header + particles in one object.

    ``petar.Snapshot('data.0')`` reads both with the auto-detected layout
    (sidecar first, structure sniff fallback; explicit kwargs win) - this is
    the dedicated home of the headered-snapshot reading path; header format
    knowledge lives in PeTarDataHeader, layout knowledge in the particle
    classes' LAYOUT_PARAMS.

    members:
        particle (petar.Particle): particle records
        header (petar.PeTarDataHeader): snapshot header (file id, N, time,
               and the c.m. offsets when the producing run used an external
               potential)
    """

    def __init__(self, path=None, **kwargs):
        if path is None:
            self.particle = Particle(**kwargs)
            self.header = PeTarDataHeader()
            return
        path = str(path)
        self.particle = Particle(path, **kwargs)
        external = self.particle.initargs.get('external_mode', 'none')
        if external == 'none':
            external = 'none'
        self.header = PeTarDataHeader(path,
                                      snapshot_format='ascii' if detectAsciiFile(path) else 'binary',
                                      external_mode=external,
                                      float_type=self.particle.initargs.get('float_type', np.float64))

    @classmethod
    def fromRecord(cls, particle, header):
        """ Construct a Snapshot from an already-read particle array and
        header - the entry point used by the packed-container reader
        (``PackedSnapshots.readSnapshot``), which seeks one record out of a
        volume file instead of reading a standalone snapshot file. The
        result is indistinguishable from ``Snapshot(path)`` of the
        extracted record. """
        snap = cls()
        snap.particle = particle
        snap.header = header
        return snap


def packVolumePath(prefix, volume):
    """ Volume filename of the packed snapshot container: ``[prefix].pack.<id>``
    (solver option ``-w 4``). """
    return '%s.pack.%04d' % (str(prefix), int(volume))


def readPackIndex(path):
    """ Read the packed-container index (``[prefix].pack.idx``) into a numpy
    structured array: one row per packed snapshot record with columns
    (nfile, time, n_body, volume, nbytes, offset), ordered by nfile. """
    path = str(path)
    if not path.endswith('.pack.idx'):
        path = path + '.pack.idx'
    dtype = [('nfile', np.int64), ('time', np.float64), ('n_body', np.int64),
             ('volume', np.int64), ('nbytes', np.int64), ('offset', np.int64)]
    rows = []
    with open(path) as fp:
        for lineno, line in enumerate(fp, 1):
            line = line.strip()
            if (not line) or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 6:
                raise ValueError('%s line %d: expected 6 columns (nfile time n_body volume nbytes offset), got %d: %r'
                                 % (path, lineno, len(parts), line))
            try:
                rows.append((int(parts[0]), float(parts[1]), int(parts[2]),
                             int(parts[3]), int(parts[4]), int(parts[5])))
            except ValueError:
                raise ValueError('%s line %d: non-numeric field in %r' % (path, lineno, line)) from None
    if not rows:
        return np.array([], dtype=dtype)
    return np.array(rows, dtype=dtype)


class PackedSnapshots:
    """ Reader front-end for the packed snapshot container (solver option
    ``-w 4``): append-only volume files ``[prefix].pack.<id>`` holding
    standalone BINARY snapshot records, plus the index
    ``[prefix].pack.idx`` (see ``readPackIndex``).

    The reader layout is resolved from the run sidecar (``[prefix].meta``)
    via ``petar.getDataKwargs`` or from explicit keyword arguments, then
    snapshots are read by seek without unpacking::

        pack = petar.PackedSnapshots('data')
        snap = pack.readSnapshot(-1)       # last snapshot (petar.Snapshot)
        snap = pack.readSnapshot(nfile=3)  # by snapshot id
        snap = pack[0]                     # by row index
        snap.particle, snap.header         # like petar.Snapshot(path)
    """

    IDX_DTYPE = [('nfile', np.int64), ('time', np.float64), ('n_body', np.int64),
                 ('volume', np.int64), ('nbytes', np.int64), ('offset', np.int64)]

    def __init__(self, prefix, **kwargs):
        prefix = str(prefix)
        if prefix.endswith('.pack.idx'):
            prefix = prefix[:-len('.pack.idx')]
        self.prefix = prefix
        self.index = readPackIndex(prefix)
        merged = dict(kwargs)
        if not all(p in merged for p in ('interrupt_mode', 'external_mode', 'float_type')):
            # Resolve the layout from the run sidecar. getDataKwargs needs a
            # real file path and the bare prefix is not one, so read the
            # sidecar directly; fall back to the index file path (whose
            # generic prefix walk locates the sidecar).
            meta = None
            try:
                meta = readMetaData(prefix + '.meta')
            except Exception:
                meta = None
            if meta:
                if meta.get('layout.interrupt_mode'):
                    merged.setdefault('interrupt_mode', meta['layout.interrupt_mode'])
                if meta.get('layout.external_mode'):
                    merged.setdefault('external_mode', meta['layout.external_mode'])
                for key in ('use_mpfrc', 'collect_sp_acc', 'spin_3d'):
                    tag = 'layout.' + key
                    if tag in meta:
                        merged.setdefault(key, meta[tag] in ('1', 'true', 'True'))
                if meta.get('layout.float_type') == 'float128':
                    merged.setdefault('float_type', np.float128)
            else:
                try:
                    detected = getDataKwargs(prefix + '.pack.idx')
                    for k, v in detected.items():
                        merged.setdefault(k, v)
                except Exception:
                    pass
        self.kwargs = merged
        self.header_size, self.header_dtype = self._headerLayout(merged)

    @staticmethod
    def _headerLayout(kwargs):
        """ (header bytes, header dtype) for the record header variant pinned
        by the layout kwargs (c.m. offsets present iff external_mode is not
        none; float128 time when use_mpfrc). """
        ft = kwargs.get('float_type', np.float64)
        cm = kwargs.get('external_mode', 'none') != 'none'
        if ft == np.float128:
            hsize = HEADER_OFFSET_WITH_CM_F128 if cm else HEADER_OFFSET_F128
        else:
            hsize = HEADER_OFFSET_WITH_CM if cm else HEADER_OFFSET
        fields = [('nfile', np.int64), ('n_body', np.int64), ('time', ft)]
        if cm:
            fields += [('pos_offset', ft, (3,)), ('vel_offset', ft, (3,))]
        return hsize, np.dtype(fields)

    def __len__(self):
        return len(self.index)

    def volumePath(self, volume):
        return packVolumePath(self.prefix, volume)

    def _rowIndex(self, row):
        row = int(row)
        if row < 0:
            row += len(self.index)
        if not (0 <= row < len(self.index)):
            raise IndexError('packed index row %d out of range (0..%d)' % (row, len(self.index) - 1))
        return row

    def readHeader(self, row):
        """ Snapshot header of index row ``row`` (petar.PeTarDataHeader). """
        row = self._rowIndex(row)
        rec = self.index[row]
        with open(self.volumePath(rec['volume']), 'rb') as fp:
            fp.seek(int(rec['offset']))
            buf = fp.read(self.header_size)
        if len(buf) != self.header_size:
            raise IOError('truncated packed record header: volume %s offset %d'
                          % (rec['volume'], rec['offset']))
        vals = np.frombuffer(buf, dtype=self.header_dtype)[0]
        header = PeTarDataHeader()
        header.file_id = int(vals['nfile'])
        header.n = int(vals['n_body'])
        header.time = float(vals['time'])
        if 'pos_offset' in self.header_dtype.names:
            header.pos_offset = [float(x) for x in vals['pos_offset']]
            header.vel_offset = [float(x) for x in vals['vel_offset']]
            header.offset_flag = True
        return header

    def _selectRow(self, row, nfile, time):
        """ Resolve the selection arguments (row index / snapshot id /
        first row at-or-after a time) to an index row. """
        if len(self.index) == 0:
            raise ValueError('the packed container index is empty')
        if nfile is not None:
            sel = np.where(self.index['nfile'] == int(nfile))[0]
            if len(sel) == 0:
                raise ValueError('nfile %d not found in the packed index' % nfile)
            row = int(sel[0])
        elif time is not None:
            sel = np.where(self.index['time'] >= float(time))[0]
            row = int(sel[0]) if len(sel) > 0 else len(self.index) - 1
        elif row is None:
            raise ValueError('one of row/nfile/time is required')
        return self._rowIndex(row)

    def readSnapshot(self, row=None, nfile=None, time=None):
        """ Read one packed snapshot as a ``petar.Snapshot`` object - the
        same type ``petar.Snapshot(path)`` returns for a standalone
        snapshot file. The record header and the particle block are seeked
        out of the volume file by the index (layout resolved from the run
        sidecar; no unpacking). Selection: row index (negative counts from
        the end), snapshot id (``nfile``), or the first row with ``time``
        >= the given time (falls back to the last row). """
        r = self._selectRow(row, nfile, time)
        rec = self.index[r]
        header = self.readHeader(r)
        particle = Particle(**self.kwargs)
        # the volume tail beyond this record consists of LATER records, so
        # the file tail past offset+count is never itemsize-aligned for
        # non-last records; the index pins the exact byte range
        # (strict_mismatch=False), and the associated generic
        # misalignment/truncation notice is suppressed for this call only
        import warnings as _w
        with _w.catch_warnings():
            _w.filterwarnings('ignore', message=r'Reading \S+: binary file size .* not aligned.*')
            particle.fromfile(self.volumePath(rec['volume']),
                              offset=int(rec['offset']) + self.header_size,
                              count=int(rec['n_body']),
                              strict_mismatch=False)
        return Snapshot.fromRecord(particle, header)

    def read(self, row=None, nfile=None, time=None):
        """ Tuple form returning (particle, header); prefer ``readSnapshot``
        (or indexing ``pack[i]``), which returns the ``petar.Snapshot``
        object. """
        snap = self.readSnapshot(row=row, nfile=nfile, time=time)
        return snap.particle, snap.header

    def __getitem__(self, row):
        return self.readSnapshot(row=row)


def getDataKwargs(path, return_provenance=False, meta='auto'):
    """ Layout kwargs for a snapshot-family file via the Particle concise
    machinery (sidecar first, structure sniff fallback). Raises
    AmbiguousSchemaError when neither source resolves the layout. """
    kwargs = {'meta': meta}
    _dat, kwargs, _path = Particle.resolvePathInput(str(path), kwargs)
    if return_provenance:
        import os
        if (meta is not False) and (meta != 'none'):
            mf = str(meta) if (meta is not None and meta != 'auto') else Particle.findLayoutMeta(str(path))
            if mf is not None and os.path.isfile(mf):
                return kwargs, 'sidecar %s' % mf
        return kwargs, 'structure'
    return kwargs


def readMetaData(path):
    """ Parse an output-meta sidecar into a plain all-keys dict (producer
    information included); accepts the shared output-meta/1 and the
    pre-release petar-meta/1 id. Accepts either the sidecar path itself or
    any run-family data file (the sidecar is located by the generic prefix
    walk); returns None when no sidecar exists. Raises ValueError for other
    formats. """
    import warnings as _w
    p = str(path)
    if not p.endswith('.meta'):
        from sdar.base import DictNpArrayMix
        p = DictNpArrayMix.findLayoutMeta(p)
        if p is None:
            return None
    meta = {}
    with open(p) as fp:
        for line in fp:
            line = line.strip()
            if (not line) or line.startswith('#') or ('=' not in line):
                continue
            key, value = (part.strip() for part in line.split('=', 1))
            meta[key] = value
    fmt = meta.get('format', '')
    if not (fmt.startswith('output-meta/') or fmt.startswith('petar-meta/')):
        raise ValueError('%s is not an output-meta sidecar (format = %r)' % (path, fmt))
    if fmt not in ('output-meta/1', 'petar-meta/1'):
        _w.warn('%s uses sidecar format %r, but this reader implements output-meta/1; '
                'parsing may be incomplete or wrong - update the petar tools.' % (path, fmt))
    return meta


def writeMetaData(fname, interrupt_mode, external_mode='none', use_mpfrc=False,
                  collect_sp_acc=False, spin_3d=None, float_type=np.float64,
                  producer='petar', petar_version=None, binary_family=None,
                  snapshot_format=None):
    """ Write a petar output-meta/1 sidecar (layout.<param> keys, matching
    the solver's writer; used by producing tools to propagate the schema
    next to outputs written outside the run directory). """
    lines = ['# output-meta/1: output schema tag; see doc/output-schema-tagging-plan.md',
             'format = output-meta/1',
             'producer = %s' % producer]
    if petar_version is not None:
        lines.append('petar_version = %s' % petar_version)
    if binary_family is not None:
        lines.append('binary_family = %s' % binary_family)
    if snapshot_format is not None:
        lines.append('snapshot_format = %s' % snapshot_format)
    lines.append('layout.interrupt_mode = %s' % interrupt_mode)
    lines.append('layout.external_mode = %s' % external_mode)
    lines.append('layout.use_mpfrc = %s' % ('true' if use_mpfrc else 'false'))
    lines.append('layout.collect_sp_acc = %s' % ('true' if collect_sp_acc else 'false'))
    if spin_3d is not None:
        lines.append('layout.spin_3d = %s' % ('true' if spin_3d else 'false'))
    lines.append('layout.float_type = %s' % ('float128' if float_type == np.float128 else 'float64'))
    with open(fname, 'w') as fp:
        fp.write('\n'.join(lines) + '\n')


def autoFillKwargs(sample_path, kwargs=None, verbose=True, label=None):
    """ CLI helper: fill missing layout kwargs from a sample file via
    Particle.resolvePathInput (detection without reading). """
    return Particle.autoFillKwargs(sample_path, kwargs, verbose=verbose, label=label)
