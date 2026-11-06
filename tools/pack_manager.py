#!/usr/bin/env python3
"""petar.pack - manage the packed snapshot container produced by the solver
option ``-w 4``: append-only volume files ``[prefix].pack.<id>`` holding
standalone BINARY snapshot records, plus the index ``[prefix].pack.idx``
(one line per record: nfile time n_body volume nbytes offset).

Subcommands:
  list     print the index table
  extract  materialize packed snapshots as standalone snapshot files
           (``<prefix>.<N>``, byte-identical to the -w 1 output)
  trim     drop packed records after a given time/snapshot id
           (restart preparation, like petar.data.clear for event streams)
  rebuild  regenerate the index by scanning the volume files
           (crash recovery; the index is derivable from the volumes)
  hdf5     offline conversion of the packed snapshots into a single HDF5
           file (requires h5py; used for distribution/archiving/interop)

See doc/packed-output_notes.md for the container layout.
"""
import argparse
import os
import sys

import numpy as np

import petar


def _load_pack(args):
    """ PackedSnapshots for the given prefix with optional explicit layout
    kwargs (falling back to the [prefix].meta sidecar). """
    kwargs = {}
    if getattr(args, 'interrupt_mode', None):
        kwargs['interrupt_mode'] = args.interrupt_mode
    if getattr(args, 'external_mode', None):
        kwargs['external_mode'] = args.external_mode
    if getattr(args, 'use_mpfrc', None):
        kwargs['use_mpfrc'] = True
        kwargs['float_type'] = np.float128
    return petar.PackedSnapshots(args.prefix, **kwargs)


def _discover_volumes(prefix):
    """ Existing volume ids of the container, sorted ascending (supports a
    directory component in the prefix). """
    d, base = os.path.split(str(prefix))
    head = base + '.pack.'
    vols = []
    for name in os.listdir(d or '.'):
        if not name.startswith(head):
            continue
        tail = name[len(head):]
        if tail.isdigit():
            vols.append(int(tail))
    return sorted(vols)


def cmd_list(args):
    pack = _load_pack(args)
    idx = pack.index
    print('# nfile time n_body volume nbytes offset')
    for rec in idx:
        print('%d %.17g %d %d %d %d' % (rec['nfile'], rec['time'], rec['n_body'],
                                        rec['volume'], rec['nbytes'], rec['offset']))
    print('# %d record(s), %d volume(s)' % (len(idx), len(set(idx['volume'])) if len(idx) else 0))
    return 0


def _extract_record(pack, rec, out_path):
    """ Byte-copy one packed record to a standalone snapshot file. """
    src = pack.volumePath(rec['volume'])
    nbytes = int(rec['nbytes'])
    with open(src, 'rb') as fin:
        fin.seek(int(rec['offset']))
        buf = fin.read(nbytes)
    if len(buf) != nbytes:
        raise IOError('truncated packed record: volume %s offset %d expects %d bytes, got %d'
                      % (src, rec['offset'], nbytes, len(buf)))
    with open(out_path, 'wb') as fout:
        fout.write(buf)
    # consistency: the record header must carry the index nfile/n_body
    head = np.frombuffer(buf[:pack.header_size], dtype=pack.header_dtype)[0]
    if int(head['nfile']) != int(rec['nfile']) or int(head['n_body']) != int(rec['n_body']):
        raise IOError('packed record header mismatch at volume %s offset %d: header (nfile %d, N %d) vs index (nfile %d, N %d)'
                      % (src, rec['offset'], head['nfile'], head['n_body'], rec['nfile'], rec['n_body']))
    return out_path


def cmd_extract(args):
    pack = _load_pack(args)
    idx = pack.index
    if len(idx) == 0:
        print('Error: the packed container index is empty.')
        return 1
    if args.snap_lst and not args.all:
        # validate before touching any output file
        print('Error: --snap-lst requires --all (the list is rewritten from the full container index).')
        return 1
    if args.all:
        rows = range(len(idx))
    elif args.nfile is not None:
        sel = np.where(idx['nfile'] == int(args.nfile))[0]
        if len(sel) == 0:
            print('Error: nfile %d not found in the packed index.' % args.nfile)
            return 1
        rows = [int(sel[0])]
    elif args.time is not None:
        sel = np.where(idx['time'] <= float(args.time))[0]
        if len(sel) == 0:
            print('Error: no packed record at or before time %g.' % args.time)
            return 1
        rows = [int(sel[-1])]
    else:
        rows = [len(idx) - 1]
    out_prefix = args.output if args.output else pack.prefix
    names = []
    for row in rows:
        rec = idx[row]
        out_path = '%s.%d' % (out_prefix, int(rec['nfile']))
        _extract_record(pack, rec, out_path)
        names.append(os.path.basename(out_path))
        print('extracted %s (time %.17g, N %d) -> %s'
              % (pack.volumePath(rec['volume']), rec['time'], rec['n_body'], out_path))
    if args.snap_lst:
        if not args.all:
            print('Error: --snap-lst requires --all (the list is rewritten from the full container index).')
            return 1
        # rewrite from the full index: idempotent across repeated extract runs
        lst_path = out_prefix + '.snap.lst'
        with open(lst_path, 'w') as fp:
            for rec in idx:
                fp.write('%s.%d\n' % (out_prefix, int(rec['nfile'])))
        print('snapshot list rewritten: %s (%d entries)' % (lst_path, len(idx)))
    return 0

def cmd_trim(args):
    if args.time is None and args.nfile is None:
        print('Error: trim requires --time or --nfile.')
        return 1
    pack = _load_pack(args)
    idx = pack.index
    if len(idx) == 0:
        print('the packed container index is empty; nothing to trim.')
        return 0
    if args.time is not None:
        keep = np.where(idx['time'] <= float(args.time) + 1e-12 * max(1.0, abs(float(args.time))))[0]
    else:
        keep = np.where(idx['nfile'] <= int(args.nfile))[0]
    if len(keep) == 0:
        print('Error: no packed record is at or before the trim point; refusing to trim everything (remove the container files to discard the whole container).')
        return 1
    if len(keep) == len(idx):
        print('nothing to trim: all %d record(s) are at or before the trim point.' % len(idx))
        return 0
    last = idx[int(keep[-1])]
    drop = len(idx) - len(keep)
    print('trimming %d record(s): keeping nfile %d (time %.17g) as the last record'
          % (drop, last['nfile'], last['time']))
    if args.dry_run:
        for rec in idx[int(keep[-1]) + 1:]:
            print('  would drop nfile %d (time %.17g, volume %d)' % (rec['nfile'], rec['time'], rec['volume']))
        print('  would truncate volume %s to %d bytes' % (pack.volumePath(last['volume']), int(last['offset']) + int(last['nbytes'])))
        return 0
    # drop whole volumes above the last kept record's volume
    for vol in _discover_volumes(pack.prefix):
        if vol > int(last['volume']):
            path = pack.volumePath(vol)
            os.remove(path)
            print('removed volume %s' % path)
    # truncate the last kept record's volume to its end
    vol_path = pack.volumePath(last['volume'])
    end = int(last['offset']) + int(last['nbytes'])
    if os.path.getsize(vol_path) > end:
        with open(vol_path, 'r+b') as fp:
            fp.truncate(end)
        print('truncated volume %s to %d bytes' % (vol_path, end))
    # rewrite the index with the kept rows
    idx_path = pack.prefix + '.pack.idx'
    with open(idx_path, 'w') as fp:
        for row in keep:
            rec = idx[int(row)]
            fp.write('%d %.17g %d %d %d %d\n' % (rec['nfile'], rec['time'], rec['n_body'],
                                                 rec['volume'], rec['nbytes'], rec['offset']))
    print('index rewritten: %s (%d record(s))' % (idx_path, len(keep)))
    return 0


def cmd_rebuild(args):
    pack = _load_pack(args)
    vols = _discover_volumes(pack.prefix)
    if not vols:
        print('Error: no packed volume files found for prefix %s.' % pack.prefix)
        return 1
    itemsize = np.dtype(petar.Particle(**pack.kwargs).collectDtype()).itemsize
    hsize, hdtype = pack.header_size, pack.header_dtype
    rows = []
    expect_nfile = None
    for vol in vols:
        path = pack.volumePath(vol)
        size = os.path.getsize(path)
        pos = 0
        while pos + hsize <= size:
            with open(path, 'rb') as fp:
                fp.seek(pos)
                head = np.frombuffer(fp.read(hsize), dtype=hdtype)[0]
            nfile = int(head['nfile'])
            n_body = int(head['n_body'])
            if n_body <= 0:
                raise IOError('%s offset %d: non-positive particle count %d - wrong header layout? (check reader kwargs)'
                              % (path, pos, n_body))
            if expect_nfile is not None and nfile != expect_nfile:
                raise IOError('%s offset %d: nfile %d breaks the sequence (expected %d) - wrong header layout? (check reader kwargs); a stale tail from restarting an older snapshot carries duplicate/decreasing nfile - trim first (petar.pack trim)'
                              % (path, pos, nfile, expect_nfile))
            nbytes = hsize + n_body * itemsize
            if pos + nbytes > size:
                raise IOError('%s offset %d: record needs %d bytes but only %d available - truncated tail record'
                              % (path, pos, nbytes, size - pos))
            rows.append((nfile, float(head['time']), n_body, vol, nbytes, pos))
            expect_nfile = nfile + 1
            pos += nbytes
        if pos != size:
            raise IOError('%s: %d trailing bytes after the last complete record - the volume is truncated or the layout kwargs are wrong'
                          % (path, size - pos))
    idx_path = pack.prefix + '.pack.idx'
    if os.path.exists(idx_path) and not args.force:
        old = petar.readPackIndex(pack.prefix)
        new = np.array(rows, dtype=petar.PackedSnapshots.IDX_DTYPE)
        if len(old) == len(new) and np.all(old['nfile'] == new['nfile']) and np.all(old['offset'] == new['offset']):
            print('existing index %s matches the scan result (%d record(s)); nothing to do (use --force to rewrite).' % (idx_path, len(new)))
            return 0
        print('Warning: existing index %s differs from the scan result; rewriting.' % idx_path)
    with open(idx_path, 'w') as fp:
        for rec in rows:
            fp.write('%d %.17g %d %d %d %d\n' % rec)
    print('index rebuilt: %s (%d record(s) from %d volume(s))' % (idx_path, len(rows), len(vols)))
    return 0


def _write_member(group, name, value, compression):
    """ Write one particle member (recursing into nested containers) as a
    dataset; vector members keep their (N, k) shape. """
    if hasattr(value, 'keys') and not isinstance(value, np.ndarray):
        for entry in value.keys:
            sub_name, _sub_type = entry[0], entry[1]
            _write_member(group, '%s.%s' % (name, sub_name), getattr(value, sub_name), compression)
        return
    arr = np.asarray(value)
    if arr.dtype.kind not in ('f', 'i', 'u'):
        return
    kwargs = {}
    if compression:
        kwargs['compression'] = compression
    group.create_dataset(name, data=arr, **kwargs)


def cmd_hdf5(args):
    try:
        import h5py
    except ImportError:
        print('Error: h5py is required for the hdf5 subcommand (pip install h5py).')
        return 1
    pack = _load_pack(args)
    idx = pack.index
    if len(idx) == 0:
        print('Error: the packed container index is empty.')
        return 1
    out = args.output if args.output else pack.prefix + '.h5'
    n_written = 0
    compression = args.compression if args.compression and args.compression != 'none' else None
    with h5py.File(out, 'w') as h5:
        try:
            meta = petar.readMetaData(pack.prefix)
        except Exception:
            meta = None
        if meta:
            for key, value in meta.items():
                h5.attrs[key] = value
        h5.attrs['source'] = 'petar packed snapshot container'
        for row in range(len(idx)):
            rec = idx[row]
            particle, header = pack.read(row=row)
            group = h5.create_group('snap_%d' % int(rec['nfile']))
            group.attrs['nfile'] = int(rec['nfile'])
            group.attrs['time'] = float(rec['time'])
            group.attrs['n_body'] = int(rec['n_body'])
            group.attrs['pos_offset'] = np.asarray(header.pos_offset, dtype=np.float64)
            group.attrs['vel_offset'] = np.asarray(header.vel_offset, dtype=np.float64)
            for entry in particle.keys:
                name = entry[0]
                _write_member(group, name, getattr(particle, name), compression)
            n_written += 1
        print('wrote %d snapshot group(s) -> %s%s'
              % (n_written, out, (' (compression: %s)' % compression) if compression else ''))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=__doc__)
    sub = parser.add_subparsers(dest='command')

    def add_layout_args(p):
        p.add_argument('prefix', help='run output prefix (e.g. data); the container is [prefix].pack.<id> + [prefix].pack.idx')
        p.add_argument('-i', '--interrupt-mode', default=None,
                       help='reader interrupt mode (default: from [prefix].meta or auto-detect)')
        p.add_argument('-e', '--external-mode', default=None,
                       help='reader external-potential mode (default: from [prefix].meta or auto-detect)')
        p.add_argument('--use-mpfrc', action='store_true',
                       help='reader float128 layout (default: from [prefix].meta or auto-detect)')

    p_list = sub.add_parser('list', help='print the index table')
    add_layout_args(p_list)
    p_list.set_defaults(func=cmd_list)

    p_ext = sub.add_parser('extract', help='materialize packed snapshots as standalone snapshot files')
    add_layout_args(p_ext)
    p_ext.add_argument('-o', '--output', default=None,
                       help='output filename prefix (default: the container prefix; files <prefix>.<N>)')
    sel = p_ext.add_mutually_exclusive_group()
    sel.add_argument('--nfile', type=int, default=None, help='extract the snapshot with this id')
    sel.add_argument('--time', type=float, default=None, help='extract the last snapshot at or before this time')
    sel.add_argument('--all', action='store_true', help='extract all packed snapshots')
    sel.add_argument('--last', action='store_true', help='extract the last packed snapshot (default)')
    p_ext.add_argument('--snap-lst', action='store_true',
                       help='append the extracted filenames to [prefix].snap.lst (petar.movie/data.process compatible)')
    p_ext.set_defaults(func=cmd_extract, all=False, nfile=None, time=None)

    p_trim = sub.add_parser('trim', help='drop packed records after a time/snapshot id (restart preparation)')
    add_layout_args(p_trim)
    p_trim.add_argument('-t', '--time', type=float, default=None, help='keep records with time <= t')
    p_trim.add_argument('--nfile', type=int, default=None, help='keep records with nfile <= n')
    p_trim.add_argument('--dry-run', action='store_true', help='print the actions without performing them')
    p_trim.set_defaults(func=cmd_trim)

    p_rb = sub.add_parser('rebuild', help='regenerate the index by scanning the volume files (crash recovery)')
    add_layout_args(p_rb)
    p_rb.add_argument('--force', action='store_true', help='rewrite the index even when it matches the scan')
    p_rb.set_defaults(func=cmd_rebuild)

    p_h5 = sub.add_parser('hdf5', help='offline conversion of the packed snapshots into a single HDF5 file')
    add_layout_args(p_h5)
    p_h5.add_argument('-o', '--output', default=None, help='output filename (default: [prefix].h5)')
    p_h5.add_argument('--compression', default=None,
                      help='h5py compression filter for datasets (e.g. gzip, lzf; default: none)')
    p_h5.set_defaults(func=cmd_hdf5)

    args = parser.parse_args(argv)
    if not getattr(args, 'command', None):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
