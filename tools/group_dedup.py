#!/usr/bin/env python3
"""Deduplicate and repair PeTar binary group output files ([prefix].group.n[N]).

Two pathologies are handled:

1. Episode (oscillation) duplication: the same particle pair repeatedly
   form/end within a short time window. Causes include group-detection
   boundary flapping of tight eccentric binaries and crash-restart loops
   that re-integrate the same time window (events are re-generated, not
   byte-identical, so hash-based dedup does not work). Each episode is
   collapsed to its first record (initial form) and last record (final
   state); episodes of a pair are separated when the time gap exceeds
   --gap (in the time unit of the file).

2. Phase-shifted appends: when a run crashes in the middle of writing a
   record, the leftover fragment shifts the phase of all subsequent
   appended records. Valid records inside misaligned regions are detected
   (header scan at all 8-byte phases), recovered and realigned; corrupt
   fragments are dropped and reported.

The binary record layout (record size, member-id byte offsets) is derived
from the petar package (GroupInfo.collectDtype) for any member count N and
interrupt/external/spin/mpfrc/float-type options; when petar or numpy is
not importable, an empirical fallback detects the layout from the data.
The external_mode option is accepted for interface parity with other petar
readers but does not change the group-record layout.

Only BINARY group output (write-group-info 2) is supported.

Usage: petar.group.dedup [options] data.group.n2
"""

import getopt
import hashlib
import math
import os
import re
import struct
import subprocess
import sys
import tempfile

S_HDR = struct.Struct('<iid')
S_Q = struct.Struct('<q')
ID_LIMIT = 1 << 62          # plausible particle id upper bound
ID_SMALL = 10 ** 9          # real (non-artificial) ids are below this
# c.m. members of hierarchical groups are recorded with negated ids
ID_OK = lambda v: -ID_LIMIT < v < ID_LIMIT
MIN_RECOVER_RUN = 8         # shortest shifted run worth recovering


def usage():
    print(__doc__)
    print("Options:")
    print("  -h(--help)            Display help information.")
    print("  -o(--output) FILE     Output file (default: <input>.dedup).")
    print("  --in-place            Replace the input file with the result after verified write.")
    print("  --gap X               Episode time gap threshold in file time units (default 0.1).")
    print("                        Same-pair records closer than this collapse to first+last; 0 disables.")
    print("  --n-proc N            Number of parallel scan processes (default 1).")
    print("  --n-members N         Group member count N (default: from .nN filename suffix or first record).")
    print("  --interrupt-mode M    none|base|merger|bse|bseEmp|mobse|dsm|sevn (default: auto-detect).")
    print("  --external-mode M     Accepted for parity; does not affect group layout.")
    print("  --spin-1d             Use 1D spin column (outputs produced before Dec 2024).")
    print("  --use-mpfrc           Add high-precision position columns.")
    print("  --float-type T        float64 (default) or float32.")
    print("  --record-size B       Manual record size override (skips layout detection).")
    print("  --id-offset-1 B       Manual byte offset of member-1 id within a record.")
    print("  --id-offset-2 B       Manual byte offset of member-2 id within a record.")
    print("  --no-recover          Do not attempt shifted-region recovery; drop invalid bytes.")
    print("  --report FILE         Write the anomalous-pair report to FILE.")
    print("  --work-dir DIR        Directory for intermediate files (default: beside output).")


# ---------------- layout detection ----------------

def header_ok(buf, off, n_members):
    t, n, tm = S_HDR.unpack_from(buf, off)
    if n != n_members or t not in (0, 1):
        return False
    if not math.isfinite(tm) or tm < 0.0 or tm > 1e30:
        return False
    for x in struct.unpack_from('<6d', buf, off + 16):
        if not math.isfinite(x) or abs(x) > 1e30:
            return False
    return True


def verify_stride(path, R, n_members, count=200):
    if R < 64 or R % 8:
        return False
    with open(path, 'rb') as f:
        buf = f.read((count + 2) * R)
    good = 0
    for k in range(count):
        off = k * R
        if off + R > len(buf):
            break
        if header_ok(buf, off, n_members):
            good += 1
        elif good > 0:
            break
    return good >= min(count, 50)


def layout_from_petar(n_members, interrupt_mode, spin_3d, use_mpfrc, float_type):
    import petar
    import numpy  # noqa: F401  (collectDtype returns dt spec; numpy does the packing)
    ft = numpy.float32 if float_type == 'float32' else numpy.float64
    kwargs = dict(interrupt_mode=interrupt_mode, use_mpfrc=use_mpfrc,
                  spin_3d=spin_3d, float_type=ft)
    gi = petar.GroupInfo(N=n_members, **kwargs)
    dt = numpy.dtype(gi.collectDtype())
    layout = {'R': dt.itemsize}
    for name in ('bin0.p1.id', 'bin0.p2.id', 'bin0.semi', 'bin0.ecc',
                 'bin0.p1.mass', 'bin0.p2.mass', 'bin0.p1.star.type', 'bin0.p2.star.type'):
        if name in dt.fields:
            layout[name.split('.', 1)[1].replace('.', '_')] = dt.fields[name][1]
    return layout


def layout_empirical(path, n_members, samples=500):
    with open(path, 'rb') as f:
        head = f.read(min(os.path.getsize(path), 4 * 1024 * 1024))
    # record size: longest run of valid headers at a fixed stride
    best = None
    for R in range(64 + (n_members - 1) * 368, min(len(head) // 100, 16384) + 1, 8):
        run = 0
        for k in range(300):
            off = k * R
            if off + R > len(head):
                break
            if header_ok(head, off, n_members):
                run += 1
            else:
                break
        if run >= 50 and (best is None or R < best):
            best = R
            break
    if best is None:
        raise ValueError('cannot detect record size; specify --record-size/--id-offset-*')
    R = best
    P = ((R - 64) // (n_members - 1) - 144) // 2
    # member id offsets: int64 fields with many small distinct values
    cands = list(range(64 + 144, 64 + 144 + P - 8, 8))
    vals = {c: [] for c in cands}
    for k in range(min(samples, len(head) // R)):
        if not header_ok(head, k * R, n_members):
            continue  # skip grid positions inside shifted/garbage zones
        for c in cands:
            if k * R + c + P + 8 > len(head):
                continue
            v1 = struct.unpack_from('<q', head, k * R + c)[0]
            v2 = struct.unpack_from('<q', head, k * R + c + P)[0]
            vals[c].append(v1)
            vals[c].append(v2)
    for c in cands:
        vs = vals[c]
        if not vs:
            continue
        distinct = len(set(vs))
        if distinct < 20 or 0 in vs:
            continue
        if sum(1 for v in vs if abs(v) < ID_SMALL) / len(vs) < 0.9:
            continue
        if all(ID_OK(v) for v in vs):
            return {'R': R, 'p1_id': c, 'p2_id': c + P,
                    'semi': 64, 'ecc': 72, 'p1_mass': 64 + 144, 'p2_mass': 64 + 144 + P}
    raise ValueError('cannot detect member-id offsets; specify --id-offset-1/--id-offset-2')


def detect_layout(path, n_members, opt):
    if opt['record_size'] and opt['id1'] is not None and opt['id2'] is not None:
        R = opt['record_size']
        if not verify_stride(path, R, n_members):
            print('warning: manual record size %d does not validate on file header' % R)
        return {'R': R, 'p1_id': opt['id1'], 'p2_id': opt['id2'],
                'semi': 64, 'ecc': 72, 'p1_mass': 64 + 144}
    if opt['record_size'] is None and opt['interrupt_mode'] is not None:
        try:
            lay = layout_from_petar(n_members, opt['interrupt_mode'], not opt['spin_1d'],
                                    opt['use_mpfrc'], opt['float_type'])
            if verify_stride(path, lay['R'], n_members):
                lay['source'] = 'petar'
                return lay
            print('warning: petar layout for interrupt_mode=%s (R=%d) does not validate'
                  % (opt['interrupt_mode'], lay['R']))
        except ImportError:
            print('note: petar/numpy not importable; using empirical layout detection')
    if opt['record_size'] is None:
        candidates = ['none', 'base', 'merger', 'bse', 'bseEmp', 'mobse', 'dsm', 'sevn']
        if opt['interrupt_mode'] is not None and opt['interrupt_mode'] in candidates:
            candidates = [opt['interrupt_mode']]
        for mode in candidates:
            try:
                lay = layout_from_petar(n_members, mode, not opt['spin_1d'],
                                        opt['use_mpfrc'], opt['float_type'])
            except ImportError:
                break
            except Exception:
                continue
            if verify_stride(path, lay['R'], n_members):
                lay['source'] = 'petar(auto,%s)' % mode
                return lay
    lay = layout_empirical(path, n_members)
    if not verify_stride(path, lay['R'], n_members):
        raise ValueError('empirical record size %d does not validate' % lay['R'])
    lay['source'] = 'empirical'
    return lay


# ---------------- stripe worker (also used as child process) ----------------

def run_worker(wid, start, end, path, outdir, lay, n_members, gap, recover):
    R = lay['R']
    id1off, id2off = lay['p1_id'], lay['p2_id']
    size = os.path.getsize(path)
    cstart = max(0, start - R)
    cend = min(size, end + R)  # partial record at tail handled by loop bounds
    CHUNK = (16 * 1024 * 1024 // R) * R

    def rec_ok(buf, o):
        if not header_ok(buf, o, n_members):
            return False
        a = S_Q.unpack_from(buf, o + id1off)[0]
        b = S_Q.unpack_from(buf, o + id2off)[0]
        return ID_OK(a) and ID_OK(b)

    fout = open(os.path.join(outdir, 'keep_%s.bin' % wid), 'wb')
    fz = open(os.path.join(outdir, 'zones_%s.txt' % wid), 'w')
    fh2 = open(path, 'rb')
    outbuf = bytearray()
    total = valid_cnt = kept = recovered = 0
    state = {}       # pair -> [t_last, rec_last]; episode first is emitted at open
    zone_start = None

    def emit(bs):
        outbuf.extend(bs)
        if len(outbuf) >= 8 * 1024 * 1024:
            fout.write(outbuf)
            del outbuf[:]

    def feed(buf, o, off, inr):
        nonlocal valid_cnt, kept
        valid_cnt += 1
        tm = S_HDR.unpack_from(buf, o)[2]
        a = S_Q.unpack_from(buf, o + id1off)[0]
        b = S_Q.unpack_from(buf, o + id2off)[0]
        p = (a, b) if a <= b else (b, a)
        st = state.get(p)
        if gap <= 0 or st is None or tm - st[0] > gap:
            if st is not None and st[1] is not None:
                emit(st[1])
                kept += 1
            if inr:
                state[p] = [tm, None]
                emit(buf[o:o + R])
                kept += 1
            else:
                state[p] = [tm, None]
        else:
            st[0] = tm
            if inr:
                st[1] = bytes(buf[o:o + R])

    def finish_zone(z0, z1):
        # recover phase-shifted valid records inside [z0,z1); feed them in order
        nonlocal recovered, total
        if not recover or z1 - z0 < MIN_RECOVER_RUN * R:
            fz.write('%d %d\n' % (z0, z1))
            return
        pos = z0
        nrec = 0
        while z1 - pos >= MIN_RECOVER_RUN * R:
            fh2.seek(pos)
            win = fh2.read(min((MIN_RECOVER_RUN + 2) * R, z1 - pos))
            phase = None
            # fragment length is arbitrary: scan every byte phase
            for ph in range(0, len(win) - MIN_RECOVER_RUN * R + 1):
                run = 0
                q = ph
                while q + R <= len(win) and rec_ok(win, q):
                    run += 1
                    q += R
                if run >= MIN_RECOVER_RUN:
                    phase = ph
                    break
            if phase is None:
                # slide the scan window forward; stop when it cannot advance
                adv = len(win) - MIN_RECOVER_RUN * R + 1
                if adv <= 0:
                    break
                pos += adv
                continue
            seg_start = pos + phase
            fh2.seek(seg_start)
            remaining = z1 - seg_start
            q = seg_start
            seg_nrec = 0
            while remaining > 0:
                b = fh2.read(min(CHUNK, remaining))
                if not b:
                    break
                remaining -= len(b)
                o = 0
                while o + R <= len(b):
                    if not rec_ok(b, o):
                        break
                    feed(b, o, q, start <= q < end)
                    if start <= q < end:
                        recovered += 1
                        total += 1
                    seg_nrec += 1
                    q += R
                    o += R
                # inner loop stopped: either bad record (segment ends at q)
                # or a chunk tail shorter than R (zone tail) -> stop either way
                if o < len(b):
                    break
            nrec += seg_nrec
            pos = q
            fz.write('#recovered %d %d run %d phase_shift %d\n'
                     % (seg_start, q, seg_nrec, phase))
        fz.write('%d %d\n' % (z0, z1))
        if nrec:
            print('worker %s: recovered %d shifted records in zone [%d,%d)' % (wid, nrec, z0, z1))

    fh = open(path, 'rb')
    pos = cstart
    while pos < cend:
        fh.seek(pos)
        buf = fh.read(min(CHUNK, cend - pos))
        try:
            os.posix_fadvise(fh.fileno(), pos, len(buf), os.POSIX_FADV_DONTNEED)
        except OSError:
            pass
        for o in range(0, len(buf) - R + 1, R):
            off = pos + o
            inr = start <= off < end
            if inr:
                total += 1
            if rec_ok(buf, o):
                if zone_start is not None:
                    finish_zone(zone_start, off)
                    zone_start = None
                feed(buf, o, off, inr)
            elif zone_start is None:
                zone_start = off
        pos += len(buf)

    if zone_start is not None:
        finish_zone(zone_start, cend)
    for st in state.values():
        if st[1] is not None:
            emit(st[1])
            kept += 1
    fout.write(outbuf)
    fout.close()
    fz.close()
    fh.close()
    fh2.close()
    with open(os.path.join(outdir, 'status_%s.txt' % wid), 'w') as f:
        f.write('%d %d %d %d %d\n' % (total, valid_cnt, kept, total - kept, recovered))
    print('worker %s: total=%d valid=%d kept=%d recovered=%d' % (wid, total, valid_cnt, kept, recovered))


def spawn_workers(path, outdir, lay, n_members, gap, recover, n_proc):
    size = os.path.getsize(path)
    main_end = size // lay['R'] * lay['R']
    bounds = [int(main_end * i / n_proc) // lay['R'] * lay['R'] for i in range(n_proc + 1)]
    # last stripe extends to the raw file tail so shifted records beyond the
    # main grid (e.g. after a crash fragment) are scanned and recovered
    bounds[-1] = size
    args = [str(lay['R']), str(lay['p1_id']), str(lay['p2_id']), str(n_members), str(gap),
            '1' if recover else '0']
    procs = []
    for i in range(n_proc):
        cmd = [sys.executable, os.path.abspath(__file__), '--worker', str(i),
               str(bounds[i]), str(bounds[i + 1]), path, outdir] + args
        log = open(os.path.join(outdir, 'wk_%d.log' % i), 'w')
        procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT))
        log.close()
    fail = 0
    for p in procs:
        if p.wait() != 0:
            fail += 1
    if fail:
        raise RuntimeError('%d worker(s) failed; see wk_*.log in %s' % (fail, outdir))
    return main_end


# ---------------- join + final stitch pass ----------------

def stitch_pass(joined, out_path, lay, n_members, gap):
    R = lay['R']
    id1off, id2off = lay['p1_id'], lay['p2_id']
    size = os.path.getsize(joined)
    state = {}
    episodes = []   # (t_first, off_first, off_last)
    fh = open(joined, 'rb')
    pos = 0
    CHUNK = (16 * 1024 * 1024 // R) * R
    total = 0
    while pos < size:
        fh.seek(pos)
        buf = fh.read(min(CHUNK, size - pos))
        try:
            os.posix_fadvise(fh.fileno(), pos, len(buf), os.POSIX_FADV_DONTNEED)
        except OSError:
            pass
        for o in range(0, len(buf) - R + 1, R):
            tm = S_HDR.unpack_from(buf, o)[2]
            a = S_Q.unpack_from(buf, o + id1off)[0]
            b = S_Q.unpack_from(buf, o + id2off)[0]
            p = (a, b) if a <= b else (b, a)
            total += 1
            # state: [off_first, t_first, t_last, off_last]
            st = state.get(p)
            if st is None:
                state[p] = [pos + o, tm, tm, pos + o]
            elif gap <= 0 or tm - st[2] > gap:
                episodes.append((st[1], st[0], st[3]))
                state[p] = [pos + o, tm, tm, pos + o]
            else:
                st[2] = tm
                st[3] = pos + o
        pos += len(buf)
    fh.close()
    for st in state.values():
        episodes.append((st[1], st[0], st[3]))

    need = sorted(set([e[1] for e in episodes]) | set([e[2] for e in episodes]))
    recs = {}
    idx = 0
    pos = 0
    fh = open(joined, 'rb')
    while pos < size and idx < len(need):
        if need[idx] < pos:
            pos = (need[idx] // CHUNK) * CHUNK
        fh.seek(pos)
        buf = fh.read(min(CHUNK, size - pos))
        if not buf:
            break
        base = pos
        pos += len(buf)
        while idx < len(need) and base <= need[idx] < base + len(buf):
            recs[need[idx]] = buf[need[idx] - base:need[idx] - base + R]
            idx += 1
    fh.close()

    episodes.sort()
    h = hashlib.sha256()
    outbuf = bytearray()
    nkept = 0
    fout = open(out_path, 'wb')
    for _, off_a, off_b in episodes:
        ra = recs[off_a]
        pair_out = (ra,) if off_b == off_a else (ra, recs[off_b])
        for rec in pair_out:
            outbuf.extend(rec)
            h.update(rec)
            nkept += 1
            if len(outbuf) >= 8 * 1024 * 1024:
                fout.write(outbuf)
                del outbuf[:]
    fout.write(outbuf)
    fout.close()
    return total, nkept, h.hexdigest()


# ---------------- report ----------------

def write_report(path, out, lay, n_members, top=25):
    R = lay['R']
    info = {}
    fh = open(path, 'rb')
    pos = 0
    size = os.path.getsize(path)
    while pos < size:
        buf = fh.read(min(64 * 1024 * 1024, size - pos))
        for o in range(0, len(buf) - R + 1, R):
            tm = S_HDR.unpack_from(buf, o)[2]
            a = S_Q.unpack_from(buf, o + lay['p1_id'])[0]
            b = S_Q.unpack_from(buf, o + lay['p2_id'])[0]
            p = (a, b) if a <= b else (b, a)
            e = info.get(p)
            if e is None:
                info[p] = [1, tm, pos + o, tm]
            else:
                e[0] += 1
                e[3] = tm
        pos += len(buf)
    ranked = sorted(info.items(), key=lambda x: -x[1][0])
    out.write('group file: %s (N=%d)\ndistinct pairs: %d\n\n' % (path, n_members, len(info)))
    out.write('Top pairs by record count (oscillation/anomaly candidates):\n')
    out.write('%-26s %-8s %-14s %-14s\n' % ('pair', 'nrec', 't_first', 't_last'))
    for p, e in ranked[:top]:
        out.write('%-26s %-8d %-14.6f %-14.6f\n' % (str(p), e[0], e[1], e[3]))
    if 'semi' in lay:
        out.write('\nPhysics of the first record of top pairs:\n')
        for p, e in ranked[:min(top, 15)]:
            fh.seek(e[2])
            rec = fh.read(R)
            semi, ecc = struct.unpack_from('<2d', rec, lay['semi'])
            m1 = struct.unpack_from('<d', rec, lay['p1_mass'])[0]
            m2 = struct.unpack_from('<d', rec, lay['p2_mass'])[0] if 'p2_mass' in lay else float('nan')
            extra = ''
            if 'p1_star_type' in lay:
                t1 = S_Q.unpack_from(rec, lay['p1_star_type'])[0]
                t2 = S_Q.unpack_from(rec, lay['p2_star_type'])[0]
                extra = ' ssetype=(%d,%d)' % (t1, t2)
            out.write('pair %s: nrec=%d t_first=%.6f semi=%.6g ecc=%.4f m=(%.3f,%.3f)%s\n'
                      % (p, e[0], e[1], semi, ecc, m1, m2, extra))
    fh.close()
    return len(info)


# ---------------- main ----------------

def main():
    opt = {'output': None, 'in_place': False, 'gap': 0.1, 'n_proc': 1, 'n_members': None,
           'interrupt_mode': None, 'external_mode': None, 'spin_1d': False, 'use_mpfrc': False,
           'float_type': 'float64', 'record_size': None, 'id1': None, 'id2': None,
           'recover': True, 'report': None, 'work_dir': None}
    longargs = ['help', 'output=', 'in-place', 'gap=', 'n-proc=', 'n-members=', 'interrupt-mode=',
                'external-mode=', 'spin-1d', 'use-mpfrc', 'float-type=', 'record-size=',
                'id-offset-1=', 'id-offset-2=', 'no-recover', 'report=', 'work-dir=']
    try:
        opts, remainder = getopt.getopt(sys.argv[1:], 'ho:', longargs)
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for k, v in opts:
        if k in ('-h', '--help'):
            usage()
            sys.exit(0)
        elif k in ('-o', '--output'):
            opt['output'] = v
        elif k == '--in-place':
            opt['in_place'] = True
        elif k == '--gap':
            opt['gap'] = float(v)
        elif k == '--n-proc':
            opt['n_proc'] = int(v)
        elif k == '--n-members':
            opt['n_members'] = int(v)
        elif k == '--interrupt-mode':
            opt['interrupt_mode'] = v
        elif k == '--external-mode':
            opt['external_mode'] = v
        elif k == '--spin-1d':
            opt['spin_1d'] = True
        elif k == '--use-mpfrc':
            opt['use_mpfrc'] = True
        elif k == '--float-type':
            opt['float_type'] = v
        elif k == '--record-size':
            opt['record_size'] = int(v)
        elif k == '--id-offset-1':
            opt['id1'] = int(v)
        elif k == '--id-offset-2':
            opt['id2'] = int(v)
        elif k == '--no-recover':
            opt['recover'] = False
        elif k == '--report':
            opt['report'] = v
        elif k == '--work-dir':
            opt['work_dir'] = v
    if len(remainder) != 1:
        usage()
        sys.exit(1)
    path = remainder[0]
    if not os.path.exists(path):
        print('input file not found: %s' % path)
        sys.exit(1)

    n_members = opt['n_members']
    if n_members is None:
        m = re.search(r'\.n(\d+)$', path)
        if m:
            n_members = int(m.group(1))
        else:
            with open(path, 'rb') as f:
                head = f.read(16)
            n_members = S_HDR.unpack(head)[1] if len(head) >= 8 else 0
        if n_members < 2:
            print('cannot determine member count; use --n-members')
            sys.exit(1)

    lay = detect_layout(path, n_members, opt)
    print('layout: R=%d id1@%d id2@%d (source: %s, N=%d)'
          % (lay['R'], lay['p1_id'], lay['p2_id'], lay.get('source', 'manual'), n_members))

    out_path = opt['output'] or (path + '.dedup')
    work_dir = opt['work_dir'] or (os.path.dirname(os.path.abspath(out_path)))
    os.makedirs(work_dir, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix='group_dedup_', dir=work_dir)

    try:
        spawn_workers(path, work_dir, lay, n_members, opt['gap'], opt['recover'], opt['n_proc'])
        total = valid = kept = dropped = recovered = 0
        joined = os.path.join(work_dir, 'joined.bin')
        with open(joined, 'wb') as fo:
            for i in range(opt['n_proc']):
                with open(os.path.join(work_dir, 'keep_%d.bin' % i), 'rb') as fi:
                    while True:
                        b = fi.read(64 * 1024 * 1024)
                        if not b:
                            break
                        fo.write(b)
                with open(os.path.join(work_dir, 'status_%d.txt' % i)) as fs:
                    t, v, k, d, r = map(int, fs.read().split())
                total += t
                valid += v
                kept += k
                dropped += d
                recovered += r
        print('scan: total=%d valid=%d invalid-zone records dropped=%d recovered(shifted)=%d'
              % (total, valid, total - valid, recovered))
        for i in range(opt['n_proc']):
            zp = os.path.join(work_dir, 'zones_%d.txt' % i)
            if os.path.exists(zp):
                with open(zp) as f:
                    lines = [ln for ln in f if ln.strip()]
                zones = [ln for ln in lines if not ln.startswith('#recovered')]
                recs = [ln for ln in lines if ln.startswith('#recovered')]
                if zones:
                    print('  stripe %d: %d invalid zone(s), %d recovered region(s)' % (i, len(zones), len(recs)))

        # stitch to a fixpoint: time-inverted episodes (restart islands with
        # backward time jumps) can re-chain after sorting, so iterate until
        # the record count stops decreasing
        cur = joined
        nkept = None
        for it in range(10):
            nxt = out_path + ('.part' if it == 0 else '.part%d' % it)
            jtotal, nkept_new, sha = stitch_pass(cur, nxt, lay, n_members, opt['gap'])
            if cur != joined:
                os.remove(cur)
            if it:
                print('stitch pass %d: records=%d' % (it + 1, nkept_new))
            if nkept_new == nkept:
                nkept = nkept_new
                os.replace(nxt, out_path + '.part')
                break
            nkept = nkept_new
            cur = nxt
        else:
            print('warning: stitch did not reach a fixpoint in 10 passes')
            os.replace(cur, out_path + '.part')
        # strict verification of the final output
        bad = 0
        with open(out_path + '.part', 'rb') as f:
            while True:
                b = f.read(64 * 1024 * 1024)
                if not b:
                    break
                for o in range(0, len(b) - lay['R'] + 1, lay['R']):
                    if not header_ok(b, o, n_members):
                        bad += 1
                        continue
                    a = S_Q.unpack_from(b, o + lay['p1_id'])[0]
                    c = S_Q.unpack_from(b, o + lay['p2_id'])[0]
                    if not (ID_OK(a) and ID_OK(c)):
                        bad += 1
        if bad:
            print('warning: %d record(s) in output failed strict validation' % bad)
        os.replace(out_path + '.part', out_path)
        if opt['in_place']:
            os.replace(out_path, path)
            out_path = path
        osize = os.path.getsize(out_path)
        print('output: %s records=%d size=%d (%.4f GB) sha256=%s'
              % (out_path, nkept, osize, osize / 1e9, sha))
        if opt['report']:
            with open(opt['report'], 'w') as f:
                write_report(out_path, f, lay, n_members)
            print('report written: %s' % opt['report'])
    finally:
        for name in os.listdir(work_dir):
            try:
                os.remove(os.path.join(work_dir, name))
            except OSError:
                pass
        try:
            os.rmdir(work_dir)
        except OSError:
            pass


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        _wid, _start, _end, _path, _outdir = sys.argv[2:7]
        _R, _id1, _id2, _n, _gap, _rec = sys.argv[7:13]
        run_worker(_wid, int(_start), int(_end), _path, _outdir,
                   {'R': int(_R), 'p1_id': int(_id1), 'p2_id': int(_id2)},
                   int(_n), float(_gap), _rec == '1')
    else:
        main()
