#!/usr/bin/env python3
"""release.note maintenance tool.

Subcommands:
  check                 verify the top list: numbering continuity, sort order,
                        entry format (verb prefix, <=150 chars, [major] tag,
                        version prefix); print a summary
  insert "<entry>" ...  insert one or more entries ("(V) Text..." or
                        "N. (V) Text...") into the top list at their sorted
                        position and renumber; use --dry-run to preview
  sort                  re-sort the whole top list by version and renumber

Sorting key is hardened for BOTH "N. (V) ..." and bare "(V) ..." forms —
a bare-prefix entry must never fall to the tail (the recurring pit this
tool exists to prevent).

Versionless entries keep their position via the PLACEMENT map below (distinctive
substring -> sort number). When a versionless entry gets its version annotated,
remove it from the map.
"""
import re
import sys

PATH_DEFAULT = "release.note"
TOP_HDR = "**experiment update referring to master**"
ISS_HDR = "**Existing issues"

# versionless entries: distinctive substring -> sort number (era position)
PLACEMENT = {
    "point sizes": 1270,
    "nested OMP": 1372,
}

VERBS = ("Add", "Fix", "Change", "Improve", "Remove")


def parse(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    top_start = lines.index(TOP_HDR) + 1
    iss = next(i for i in range(top_start, len(lines)) if lines[i].startswith(ISS_HDR))
    top = [l for l in lines[top_start:iss] if re.match(r"^\d+\. ", l)]
    return lines, top_start, iss, top


def sort_key(line):
    # hardened: matches both "N. (V) ..." and bare "(V) ..."
    m = re.match(r"^(?:\d+\. )?\((\d+)", line)
    if m:
        return (int(m.group(1)), 0)
    for sub, num in PLACEMENT.items():
        if sub in line:
            return (num, 0)
    return (9999, 9)


def renumber(entries):
    return [re.sub(r"^\d+\. ", f"{i+1}. ", l) for i, l in enumerate(entries)]


def check(top):
    problems = []
    prev, prevn = 0, 0
    for i, l in enumerate(top, 1):
        m = re.match(r"^(\d+)\. ", l)
        if not m:
            problems.append(f"line {i}: no entry number: {l[:60]}")
            continue
        n = int(m.group(1))
        if n != prevn + 1:
            problems.append(f"line {i}: numbering gap/jump: {prevn} -> {n}")
        prevn = n
        vm = re.match(r"^\d+\. \(([^)]+)\) ", l)
        if not vm:
            problems.append(f"#{n}: no version prefix: {l[:60]}")
        else:
            v = int(re.match(r"\d+", vm.group(1)).group(0))
            if v < prev:
                problems.append(f"#{n}: sort break: {prev} -> {v}: {l[:60]}")
            prev = v
            # length is enforced only for entries after the rule took effect
            # (2026-10-24, ~1848e); older ones are grandfathered
            if len(l) > 150 and v >= 1848:
                problems.append(f"#{n}: over 150 chars ({len(l)})")
        vm2 = re.match(r"^\d+\. (?:\([^)]*\) )?(\w+)", l)
        if not vm2 or vm2.group(1) not in VERBS:
            problems.append(f"#{n}: non-standard verb: {l[:60]}")
        if l.startswith("*"):
            problems.append(f"#{n}: leading star marker (use trailing [major])")
    fallback = [l for l in top if sort_key(l)[0] == 9999]
    for l in fallback:
        problems.append(f"unsorted-fallback (no version, not in PLACEMENT): {l[:60]}")
    print(f"entries: {len(top)}; problems: {len(problems)}")
    for p in problems:
        print("  " + p)
    return len(problems) == 0


def insert(top, new_entries, dry):
    all_e = top + [re.sub(r"^\d+\. ", "", e) for e in new_entries]
    all_e.sort(key=sort_key)
    return renumber(all_e)


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return 1
    cmd, rest = args[0], args[1:]
    dry = "--dry-run" in rest
    rest = [a for a in rest if a != "--dry-run"]
    path = rest[0] if rest and (cmd != "insert" or rest[0].endswith(".note")) else PATH_DEFAULT
    if cmd == "insert":
        entries = [a for a in rest if not a.endswith(".note")]
        lines, ts, iss, top = parse(path)
        new = insert(top, entries, dry)
        if dry:
            print("\n".join(new))
        else:
            open(path, "w", encoding="utf-8").write(
                "\n".join(lines[:ts] + new + [""] + lines[iss:]))
            print(f"inserted {len(entries)}; total {len(new)}")
    elif cmd == "sort":
        lines, ts, iss, top = parse(path)
        new = renumber(sorted(top, key=sort_key))
        open(path, "w", encoding="utf-8").write(
            "\n".join(lines[:ts] + new + [""] + lines[iss:]))
        print(f"sorted {len(new)} entries")
    elif cmd == "check":
        sys.exit(0 if check(parse(path)[3]) else 1)
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
