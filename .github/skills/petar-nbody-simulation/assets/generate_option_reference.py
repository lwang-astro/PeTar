#!/usr/bin/env python3
"""
Generate a complete CLI option reference from PeTar source headers.

Parses header files for `IOParams` declarations and their initializers
to produce a catalog of ALL possible PeTar command-line options,
grouped by configure flag. No binary compilation required.

Usage: python3 generate_option_reference.py
"""

import re
import os
import sys
from collections import OrderedDict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# assets/ → petar-nbody-simulation/ → skills/ → .github/ → PeTar root
REPO_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", ".."))
ASSET_DIR = SCRIPT_DIR

# ---------------------------------------------------------------------------
HEADER_FILES = OrderedDict([
    ("src/petar.hpp", {
        "label": "Core Simulation",
        "guard_group": "core",
        "conditional_guards": {
            "PARTICLE_SIMULATOR_MPI_PARALLEL": "MPI build",
            "KDKDK_4TH": "4th-order KDKDK step mode",
            "USE_GPU": "--enable-gpu",
            "USE_FUGAKU": "--enable-fugaku",
            "PROFILE": "--enable-profile",
            "GALPY": "--with-external=galpy",
            "AGAMA": "--with-external=agama",
            "BSE_BASE": "--with-interrupt",
            "DISK_STAR_MERGER": "--with-interrupt=dsm",
            "MPI_DEBUG": "--with-debug=g (MPI debug)",
        },
    }),
    ("src/hard.hpp", {
        "label": "Hard Integrator",
        "guard_group": "core",
        "conditional_guards": {
            "HARD_CHECK_ENERGY": "development debug",
            "ORBIT_SAMPLING": "--enable-orbit-sampling",
            "STELLAR_EVOLUTION": "--with-interrupt",
            "BSE_BASE": "--with-interrupt=base",
            "ADJUST_GROUP_PRINT": "--enable-adjust-group-print",
            "HERMITE_PN": "--with-pn (Hermite)",
            "SDAR_PN": "--with-pn (SDAR)",
            "HERMITE_ONLY_CALC_NEIGHBOR_FORCE": "--enable-hermite-only-calc-neighbor-force",
        },
    }),
    ("bse-interface/bse_interface.h", {
        "label": "BSE Stellar Evolution",
        "guard_group": "bse",
        "conditional_guards": {
            "BSEBBF": "--with-interrupt=bse",
            "BSEEMP": "--with-interrupt=bseEmp",
            "MOBSE": "--with-interrupt=mobse",
        },
    }),
    ("src/disk_star_merger.hpp", {
        "label": "DSM (Disk Star Merger)",
        "guard_group": "dsm",
    }),
    ("src/gas_drag.hpp", {
        "label": "Gas Drag",
        "guard_group": "gasdrag",
    }),
    ("src/external_hard.hpp", {
        "label": "External Hard",
        "guard_group": "external-hard",
    }),
    ("galpy-interface/galpy_interface.h", {
        "label": "Galpy External Potential",
        "guard_group": "galpy",
    }),
    ("agama-interface/agama_interface.h", {
        "label": "Agama External Potential",
        "guard_group": "agama",
    }),
    ("parallel-random/rand_io.hpp", {
        "label": "Random Number Generator",
        "guard_group": "core",
    }),
])

GUARD_GROUP_TABLE = OrderedDict([
    ("core",          "Core (always available)"),
    ("bse",           "BSE Interrupt (`--with-interrupt=bse|mobse|bseEmp`)"),
    ("dsm",           "DSM Interrupt (`--with-interrupt=dsm`)"),
    ("galpy",         "Galpy External (`--with-external=galpy`)"),
    ("agama",         "Agama External (`--with-external=agama`)"),
    ("gasdrag",       "Gas Drag (`--with-external-hard=gasdrag`)"),

    ("external-hard", "External Hard (`--with-external-hard`)"),
])


def read_file(path):
    full = os.path.join(REPO_ROOT, path)
    if not os.path.isfile(full):
        return ""
    with open(full, "r") as f:
        return f.read()


def extract_guard_per_line(lines):
    """List of (1-based lineno, set of active guard names)."""
    result = []
    active = []
    for i, line in enumerate(lines, start=1):
        m = re.match(r'^\s*#\s*(ifdef|ifndef)\s+(\w+)', line)
        if m:
            d, name = m.group(1), m.group(2)
            active.append(name if d == "ifdef" else "!" + name)
            result.append((i, set(active)))
            continue
        if re.match(r'^\s*#\s*else', line) and active:
            last = active[-1]
            active[-1] = "!" + last if not last.startswith("!") else last[1:]
            result.append((i, set(active)))
            continue
        if re.match(r'^\s*#\s*endif', line) and active:
            active.pop()
            result.append((i, set(active)))
            continue
        result.append((i, set(active)))
    return result


def extract_options(lines, guard_per_line, source_file, module_label):
    """
    Two-pass extraction.
    Pass 1: `IOParams<TYPE> VAR_NAME;`  →  record VAR_NAME → TYPE
    Pass 2: `VAR_NAME (input_par_store, DEFAULT, "OPT", "DESC")`  →  combine
    """
    # Pass 1 — declarations
    decl_re = re.compile(r'IOParams\s*<\s*([^>]+?)\s*>\s+(\w+)\s*;')
    var_types = {}
    for i, (line_no, guards) in enumerate(guard_per_line):
        if i >= len(lines):
            break
        m = decl_re.search(lines[i])
        if m:
            var_types[m.group(2)] = {"type": m.group(1).strip(), "guards": set(guards)}

    # Pass 2 — initializers
    #   VAR_NAME (input_par_store, DEFAULT, "OPT-NAME", "DESC"
    # The description may contain commas inside quotes; we match the closing ") pattern.
    init_re = re.compile(
        r'(\w+)\s*\(\s*input_par_store\s*,\s*'
        r'([^,]*?)\s*,\s*'
        r'"([^"\\]*(?:\\.[^"\\]*)*)"\s*,\s*'
        r'"([^"\\]*(?:\\.[^"\\]*)*)"'
    )

    options = []
    for i, (line_no, guards) in enumerate(guard_per_line):
        if i >= len(lines):
            break
        m = init_re.search(lines[i])
        if not m:
            continue
        var_name = m.group(1)
        default_val = m.group(2).strip().rstrip(",").rstrip(")")
        opt_name = m.group(3).strip()
        description = m.group(4).strip()

        if not opt_name or var_name not in var_types:
            continue

        vinfo = var_types[var_name]
        options.append({
            "option_name": opt_name,
            "description": description,
            "default": default_val,
            "type": vinfo["type"],
            "line": line_no,
            "guards": vinfo["guards"] | set(guards),
            "source_file": source_file,
            "source_module": module_label,
        })
    return options


def build_guard_group_map():
    m = {}
    for fp, info in HEADER_FILES.items():
        gg = info["guard_group"]
        if "conditional_guards" in info:
            for cpp_name, expl in info["conditional_guards"].items():
                m[cpp_name] = (gg, expl)
    return m


def main():
    print("Parsing PeTar source headers for IOParams declarations...", file=sys.stderr)

    guard_group_map = build_guard_group_map()
    all_options = []

    for fp, info in HEADER_FILES.items():
        content = read_file(fp)
        if not content:
            print(f"  SKIP: {fp}", file=sys.stderr)
            continue
        lines = content.split("\n")
        gpl = extract_guard_per_line(lines)
        opts = extract_options(lines, gpl, fp, info["label"])
        all_options.extend(opts)
        print(f"  {fp}: {len(opts)} options", file=sys.stderr)

    # Group by guard_group
    groups = OrderedDict((k, {"title": v, "options": []}) for k, v in GUARD_GROUP_TABLE.items())

    # build reverse: filepath → guard_group
    file_group = {}
    for fp, info in HEADER_FILES.items():
        file_group[fp] = info["guard_group"]

    for opt in all_options:
        assigned = False
        # First try: match by individual option guards
        for g in opt["guards"]:
            gn = g.lstrip("!")
            if gn in guard_group_map:
                gg, _ = guard_group_map[gn]
                if gg in groups:
                    groups[gg]["options"].append(opt)
                    assigned = True
                    break
        # Second try: fall back to file-level group
        if not assigned:
            fg = file_group.get(opt["source_file"], "core")
            if fg in groups:
                groups[fg]["options"].append(opt)
            else:
                groups["core"]["options"].append(opt)

    # ---- Write option-reference.md ----
    from datetime import datetime, timezone
    output_path = os.path.join(ASSET_DIR, "option-reference.md")
    total = sum(len(g["options"]) for g in groups.values())
    with open(output_path, "w") as f:
        f.write("# PeTar CLI Option Reference (Source-Generated)\n\n")
        f.write("All possible PeTar command-line options across every configure variant.\n")
        f.write("Generated from source headers — no binary compilation needed.\n\n")
        f.write(f"> Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n")
        f.write(f"> Regenerate: `python3 .github/skills/petar-nbody-simulation/assets/generate_option_reference.py`\n\n")
        f.write("See [`option-matrix.md`](option-matrix.md) for currently installed binaries.\n")
        f.write("Use `<binary> -h` for exact runtime option validation.\n\n---\n\n")

        f.write("## Summary\n\n")
        f.write("| Configure Flag Group | Option Count |\n")
        f.write("|-----------------------|-------------|\n")
        for gg, gdata in groups.items():
            f.write(f"| {gdata['title']} | {len(gdata['options'])} |\n")
        f.write(f"| **Total** | **{total}** |\n\n---\n\n")

        for gg, gdata in groups.items():
            if not gdata["options"]:
                continue
            f.write(f"## {gdata['title']}\n\n")

            seen = {}
            for opt in gdata["options"]:
                if opt["option_name"] not in seen:
                    seen[opt["option_name"]] = opt

            for name in sorted(seen):
                opt = seen[name]
                f.write(f"### `--{opt['option_name']}`\n\n")
                f.write(f"- **Type**: `{opt['type']}`\n")
                f.write(f"- **Default**: `{opt['default']}`\n")
                if opt["description"]:
                    f.write(f"- **Description**: {opt['description']}\n")
                f.write(f"- **Source**: `{opt['source_file']}`:{opt['line']}\n")
                if opt["guards"]:
                    human = []
                    for g in sorted(opt["guards"]):
                        gn = g.lstrip("!")
                        neg = " (negated)" if g.startswith("!") else ""
                        if gn in guard_group_map:
                            _, expl = guard_group_map[gn]
                            human.append(f"`{gn}` → {expl}{neg}")
                        else:
                            human.append(f"`{g}`{neg}")
                    f.write(f"- **Guards**: {', '.join(human)}\n")
                f.write("\n")

    print(f"Wrote {output_path} ({total} options)", file=sys.stderr)

    # ---- option-list.txt ----
    list_path = os.path.join(ASSET_DIR, "option-list.txt")
    all_names = sorted(set(o["option_name"] for o in all_options if o["option_name"]))
    with open(list_path, "w") as f:
        for n in all_names:
            f.write(f"{n}\n")
    print(f"Wrote {list_path} ({len(all_names)} option names)", file=sys.stderr)


if __name__ == "__main__":
    main()
