#!/usr/bin/env python3
"""Regenerate compile_commands.json for clangd from PeTar's live build configuration.

PeTar's active feature macros depend on ``./configure`` options and on the binary
family selected by ``petar.select`` (for example ``--with-se=dsm``,
``--with-ext=galpy``, ``--with-ext-hard=gasdrag``).  ``compile_commands.json`` is
a checked-in file, so it silently goes stale whenever the build is reconfigured:
clangd then greys out code guarded by macros that are actually enabled
(``STELLAR_EVOLUTION``, ``GALPY``, ``GAS_DRAG``, ``DISK_STAR_MERGER``,
``EXTERNAL_POT_IN_PTCL``, ...).

This script asks the Makefile for the *currently expanded* flag variables and
writes a ``compile_commands.json`` matching what the compiler really sees, so the
file can be refreshed with a single command after reconfiguring:

    python3 tools/gen_compile_commands.py
    python3 tools/gen_compile_commands.py --print   # inspect the recipe only

The compile command is modelled on the ``*.hard.debug`` target recipes, i.e. it
uses ``DEBUG_OPT_FLAGS`` (-g -O0) plus ``HARD_DEBFLAGS``, so that debug-only code
paths stay visible while browsing.  All macros from the real build
(``MT_FLAGS``, ``CXXFLAGS``, ``OMPFLAGS``, ``FDPSFLAGS``, ...) are included.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Makefile variables that together describe a real compilation of a PeTar source.
MAKE_VARS = (
    "PETAR_INCLUDE",
    "OPTFLAGS",
    "DEBUG_OPT_FLAGS",
    "CXXFLAGS",
    "OMPFLAGS",
    "FDPSFLAGS",
    "MT_FLAGS",
    "DEBFLAGS",
    "HARD_DEBFLAGS",
)

# Only translation units of the PeTar build itself are listed; clangd infers the
# same command for the headers they pull in.  Other trees (test/, sample/, the
# interface directories) are deliberately skipped.
SOURCE_DIRS = ("src",)
SOURCE_SUFFIXES = (".cxx", ".cc")

COMPILER = "/usr/bin/g++"
DRIVER = "/usr/bin/mpic++"  # only consulted to extract the MPI include dirs


def make_var(name: str) -> list[str]:
    """Return the expanded value of a Makefile variable as a list of words."""
    # GNU make expands $(name) at recipe time, i.e. with all conditionals and
    # += accumulations from the configured Makefile already applied.
    recipe = "petar-print-var: ; @printf '%s' '$(" + name + ")'"
    result = subprocess.run(
        ["make", "--no-print-directory", "--eval", recipe, "petar-print-var"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.exit(f"error: could not read Makefile variable {name}:\n{result.stderr}")
    return result.stdout.split()


def mpi_compile_flags() -> list[str]:
    """MPI include/link flags the mpic++ wrapper would add."""
    try:
        out = subprocess.run(
            [DRIVER, "--showme:compile"], capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"warning: cannot query {DRIVER} ({exc}); no MPI includes added", file=sys.stderr)
        return []
    return shlex.split(out)


def build_flags(variables: dict[str, list[str]]) -> list[str]:
    return normalize_includes(
        mpi_compile_flags()
        + ["-pthread"]
        + variables["PETAR_INCLUDE"]
        + variables["DEBUG_OPT_FLAGS"]
        + variables["CXXFLAGS"]
        + variables["OMPFLAGS"]
        + variables["FDPSFLAGS"]
        + variables["MT_FLAGS"]
        + variables["HARD_DEBFLAGS"]
        + ["-D", "STABLE_CHECK_DEBUG"]
    )


def find_sources() -> list[Path]:
    sources: list[Path] = []
    for directory in SOURCE_DIRS:
        base = ROOT / directory
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if path.suffix in SOURCE_SUFFIXES and path.is_file():
                sources.append(path)
    return sources


def normalize_includes(flags: list[str]) -> list[str]:
    """Collapse ``..`` segments in -I paths so the file stays readable/diffable."""
    result: list[str] = []
    index = 0
    while index < len(flags):
        flag = flags[index]
        if flag == "-I" and index + 1 < len(flags):
            result.append("-I" + os.path.normpath(flags[index + 1]))
            index += 2
            continue
        if flag.startswith("-I") and len(flag) > 2 and not flag.startswith("-I-"):
            result.append("-I" + os.path.normpath(flag[2:]))
            index += 1
            continue
        result.append(flag)
        index += 1
    return result


def macro_names(flags: list[str]) -> list[str]:
    """Extract the macros defined by the flag list (both -DNAME and -D NAME)."""
    names: list[str] = []
    index = 0
    while index < len(flags):
        flag = flags[index]
        if flag == "-D" and index + 1 < len(flags):
            names.append(flags[index + 1])
            index += 2
            continue
        if flag.startswith("-D") and len(flag) > 2:
            names.append(flag[2:])
        index += 1
    return names


def main() -> int:
    parser = argparse.ArgumentParser(description="generate compile_commands.json for clangd")
    parser.add_argument(
        "--print",
        dest="print_only",
        action="store_true",
        help="print the compile recipe instead of writing compile_commands.json",
    )
    parser.add_argument(
        "--output",
        default=str(ROOT / "compile_commands.json"),
        help="path of the file to write (default: compile_commands.json)",
    )
    args = parser.parse_args()

    variables = {name: make_var(name) for name in MAKE_VARS}
    flags = build_flags(variables)

    if args.print_only:
        print(" ".join([COMPILER, "-c"] + flags + ["src/main.cc"]))
        return 0

    entries = [
        {
            "directory": str(ROOT),
            "file": source.relative_to(ROOT).as_posix(),
            "arguments": [COMPILER, "-c"] + flags + [source.relative_to(ROOT).as_posix()],
        }
        for source in find_sources()
    ]
    if not entries:
        sys.exit("error: no translation units found under " + ", ".join(SOURCE_DIRS))

    out = Path(args.output)
    out.write_text(json.dumps(entries, indent=1) + "\n", encoding="utf-8")

    macros = sorted(set(macro_names(flags)))
    print(f"wrote {out} with {len(entries)} entries")
    print(f"  {len(macros)} macros: {' '.join(macros)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
