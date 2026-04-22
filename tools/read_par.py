#!/usr/bin/env python3

import getopt
import sys


def usage():
    print("A tool to read PeTar input parameter files and display floating-point values in decimal.")
    print("Usage: petar.read.par [options] [input parameter filename]")
    print("Options:")
    print("  -h(--help)           Display help information.")
    print("  -o(--output) [FILE]  Write converted content to FILE instead of standard output.")


def parse_float(value):
    try:
        return float(value)
    except ValueError:
        return float.fromhex(value)


def convert_line(line):
    stripped = line.strip()
    if stripped == '' or stripped.startswith('#'):
        return line

    columns = stripped.split(None, 2)
    if len(columns) != 3:
        raise ValueError("Line does not match the expected 3-column input.par format")

    value_type, key, value = columns
    if value_type == 'F':
        value = repr(parse_float(value))

    return f"{value_type} {key} {value}\n"


if __name__ == '__main__':
    output_path = None

    try:
        shortargs = 'ho:'
        longargs = ['help', 'output=']
        opts, remainder = getopt.getopt(sys.argv[1:], shortargs, longargs)

        for opt, arg in opts:
            if opt in ('-h', '--help'):
                usage()
                sys.exit(0)
            elif opt in ('-o', '--output'):
                output_path = arg
            else:
                assert False, 'unhandled option'

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    if len(remainder) != 1:
        usage()
        sys.exit(1)

    input_path = remainder[0]
    with open(input_path, 'r') as source:
        converted = ''.join(convert_line(line) for line in source)

    if output_path is None:
        sys.stdout.write(converted)
    else:
        with open(output_path, 'w') as destination:
            destination.write(converted)