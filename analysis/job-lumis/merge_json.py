#!/usr/bin/python

import sys
import json


def merge_file(filename):
    """merge_file(filename)

    filename -- the path to the concatenated json certificates file
    format:
        [ ["run_num", [lumisections]] ]

    (should be prepared from `cat *json > all.json` in advance)
    """

    with open(filename, 'r') as f:
        d = json.load(f)

    # d is [ ["run_num", [lumisections]] ]

    m = {}

    for i, record in enumerate(d):
        if len(record) != 2: raise(ValueError('Line %s is not ["run_num", [lumisecs]]' % i))
        run_num, lumies = record
        m.setdefault(run_num, []).extend(lumies)

    return m

if __name__ == '__main__':
    print("argv", sys.argv)
    if len(sys.argv) != 3:
        print("Usage:\nmerge_json.py filename out_filename")
        print("Input file format:")
        print('    [ ["run_num", [lumisections]] ]')
        print('    (should be prepared in advance from regular `cat *json > all.json`)')
        print('Output file format:')
        print('    {"run_num": [lumies]}')
        exit(1)

    filename, out_filename = sys.argv[1:3]
    print("Merging %s" % filename)
    res = merge_file(filename)
    print("Dumping to %s" % out_filename)
    with open(out_filename, 'w') as f:
        json.dump(res, f)
    print("Done.")
    exit()

