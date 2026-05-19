#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
   rebuild a table with the usual 13 fields of metadata
"""

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2020/12/22"
__version__ = "$Revision: 1.0"

import sys
import argparse


# *************************************************************************** #
#                                                                             #
#                                  Functions                                  #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':
    """
    Parse arguments from command line.
    """
    parser = argparse.ArgumentParser(
        description="rebuild an occurrence table with 13 fields of metadata.")

    parser.add_argument("--mumu_table",
                        dest="mumu_table",
                        required=True,
                        help="occurrence table produced by mumu")

    parser.add_argument("--old_table",
                        dest="old_table",
                        required=True,
                        help="occurrence table pre-mumu")

    ARGS = parser.parse_args()


def parse_old_table(old_table):
    """
    Map OTUs and metadata.
    """
    separator = '\t'
    amplicons = [dict() for i in range(0, 256)]

    with open(old_table, 'r') as old_table:
        print("PROGRESS: parsing old table", file=sys.stderr)
        is_first_line = True
        for line in old_table:
            
            # assume sample columns in both tables are in the same order
            if is_first_line:
                print(line.strip(), file=sys.stdout)
                is_first_line = False
                continue

            # index metadata (assume chimera status is always "N")
            line = line.strip().split(separator)
            amplicon = line[3]
            index = int(amplicon[0:2], 16)
            amplicons[index][amplicon] = {"length": line[4],
                                          "abundance": line[5],
                                          "quality": line[8],
                                          "sequence": line[9],
                                          "identity": line[10],
                                          "taxonomy": line[11],
                                          "references": line[12]}
            
    return amplicons


def parse_mumu_table(mumu_table, amplicons):
    """
    Compute and merge metadata, make a new table.
    """
    separator = '\t'
    cloud = "NA"  # that data is lost during the mumu step
    chimera = "N"  # assume chimera status is always "N"
    with open(mumu_table, 'r') as mumu_table:
        print("PROGRESS: parsing mumu table", file=sys.stderr)
        line_counter = 0
        for line in mumu_table:
            # discard mumu's header
            if line_counter == 0:
                line_counter += 1
                continue

            # index metadata
            line = line.strip().split(separator)
            amplicon = line[0]
            index = int(amplicon[0:2], 16)

            # print the new line, add old fields
            print(line_counter,
                  sum([int(i) for i in line[1:]]),  # new total
                  cloud,
                  amplicon,
                  amplicons[index][amplicon]["length"],
                  amplicons[index][amplicon]["abundance"],
                  chimera,
                  len([e for e in line[1:] if e != "0"]),  # new spread
                  amplicons[index][amplicon]["quality"],
                  amplicons[index][amplicon]["sequence"],
                  amplicons[index][amplicon]["identity"],
                  amplicons[index][amplicon]["taxonomy"],
                  amplicons[index][amplicon]["references"],
                  '\t'.join(line[1:]),  # occurrence values
                  sep='\t')
            line_counter += 1
    
    return None


def main():
    """
    rebuild an occurrence table with 13 fields of metadata.
    """
    # capture arguments
    mumu_table = ARGS.mumu_table
    old_table = ARGS.old_table

    # parse and merge
    amplicons = parse_old_table(old_table)
    parse_mumu_table(mumu_table, amplicons)

    return None


# *************************************************************************** #
#                                                                             #
#                                     Body                                    #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':

    main()

sys.exit(0)
