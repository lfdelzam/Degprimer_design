#!/usr/bin/env python3

import os
import argparse
import string

usage = 'Primer_screen.py -d -o -c -g'
description = 'This program selects primers suggested by degeprime above user-defined coverage and GC content'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-d',
    dest='d',
    help='input directory, degeprime output *_selected',
    required=True)
parser.add_argument(
    '-o',
    dest='o',
    help='output file, list of selected primers',
    required=True)
parser.add_argument(
    '-c',
    dest='c',
    help='coverage threshold, fraction, default 0.99',
    default=0.99)
parser.add_argument(
    '-g',
    dest='g',
    help='gc content threshold, in percentage, default 33',
    default=33)

args = parser.parse_args()

degen = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
         "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"],
         "V": ["A", "G", "C"], "D": ["A", "G", "T"],
         "H": ["A", "T", "C"], "B": ["C", "G", "T"],
         #       "I": "N" }
         "I": ["A", "G", "C", "T"], "N": ["A", "G", "C", "T"]}


def undegenerating(primer):
    total_deg = 1
    changes = []  # List of changes per position
    primers = {}  # Dictionary storing the primers
    for item in primer:
        if item in degen.keys():
            total_deg = total_deg * len(degen[item])  # Total number of primers
            changes.append(total_deg)
        else:
            changes.append(1)

    for i in range(0, total_deg):
        primers[i] = ["N" for n in range(0, len(primer))]  # initializing

    for it in range(0, len(primer)):
        deg = primer[it]
        pr = 0
        while pr < total_deg:  # for all the primers
            if deg in degen.keys():
                Nn = len(degen[deg])  # number of nucleotides
                Ntc = changes[it]  # number of changes
                # number of times the same nucleotide should be used for if
                # window change
                window = total_deg / Ntc
                for ch in range(
                        0, Nn):  # for each nucleotide in degenerated nucleotide
                    counter = 0
                    while counter < window:
                        primers[pr][it] = degen[deg][ch]
                        counter += 1
                        pr += 1
            else:
                primers[pr][it] = deg
                pr += 1
    return primers


def gc_cont(new):
    G = new.count("G")
    C = new.count("C")
    S = new.count("S")
    GC = round((G + C + S) * 100 / len(new), 2)
    return GC


def revcomp(line):
    old_chars = "ACGTRYMKSWVBDHIN"
    replace_chars = "TGCAYRKMSWBVHDNN"
    new = line.translate(line.maketrans(old_chars, replace_chars))[::-1]
    return new


List_of_files = [f for f in os.listdir(args.d) if f.endswith("selected")]
PRIM = set()
selected_lines = []
firsttime = True
# printing out selected primers
with open(args.o, "w") as fout:
    for file in List_of_files:
        with open(os.path.join(args.d, file), "r") as fin:
            for line in fin:
                line = line.rstrip()
                if line.startswith("Pos"):
                    if firsttime:
                        print(line, file=fout)
                        firsttime = False
                else:
                    line = line.split()
                    # Pos     TotalSeq        UniqueMers      Entropy PrimerDeg
                    # PrimerMatching  PrimerSeq
                    cov = float(line[8])
                    GC = float(line[7])
                    Posit = line[0]
                    sequ = line[6]
                    if (GC >= float(args.g) and cov >= float(args.c)) and (
                            line not in selected_lines) and ((Posit, sequ) not in PRIM):
                        # if GC and Cov criteria are fullfil                  #
                        # selected_lines and PRIM condition are used to avoid
                        # printing the same primer more then onece
                        selected_lines.append(line)
                        PRIM.add((Posit, sequ))

    selected_lines = sorted(
        selected_lines,
        reverse=False,
        key=lambda x: int(
            x[0]))  # order by position
    for item in selected_lines:
        print("\t".join(map(str, item)), file=fout)
