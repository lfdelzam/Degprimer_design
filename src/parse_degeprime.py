#!/usr/bin/env python3

import os
import argparse
import string

usage = 'parse_degeprime.py -i -o -c -g -s'
description = 'This program selects primers suggested by degeprime above user-defined coverage and GC content'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i', dest='i', help='input file, degeprime output', required=True)
parser.add_argument(
    '-o',
    dest='o',
    help='output file, list of selected primers',
    required=True)
parser.add_argument(
    '-d',
    dest='d',
    help='output directory, default (current directory)',
    default=os.getcwd())
parser.add_argument(
    '-c',
    dest='c',
    help='coverage threshold, fraction, default 0.98',
    default=0.98)
parser.add_argument(
    '-g',
    dest='g',
    help='gc content threshold, in percentage, default 30',
    default=30)
# parser.add_argument(
#    '-t', dest='t', help='tm threshold, in C, default 30', default=30)
parser.add_argument(
    '-s', dest='s', help='suffix for figures', required=True)

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

# creating output directory


if not os.path.exists(args.d):
    os.makedirs(args.d)
if not os.path.exists(os.path.join(args.d, "Figures")):
    os.makedirs(os.path.join(args.d, "Figures"))

outfile = os.path.join(args.d, args.o)
with open(args.i, "r") as fin, open(outfile, "w") as fout:
    for line in fin:
        line = line.rstrip()
        if line.startswith("Pos"):
            line += "\tminGC\tCOV\tREVCOMP\tObservations"
            print(line, file=fout)
        else:
            line = line.split()
            # Pos     TotalSeq        UniqueMers      Entropy PrimerDeg
            # PrimerMatching  PrimerSeq
            cov = round(int(line[5]) / int(line[1]), 3)
            seq = line[6]
            all_primers = undegenerating(seq)
            GCrange = []
            for p in all_primers.values():
                GCrange.append(gc_cont(p))

            GC = min(GCrange)

            if (GC >= float(args.g) and cov >= float(args.c)):
                line.append(GC)
                line.append(cov)
                line.append(revcomp(seq))
                text = "GC range % (" + str(min(GCrange)) + \
                    "," + str(max(GCrange)) + ")"
                line.append(text)
                print("\t".join(map(str, line)), file=fout)
