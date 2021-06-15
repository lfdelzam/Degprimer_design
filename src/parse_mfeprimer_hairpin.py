#!/usr/bin/env python3

import os
import argparse
import string
import statistics

usage = 'parse_mfeprimer_hairpin.py -i -o -c -f'
description = 'This program selects primers according to hairpin formation'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i', dest='i', help='input file mfeprimer3', required=True)
parser.add_argument(
    '-f', dest='f', help='input file .fasta, primers', required=True)
parser.add_argument(
    '-o', dest='o', help='output file .fasta, selected primers', required=True)

args = parser.parse_args()

problem="Hairpin"
problem_list=problem+" List"
PCR_problem={}

with open(args.i, "r") as fin:
    for line in fin:
        line=line.rstrip()
        if line.startswith("Primer_"):
            prim=line.split()[0]
            prim=prim.split("_at_")[0]
            if not prim in PCR_problem:
                PCR_problem[prim]=0

        if line.startswith("Hairpin ") and not line.startswith("Hairpin List"):
#            Hairpin 1: Primer_145_at_481.53
            p=line.split(":")[1]
            primer=p.split("_at_")[0]
            primer=primer.replace(" ","")
            PCR_problem[primer]+=1

#print(PCR_problem)
selected_primers=set()


for primer, prob in PCR_problem.items():
    if prob == 0 :
        selected_primers.add(primer)


print("Total number of hairpin primers (excluded):", len(PCR_problem.keys())-len(selected_primers))
print("Number of selected degenerated primers:", len(selected_primers))

with open(args.f, "r") as fp, open(args.o, "w") as fout:
    copy=False
    for line in fp:
        line=line.rstrip()
        if line.startswith(">"):
            primer_name=line.split("_at_")[0][1:]
            if primer_name in selected_primers:
                header=line
                copy=True
        else:
            if copy:
                print(header, file=fout)
                print(line,file=fout)
                copy=False
