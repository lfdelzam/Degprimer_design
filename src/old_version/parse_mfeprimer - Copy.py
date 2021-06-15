#!/usr/bin/env python3

import os
import argparse
import string
import statistics

usage = 'parse_mfeprimer.py -i -o -c -t -f'
description = 'This program selects primers according to dimer or hairpin formation'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i', dest='i', help='input file mfeprimer3', required=True)
parser.add_argument(
    '-f', dest='f', help='input file .fasta, primers', required=True)
parser.add_argument(
    '-o', dest='o', help='output file .fasta, selected primers', required=True)
parser.add_argument(
    '-c', dest='c', help='Hairpin or Dimer', default='Dimer')
parser.add_argument(
    '-t', dest='t', help='Hairpin/Dimer content threshold', default="cal")
parser.add_argument(
    '-p', dest='p', help='print table option', default=False)

args = parser.parse_args()
#Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg
#                                                                      (bp)      (%)    (°C)  (kcal/mol)
#
#Primer_74_at_326.1             GGGTATGAATCCGATGG                         17   52.94   51.43    -16.94

#remover primer que hacen dimer con ellos mismos antes de remover otros primers


degen = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
         "S": ["G", "C"], "W": ["A", "T"],  "K": ["G", "T"],
         "V": ["A", "G", "C"], "D": ["A", "G", "T"],
         "H": ["A", "T", "C"], "B": ["C", "G", "T"],
         #       "I": "N" }
         "I": ["A", "G", "C", "T"], "N": ["A", "G", "C", "T"]}

#RYMKSWVBDHIN / YRKMSWBVHDNN

def degenerecy(primer):
    total_deg = 1
    for item in primer:
        if item in degen.keys():
            total_deg = total_deg*len(degen[item])  # Total number of primers
    return total_deg

def revcomp(line):
    old_chars = "ACGTRYMKSWVBDHIN"
    replace_chars = "TGCAYRKMSWBVHDNN"
    new = line.translate(line.maketrans(old_chars, replace_chars))[::-1]
    return new

problem=args.c+" "
problem_list=problem+"List"
PCR_problem={}
number_of_primers=0
bad_primers=set()
if args.p == "True":
    primers_seq={}
with open(args.i, "r") as fin:
    for line in fin:
        line=line.rstrip()
        if line.startswith("Primer_"):
            number_of_primers += 1
            prim=line.split()[0]
            prim=prim.split("_at_")[0]
            if not prim in PCR_problem:
                PCR_problem[prim]=0
            if args.p == "True":
                if prim in primers_seq:
                    primers_seq[prim]["GC"].append(float(line.split()[3]))
                    primers_seq[prim]["TM"].append(float(line.split()[4]))
                else:
                    primers_seq[prim]={"Length": int(line.split()[2]), "GC": [float(line.split()[3])], "TM": [float(line.split()[4])]}

        if line.startswith(problem) and not line.startswith(problem_list):
            #Dimer 1: Primer_78_at_327.1 x Primer_78_at_327.9
            dim=line.split(":")[1]
            dim=dim.split(" x ")
            if dim[0].split("_at_")[0] == dim[1].split("_at_")[0]:
                badprimer=dim[0].split("_at_")[0]
                badprimer=badprimer.replace(" ","")
                bad_primers.add(badprimer)

            for p in dim:
                primer=p.split("_at_")[0]
                primer=primer.replace(" ","")
                PCR_problem[primer]+=1

selected_primers=set()

#print(len(PCR_problem.keys()))
#print(int(float(sum(PCR_problem.values())/len(PCR_problem.keys()))*0.95))
#print(min(PCR_problem.values()))
#print(max(PCR_problem.values()))
#print(statistics.median(PCR_problem.values()))
print("Number of degenerated primers:", len(PCR_problem.keys()))
print("     Number of undegenerated primers:", number_of_primers)
print("     Total {}: {}".format(args.c,int(sum(PCR_problem.values()))))
print("     Average {} per undegenerated primer: {}".format(args.c, round(float(sum(PCR_problem.values())/number_of_primers),2)))
print("     Average {} per degenerated primers: {}".format(args.c, round(float(sum(PCR_problem.values())/len(PCR_problem.keys())),2)))
print("     Median:", int(statistics.median(PCR_problem.values())))
dict={key: value for key, value in sorted(PCR_problem.items(), key=lambda x: x[1], reverse=False)}
top=list(dict.items())[int(len(PCR_problem.keys())*0.25)] #0.25

if args.t == "cal" :
#    cutoff_dimers=int(float(sum(PCR_problem.values())/len(PCR_problem.keys()))*0.95)
#    cutoff_dimers=int(statistics.median(PCR_problem.values())*0.15)
#    cutoff_dimers=int(statistics.median(PCR_problem.values())*0.5)
#    cutoff_dimers=int(statistics.median(PCR_problem.values()))
    cutoff_dimers=int(top[1])


else:
    cutoff_dimers=int(args.t)

print("     {} cutoff: {}".format(args.c, cutoff_dimers))

for primer, prob in PCR_problem.items():
    if (prob <= cutoff_dimers) and (primer not in bad_primers):
        selected_primers.add(primer)

print("Number of selected degenerated primers:", len(selected_primers))

if args.p == "True" :
    ft=open(args.o.replace(".fasta", "_table"), "w")
    print("Primer\tSequence\tRevComp\tDegenerecy\tLength\tGC range\tTM range", file=ft)

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

                if args.p == "True" :
                    print("{}\t{}\t{}\t{}\t{}\t({}-{})\t({}-{})".format(header[1:], line, revcomp(line), degenerecy(line), primers_seq[primer_name]["Length"], min(primers_seq[primer_name]["GC"]),  max(primers_seq[primer_name]["GC"]), min(primers_seq[primer_name]["TM"]), max(primers_seq[primer_name]["TM"]) ), file=ft )

if args.p == "True" :
    ft.close()
