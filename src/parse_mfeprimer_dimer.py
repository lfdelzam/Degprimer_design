#!/usr/bin/env python3

import os
import argparse
import string
import statistics

usage = 'parse_mfeprimer_dimer.py -i -f -o -t -m -M -s'
description = 'This program selects primers according to dimer formation'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i', dest='i', help='input file mfeprimer3', required=True)
parser.add_argument(
    '-f', dest='f', help='input file .fasta, primers', required=True)
parser.add_argument(
    '-o', dest='o', help='output file .txt, selected primers', required=True)
parser.add_argument(
    '-t', dest='t', help='output file table, selected primers', required=True)
parser.add_argument(
    '-m', dest='m', help='minimum amplicon size', default=300)
parser.add_argument(
    '-M', dest='M', help='maximum amplicon size', default=900)
parser.add_argument(
    '-s', dest='s', help='primer gene target', default='GroEl')

args = parser.parse_args()

degen = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
         "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"],
         "V": ["A", "G", "C"], "D": ["A", "G", "T"],
         "H": ["A", "T", "C"], "B": ["C", "G", "T"],
         #       "I": "N" }
         "I": ["A", "G", "C", "T"], "N": ["A", "G", "C", "T"]}


def degenerecy(primer):
    total_deg = 1
    for item in primer:
        if item in degen.keys():
            total_deg = total_deg * len(degen[item])  # Total number of primers
    return total_deg


def revcomp(line):
    old_chars = "ACGTRYMKSWVBDHIN"
    replace_chars = "TGCAYRKMSWBVHDNN"
    new = line.translate(line.maketrans(old_chars, replace_chars))[::-1]
    return new


def bad_primer_list(file, problem, problem_list):
    bad_primers = set()
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(problem) and not line.startswith(problem_list):
                # Dimer 1: Primer_78_at_327.1 x Primer_78_at_327.9
                dim = line.split(":")[1]
                dim = dim.split(" x ")
                opt1 = dim[0].split("_at_")[0]
                opt1 = opt1.replace(" ", "")
                opt2 = dim[1].split("_at_")[0]
                opt2 = opt2.replace(" ", "")
                if opt1 == opt2:
                    bad_primers.add(opt1)

    return bad_primers


def get_sequences(file, list_of_primers):
    Sequences = {}
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                primer_name = line.split("_at_")[0][1:]
            else:
                if primer_name in list_of_primers:
                    Sequences[primer_name] = line
    return Sequences


def p_to_p_dimers(file, problem, problem_list, bad_primers):
    PCR_problem = {}
    primers_info = {}

    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith("Primer_"):
                prim = line.split()[0]
                pos = prim.split("_at_")[1]
                pos = pos.split(".")[0]
                prim = prim.split("_at_")[0]

                if prim not in bad_primers:
                    if prim not in PCR_problem:  # initializing primer count
                        PCR_problem[prim] = set()
                    if prim in primers_info:
                        primers_info[prim]["GC"].append(float(line.split()[3]))
                        primers_info[prim]["TM"].append(float(line.split()[4]))
                    else:
                        primers_info[prim] = {
                            "Length": int(
                                line.split()[2]), "GC": [
                                float(
                                    line.split()[3])], "TM": [
                                float(
                                    line.split()[4])], "pos": pos}

            if line.startswith(problem) and not line.startswith(problem_list):
                # Dimer 1: Primer_78_at_327.1 x Primer_78_at_327.9
                dim = line.split(":")[1]
                dim = dim.split(" x ")
                primer1 = dim[0].split("_at_")[0]
                primer1 = primer1.replace(" ", "")
                primer2 = dim[1].split("_at_")[0]
                primer2 = primer2.replace(" ", "")
                if primer1 not in bad_primers:
                    PCR_problem[primer1].add(primer2)
                if primer2 not in bad_primers:
                    PCR_problem[primer2].add(primer1)

    return PCR_problem, primers_info


bad_primers = bad_primer_list(args.i, "Dimer ", "Dimer List")
PCR_problem, primers_info = p_to_p_dimers(
    args.i, "Dimer ", "Dimer List", bad_primers)

print("Total number of self-dimer primers (excluded):", len(bad_primers))
print("Number of selected degenerated primers:", len(PCR_problem.keys()))

list_of_primers = list(PCR_problem.keys())
Sequences = get_sequences(args.f, list_of_primers)
counter = 1
# Printing out selected primers
with open(args.t, "w") as fo, open(args.o, "w") as fout:
    print("#{}\t{}\t{}\t{}\t{}\t{}\t{})".format("primer", "sequence 5'-> 3'", "revcomp", "degenerecy", "Length", "GC range", "TM range C"), file=fo)
    print("{} {} {} {} {}".format("#Gene_target", "Forward_primer", "sequence_(5'-> 3')", "Reverse_primer", "sequence_(5'-> 3')"), file=fout)
    NO_pirmer=True
    for i in range(0, len(list_of_primers)):
        primer = list_of_primers[i]
        no_pair = PCR_problem[primer]
        pos_ini = primers_info[primer]["pos"]
        TM_min = min(primers_info[primer]["TM"])
        TM_max = max(primers_info[primer]["TM"])
        for o in range(i + 1, len(list_of_primers)):
            primer_pair = list_of_primers[o]
            if (len(primer_pair) > 0) and (primer_pair not in no_pair):
                pos_final = primers_info[primer_pair]["pos"]
                TM_min2 = min(primers_info[primer_pair]["TM"])
                TM_max2 = max(primers_info[primer_pair]["TM"])
                size = int(pos_final) - int(pos_ini)
                delta_TMmin = int(abs(TM_min2 - TM_min))
                delta_TMmax = int(abs(TM_max2 - TM_max))

                if size >= int(
                        args.m) and size <= int(
                        args.M) and delta_TMmax <= 5 and delta_TMmin <= 5:
                    NO_pirmer=False
                    print("{} {} {} {} {}".format(args.s + "_" + str(counter),
                                                  primer + "_at_" + primers_info[primer]["pos"],
                                                  Sequences[primer],
                                                  primer_pair + "_at_" + primers_info[primer_pair]["pos"],
                                                  revcomp(Sequences[primer_pair])),
                                                  file=fout)
                    counter += 1
            print("{}\t{}\t{}\t{}\t{}\t({}-{})\t({}-{})".format(primer + "_at_" + pos_ini,
                                                            Sequences[primer],
                                                            revcomp(Sequences[primer_pair]),
                                                            degenerecy(Sequences[primer]),
                                                            primers_info[primer]["Length"],
                                                            min(primers_info[primer]["GC"]),
                                                            max(primers_info[primer]["GC"]),
                                                            TM_min, TM_max),
                                                            file=fo)
    if NO_pirmer:
        print("There is no primer that fullfil the amplicon size criteria, or deltaTM is > 5", file=fout)                                                    
