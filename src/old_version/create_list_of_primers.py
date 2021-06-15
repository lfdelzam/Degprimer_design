#!/usr/bin/env python3

import os
import argparse
import string
import statistics

usage = 'create_list_of_primers.py -i -o -m -M -s'
description = 'This program prints out the list of select primers for in silico primer evaluation'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i', dest='i', help='input file *_table', required=True)
parser.add_argument(
    '-o', dest='o', help='output file list_of_primers.txt, selected primers', required=True)
parser.add_argument(
    '-m', dest='m', help='minimum amplicon size', default=200)
parser.add_argument(
    '-M', dest='M', help='maximum amplicon size', default=1000)
parser.add_argument(
    '-s', dest='s', help='primer gene target', default='gene')
args = parser.parse_args()

primers={}
with open(args.i, "r") as fin:
    for line in fin:
        line=line.rstrip()
        if line.startswith("Primer_"):
            #Primer_168_at_528       GTNGAAGAAGGYCARGS       SCYTGRCCTTCTTCNAC       17      (47.06-64.71)   (48.18-58.1)
            primer=line.split()[0]
            pos=primer.split("_at_")[1]
            primers[primer]={ "pos": pos, "fw": line.split()[1], "rv": line.split()[2] }

counter=1
primer_names=list(primers.keys())
with open(args.o, "w") as fout:
    for i in range(0,len(primers.keys())):
        for l in range(i+1,len(primers.keys())):
            size=int(primers[primer_names[l]]["pos"])-int(primers[primer_names[i]]["pos"])
            if size >= int(args.m) and size <= int(args.M):
                print("{}_{} {} {} {} {}".format(args.s, counter, primer_names[i], primers[primer_names[i]]["fw"], primer_names[l], primers[primer_names[l]]["rv"]), file=fout)
                counter += 1
