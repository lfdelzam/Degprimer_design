#!/usr/bin/env python3

import os
import argparse
import re
import string
import io
import gzip


usage = 'from_refseq_gff_to_fasta.py -i -d -e -o'
description = 'This program retrieves hsp60 sequences, using gff and genome files from refseq NCBI, printing out in 5 -> 3'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='inf', help='selected lines from gff refseq files ', required=True)
parser.add_argument('-d', dest='d', help='Directory - input file, genomes.fna', required=True)
parser.add_argument('-e', dest='e', help='file extention input files', default=".fna.gz")
parser.add_argument('-o', dest='out', help='output file', required=True)

args = parser.parse_args()

old_chars = "ACGTacgt"
replace_chars = "TGCAtgca"

species={}
loc={}
ext=args.e
#GCF_003612695.1_ASM361269v1_genomic.fna.gz

set_of_genomes=set([str(f)[:-(len(ext))] for f in os.listdir(args.d) if f.endswith(ext)])

with open(args.inf, "r") as fin:
  for line in fin:
    line = line.rstrip()
#GCF_000006745.1_ASM674v1_genomic.gff.gz	NC_002505.1	Protein Homology	CDS	2832620	2834254	.	-	0	ID=cds-WP_000729136.1;Parent=gene-VC_RS12830;Dbxref=Genbank:WP_000729136.1;Name=WP_000729136.1;Note=60 kDa chaperone family%3B promotes refolding of misfolded polypeptides especially under stressful conditions%3B forms two stacked rings of heptamers to form a barrel-shaped 14mer%3B ends can be capped by GroES%3B misfolded proteins enter the barrel where they are refolded when GroES
    if line.startswith("GCF"):
      line=line.split("\t")
      file=line[0][:-7] #GCF_000006745.1_ASM674v1_genomic
      in_genome=line[1] #NC_002505.1
      start=line[4]
      stop=line[5]
      dir=line[7]
      loc[file]=(in_genome, dir, start, stop)
    else:
      line=line.split("\t")
      in_genome=line[0] #NC_002505.1
      start=line[3]
      stop=line[4]
      dir=line[6]
      if (int(loc[file][3])-int(loc[file][2])) < (int(stop)-int(start)):
                loc[file]=(in_genome, dir, start, stop)


for g in loc.keys():
    if g in set_of_genomes:
        with io.TextIOWrapper(gzip.open(os.path.join(args.d, g+ext), "rb"), encoding='utf-8') as fg, open(args.out, "a") as fout:
            #List of locis where gene GroEL is present:
            copy=False
            for line in fg:
                line=line.rstrip()
                if line.startswith(">"):
        #>NZ_CP046820.1 Vibrio metoecus strain 2011V-1015 chromosome 1, complete sequence
        #GACTTTATAAAAAAGTGGCATTTTCAACGTATAACTATACTTCAATTAACAATTCAATGACAGTGTAGTAACCTTACACT - 80 (start 30) (51 copied)
        #CCATCAGGATTTGGCCTTTTATTTCAAAAATTTGAAGTTTGTTACACGAGGAACATGCATTCATGGCTGCAGTAAAACTA - 160 (stop 100) (20 copied)
                    id=line.split()[0][1:]

                    if id == loc[g][0]:
                        header=line.replace(" ","_")
                        header=header.replace(",", " ")
                        header=re.sub(".chromosome.*$", "",header,count=1)
                        header=re.sub(".genome.*$", "",header,count=1)
                        header=re.sub(".plasmid.*$", "",header,count=1)
                        header=re.sub(".DNA.*$", "",header,count=1)
                        header=re.sub(".complete.*$", "",header,count=1)
                        header=re.sub("sp.", "sp",header,count=1)
                        header+=" [direction=+][gene=GroEL]"
                        copy=True
                        counter=0
                        tocopy=""
                        firsttime=True
                        #[('NZ_CP046820.1', '-', '1230695', '1232329')]
                        st=int(loc[g][2])
                        stp=int(loc[g][3])
                        direct=str(loc[g][1])
                else:
                    if copy:
                        counter +=len(line)
                        if counter >= st and counter <= stp:

                            if firsttime:
                                initial=counter-st + 1
                                tocopy +=line[-(initial):] #the N last characters
                                firsttime=False
                            else:
                                tocopy +=line
                        if  counter > stp:
                            ###
                            #extra for interleaved files
                            if firsttime:
                                initial=st - 1
                                firsttime=False
                                tocopy += line[initial:stp]
                            else:
                                last=counter-stp
                                tocopy += line[:-(last)] # copy line but the last characters

                            ###

#                            last=counter-stp
#                            tocopy += line[:-(last)] # copy line but the last characters
                            print(header + "[length={}]".format(len(tocopy)), file=fout)
                            s=re.sub("^.*_Vibrio", "Vibrio",header,count=1)
                            s=s.split("_")[1]
                            s=s.split()[0]
                            if s in species:
                                species[s]+=1
                            else:
                                species[s]=1

                            if direct == "-":
                                print(tocopy.translate(tocopy.maketrans(old_chars,replace_chars))[::-1], file=fout)
                            else:
                                print(tocopy, file=fout)
                            copy=False

print("# Total {} Hsp60-GroEL unique sequences from {} Vibrio species".format(sum(species.values()), len(species.keys())))
