#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC

parser=argparse.ArgumentParser(description='Filter SPAdes contigs according to length, GC content, coverage')
parser.add_argument('-i', help='multifasta of SPAdes contigs',type=str)
parser.add_argument('-l', help='minimum contig length, set to 0 to skip',type=int, default='0')
parser.add_argument('-min_cov', help='minimum coverage, set to 0 to skip',type=float, default='0')
parser.add_argument('-max_cov', help='maximum coverage, set to 0 to skip',type=float, default='0')
parser.add_argument('-min_gc', help='minimum GC, set to 0 to skip',type=float, default='0')
parser.add_argument('-max_gc', help='maximum GC, set to 0 to skip',type=float, default='0')
parser.add_argument('-o', help='Output file name',type=str)

args = parser.parse_args()

assembly = SeqIO.parse(args.i, "fasta")

sequences=[]

for contig in assembly:
	lencov=[int(contig.id.split("_")[3]),float(contig.id.split("_")[5])]
	if lencov[0] >= args.l or args.l==0:
		if lencov[1] >= args.min_cov or args.min_cov==0:
			if lencov[1] <= args.max_cov or args.max_cov==0:
				if GC(contig.seq) >= args.min_gc or args.min_gc==0:
					if GC(contig.seq) <= args.max_gc or args.max_gc==0:
						sequences.append(contig)

SeqIO.write(sequences, args.o, "fasta")
