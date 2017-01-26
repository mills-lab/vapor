#!/usr/bin/python
import argparse
import sys
import os
import random
import math
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Simulate structural variants in the human reference genome')
parser.add_argument('-r', '--reference', help='reference genome in fasta format')
parser.add_argument('-s', '--structuralVariants', help='file containing structual variants to simulate')
parser.add_argument('-o', '--outdir', help='directory to save breakpoints files')
args = parser.parse_args()

## Read in reference file
reference = SeqIO.to_dict(SeqIO.parse(args.reference, 'fasta'))

## generate chromosome length dictionary
chrWeight = []
for seq in reference:
	os.system("touch " + args.outdir + "/" + reference[seq].id + "_sv.txt")
	extString = [reference[seq].id] * max(1, (len(reference[seq].seq) // 1000000))
	chrWeight.extend(extString)

## file to keep track of number of variants simulated 
numSimFile = open(args.outdir + '/svCounts.txt', 'w')
## open structural variants file 
svFile = open(args.structuralVariants,'r')
## remove header
svFile.readline()
## generate simulated genome from specified variants
for line in svFile:
	print((line.strip()))
	splitLine = line.strip().split('\t')
	variantType = splitLine[2]
	refAllele = splitLine[0]
	altAllele = splitLine[1]
	variantSize = splitLine[3].split('-')
	simNumberEst = int(splitLine[4])
	
	simNumber = random.randrange(math.floor(simNumberEst * 0.9), math.ceil(simNumberEst * 1.1))
	print(simNumber)
	numSimFile.write('\t'.join([variantType, refAllele, altAllele, str(simNumber)]) + '\n')

	## generate simulated breakpoints
	for i in range(simNumber):
		print(("Chromosome " + str(chr) + " has been selected"))

		## select a chromosome from the chromosome list (weighted random choice)
		chr = random.choice(chrWeight)

		## select a size for the SV within the given range
		svSize = random.randint(int(variantSize[0]),int(variantSize[1]))

		## write SV to appropriate chromosome file
		outfile = open(args.outdir + "/" + chr + "_sv.txt", 'a')
		outfile.write('\t'.join([refAllele, altAllele, variantType, str(svSize)]) + '\n')
		outfile.close()
