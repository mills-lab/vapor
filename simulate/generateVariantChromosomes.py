#!/usr/bin/python

import argparse
import sys
import os
import random
import math
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

#*** Function: possibleChromePos(chrom, variantSize)
#	chrom - chomosome character (1-22,X,Y,M)
#	variantSize - total size of structural variant
#	
#	This function determines the chromosome positions that can not be selected
#	as at the start site of the strucal variant and returning a list of the positions
#***
def possibleChromPos(size):
	sites = posStartSites[:]
	sites[len(sites) - size : len(sites) + 1] = '-' * size
	
	for i in range(len(svLocations)):
		sites[max(0,svLocations[i] - size) : svLocations[i]] = '-' * min(size, svLocations[i])
	sites = [x for x in sites if x != '-']
	return sites

#*** Function: generate_break_points(numBlocks, size):
#	numBlocks - number of blocks present in the structural variant
#	size - a list containing the upper and lower bounds for the  total SV size
#
#	This function randomly selects breakpoint locations for the SV and returns it 
#	in the form of a string (chr:bp1:bp2:..:bpN)
#***
def generate_break_points(numBlocks, svSize):
	
	## list to hold the size of each block
	sizeList = random.sample(list(range(5,svSize-4)), numBlocks - 1)
	sizeList.sort()
	sizeList.append(svSize)

	## remaining usable sites for selected chromosome
	startSiteOptions = possibleChromPos(svSize)

	## randomly select start site
	try:
		variantStartSite = random.choice(startSiteOptions)
	except IndexError:
		print("Error: Your have run out of room to simulate all the variants for this chromosome.") 
		sys.exit()
	print((str(variantStartSite) + " has been selected as the start site."))

	## concatenate breakpoints into string
	bp = [variantStartSite]
	for i in range(len(sizeList)):
		bp.append(variantStartSite + sizeList[i])

	return bp
## Create a set of sequences to be used as micro insertions at breakpoints
def micro_ins_seq(seq):
	size = random.randrange(1,11)
	if len(seq) - 10000 > 10001:
		seq_start = random.randrange(10000,len(seq)-10000)
	else:
		seq_start = random.randrange(1, len(seq) - size + 1)

	if 'N' not in seq[seq_start : seq_start + size]:
		micros.append(seq[seq_start:seq_start + size])

## Add micro insertions/deletions at one or both breakpoints of SV
def micro_insertions(chrSeq, bpStart, bpStop):
	
	types = ['ins','del']
	loc = [0,0]
	## Determine if there will be an insertion/deletion at the first, second, or both breakpoints
	andOr = random.randrange(3)
	## Determine if breakpoint will have insertion or deletion
	if andOr == 0 or andOr == 2:
		loc[0] = random.choice(types)
	if andOr == 1 or andOr == 2:
		loc[1] = random.choice(types)

	microInsertSeq = ''
	if loc[0] == 'ins':
		microInsertSeq = random.choice(micros)
		duplications.append([bpStart - 1, microInsertSeq])
	if loc[1] == 'ins':
		microInsertSeq = random.choice(micros)
		duplications.append([bpStop, microInsertSeq])

	if loc[0] == 'del':
		delete_sequence(chrSeq, bpStart - random.randrange(1,11), bpStart - 1)
	if loc[1] == 'del':
		delete_sequence(chrSeq, bpStop + 1, bpStop + random.randrange(1,11))


## inverted sequence includes the breakpoints
def invert_sequence(chrSeq, start, stop):
	
	invert = chrSeq[int(start)-1:int(stop)]
	MutableSeq.reverse_complement(invert)
	begin = MutableSeq.__add__(chrSeq[:int(start)-1], invert)
	chrSeq = MutableSeq.__add__(begin, chrSeq[int(stop):])
	return chrSeq

## duplication will occur after the insert location
def duplicate_sequence(chrSeq, dupStart, dupStop, insertLoc, numDup = 1, invert = False):

	duplication = str(chrSeq[int(dupStart)-1:int(dupStop)]) * int(numDup)
	if invert == True:
		MutableSeq.reverse(duplication)
	begin = MutableSeq.__add__(chrSeq[:int(insertLoc)], duplication)
	chrSeq = MutableSeq.__add__(begin, chrSeq[int(insertLoc):])
	return chrSeq
	
## deletion includes the breakpoints
def delete_sequence(chrSeq, start, stop):

	chrSeq[int(start)-1:int(stop)] = "-" * (int(stop) - int(start) + 1)
	return chrSeq	

## insertion will occur after the insert location
def insert_sequence(chrSeq, insertLoc, insertSeq):
	
	begin = MutableSeq.__add__(chrSeq[:int(insertLoc)], insertSeq)
	chrSeq = MutableSeq.__add__(begin, chrSeq[int(insertLoc):])
	return chrSeq

parser = argparse.ArgumentParser(description='Simulate structural variants in the human reference genome')
parser.add_argument('-r', '--reference', help='reference genome in fasta format')
parser.add_argument('-s', '--structuralVariants', help='file containing structual variants to simulate')
parser.add_argument('-o', '--outdir', help='directory to save breakpoints files')
parser.add_argument('-f', '--fastaOutfile', help='file to write modififed reference')
parser.add_argument('-b', '--blacklist', help='bedfile containing blacklisted regions for chromosome')
parser.add_argument('-l', '--lineFile', help='fasta file containing examples of LINEs')
parser.add_argument('-a', '--aluFile', help='fasta file containing examples of ALUs')
parser.add_argument('-v', '--svaFile', help='fasta file containing examples of SVAs')
args = parser.parse_args()

## Set size of buffer between SVs
bufferSize = 3000

## Read in reference file
reference = list(SeqIO.parse(args.reference, 'fasta'))[0]

## Determine the length of the chromosome
chrLen = len(reference.seq)

## Create list of base pair positions for start sites
posStartSites = list(range(1, chrLen + 1))

## Read in insertion examples
ALUs = list(SeqIO.parse(args.aluFile, 'fasta'))
LINEs = list(SeqIO.parse(args.lineFile, 'fasta'))
SVAs = list(SeqIO.parse(args.svaFile, 'fasta'))

## Read in the black list
blacklist = open(args.blacklist, 'r')
for line in blacklist:
	start = int(line.strip().split("\t")[1])
	stop = int(line.strip().split("\t")[2])
	posStartSites[start - 1:stop] = "-" * (stop - start + 1)

## store duplications
duplications = []

## keep track of simulated variant locations
svLocations = [] 

## Breakpoints output file
bpOut = open(args.outdir + "/" + reference.id + "_svBreakpoints.txt", 'w')

## open structural variants file 
svFile = open(args.structuralVariants,'r')
svs = svFile.readlines()
svFile.close()

## generate micro insertion sequences
micros = []
numSVs = len(svs)
for i in range(max(2, int(math.ceil(numSVs * 0.20)))):
	micro_ins_seq(reference.seq)

## generate simulated genome from specified variants
for sv in svs:
	splitLine = sv.strip().split('\t')
	variantType = splitLine[2]
	refAllele = splitLine[0]
	altAllele = splitLine[1]
	variantSize = splitLine[3]
	
	## determine number of blocks needed for structural variant
	numBlocks = len(refAllele.split("/")[0])
	
	## generate the breakpoint
	breakpoint = generate_break_points(numBlocks, int(variantSize))

	## add breakpoint location (with buffer) to svLocations
	svLocations.append(breakpoint[0] - bufferSize)

	## remove variant location from possible start sites
	posStartSites[(breakpoint[0]-1 - bufferSize) : (breakpoint[-1] + bufferSize)] = '-' * (6000 + int(variantSize)) 

	## add variants into chromosome
	chrSeq = reference.seq.tomutable()

	insertSeq = ''
	
	## call appropriate modification function for variant type
	if variantType == "del":
		chrSeq = delete_sequence(chrSeq, breakpoint[0], breakpoint[1])
	elif variantType == "inv":
		chrSeq = invert_sequence(chrSeq, breakpoint[0], breakpoint[1])
	elif variantType == "ins":
#		insertSeq = ''
		if altAllele == 'ALU':
			insertSeq = random.choice(ALUs).seq
		if altAllele == 'LINE1':
			insertSeq = random.choice(LINEs).seq
		if altAllele == 'SVA':
			insertSeq = random.choice(SVAs).seq
		duplications.append([breakpoint[0], insertSeq])
	elif variantType == "tan_dup":
		numDup = random.randrange(50)
		duplications.append([breakpoint[1], breakpoint[0], breakpoint[1], numDup, 0])
	elif variantType == "dis_dup":
		if altAllele == 'ab/aba':
			duplications.append([breakpoint[2], breakpoint[0], breakpoint[1], 0]) 	
		elif altAllele == 'ab/bab':
			duplications.append([int(breakpoint[0])-1, int(breakpoint[1])+1, breakpoint[2], 0])
	elif variantType == "del_inv":
		if altAllele == 'b^/ab':
			chrSeq = delete_sequence(chrSeq, breakpoint[0], breakpoint[1])
			chrSeq = invert_sequence(chrSeq, int(breakpoint[1])+1, breakpoint[2])
		elif altAllele == 'a^/ab':
			chrSeq = delete_sequence(chrSeq, int(breakpoint[1])+1, breakpoint[2])
			chrSeq = invert_sequence(chrSeq, breakpoint[0], breakpoint[1])
		elif altAllele == 'b^/abc':
			chrSeq = delete_sequence(chrSeq, breakpoint[0], breakpoint[1])
			chrSeq = delete_sequence(chrSeq, int(breakpoint[2])+1, breakpoint[3])
			chrSeq = invert_sequence(chrSeq, int(breakpoint[1])+1, breakpoint[2])
	elif variantType == "dup_inv_ins":
		if altAllele == "ab/aba^":
			duplications.append([breakpoint[2], breakpoint[0], breakpoint[1], 1]) 	
		elif altAllele == "ab/b^ab":
			duplications.append([int(breakpoint[0])-1, int(breakpoint[1])+1, breakpoint[2], 1])
	elif variantType == "del_dup":
		if altAllele == "cbc/abc":
			chrSeq = delete_sequence(chrSeq, breakpoint[0], breakpoint[1])
			duplications.append([breakpoint[1], int(breakpoint[2])+1, breakpoint[3], 0])	
		elif altAllele == "aba/abc":
			chrSeq = delete_sequence(chrSeq, int(breakpoint[2])+1, breakpoint[3])
			duplications.append([breakpoint[2], breakpoint[0], breakpoint[1], 0])	
	elif variantType == "del_dup_inv":
		if altAllele == "c^bc/abc":
			chrSeq = delete_sequence(chrSeq, breakpoint[0], breakpoint[1])
			duplications.append([breakpoint[1], int(breakpoint[2])+1, breakpoint[3], 1]) 	
		elif altAllele == "aba^/abc":
			chrSeq = delete_sequence(chrSeq, int(breakpoint[2])+1, breakpoint[3])
			duplications.append([breakpoint[2], breakpoint[0], breakpoint[1], 1])
	else:
		print((variantType + ": variant type not accecpted. This variant has been skipped"))

	if random.random() <= 0.12:
		micro_insertions(chrSeq, breakpoint[0],breakpoint[-1])


	reference.seq = chrSeq.toseq()

	## add breakpoint to breakpoint file
	if variantType == 'ins':
		bpOut.write('\t'.join([':'.join([reference.id,':'.join([str(i) for i in breakpoint])]), variantType, refAllele, altAllele, insertSeq.tostring()]) + "\n")
	else:
		bpOut.write('\t'.join([':'.join([reference.id,':'.join([str(i) for i in breakpoint])]), variantType, refAllele, altAllele]) + "\n")

bpOut.close()

## any insertion (insertions and duplications) are added from the end of the chromosome
## to the begining so as to not change the indicies

## sort insertions by insertion point
duplications.sort(reverse = True)
for dup in duplications:
	## single duplication
	if len(dup) == 4:
		invert = False
		if dup[3] == 1:
			invert = True
		reference.seq = duplicate_sequence(reference.seq.tomutable(), dup[1], dup[2], dup[0], invert).toseq()
	## multiple duplications (tandum duplications)
	if len(dup) == 5:
		invert = False
		if dup[4] == 1:
			invert = True
		reference.seq = duplicate_sequence(reference.seq.tomutable(), dup[1], dup[2], dup[0], dup[3], invert).toseq()
	## insertions
	if len(dup) == 2:
		reference.seq = insert_sequence(reference.seq.tomutable(), dup[0], dup[1]).toseq()

reference.seq = reference.seq.ungap("-")

with open(args.fastaOutfile, "w") as output_handle:
	SeqIO.write(reference, output_handle, "fasta")

