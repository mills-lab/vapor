# VaLoR
##Description
*VaLoR* is a structural variants (SVs) validator based on long reads.  *VaLoR* could read in VCF files in format 4.1 or newer, and the corresponding long reads in bam format to sore each prediction. 

##Required third-party resources
```
python:   https://www.python.org/ 
samtools: http://samtools.sourceforge.net/
```
##Usage:
Usage: VaLoR [Options] [Parameters]

Options: 

	svelter
	vcf
	bed
	
Parameters:

	--sv-input:
	--output-path:
	--pacbio-input:
	--reference:
	
Optional Parameters:

	--window-size
	--sv-type
	

##Quick Start
Download and Install
```
git clone https://github.com/mills-lab/valor.git
cd valor
python setup install --user
```

##Requirement for input file:
###requirement for SV input in bed format:
SV input in bed format should be consist of 4 columns, specifying the chromosome , start, end and description of the variance. The description column should specify the type of the variance (eg. DEL, DUP, INV, INS). For insertions, the description column should have either the sequence of the length of insertion, or both, separated by '_' 

Here's an example of the bed input:
```
chr1	831220	833732	DUP
chr1	2093108	2095944	DEL
chr1	2994972	2995286	DEL
chr1	3935209	3935209	INS_191
```
