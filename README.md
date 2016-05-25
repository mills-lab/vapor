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
git clone https://github.com/xuefzhao/VaLoR.git
cd VaLoR
python setup install --user
```
