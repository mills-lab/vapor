# VaPoR

## Description

*VaPoR* is a structural variants (SVs) validator based on long reads from sequencing technology represented by PacBio.  *VaPoR* reads in BED files or VCF files in format 4.1 or newer, and the corresponding long reads in bam format to score each prediction. 

## Required third-party resources
```
python:   https://www.python.org/ 
samtools: http://samtools.sourceforge.net/
```
## Usage:
Usage: vapor [Options] [Parameters]

Options: 

	vcf - if input structural variants are encoded in vcf format
	bed - if input structural variants are encoded in bed format
	
Parameters:

	--sv-input:		input file in bed or vcf format
	--output-path:		folder where the recurrence plot will be kept
	--reference:		reference genome that pacbio files are aligned against
	--pacbio-input:		absolute path of input pacbio file



## Quick Start
Download and Install

```
git clone https://github.com/mills-lab/vapor.git
cd vapor
python setup.py install --user
export PATH=$PATH:$HOME/.local/bin
```
please add 'export PATH=$PATH:$HOME/.local/bin' to your bashrc file


Run VaPoR on a bed file
```
vapor bed --sv-input ../input.bed --output-path ../vapor_result/ --reference ../reference.fa --pacbio-input ../sample.bam
```

Run VaPoR on a vcf file
```
vapor vcf --sv-input ../input.vcf --output-path ../vapor_result/ --reference ../reference.fa --pacbio-input ../sample.bam
```


## Instructions for input file:
### SV input in bed format:
SV input in bed format should be consist of 4 columns, specifying the chromosome , start, end and description of the variance. The description column should specify the type of the variance (eg. DEL, DUP, INV, INS). For insertions, the description column should have either the sequence of the length of insertion, or both, separated by '_' 

Here's an example of the bed input:
```
chr1	831220	833732	DUP
chr1	2093108	2095944	DEL
chr1	2994972	2995286	DEL
chr1	3935209	3935209	INS_191
chr12   28226407        28226407        INS_GTAAGTTCCAGTGATCTATTGTATAGCAGGATGACTGTAGTTAATGTACTGTATTCTTTTTTTTTTTTTTCTTTTTTTTTATTTTTTTATTATACTCTAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAACGTGTCATCTAGCATTAGGTATATCTCCCAATGCTATCCCTCCCCCTCCCCCGACCCCACCACAGTCCCAGAGTGTGATATTCCCTTCCTGTGTCCAAGTGTTCTCATTGTTCAATTCCCACCTATGAGTGAGAATATGCGGTGTTTGTTTTTTGTTCTTGCGATAGTTTACTGAGAATGATGGTTTCCAATTTCATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATGTGCCACATTTTCTTAATTCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTTTCTTTATAGCAGCATGATTTATAGTCCTTTGGGTATATACCCAGTAATGGGATAGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACAGTCCCACCAACAGTGTAAAAGTGTTCCTATTTCTCCACATCCTCTCCAGCACCTGTTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGATATCTCATAGTGGTTTTGATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATTTCTTCATGTGTTTTTTGGCTGCATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATGTCCTTCGCCCACTTTTTGATGGGGTTTTTGTTTTTTTCTTGTAAATTTGTTTGAGTTCATTGTAGATTCTGGATATTAGCCCTTTGTCAGATGAGTAGGTTGCGAAAATTTTTTCCCATGTTGTAGGTTGCCTGTTCACTCTGATGGTAGTTTCTTTTGCTGTGCAGAAGCTCTTTAGTTTAATTAGATCCCATTTGTCAATTTTGGCTTTTGTTGCAATTGCTTTTGGTGTTTTGGACATGAAGTCCTTGCCCACGCCTATGTCCTGAATGGTAATGCCTAGGTTTTCTTCTAGGGTTTTTTATGGTTTTAGGTCTAACGTTTAAATCTTTAATCCATCTTGAATTGATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCCAGCACATTTATTAAATAGGGAATCCTTTCCCCATTGCTTGTTTTTCTCAGGTTTGTCAAAGATCAGATAGTTGTAGATATGTGGCATTATTTCTGAGGGCTCTGTTCTGTTCCATTGATCTATATCTCTGTTTTGGTACCAGTACCATGCTGTTTTGGTTACTGTAGCCTTGTAGTATAGTTTGAAGTCAGGTAGTGTGATGCCTCCAGCTTTGTTCTTTTGGCTTAGGATTGACTTGGCGATGCGGGCTCTTTTTTGGTTCCATATGAACTTTAAAGTAGTTTTTTCCAATTCTGTGAAGAAAGTCATTGGTAGCTTGATGGGGATGGCATTGAATCTGTAAATTACCTTGGCAGTATGGCCATTTTCACGATATTGATTCTTCCTACCCATGAGCATGGAATGTTCTTCCATTTGTTTGTGTCCTCTTTTATTTCCTTGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCACATCCCTTGTAAGTTGGATTCCTAGCTATTTTATTCTCTTTGAAGCAATTGTGAATGGGAGTTCACTCATGATTTGGCTCTCTGTTTGTCTGTTGTTGGTGTATAAGAATGCTTTGTGATTTTTGTACATTGATTTTGTATCCTGAGACTTTGCTGAAGTTGCTTATCAGCTTAAGGAGATTTTGGGCTGAGACGATGGGGTTTTCTAGATAAACAATCATGTCGTCTGCAAACAGGGACAATTTGACTTCCTCTTTCCTAATTGAATACCCTTTATTTCCTTCTCCTGCCTGATTGCCCTGGCCAGAACTTCCAACACTATGTTGAATAGTAGCGGTGAGAGAGGGCATCCCTGTCTTGTG
```
### SV output for input in bed format:
additional columns would be added to the bed format:

	VaPoR_qs: highest score of assessed read
	VaPoR_gs: proportion of reads supportive of predicted structure versus all reads assessed
	VaPoR_GT: genotype of the proposed SV as assessed by VaPoR
	VaPoR_GQ: genotype quality of the proposed SV as assessed by VaPoR
	VaPoR_Rec: list of scores for all reads assessed




### SV input in vcf format:
If complex structural variants are encoded within the vcf input, the INFO column should describe the complex event in the similar format as the examples listed below:
```
chr6	119011725	chr6:119011725:119012294:119012497:119013928	t	<DISDUP>	92	PASSSVTYPE=disdup;END=119012294;insert_point=chr6:119013928;Other=abc/abc_abc/abca_chr6:119011725:119012294:119012497:119013928	GT	0/1
chr9	25594104	chr9:25594104:25594327:25594448	c	<DEL_INV>	78	PASS	SVTYPE=del_inv;END=25594448;del=chr9:25594104-25594327;inv=chr9:25594327-25594448;Other=ab/ab_ab/b^_chr9:25594104:25594327:25594448	GT	0/1
chr3	12913583	chr3:11956093:12913583:12913799	G	<DEL_DUP_INV>	-24	LowQual	SVTYPE=del_dup_inv;END=12913799;del=chr3:11956093-12913583;dup_inv=chr3:12913583-12913799;insert_point=chr3:12913799;Other=ab/ab_a/bb^_chr3:11956093:12913583:12913799	GT	0/1
chr7	40412883	chr7:40412883:40412973:40413111	A	<DUP_INV>	98	PASS	SVTYPE=dup_inv;END=40412973;insert_point=chr7:40412883;Other=ab/ab_a^ab/abb_chr7:40412883:40412973:40413111	GT	1/0
```


### SV output for input in vcf format:
additional information will be added to the info column of the VCF:

	##INFO=<ID=VaPoR_GS,Number=1,Type=Float,Description="VaPoR Score, representing the percentage of transverse long reads that support the prediction">'
	##INFO=<ID=VaPoR_GT,Number=1,Type=String,Description="Genotype with the highest likelihood as estimated by VaPoR">'
	##INFO=<ID=VaPoR_GQ,Number=1,Type=Float,Description="Genotype quality score - likelihood of the second most likely genotype on a -log10 normalized scale"'
	##INFO=<ID=VaPoR_REC,Number=.,Type=Float,Description="Similarity scores assigned to each of the reads traversings the predicted SV">'
 
