# VaLoR
##Description
*VaLoR* is a structural variants (SVs) validator based on long reads.  *VaLoR* could read in VCF files in format 4.1 or newer, and the corresponding long reads in bam format to sore each prediction. 

##Required third-party resources
```
python:   https://www.python.org/ 
samtools: http://samtools.sourceforge.net/
```
##Usage:
Usage: valor [Options] [Parameters]

Options: 

	vcf - if input structural variants are encoded in vcf format
	bed - if input structural variants are encoded in bed  format
	
Parameters:

	--sv-input:		input file in bed or vcf format
	--output-path:		folder where the recurrence plot will be kept
	--reference:		reference genome that pacbio files are aligned against
	--pacbio-input:		absolute path of input pacbio file

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

Run VaLoR on a bed file
```
valor bed --sv-input ../input.bed --output-path ../valor_result/ --reference ../reference.fa --pacbio-input sample.bam
```

Run VaLoR on a vcf file
```
valor vcf --sv-input ../input.vcf --output-path ../valor_result/ --reference ../reference.fa --pacbio-input sample.bam
```


##Instructions for input file:
###SV input in bed format:
SV input in bed format should be consist of 4 columns, specifying the chromosome , start, end and description of the variance. The description column should specify the type of the variance (eg. DEL, DUP, INV, INS). For insertions, the description column should have either the sequence of the length of insertion, or both, separated by '_' 

Here's an example of the bed input:
```
chr1	831220	833732	DUP
chr1	2093108	2095944	DEL
chr1	2994972	2995286	DEL
chr1	3935209	3935209	INS_191
```


###SV input in valor format:
If complex structural variants are encoded within the vcf input, the INFO column should describe the complex event in the similar format as the examples listed below:
```
chr6	119011725	chr6:119011725:119012294:119012497:119013928	t	<DISDUP>	92	PASSSVTYPE=disdup;END=119012294;insert_point=chr6:119013928;Other=abc/abc_abc/abca_chr6:119011725:119012294:119012497:119013928	GT	0/1
chr9	25594104	chr9:25594104:25594327:25594448	c	<DEL_INV>	78	PASS	SVTYPE=del_inv;END=25594448;del=chr9:25594104-25594327;inv=chr9:25594327-25594448;Other=ab/ab_ab/b^_chr9:25594104:25594327:25594448	GT	0/1
chr3	12913583	chr3:11956093:12913583:12913799	G	<DEL_DUP_INV>	-24	LowQual	SVTYPE=del_dup_inv;END=12913799;del=chr3:11956093-12913583;dup_inv=chr3:12913583-12913799;insert_point=chr3:12913799;Other=ab/ab_a/bb^_chr3:11956093:12913583:12913799	GT	0/1
chr7	40412883	chr7:40412883:40412973:40413111	A	<DUP_INV>	98	PASS	SVTYPE=dup_inv;END=40412973;insert_point=chr7:40412883;Other=ab/ab_a^ab/abb_chr7:40412883:40412973:40413111	GT	1/0
```



