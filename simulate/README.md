# Simulate Structural Variants in Genome  

## Description

__Summary__: Together, these scripts can be used to simulate non-overlapping structural variants in the reference genome.  
  
__Details__: Two scripts are to be run in succession. The first takes in the  reference genome as well as a file describing the varaints to be simulated (see below for more info). The purpose of this first script is to distribute the variants to  be simulated across the chromosomes allowing the actual simulation of variants on each chromosome to run independently of one another. Using the generated variant counts for a given chromosome, the second script will randomly select locations for the SVs on the chromosome, avoiding regions specified in the blacklist and ensuring there is no overlap between simulated variants. 

## Output Files
Sequences in bam format can be downloaded [here] (https://umich.box.com/v/vapor) 

## Required third party resources
```
python:      https://www.python.org/ 
biopython:   https://github.com/biopython/biopython.github.io/
```

## Details  

### Supported Structural Variants

| Variant Type                    | Abbreviation | Reference Allele | Alternate Allele     |
|---------------------------------|--------------|------------------|----------------------|
| Deletions                       | del          | a/a              | /a                   |
| Insertions                      | ins          | a/a              | ALU, LINE1, SVA      |
| Inversion                       | inv          | a/a              | a/a^                 |
| Tandem Duplications             | tan_dup      | a/a              | a/aa                 |
| Dispersed Duplications          | dis_dup      | ab/ab            | ab/aba, ab/bab       |
| Deletion-Inversion              | del_inv      | ab/ab, abc/abc   | b^/ab, a^/ab, b^/abc |
| Duplication-Inversion-Insertion | del_inv_ins  | ab/ab            | ab/aba^, ab/b^ab     |
| Deletion-Duplication            | del_dup      | abc/abc          | cbc/abc, aba/abc     |
| Deletion-Duplication-Inversion  | del_dup_inv  | abc/abc          | c^bc/abc, aba^/abc   |

### Blacklisted regions

Regions within the genome where simulated SVs are not desired can be specified in a bed file and input into the script using the -b BLACKLIST flag. 

### Insertion sequences

Three types of insertions are supported: LINES, SVA, and ALU. Examples of each of these type of insertions must be provided to the script.

## Usage:  
  
### 1) selectVariantChromosomes.py     
  
python selectVariantChromosome.py  [parameters]  
  
Parameters:  
```
  -r REFERENCE, --reference REFERENCE:  
  		reference genome in fasta format (all chromosomes 1-22,X,Y)
  -s STRUCTURALVARIANTS, --structuralVariants STRUCTURALVARIANTS:		
  		file containing structual variants to simulate  
  -o OUTDIR, --outdir OUTDIR: 
  		directory to save chromosome files SV simulation counts 
``` 

#### Input Files
##### Structural variants file format:  

This is a tab separated text file with a header row and the following columns: ref,alt,sv,size,num

```
ref:   
   reference genotype structure
alt:   
   alternate genotype structure
sv:    
   structural varaiant description (from table above)
sive:  
   range of bp size for SV separated by '-' (a random length will be selected from this range)
num:   
   number of SVs of this type/size to be simulated (a random number will be selected between [0.9 * num, 1.1 * num]
```

Example:  
```
ref   	alt  	sv    	size	num
a/a     /a      del     50-5000 4000
a/a     a/a^    inv     50-5000 1000
a/a     LINE1   ins     0-0     200
a/a     a/aa    tan_dup 50-5000 1000
ab/ab   ab/aba  dis_dup 5000-100000     10
ab/ab   ab/bab  dis_dup 50-5000 100
ab/ab   a^/ab   del_inv 5000-100000     10
abc/abc b^/abc  del_inv 50-5000 100
ab/ab   ab/aba^ dup_inv_ins     50-5000 1000
abc/abc aba^/abc        del_dup_inv     50-5000 1000
```

#### Output Files  
One file will be created for each chromosome containing the SVs that should be simulated on the chromosome. These files are the input for the next script. An additional file, svCounts.txt, contains the number of each type of variant will be simulated.  

### 2) generateVariantChromosomes.py

python generateVariantChromosomes.py [parameters]  
  
Parameters:  

```
  -r REFERENCE, --reference REFERENCE  
             reference genome in fasta format  (single chromosome)
  -s STRUCTURALVARIANTS, --structuralVariants STRUCTURALVARIANTS  
             file containing structual variants to simulate  (generated from selectVariantChromosomes.py)
  -o OUTDIR, --outdir OUTDIR  
             directory to save breakpoints file  
  -f FASTAOUTFILE, --fastaOutfile FASTAOUTFILE  
             file path to write modififed reference  
  -b BLACKLIST, --blacklist BLACKLIST  
             bedfile containing blacklisted regions for chromosome  
  -l LINEFILE, --lineFile LINEFILE  
             fasta file containing examples of LINEs  
  -a ALUFILE, --aluFile ALUFILE  
             fasta file containing examples of ALUs  
  -v SVAFILE, --svaFile SVAFILE  
             fasta file containing examples of SVAs  
```

#### Output Files

##### Breakpoint file
File containing the breakpoints and types of structural varinat contained in the chromosome.

##### Fasta file
File containing chromosome with simulated structural variants in fasta format.

#### Important Note
If no fasta file is generated and you see the following error message: 
```
"Error: Your have run out of room to simulate all the variants for this chromosome."
```

Ensure the combined length of your SV's and [buffer size (3000 bp) * number of SVs], is less than the length of the chromosome. Running the simulation again with the same number of SVs to simulate may result in an alternative arrangement of SVs allowing all variants to fit in the chromosome. 
