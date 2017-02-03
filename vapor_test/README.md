This repository provides PacBio sequences from chr10 of NA12878 as a test case. 

Sequences in bam format can be downloaded [here](https://umich.box.com/v/vapor)

To run vapor on the test:
```
vapor vcf --sv-input vapor_test.vcf --output-path ./vapor_vcf_test --pacbio-input ./vapor_test_bam/NA12878.hg19.Pacbio.chr10_127191658_134925265.bam --reference hg19.fa

vapor bed --sv-input vapor_test.bed --output-path ./vapor_bed_test --pacbio-input ./vapor_test_bam/NA12878.hg19.Pacbio.chr10_127191658_134925265.bam --reference hg19.fa
```

Summary of vapor results can be found at *vapor_test.vcf.vapor* for vcf input and *vapor_test.bed.vapor* for bed input. Recurrence plots for each SV are kept under the folder *vapor_vcf_test* and *vapor_bed_test*