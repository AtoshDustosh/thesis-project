# ***Thesis Project Bioinfo [developing]***

**note** *: This is a repository for code storage and communication when working on the thesis project.*

sub-repo: *https://github.com/AtoshDustosh/genome-file-handler*

***  

## *Project of stage 1*

- [x] Objective: able to read <\*.sam/bam>|<\*.vcf/bcf> files and get detailed data of every record. (study how to use htslib)

- Directory path: `~/htslib_Study`

- Console compile: `~/htslib_Study $ gcc *.c *.h -o test -lhts`

Please note that `'<samtools/htslib>'` is required.

- Usage: `$ ./test [<*.sam/bam>|<*.vcf/bcf>]`

Options for command lines are not realized yet.


## *Project of stage 2*

- [x] Load and index all reference genome files (\*.fa/\*.fna)

It takes about 50 seconds to load and index all bases in GRCh38_full_analysis_set_plus_decoy_hla.fa (the complete genome of human). Not bad! :)

- [x] Load and index all files of reads (\*.sam)

2.5 seconds for about 80,000 records. The time will explode if there are too many records. 

- [x] Load and index all files of variations (\*.vcf)

The time will explode if there are too many records. Almost unable to process ori.vcf

- [ ] Optimize loading processs and algorithms (put aside for the moment)

## *Project of stage 3* 

- [ ] Able to integrate any selected SNP into the reference genome

- [ ] Able to integrate any selected SV into the reference genome
