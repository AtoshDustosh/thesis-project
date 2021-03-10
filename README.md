# ***Thesis Project Bioinfo***

**note** *: This is a repository for code storage and communication when working on the thesis project.*

sub-repo: 
	https://github.com/AtoshDustosh/genome-file-handler

***  

## *Project of stage 1*

- [x] Objective: able to read <\*.sam/bam>|<\*.vcf/bcf> files and get detailed data of every record. 

- Directory path: `~/htslib_Study`

- Console compile: `~/htslib_Study $ gcc *.c *.h -o test -lhts`

Please note that `'<samtools/htslib>'` is required.

- Usage: `$ ./test [<*.sam/bam>|<*.vcf/bcf>]`

Options for command lines are not realized yet.


## *Project of stage 2*

- [x] Load and index all reference genome files (\*.fa/\*.fna)

- [x] Load and index all reads files (\*.sam)

- [x] Load and index all variation files (\*.vcf)

- [ ] Optimize loading processs (put aside for the moment)

## *Project of stage 3* 

- [ ] 

- [ ] Integrate SNP into the reference genome

- [ ] Integrate SV into the reference genome
