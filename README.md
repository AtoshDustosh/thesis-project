# ***Thesis Project Bioinfo [developing]***

**note** *: This is a repository for code storage and communication when working on the thesis project.*

***  

## Compile and Run

- OS requirement: Linux

- compile:
<font color=#006400 size=4>.../THE_project/Genome_Realignment_Based_on_Variants $</font> gcc \*.c -o main -lhts -lpthread

Required additional libraries and tools:

- for compilation, htslib
- for data process, samtools & bcftools
- for simulation: varsim, art_illumina

---

## Data Load & In-memory Index

### Reference Genome (\*.fa/fna/fasta)

- source code files: genomeFa.h, genomeFa.c, genomeFaMacros.h, genomeFaMacros.c

- functions provided:

  - load \*.fa/fna/fasta files into memory
  - randomly access any base with specified position
  - randomly access any sequence (length >= 0) with specified start pos and end pos.

- strategy:
  - use binary coding for compressed storage in memory
  - use simple arrays of 64-bit to store coded bases

### Alignment Results (*.sam/*.bam)

- source code files: genomeSam.h, genomeSam.c

- functions provided:

  - load \*.sam/bam files into memory
  - sequential access and traverse all sam records
  - methods for accessing data fields of a sam record
  - iterator for all records

- strategy:
  - use singly linked list to store all records and traverse

### Variants (*.vcf/bcf)

- source code files: genomeVcf_bPlus.h, genomeVcf_bPlus.c

- functions provided:

  - load \*.vcf/bcf files into memory
  - given POS of a vcf record, return pointer to the data, or if it doesn't exist, return pointer to another data whose POS is next to but bigger than the  given POS.
  - traverse of all vcf records
    - methods for accessing data fields of a vcf record

- strategy:
  - (old version): use singly linked list to store records and traverse
  - (updated version): use b+ tree to store, index and access records, which greatly improved performance. :)

## Libraries Used

- ksw2 from GitHub. (selected)
- edlib from GitHub. (c++, incompatible with this project)
- Complete-Striped-Smith-Waterman from GitHub. (tried, but doesn't satisfy the requirements of the project)

---

## Integration of Variants (mapped reads)

Suppose you need to process a sam record
<font color=#1E90FF>

(QNAME | FLAG | RNAME | POS | MAPQ | CIGAR | RNEXT | PNEXT | TLEN | SEQ | QUAL)

</font>

### Partition of Reference Sequence

Suppose the sam record under processing has a CIGAR like: **30M 2I 4D 60M 2D 10M 2I**

We have a picture like this:

    reference: ... AAAA...AAAA -- CCCCC.....CCCCCC GG TTT....TTTT -- ...
    read seq:      AAAA...AAAA AA CCCCC.....CCCCCC -- TTT....TTTT AA
    read cigar:    <-- 30M --> II <---- 60M -----> DD <-- 10M --> II
                   <- part1 ->    <--- M_area ---> <-- part2 --->

We call the area that has the maximum count of cigar op 'M' as "M_area". And in the case above it's the area corresponding to "60M".

Extract the sequences towards boths sides of the area, but do not include any bases of the M_area.

Integration will be performed separately on these 2 parts.

### Description of Process for Integration

Extract RNAME, POS, SEQ, and locate the interval "old_ref" on reference sequence to which this SEQ is aligned.

Use accessor for variants to extract all variants that may influence the "old_ref".

And now we have the following picuture:

variants ("1.1------" indicates the 1st allele of the 1st variant. And "-----" indicates length of the allele's REF field if it's a DEL)

    old_ref:         ...AAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXBBBBBBBCCCCC...
    variants:
                |var1     |var2   |var3   |var4   |var5   |var6    |var7
    alleles of variants:
                |1.1--------------|3.1-   |4.1-   |5.1----|6.1-    |7.1----
                |1.2----- |2.1-   |3.2-------------------------    |7.2-
                                          |4.2-
                                          |4.3-

In the situation above, 1.1 and 2.1 cannot be integrated at the same time. Otherwise conflicts will occur. So does 3.2 and any of 4.1, 4.2, 4.3, 5.1, 6.1. Because when 3.2 is integrated, the bases at the position of 4.1 is deleted.

In order to prevent the reference sequence becoming too short after integrating some DELs, we expand the "old_ref" towards both sides.

And in the situation above, if we select 2.1 and 3.2 for integration, the reference sequence after integration will look like

    old_ref:         ...AAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXBBBBBBBCCCCC...
    variants:             |       |                           |
                |var1     |var2   |var3   |var4   |var5   |var6    |var7
    alleles of variants:  |       |                           |
                          |2.1-   |3.2------------------------- 
    ref after interation: |       |                           |
                     ...AADAAAAAAABBBBBBBCCCCCCCCCCCCCCCCCCCCCCCCC...

And now we have a reference sequence with variants integrated.

But this is only one instance corresponding to one combination of variants. And  $ C_n^m = C_n^{n-m}$

Actually if there are N variants that intersect with the interval, we need to filter and process cases of more than this number:

$ \sum_{i=1}^N{ C_N^i } $

### Strategy of Integration

There are 3 strategies - 1) SNP only; 2)SV only; 3) SNP and SV.

This is implemented by filtering variants by their lengths when selecting them. And this step is done before calculating the combinations of variants.

### Realignment

After integration, suppose we have the following picture:

    old reference: ... AAAA...AAAA -- CCCCC.....CCCCCC GG TTT....TTTT -- ...
    integrated:    ... ACCC...GGTA -- CCCCC.....CCCCCC -- TTT....TTTT AA ...
    read seq:          AAAA...AAAA AA CCCCC.....CCCCCC -- TTT....TTTT AA
    read cigar:        <-- 30M --> II <---- 60M -----> DD <-- 10M --> II
                       <- part1 ->    <--- M_area ---> <-- part2 --->

Use ksw2 to align integrated sequence with the read seq. Align part2 of integrated seq with part 2 of read seq. Do the same to part1.

When aligning part1, sequences should be reversed. This is because we need to make ksw2 align from left to right.

As ksw2 is based on global alignment, we need to cut off the 'D's at the end of cigar in the alignment result and create an alignment result similar to local alignment.

Then merge the alignment result of part1 and part2, together with the M_area.

### Output

Results of realignment (both CIGAR and updated POS) and integrated variants will be output as another sam record.

Integrated variants will be output using optional fields of according to SAM format.

      XV:Z:value
      
      Suppose 3 varaints are selected.
      Their IDs are (v1, v2, v3). 
      Their cnt_alleles are (2, 4, 1)
      Their POSs are (p1, p2, p3)
      Integrated alleles are (v1.1, v2.3, v3.0). (0-based index)

      The optional fields will be "XV:Z:v1;p1;1 v2;p2;3 v3;p3;0"

      As whitespaces and semicolons are not allowed in the ID of a vcf
      record, we use whitespace and semicolon to separate information. 

---

## Integration and Kmer Generation (unmapped reads)

### Input & Output

Reference Genome will be partitioned first according to input data, which follows a format like: [id_chrom,pos_start,pos_end] no whitespace or tab allowed.

    [1,123211242,123214242]
    [1,123224242,123227242]
    [1,123247242,123249242]
    ...
    [5,444249242,444249242]

id_chrom, pos_start and pos_end are all 1-based. And pos_start and pos_end indicate the positions on the chromosome instead of the absolute positions on the whole genome. (absolute position of a kmer = absolute offset of the chrom + kmer's position on the chrom)

The program will load these intervals on the reference genome and then extract all variants WITHIN the interval.

Kmers from the original ref sequence will be output first and kmers after integration with variants will be output later.

And kmers containing bases 'N', 'M' or 'R' will be ignored.

The format of output file is as follows. Whitespaces are used for "fscanf". This function doesn't recognize "," as a delimiter when there exists string input.

    # [id_chrom pos_start pos_end char_input char_output ]  (this is the comment line)
    [4 72 ACGTACGTAAGGGTCTAACCAA C G ]
    [4 77 CGTAAGGGTCTAACCAAGTAAC C C ]
    [4 88 AACCAAGTAACAAACAAAAAAA G T ]
    [4 96 AACAAACAAAAAAACCCCCGGA A T ]
    [4 102 CAAAAAAACCCCCGGATTTTAA C A ]

id_chrom, pos_start and pos_end are all 1-based like in the input file. But pos_start and poss_end indicate the absolute positions on the whole genome.

The program will try to fiter duplicated kmers from the same interval, but doesn't promise of uniqueness of every kmer. And the program will not consider filtering duplicated kmers among different intervals.

- About char_input and char_output:

      reference  ... AAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCGGGGGGGGGGGGG ...
      kmer_i                          CCCCCCCCCCCCCCCCC
      char_input(kmer_i) = A
      char_input(kmer_i) = G

  They are not the first and the last char of the kmer, but the first chars beside the kmer on the reference genome.

  When we extract the char_input and char_output of a kmer, we cannot simply extract the (k+2)mer and select the first and last bases. That's not valid. Because when you do this, you need to expand the interval towards both sides by 1 base. Afte that, another variant might be integrated in the process, and results will be different.

### Description of Process

We suppose kmers at the end of the chromosome or at the beginning will be ignored.

#### Extraction of Kmers from Original Sequence

Suppose length of kmers is 6.

                |<------------- interval -------------->|
    seq:  AAAAA AAAAAAAAACCCCCCCCCCCCCGGGGGGGGGGGGGTTTTTT TTTTTT
    kmers:
          AAAAA A
              ...
              A AAAAA
                AAAAAA
                   ... ... ... ... ... ... ... ... ... 
                                                   TTTTTT
                                                    TTTTT T
                                                        ...
                                                        T TTTTT

Positions of kmers keep their positions on the original sequence
char_input and char_output will be extracted as well.

// TODO

#### Extraction of Kmers after Integration

Suppose length of kmers is 6.

                    |<------------- interval -------------->|
    seq:  AAAAAAAGG AAAAAAAAACCCCCCCCCCCCCGGGGGGGGGGGGGTTTTTT TTTTTT
    variants:       |    |   |   |  |      |      |   |
                 <- ----->   |   |  <------>      ^   |
                     DEL    SNP SNP    DEL       INS SNP
    kmers after integration:
                    |    |   |   |  |      |      |   |
              AAAAA ------A   ...            GGGGG|..
               AAAA ------AA   ...         [kmers in INS]
                AAA ------AAA   ...             ..|GGGGG
                 AA ------AAAA   ...
                  A ------AAAAA   ...

##### Special cases

###### Long DEL

There may exist some long DEL in the interval, and in such cases, it is very likely we need to expand the reference sequences to the right. Not by (kmerLength - 1), but some other calculated value. And for simplification in implementation, we just make the expanded rbound_ref = rbound_interval + kmerLength - 1 + length_DEL_REF (length of the long DEL's REF field).

When expanding, we need to be careful not making the rbound_ref out of the chromosome's boundary.

Suppose length of kmers is 5.

                                                "expanded"
                     |<--- |<----interval---->| -- ->|------------>|
    ref seq: AAAACCCCAAGCC CCCCC CGGGGGGGGGGGTT TT TTTTTTTAAAAAAACCCCCCCCGGGGGGG
    variant:                     |<-----DEL----->|     |
    kmers             AGCC C     |               |     |
                       ... ...   |               |     |
                           CCCCC |               |     |
                             ... |               |     |
                              CC |               | TTT |
                               C |               | TTTT 

###### Intervals at both ends

For example, when length of kmers is 22.

    chromosome:
      NNNNNNN ... NNNNNNNN ... ... NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ... NNNNNN
    intervals:
         |<--interval1-->|                              |<--interval2-->|
      ----1st kmer----
       ----2nd kmer----

Directly expand (kmerLength - 1) towards both sides may result in an out of boundary error. And in such cases, we need to check the expanded boundaries of the reference sequence.

##### Calculation of positions for kmers

- Positions of kmers near a DEL

  Keep the position of the kmer's input char as its position

- Positions of kmers near a INS

  From the pos where INS is, every kmer that intersect with the inserted sequence will have a ascending offset on the position.

- Positions of kmers near a SNP

  Same as the position on the original sequence

---

## Usage

    Usage: ./main [commands] [arguments]
    Run one task at a time. Some commands conflict with each other.

    Commands:
    -- Set files. Do this first!
      outputFile [filepath] set output file
      faFile [filepath] set reference genome file
      fastqFile [filepath] set fastq file
      samFile [filepath] set sam file
      vcfFile [filepath] set vcf file
      sv_min_len [length] set minimal length for a SV. Designed for integration.
      sv_max_len [length] set maximal length for a SV. Designed for integration. 
                            Do not set this parameter too big. That may cause the
                            program running for decades! (combinations of too many 
                            variants generated) Recommended value: 300
      match [score] set score for match
      mismatch [score] set score for mismatch
      gapOpen [score] set score for gapOpen
      gapExtension [score] set score for gapExtension

    -- Program infos
      verbose verbose mode. Private usage only. Do no choose this.

    -- Simple operations
      countRec count records for all input files. Execute successfully only when 
                    the files' formats are correct
      firstLines [number] print the first [number] lines for all files to console
      extractChrom [chrom_idx 1-based] extract bases of the selected chromosome 
                                            and write into designated output file 
                                            together with the chromosome's info field
      statistics_vcf  collect statistics from a vcf file. Statistics includes number 
                      of snp, small_ins, small_del, mnp, sv_ins, sv_del and other types 
                      of variants.  Variants using tags like <INV> will be  classified 
                      separately

    -- GRBV operations
      threads [NUM_threads] use multi-threads methods to run the program. This only 
                            works for integrateVcfToSam.
      selectBadReads [MAPQ_threshold] select reads with MAPQ lower than 
                                            MAPQ_threshold from previously set sam 
                                            file and then output them into the 
                                            previously specified output file
      integrateVcfToSam [integration_strategy] integrate variants from *.vcf file 
                                                    with *.sam file. This will perform 
                                                    realignment for all reads in the 
                                                    *.sam file with new created 
                                                    reference genome. It's actually 
                                                    the main purpose of the project. 
                        [integration_strategy]: [1] SNP and small indel only; 
                                                [2] SV only; 
                                                [3] SNP, small indel and SV
