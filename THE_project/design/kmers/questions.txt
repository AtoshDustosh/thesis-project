Suppose length of kmers is 6.
(note that most "|" are just marks used for better reading)

            "expanded"                                    "expanded"
              |<--- |<------------- interval -------------->| --->|
    seq:  AAAAAAAGG AAAAAAACCCCCAACCCCCCCCCCCCCGGGGGGGGGGGGGT TTTTT TTTTTT
    variants: |  |  |    |        |   |  |      |      |   |
              |  <- ----->        |   |  <------>      ^   |
              |      DEL         SNP SNP    DEL       INS SNP
              |                  [T] [T]                  [C]
    kmers after integration:
                    |    |        |   |  |      |      |   |
              AAAGG ------A       |   |    ...    GGGGG|..
               AAGG ------AC      |   |    ...   [kmers in INS]
                AGG ------ACC     |   |    ...       ..|GGGGG
                 GG ------ACCC    |   |    ...
                  G ------ACCCC   |   |    ...
                             CCCAAT   |   
                                ...   |
                                  TCCCCC
                                  |   |
                                 AACCCT
                                  | ...
                                  |   TCCCCC
                                 ATCCCT   [combinations of 2 SNP]
                                  TCCCTC   [combinations of 2 SNP]
                                  ... [other combinations]

问题：variants中的第一个变异DEL
    给定区间为interval，向两侧延伸（kmerLength-1）的部分为“expanded”
    第一个DEL的前半部分在expanded里，后半部分在interval里
    这种情况下，是只需要把DEL的后半部分在interval内的整合，expanded的部分不整合，
    （得到DEL附近的kmer为AAAGGA,AAGGAC,...,GACCCC）
    还是把DEL的所有部分都整合？
    （得到DEL附近的kmer为AAAACC,AAACCC,AACCCC)

    另外，之前提到的
      “先根据变异和区间确定新的区间，然后提取kmer”
      “保证新区间的长度=原区间长度+2*expanded长度”

    那么上面的这个情况下，整合第一个DEL之后，区间后面的部分也需要延伸对应被删除的长度吗？
    （比如DEL删了6个碱基，区间就向后延伸6个碱基）