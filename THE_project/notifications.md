# Notifications about the Program

## Vcf Exceptions

### Mixed records

For records as the following, bug will occur during integration between PTV042 and PTV043. 

When the SNP of PTV042 is integrated, PTV043 will be ignored, which is in fact wrong.

```
chr4	80	PTV042	AAGGGTTTAACCAAGTAACAAAAAAAAAAACCCCCGGATT	A,CAGGGTTTAACCAAGTAACAAAAAAAAAAACCCCCGGATT	200	PASS	.
chr4	82	PTV043	G	C,T	200	PASS	.

```

But if you split PTV042 into 2 records, the error will be fixed.

Note that in VCFv4.3 format, multiple records with the sam POS is permitted.

```
chr4	80	PTV04X	A	C	200	PASS	.
chr4	80	PTV042	AAGGGTTTAACCAAGTAACAAAAAAAAAAACCCCCGGATT	A	200	PASS	.
chr4	82	PTV043	G	C,T	200	PASS	.
```

