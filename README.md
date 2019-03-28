# filterSomatic
Filter SNVs or Indels based on non genotype specific quality parameters.

Takes a vcf as its input and generates filtered vcfs as its output. Filtering can greatly reduce a vcf in size. This allows a vcf to be loaded into R or Python with reasonable memory usage, were complex genotype filtering can be performed.

Filtering is performed on several parameters:
- FILTER field in vcf
- Remove all SNVs or Indels
- QUAL
- Max alleles
- Chromosomes
- MQ
- Blacklists

## USAGE:
```bash
python filterSomatic.py -i filterSomatic.ini
```
An example ini file is provided.

## Dependencies
- Python >= 3.6.5
- vcftools >= 0.1.15
- snpeff >= 4.1
- coreutils >= 8.30 (Standard unix tools)
