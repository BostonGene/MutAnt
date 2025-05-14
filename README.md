# MutAnt
Mutation Annotation tool for mutation pathogenicity classification
To work with notebooks please use jupyter.
First install dbnsfp_annotation package, then mutant.
Models are attached to the repository.

## Installation

```bash
python setup.py install
```
OR
```bash
pip install -e .
```
## Usage

```python
from mutant.classifier import add_mutant_annotation

# example, annotate vcf, genome version hg19 (this vcf file are not attached)
add_mutant_annotation('TCGA-A1-A0SB.vcf', file_type='vcf', genome_version = 'hg19', development_mode=True, method='predict_proba')

# example, annotate maf, genome version hg38
add_mutant_annotation('TCGA-A1-A0SB.maf', file_type='maf', genome_version = 'hg38', development_mode=True, method='predict_proba')

# convert maf DNP to SNP and annotate, genome version hg38
add_mutant_annotation('TCGA-A1-A0SB.maf', file_type='maf', genome_version = 'hg38', dnp_convert_info=True, method='predict_proba')
# After annotation you will get «DNP_Converted» column in your maf with True or False values,
# where True means that DNP was converted to SNP while annotation,
# False - DNP could not be converted to SNP because of for amino acid change it need in two or more changes in codon.
```