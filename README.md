# MutAnt
Mutation Annotation tool for mutation pathogenicity classification
To work with notebooks please use jupyter.
First install dbnsfp_annotation package, then mutant

## Installation

```bash
python setup.py install
```
OR
```bash
pip install .
```
## Usage

```python
from mutant.classifier import add_mutant_annotation

# annotate vcf, genome version hg19
add_mutant_annotation('/uftp/projects/alterationsDB/annotation_output/LUAD/TCGA-91-7771/TCGA-91-7771.vcf', file_type='vcf', genome_version = 'hg19', development_mode=True, method='predict_proba')

# annotate maf, genome version hg38
add_mutant_annotation('/uftp/mvp-data/patients/Early_Adopters/BG000001/somatic.maf', file_type='maf', genome_version = 'hg38', development_mode=True, method='predict_proba')

# convert maf DNP to SNP and annotate, genome version hg38
add_mutant_annotation('/uftp/mvp-data/patients/Early_Adopters/BG001542/somatic.maf',
file_type='maf', genome_version = 'hg38', dnp_convert_info=True, method='predict_proba')
# After annotation you will get «DNP_Converted» column in your maf with True or False values,
# where True means that DNP was converted to SNP while annotation,
# False - DNP could not be converted to SNP because of for amino acid change it need in two or more changes in codon.
```