# Changelog
## 0.4.0:
### Features:
- Added support for sex-assignation
- Added `01-preprocess/00-quality-control` workflow, with `mosdepth` and `bedtools` sequencing depth / coverage calculation.
### Bugfixes:
- READ now automatically pre-filters samples having an SNP-panel coverage below the requested threshold (see. config['kinship']['READ']['coverage-threshold'] parameter). This prevents the method from bugging out when the overlap between a pair of individuals reaches 0.

### Bugfixes:
- Fixes buggy ANGSD rule definition.
## 0.3.0:
### Features:  
- Added support for pmd-mask 
## 0.2.0:
### Features:
- Added support for KIN kinship estimation method
