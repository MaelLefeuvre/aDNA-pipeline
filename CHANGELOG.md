# Changelog
## 0.4.3:
### Breaking change(s)
- Minimum supported snakemake version is now snakemake>=8

### Features
- (Breaking) added `path` entry to samples.yml. This allows users to specify custom paths for their R1 and R2 files, but changes the structure of the yaml, and forces the use of Snakemake
## 0.4.2:
### Features
- Bump to snakemake-8+ (8.30.0). (Pipeline is still retrocompatible with snakemake>=7.20.0)
### Fixes
- Fix `pmd-mask-0.3.2` and `grups-rs-0.3.2` conda environments
- Explicitly specify python3 instead of `python` in `get_bamlist_panel_coverage`

## 0.4.1:

### Features:
- Simplify reference genome selection and download on the user end
- Update to [pmd-mask v0.3.2](https://github.com/MaelLefeuvre/pmd-mask/releases/tag/v0.3.2)
- Update to [grups-rs v0.3.2](https://github.com/MaelLefeuvre/grups-rs/releases/tag/v0.3.2)

### Bugfixes
- Freeze `READ` to version [`v1.0`](https://bitbucket.org/tguenther/read/src/v1.0/), as breaking changes were introduced within subsequent versions (i.e.: Removal of `Read_output_ordered`)
- Fix rule `filter_low_coverage_samples` reaching an error state when no low-coverage samples are to be found.
- Force rule `download_HapMapII_recombination_map` to ignore and delete the downloaded `README.md` file.
- Fix bug when parsing ReichLab's AADR Compendium version tag: `resolve_aadr_version()` now accepts patch version identifiers
- Added temporary workaround to `ftputil.error.FTPOSError: [Errno 104]` arising when solving the DAG takes too long.
- Fixed redundand logging when solving DAG of jobs.
- Fixed `SSL_VERIFY` error in call to `HTTPRemoveProvider()`, when downloading some datasets.
- Fixed strict channel priority of multiple conda environments: [samtools-1.15.yml](/workflow/envs/samtools-1.15.yml), [mapdamage-2.2.1](/workflow/envs/mapdamage-2.2.1.yml), [kin-3.1.3.yml](/workflow/envs/kin-3.1.3.yml), [grups-rs-0.3.2.yml](/workflow/envs/grups-rs-0.3.2.yml), [coreutils-9.1.yml](/workflow/envs/coreutils-9.1.yml), [pmd-mask-0.3.2.yml](/workflow/envs/pmd-mask-0.3.2.yml)
- Remove dependency on gitmodules [`grups-rs`](https://github.com/MaelLefeuvre/grups-rs) and [`pmd-mask`](https://github.com/MaelLefeuvre/pmd-mask). These dependencies are now directly fetched from github when building their corresponding conda environments.

---

## 0.4.0:

### Features:
- Added support for sex-assignation
- Added `01-preprocess/00-quality-control` workflow, with `mosdepth` and `bedtools` sequencing depth / coverage calculation.

### Bugfixes:
- READ now automatically pre-filters samples having an SNP-panel coverage below the requested threshold (see. `config["kinship"]["READ"]["coverage-threshold"]` parameter). This prevents the method from bugging out when the overlap between a pair of individuals reaches 0.
- Fixes buggy ANGSD rule definition.

---

## 0.3.0:

### Features:  
- Added support for [pmd-mask](https://github.com/MaelLefeuvre/pmd-mask)

---

## 0.2.0:

### Features:
- Added support for KIN kinship estimation method
