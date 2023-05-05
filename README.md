# A basic ancient DNA kinship estimation pipeline

## Installation

0. Clone this repository

  ```Bash
user@desktop:~$ git clone --recursive git@github.com:MaelLefeuvre/aDNA-pipeline.git
  ```

1. Install [Conda](https://docs.conda.io/en/latest/)

  - Check the documentation of [miniconda3](https://docs.conda.io/en/latest/miniconda.html) and review the detailled [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) beforehand.

  - On a x86_64 bits Linux architecture:

    ```Bash
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    wget $MINICONDA_URL && bash Miniconda3-latest-Linux-x86_64.sh
    ```

2. [Optional] Install [Mamba](https://github.com/mamba-org/mamba)

  - Integrates seamlessly with conda, and greatly speeds-up the environment solver.

  - On a x86_64 bit Linux architecture:

    ```Bash
    conda install -n base -c conda-forge mamba
    ```

3. Install [Snakemake](https://snakemake.github.io/)

  - Check the [documentation](https://snakemake.readthedocs.io/en/stable/) and review the detailled [installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

  - A dedicated environment is available within this repository. On a x86_64 bit Linux architecture:

    ```Bash
    mamba env install -f ./workflow/envs/snakemake-7.12.0.yml
    ```

## Execution

1. Activate snakemake

  ```Bash
  mamba activate snakemake-7.12.0.yml
  ```

2. Place your raw FASTQ files in the `original-data/samples` directory.
  - This pipeline currently only allows for PE data.
  - Currently, filenames and directory structure must **strictly** follow the current nomenclature: 
    ```original-data/samples/{sample-name}/{run-id}/{sample-name}_R{1,2}.fastq.gz```
    i.e.:

  ```Bash
  original-data
  └── samples
      ├── MT23
      │   ├── 2020_run0
      │   │   ├── MT23_R1.fastq.gz
      │   │   └── MT23_R2.fastq.gz
      │   ├── 2020_run1
      │   │   ├── MT23_R1.fastq.gz
      │   │   └── MT23_R2.fastq.gz
      │   └── 2021_run1
      │       ├── MT23_R1.fastq.gz
      │       └── MT23_R2.fastq.gz
      ├── MT26
      │   ├── 2020_run0
      │   │   ├── MT26_R1.fastq.gz
      │   │   └── MT26_R2.fastq.gz
      │   ├── 2020_run1
      │   │   ├── MT26_R1.fastq.gz
      │   │   └── MT26_R2.fastq.gz
      │   └── 2021_run1
      │       ├── MT26_R1.fastq.gz
      │       └── MT26_R2.fastq.gz
      └── SK27300
          ├── 2020_run0
          │   ├── SK27300_R1.fastq.gz
          │   └── SK27300_R2.fastq.gz
          ├── 2020_run1
          │   ├── SK27300_R1.fastq.gz
          │   └── SK27300_R2.fastq.gz
          ├── 2021_run1
          │   ├── SK27300_R1.fastq.gz
          │   └── SK27300_R2.fastq.gz
          └── 2021_run2
              ├── SK27300_R1.fastq.gz
              └── SK27300_R2.fastq.gz
  ```

3. Specify the targeted sample names into the `config/config.yml`, and modify parameters at leasure.

4. Run the pipeline (using all cores + 64GB of RAM).

  ```Bash
  snakemake --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete --cores `proc` --resources mem_mb=64000
  ```

## Workflow overview (Simplified Rulegraph)
![workflow](dags/simplified-rulegraph.svg)
