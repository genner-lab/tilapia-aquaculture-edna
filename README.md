[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16699183.svg)](https://doi.org/10.5281/zenodo.16699183)

# Tilapia Aquaculture eDNA

Collins RA, Saxon AD, Shechonge AH, Kishe MA, Ngatunga BP & Genner MJ. (2026). Environmental DNA-based quantification of an invasive tilapia species in Tanzanian inland aquaculture. _Aquaculture, Fish and Fisheries_ [https://doi.org/10.1002/aff2.70237](https://doi.org/10.1002/aff2.70237).

Code and data for article on monitoring Tanzania tilapia aquaculture using eDNA.

### Download data and install R packages

```bash
# clone the repository onto your local machine
git clone https://github.com/genner-lab/tilapia-aquaculture-edna.git
cd tilapia-aquaculture-edna
mkdir temp
# install R packages - requires R v4.4.1
Rscript -e "renv::restore()"
```

### Process qPCR 

```bash
# run scripts to generate tables and figures
scripts/qpcr.R 
scripts/models.R
scripts/haplotyping.R
scripts/primer-efficiency.sh
```

### Repository contents

* **`data/`** - Raw and processed data used in analyses.
    - `Ciezarek_MtDNA_GenomicAncestry_Correspondence.txt` - data from Ciezarek et al. studies
    - `events-master.csv` - table of eDNA sampling events
    - `extractions-master.csv` - table of eDNA extractions and metadata
    - `nd1-metadata.csv` - metadata describing all ND1 sequences used
    - `qpcr-results.csv` - raw data from qPCR analyses
    - `sra-nd1-references.fasta` - ND1 reference sequences obtained from SRA
    - `tilapia-ml-tree.nwk` - phylogenetic tree used for ND1 haplotyping
    - `tissues.fasta` - FASTA file containing nucleotide data from sequenced tissue samples
    - `tissues-master.csv` - table containing metadata from tissue samples
* **`renv/`** - Settings for the R environment.
* **`scripts/`** - R scripts to run analyses.
    - `genbank-submit.R` - script to generate and check files to submit to GenBank
    - `haplotyping.R`  -script to identify haplotyped individuals
    - `load-libs.R` - script to load package libraries and custom functions
    - `models.R` - script to run site occupancy models
    - `primer-efficiency.R` - script to estimate potential off-target primer binding
    - `primer-efficiency.sh` - script to estimate potential off-target primer binding
    - `qpcr.R` - script to process raw qPCR data
* **`temp/`** - Temporary file directory ignored by git.
* `LICENSE` - Legal stuff
* `README.md` - This file
* `renv.lock` - R packages required and managed by renv
* `.gitignore` - files and directories ignored by git
* `.python-version` - python version required
* `.Rprofile` - activates renv
* `.R-version` - required version of R (used by renv-installer)
* `Collins_eDNAAssay_SuppFigures_13Jan2026.pdf` - Supporting figures for publication
* `Collins_eDNAAssay_SuppTables_13Jan2026.xlsx` - Supporting tables for publication
