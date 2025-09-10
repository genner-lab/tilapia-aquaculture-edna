[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16699183.svg)](https://doi.org/10.5281/zenodo.16699183)

# Tilapia Aquaculture eDNA

Collins RA, Saxon AD, Shechonge AH, Kishe MA, Ngatunga BP & Genner MJ. Environmental DNA-based quantification of an invasive tilapia species in Tanzanian inland aquaculture.

Code and data for manuscript on monitoring Tanzania tilapia aquaculture using eDNA.

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
```

### Repository contents

* **`data/`** - Raw and processed data used in analyses.
    - `events-master.csv` - table of eDNA sampling events
    - `extractions-master.csv` - table of eDNA extractions and metadata
    - `qpcr-results.csv` - raw data from qPCR analyses
    - `tissues.fasta` - FASTA file containing nucleotide data from sequenced tissue samples
    - `tissues-master.csv` - table containing metadata from tissue samples
* **`renv/`** - Settings for the R environment.
* **`scripts/`** - R scripts to run analyses.
    - `genbank-submit.R` - script to generate and check files to submit to GenBank
    - `load-libs.R` - script to load package libraries and custom functions
    - `models.R` - script to run site occupancy models
    - `qpcr.R` - script to process raw qPCR data
* **`temp/`** - Temporary file directory ignored by git.
* `LICENSE` - Legal stuff
* `README.md` - This file
* `renv.lock` - R packages required and managed by renv
* `.gitignore` - files and directories ignored by git
* `.Rprofile` - activates renv
* `.R-version` - required version of R (used by renv-installer)
* `Collins_eDNAAssay_9June2025_SupportingInformation.xlsx` - Supporting Information for publication
