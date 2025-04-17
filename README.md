[![DOI](https://zenodo.org/badge/x.svg)](https://zenodo.org/badge/latestdoi/x)

# Tilapia Aquaculture eDNA

Collins RA, Saxon AD, Shechonge AH, Kishe MA, Ngatunga BP & Genner MJ. Environmental DNA-based quantification of an invasive tilapia species in Tanzanian inland aquaculture. _Journal_ [https://doi.org/x](https://doi.org/x).

Code and data for manuscript on monitoring Tanzania tilapia aquaculture using eDNA.

### Repository contents (A-Z)

* **`data/`** - Raw and processed data used in analyses.
    - `events-master.csv` - table of eDNA sampling events
    - `extractions-master.csv` - table of eDNA extractions and metadata
    - `qpcr-results.csv` - raw data from qPCR analyses
    - `tissues.fasta` - FASTA file containing nucleotide data from sequenced tissue samples
    - `tissues-master.csv` - table containing metadata from tissue samples
* **`renv/`** - Settings for the R environment.
* **`scripts/`** - R scripts to run analyses.
    - `models.R` - script to run site occupancy models
    - `qpcr.R` - script to process raw qPCR data
* **`temp/`** - Temporary file directory ignored by git.
* `LICENSE` - Legal stuff
* `README.md` - This file
* `renv.lock` - R packages required and managed by renv
* `.gitignore` - files and directories ignored by git
* `.Renviron` - ??
* `.Rprofile` - activates renv
* `.R-version` - required version of R (used by renv-installer)
