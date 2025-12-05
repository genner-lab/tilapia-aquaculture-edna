#!/usr/bin/env sh

# INSTALL CRABS
# see instructions in dyfi-edna project

# set up
mkdir crabs


###################
#### TAXONOMY ####
###################

# download ncbi taxonomy 
crabs --download-taxonomy --output crabs
# or symlink from previous


##########################
#### MIDORI ND5 LEUCO ####
##########################

# download midori nd5
#crabs --download-midori --output crabs/leuco.fasta --gb-number 268_2025-08-14 --gene ND5 --gb-type uniq
wget https://www.reference-midori.info/download/Databases/GenBank268_2025-08-14/RAW_sp/uniq/MIDORI2_UNIQ_SP_NUC_GB268_ND5_RAW.fasta.gz -O crabs/leuco.fasta.gz
gzip -d crabs/leuco.fasta.gz

# import midori
crabs --import --import-format midori --input crabs/leuco.fasta --names crabs/names.dmp --nodes crabs/nodes.dmp --acc2tax crabs/nucl_gb.accession2taxid --output crabs/leuco.tsv --ranks 'kingdom;phylum;class;order;family;genus;species'

# in silico pcr
crabs --in-silico-pcr --input crabs/leuco.tsv --output crabs/leuco.hits.tsv --forward CTCTGCCCTACTGCACTCG --reverse CATGGGGCTTATACGGATGAGA --threads 8 --untrimmed crabs/leuco.nohits.tsv --mismatch 5

# derep
crabs --dereplicate --input crabs/leuco.hits.tsv --output crabs/leuco.hits.derep.tsv --dereplication-method 'unique_species'

# run get primer eff
scripts/primer-efficiency.R -a leuco -t 8


#########################
#### MIDORI ND1 NILO ####
##########################

# download midori nd1
#crabs --download-midori --output crabs/nilo.fasta --gb-number 268_2025-08-14 --gene ND1 --gb-type uniq
wget https://www.reference-midori.info/download/Databases/GenBank268_2025-08-14/RAW_sp/uniq/MIDORI2_UNIQ_SP_NUC_GB268_ND1_RAW.fasta.gz -O crabs/nilo.fasta.gz
gzip -d crabs/nilo.fasta.gz

# import midori
crabs --import --import-format midori --input crabs/nilo.fasta --names crabs/names.dmp --nodes crabs/nodes.dmp --acc2tax crabs/nucl_gb.accession2taxid --output crabs/nilo.tsv --ranks 'kingdom;phylum;class;order;family;genus;species'

# in silico pcr
crabs --in-silico-pcr --input crabs/nilo.tsv --output crabs/nilo.hits.tsv --forward TGGAGGTTTTACCCTACAGACC --reverse GTCGAAGGGAGCTCGGTTA --threads 8 --untrimmed crabs/nilo.nohits.tsv --mismatch 5

# derep
crabs --dereplicate --input crabs/nilo.hits.tsv --output crabs/nilo.hits.derep.tsv --dereplication-method 'unique_species'

# run get primer eff
scripts/primer-efficiency.R -a nilo -t 8


########################
#### MIDORI ND6 URO ####
########################

# download midori nd6
#crabs --download-midori --output crabs/uro.fasta --gb-number 268_2025-08-14 --gene ND6 --gb-type uniq
wget https://www.reference-midori.info/download/Databases/GenBank268_2025-08-14/RAW_sp/uniq/MIDORI2_UNIQ_SP_NUC_GB268_ND6_RAW.fasta.gz -O crabs/uro.fasta.gz
gzip -d crabs/uro.fasta.gz

# import midori
crabs --import --import-format midori --input crabs/uro.fasta --names crabs/names.dmp --nodes crabs/nodes.dmp --acc2tax crabs/nucl_gb.accession2taxid --output crabs/uro.tsv --ranks 'kingdom;phylum;class;order;family;genus;species'

# in silico pcr
crabs --in-silico-pcr --input crabs/uro.tsv --output crabs/uro.hits.tsv --forward CTAAGCCTCGTGTTAACTCCAG --reverse TTTGATGTATTCGGCAGGTGGA --threads 8 --untrimmed crabs/uro.nohits.tsv --mismatch 5

# derep
crabs --dereplicate --input crabs/uro.hits.tsv --output crabs/uro.hits.derep.tsv --dereplication-method 'unique_species'

# run get primer eff
scripts/primer-efficiency.R -a uro -t 8


#############
#### FIN ####
#############

printf "...\nDone\n"
