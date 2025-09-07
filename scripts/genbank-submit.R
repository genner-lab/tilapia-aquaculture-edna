#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################
#rm(list=ls())
source(here::here("scripts/load-libs.R"))
report("R packages loaded",type="s") 


#####################
##### LOAD DATA #####
#####################

# load data
tissues.meta <- readr::read_csv(here("data/tissues-master.csv"),show_col_types=FALSE)
seqs <- ape::read.FASTA(here("data/tissues.fasta"))
tissues.suppl <- readxl::read_xlsx(here("Collins_eDNAAssay_9June2025_SupportingInformation.xlsx"),sheet=1,skip=2) |> dplyr::rename_with(~ stringr::str_replace_all(.x," ","_"))


######################
##### CHECK DATA #####
######################

# subset this study
tissues.suppl <- tissues.suppl |> dplyr::filter(grepl("This study",Source_database)) 

# check against tissues metadata
tissues.suppl.list <- tissues.suppl |> 
    dplyr::distinct(Tissue_code) |> 
    dplyr::pull(Tissue_code)
missing <- setdiff(tissues.suppl.list,dplyr::pull(tissues.meta,otherCatalogNumbers))
if(rlang::is_empty(missing)) {report("All codes found in Metadata",type="s")} else {report("Codes missing from metadata",type="f")}

# subset cols in suppl
tissues.suppl.names <- tissues.suppl |> dplyr::select(Tissue_code,Inferred_mtDNA_lineage)

# make sci names in metadata
tissues.meta <- tissues.meta |> dplyr::mutate(scientificName=
    dplyr::case_when(
        taxonRank=="species" ~ glue::glue("{genus} {specificEpithet}"),
        taxonRank=="genus" ~ glue::glue("{genus} {identificationQualifier}")
    ))

# subset labels 
tissues.meta.names <- tissues.meta |> 
    dplyr::rename(Tissue_code=otherCatalogNumbers) |>
    dplyr::mutate(tiplabel=glue::glue("{Tissue_code} | {scientificName}")) |>
    dplyr::select(Tissue_code,scientificName,tiplabel)

# make genbank source mod table
tissues.meta.sm <- tissues.meta |> 
    dplyr::select(otherCatalogNumbers,scientificName,recordedBy,eventDate,country,stateProvince,waterBody,verbatimLocalityClean,decimalLatitude,decimalLongitude) |>
    dplyr::mutate(recordedBy=stringr::str_replace_all(recordedBy,";",",")) |>
    dplyr::mutate(decimalLatitude=convert_latitude(decimalLatitude),decimalLongitude=convert_longitude(decimalLongitude),Lat_Lon=glue::glue("{decimalLatitude} {decimalLongitude}")) |> 
    dplyr::mutate(Country=
        dplyr::case_when(
            !is.na(waterBody) ~ glue::glue("{country}: {stateProvince}, {verbatimLocalityClean} ({waterBody})"),
            is.na(waterBody) ~ glue::glue("{country}: {stateProvince}, {verbatimLocalityClean}")
        )) |>
    dplyr::mutate(Collection_date=lubridate::ymd(eventDate), Collection_date=format(Collection_date,"%d-%b-%Y")) |>
    dplyr::rename(Sequence_ID=otherCatalogNumbers,Collected_by=recordedBy) |>
    dplyr::select(Sequence_ID,Collected_by,Collection_date,Country,Lat_Lon,scientificName)
tissues.meta.sm |> dplyr::glimpse()

# join and check
tissues.suppl.names |> 
    dplyr::left_join(tissues.meta.names,by=dplyr::join_by(Tissue_code)) |> 
    dplyr::select(-tiplabel) |>
    dplyr::filter(Inferred_mtDNA_lineage!=scientificName) |> 
    dplyr::arrange(Tissue_code)


#####################
######## ND1 ########
#####################

# subset nd1 sequences from fasta
seqs.nd1 <- seqs[grepl("leu\\+ile",names(seqs))] 
names(seqs.nd1) <- names(seqs.nd1) |> stringr::str_replace_all("_leu\\+ile","")

# get nd1 codes in suppl
tissues.suppl.nd1 <- tissues.suppl |> 
    dplyr::filter(grepl("leu\\+ile",Source_database)) |> 
    dplyr::distinct(Tissue_code) |> 
    dplyr::pull(Tissue_code)

# check both are same
if(setequal(names(seqs.nd1),tissues.suppl.nd1)) {report("All codes match",type="s")} else {report("Codes DO NOT match",type="f")}

# make phylogenetic tree and root
tr.nd1 <- ape::nj(ape::dist.dna(as.matrix(seqs.nd1),model="TN93"))
tr.nd1$edge.length[tr.nd1$edge.length < 0] <- 0
tr.nd1.root <- tr.nd1 |> castor::root_in_edge(which.max(tr.nd1$edge.length))

# generate random colour pal
nsp.nd1 <- tissues.meta.names |> 
    dplyr::filter(Tissue_code %in% names(seqs.nd1)) |> 
    dplyr::distinct(scientificName) |> 
    nrow()
cpal <- withr::with_seed(seed=6,code=randomcoloR::distinctColorPalette(k=nsp.nd1))

# plot tree with ggtree
p <- ggtree::ggtree(tr.nd1.root, ladderize=TRUE,right=TRUE,size=0.5) %<+% tissues.meta.names
pp <- p + ggtree::geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE,size=2) +
        ggtree::geom_tippoint(aes(color=scientificName),size=1.5) +
        ggtree::scale_color_manual(values=cpal) +
        ggtree::theme(legend.position="none") +
        ggtree::xlim(0,0.2)
ggtree::ggsave(filename=here("temp/nd1.genbank.pdf"),plot=pp,limitsize=FALSE,width=210,height=297,units="mm")

# write out
tissues.meta.sm |> 
    dplyr::filter(Sequence_ID %in% names(seqs.nd1)) |>
    dplyr::mutate(Sequence_ID=glue::glue("{Sequence_ID}_nd1")) |> 
    dplyr::select(-scientificName) |> 
    readr::write_tsv(here("temp/nd1.genbank.tsv"))

# format fasta and write out
seqs.nd1.fa <- fas2tab(seqs.nd1) |> 
    dplyr::rename(Sequence_ID=label) |> 
    dplyr::left_join(tissues.meta.sm,by=dplyr::join_by(Sequence_ID)) |> 
    dplyr::mutate(Sequence_ID_iso=Sequence_ID) |>
    dplyr::mutate(Sequence_ID=glue::glue("{Sequence_ID}_nd1")) |>
    dplyr::mutate(label=glue::glue("{Sequence_ID} [organism={scientificName}] [isolate={Sequence_ID_iso}] NADH dehydrogenase subunit 1 (ND1) gene, complete cds")) |> 
    tab2fas(seqcol="nucleotides",namecol="label")
seqs.nd1.fa |> ape::write.FASTA(here::here("temp/nd1.genbank.fasta"))


###############
##### ND5 #####
###############

# subset nd5 sequences from fasta
seqs.nd5 <- seqs[grepl("nd5\\+glu",names(seqs))] 
names(seqs.nd5) <- names(seqs.nd5) |> stringr::str_replace_all("_nd5\\+glu","")

# get nd5 codes in suppl
tissues.suppl.nd5 <- tissues.suppl |> 
    dplyr::filter(grepl("nd5\\+glu",Source_database)) |> 
    dplyr::distinct(Tissue_code) |> 
    dplyr::pull(Tissue_code)

# check both are same
if(setequal(names(seqs.nd5),tissues.suppl.nd5)) {report("All codes match",type="s")} else {report("Codes DO NOT match",type="f")}

# make phylogenetic tree and root
seqs.nd5.al <- ips::mafft(seqs.nd5,exec="mafft")
tr.nd5 <- ape::nj(ape::dist.dna(seqs.nd5.al,model="TN93"))
tr.nd5$edge.length[tr.nd5$edge.length < 0] <- 0
tr.nd5.root <- tr.nd5 |> castor::root_in_edge(which.max(tr.nd5$edge.length))

# generate random colour pal
nsp.nd5 <- tissues.meta.names |> 
    dplyr::filter(Tissue_code %in% names(seqs.nd5)) |> 
    dplyr::distinct(scientificName) |> 
    nrow()
cpal <- withr::with_seed(seed=8,code=randomcoloR::distinctColorPalette(k=nsp.nd5))

# plot tree with ggtree
p <- ggtree::ggtree(tr.nd5.root, ladderize=TRUE,right=TRUE,size=0.5) %<+% tissues.meta.names
pp <- p + ggtree::geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE,size=2) +
        ggtree::geom_tippoint(aes(color=scientificName),size=1.5) +
        ggtree::scale_color_manual(values=cpal) +
        ggtree::theme(legend.position="none") +
        ggtree::xlim(0,0.2)
ggtree::ggsave(filename=here::here("temp/nd5.genbank.pdf"),plot=pp,limitsize=FALSE,width=105,height=105,units="mm")

# write out
tissues.meta.sm |> 
    dplyr::filter(Sequence_ID %in% names(seqs.nd5)) |> 
    dplyr::mutate(Sequence_ID=glue::glue("{Sequence_ID}_nd5")) |> 
    dplyr::select(-scientificName) |> 
    readr::write_tsv(here("temp/nd5.genbank.tsv"))

# format fasta and write out
seqs.nd5.fa <- fas2tab(seqs.nd5) |> 
    dplyr::rename(Sequence_ID=label) |> 
    dplyr::left_join(tissues.meta.sm,by=dplyr::join_by(Sequence_ID)) |> 
    dplyr::mutate(Sequence_ID_iso=Sequence_ID) |>
    dplyr::mutate(Sequence_ID=glue::glue("{Sequence_ID}_nd5")) |>
    dplyr::mutate(length=stringr::str_length(nucleotides)) |> # 1057
    dplyr::mutate(label=dplyr::case_when(
        length < 1057 ~ glue::glue("{Sequence_ID} [organism={scientificName}] [isolate={Sequence_ID_iso}] NADH dehydrogenase subunit 5 (ND5), partial cds"),
        length > 1057 ~ glue::glue("{Sequence_ID} [organism={scientificName}] [isolate={Sequence_ID_iso}] mitochondrion ND5 and ND6 genes, partial cds")
        )) |> 
    tab2fas(seqcol="nucleotides",namecol="label")
seqs.nd5.fa |> ape::write.FASTA(here::here("temp/nd5.genbank.fasta"))


###############
##### END #####
###############

report("Done!",type="s")
