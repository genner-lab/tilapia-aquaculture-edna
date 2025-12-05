#!/usr/bin/env Rscript

###########################
#### LIBRARIES ####
###########################

suppressPackageStartupMessages({
    library("here")
    library("glue")
    library("tidyverse")
    library("optparse")
    library("DECIPHER")
    library("Biostrings")
    library("parallel")
})

###################
#### ARGUMENTS ####
###################

# get args
option_list <- list( 
    make_option(c("-a","--assay"), type="character"),
    make_option(c("-t","--threads"), type="numeric")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#opt$assay <- "uro"
#opt$threads <- 8

###################
#### FUNCTIONS ####
###################

# function to run AmplifyDNA to get primer efficiencies
extract_eff <- function(qseq,ppair,temp,label){
    res <- DECIPHER::AmplifyDNA(primers=Biostrings::DNAStringSet(ppair),myDNAStringSet=Biostrings::DNAStringSet(qseq),maxProductSize=16000,annealingTemp=temp,P=8e-7,minEfficiency=0,taqEfficiency=TRUE,maxDistance=0.4,maxGaps=2,processors=1)
    out <- vector("character", 2)
    if(length(res)==0) {out[] <- NA_character_}
    else {
        out[1] <- as.character(as.numeric(stringr::str_replace_all(res@ranges@NAMES[1],"%.*",""))/100)
        out[2] <- as.character(res@ranges@width[1])
    }
    return(tibble::tibble(
        accession=label,
        primerEfficiency=out[1],
        ampliconLength=out[2]
    ))
}

# FUN TO CALC PROBE EFFICIENCY
# ions ≈ 0.2 M (i.e. 0.2 mol/L sodium-equivalent) — this is often used in DECIPHER’s own primer/amplicon-simulation examples. 
# P ≈ 4e-7 (i.e. 4×10⁻⁷ M, which is 0.4 µM) — used in DECIPHER vignettes/examples. 
probe_eff <- function(template,probe) {
    template.bs <- Biostrings::DNAString(template)
    probe <- Biostrings::DNAString(probe)
    loc <- Biostrings::matchPattern(probe,template.bs,max.mismatch=6,with.indels=FALSE)
    probe.loc <- Biostrings::DNAStringSet(Biostrings::Views(template.bs,start(loc),end(loc)))
    if(length(probe.loc)==0) {out <- 0} 
    else {
        out <- DECIPHER::CalculateEfficiencyPCR(Biostrings::DNAStringSet(probe), Biostrings::reverseComplement(probe.loc),temp=60,ions=0.2,P=4e-7)
    }
    names(out) <- names(template)
    return(out)
}


###################
#### LOAD DATA ####
###################

# read in hits
cc <- c("accession","scientificName","taxid","kingdom","phylum","class","order","family","genus","species","nucleotides")
hits.tbl <- readr::read_tsv(here::here(glue::glue("crabs/{opt$assay}.hits.derep.tsv")),col_types=cols(.default="c"),col_names=cc)
#glimpse(hits.tbl)

# read in all
all.tbl <- readr::read_tsv(here::here(glue::glue("crabs/{opt$assay}.tsv")),col_types=cols(.default="c"),col_names=cc)
#glimpse(all.tbl)

# subset
hits.full <- all.tbl |> dplyr::filter(accession %in% dplyr::pull(hits.tbl,accession))

# set primers
if(opt$assay=="leuco") {
    primers <- c("CTCTGCCCTACTGCACTCG","CATGGGGCTTATACGGATGAGA")
    probe <- "AGCACCATAGTCGTAGCCGGCATCT"
} else if(opt$assay=="nilo") {
    primers <- c("TGGAGGTTTTACCCTACAGACC","GTCGAAGGGAGCTCGGTTA")
    probe <- "AGTGTCTGACTAATCCTTCCCGCCTGAC"
} else if(opt$assay=="uro") {
    primers <- c("CTAAGCCTCGTGTTAACTCCAG","TTTGATGTATTCGGCAGGTGGA") #c("TCCACCTGCCGAATACATCAAA","CTGGAGTTAACACGAGGCTTAG")
    probe <- "GGGCGTGAGTGTTGTTGTTGACGTTG"#"CAACGTCAACAACAACACTCACGCCC"#
} else {
    stop("Assay must be 'leuco', 'nilo', 'uro'.")
}

# make list of named nucleotides
hits.nuc <- dplyr::pull(hits.full,nucleotides)
names(hits.nuc) <- dplyr::pull(hits.full,accession)


###########################
#### PRIMER EFFICIENCY ####
###########################

writeLines("\n\nCalculating primer and probe efficiencies ...")

# get eff for all primers 
primers.res <- mcmapply(function(x,y) extract_eff(x,ppair=primers,temp=60,label=y), x=hits.nuc, y=names(hits.nuc), mc.cores=opt$threads, SIMPLIFY=FALSE, USE.NAMES=TRUE)

# tabulate
primers.res.tbl <- dplyr::bind_rows(primers.res)
#primers.res.tbl |>  filter(primerEfficiency!=0) |> arrange(primerEfficiency) |> print(n=Inf)

##########################
#### PROBE EFFICIENCY ####
##########################

# calc probe eff on all targets
probes.res <- mcmapply(function(x) probe_eff(template=x,probe=probe), x=hits.nuc, mc.cores=opt$threads, SIMPLIFY=TRUE, USE.NAMES=TRUE)

# tabulate
probes.res.tbl <- probes.res |> 
    tibble::enframe(name="accession",value="probeEfficiency") |> 
    dplyr::mutate(probeEfficiency=as.character(round(probeEfficiency,3)))
#probes.res.tbl |> filter(probeEfficiency!=0) |> arrange(probeEfficiency) |> print(n=Inf)

#######################
#### COMBINE WRITE ####
#######################

# join and tidy
hits.eff <- hits.full |> 
    dplyr::left_join(primers.res.tbl,by=join_by(accession)) |> 
    dplyr::left_join(probes.res.tbl,by=join_by(accession)) |>
    dplyr::mutate(assay=opt$assay) |>
    dplyr::select(assay,accession,phylum,class,order,family,species,primerEfficiency,probeEfficiency,ampliconLength) |>
    dplyr::filter(primerEfficiency!=0 | probeEfficiency!=0) |> 
    dplyr::arrange(desc(primerEfficiency),desc(probeEfficiency),species)

# write out
hits.eff |> readr::write_csv(here::here(glue::glue("crabs/{opt$assay}.efficiencies.csv")))

writeLines("\nDone\n")
