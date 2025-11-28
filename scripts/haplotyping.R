#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################

source(here::here("scripts/load-libs.R"))
report("R packages loaded",type="s") 

#####################
##### LOAD DATA #####
#####################

# load nd1 data
nd1.table <- readr::read_csv(here::here("data/nd1-metadata.csv"),show_col_types=FALSE) |> dplyr::relocate(accession,.before="sourceDatabase")

# load sra
nd1.sra <- ape::read.FASTA(here::here("data/sra-nd1-references.fasta")) |> as.matrix()

# report
report("Data loaded",type="s")


##########################
##### DOWNLOAD FASTA #####
##########################

# filter genbank accs
nd1.table.gb <- nd1.table |> dplyr::filter(sourceDatabase=="Nucleotide" | grepl("This study",sourceDatabase))

# download from genbank
nd1.fas <- ape::read.GenBank(access.nb=dplyr::pull(nd1.table.gb,accession))

# report
report("Sequences obtained from GenBank",type="s")


###############################
##### ALIGN AND MAKE TREE #####
###############################

# align with mafft
nd1.fas.ali <- ips::mafft(nd1.fas,exec="mafft",method="auto",thread=4)

# subset coding sequence
nd1.fas.ali.sub <- nd1.fas.ali[,2884:3858]

# join with SRA data
nd1 <- rbind(nd1.fas.ali.sub,nd1.sra)

# make ml tree
nd1.tr <- phangorn::pml_bb(nd1,model="TrN+G(4)",rearrangement="stochastic",method="unrooted")
#plot(nd1.tr$tree)

# root 
nd1.tr.root <- castor::root_in_edge(nd1.tr$tree,root_edge=which.max(nd1.tr$tree$edge.length))
# nd1.tr.root |> ape::ladderize() |> plot()

# report
report("Phylogenetic tree created",type="s")


#####################
##### PLOT TREE #####
#####################

# make tip labels
nd1.table.tips <- nd1.table |> dplyr::mutate(
    tiplabel=dplyr::case_when(
        is.na(tissueCode) ~ glue::glue("{accession} | {speciesInferred}"),
        !is.na(tissueCode) ~ glue::glue("{accession} | {tissueCode} | {speciesInferred}")
    )
)

# generate random colour pal
cpal <- withr::with_seed(seed=42,code=randomcoloR::distinctColorPalette(k=3))

# plot tree with ggtree
p <- treeio::as.treedata(nd1.tr.root) |> ggtree::ggtree(ladderize=TRUE,right=TRUE,size=0.5) %<+% nd1.table.tips
pp <- p + ggtree::geom_tiplab(offset=0.001,aes(label=tiplabel),align=FALSE,size=2) +
        #ggtree::geom_tippoint(aes(color=scientificName),size=1.5) +
        #ggtree::geom_tippoint(aes(color=speciesInferred),size=1.5) +
        ggtree::geom_tippoint(aes(color=role),size=1.5) +
        ggtree::scale_color_manual(values=cpal) +
        ggtree::theme(legend.position="none") +
        ggtree::xlim(0,0.37)

# export tree as pdf
ggtree::ggsave(filename=here::here("temp/nd1.haplotyping.tree.pdf"),plot=pp,limitsize=FALSE,width=210,height=390,units="mm")


###############
##### END #####
###############

report("All steps complete; tree plot written to {.file temp/nd1.haplotyping.tree.pdf}",type="s")
