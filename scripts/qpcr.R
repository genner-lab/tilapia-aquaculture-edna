#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################

source(here::here("scripts/load-libs.R"))
report("R packages loaded",type="s") 

#####################
##### LOAD DATA #####
#####################

# load qpcr data
qcr.results <- readr::read_csv(here("data/qpcr-results.csv"),show_col_types=FALSE)

# load extraction data (sample info)
extractions.master <- readr::read_csv(here("data/extractions-master.csv"),show_col_types=FALSE)

# load events data (event info)
events.master <- readr::read_csv(here("data/events-master.csv"),show_col_types=FALSE)

# load tissue data (species in ponds)
tissues.master <- readr::read_csv(here("data/tissues-master.csv"),show_col_types=FALSE)

# report
report("Data loaded",type="s")


################################
##### NO TEMPLATE CONTROLS #####
################################

# report
report("/// NO TEMPLATE CONTROLS (NTC) ///",type="i")

# summarise no template controls
qcr.results.neg <- qcr.results |> 
    dplyr::filter(role=="NTC") |> 
    dplyr::arrange(desc(copies))

# summarise negative and positive NTCs
qcr.results.neg |> 
    dplyr::mutate(positive=if_else(is.na(cq),0,1),negative=if_else(is.na(cq),1,0)) |> 
    dplyr::group_by(assay) |>
    dplyr::summarise(negative=sum(negative),positive=sum(positive)) |> 
    knitr::kable(caption="Sum of positive and negative NTCs",format="simple")

# report positive NTC
qcr.results.neg |>
    dplyr::filter(!is.na(copies)) |>
    dplyr::select(experiment,well,assay,role,cq,copies) |>
    knitr::kable(digits=2,caption="Positive NTCs",format="simple")


###########################
##### QUANT STANDARDS #####
###########################

# report
report("/// QUANTIFICATION STANDARDS ///",type="i")

# check standards
qcr.results |> 
    dplyr::filter(is.na(exclude)) |> 
    dplyr::filter(role=="Standard") |> 
    dplyr::group_by(assay,copies) |>
    dplyr::summarise(meancq=mean(cq),npositive=n(),.groups="drop") |>
    dplyr::arrange(assay,desc(copies)) |>
    knitr::kable(format.args=list(scientific=FALSE,big.mark=","),digits=2,caption="Quantification standards",format="simple")


######################
##### CLEAN DATA #####
######################

# remove positive controls and standards
# change copies to zero for excluded samples
# split the sample and dilution factor cols
qcr.results.clean <- qcr.results |> 
    dplyr::filter(role=="Unknown") |>
    dplyr::select(-role,-thresholdSetting,-thresholdCq,-baselineSetting,-baselineStart,-baselineEnd) |>
    dplyr::mutate(copies=dplyr::if_else(is.na(copies),0,copies)) |>
    tidyr::separate(col=sample,into=c("sample","dilutionFactor"),sep=" \\(",remove=FALSE) |>
    dplyr::mutate(dilutionFactor=as.numeric(str_split_fixed(dilutionFactor,"",4)[,3])) |>
    dplyr::rename(extractionTubeLabel=sample)

# check extractions sample names match
dont_print({
    # (should be empty)
    qcr.results.clean |> dplyr::filter(!extractionTubeLabel %in% dplyr::pull(extractions.master,extractionTubeLabel))
    # (should all have "no sample" in notes)
    extractions.master |> dplyr::filter(!extractionTubeLabel %in% dplyr::pull(qcr.results.clean,extractionTubeLabel))
})

# join with the sample info
qcr.results.ext <- qcr.results.clean |> dplyr::left_join(extractions.master,by="extractionTubeLabel")

# check all blanks are present (don't print)
dont_print({
    # blanks (should all be zero copies)
    qcr.results.ext |> dplyr::filter(!is.na(blankType)) |> dplyr::filter(copies>0)
    # check event ids (should all be blanks)
    qcr.results.ext |> dplyr::filter(!eventID %in% dplyr::pull(events.master,eventID))# |> print(n=Inf)
    # check event ids (should all be KIL)
    events.master |> dplyr::filter(!eventID %in% dplyr::pull(qcr.results.ext,eventID))# |> print(n=Inf)
})


##################################
##### FIELD AND LAB CONTROLS #####
##################################

# get numbers for field/lab blank controls
report("/// FIELD AND LAB BLANK CONTROLS ///",type="i")

# numbers of blank extractions
qcr.results.ext |> 
    dplyr::filter(!is.na(blankType)) |> 
    dplyr::distinct(extractionTubeLabel,blankType) |> 
    dplyr::count(blankType) |>
    knitr::kable(caption="Number of field and lab blank extractions",format="simple")

# numbers of blank PCRs
qcr.results.ext |> 
    dplyr::filter(!is.na(blankType)) |> 
    dplyr::count(blankType) |>
    knitr::kable(caption="Number of field and lab blank PCRs",format="simple")


##################
##### EVENTS #####
##################

# add events data
# clean up - change names, remove blanks and superfluous columns
qcr.results.joined <- qcr.results.ext |> 
    dplyr::left_join(events.master,by="eventID") |>
    dplyr::filter(is.na(blankType)) |>
    dplyr::select(eventID,localityID,localitySite,,nickname,assay,extractionTubeLabel,extractionDate,dilutionFactor,repVol,copies) |>
    dplyr::mutate(localityID=stringr::str_replace_all(localityID,"MG","Site"),localitySite=stringr::str_replace_all(localitySite,"MG","Site")) |>
    dplyr::mutate(localityID=glue::glue("{localityID}-{nickname}"), localitySite=stringr::str_replace_all(localitySite,"-",glue::glue("-{nickname}-")))
# qcr.results.joined |> filter(extractionTubeLabel=="10.06.19/01")

# get df of names 
dont_print({
    qcr.results.names <- qcr.results.joined |> dplyr::distinct(localityID,localitySite,extractionTubeLabel) 
    qcr.results.names |> dplyr::distinct(localityID) |> dplyr::arrange(localityID) #|> print(n=Inf) # Site05 and Site13 has two names
    qcr.results.names |> dplyr::distinct(localitySite) |> dplyr::arrange(localitySite) #|> print(n=Inf) 
    qcr.results.names |> dplyr::distinct(extractionTubeLabel) |> dplyr::arrange(extractionTubeLabel) #|> print(n=Inf) 
})

# get numbers for field/lab blank controls
report("/// EDNA SAMPLING EVENTS ///",type="i")

# count number localities etc
qcr.results.names |> dplyr::summarise(
        site=length(unique(localityID)),
        pond=length(unique(localitySite)),
        waterSample=length(unique(extractionTubeLabel))
        ) |>
    tidyr::pivot_longer(tidyr::everything(),names_to="samplingLevel",values_to="n") |>
    knitr::kable(caption="Number of eDNA samples from sites, ponds, and water sample replicates",format="simple")

# add a new column for copies/L of water given x dilutionFactor, 1 uL in each PCR reaction, 105 uL extraction volume, and x mL filtration volume (repVol)
# sum total copies in a filter extraction
# get number PCR replicates and postives
qcr.results.by.extraction <- qcr.results.joined |> 
    dplyr::mutate(copiesLitre=((((copies*dilutionFactor)*105)/repVol)*1000)) |>
    dplyr::select(localityID,localitySite,extractionTubeLabel,assay,copiesLitre) |>
    dplyr::group_by(extractionTubeLabel,assay) |> 
    dplyr::mutate(copiesLitreExtraction=sum(copiesLitre),nPcrReps=n(),nPcrRepsPositive=length(which(copiesLitre>0)),propPcrRepsPostive=nPcrRepsPositive/nPcrReps) |> 
    dplyr::ungroup() |>
    dplyr::arrange(localitySite,extractionTubeLabel,assay,desc(copiesLitre))
# qcr.results.by.extraction |> filter(localitySite=="Site02-Msauzi-A" & assay=="leucostictus") |> print(n=Inf)
# qcr.results.by.extraction 


##########################
##### SUMMARISE QPCR #####
##########################

# summarise by groups
dont_print({
    # summarise by 'localityID'
    qcr.results.by.extraction |> summarise_qpcr(grouping=localityID)
    # summarise by 'localitySite'
    qcr.results.by.extraction |> summarise_qpcr(grouping=localitySite)
    # summarise by 'extractionTubeLabel'
    qcr.results.by.extraction |> summarise_qpcr(grouping=extractionTubeLabel)
})

# summarise by 'localitySite' (pond)
qcr.results.summarised <- qcr.results.by.extraction |> summarise_qpcr(grouping=localitySite)


##########################
##### TISSUE SAMPLES #####
##########################

# Summarising tissue sampling
report("/// TISSUE SAMPLES ///",type="i")

# format the tissues
tissues.master.format <- tissues.master |> 
    dplyr::filter(eventID %in% dplyr::pull(qcr.results.ext,eventID),genus=="Oreochromis",identifiedBy=="mtDNA",!is.na(specificEpithet)) |>
    dplyr::left_join(distinct(events.master,eventID,localityID,localitySite,nickname),by="eventID") |>
    dplyr::mutate(localityID=str_replace_all(localityID,"MG","Site"),localitySite=stringr::str_replace_all(localitySite,"MG","Site")) |>
    dplyr::mutate(localityID=paste(localityID,nickname,sep="-"),localitySite=stringr::str_replace_all(localitySite,"-",glue::glue("-{nickname}-")))

# tissues per pond
tissues.master.format |> 
    dplyr::count(localitySite) |> 
    knitr::kable(caption="Number of tissue samples from each pond",format="simple")
#tissues.master.format |> dplyr::group_by(localitySite) |> dplyr::summarise(n=n()) |> dplyr::pull(n) |> mean()

# tissues per species
tissues.master.format |> 
    dplyr::distinct(specificEpithet,localitySite) |> 
    dplyr::count(specificEpithet) |> 
    knitr::kable(caption="Number of tissue samples by species",format="simple")

# by site
tissues.master.format |> 
    summarise_tissues(grouping=localityID) |>
    knitr::kable(caption="Presence/absence of species by site",format="simple")

# by pond
tissues.master.format |> 
    summarise_tissues(grouping=localitySite) |>
    knitr::kable(caption="Presence/absence of species by pond",format="simple")

# by pond
tissues.binary <- tissues.master.format |> summarise_tissues(grouping=localitySite)

# add tissues data to SUMMARISED
qcr.results.summarised.tissues <- qcr.results.summarised |> 
    dplyr::left_join(pivot_longer(tissues.binary,cols=c(niloticus,urolepis,leucostictus),names_to="species",values_to="presence"),by=c("group","species")) |>
    dplyr::group_by(group) |>
    dplyr::mutate(nPcrRepsPositiveByGroup=sum(nPcrRepsPositive)) |> 
    dplyr::ungroup()

# add tissues data to BY EXTRACTION
qcr.results.by.extraction.tissues <- qcr.results.by.extraction |> 
    select(localityID,localitySite,assay,extractionTubeLabel,copiesLitre) |> 
    rename(group=localitySite,site=localityID,species=assay) |>
    left_join(tidyr::pivot_longer(tissues.binary,cols=c(niloticus,urolepis,leucostictus),names_to="species",values_to="presence"),by=c("group","species")) |>
    rename(pond=group,assay=species,fishPresence=presence)

# write out
qcr.results.by.extraction.tissues |> readr::write_csv(here::here("temp/qpcr-results-processed.csv"))

###########################
##### SAMPLING EVENTS #####
###########################

report("/// SAMPLING SITES ///",type="i")

# print out events
events.master |> 
    dplyr::distinct(localitySite,localityID,nickname,fieldID,decimalLatitude,decimalLongitude,year,month,day,waterBody) |>
    dplyr::mutate(group=stringr::str_replace_all(localitySite,"MG","Site"),group=stringr::str_replace_all(group,"-",glue::glue("-{nickname}-"))) |>
    dplyr::filter(group %in% pull(qcr.results.names,localitySite)) |>
    tidyr::separate(col=group,into=c("site","name","pond"),sep="-",remove=TRUE) |>
    dplyr::mutate(samplingDate=ymd(paste(year,month,day,sep="-"))) |>
    dplyr::mutate(decimalLatLon=paste(round(decimalLatitude,digits=3),round(decimalLongitude,digits=3),sep=", ")) |>
    dplyr::arrange(site,name,pond) |>
    dplyr::select(site,name,pond,fieldID,waterBody,samplingDate,decimalLatLon) |>
    dplyr::mutate(samplingDate=as.character(samplingDate)) |>
    dplyr::rename(Site=site,Name=name,Pond=pond,`Field ID code`=fieldID,Drainage=waterBody,`Sampling date (ymd)`=samplingDate,`Location (decimal lat/lon)`=decimalLatLon) |>
    knitr::kable(caption="Summary of sampling sites",format="simple")

# get status
#qcr.results.summarised.tissues |> dplyr::group_by(species) |> dplyr::summarise(n=length(which(copiesTotal>0)),tot=n(),prop=n/tot)
#qcr.results.summarised.tissues |> dplyr::group_by(species) |> dplyr::summarise(n=length(which(proportion>0.1)),tot=n(),prop=n/tot)
#qcr.results.summarised.tissues |> dplyr::group_by(species) |> dplyr::summarise(n=length(which(proportion>0.9)),tot=n(),prop=n/tot)
#qcr.results.summarised.tissues |> dplyr::group_by(group) |> dplyr::summarise(n=length(which(copiesTotal>0))) |> dplyr::group_by(n) |> dplyr::summarise(count=n()) |> dplyr::mutate(tot=sum(count),prop=count/tot)
#qcr.results.summarised.tissues |> dplyr::filter(presence==1 & copiesTotal==0)
#qcr.results.summarised.tissues |> dplyr::filter(presence==1 & copiesTotal>0) |> print(n=Inf)
#

# report
report("/// RESULTS SUMMARY TABLES ///",type="i")

# print results for each pond
qcr.results.summarised.tissues |> 
    dplyr::mutate(filterPositive=paste(nFilterRepsPositive,nFilterReps,sep="/"),pcrPositive=paste(nPcrRepsPositive,nPcrReps,sep="/")) |>
    dplyr::group_by(group) |>
        dplyr::mutate(maxProportion=max(proportion),dominantSpecies=if_else(maxProportion>0,species[which(proportion==maxProportion)],"NA")) |> #
        dplyr::ungroup() |>
    dplyr::mutate(copiesTotal=round(copiesTotal,digits=0),proportion=round(proportion,digits=4)) |>
    tidyr::separate(col=group,into=c("site","name","pond"),sep="-",remove=TRUE) |>
    dplyr::select(-dominantSpecies,-maxProportion,-nPcrRepsPositiveByGroup,-nFilterReps,-nFilterRepsPositive,-nPcrReps,-nPcrRepsPositive) |>
    dplyr::arrange(site,name,pond,desc(proportion)) |>
    #write_csv(here("temp/qcr-results-summarised.csv"))
    dplyr::mutate(proportion=round(proportion*100,digits=1)) |>
    dplyr::mutate(presence=dplyr::if_else(is.na(presence),"",as.character(presence))) |>
    #mutate(presence=as.character(presence),presence=str_replace_all(presence,NA,"")) |>
    dplyr::rename(Site=site,Name=name,Pond=pond,Species=species,`Total copies`=copiesTotal,`Proportion (%)`=proportion,`Species presence`=presence,`Positive filters (n/total)`=filterPositive,`Positive PCRs (n/total)`=pcrPositive) |>
    knitr::kable(caption="Summary of results by pond",format="simple",format.args=list(big.mark=","),align="llllrrrrr")


###############
##### END #####
###############

report("All steps complete; processed qPCR results written to {.file temp/qpcr-results-processed.csv}",type="s")
