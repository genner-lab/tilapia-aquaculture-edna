#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################

library("here")
library("glue")
library("knitr")
library("tidyverse")
library("lubridate")
library("xtable")
# plot mixed models 
library("lme4")
library("glmmTMB")
library("sjPlot")
library("emmeans")
library("performance")
library("ggthemes")
#rm(list=ls())
options(width=180)

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


######################
##### CHECK DATA #####
######################

# check no template controls
qcr.results.ntc <- qcr.results |> 
    dplyr::filter(role=="NTC") |> 
    dplyr::arrange(desc(copies))
qcr.results.ntc |> dplyr::glimpse()
qcr.results.ntc |> dplyr::filter(!is.na(copies)) |> dplyr::glimpse()

# check standards
qcr.results |> 
    dplyr::filter(is.na(exclude)) |> 
    dplyr::filter(role=="Standard") |> 
    dplyr::group_by(assay,copies) |>
    dplyr::summarise(meancq=mean(cq),n=n(),.groups="drop") |>
    dplyr::arrange(assay,desc(copies)) |>
    print(n=Inf)
    #knitr::kable()


######################
##### CLEAN DATA #####
######################

# remove positive controls and standards
# change copies to zero for excluded
# split into dilution factor
qcr.results.clean <- qcr.results |> 
    dplyr::filter(role=="Unknown") |>
    dplyr::select(-role,-thresholdSetting,-thresholdCq,-baselineSetting,-baselineStart,-baselineEnd) |>
    dplyr::mutate(copies=dplyr::if_else(is.na(copies),0,copies)) |>
    tidyr::separate(col=sample,into=c("sample","dilutionFactor"),sep=" \\(",remove=FALSE) |>
    dplyr::mutate(dilutionFactor=as.numeric(str_split_fixed(dilutionFactor,"",4)[,3])) |>
    dplyr::rename(extractionTubeLabel=sample)

# check for missing samples
# (should be empty)
qcr.results.clean |> dplyr::filter(!extractionTubeLabel %in% dplyr::pull(extractions.master,extractionTubeLabel))
# (should all have "no sample" in notes)
extractions.master |> dplyr::filter(!extractionTubeLabel %in% dplyr::pull(qcr.results.clean,extractionTubeLabel)) |> print(width=200)

# join with the sample info
qcr.results.ext <- qcr.results.clean |> dplyr::left_join(extractions.master,by="extractionTubeLabel")

# blanks (should all be zero copies)
qcr.results.ext |> 
    dplyr::filter(!is.na(blankType)) |> 
    dplyr::filter(copies>0) 

# check event ids (should all be blanks)
qcr.results.ext |> dplyr::filter(!eventID %in% dplyr::pull(events.master,eventID)) |> print(n=Inf)
# check event ids (should all be KIL)
events.master |> dplyr::filter(!eventID %in% dplyr::pull(qcr.results.ext,eventID)) |> print(n=Inf)

# get numbers for controls
qcr.results.ext |> dplyr::filter(blankType=="FIELD") |> dplyr::distinct(extractionTubeLabel)
qcr.results.ext |> dplyr::filter(blankType=="FIELD") |> dplyr::glimpse()
qcr.results.ext |> dplyr::filter(blankType=="LAB") |> dplyr::distinct(extractionTubeLabel)
qcr.results.ext |> dplyr::filter(blankType=="LAB") |> dplyr::glimpse()

# add events data
# clean up - change names, remove blanks and superfluous columns
qcr.results.joined <- qcr.results.ext |> 
    dplyr::left_join(events.master,by="eventID") |>
    dplyr::filter(is.na(blankType)) |>
    dplyr::select(eventID,localityID,localitySite,,nickname,assay,extractionTubeLabel,extractionDate,dilutionFactor,repVol,copies) |>
    dplyr::mutate(localityID=stringr::str_replace_all(localityID,"MG","Site"),localitySite=stringr::str_replace_all(localitySite,"MG","Site")) |>
    dplyr::mutate(localityID=glue::glue("{localityID}-{nickname}"), localitySite=stringr::str_replace_all(localitySite,"-",glue::glue("-{nickname}-")))

qcr.results.joined |> filter(extractionTubeLabel=="10.06.19/01")

# get df of names 
qcr.results.names <- qcr.results.joined |> dplyr::distinct(localityID,localitySite,extractionTubeLabel) 
qcr.results.names |> dplyr::distinct(localityID) |> dplyr::arrange(localityID) |> print(n=Inf) # Site05 and Site13 has two names
qcr.results.names |> dplyr::distinct(localitySite) |> dplyr::arrange(localitySite) |> print(n=Inf) 
qcr.results.names |> dplyr::distinct(extractionTubeLabel) |> dplyr::arrange(extractionTubeLabel) |> print(n=Inf) 


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
qcr.results.by.extraction 

#########################
##### SUMMARY STATS #####
#########################


# function to summarise by grouping var ()
summarise_qpcr <- function(df,grouping) {
    df.prop <- df |> 
        dplyr::group_by({{grouping}}) |>
        dplyr::mutate(copiesLitreTotal=sum(copiesLitre)) |>
        dplyr::ungroup() |>
        dplyr::group_by(assay,{{grouping}},copiesLitreTotal) |>
        dplyr::summarise(copiesLitreSpecies=sum(copiesLitre),
            nFilterReps=length(unique(extractionTubeLabel)),
            nFilterRepsPositive=length(unique(extractionTubeLabel[nPcrRepsPositive>0])),
            nPcrReps=n(),
            nPcrRepsPositive=length(which(copiesLitre>0)),
            .groups="drop") |>
        dplyr::mutate(proportion=copiesLitreSpecies/copiesLitreTotal) |> 
        dplyr::arrange({{grouping}},assay) |>
        dplyr::select({{grouping}},assay,copiesLitreSpecies,proportion,nFilterReps,nFilterRepsPositive,nPcrReps,nPcrRepsPositive) |>
        dplyr::rename(group={{grouping}},species=assay,copiesTotal=copiesLitreSpecies) |>
        dplyr::mutate(proportion=if_else(is.nan(proportion),0,proportion))
return(df.prop)
} 


# summarise by 'localityID'
#qcr.results.summarised <- qcr.results.by.extraction |> summarise_qpcr(grouping=localityID)#  |> filter(group=="Site02-Msauzi")
# summarise by 'localitySite'
qcr.results.summarised <- qcr.results.by.extraction |> summarise_qpcr(grouping=localitySite)# |> filter(group=="Site02-Msauzi-A")
# summarise by 'extractionTubeLabel'
#qcr.results.summarised <- qcr.results.by.extraction |> summarise_qpcr(grouping=extractionTubeLabel) |> 
#    left_join(rename(qcr.results.names,group=extractionTubeLabel)) |> 
#    mutate(group=localityID) |>
#    #mutate(group=localitySite) |>
#    select(-localityID,-localitySite)


# format the tissues
tissues.master.format <- tissues.master |> 
    dplyr::filter(eventID %in% dplyr::pull(qcr.results.ext,eventID),genus=="Oreochromis",identifiedBy=="mtDNA",!is.na(specificEpithet)) |>
    dplyr::left_join(distinct(events.master,eventID,localityID,localitySite,nickname),by="eventID") |>
    dplyr::mutate(localityID=str_replace_all(localityID,"MG","Site"),localitySite=stringr::str_replace_all(localitySite,"MG","Site")) |>
    dplyr::mutate(localityID=paste(localityID,nickname,sep="-"),localitySite=stringr::str_replace_all(localitySite,"-",glue::glue("-{nickname}-")))

# numbers
tissues.master.format |> dplyr::glimpse()
tissues.master.format |> dplyr::group_by(localitySite) |> dplyr::summarise(n=n()) 
tissues.master.format |> dplyr::group_by(localitySite) |> dplyr::summarise(n=n()) |> dplyr::pull(n) |> mean()
tissues.master.format |> dplyr::distinct(specificEpithet,localitySite) |> dplyr::group_by(specificEpithet) |> dplyr::summarise(n=n())
tissues.master.format |> dplyr::distinct(localitySite) 


# summarise tissues
summarise_tissues <- function(df,grouping) {
    df.bin <- df |> 
        dplyr::select({{grouping}},specificEpithet) |>
        dplyr::group_by({{grouping}},specificEpithet) |>
        dplyr::count() |> 
        tidyr::pivot_wider(names_from=specificEpithet,values_from=n,values_fill=0) |>
        dplyr::mutate(across(c(niloticus,urolepis,leucostictus),~if_else(.x>0,1,0))) |>
        dplyr::ungroup() |>
        dplyr::rename(group={{grouping}})
return(df.bin)
}

# run for groupvar
#tissues.binary <- tissues.master.format |> summarise_tissues(grouping=localityID)
tissues.binary <- tissues.master.format |> summarise_tissues(grouping=localitySite)


# add tissues data to SUMMARISED
qcr.results.summarised.tissues <- qcr.results.summarised |> 
    dplyr::left_join(pivot_longer(tissues.binary,cols=c(niloticus,urolepis,leucostictus),names_to="species",values_to="presence"),by=c("group","species")) |>
    dplyr::group_by(group) |>
    dplyr::mutate(nPcrRepsPositiveByGroup=sum(nPcrRepsPositive)) |> 
    dplyr::ungroup()

# add tissues data to BY EXTRACTION
#qcr.results.summarised.tissues <- qcr.results.by.extraction |> 
#    select(localityID,localitySite,extractionTubeLabel,assay,copiesLitre,nPcrRepsPositive) |> 
#    rename(group=localitySite,site=localityID,species=assay) |>
#    left_join(pivot_longer(tissues.binary,cols=c(niloticus,urolepis,leucostictus),names_to="species",values_to="presence"),by=c("group","species")) |>
#    group_by(group) |>
#    mutate(nPcrRepsPositiveByGroup=sum(nPcrRepsPositive)) |> 
#    ungroup()


# format the events
events.master |> 
    dplyr::distinct(localitySite,localityID,nickname,fieldID,decimalLatitude,decimalLongitude,year,month,day,waterBody) |>
    dplyr::mutate(group=str_replace_all(localitySite,"MG","Site"),group=str_replace_all(group,"-",glue::glue("-{nickname}-"))) |>
    dplyr::filter(group %in% pull(qcr.results.names,localitySite)) |>
    tidyr::separate(col=group,into=c("site","name","pond"),sep="-",remove=TRUE) |>
    dplyr::mutate(samplingDate=ymd(paste(year,month,day,sep="-"))) |>
    dplyr::mutate(decimalLatLon=paste(round(decimalLatitude,digits=3),round(decimalLongitude,digits=3),sep=", ")) |>
    #group_by(site,name) |>
    #    mutate(numberPonds=n()) |>
    #    ungroup() |>
    dplyr::arrange(site,name,pond) |>
    dplyr::select(site,name,pond,fieldID,waterBody,samplingDate,decimalLatLon) |>
    #write_csv(here("temp/sampling-sites-summarised.csv"))
    dplyr::mutate(samplingDate=as.character(samplingDate)) |>
    dplyr::rename(Site=site,Name=name,Pond=pond,`Field ID code`=fieldID,Drainage=waterBody,`Sampling date (ymd)`=samplingDate,`Location (decimal lat/lon)`=decimalLatLon) |>
    xtable::xtable(caption="blahhh", digits=rep(0,8)) |>
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="tiny")

# get status
qcr.results.summarised.tissues |> dplyr::group_by(species) |> dplyr::summarise(n=length(which(copiesTotal>0)),tot=n(),prop=n/tot)
qcr.results.summarised.tissues |> dplyr::group_by(species) |> dplyr::summarise(n=length(which(proportion>0.1)),tot=n(),prop=n/tot)
qcr.results.summarised.tissues |> dplyr::group_by(species) |> dplyr::summarise(n=length(which(proportion>0.9)),tot=n(),prop=n/tot)
qcr.results.summarised.tissues |> dplyr::group_by(group) |> dplyr::summarise(n=length(which(copiesTotal>0))) |> dplyr::group_by(n) |> dplyr::summarise(count=n()) |> dplyr::mutate(tot=sum(count),prop=count/tot)
qcr.results.summarised.tissues |> dplyr::filter(presence==1 & copiesTotal==0)
qcr.results.summarised.tissues |> dplyr::filter(presence==1 & copiesTotal>0) |> print(n=Inf)


# write out the results for Appendix
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
    dplyr::mutate(across(where(is.numeric),~prettyNum(.x,big.mark=","))) |>
    xtable::xtable(caption="blahhh", digits=rep(0,10)) |>
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="tiny")

############################################### 
###############################################

# plot logistic glm
qcr.results.summarised.tissues |> 
    dplyr::filter(nPcrRepsPositiveByGroup>0 & !is.na(presence)) |>
    #ggplot(aes(x=proportion,y=presence)) + 
    ggplot2::ggplot(aes(x=log(copiesTotal+1),y=presence)) + 
        geom_point() + 
        stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE) +
        facet_wrap(vars(species),scales="free")

# check outliers
qcr.results.summarised.tissues |> filter(presence==1 & proportion<0.05)
qcr.results.summarised.tissues |> filter(presence==0 & proportion>0.25)


# format for mixed model 
qcr.results.summarised.tissues.reduced <- qcr.results.summarised.tissues |> 
    filter(nPcrRepsPositiveByGroup>0 & !is.na(presence)) |> 
    mutate(presence=as.factor(presence)) |>
    mutate(site=str_split_fixed(group,"-",3)[,1])

# run a glmm
m0 <- glmmTMB(
    formula=log10(copiesTotal+1) ~ presence + species + (1|site/group), 
    data=qcr.results.summarised.tissues.reduced, 
    family=gaussian(link="identity"),
    REML=TRUE,
    na.action=na.fail)

# get summary
summary(m0)

# run lmm
m0 <- lmer(log10(copiesTotal+1) ~ presence + species + (1|site/group), data = qcr.results.summarised.tissues.reduced, na.action="na.omit")
summary(m0)
sjPlot::plot_model(m0)
sjPlot::tab_model(m0)

# check model stats
performance::check_model(m0)

# estimated marginal means
m0.means <- emmeans(m0,specs="presence")
#m0.means |> plot()
pairs(m0.means,type="response",reverse=TRUE)
# contrast ratio    SE df null t.ratio p.value
# 1 / 0     8484 13401 39    1   5.727  <.0001

# plot with emmeans
p <- qcr.results.summarised.tissues.reduced |> 
    ggplot(aes(y=log10(copiesTotal+1),x=presence)) + 
    #geom_point(alpha=0.5,pch=21) + 
    geom_jitter(alpha=0.5,pch=21,width=0.05,height=0) + 
    geom_pointrange(data=as_tibble(m0.means),mapping=aes(x=presence,y=emmean, ymin=lower.CL, ymax=upper.CL),color="orange",size=1) +
    theme_clean() +
    labs(x="Species presence",y="Total copies (log10 transformed)")
plot(p)
# save
ggsave(filename=here("temp/maps/copies-by-presence.svg"),plot=p,width=120,height=120,units="mm")
