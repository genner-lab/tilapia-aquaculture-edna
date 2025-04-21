#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################

# load libs and funs
source(here::here("scripts/load-libs.R"))
report("R packages loaded",type="s")

# load data
model.data <- readr::read_csv(here::here("temp/qpcr-results-processed.csv"),show_col_types=FALSE)
report("Data loaded from {.file temp/qpcr-results-processed.csv}",type="s")

### summarise, format and transform data

model.data.summary <- model.data |> 
    dplyr::group_by(pond,assay,extractionTubeLabel) |>
    dplyr::summarise(copiesLitre=mean(copiesLitre),fishPresence=max(fishPresence),.groups="drop") |> 
    dplyr::group_by(extractionTubeLabel) |> 
    dplyr::mutate(copiesLitreExtraction=sum(copiesLitre)) |>
    dplyr::ungroup() |>
    dplyr::mutate(copiesLitreTrans=sqrt(sqrt(copiesLitre))) |>
    dplyr::mutate(proportion=copiesLitre/copiesLitreExtraction) |>
    dplyr::mutate(proportionTrans=asin(sqrt(proportion)))

### filter to limit only to the 15 fish collected sites
model.data.summary.proportion <- model.data.summary |> dplyr::filter(!is.na(fishPresence))


#####################################################
report("/// LINEAR MIXED EFFECTS MODEL ///",type="i")
#####################################################

### model the quantities of DNA across species with linear mixed effects model

mod.all <- lmerTest::lmer(copiesLitreTrans ~ factor(fishPresence) + factor(assay) + (1|pond), data = model.data.summary.proportion)

# anova table
mod.all |> anova() |> parameters::model_parameters()

# check model
dont_print({
    mod.all |> performance::check_model()
})

# model params
mod.all |> parameters::model_parameters()

# get marginal means
mod.all |> modelbased::estimate_means(c("assay","fishPresence"),transform=power4)
mod.all |> modelbased::estimate_means("fishPresence",transform=power4)
mod.all |> modelbased::estimate_contrasts("fishPresence",transform=power4)


### model the quantities of DNA and proportion of DNA comprised by each species

####################################################
report("/// LEUCOSTICTUS ASSAY COPIES ///",type="i")
####################################################

# leucostictus
model.data.summary.proportion.leuco <- model.data.summary.proportion |> dplyr::filter(assay=="leucostictus")

# copies
mod.leuco <- lmerTest::lmer(copiesLitreTrans ~ factor(fishPresence) + (1|pond), data = model.data.summary.proportion.leuco)
mod.leuco |> anova() |> parameters::model_parameters()
mod.leuco |> modelbased::estimate_means("fishPresence",transform=power4)


#########################################################
report("/// LEUCOSTICTUS ASSAY PROPORTIONS ///",type="i")
#########################################################

# leuco proportions
mod.leuco.prop <- lmerTest::lmer(proportionTrans ~ factor(fishPresence) + (1|pond), data = model.data.summary.proportion.leuco)
mod.leuco.prop |> anova() |> parameters::model_parameters()
mod.leuco.prop |> modelbased::estimate_means("fishPresence",transform=sqsine)


#################################################
report("/// NILOTICUS ASSAY COPIES ///",type="i")
#################################################

# niloticus
model.data.summary.proportion.nilo <- model.data.summary.proportion |> dplyr::filter(assay=="niloticus")

# copies
mod.nilo <- lmerTest::lmer(copiesLitreTrans ~ factor(fishPresence) + (1|pond), data = model.data.summary.proportion.nilo)
mod.nilo |> anova() |> parameters::model_parameters()
mod.nilo |> modelbased::estimate_means("fishPresence",transform=power4)


######################################################
report("/// NILOTICUS ASSAY PROPORTIONS ///",type="i")
######################################################

# nilo proportions
mod.nilo.prop <- lmerTest::lmer(proportionTrans ~ factor(fishPresence) + (1|pond), data = model.data.summary.proportion.nilo)
mod.nilo.prop |> anova() |> parameters::model_parameters()
mod.nilo.prop |> modelbased::estimate_means("fishPresence",transform=sqsine)


################################################
report("/// UROLEPIS ASSAY COPIES ///",type="i")
################################################

# urolepis
model.data.summary.proportion.uro <- model.data.summary.proportion |> dplyr::filter(assay=="urolepis")

# copies
mod.uro <- lmerTest::lmer(copiesLitreTrans ~ factor(fishPresence) + (1|pond), data = model.data.summary.proportion.uro)
mod.uro |> anova() |> parameters::model_parameters()
mod.uro |> modelbased::estimate_means("fishPresence",transform=power4)


#####################################################
report("/// UROLEPIS ASSAY PROPORTIONS ///",type="i")
#####################################################

# uro proportions
mod.uro.prop <- lmerTest::lmer(proportionTrans ~ factor(fishPresence) + (1|pond), data = model.data.summary.proportion.uro)
mod.uro.prop |> anova() |> parameters::model_parameters()
mod.uro.prop |> modelbased::estimate_means("fishPresence",transform=sqsine)


############################################
report("/// PLOTTING BOXPLOTS ///",type="i")
############################################

# plot 1
plot1 <- model.data.summary.proportion |> 
    ggplot2::ggplot(aes(y=copiesLitreTrans, x=factor(fishPresence))) + 
        ggplot2::geom_boxplot(aes(fill=assay),outliers=FALSE) + 
        ggplot2::geom_jitter(position=ggplot2::position_jitter(height=0,width=0.25,seed=42),shape=1,alpha=0.5) +
        ggplot2::facet_grid(~factor(assay,levels=c("niloticus","urolepis","leucostictus"))) +
        ggplot2::scale_fill_manual(values=c("gold","steelblue2","red")) +
        ggthemes::theme_clean() +
        ggplot2::theme(legend.position="none") +
        ggplot2::labs(x="Confirmed species presence in water body", y="eDNA copies (4th root transformed)")
#print(plot1)
ggplot2::ggsave(plot=plot1,filename=here::here("temp/plot1.svg"),device=svglite::svglite,width=8,height=3,system_fonts=list(sans="Arial"))
# to fix text editing in Inkscape
# sed -i -E "s/ (textLength|lengthAdjust)='[^']*'//g" plot1.svg
report("Plot saved to {.file temp/plot1.svg}",type="s")

# plot2
plot2 <- model.data.summary.proportion |> 
    dplyr::filter(!is.nan(proportion)) |> 
    ggplot2::ggplot(aes(y=proportion, x=factor(fishPresence))) + 
        ggplot2::geom_boxplot(aes(fill=assay),outliers=FALSE) + 
        ggplot2::geom_jitter(position=ggplot2::position_jitter(height=0,width=0.25,seed=4242),shape=1,alpha=0.5) +
        ggplot2::facet_grid(~factor(assay,levels=c("niloticus","urolepis","leucostictus"))) +
        ggplot2::scale_fill_manual(values=c("gold","steelblue2","red")) +
        ggthemes::theme_clean() +
        ggplot2::theme(legend.position="none") +
        ggplot2::labs(x="Confirmed species presence in water body", y="Proportion of eDNA in water body")
#print(plot2)
ggplot2::ggsave(plot=plot2,filename=here::here("temp/plot2.svg"),device=svglite::svglite,width=8,height=3,system_fonts=list(sans="Arial"))
# to fix text editing in Inkscape
# sed -i -E "s/ (textLength|lengthAdjust)='[^']*'//g" plot2.svg
report("Plot saved to {.file temp/plot2.svg}",type="s")

# fin
report("All steps complete",type="s")
