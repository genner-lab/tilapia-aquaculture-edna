#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################

library("ggplot2")
library("ggpubr")
library("dplyr")
library("car")
library("emmeans")
library("lme4")
library("lmerTest")
library("DHARMa")

ModelData_All <- read.table(here::here("temp/ModelData_19Oct2024.tsv"),header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

### First estimate proportions of each species across all collected sites

ModelData_All_FC_Summary <- ModelData_All %>%
  group_by(eventID, assay, extractionTubeLabel) %>%
  reframe(copies_l_water = mean(copies_l_water), Target_present_in_fish = max(Target_present_in_fish),fish_collected = max(fish_collected),n = n())

ModelData_All_FC_Proportion <- ModelData_All_FC_Summary %>%
  group_by(extractionTubeLabel) %>%
  mutate(proportion =  copies_l_water/sum(copies_l_water))

ModelData_All_FC_Proportion_Site <- ModelData_All_FC_Proportion %>%
  group_by(eventID, assay) %>%
  reframe(copies_l_water = mean(copies_l_water))

write.table (ModelData_All_FC_Proportion_Site, file="CopyNumberPerSite.txt")

### Next filter to limit only to the 15 fish collected sites

ModelData_All_FC_Proportion <- filter(ModelData_All_FC_Proportion, fish_collected == "Yes")

### Clean house

rm(ModelData_All)
rm(ModelData_All_FC_Proportion_Site)
rm(ModelData_All_FC_Summary)

### Model the quantities of DNA across species

ModelData_All_FC_Proportion$copies_l_water_4th <- sqrt(sqrt(ModelData_All_FC_Proportion$copies_l_water))

ModelAll1 <- lmer(copies_l_water_4th ~ Target_present_in_fish + assay + (1|eventID),
                      data = ModelData_All_FC_Proportion)
summary(ModelAll1)
anova(ModelAll1)

emmeans(ModelAll1, specs = pairwise ~ Target_present_in_fish)

(-0.611^4)
(0.912^4)
(2.44^4)

(3.628^4)
(5.149^4)
(6.67^4)

### Model the quantities of DNA and proportion of DNA comprised by each species

ModelData_Leuco_Proportion <- ModelData_All_FC_Proportion[which(ModelData_All_FC_Proportion$assay=='leucostictus'), ]

ModelLeuco1 <- lmer(copies_l_water_4th ~ Target_present_in_fish + (1|eventID),
                   data = ModelData_Leuco_Proportion)
summary(ModelLeuco1)
anova(ModelLeuco1)

ModelLeuco2 <- lmer(asin(proportion) ~ Target_present_in_fish + (1|eventID),
                        data = ModelData_Leuco_Proportion)
summary(ModelLeuco2)
anova(ModelLeuco2)

###

ModelData_Nilo_Proportion <- ModelData_All_FC_Proportion[which(ModelData_All_FC_Proportion$assay=='niloticus'), ]

ModelNilo1 <- lmer(copies_l_water_4th ~ Target_present_in_fish + (1|eventID),
                   data = ModelData_Nilo_Proportion)
summary(ModelNilo1)
anova(ModelNilo1)

ModelNilo2 <- lmer(asin(proportion) ~ Target_present_in_fish + (1|eventID),
                        data = ModelData_Nilo_Proportion)
summary(ModelNilo2)
anova(ModelNilo2)

###

ModelData_Uro_Proportion <- ModelData_All_FC_Proportion[which(ModelData_All_FC_Proportion$assay=='urolepis'), ]

ModelUro1 <- lmer(copies_l_water_4th ~ Target_present_in_fish + (1|eventID),
                  data = ModelData_Uro_Proportion)
summary(ModelUro1)
anova(ModelUro1)

ModelUro2 <- lmer(asin(proportion) ~ Target_present_in_fish + (1|eventID),
                       data = ModelData_Uro_Proportion)
summary(ModelUro2)
anova(ModelUro2)

### Plot boxplots

ModelData_All_FC_Proportion$Target_present_in_fish<-as.factor(ModelData_All_FC_Proportion$Target_present_in_fish)

plot1 <- ggplot(ModelData_All_FC_Proportion, aes(y=copies_l_water_4th, x=Target_present_in_fish)) + 
  geom_boxplot(aes(fill = assay),outliers = FALSE) + 
  geom_jitter(height=0,width=0.3,shape=1) +
  facet_grid(~factor(assay, levels=c('niloticus', 'urolepis', 'leucostictus'))) +
  theme_classic() +
  scale_fill_manual(values=c("gold", "steelblue2", "red")) + # set the color of the values manualy
  labs(x ="Confirmed species presence in water body", y =  bquote('eDNA copies (4th root transformed)'),) +
  theme(legend.position = "right")
plot1

ModelData_All_FC_Proportion2 <- ModelData_All_FC_Proportion[complete.cases(ModelData_All_FC_Proportion), ]

plot2 <- ggplot(ModelData_All_FC_Proportion2, aes(y=proportion, x=Target_present_in_fish)) + 
  geom_boxplot(aes(fill = assay),outliers = FALSE) + 
  geom_jitter(height=0,width=0.3,shape=1) +
  facet_grid(~factor(assay, levels=c('niloticus', 'urolepis', 'leucostictus'))) +
  theme_classic() +
  scale_fill_manual(values=c("gold", "steelblue2", "red")) + # set the color of the values manualy
  labs(x ="Confirmed species presence in water body", y =  bquote('Proportion of eDNA in water body'),) +
  theme(legend.position = "right")
plot2

figure <- ggarrange (plot1, plot2, ncol = 1, nrow =2, legend = "none")
figure

#save as 10x6 landscape

#######

