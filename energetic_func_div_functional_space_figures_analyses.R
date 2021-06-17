#Energetic Functional Space Expansion
#Manuscript: "Functional space expansion driven by transitions between advantageous traits across a deep-sea energetic gradient"
#Authors: S. River D. Bryant & Craig R. McClain
##Functional Diversity Metric Calculations, Plots, & Analyses
##Code Author: River Dixon Bryant

################### SETUP ###################

#LOAD PACKAGES
library(dplyr)      #for data manipulation
library(tidyr)      #for data manipulation
library(ggplot2)    #for plotting
library(gridExtra)  #for arranging plots into figure pdfs
library(multirich)  #for UTC, sUTC, raw overlap, simple overlap, coverage overlap, mean, median, max, and min overlap 
library(vegan)      #for diversity metrics
library(FD)         #for functional diversity metrics
library(corrplot)   #for creating correlation plot
library(stargazer)  #for creating output tables
library(car)        #for ANOVAs
library(rcompanion) #for pseudo r-square
library(voxel)      #for plotgam
library(visreg)     #for visualizaton of models
library(mgcv)       #for general additive model (GAMM)
library(mgcViz)     #for GAM Plots
library(grid)       #for GAM Plots
library(ggplotify)  #for GAM Plots
library(cowplot)    #for GAM Plots
library(ggpubr)

##ggplot stuff
theme_craig <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      # change stuff here
      axis.line = element_line(colour = "#6B7D83"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      legend.position="none",
      plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
}


################### DATA MANIPULATION FOR ANALYSES ##############################

#CREATE DATA HOLDER FOR SPECIES X STATION MATRIX
biv_taxa = read.csv("allen_bivalve_sp.station.csv", header=TRUE)
row.names(biv_taxa) = biv_taxa$Station


#MAKE IT PRESENCE/ABSENCE
biv_taxa_binary = +(biv_taxa >= 1)

#FIRST ROW IS REGIONAL POOL
regional = matrix(1, ncol = ncol(biv_taxa_binary))

biv_taxa_binary.regional = rbind(regional, biv_taxa_binary)
row.names(biv_taxa_binary.regional)[1] = "Regional"

#

#CREATE MAIN STATION ENVIRONMENTAL DATA HOLDER
biv_enviro = read.csv("allen_bivalve_station.csv")

row.names(biv_enviro) = biv_enviro$Station2
biv_enviro = biv_enviro[,-c(10)]
glimpse(biv_enviro)

#

#CREATE MAIN TRAIT DATA HOLDER
biv_trait_orig = read.csv("allen_bivalve_species.traits.csv", header=TRUE)
glimpse(biv_trait_orig)

biv_trait = biv_trait_orig %>% 
  select(tiering:ecospace)
row.names(biv_trait) = biv_trait_orig$Species_ID

#

#CULL TRAIT & BINARY TAXA DATASETS DOWN
biv_trait2 = na.omit(biv_trait)

keepspecies = colnames(biv_taxa_binary.regional) [(colnames(biv_taxa_binary.regional) 
                                                   %in% row.names(biv_trait2))] 

biv_taxa_binary.regional2 = biv_taxa_binary.regional[, keepspecies]

dim(biv_trait2)
dim(biv_taxa_binary.regional2)

#

#DELETE STATIONS WITH < 3 SPECIES
#CREATE THRESHOLD FOR STATIONS WITH MORE THAN THREE SPECIES PRESENT
threshhold = which(rowSums(biv_taxa_binary.regional2) > 3)

#cREATE NEW MATRIX WITH ONLY THOSE STATIONS WHICH MEET THRESHOLD
biv_taxa_binary.regional3 = biv_taxa_binary.regional2[threshhold, ]

#CHECK IT WORKS
dim(biv_taxa_binary.regional2)
dim(biv_taxa_binary.regional3)

################### Calculating UTC and Functional Overlap #########################

#UTC, sUTC, raw overlap, simple overlap, coverage overlap, mean, median, max, and min overlap

mvr.biv = mvfd(as.matrix(biv_trait2), as.matrix(biv_taxa_binary.regional3))  
summary(mvr.biv)

#PUT THE RESULTS INTO A NEW DATA FRAME
biv_results = data.frame(rownames(biv_taxa_binary.regional3)) 
colnames(biv_results)[1] <- "Station"

biv_results$utc = mvr.biv$utc
biv_results$sutc = mvr.biv$sutc
biv_results$rmvo = mvr.biv$rmvo
biv_results$smvo = mvr.biv$smvo
biv_results$cmvo = mvr.biv$cmvo
biv_results$meanmvo = mvr.biv$meanmvo
biv_results$medmvo = mvr.biv$medmvo
biv_results$maxmvo = mvr.biv$maxmvo
biv_results$minmvo = mvr.biv$minmvo

biv_results = biv_results [2:224,] #delete regional row

#CULL ORIGINAL TAXA DATA DOWN
biv_taxa2 = biv_taxa [, keepspecies] #take out species we don't have trait info for
biv_taxa2_alpha = biv_taxa2 [, order(names(biv_taxa2))] #alphabetize columns so they match order of traits

#

#DROP SPECIES THAT DO NOT OCCUR AT ANY STATION
biv_taxa2_alpha$Arc48 = NULL #drop Arc48 because does not occur at any station
biv_trait2 = biv_trait2[-c(5),] #drop Arc48 from the trait data too

#get information about how many species with each trait
taxa_summary_tier = biv_trait2 %>% group_by(tiering) %>% summarise(n = n())
taxa_summary_mot = biv_trait2 %>% group_by(motility) %>% summarise(n = n())
taxa_summary_feed = biv_trait2 %>% group_by(feeding) %>% summarise(n = n())


#################### CALCULATING FD METRICS ###################################

#Note: Have to set working directory to local computer for this chunk

#CALCULATE FD METRICS
biv_FD = dbFD(biv_trait2, biv_taxa2_alpha, corr = c("sqrt"), m = "max")
summary(biv_FD)
biv_FD_dataframe = as.data.frame(biv_FD)

#CULL TO MAKE SAME LENGTH AS RESULTS DATA FRAME
row.names(biv_results) = biv_results$Station
keepsites = rownames(biv_FD_dataframe) [(rownames(biv_FD_dataframe) %in% row.names(biv_results))] 
biv_FD_dataframe2 = biv_FD_dataframe [keepsites, ]

#ADD RESULTS TO MAIN DATA FRAME
biv_results$FRic = biv_FD_dataframe2$FRic #functional richness
biv_results$FEve = biv_FD_dataframe2$FEve #functional eveness
biv_results$FDiv = biv_FD_dataframe2$FDiv #functional divergence
biv_results$FDis = biv_FD_dataframe2$FDis #functional dispersion
biv_results$RaoQ = biv_FD_dataframe2$RaoQ #Rao's quadratic entropy (Q)


############## CALCULATE SPP RICHNESS AND FINISH FINAL DATA FRAME ###########################

#CALCULATE SPECIES RICHNESS AND ES(50) AND ES(20)
biv_results$Richness = specnumber(biv_taxa_binary.regional3[2:224,])

biv_results$rare20 = rarefy(biv_taxa_binary.regional3[2:224,], 20)
biv_results$rare50 = rarefy(biv_taxa_binary.regional3[2:224,], 50)

#

#CREATE FINAL DATA FRAME WITH ALL RESULTS
names(biv_results)[1] = "Station2"
biv_final = left_join(biv_results, biv_enviro, by="Station2")
dim(biv_final)

#

#ADD LOG(FLUX) AND LOG(FLUX^2)
biv_final$logFlux = log10(biv_final$Flux)
biv_final$logFlux2 = log10(biv_final$Flux^2)

#ADD LOG(TEMP)
biv_final$logTemp = log10(biv_final$Temp)


###################### GAMMS: ANALYSES FOR TABLE 1 #######################################

##Diversity

#GAM1 : Species Richness
gam_1 = mgcv::gam(Richness ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_1)


#GAM9 : ES(20)
gam_9 = mgcv::gam(rare20 ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_9)


##Expansion Metrics

#GAM2 : UTC
gam_2 = mgcv::gam(utc ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_2)


#GAM3 : FRic
gam_3 = mgcv::gam(FRic ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_3)


##Packing Metrics

#GAM4 : meanMVO (Functional Overlap)
gam_4 = mgcv::gam(meanmvo ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_4)


#GAM5 : FDis
gam_5 = mgcv::gam(FDis ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_5)


#GAM6 : RaoQ
gam_6 = mgcv::gam(RaoQ ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_6)

##Evenness of Packing

#GAM7 : FEve
gam_7 = mgcv::gam(FEve ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_7)


#GAM8 : FDiv
gam_8 = mgcv::gam(FDiv ~ s(logFlux) + s(logTemp) + Basin,
                  data = biv_final, 
                  method = "REML", select = TRUE)
summary(gam_8)

############### FIGURES 1 & 2 ##############################

#First have to convert each gam to a new object for use with mgcViz package

viz_1 = getViz(gam_1)
viz_2 = getViz(gam_2)
viz_3 = getViz(gam_3)
viz_4 = getViz(gam_4)
viz_5 = getViz(gam_5)
viz_6 = getViz(gam_6)
viz_7 = getViz(gam_7)
viz_8 = getViz(gam_8)
viz_9 = getViz(gam_9)


#Species Richness 
p1 = plot( sm(viz_1, 1) ) +
  l_fitLine(colour = "steelblue1", size = 1.5) +
  l_ciLine(mul = 5, colour = "steelblue4", linetype = 2, size = 1) + 
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Species Richness",
       tag = "(a)") +
  theme_craig()

#ES(20)
p9 = plot( sm(viz_9, 1) )  +
  l_fitLine(colour = "steelblue1", size = 1.5) +
  l_ciLine(mul = 5, colour = "steelblue4", linetype = 2, size = 1) +
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "ES(20)",
       tag = "(b)") +
  theme_craig()

#Unique Trait Combinations
p2 = plot( sm(viz_2, 1) ) +
  l_fitLine(colour = "steelblue1", size = 1.5) +
  l_ciLine(mul = 5, colour = "steelblue4", linetype = 2, size = 1) + 
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Unique Trait Combinations",
       tag = "(b)") +
  theme_craig()

#Functional Richness
p3 = plot( sm(viz_3, 1) ) +
  l_fitLine(colour = "steelblue1", size = 1.5) +
  l_ciLine(mul = 5, colour = "steelblue4", linetype = 2, size = 1) +
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Functional Richness",
       tag = "(a)") +
  theme_craig()

#Functional Overlap
p4 = plot( sm(viz_4, 1) ) +
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Functional Overlap",
       tag = "(c)") +
  theme_craig()

#Functional Dispersion
p5 = plot( sm(viz_5, 1) ) +
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Functional Dispersion",
       tag = "(d)") +
  theme_craig()

#Rao's Q
p6 = plot( sm(viz_6, 1) ) +
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Rao's Quadratic Entropy",
       tag = "(d)") +
  theme_craig()

#Functional Evenness
p7 = plot( sm(viz_7, 1) ) +
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Functional Evenness",
       tag = "(e)") +
  theme_craig()

#Functional Divergence
p8 = plot( sm(viz_8, 1) ) +
  l_fitLine(colour = "steelblue1", size = 1.5) +
  l_ciLine(mul = 5, colour = "steelblue4", linetype = 2, size = 1) + 
  l_points(shape = 16, size = 2, alpha = 0.5) + 
  labs(x = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}),
       y= "Partial Residuals",
       subtitle = "Functional Divergence",
       tag = "(f)") +
  theme_craig()






