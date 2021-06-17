#Energetic Functional Space Expansion
#Manuscript: "Functional space expansion driven by transitions between advantageous traits across a deep-sea energetic gradient"
#Authors: S. River D. Bryant & Craig R. McClain
##Individual Functional Traits Across Flux: Plots & Analyses
##Code Author: River Bryant


################################ Set-Up ########################################
##Load Packages
require(ggplot2)      #For plotting
require(reshape2)     #For fixing trait data
library(stringr)      #For data wrangling  
library(dplyr)        #For data wrangling
library(tidyr)        #For data wrangling
library(gridExtra)    #For the Figure
library(grid)         #For the figure
library(lme4)         #For glm with random effect
library(data.table)   #For casting with multiple "value" variables
library(forcats)      #For re-ordering factors
library(ggExtra)      #For pretty plots
library (mgcv)        #For GAMMs
library(broom)

##ggplot stuff
theme_craig <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      # change stuff here
      axis.line = element_line(colour = "darkgrey"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      legend.position="none",
      plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
}

##Read in data
traits    = read.csv(file = "allen_bivalve_species.traits.csv", header = T)
spStation = read.csv(file = "allen_bivalve_sp.station.csv", header = T)
stations  = read.csv(file = "allen_bivalve_station.csv", header = T)

#Create triplet from matrix
ecoStation.triplet = spStation %>%
  gather(Species_ID, Abundance, -Station)

#Calculate log10Flux
stations$log10Flux = log10(stations$Flux)

#Add log10Flux, tiering, motility, feeding to triplet
ecoStation.triplet$tiering = traits$tiering[match(ecoStation.triplet$Species_ID,traits$Species_ID)]
ecoStation.triplet$motility = traits$motility[match(ecoStation.triplet$Species_ID,traits$Species_ID)]
ecoStation.triplet$feeding = traits$feeding[match(ecoStation.triplet$Species_ID,traits$Species_ID)]
ecoStation.triplet$log10Flux = stations$log10Flux[match(ecoStation.triplet$Station,stations$Station2)]


#######################Proportion of Species#########################
sprows = spStation[,1]
spStation2 = spStation[,2:529]
spStation_binary = as.data.frame(+(spStation2 >= 1))
row.names(spStation_binary) = sprows
spStation_binary$Station = sprows

binary_triplet = reshape2::melt(spStation_binary, id.vars=c("Station"), variable.name = "Species_ID", value.name = "PA") 

#Add log10Flux, tiering, motility, feeding to triplet
binary_triplet$tiering = traits$tiering[match(binary_triplet$Species_ID,traits$Species_ID)]
binary_triplet$motility = traits$motility[match(binary_triplet$Species_ID,traits$Species_ID)]
binary_triplet$feeding = traits$feeding[match(binary_triplet$Species_ID,traits$Species_ID)]
binary_triplet$log10Flux = stations$log10Flux[match(binary_triplet$Station,stations$Station2)]
binary_triplet$Temp = stations$Temp[match(binary_triplet$Station,stations$Station2)]
binary_triplet$Basin = stations$Basin[match(binary_triplet$Station,stations$Station2)]


prop_feed_data = binary_triplet %>% select("log10Flux", "PA", "feeding", "Basin", "Temp") %>%
  group_by(Basin, Temp, log10Flux, feeding) %>% 
  summarize(Richness = sum(PA))

prop_tier_data = binary_triplet %>% select("log10Flux", "PA", "tiering", "Basin", "Temp") %>%
  group_by(Basin, Temp, log10Flux, tiering) %>% 
  summarize(Richness = sum(PA))

prop_mot_data = binary_triplet %>% select("log10Flux", "PA", "motility", "Basin", "Temp") %>%
  group_by(Basin, Temp, log10Flux, motility) %>% 
  summarize(Richness = sum(PA))


######################FIGURE 5#####################
###Feeding

#Suspension Feeder
feed1 = prop_feed_data %>% filter(feeding == 1)
total_richness = prop_feed_data %>% group_by(log10Flux) %>% summarize (totalRichness = sum(Richness))
susp_figure2 = merge(x = feed1, y = total_richness, by.x = "log10Flux")

susp_figure2$rel_richness = susp_figure2$Richness / susp_figure2$totalRichness

p11 = ggplot(aes(x=log10Flux, y = rel_richness), data = susp_figure2)+
  geom_point(alpha = 0.35, fill = "black", size = 2)+
  geom_smooth(method = "gam",color = "steelblue1",  se=F)+
  geom_smooth(method = "loess",color = "steelblue4", linetype = "dashed",  se=F)+
  ylim(0,1)+
  theme_craig()+
  labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), 
       y="Proportion of Species",
       tag = "(a)")+
  ggtitle("Suspension Feeders")



#Deposit Feeder
feed2 = prop_feed_data %>% filter(feeding == 2)
dep_figure2 = merge(x = feed2, y = total_richness, by.x = "log10Flux")

dep_figure2$rel_richness = dep_figure2$Richness / dep_figure2$totalRichness

p12 = ggplot(aes(x=log10Flux, y = rel_richness), data = dep_figure2)+
  geom_point(alpha = 0.35, fill = "black", size = 2)+
  geom_smooth(method = "gam",color = "steelblue1",  se=F)+
  geom_smooth(method = "loess",color = "steelblue4", linetype = "dashed",  se=F)+
  ylim(0,1)+
  theme_craig()+
  labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), 
       y=" ",
       tag = "(b)")+
  ggtitle("Deposit Feeders")


#####Motility

#Freely, Slow
mot2 = prop_mot_data %>% filter(motility == 2)
total_richness_mot = prop_mot_data %>% group_by(log10Flux) %>% summarize (totalRichness = sum(Richness))
free_figure2 = merge(x = mot2, y = total_richness_mot, by.x = "log10Flux")

free_figure2$rel_richness = free_figure2$Richness / free_figure2$totalRichness

p19 = ggplot(aes(x=log10Flux, y = rel_richness), data = free_figure2)+
  geom_point(alpha = 0.35, fill = "black", size = 2)+
  geom_smooth(method = "gam",color = "steelblue1",  se=F)+
  geom_smooth(method = "loess",color = "steelblue4", linetype = "dashed",  se=F)+
  ylim(0,1)+
  theme_craig()+
  labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), 
       y=" ",
       tag = "(d)")+
  ggtitle("Freely, Slow")

#Faculative unattached
mot3 = prop_mot_data %>% filter(motility == 3)
facun_figure2 = merge(x = mot3, y = total_richness_mot, by.x = "log10Flux")

facun_figure2$rel_richness = facun_figure2$Richness / facun_figure2$totalRichness

p20 = ggplot(aes(x=log10Flux, y = rel_richness), data = facun_figure2)+
  geom_point(alpha = 0.35, fill = "black", size = 2)+
  geom_smooth(method = "gam",color = "steelblue1",  se=F)+
  geom_smooth(method = "loess",color = "steelblue4", linetype = "dashed",  se=F)+
  ylim(0,1)+
  theme_craig()+
  labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), 
       y="Proportion of Species",
       tag = "(c)")+
  ggtitle("Facultative, Unattached")

#####Tiering

#Surficial
tier3 = prop_tier_data %>% filter(tiering == 3)
total_richness_tier = prop_tier_data %>% group_by(log10Flux) %>% summarize (totalRichness = sum(Richness))
surf_figure2 = merge(x = tier3, y = total_richness_tier, by.x = "log10Flux")

surf_figure2$rel_richness = surf_figure2$Richness / surf_figure2$totalRichness

p24 = ggplot(aes(x=log10Flux, y = rel_richness), data = surf_figure2)+
  geom_point(alpha = 0.35, fill = "black", size = 2)+
  geom_smooth(method = "gam",color = "steelblue1",  se=F)+
  geom_smooth(method = "loess",color = "steelblue4", linetype = "dashed",  se=F)+
  ylim(0,1)+
  theme_craig()+
  labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), 
       y="Proportion of Species",
       tag = "(e)")+
  ggtitle("Surficial")

#Shallow Infaunal
tier5 = prop_tier_data %>% filter(tiering == 5)
shin_figure2 = merge(x = tier5, y = total_richness_tier, by.x = "log10Flux")

shin_figure2$rel_richness = shin_figure2$Richness / shin_figure2$totalRichness

p25 = ggplot(aes(x=log10Flux, y = rel_richness), data = shin_figure2)+
  geom_point(alpha = 0.35, fill = "black", size = 2)+
  ylim(0,1)+
  theme_craig()+
  labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), 
       y=" ",
       tag = "(f)")+
  ggtitle("Shallow Infaunal")

#Final Figure
grid.arrange(p11, p12, p20, p19, p24, p25, nrow = 3, respect = T)

###################Stats: GAMMS########################
susp_figure2$log10Temp = log10(susp_figure2$Temp)
dep_figure2$log10Temp = log10(dep_figure2$Temp)

#Suspension
gam_susp = mgcv::gam(rel_richness ~ s(log10Flux) +  Basin,
                     data = susp_figure2, 
                     method = "REML", select = TRUE)
summary(gam_susp)

#Deposit
gam_dep = mgcv::gam(rel_richness ~ s(log10Flux)  + Basin,
                    data = dep_figure2, 
                    method = "REML", select = TRUE)
summary(gam_dep)

#Facultative, unattached
gam_facun = mgcv::gam(rel_richness ~ s(log10Flux)  + Basin,
                    data = facun_figure2, 
                    method = "REML", select = TRUE)
summary(gam_facun)

#Freely, Slow
gam_freesl = mgcv::gam(rel_richness ~ s(log10Flux)  + Basin,
                      data = free_figure2, 
                      method = "REML", select = TRUE)
summary(gam_freesl)

#Surficial
gam_surf = mgcv::gam(rel_richness ~ s(log10Flux)  + Basin,
                       data = surf_figure2, 
                       method = "REML", select = TRUE)
summary(gam_surf)

#Shallow infaunal
gam_shin = mgcv::gam(rel_richness ~ s(log10Flux)  + Basin,
                       data = shin_figure2, 
                       method = "REML", select = TRUE)
summary(gam_shin)

