#Energetic Functional Space Expansion
#Manuscript: "Functional space expansion driven by transitions between advantageous traits across a deep-sea energetic gradient"
#Authors: S. River D. Bryant & Craig R. McClain
##Functional Types Across Flux: Figures & Analyses
##Code Author: River Dixon Bryant

#NOTE: "Ecospace" is synonymous with "Functional Type"


##########SETUP##############

##Load Packages
require(ggplot2)      #For plotting
require(reshape2)     #For fixing trait data
library(stringr)      #For data wrangling  
library(dplyr)        #For data wrangling
library(tidyr)        #For data wrangling
library(plotly)       #For 3D Plots
library(gridExtra)    #For Figure
library(grid)         #For figure
library(lme4)         #For glm with random effect
library(data.table)   #For casting with multiple "value" variables
library(forcats)      #For re-ordering factors
library(ggExtra)      #For plots
library(sjstats)      #For ANOVA effect sizes
library(ggridges)     #For plots

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

###################DATA MANIPULATION#############
#Getting Data into form for ecospace ~ flux analysis

#Create triplet from matrix
ecoStation.triplet = spStation %>%
  gather(Species_ID, Abundance, -Station)

#add ecospace and binary presence/absence to triplet
ecoStation.triplet$ecospace = traits$ecospace[match(ecoStation.triplet$Species_ID,traits$Species_ID)]
ecoStation.triplet$PA = as.numeric(ifelse(ecoStation.triplet$Abundance > 0, "1", "0"))

#Calculate abundance and richness per ecospace
ecodata1 = ecoStation.triplet %>% 
  group_by(Station, ecospace) %>%
  summarize(
    I = sum(Abundance),
    R = sum(PA)) 

#Add Presence/absence binary back into dataframe
ecodata1$PA = as.numeric(ifelse(ecodata1$I > 0, "1", "0"))
ecodata1 = as.data.frame(ecodata1)

#cast with multiple "value" variables
ecodata2 = dcast(setDT(ecodata1), Station~ecospace, value.var=c('I', 'R', 'PA'))

#delete data for incomplete ecospace "?61" and NAs
ecodata = ecodata2 %>% 
  select(-c("I_61", "R_61", "PA_61"))
ecodata = ecodata %>% mutate (I_NA = NULL,
                              R_NA = NULL,
                              PA_NA = NULL)

#add environmental data back in
ecodata$Station = as.character(ecodata$Station)
stations$Station2 = as.character(stations$Station2)

ecodata_env = left_join(ecodata, stations, b=c("Station" = "Station2"))

ecodata_env = ecodata_env %>%
  select(-c("Individuals", "Lat", "Long", "Depth",
            "Salinity", "DissOx", "Station.y", "Station3"))

ecodata_env$log10Flux = log10(ecodata_env$Flux)


#PA: Presence/Absence
PAdata1 = ecodata %>%
  select(c(1, 42:61))

PAdata1 = as.data.frame(PAdata1)
PAdata1$Station = as.character(PAdata1$Station)

PAdata2 = left_join(PAdata1, stations, b=c("Station" = "Station2"))

PAdata = PAdata2 %>%
  select(-c("Individuals", "Lat", "Long", "Depth",
            "Salinity", "DissOx", "Station.y", "Station3"))

PAdata$log10Flux = log10(PAdata$Flux)

PAdata_long = PAdata %>%
  gather(Ecospace, PA, PA_321:PA_651)
PAdata_long$eco1 = str_sub(PAdata_long$Ecospace, start= -3)

#Abundance (I)
Idata1 = ecodata %>%
  select(c(1:21))

Idata1 = as.data.frame(Idata1)
Idata1$Station = as.character(Idata1$Station)

Idata2 = left_join(Idata1, stations, b=c("Station" = "Station2"))

Idata = Idata2 %>%
  select(-c("Individuals", "Lat", "Long", "Depth",
            "Salinity", "DissOx", "Station.y", "Station3"))

Idata$log10Flux = log10(Idata$Flux)

Idata_long = Idata %>%
  gather(Ecospace, I, I_321:I_651)
Idata_long$eco1 = str_sub(Idata_long$Ecospace, start= -2)

##Reorder so that it plots with descending mean
Idata_long$eco1 = as.factor(Idata_long$eco1)
Idata_long$eco1 = with(Idata_long, fct_reorder(eco1, log10Flux, mean, .desc = T))

Idata_long$eco2 = as.factor(ifelse( Idata_long$eco1 == "321", "1",
                                     ifelse(Idata_long$eco1 == "331", "2",
                                            ifelse(Idata_long$eco1 == "335", "3",
                                                   ifelse(Idata_long$eco1 == "341", "4", 
                                                          ifelse(Idata_long$eco1 == "345", "5", 
                                                                 ifelse(Idata_long$eco1 == "351", "6",
                                                                        ifelse(Idata_long$eco1 == "361", "7",
                                                                               ifelse(Idata_long$eco1 == "521", "8",
                                                                                      ifelse(Idata_long$eco1 == "522", "9",
                                                                                             ifelse(Idata_long$eco1 == "525", "10",
                                                                                                    ifelse(Idata_long$eco1 == "526", "11",      
                                                                                                           ifelse(Idata_long$eco1 == "531", "12",
                                                                                                                  ifelse(Idata_long$eco1 == "532", "13",
                                                                                                                         ifelse(Idata_long$eco1 == "533", "14",
                                                                                                                                ifelse(Idata_long$eco1 == "535", "15",
                                                                                                                                       ifelse(Idata_long$eco1 == "536", "16",      
                                                                                                                                              ifelse(Idata_long$eco1 == "551", "17",
                                                                                                                                                     ifelse(Idata_long$eco1 == "561", "18",
                                                                                                                                                            ifelse(Idata_long$eco1 == "631", "19",
                                                                                                                                                                   ifelse(Idata_long$eco1 == "651", "20", NA
                                                                                                                                                                   )))))))))))))))))))))

##Reorder so that it plots with descending mean
Idata_long$eco2 = with(Idata_long, fct_reorder(eco2, log10Flux, mean, .desc = T))









######### FIGURE 3 ############
##Reorder so that it plots with descending mean
PAdata_long$eco1 = as.factor(PAdata_long$eco1)
PAdata_long$eco1 = with(PAdata_long, fct_reorder(eco1, log10Flux, mean, .desc = T))

## Set-up ecospaces as 1-18 to plot correctly
PAdata_long$eco2 = as.factor(ifelse( PAdata_long$eco1 == "321", "1",
                                     ifelse(PAdata_long$eco1 == "331", "2",
                                            ifelse(PAdata_long$eco1 == "335", "3",
                                                   ifelse(PAdata_long$eco1 == "341", "4", 
                                                          ifelse(PAdata_long$eco1 == "345", "5", 
                                                                 ifelse(PAdata_long$eco1 == "351", "6",
                                                                        ifelse(PAdata_long$eco1 == "361", "7",
                                                                               ifelse(PAdata_long$eco1 == "521", "8",
                                                                                      ifelse(PAdata_long$eco1 == "522", "9",
                                                                                             ifelse(PAdata_long$eco1 == "525", "10",
                                                                                                    ifelse(PAdata_long$eco1 == "526", "11",      
                                                                                                           ifelse(PAdata_long$eco1 == "531", "12",
                                                                                                                  ifelse(PAdata_long$eco1 == "532", "13",
                                                                                                                         ifelse(PAdata_long$eco1 == "533", "14",
                                                                                                                                ifelse(PAdata_long$eco1 == "535", "15",
                                                                                                                                       ifelse(PAdata_long$eco1 == "536", "16",      
                                                                                                                                              ifelse(PAdata_long$eco1 == "551", "17",
                                                                                                                                                     ifelse(PAdata_long$eco1 == "561", "18",
                                                                                                                                                            ifelse(PAdata_long$eco1 == "631", "19",
                                                                                                                                                                   ifelse(PAdata_long$eco1 == "651", "20", NA
                                                                                                                                                                   )))))))))))))))))))))

##Reorder so that it plots with descending mean
PAdata_long$eco2 = with(PAdata_long, fct_reorder(eco2, log10Flux, mean, .desc = T))


##Filter out only presence values
PAdata_long_P = PAdata_long %>%
  filter(PA == "1")

PAdata_long_P$eco2 = with(PAdata_long_P, fct_reorder(eco2, log10Flux, mean, .desc = T))
PAdata_long_P$eco1 = with(PAdata_long_P, fct_reorder(eco1, log10Flux, mean, .desc = T))


##Boxplot
#breaks are along quartiles
p1 = ggplot(PAdata_long_P, aes(x =eco1, y = log10Flux)) +
  geom_boxplot(alpha = 0.5, fill = "black") +
  geom_hline(yintercept = c(0.56,0.74, 1.01, 1.94), linetype="dashed", 
             color = "steelblue4")+
  annotate(geom = "text", x=as.factor(536), y=0.35, label="Lower 25% Flux", color = "steelblue4", fontface = "bold")+
  annotate(geom = "text", x=as.factor(536), y=0.8, label="Middle 50% Flux", color = "steelblue4", fontface = "bold")+
  annotate(geom = "text", x=as.factor(536), y=1.5, label="Upper 25% Flux", color = "steelblue4", fontface = "bold")+
  coord_flip() +
  labs(x = "Ecospace",
       y = expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2})) + 
  theme_craig()

#################### STATS: ANOVA & EMMEANS #######################
PAdata_long_P$scaleFlux = scale(PAdata_long_P$Flux, center = T, scale = T)

##Testing Assumptions
#Q-Q plot
qqnorm(PAdata_long_P$log10Flux)
qqline(PAdata_long_P$log10Flux)

qqnorm(PAdata_long_P$Flux)
qqline(PAdata_long_P$Flux)

qqnorm(PAdata_long_P$scaleFlux)
qqline(PAdata_long_P$scaleFlux)

#   log-transformed definitely best

#Shapiro-Wilke Test
shapiro.test(PAdata_long_P$log10Flux)
shapiro.test(PAdata_long_P$Flux)
shapiro.test(PAdata_long_P$scaleFlux)

#visualize with histogram
hist(PAdata_long_P$log10Flux)
hist(PAdata_long_P$Flux)
hist(PAdata_long_P$scaleFlux)

#Bartlett Test
#requires 2+ observations per group so must 
#  drop ecospaces "536", "631", "351"
PAdata_bartlett = PAdata_long_P %>% 
                  filter(eco1 != "536", 
                         eco1 != "631", 
                         eco1 != "351")
bartlett.test(log10Flux ~ eco1, data = PAdata_bartlett) 
#    variance is different among groups (p=0.001)



##ANOVA for Functional Types
ecospace.aov = aov(log10Flux ~ eco1, data = PAdata_long_P)
summary(ecospace.aov)

#posthoc
TukeyHSD(ecospace.aov)


##Estimated Marginal Means (emmeans) analyses
library(emmeans)

spp_per_eco = traits %>% group_by(ecospace) %>% summarize(eco_n = n())
spp_per_eco1 = spp_per_eco[2:21,]
spp_per_eco2 = transform(spp_per_eco1, tiering = substr(ecospace, 1, 1), 
                          motility = substr(ecospace, 2, 2),
                          feeding = substr(ecospace, 3, 3))

richness_trait.aov = aov(eco_n ~ tiering + motility + feeding, data = spp_per_eco2)
summary(richness_trait.aov)
TukeyHSD(richness_trait.aov)

################### Figure 4 ####################################
emm1 = emmeans(ecospace.aov, ~ eco1)

pw <- pwpp(emm1)+
  geom_ribbon(aes(xmin=0.05, xmax=1), alpha=0.01)+
  annotate("text", x = 1.1, y = 1:20, label = c("2", "2", "26", "87", "6", "2", "93", "105", "32", "8", "25", "23", "9", "1", "7", "29", "8", "2", "1", "1"), alpha = 0.8, size = 3)+
  theme_craig()+
  xlab("Tukey-Adjusted P-Value")+
  ylab("Ecospace")+
  ggtitle("Log-Transformed Carbon Flux")+
  theme(plot.title = element_text(size = 10, face = "bold"))

