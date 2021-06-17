#Energetic Functional Space Expansion
#Manuscript: "Functional space expansion driven by transitions between advantageous traits across a deep-sea energetic gradient"
#Authors: S. River D. Bryant & Craig R. McClain
##Null Model
##Code Author: River Bryant

##################### Set-Up #################################
#Load Packages
library(dplyr)         #for data manipulation
library(vegan)         #for diversity calculations
library(multirich)     #for functional diversity metrics
library(FD)            #for functional diversity metrics
library(ggplot2)       #for plotting
library(gridExtra)     #for plotting
library(reshape)       #for data manipulation

#Create ggplot theme
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

#Load Data
taxa = read.csv(file = "allen_bivalve_species.traits.csv", header = T)

#Delete taxa without trait info
taxa = na.omit(taxa)




################# Null Model Loop  ###############################

#NOTE: COMPUTATIONALLY INTENSIVE***************

#initialize a df for the results to feed into
results = data.frame(n=integer(),
                    utc=integer())
# 
#nested for loop - basic logic is it runs i times (1000) for each value of k ("spp richness")
 for(k in c(3:63)){
 for(i in c(1:1000)){
   
   n_result = sample_n(taxa, k, replace = FALSE) %>%    #sample k number of random species from the data without replacement
   distinct(ecospace) %>%                             #take only the unique ecospaces
     summarize(n = k, utc = n())                        #make a df with k and utc is the length of the df comprising unique ecospaces
   
   results = rbind(results, n_result)                   #add the results to the results df
   
 }
 }
# 
#save it off 
write.csv(results, "null_model_results.csv", row.names = FALSE)



############################## Empirical Comparison ##########################

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

#

#CREATE MAIN TRAIT DATA HOLDER
biv_trait_orig = read.csv("allen_bivalve_species.traits.csv", header=TRUE)

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


#DROP SPECIES THAT DO NOT OCCUR AT ANY STATION
biv_taxa2_alpha$Arc48 = NULL #drop Arc48 because does not occur at any station
biv_trait2 = biv_trait2[-c(5),] #drop Arc48 from the trait data too

#CALCULATE SPECIES RICHNESS AND ES(50) AND ES(20)
biv_results$Richness = specnumber(biv_taxa_binary.regional3[2:224,])

biv_results$rare20 = rarefy(biv_taxa_binary.regional3[2:224,], 20)
biv_results$rare50 = rarefy(biv_taxa_binary.regional3[2:224,], 50)


#CREATE FINAL DATA FRAME WITH ALL RESULTS
names(biv_results)[1] = "Station2"
biv_final = left_join(biv_results, biv_enviro, by="Station2")
dim(biv_final)

#ADD LOG(FLUX) AND LOG(FLUX^2)
biv_final$logFlux = log10(biv_final$Flux)
biv_final$logFlux2 = log10(biv_final$Flux^2)

#ADD LOG(TEMP)
biv_final$logTemp = log10(biv_final$Temp)



########################## SUPPLEMENTAL FIGURE 1 ##############################
#Read in null model data for the plot
plot_data = read.csv("null_model_results.csv", header = TRUE)

ggplot()+
  geom_smooth(data = plot_data, aes(x = n, y = utc), method = "gam", color = "blue", linetype = "dashed")+
  geom_smooth(data = biv_final, aes(x = Richness, y = utc), method = "gam", color = "red")+
  labs(x = "Species Richness",
       y = "UTC")+
  theme_craig()+
  theme(legend.position = "right")


