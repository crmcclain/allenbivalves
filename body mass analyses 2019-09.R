# Deep sea bivalve diveristy, energy, metabolism
# Aug 2019 version

###############################load packages###################################################
#' Load required packages
  library(tidyverse)
  library(vegan)
  library(lme4)
  library(MuMIn)
  library(ggplot2)
  library(RColorBrewer)
  library(quantreg)
  library(broom)
  library(lmerTest)
  library(emmeans)
  library(gridExtra)
  library(grid)
  library(ggpubr)

############################### load data###################################################
#' BIVALVES
#' 
##' Read in and tidy data
  setwd("~/Dropbox/DS Metabolic Diversity/Data/")
  b.speciesstation = read_csv("allen_bivalve.csv") %>% rename(station = X1)
  b.stations = read_csv("bivalve_stations.csv")
  b.species = read_csv("bivalve_species.csv")

############################### custom theme###################################################
  theme_craig <- function () { 
    theme_bw(base_size=12) %+replace% 
      theme(
        # change stuff here
        axis.line = element_line(colour = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0))}

############################### data manipulation###################################################
# Calculate station environmental and species trait metrics
    b.species <- b.species %>%
      mutate(
        log2BV = log2(Biovolume_mm3),
        log10Mass = (0.9575 * log10(Biovolume_mm3)) - 4.8939,
        Mass = 10^log10Mass,
        log10Mtot = log10(`B0 Order` * Mass^.75)
      )

    b.stations <- b.stations %>%
      mutate(
        log10Flux = log10(Flux),
        log10Temp = log10(Temp),
        standFlux = scale(log10Flux)[,1],
        standTemp = scale(log10Temp)[,1]
      )

#' Create a `stations` variable in `b.stations` to match station codes as listed in b.speciestation.
#' There are 5 anomalies that need dealing with separately
    b.stations <- b.stations %>%
      mutate(
        stat_code = str_replace_all(Station2, " ", "")
      ) %>%
      mutate(
        stat_code = case_when(
          stat_code == "st1984" ~ "st198",
          stat_code == "st3063" ~ "st306",
          stat_code == "stCh.10" ~ "stCh10",
          stat_code == "stCh.6" ~ "stCh6",
          stat_code == "stPolygasDS234" ~ "stPolygasDS23",
          TRUE ~ stat_code
        )
      )

#' Create tidy version of species x station
    b_spxst <- b.speciesstation %>%
      gather(species, abundance, -station) %>%
      filter(abundance > 0)

#' Note some species not present in main species list - are these important?
    sort(unique(b_spxst$species[!(b_spxst$species %in% b.species$Species_ID)]))

#' response from Craig: checked the absent data on the bivalves, I left those species out because
#' I had concerns about their IDs, i.e. I could not confirm a valid species name for them.
#' The one exception was Serp1, which was in the dataset as SERP1. So:
    b.species$Species_ID[b.species$Species_ID == "SERP1"] <- "Serp1"


#' Add in some relevant species- and site- specific parameters
    b_spxst <- left_join(b_spxst, dplyr::select(b.stations, Basin:DissOx, stat_code), by = c("station" = "stat_code"))
    b_spxst <- left_join(b_spxst, dplyr::select(b.species, Species_ID, Biovolume_mm3, Order, `B0 Order`, log10Mass, Mass, log10Mtot),
                         by = c("species" = "Species_ID")) %>% rename(B0 = `B0 Order`)

#' For analysis, remove rows with no mass info:
    b_spxst <- b_spxst %>% filter(!is.na(Mass))

## Body mass groups based on quantiles
#' Decide on quantiles to use - here I use 4 groups (i.e. 25% quantiles) but this can be set to what you want
    n_m_qs <- 4

#' Use this to create a new categorical variable,
#' `M_quant`, classifying each species according to its `Mass`:
    b.species <- b.species %>%
      mutate(M_quant = cut(log10(Mass),
                           breaks = quantile(log10(Mass), 0:n_m_qs/n_m_qs, na.rm = TRUE),
                           labels = round(0:(n_m_qs - 1)/n_m_qs, 2), include.lowest = TRUE)
      )
    
    #alternative cut
    sizerange <- range(log10(b.species$Mass))
    size_by <- abs(sizerange[1]-sizerange[2])/4
    size_cuts <- seq(sizerange[1], sizerange[2], by=size_by)
    size_cuts
    
    b.species <- b.species %>%
      mutate(M_bins = cut(log10(Mass),
                           breaks = size_cuts,
                           labels = round(0:(n_m_qs - 1)/n_m_qs, 2), include.lowest = TRUE)
      )

#' Now, join the quantile info to the species x station table:
    b_spxst <- left_join(b_spxst, select(b.species, Species_ID, M_quant), by = c("species" = "Species_ID"))
    
    b_spxst_bin <- left_join(b_spxst, select(b.species, Species_ID, M_bins), by = c("species" = "Species_ID"))
    
#' In order to get diversity per quantile per station, we need a function that can take a tidy data frame
#' and apply the `vegan` diveristy functions (this also returns number of individuals and species):
    tidy_div <- function(df, site = "station", sp = "species", ab = "abundance", d_index = "shannon"){
      x <- spread(select(df, site, sp, ab), sp, ab) %>%
        replace(., is.na(.), 0) %>%
        as.data.frame()
      
      rownames(x) <- select(x, site) %>% pull()
      x <- x[, -1]
      d <- vegan::diversity(x, index = d_index)
      d <- tibble(site = names(d), d = as.vector(d))
      d <- d %>% mutate(n_ind = rowSums(x), n_sp = rowSums(x > 0))
      
      d
    }

#' To run this over all M_quants to get diversity per metabolic quantile per station, and tidy up the output:
    b_spxst_q <- b_spxst %>%
      group_by(M_quant) %>% 
      do(., H = tidy_div(df = .)) %>% 
      unnest() %>% 
      filter(!is.na(M_quant))
    
    b_spxst_q_bin <- b_spxst_bin %>%
      group_by(M_bins) %>% 
      do(., H = tidy_div(df = .)) %>% 
      unnest() %>% 
      filter(!is.na(M_bins))

#' Add the station data to this, and rename the generic diversity variable 'd' to H
#' (which is what we calculated in this instance):
    b_spxst_q <- left_join(b_spxst_q, b.stations, by = c("site" = "stat_code")) %>% 
      rename(H = d)

    b_spxst_q_bin <- left_join(b_spxst_q_bin, b.stations, by = c("site" = "stat_code")) %>% 
      rename(H = d)

############################### correlation analysis###################################################
## Correlation between mass and metabolism
    b.species2 <- b.species %>% filter(!is.na(log10Mtot))
    cor(b.species2$log10Mass, b.species2$log10Mtot) 
    
    ggplot(data=b.species2, aes(log10Mass, log10Mtot))+
      geom_point(cex=3, alpha=0.1)+
      labs(x=expression(log["10"]*Mass*" "*(g)), y=expression(log["10"]*Met["tot"]*" "*(mW)))+
      theme_craig()
    
############################### plot of body sizes###################################################
ggplot(data=b_spxst, aes(log10Mass, group=M_quant, fill=M_quant))+
      geom_histogram(bins=200)+
      labs(x=expression(log["10"]*Mass*" "*(g)), y="Count", fill="Mass Quantiles")+
      theme_craig()+
      theme(legend.position = c(0.8, 0.8))+
      scale_fill_brewer(palette = "YlOrRd")
  
ggsave("figure1.pdf")

ggplot(data=b_spxst, aes(log10Mass, group=M_quant, fill=M_quant))+
  geom_histogram(bins=200)+
  labs(x=expression(log["10"]*Mass*" "*(g)), y="Count", fill="Mass Quantiles")+
  theme_craig()+
  theme(legend.position = c(0.8, 0.8))+
  scale_fill_brewer(palette = "YlOrRd")

############################### Species v Individuals###################################################
#' Look at number of species v number of individuals for each Mass quantile:

b_spxst_q$M_quant <- as.factor(b_spxst_q$M_quant)
b_spxst_q$logI <- log10(b_spxst_q$n_ind)
b_spxst_q$logS <- log10(b_spxst_q$n_sp)

#plot
    (sp_v_ind_lm <- ggplot(b_spxst_q, aes(x = logI, y = logS)) +
        geom_smooth(aes(colour = M_quant, fill = M_quant), method = "lm", )+
        labs(x=expression(log["10"]*Individuals), y=expression(log["10"]*Species), color="Mass Quantiles")+
        scale_color_brewer(palette = "YlOrRd")+
        scale_fill_brewer(palette = "YlOrRd")+
        guides(fill = FALSE, color=guide_legend(override.aes=list(fill=NA)))+
        theme_craig())


  ggsave("figure2.pdf")
#analysis
    
    sp_v_ind_analyses2 <- lm(data=b_spxst_q, logS~logI+M_quant+logI*M_quant)
    anova(sp_v_ind_analyses2)
    summary(sp_v_ind_analyses2)
    pvalues <- summary(sp_v_ind_analyses2)$coefficients[,4]  
 
    #posthoc tests
    #slopes
    emtrends(sp_v_ind_analyses2, pairwise ~ M_quant, var = "logI",adjust = "holm" )
    #intercepts
    emmeans(sp_v_ind_analyses2, pairwise ~ M_quant | logI, adjust = "holm" )
    
    #independent regressions for each M_quant
    sp_v_ind_analyses <- lmList(data=b_spxst_q, logS~logI | M_quant)
    summary(sp_v_ind_analyses)
    sapply(sp_v_ind_analyses,function(x) summary(x)$r.squared)
    
    
    
############################### Basic Plotting Function################################################### 
#' Now dig in to the relationships between flux and diversity etc. First, produce some plots
#' These ignore basin for now.
#' This is a generic function that allows you to specify a diversity variable for the y axis
#' defaults to (H') against an environmental variable on the x axis (defaults to flux),
#' separately for each Mass group, fitting a second order polynomial overall as well as for the
#' 0.1 and 0.9 quantiles in each case. In other words, for each mass group we are examining whether
#' there is overall a relationship between flux and diversity, and whether this is driven more
#' by changes in maximum or minimum diversity with env.

div_v_env_plot <- function(dat = b_spxst_q, div_var = "H", env_var = "log10Flux",
                           w_var = NA, facet_var = "M_quant", d_lims = c(0, NA)){

  #dat = b_spxst_q; div_var = "H"; env_var = "log10Flux"
  #w_var = NA; facet_var = "M_quant"; d_lims = c(0, NA)
  #rm(dat, div_var, env_var, w_var, facet_var, d_lims)
  
  div_var <- sym(div_var)
  env_var <- sym(env_var)
  
  if(is.na(w_var)){
    dplot <-
      ggplot(dat, aes(x = !!env_var, y = !!div_var)) +
      geom_point(alpha = .7, cex=3, aes(color=M_quant)) +
      geom_smooth(method = "lm", formula = y ~ poly(x, 2), colour = "black") +
      geom_quantile(quantiles = c(0.1, 0.9), formula = y ~ poly(x, 2), colour = "grey50") +
      ylim(d_lims) +
      theme_craig()+
      labs(x=expression(log["10"]*Flux*" "*mgCDay^{-1}*M^{-2}), y="H'", color="Mass Quantiles")+
      scale_color_brewer(palette = "YlOrRd")+
      facet_wrap(~ dat[[facet_var]])
  } else {
    w_var <- sym(w_var)
    dplot <-
      ggplot(dat, aes(x = !!env_var, y = !!div_var, weight = !!w_var)) +
      geom_point(alpha = .7, cex=3, aes(color=M_quant)) +
      geom_smooth(method = "glm", method.args = list(family = binomial), colour = "black") +
      geom_quantile(quantiles = c(0.1, 0.9), formula = y ~ poly(x, 2), colour = "grey50") +
      ylim(d_lims) +
      theme_craig()+
      labs(x=expression(log["10"]*Flux*mgCDay^{-1}*M^{-2}), y="H'", color="Mass Quantiles")+
      scale_color_brewer(palette = "YlOrRd")+
      facet_wrap(~ dat[[facet_var]])
    }
  
  dplot
  
  }


############################### H' Plots and Analyses ###############################
#' so to plot H against flux (removing mass-group/site combos with only one species
  #' which have H = 0 by definition):
  div_v_env_plot(dat = filter(b_spxst_q, n_sp > 1))
  ggsave("figure3.pdf")

#' From this - H' has a hump-shaped relationship with flux (peaking at intermediate fluxes) for species
  #' of low-to-medium size, this changes to a largely increasing trend in the largest species.
  #' These trends are driven by changes in maximum diversity (top orange line, 90% quantile regression),
  #' as very low diversity is observed in all metabolic groups at all flux levels
  #' (bottom orange line, 10% quantile regression).

#' H against temperature:
  div_v_env_plot(env_var = "log10Temp", dat = filter(b_spxst_q, n_sp > 1))
#' H against depth:
  div_v_env_plot(env_var = "Depth", dat = filter(b_spxst_q, n_sp > 1))

#' Fit some models:
  #' This is the 'full' model, which includes flux (2nd order poly), temperature (2nd order poly), mass group
  #' (as a 4-level factor), Basin, and all interactions, removing site-M_quant combos with only one species:
  b_spxst_q_sub <- subset(b_spxst_q, n_sp > 1)
  
  ##full model
      b_h_flux_temp_mod_global <- lm(H ~ (poly(log10Flux, 2, raw = TRUE) + log10Temp +
                                M_quant + Basin)^2,
                              data = b_spxst_q_sub,
                              na.action = "na.fail")
      b_h_flux_temp_mod_global.aov <- aov(b_h_flux_temp_mod_global)
      
      b_h_flux_temp_mod_global.aov.tidy <- tidy(b_h_flux_temp_mod_global.aov)
      write.csv(b_h_flux_temp_mod_global.aov.tidy, file = "b_h_flux_temp_mod_global_aov_tidy.csv")
      
      b_h_flux_temp_mod_global.aov.glance <- glance(b_h_flux_temp_mod_global.aov)
      write.csv(b_h_flux_temp_mod_global.aov.glance, file = "b_h_flux_temp_mod_global_aov_glance.csv")
      

  #drop temp because not signficant  
      b_h_flux_mod_global <- lm(H ~ (poly(log10Flux, 2, raw = TRUE)+
                                       M_quant + Basin)^2,
                                     data = b_spxst_q_sub,
                                     na.action = "na.fail")
      
      b_h_flux_mod_global.aov <- aov(b_h_flux_mod_global)
    
      b_h_flux_mod_global.aov.tidy <- tidy(b_h_flux_mod_global.aov)
      write.csv(b_h_flux_mod_global.aov.tidy, file = "b_h_flux_mod_global_aov_tidy.csv")
      
      b_h_flux_mod_global.aov.glance <- glance(b_h_flux_mod_global.aov)
      write.csv(b_h_flux_mod_global.aov.glance, file = "b_h_flux_mod_global_aov_glance.csv")
      
  #posthoc tests
      #slopes
      emtrends(b_h_flux_mod_global, pairwise ~ M_quant, var = "poly(log10Flux, 2, raw = TRUE)",adjust = "holm" )
      #intercepts
      emmeans(b_h_flux_mod_global, pairwise ~ M_quant | poly(log10Flux, 2, raw = TRUE), adjust = "holm" )
      
  #' Although Basin is not in the top vars, probably a good idea to fit it as a random effect:
  b_h_flux_temp_mod_lme <- lmer(H ~ poly(log10Flux, 2, raw = TRUE) +
                                  M_quant +
                                  M_quant:poly(log10Flux, 2, raw = TRUE) +
                                  (1|Basin),
                                 data = filter(b_spxst_q, n_sp > 1),
                                 na.action = "na.fail")
 
    summary(b_h_flux_temp_mod_lme)
    anova(b_h_flux_temp_mod_lme)
    glance(b_h_flux_temp_mod_lme)
    
    
  #Cook's D
    
    cooksd <- cooks.distance(b_h_flux_mod_global)
    # Plot the Cook's Distance using the traditional 4/n criterion
    sample_size <- nrow(bb_spxst_q_sub)
    plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
    abline(h = 4/sample_size, col="red")  # add cutoff line
    text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
    
    # Removing Outliers
    # influential row numbers
    influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
    influential <- na.omit(influential)
    
    b_spxst_q_cook <- b_spxst_q_sub[-influential, ]
  
    #new plots
    div_v_env_plot(dat = filter(b_spxst_q_sub, n_sp > 1)) + ggtitle("Original")
    div_v_env_plot(dat = filter(b_spxst_q_cook, n_sp > 1)) + ggtitle("Filtered")
    
    b_h_flux_mod_global2 <- lm(H ~ (poly(log10Flux, 2, raw = TRUE)+
                                     M_quant + Basin)^2,
                              data = b_spxst_q_cook,
                              na.action = "na.fail")
    
    b_h_flux_mod_global.aov2 <- aov(b_h_flux_mod_global2)
    summary(b_h_flux_mod_global2)
    summary(b_h_flux_mod_global.aov2)
    
    
    #slopes
    emtrends(b_h_flux_mod_global2, pairwise ~ M_quant, var = "poly(log10Flux, 2, raw = TRUE)",adjust = "holm" )
    
    #Basin Differences?
    emtrends(b_h_flux_mod_global2, pairwise ~ Basin, var = "poly(log10Flux, 2, raw = TRUE)",adjust = "holm" )
    emmeans(b_h_flux_mod_global, pairwise ~ Basin | poly(log10Flux, 2, raw = TRUE), adjust = "holm" )
    
  #' Finally, as most of the action seems to be in the biggest mass group, collapse the smaller levels
  #' So that we have a contrast between small (first 3 mass quartiles) and big (top mass quartile):
      b_spxst_q  <- b_spxst_q %>%
        mutate(mass_cat = ifelse(M_quant == "0.75", "top", "small"))
      # 're-fit the model
      b_h_flux_temp_mod_lme2 <- lmer(H ~ poly(log10Flux, 2, raw = TRUE) +
                                       mass_cat +
                                       mass_cat:poly(log10Flux, 2, raw = TRUE) +
                                       (1|Basin),
                                     data = filter(b_spxst_q, n_sp > 1),
                                     na.action = "na.fail")
      summary(b_h_flux_temp_mod_lme2)
      
      #' look at on plot:
      div_v_env_plot(dat = filter(b_spxst_q, n_sp > 1), facet_var = "mass_cat")
      #' Hump-shaped ralationship for smaller species turns to ~flat (possible U-shaped) for biggest species.
      ggsave("Figure4.pdf")
      
      
   
    
############################### Proportion of Species Plots and Analyses ###############################

#' We can also consider what proportion of individuals or species are from each mass group,
#' and how this varies with flux. First, get total number of individuals and species per site,
#' add to our quantile data, and calculate relevant proportions:
    total_ind_stat <- b_spxst %>% group_by(station) %>%
      summarise(total_ind = sum(abundance),
                total_sp = n_distinct(species))
    b_spxst_q <- left_join(b_spxst_q, total_ind_stat, by = c("site" = "station"))
    b_spxst_q <- b_spxst_q %>% mutate(p_ind = n_ind / total_ind, p_sp = n_sp / total_sp)


#Plots for proportion of species:
    div_v_env_plot(div_var = "p_sp", w_var = "total_sp", d_lims = c(0, 1))+
      ylab("Proportion of Total Species")
    ggsave("Figure 5.pdf")


    #' Try the same for P(species) - here just model the proportion 0.75
    b_psp_flux_temp_mod <- glm(p_sp ~ poly(log10Flux, 2, raw = TRUE) + Basin,
                               weights = total_sp,
                               data = filter(b_spxst_q, M_quant == "0.75"),
                               family = binomial,
                               na.action = "na.fail")
    
    anova(b_psp_flux_temp_mod, test = "Chisq")
    summary(b_psp_flux_temp_mod, test = "Chisq")
    
    #' Try the same for P(species) - here just model the proportion 0.5
    b_psp_flux_temp_mod <- glm(p_sp ~ poly(log10Flux, 2, raw = TRUE) + Basin,
                               weights = total_sp,
                               data = filter(b_spxst_q, M_quant == "0.5"),
                               family = binomial,
                               na.action = "na.fail")
    
    anova(b_psp_flux_temp_mod, test = "Chisq")
    summary(b_psp_flux_temp_mod, test = "Chisq")
    
    #' Try the same for P(species) - here just model the proportion 0.25
    b_psp_flux_temp_mod <- glm(p_sp ~ poly(log10Flux, 2, raw = TRUE) + Basin,
                               weights = total_sp,
                               data = filter(b_spxst_q, M_quant == "0.25"),
                               family = binomial,
                               na.action = "na.fail")
    
    anova(b_psp_flux_temp_mod, test = "Chisq")
    summary(b_psp_flux_temp_mod, test = "Chisq")
    
    #' Try the same for P(species) - here just model the proportion 0
    b_psp_flux_temp_mod <- glm(p_sp ~ poly(log10Flux, 2, raw = TRUE) + Basin,
                               weights = total_sp,
                               data = filter(b_spxst_q, M_quant == "0"),
                               family = binomial,
                               na.action = "na.fail")
    
    anova(b_psp_flux_temp_mod, test = "Chisq")
    summary(b_psp_flux_temp_mod, test = "Chisq")
    
#comparing proportion of species 
    b_spxst_q_spread <- b_spxst_q %>%
      group_split(M_quant)
    
    b_spxst_q_spread2 <- b_spxst_q_spread[[1]] %>%
      full_join(b_spxst_q_spread[[2]], by='Station') %>%
      full_join(b_spxst_q_spread[[3]], by='Station') %>%
      full_join(b_spxst_q_spread[[4]], by='Station')

    #plots
      p0_25 <- ggplot(data=b_spxst_q_spread2, aes(x=p_sp.x, y=p_sp.y))+
        geom_point(alpha = .4, cex=3)+
        geom_smooth(method = "lm", color="grey20")+
        xlab("Mass Quantile 0")+
        ylab("Mass Quantile 0.25")+
        theme_craig()+
        stat_cor(method = "pearson", label.x = .1, label.y = .9)+
        xlim(c(0,1))+
        ylim(c(0,1))
      
      p0_50 <- ggplot(data=b_spxst_q_spread2, aes(x=p_sp.x, y=p_sp.x.x))+
        geom_point(alpha = .4, cex=3)+
        geom_smooth(method = "lm", color="grey20")+
        xlab("Mass Quantile 0")+
        ylab("Mass Quantile 0.50")+
        theme_craig()+
        stat_cor(method = "pearson", label.x = .1, label.y = .9)+
        xlim(c(0,1))+
        ylim(c(0,1))
      
      p0_75 <- ggplot(data=b_spxst_q_spread2, aes(x=p_sp.x, y=p_sp.y.y))+
        geom_point(alpha = .4, cex=3)+
        geom_smooth(method = "lm", color="grey20")+
        xlab("Mass Quantile 0")+
        ylab("Mass Quantile 0.75")+
        theme_craig()+
        stat_cor(method = "pearson", label.x = .1, label.y = .9)+
        xlim(c(0,1))+
        ylim(c(0,1))
      
      p25_50 <- ggplot(data=b_spxst_q_spread2, aes(x=p_sp.y, y=p_sp.x.x))+
        geom_point(alpha = .4, cex=3)+
        geom_smooth(method = "lm", color="grey20")+
        xlab("Mass Quantile 0.25")+
        ylab("Mass Quantile 0.50")+
        theme_craig()+
        stat_cor(method = "pearson", label.x = .1, label.y = .9)+
        xlim(c(0,1))+
        ylim(c(0,1))
      
      p25_75 <- ggplot(data=b_spxst_q_spread2, aes(x=p_sp.y, y=p_sp.y.y))+
        geom_point(alpha = .4, cex=3)+
        xlab("Mass Quantile 0.25")+
        ylab("Mass Quantile 0.75")+
        theme_craig()+
        stat_cor(method = "pearson", label.x = .1, label.y = .9)+
        xlim(c(0,1))+
        ylim(c(0,1))
      
      p50_75 <- ggplot(data=b_spxst_q_spread2, aes(x=p_sp.x.x, y=p_sp.y.y))+
        geom_point(alpha = .4, cex=3)+
        xlab("Mass Quantile 0.50")+
        ylab("Mass Quantile 0.75")+
        theme_craig()+
        stat_cor(method = "pearson", label.x = .1, label.y = .9)+
        xlim(c(0,1))+
        ylim(c(0,1))
      
      blank <- grid.rect(gp=gpar(col="white"))
      
      pdf(file="figure6.pdf",width=8.5, height=11, useDingbats=FALSE)
        par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
        
        grid.arrange(p0_25, p0_50, p0_75,
                   blank, p25_50, p25_75,
                   blank, blank, p50_75)
      dev.off()
      
      

    
    
      
 