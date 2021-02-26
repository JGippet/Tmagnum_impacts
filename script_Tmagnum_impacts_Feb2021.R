#Script datant du 1 janvier 2020, modifié pour l'article

library(lme4)
library(car)
library(glmmTMB)
library(effects)
library(ggeffects)
library(multcomp)
library(DHARMa)
library(sjstats)
library(bbmle) #for AICtab
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(lmerTest)
library(emmeans)



                                        ## colors ##

colorCardinals <- c("#440154ff", "#35b779ff", "#fde725ff", "#31688eff")
# North: "#440154ff"
# East: "#35b779ff"
# South: "#fde725ff"
# West: "#31688eff"

colorInvasion <- c("#b3bde2ff", "#8e453cff")
# Invaded: "#8e453cff"
# Non-invaded: "#b3bde2ff"




