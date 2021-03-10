#Script datant du 1 janvier 2020, modifié pour l'article

library(ggplot2)
library(glmmTMB)
library(car)
library(DHARMa)
library(ggeffects)
library(emmeans)
library(multcomp)
library(vioplot)
library(vegan)
library(gridExtra)
library(buildmer)


library(see) # if using performance::check_model


library(lme4)
library(effects)

library(sjstats)
library(bbmle) #for AICtab
library(cowplot)
library(RColorBrewer)
library(ggpubr)

library(lmerTest)




                                ######### [PART 1] #########

                                        ## colors ##

colorCardinals <- c("#440154ff", "#35b779ff", "#fde725ff", "#31688eff")
# North: "#440154ff"
# East: "#35b779ff"
# South: "#fde725ff"
# West: "#31688eff"

colorInvasion <- c("#b3bde2ff", "#8e453cff")
# Invaded: "#8e453cff"
# Non-invaded: "#b3bde2ff"



                                ######### [PART 2] #########

                                    ## Load datasets ##

# Baits' data
dataTI <- read.table("Data_impactsTapinoma_vA.txt", h=T, sep="\t") 
head(dataTI)
dim(dataTI)

# Buildings' data (at the bait level >> lot of duplicated values, that's normal)
dataBuildings <- read.table("Data_Buildings_v3.txt", h=T, sep="\t")
head(dataBuildings)
dim(dataBuildings)

                                ######### [PART 3] #########
                    ## Data manipulation and preparation for analyses ##

# merging baits and buildings data
dataTI2 <- cbind(dataTI, dataBuildings[,7:16][match(dataTI$ID, dataBuildings$ID),])
head(dataTI2)  
dim(dataTI2)   # 3840 30    ## 23 + 10 = 33


# Create new column indicating if the bait was near (A position) or far (B position) from the building
dataTI2$positionAB <- rep(c("A","B"), 1920)

for (i in 1:dim(dataTI2)[1]){
  if (grepl("A",dataTI2$position[i])){ 
    dataTI2$positionAB[i] <- "A"
  } else {dataTI2$positionAB[i] <- "B"}
}                                        

# re-ordering factors levels
dataTI2$time <- as.factor(dataTI2$time)
dataTI2$time = factor(dataTI2$time, levels(dataTI2$time)[c(2,3,1)]) 
dataTI2$Bface <- as.factor(dataTI2$Bface)
dataTI2$Bface = factor(dataTI2$Bface, levels(dataTI2$Bface)[c(2,1,3,4)]) 
dataTI2$Date <- as.factor(dataTI2$Date)
dataTI2$building <- as.factor(dataTI2$building)
dataTI2$building = factor(dataTI2$building, levels(dataTI2$building)[c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]) 

# Sampling design (number of baits by sampling events and day at which each sampling event was performed, by building)
table(dataTI2$Date, dataTI2$time, dataTI2$building)

# Total number of native ant workers recruited at each bait
dataTI2$Ab_allNatives <- rowSums(dataTI2[ ,10:23])
head(dataTI2)
dim(dataTI2)


# changing the ID column with a new IDbait column (without tuna an honey information)
dataTI2$IDbait <- as.factor(gsub('.{4}$', '', dataTI2$ID))
dataTI2 <- dataTI2[,2:36]
head(dataTI2)
dim(dataTI2)

## merging tuna and honey lines by bait to have a bait level dataset
dataTI3_part1 <- dataTI2[seq(2,3840, by=2),c(35,31,32,33,1,2,3,4,5,7,23,24,26,28,29,30)] 
head(dataTI3_part1)
dim(dataTI3_part1) #1920  16

# On sélectionne ensuite les variables quantitatives et on sum les données des deux lignes (miel et thon) 
dataTI3_part2 <- aggregate(dataTI2[,c(8:22,34)], by=list(dataTI2$IDbait), sum)  # De base: c(10:20,29,30)
head(dataTI3_part2)
dim(dataTI3_part2) #1920  18
colnames(dataTI3_part2)[1] <- "IDbait"

dataTI3 <- cbind(dataTI3_part1, dataTI3_part2[match(dataTI3_part1$IDbait, dataTI3_part2$IDbait),][,-1])
dataTI3$presence_Natives <- 0
for (i in 1:1920){if(dataTI3$Ab_allNatives[i] > 0){dataTI3$presence_Natives[i] <- 1}}
head(dataTI3)


dataTI3$positionAB <- as.factor(dataTI3$positionAB)
dataTI3$zone <- as.factor(dataTI3$zone)
dataTI3$building <- as.factor(dataTI3$building)
dataTI3$position <- as.factor(dataTI3$position)
dataTI3$V_richness <- as.factor(dataTI3$V_richness)
head(dataTI3)
dim(dataTI3)   #1920  33


                                ######### [PART 4] #########
                  ## Microclimatic variations induced by shading conditions ##


dataTI3_Temperature <- dataTI3
dataTI3_Temperature$IDface <- paste0(dataTI3_Temperature$building, dataTI3_Temperature$time, dataTI3_Temperature$Bface)
dataTI3_Temperature <-  dataTI3_Temperature[!duplicated(dataTI3_Temperature$IDface),]
dim(dataTI3_Temperature)
head(dataTI3_Temperature)
hist(dataTI3_Temperature$Tground_mean)

mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$time=="morning"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$time=="noon"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$time=="afternoon"])

mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="N"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="E"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="S"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="W"])

dataTI3_Temperature$Light_mean_sqrt <- sqrt(dataTI3_Temperature$Light_mean*10764)
plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$time=="morning",])
plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$time=="noon",])
plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$time=="afternoon",])

plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="N",])
plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="E",])
plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="S",])
plot(Tground_mean ~ Light_mean_sqrt, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="W",])


Tmean_test <- glmmTMB(Tground_mean   ~ 1 +
                        Bface + 
                        time + 
                        #zone + 
                        time:Bface + 
                        #zone:time +
                        #zone:Bface +
                        (1|building) + (1|Date), 
                      data=dataTI3_Temperature, 
                      family=gaussian)
Anova(Tmean_test, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: Tground_mean
# Chisq Df Pr(>Chisq)    
# (Intercept) 373.4206  1  < 2.2e-16 ***
#   Bface         5.0495  3     0.1682    
#   time         23.2313  2  9.024e-06 ***
#   Bface:time   27.9366  6  9.658e-05 ***

summary(Tmean_test)
# Family: gaussian  ( identity )
# Formula:          Tground_mean ~ 1 + Bface + time + time:Bface + (1 | building) +      (1 | Date)
# Data: dataTI3_Temperature
# 
# AIC      BIC   logLik deviance df.resid 
# 1042.4   1091.3   -506.2   1012.4      177 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# building (Intercept) 1.030    1.015   
# Date     (Intercept) 4.730    2.175   
# Residual             9.636    3.104   
# Number of obs: 192, groups:  building, 16; Date, 9
# 
# Dispersion estimate for gaussian family (sigma^2): 9.64 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           21.3191     1.1032  19.324  < 2e-16 ***
#   BfaceE                 1.9438     1.0975   1.771 0.076543 .  
#   BfaceS                 1.2781     1.0975   1.165 0.244193    
#   BfaceW                -0.1313     1.0975  -0.120 0.904799    
#   timenoon               4.8842     1.1084   4.407 1.05e-05 ***
#   timeafternoon          4.2562     1.0975   3.878 0.000105 ***
#   BfaceE:timenoon        4.5437     1.5521   2.927 0.003417 ** 
#   BfaceS:timenoon        5.6062     1.5521   3.612 0.000304 ***
#   BfaceW:timenoon        5.2937     1.5521   3.411 0.000648 ***
#   BfaceE:timeafternoon  -0.7938     1.5521  -0.511 0.609060    
#   BfaceS:timeafternoon   3.7969     1.5521   2.446 0.014432 *  
#   BfaceW:timeafternoon   3.0688     1.5521   1.977 0.048019 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

system.time(sr_Tmean_test <- simulateResiduals(Tmean_test, n=1000))
testDispersion(simulationOutput = sr_Tmean_test, alternative ="two.sided")
plot(sr_Tmean_test)

performance::r2_nakagawa(Tmean_test)
# Conditional R2: 0.71
# Marginal R2: 0.53

### Plotting effects
ef_Tmean_test <- ggemmeans(Tmean_test, c("time", "Bface"), type = "fe")
plot(ef_Tmean_test, line.size=1.25, col=colorCardinals) + theme_classic()


# multiple comparisons
contrasts_Tmean_test <- emmeans(Tmean_test, specs = pairwise ~ time + Bface, type = "response" )
cld(contrasts_Tmean_test$emmeans, Letters = letters)
# time      Bface emmean  SE  df lower.CL upper.CL .group 
# morning   W       21.2 1.1 177     19.0     23.4  a     
# morning   N       21.3 1.1 177     19.1     23.5  a     
# morning   S       22.6 1.1 177     20.4     24.8  ab    
# morning   E       23.3 1.1 177     21.1     25.4  abc   
# afternoon N       25.6 1.1 177     23.4     27.8   bcd  
# noon      N       26.2 1.1 177     24.0     28.4   bcd  
# afternoon E       26.7 1.1 177     24.5     28.9    cd  
# afternoon W       28.5 1.1 177     26.3     30.7     de 
# afternoon S       30.7 1.1 177     28.5     32.8      ef
# noon      W       31.4 1.1 177     29.2     33.5      ef
# noon      E       32.7 1.1 177     30.5     34.9       f
# noon      S       33.1 1.1 177     30.9     35.3       f
# 
# Confidence level used: 0.95 
# P value adjustment: tukey method for comparing a family of 12 estimates 
# significance level used: alpha = 0.05 

# zone are similar in temperature
vioplot(Tground_mean ~ zone, data=dataTI3_Temperature, col=colorInvasion)

# additional plot/statistics
vioplot(Tground_mean ~ Bface, data=dataTI3_Temperature[dataTI3_Temperature$time=="morning",], col=colorCardinals)
vioplot(Tground_mean ~ Bface, data=dataTI3_Temperature[dataTI3_Temperature$time=="noon",], col=colorCardinals)
vioplot(Tground_mean ~ Bface, data=dataTI3_Temperature[dataTI3_Temperature$time=="afternoon",], col=colorCardinals)

vioplot(Tground_mean ~ time, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="N",])
vioplot(Tground_mean ~ time, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="E",])
vioplot(Tground_mean ~ time, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="S",])
vioplot(Tground_mean ~ time, data=dataTI3_Temperature[dataTI3_Temperature$Bface=="W",])

vioplot(Tground_mean ~ time, data=dataTI3_Temperature)
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$time=="morning"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$time=="noon"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$time=="afternoon"])

mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="N"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="E"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="S"])
mean(dataTI3_Temperature$Tground_mean[dataTI3_Temperature$Bface=="W"])






                                ######### [PART 5] #########
                      ## Effect of T. magnum on native ant communities ##

# Does Tapinoma influence the abundance, richness and diversity of native ant species?
dataTI3presence <- dataTI3
for (i in 1:1920) {for (j in 17:31) {if(dataTI3presence[i,j] > 0){dataTI3presence[i,j] <- 1}}} 
head(dataTI3presence)

# aggregate at the building level
dataTI3presence_v2 <- aggregate(dataTI3presence[,17:31], 
                                by = list(dataTI3presence$zone, 
                                          dataTI3presence$building, 
                                          dataTI3presence$time,
                                          dataTI3presence$Date,
                                          dataTI3presence$Age_building), 
                                sum)
colnames(dataTI3presence_v2)[1:5] <- c("zone", "building", "time", "Date", "age_building")
head(dataTI3presence_v2)
dim(dataTI3presence_v2)  # 48 18
dataTI3presence_v2 <- dataTI3presence_v2[order(dataTI3presence_v2$zone, dataTI3presence_v2$building, dataTI3presence_v2$time),]

dataTI3presence_v2$nbBaits <- 40
dataTI3presence_v2$allNative <-  rowSums(dataTI3presence_v2[,7:20])
dataTI3presence_v2$allANTS <-  rowSums(dataTI3presence_v2[,6:20])



# Species richness by sampling event
dataTI3presence_v2$richness <- NA
for (i in 1:dim(dataTI3presence_v2)[1]){
  dataTI3presence_v2$richness[i] = 14 - length(which(dataTI3presence_v2[i,7:20]==0))
}
hist(dataTI3presence_v2$richness)

# Shannon diversity by sampling event
dataTI3presence_v2$diversity <- diversity(dataTI3presence_v2[,7:20], index = "shannon")
hist(dataTI3presence_v2$diversity)
head(dataTI3presence_v2)
dim(dataTI3presence_v2)

vioplot(richness ~ zone, data= dataTI3presence_v2, ylim=c(0,10))
vioplot(diversity ~ zone, data= dataTI3presence_v2, ylim=c(0,2))
vioplot(allNative ~ zone, data= dataTI3presence_v2, ylim=c(0,40))

library(randomcoloR)
library(viridis)
dataTI3presence_v2$building
colorBuilding <- c(rep(distinctColorPalette(8), each=3),rep(distinctColorPalette(8), each=3))
shapeBuilding <- c(rep(c(0:7), each=3), rep(c(0:7), each=3))


# native ants richness
richness_test <- glmmTMB(richness ~ 1 +
                           zone + 
                           (1|building) + 
                           (1|Date), 
                         data=dataTI3presence_v2, 
                         family=gaussian)
summary(richness_test)
Anova(richness_test, type='III')

system.time(sr_richness_test <- simulateResiduals(richness_test, n=1000))
testDispersion(simulationOutput = sr_richness_test, alternative ="two.sided")
plot(sr_richness_test) # very good...

p_richness <- ggplot(dataTI3presence_v2, aes(y=richness, x=zone, fill=building)) + 
  geom_boxplot(width=0.5, alpha=0.9, fill=colorInvasion, show.legend = FALSE) +
  geom_point(position=position_jitter(width = 0.25, height = 0, seed=3), shape=shapeBuilding, size=2.5, stroke=1.5, aes(color=zone), show.legend = FALSE) +
  geom_text(position=position_jitter(width = 0.25, height = 0,seed=3), aes(label=building), size=2.5) +
  scale_color_manual(values = colorInvasion) +
  ylim(c(0,10)) +
  scale_y_continuous( breaks=c(0,1,2,3,4,5,6,7,8,9), labels=c(0,1,2,3,4,5,6,7,8,9), limits=c(0.75,9.25)) +
  theme_classic()
p_richness


# native ants diversity
diversity_test <- glmmTMB(diversity ~ 1 +
                           zone + 
                           (1|building) + 
                           (1|Date), 
                         data=dataTI3presence_v2, 
                         family=gaussian)
summary(diversity_test)
Anova(diversity_test, type='III')

system.time(sr_diversity_test <- simulateResiduals(diversity_test, n=1000))
testDispersion(simulationOutput = sr_diversity_test, alternative ="two.sided")
plot(sr_diversity_test) # very good...

p_diversity <- ggplot(dataTI3presence_v2, aes(y=diversity, x=zone, fill=building)) + 
  geom_boxplot(width=0.5, alpha=0.9, fill=colorInvasion, show.legend = FALSE) + 
  geom_point(position=position_jitter(width = 0.25, height = 0, seed=3), shape=shapeBuilding, size=2.5, stroke=1.5, aes(color=zone), show.legend = FALSE) +
  geom_text(position=position_jitter(width = 0.25, height = 0,seed=3), aes(label=building), size=2.5) +
  scale_color_manual(values = colorInvasion) +
  ylim(c(0,10)) +
  scale_y_continuous( breaks=c(0.5,0.75,1,1.25,1.5,1.75), labels=c(0.5,0.75,1,1.25,1.5,1.75), limits=c(0.4,1.85)) +
  theme_classic()
p_diversity

# native ants abundance (proportion of baits occupied)
allNative_test <- glmmTMB(cbind(allNative, nbBaits) ~ 1 +
                            zone + 
                            (1|building) + 
                            (1|Date), 
                          data=dataTI3presence_v2, 
                          family=binomial)
summary(allNative_test)
Anova(allNative_test, type='III')

system.time(sr_allNative_test <- simulateResiduals(allNative_test, n=1000))
testDispersion(simulationOutput = sr_allNative_test, alternative ="two.sided")
plot(sr_allNative_test) # very good...

p_abundance <- ggplot(dataTI3presence_v2, aes(y=allNative, x=zone, fill=building)) + 
  geom_boxplot(width=0.5, alpha=0.9, fill=colorInvasion, show.legend = FALSE) + 
  geom_point(position=position_jitter(width = 0.25, height = 0, seed=3), shape=shapeBuilding, size=2.5, stroke=1.5, aes(color=zone), show.legend = FALSE) +
  geom_text(position=position_jitter(width = 0.25, height = 0,seed=3), aes(label=building), size=2.5) +
  scale_color_manual(values = colorInvasion) +
  ylim(c(0,42)) +
  scale_y_continuous( breaks=c(0,8,16,24,32,40), labels= round(c(0,8,16,24,32,40)/40*100)) +
  theme_classic()
p_abundance





figure2 <- function(){
  grid.arrange(p_richness, p_abundance, p_diversity, nrow = 1)
}

pdf(file = "figure2_v1.pdf", width=15, height = 5)
figure2()
dev.off()



plot(las_nig ~tap_mag, data=dataTI3presence_v2)


# native + tapinoma magnum ants abundance (proportion of baits occupied)
allANTS_test <- glmmTMB(cbind(allANTS, nbBaits) ~ 1 +
                          zone + 
                          (1|building) + 
                          (1|Date), 
                        data=dataTI3presence_v2, 
                        family=binomial)
summary(allANTS_test)
Anova(allANTS_test, type='III')

# # Lasius niger abundance (proportion of baits occupied)
# lasnig_test <- glmmTMB(cbind(las_nig, nbBaits) ~ 1 +
#                          zone + 
#                          (1|building) + 
#                          (1|Date), 
#                        data=dataTI3presence_v2, 
#                        family=binomial)
# summary(lasnig_test)
# 
# # Myrmica specioides abundance (proportion of baits occupied)
# myrspe_test <- glmmTMB(cbind(myr_spe, nbBaits) ~ 1 +
#                          zone + 
#                          (1|building) + 
#                          (1|Date), 
#                        data=dataTI3presence_v2, 
#                        family=binomial)
# summary(myrspe_test)
# 
# # Myrmica sabuleti abundance (proportion of baits occupied)
# myrsab_test <- glmmTMB(cbind(myr_sab, nbBaits) ~ 1 +
#                          zone + 
#                          (1|building) + 
#                          (1|Date), 
#                        data=dataTI3presence_v2, 
#                        family=binomial)
# summary(myrsab_test)
# 
# # Tetramorium sp. abundance (proportion of baits occupied)
# tetsp_test <- glmmTMB(cbind(tet_sp, nbBaits) ~ 1 +
#                         zone + 
#                         (1|building) + 
#                         (1|Date), 
#                       data=dataTI3presence_v2, 
#                       family=binomial)
# summary(tetsp_test)




# Plotting native species species recruitment at baits

# new dataset to add some info
dataTI3_2 <- dataTI3
for (i in 17:31) {
  for (j in 1: dim(dataTI3_2)[1]){
    if (dataTI3_2[j,i] > 0) {dataTI3_2[j,i] <- 1}
  }
}
head(dataTI3_2)
dataTI3_2$Ab_allNatives[dataTI3_2$Ab_allNatives>0] <- 1
dataTI3_2$Ab_baits_allspecies <- as.numeric(as.character(dataTI3_2$Ab_allNatives)) + dataTI3_2$tap_mag
dataTI3_2$Ab_baits_allspecies <- as.numeric(as.character(dataTI3_2$Ab_baits_allspecies))
range(dataTI3_2$Ab_baits_allspecies)
dataTI3_2$Ab_baits_allspecies[dataTI3_2$Ab_baits_allspecies>0] <- 1

# creates dataset for plotting the barplot
ALLSpecies_occBaits <- aggregate(as.numeric(as.character(dataTI3_2$Ab_baits_allspecies)), by=list(dataTI3_2$building), mean)
colnames(ALLSpecies_occBaits) <- c("building", "meanOccBaits")
ALLSpecies_occBaits$zone <- c(rep("Invaded", 8), rep("Non-Invaded", 8))
w = wilcox.test(ALLSpecies_occBaits$meanOccBaits ~ ALLSpecies_occBaits$zone)
w2 <- c(round(as.numeric(w$statistic, digits=0)), round(w$p.value, digits=4))
mean_species_occBaits <- aggregate(ALLSpecies_occBaits$meanOccBaits, by=list(ALLSpecies_occBaits$zone), mean)
colnames(mean_species_occBaits) <- c("zone", "meanOccBaits")
mean_species_occBaits$median <- aggregate(ALLSpecies_occBaits$meanOccBaits, by=list(ALLSpecies_occBaits$zone), median)[,2]
mean_species_occBaits$se <- aggregate(ALLSpecies_occBaits$meanOccBaits, by=list(ALLSpecies_occBaits$zone), sd)[,2]/sqrt(8)
mean_species_occBaits$species <- "All species"
mean_species_occBaits$statsWilcox <- w2
Stats_occBaits <- mean_species_occBaits

NativeSpecies_occBaits <- aggregate(as.numeric(as.character(dataTI3_2$Ab_allNatives)), by=list(dataTI3_2$building), mean)
colnames(NativeSpecies_occBaits) <- c("building", "meanOccBaits")
NativeSpecies_occBaits$zone <- c(rep("Invaded", 8), rep("Non-Invaded", 8))
w = wilcox.test(NativeSpecies_occBaits$meanOccBaits ~ NativeSpecies_occBaits$zone)
w2 <- c(round(as.numeric(w$statistic, digits=0)), round(w$p.value, digits=4))
mean_species_occBaits <- aggregate(NativeSpecies_occBaits$meanOccBaits, by=list(NativeSpecies_occBaits$zone), mean)
colnames(mean_species_occBaits) <- c("zone", "meanOccBaits")
mean_species_occBaits$median <- aggregate(NativeSpecies_occBaits$meanOccBaits, by=list(NativeSpecies_occBaits$zone), median)[,2]
mean_species_occBaits$se <- aggregate(NativeSpecies_occBaits$meanOccBaits, by=list(NativeSpecies_occBaits$zone), sd)[,2]/sqrt(8)
mean_species_occBaits$species <- "Native only"
mean_species_occBaits$statsWilcox <- w2
Stats_occBaits <- rbind(Stats_occBaits, mean_species_occBaits)


for (i in 17:31){
  species_occBaits <- aggregate(dataTI3_2[,i], by=list(dataTI3_2$building), mean)
  colnames(species_occBaits) <- c("building", "meanOccBaits")
  species_occBaits$zone <- c(rep("Invaded", 8), rep("Non-Invaded", 8))
  w = wilcox.test(species_occBaits$meanOccBaits ~ species_occBaits$zone)
  w2 <- c(round(as.numeric(w$statistic, digits=0)), round(w$p.value, digits=4))
  
  mean_species_occBaits <- aggregate(species_occBaits$meanOccBaits, by=list(species_occBaits$zone), mean)
  colnames(mean_species_occBaits) <- c("zone", "meanOccBaits")
  mean_species_occBaits$median <- aggregate(species_occBaits$meanOccBaits, by=list(species_occBaits$zone), median)[,2]
  mean_species_occBaits$se <- aggregate(species_occBaits$meanOccBaits, by=list(species_occBaits$zone), sd)[,2]/sqrt(8)
  mean_species_occBaits$species <- colnames(dataTI3_2)[i]
  mean_species_occBaits$statsWilcox <- w2
  
  Stats_occBaits <- rbind(Stats_occBaits, mean_species_occBaits)
}

order_NI <- Stats_occBaits[Stats_occBaits$zone=="Non-Invaded",]
order_species_NI <- as.character(order_NI$species[order(order_NI$meanOccBaits, decreasing=T) ])

# plot
Stats_occBaits$zone <- ordered(Stats_occBaits$zone, levels = c("Non-Invaded","Invaded"))
Stats_occBaits$species <- ordered(Stats_occBaits$specie, levels = order_species_NI)


figure3 <- ggplot(data=Stats_occBaits[-c(1:4),], aes(x=species, y=meanOccBaits, fill=zone)) +
  geom_bar(stat = "identity", position="dodge",lwd=0.65, col="black", width=0.75) +
  geom_errorbar( aes(x=species, ymin=(meanOccBaits-se), ymax=(meanOccBaits+se)), 
                 position=position_dodge(0.75),
                 width=0.25, colour="black", alpha=0.9, size=0.5) +
  geom_point(aes(x=species, y=median, fill=zone),
             position=position_dodge(0.75),
             color= "red", size=1.25) +
  labs(y="Proportion of baits occupied",
       x="Species") +
  scale_fill_manual(values=colorInvasion) + 
  scale_y_continuous(trans = 'pseudo_log') +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1.15),
        axis.text.y = element_text(size=14))

pdf(file = "figure3_v1.pdf", width=10, height = 7.5)
figure3
dev.off()

#



                                       ######### [PART 6] #########
            ## Effect of T. magnum and shading conditions on native ants' foraging activity ##

## Probability of presence at baits

head(dataTI3presence)
dim(dataTI3presence)

# Lasius niger
binomial_lasnig <- glmmTMB(las_nig ~ 1 +
                          zone + 
                          Bface +
                          time +
                          zone:Bface +
                          #zone:time +
                          time:Bface +
                          (1|building) +
                          (1|Date), 
                        data=dataTI3presence, 
                        family=binomial)
Anova(binomial_lasnig, type="III")
summary(binomial_lasnig)
r2_tjur(binomial_lasnig)

system.time(sr_binomial_lasnig <- simulateResiduals(binomial_lasnig, n=1000))
testDispersion(simulationOutput = sr_binomial_lasnig, alternative ="two.sided")
plot(sr_binomial_lasnig) # very good...

ef_binomial_lasnig_all <- ggemmeans(binomial_lasnig, c("time", "Bface", "zone"), type = "fe")
plot(ef_binomial_lasnig_all, col=colorCardinals, line.size=1.75, dot.size=4, dodge=0.4) 

ef_binomial_lasnig_side <- ggemmeans(binomial_lasnig, c("Bface", "zone"), type = "fe")
plot(ef_binomial_lasnig_side, col=colorInvasion, line.size=2, dot.size=5, dodge=0.25)

pdf(file = "figure4a_v1.pdf", width=10, height = 7)
plot(ef_binomial_lasnig_side, col=colorInvasion, line.size=2, dot.size=5, dodge=0.25)
dev.off()

marginal_binomial_lasnig_side <- emmeans(binomial_lasnig, specs =  pairwise ~ zone | Bface, type = "response", adjust="tukey") 
cld(marginal_binomial_lasnig_side$emmeans,
    alpha=0.05,
    Letters=letters)

marginal_binomial_lasnig_zone <- emmeans(binomial_lasnig, specs =  pairwise ~ zone, type = "response", adjust="tukey") 
marginal_binomial_lasnig_zone

binomial_lasnig_nre <- glmmTMB(las_nig ~ 1 +
                             zone + 
                             Bface +
                             time +
                             zone:Bface +
                             #zone:time +
                             time:Bface, 
                           data=dataTI3presence, 
                           family=binomial)
r2_tjur(binomial_lasnig_nre)


# Myrmica specioides
binomial_myrspe <- glmmTMB(myr_spe ~ 1 +
                             zone + 
                             #Bface +
                             time +
                             #zone:Bface +
                             #zone:time +
                             #time:Bface +
                             (1|building) +
                             (1|Date), 
                           data=dataTI3presence, 
                           family=binomial)
Anova(binomial_myrspe, type="III")
summary(binomial_myrspe)
r2_tjur(binomial_myrspe)

system.time(sr_binomial_myrspe <- simulateResiduals(binomial_myrspe, n=1000))
testDispersion(simulationOutput = sr_binomial_myrspe, alternative ="two.sided")
plot(sr_binomial_myrspe) # very good...

ef_binomial_myrspe_side <- ggemmeans(binomial_myrspe, c("time", "zone"), type = "fe")
plot(ef_binomial_myrspe_side, col=colorInvasion, line.size=2, dot.size=5, dodge=0.25)# + ylim(0.15, 0.62)

marginal_binomial_myrspe_time <- emmeans(binomial_myrspe, specs =  pairwise ~ zone | time, type = "response", adjust="tukey") 
cld(marginal_binomial_myrspe_time$emmeans,
    alpha=0.05,
    Letters=letters)

marginal_binomial_myrspe_zone <- emmeans(binomial_myrspe, specs =  pairwise ~ zone, type = "response", adjust="tukey") 
marginal_binomial_myrspe_zone

binomial_myrspe_nre <- glmmTMB(myr_spe ~ 1 +
                                 zone + 
                                 #Bface +
                                 time 
                               #zone:Bface +
                               #zone:time +
                               #time:Bface +
                               , 
                               data=dataTI3presence, 
                               family=binomial)
r2_tjur(binomial_myrspe_nre)

# Myrmica sabuleti
binomial_myrsab <- glmmTMB(myr_sab ~ 1 +
                             zone + 
                             Bface +
                             time +
                             #zone:Bface +
                             #zone:time +
                             #time:Bface +
                             (1|building) +
                             (1|Date), 
                           data=dataTI3presence, 
                           family=binomial)
Anova(binomial_myrsab, type="III")
summary(binomial_myrsab)
r2_tjur(binomial_myrsab)

system.time(sr_binomial_myrsab <- simulateResiduals(binomial_myrsab, n=1000))
testDispersion(simulationOutput = sr_binomial_myrsab, alternative ="two.sided")
plot(sr_binomial_myrsab) # very good...

ef_binomial_myrsab_all <- ggemmeans(binomial_myrsab, c("time", "Bface", "zone"), type = "fe")
plot(ef_binomial_myrsab_all, col=colorCardinals, line.size=2, dot.size=5, dodge=0.25)# + ylim(0.15, 0.62)

marginal_binomial_myrsab_zone <- emmeans(binomial_myrsab, specs =  pairwise ~ zone, type = "response", adjust="tukey") 
marginal_binomial_myrsab_zone

binomial_myrsab_nre <- glmmTMB(myr_sab ~ 1 +
                             zone + 
                             Bface +
                             time 
                             #zone:Bface +
                             #zone:time +
                             #time:Bface +
                               , 
                           data=dataTI3presence, 
                           family=binomial)
r2_tjur(binomial_myrsab_nre)

# Tetramorium sp.
binomial_tetsp <- glmmTMB(tet_sp ~ 1 +
                             #zone + 
                             Bface +
                             #time +
                             #zone:Bface +
                             #zone:time +
                             #time:Bface +
                             (1|building) +
                             (1|Date), 
                           data=dataTI3presence, 
                           family=binomial)
Anova(binomial_tetsp, type="III")
summary(binomial_tetsp)
r2_tjur(binomial_tetsp)

system.time(sr_binomial_tetsp <- simulateResiduals(binomial_tetsp, n=1000))
testDispersion(simulationOutput = sr_binomial_tetsp, alternative ="two.sided")
plot(sr_binomial_tetsp) # very good...

ef_binomial_tetsp_side <- ggemmeans(binomial_tetsp, c("Bface"), type = "fe")
plot(ef_binomial_tetsp_side, line.size=2, col=colorInvasion, dot.size=5, dodge=0.25)# + ylim(0.15, 0.62)

marginal_binomial_tetsp_side <- emmeans(binomial_tetsp, specs =  pairwise ~ Bface, type = "response", adjust="tukey") 
cld(marginal_binomial_tetsp_side$emmeans,
    alpha=0.05,
    Letters=letters)


binomial_tetsp_nre <- glmmTMB(tet_sp ~ 1 +
                            #zone + 
                            Bface 
                            #time +
                            #zone:Bface +
                            #zone:time +
                            #time:Bface +
                            , 
                          data=dataTI3presence, 
                          family=binomial)
r2_tjur(binomial_tetsp_nre)

# Tapinoma magnum
binomial_tapmag <- glmmTMB(tap_mag ~ 1 +
                             Bface +
                             #time +
                             #time:Bface +
                             (1|building) +
                             (1|Date), 
                           data=dataTI3presence[dataTI3presence$zone=="invaded",], 
                           family=binomial)
Anova(binomial_tapmag, type="III")
summary(binomial_tapmag)
r2_tjur(binomial_tapmag)

system.time(sr_binomial_tapmag <- simulateResiduals(binomial_tapmag, n=1000))
testDispersion(simulationOutput = sr_binomial_tapmag, alternative ="two.sided")
plot(sr_binomial_tapmag) # very good...

ef_binomial_tapmag_side <- ggemmeans(binomial_tapmag, c("Bface"), type = "fe")
plot(ef_binomial_tapmag_side, col=colorCardinals, line.size=2, dot.size=5, dodge=0.25)

marginal_binomial_tapmag_side <- emmeans(binomial_tapmag, specs =  pairwise ~ Bface, type = "response", adjust="tukey") 
cld(marginal_binomial_tapmag_side$emmeans,
    alpha=0.05,
    Letters=letters)

binomial_tapmag_nre <- glmmTMB(tap_mag ~ 1 +
                             Bface 
                             #time +
                             #time:Bface +
                            , 
                           data=dataTI3presence[dataTI3presence$zone=="invaded",], 
                           family=binomial)
r2_tjur(binomial_tapmag_nre)



## Number of workers recruited at baits

# Lasius niger >> DONE
poisson_lasnig <- glmmTMB(las_nig ~ 1 + 
                       zone +
                       time + 
                       Bface +
                       zone:time +
                       zone:Bface +
                       #time:Bface +
                       (1|building)  + 
                       (1|Date),
                     
                     ziformula= ~ 1 +
                       zone +
                       time +
                       Bface +
                       #zone:time +
                       zone:Bface +
                       time:Bface +
                       (1|building) + 
                       (1|Date),
                     
                     family=truncated_nbinom1, 
                     data=dataTI3)

Anova(poisson_lasnig, type="III")
summary(poisson_lasnig)
performance::r2_zeroinflated(poisson_lasnig, method="correlation")

system.time(sr_poisson_lasnig <- simulateResiduals(poisson_lasnig, n=1000))
testDispersion(simulationOutput = sr_poisson_lasnig, alternative ="two.sided")
plot(sr_poisson_lasnig) # very good...


# all effects in Lasius niger
ef_poisson_lasnig_all <- ggemmeans(poisson_lasnig, c("time", "Bface", "zone"), type = "fe")
plot(ef_poisson_lasnig_all, line.size=1.25, col=colorCardinals) 


# interaction between building side and invasion in Lasius niger
ef_poisson_lasnig_side <- ggemmeans(poisson_lasnig, c("Bface", "zone"), type = "fe", ci.lvl = 0.95)
plot(ef_poisson_lasnig_side, line.size=2, col=colorInvasion, dot.size=5, dodge=0.25 )# + ylim(0.15, 0.62)

pdf(file = "figure4b_v1.pdf", width=10, height = 7)
plot(ef_poisson_lasnig_side, line.size=2, col=colorInvasion, dot.size=5, dodge=0.25 )# + ylim(0.15, 0.62)
dev.off()

marginal_poisson_lasnig_side <- emmeans(poisson_lasnig, specs =  pairwise ~ zone | Bface, type = "response", adjust="tukey") 
cld(marginal_poisson_lasnig_side$emmeans,
    alpha=0.05,
    Letters=letters)

# interaction between time and invasion in Lasius niger
ef_poisson_lasnig_time <- ggemmeans(poisson_lasnig, c("time", "zone"), type = "fe")
plot(ef_poisson_lasnig_time, line.size=2, col=colorInvasion, dot.size=5, dodge=0.25)# + ylim(0.15, 0.62)

pdf(file = "figure4c_v1.pdf", width=10, height = 7)
plot(ef_poisson_lasnig_time, line.size=2, col=colorInvasion, dot.size=5, dodge=0.25)
dev.off()

marginal_poisson_lasnig_time <- emmeans(poisson_lasnig, specs =  pairwise ~ zone | time, type = "response", adjust="tukey") 
cld(marginal_poisson_lasnig_time$emmeans,
    alpha=0.05,
    Letters=letters)

plot(log(predict(poisson_lasnig, type="response")), log(dataTI3$las_nig))

poisson_lasnig_nre <- glmmTMB(las_nig ~ 1 +
                                zone +
                                time + 
                                Bface +
                                zone:time +
                                zone:Bface ,
                                #time:Bface +
                              
                              ziformula= ~ 1 +
                                zone +
                                time +
                                Bface +
                                #zone:time +
                                zone:Bface +
                                time:Bface,
                              family=truncated_nbinom1,
                              data=dataTI3)

performance::r2_zeroinflated(poisson_lasnig_nre, method="correlation")

# Myrmica specioides >> DONE
poisson_myrspe <- glmmTMB(myr_spe ~ 1 + 
                            #zone +
                            #time + 
                            Bface +
                            #zone:time +
                            #zone:Bface +
                            #time:Bface +
                            (1|building)  + 
                            (1|Date),
                          
                          ziformula= ~ 1 +
                            zone +
                            time +
                            #Bface +
                            #zone:time +
                            #zone:Bface +
                            #time:Bface +
                            (1|building) + 
                            (1|Date),
                          
                          family= truncated_nbinom1, 
                          data=dataTI3)

Anova(poisson_myrspe, type="III")
summary(poisson_myrspe)
performance::r2_zeroinflated(poisson_myrspe, method="correlation")

system.time(sr_poisson_myrspe <- simulateResiduals(poisson_myrspe, n=1000))
testDispersion(simulationOutput = sr_poisson_myrspe, alternative ="two.sided")
plot(sr_poisson_myrspe) # very good...

ef_poisson_myrspe_side <- ggemmeans(poisson_myrspe, c("Bface"), type = "fe")
plot(ef_poisson_myrspe_side, line.size=2, col=colorCardinals, dot.size=5, dodge=0.25)# + ylim(0.15, 0.62)

marginal_poisson_myrspe_side <- emmeans(poisson_myrspe, specs =  pairwise ~ Bface, type = "response", adjust="tukey") 
cld(marginal_poisson_myrspe_side$emmeans,
    alpha=0.05,
    Letters=letters)

plot(log(predict(poisson_myrspe, type="response")), log(dataTI3$myr_spe))

poisson_myrspe_nre <- glmmTMB(myr_spe ~ 1 +
                              #zone +
                              #time ,
                              Bface ,
                              #zone:time +
                              #zone:Bface +
                              #time:Bface +
                              
                              ziformula= ~ 1 +
                                zone +
                                time ,
                                #Bface +
                              #zone:time +
                              #zone:Bface +
                              #time:Bface +
                              family=truncated_nbinom1,
                              data=dataTI3)

performance::r2_zeroinflated(poisson_myrspe_nre, method="correlation")


# Myrmica sabuleti >> DONE
poisson_myrsab <- glmmTMB(myr_sab ~ 1 + 
                            #zone +
                            #time + 
                            #Bface +
                            #zone:time +
                            #zone:Bface +
                            #time:Bface +
                            (1|building)  + 
                            (1|Date),
                          
                          ziformula= ~ 1 +
                            zone +
                            time +
                            Bface +
                            #zone:time +
                            #zone:Bface +
                            #time:Bface +
                            (1|building) + 
                            (1|Date),
                          
                          family=truncated_nbinom1, 
                          data=dataTI3)

Anova(poisson_myrsab, type="III")
summary(poisson_myrsab)
performance::r2_zeroinflated(poisson_myrsab, method="correlation")

system.time(sr_poisson_myrsab <- simulateResiduals(poisson_myrsab, n=1000))
testDispersion(simulationOutput = sr_poisson_myrsab, alternative ="two.sided")
plot(sr_poisson_myrsab) # very good...

plot(log(predict(poisson_myrsab, type="response")), log(dataTI3$myr_sab))

poisson_myrsab_nre <- glmmTMB(myr_sab ~ 1 ,
                               #zone +
                               #time ,
                             #Bface +
                             #zone:time +
                             #zone:Bface +
                             #time:Bface +
                             
                             ziformula= ~ 1 +
                               zone +
                               time +
                               Bface,
                             #zone:time +
                             #zone:Bface +
                             #time:Bface +
                             family=truncated_nbinom1,
                             data=dataTI3)

performance::r2_zeroinflated(poisson_myrsab_nre, method="correlation")

# Tetramorium sp. >> DONE
poisson_tetsp <- glmmTMB(tet_sp ~ 1 + 
                            #zone +
                            time + 
                            #Bface +
                            #zone:time +
                            #zone:Bface +
                            #time:Bface +
                            (1|building)  + 
                            (1|Date),
                          
                          ziformula= ~ 1 +
                            #zone +
                            #time +
                            Bface +
                            #zone:time +
                            #zone:Bface +
                            #time:Bface +
                            (1|building) + 
                            (1|Date),
                          
                          family=truncated_nbinom1, 
                          data=dataTI3)

Anova(poisson_tetsp, type="III")
summary(poisson_tetsp)
performance::r2_zeroinflated(poisson_tetsp, method="correlation")

system.time(sr_poisson_tetsp <- simulateResiduals(poisson_tetsp, n=1000))
testDispersion(simulationOutput = sr_poisson_tetsp, alternative ="two.sided")
plot(sr_poisson_tetsp) # very good...

ef_poisson_tetsp_time <- ggemmeans(poisson_tetsp, c("time"), type = "fe")
plot(ef_poisson_tetsp_time, line.size=2, dot.size=5, dodge=0.25)# + ylim(0.15, 0.62)

marginal_poisson_tetsp_time <- emmeans(poisson_tetsp, specs =  pairwise ~ time, type = "response", adjust="tukey") 
cld(marginal_poisson_tetsp_time$emmeans,
    alpha=0.05,
    Letters=letters)

plot(log(predict(poisson_tetsp, type="response")), log(dataTI3$tet_sp))

poisson_tetsp_nre <- glmmTMB(tet_sp ~ 1 +
                               #zone +
                                time ,
                              #Bface +
                             #zone:time +
                             #zone:Bface +
                             #time:Bface +
                              
                              ziformula= ~ 1 +
                               #zone +
                                #time +
                                Bface,
                             #zone:time +
                             #zone:Bface +
                             #time:Bface +
                              family=truncated_nbinom1,
                              data=dataTI3)

performance::r2_zeroinflated(poisson_tetsp_nre, method="correlation")



# Tapinoma magnum
poisson_tapmag <- glmmTMB(tap_mag ~ 1 + 
                           time + 
                           #Bface +
                           #time:Bface +
                           (1|building)  + 
                           (1|Date),
                         
                         ziformula= ~ 1 +
                           #time +
                           Bface +
                           #time:Bface +
                           (1|building) + 
                           (1|Date),
                         
                         family=truncated_nbinom1,
                         data=dataTI3[dataTI3$zone=="invaded" & dataTI3$tap_mag <= max(sort(dataTI3$tap_mag[dataTI3$zone=="invaded"])[-c((960*0.95):960)]),])

Anova(poisson_tapmag, type="III")
summary(poisson_tapmag)
performance::r2_zeroinflated(poisson_tapmag, method="correlation")

system.time(sr_poisson_tapmag <- simulateResiduals(poisson_tapmag, n=1000))
testDispersion(simulationOutput = sr_poisson_tapmag, alternative ="two.sided")
plot(sr_poisson_tapmag) # very good...

ef_poisson_tapmag_time <- ggemmeans(poisson_tapmag, c("time"), type = "fe")
plot(ef_poisson_tapmag_time, line.size=2, dot.size=5, dodge=0.25)


marginal_poisson_tapmag_time <- emmeans(poisson_tapmag, specs =  pairwise ~ time , type = "response", adjust="tukey") 
cld(marginal_poisson_tapmag_time$emmeans,
    alpha=0.05,
    Letters=letters)

poisson_tapmag_nre <- glmmTMB(tap_mag ~ 1 + 
                            time ,
                            #Bface +
                            #time:Bface
                          
                          ziformula= ~ 1 +
                            #time +
                            Bface,
                            #time:Bface
                          family=truncated_nbinom1,
                          data=dataTI3[dataTI3$zone=="invaded" & dataTI3$tap_mag <= max(sort(dataTI3$tap_mag[dataTI3$zone=="invaded"])[-c((960*0.95):960)]),])

performance::r2_zeroinflated(poisson_tapmag_nre, method="correlation")

plot(log(predict(poisson_tapmag, type="response")), log(dataTI3$tap_mag[dataTI3$zone=="invaded" & dataTI3$tap_mag <= max(sort(dataTI3$tap_mag[dataTI3$zone=="invaded"])[-c((960*0.95):960)])]))



