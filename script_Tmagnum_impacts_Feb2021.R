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
head(dataTI3)
dim(dataTI3)   #1920  33






