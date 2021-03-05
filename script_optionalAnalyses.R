library(cooccur)


######### [PART 7] #########
## Co-occurrences patterns ##



head(dataTI3)
str(dataTI3)

dataTI3_free <- dataTI3[dataTI3$zone=="free",]
dataTI3_invaded <- dataTI3[dataTI3$zone=="invaded",]

table(dataTI3_invaded$presence_Natives)

# co-occurrence of Tapinoma magnum with other species
dataTI3_invaded$cooccurence_tapmag <- "alone"
for (i in 1:dim(dataTI3_invaded)[1]){
  if (dataTI3_invaded$tap_mag[i] >0){
    if(sum(dataTI3_invaded[i,17:31][-1]>0)){
      dataTI3_invaded$cooccurence_tapmag[i] <- colnames(dataTI3_invaded[i,17:31][-1])[which(dataTI3_invaded[i,17:31][-1]>0)]
    }
  }
}

table(dataTI3_invaded$cooccurence_tapmag)

# co-occurrence of Lasius niger with other species in non-invaded sites
dataTI3_free$cooccurence_lasnig <- "alone"
for (i in 1:dim(dataTI3_free)[1]){
  if (dataTI3_free$las_nig[i] >0){
    if(sum(dataTI3_free[i,17:31][-3]>0)){
      dataTI3_free$cooccurence_lasnig[i] <- colnames(dataTI3_free[i,17:31][-3])[which(dataTI3_free[i,17:31][-3]>0)]
    }
  }
}

# co-occurrence of Lasius niger with other species in invaded sites
dataTI3_invaded$cooccurence_lasnig <- "alone"
for (i in 1:dim(dataTI3_invaded)[1]){
  if (dataTI3_invaded$las_nig[i] >0){
    if(sum(dataTI3_invaded[i,17:31][-3]>0)){
      dataTI3_invaded$cooccurence_lasnig[i] <- colnames(dataTI3_invaded[i,17:31][-3])[which(dataTI3_invaded[i,17:31][-3]>0)]
    }
  }
}

# co-occurrence of Myrmica specioides with other species in non-invaded sites
dataTI3_free$cooccurence_myrspe <- "alone"
for (i in 1:dim(dataTI3_free)[1]){
  if (dataTI3_free$myr_spe[i] >0){
    if(sum(dataTI3_free[i,17:31][-5]>0)){
      dataTI3_free$cooccurence_myrspe[i] <- colnames(dataTI3_free[i,17:31][-5])[which(dataTI3_free[i,17:31][-5]>0)]
    }
  }
}

# co-occurrence of Myrmica specioides with other species in invaded sites
dataTI3_invaded$cooccurence_myrspe <- "alone"
for (i in 1:dim(dataTI3_invaded)[1]){
  if (dataTI3_invaded$myr_spe[i] >0){
    if(sum(dataTI3_invaded[i,17:31][-5]>0)){
      dataTI3_invaded$cooccurence_myrspe[i] <- colnames(dataTI3_invaded[i,17:31][-5])[which(dataTI3_invaded[i,17:31][-5]>0)]
    }
  }
}

# co-occurrence of Myrmica sabuleti with other species in non-invaded sites
dataTI3_free$cooccurence_myrsab <- "alone"
for (i in 1:dim(dataTI3_free)[1]){
  if (dataTI3_free$myr_sab[i] >0){
    if(sum(dataTI3_free[i,17:31][-6]>0)){
      dataTI3_free$cooccurence_myrsab[i] <- colnames(dataTI3_free[i,17:31][-6])[which(dataTI3_free[i,17:31][-6]>0)]
    }
  }
}

# co-occurrence of Myrmica sabuleti with other species in invaded sites
dataTI3_invaded$cooccurence_myrsab <- "alone"
for (i in 1:dim(dataTI3_invaded)[1]){
  if (dataTI3_invaded$myr_sab[i] >0){
    if(sum(dataTI3_invaded[i,17:31][-6]>0)){
      dataTI3_invaded$cooccurence_myrsab[i] <- colnames(dataTI3_invaded[i,17:31][-6])[which(dataTI3_invaded[i,17:31][-6]>0)]
    }
  }
}

# co-occurrence of Tetramorium with other species in non-invaded sites
dataTI3_free$cooccurence_tetsp <- "alone"
for (i in 1:dim(dataTI3_free)[1]){
  if (dataTI3_free$tet_sp[i] >0){
    if(sum(dataTI3_free[i,17:31][-9]>0)){
      dataTI3_free$cooccurence_tetsp[i] <- colnames(dataTI3_free[i,17:31][-9])[which(dataTI3_free[i,17:31][-9]>0)]
    }
  }
}

# co-occurrence of Tetramorium with other species in invaded sites
dataTI3_invaded$cooccurence_tetsp <- "alone"
for (i in 1:dim(dataTI3_invaded)[1]){
  if (dataTI3_invaded$tet_sp[i] >0){
    if(sum(dataTI3_invaded[i,17:31][-9]>0)){
      dataTI3_invaded$cooccurence_tetsp[i] <- colnames(dataTI3_invaded[i,17:31][-9])[which(dataTI3_invaded[i,17:31][-9]>0)]
    }
  }
}

table(dataTI3_invaded$cooccurence_tapmag)
length(dataTI3_invaded$tap_mag[dataTI3_invaded$tap_mag==0])
length(dataTI3_invaded$tap_mag[dataTI3_invaded$tap_mag!=0])

dataTI3cooccurences <- rbind(dataTI3_free, dataTI3_invaded)
head(dataTI3cooccurences)

dataTI3cooccurences$las_nig_presence <- 0
dataTI3cooccurences$myr_spe_presence <- 0
dataTI3cooccurences$myr_sab_presence <- 0
dataTI3cooccurences$tet_sp_presence <- 0
dataTI3cooccurences$tap_mag_presence <- 0
for (i in 1:dim(dataTI3cooccurences)[1]) {
  if (dataTI3cooccurences$las_nig[i]>0){dataTI3cooccurences$las_nig_presence[i] <- 1}
  if (dataTI3cooccurences$myr_spe[i]>0){dataTI3cooccurences$myr_spe_presence[i] <- 1}
  if (dataTI3cooccurences$myr_sab[i]>0){dataTI3cooccurences$myr_sab_presence[i] <- 1}
  if (dataTI3cooccurences$tet_sp[i]>0){dataTI3cooccurences$tet_sp_presence[i] <- 1}
  if (dataTI3cooccurences$tap_mag[i]>0){dataTI3cooccurences$tap_mag_presence[i] <- 1}
}

table(dataTI3cooccurences$cooccurence_lasnig, dataTI3cooccurences$zone, dataTI3cooccurences$las_nig_presence)
length(dataTI3cooccurences$las_nig[dataTI3cooccurences$las_nig>0 & dataTI3cooccurences$zone=="invaded"])


table(dataTI3cooccurences$cooccurence_myrspe, dataTI3cooccurences$zone, dataTI3cooccurences$myr_spe_presence)
length(dataTI3cooccurences$myr_spe[dataTI3cooccurences$myr_spe>0 & dataTI3cooccurences$zone=="invaded"])


table(dataTI3cooccurences$cooccurence_myrsab, dataTI3cooccurences$zone, dataTI3cooccurences$myr_sab_presence)
length(dataTI3cooccurences$myr_sab[dataTI3cooccurences$myr_sab>0 & dataTI3cooccurences$zone=="invaded"])


table(dataTI3cooccurences$cooccurence_tetsp, dataTI3cooccurences$zone, dataTI3cooccurences$tet_sp_presence)
length(dataTI3cooccurences$tet_sp[dataTI3cooccurences$tet_sp>0 & dataTI3cooccurences$zone=="invaded"])


table(dataTI3cooccurences$cooccurence_tapmag, dataTI3cooccurences$zone, dataTI3cooccurences$tap_mag_presence)
length(dataTI3cooccurences$tap_mag[dataTI3cooccurences$tap_mag>0 & dataTI3cooccurences$zone=="invaded"])



table(dataTI3_free$cooccurence_lasnig)
table(dataTI3_invaded$cooccurence_lasnig)
length(dataTI3_invaded$las_nig[dataTI3_invaded$las_nig==0])
length(dataTI3_invaded$las_nig[dataTI3_invaded$las_nig!=0])

table(dataTI3_free$cooccurence_myrspe)
table(dataTI3_invaded$cooccurence_myrspe)
length(dataTI3_invaded$myr_spe[dataTI3_invaded$myr_spe==0])
length(dataTI3_invaded$myr_spe[dataTI3_invaded$myr_spe!=0])

table(dataTI3_free$cooccurence_myrsab)
table(dataTI3_invaded$cooccurence_myrsab)
length(dataTI3_invaded$myr_sab[dataTI3_invaded$myr_sab==0])
length(dataTI3_invaded$myr_sab[dataTI3_invaded$myr_sab!=0])

table(dataTI3_free$cooccurence_tetsp)
table(dataTI3_invaded$cooccurence_tetsp)
length(dataTI3_invaded$tet_sp[dataTI3_invaded$tet_sp==0])
length(dataTI3_invaded$tet_sp[dataTI3_invaded$tet_sp!=0])


hist(dataTI3$tap_mag[dataTI3$tap_mag>0])

####

head(dataTI3presence)
dim(dataTI3presence)


dataTI3presence_free <- dataTI3presence[dataTI3presence$zone=="free", 18:31]
dataTI3presence_free_others <- rowSums(dataTI3presence_free[,-c(2,4,5,8)])
dataTI3presence_free_others[dataTI3presence_free_others>1] = 1
dataTI3presence_free_v2 <- cbind(dataTI3presence_free[,c(2,4,5,8)], dataTI3presence_free_others)
colnames(dataTI3presence_free_v2)[5] = "others"
dataTI3presence_free_v2 <- t(dataTI3presence_free_v2)

cooccur.ants_free <- cooccur(mat = dataTI3presence_free_v2, type = "spp_site", thresh = F, spp_names = TRUE)
prob.table(cooccur.ants_free)
plot(cooccur.ants_free)



dataTI3presence_invaded <- dataTI3presence[dataTI3presence$zone=="invaded", 17:31][,-c(2,7,8,15)]
dataTI3presence_invaded_others <- rowSums(dataTI3presence_invaded[,-c(1,2,4,5,6)])
dataTI3presence_invaded_others[dataTI3presence_invaded_others>1] = 1
dataTI3presence_invaded_v2 <- cbind(dataTI3presence_invaded[,c(1,2,4,5,6)], dataTI3presence_invaded_others)
colnames(dataTI3presence_invaded_v2)[6] = "others"
head(dataTI3presence_invaded_v2)
dataTI3presence_invaded_v2 <- t(dataTI3presence_invaded_v2)

cooccur.ants_invaded <- cooccur(mat = dataTI3presence_invaded_v2, type = "spp_site", thresh = F, spp_names = TRUE)
prob.table(cooccur.ants_invaded)
plot(cooccur.ants_invaded)


grid.arrange(pair.profile(cooccur.ants_invaded), obs.v.exp(cooccur.ants_invaded), ncol = 2)
