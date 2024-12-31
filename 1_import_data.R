library(tidyverse)
library(plyr)
library(plotrix)
library(phytools)
library(geiger)
library(phangorn)

##species means##
rawdata <- read.csv("rawdata.csv")
family <- rawdata$family
species <- rawdata$species

species_avg <- ddply(rawdata, .(species, family), colwise(mean, na.rm = TRUE))
rownames(species_avg) <- species_avg[,1]
data <- species_avg


multitree <- read.nexus("Upham_1k_carnivorans.nex")
multitree <- lapply(multitree,method = "extend", force.ultrametric)
multitree <- lapply(multitree, ladderize)
multitree <- lapply(multitree, multi2di)

multitree_prune <- list(length=length(multitree))
for(i in 1:length(multitree)){
  multitree_prune[[i]] <- treedata(phy = multitree[[i]], data = data, warnings=TRUE)$phy
}
class(multitree_prune) <- "multiPhylo"


tree <- maxCladeCred(multitree)
#write.tree(tree, file = "data/Upham_mcctree_carnivorans.tre")
tree_prune <- treedata(tree, species_avg, sort=TRUE, warnings = TRUE)$phy #check that all species are in tree (nothing should appear)
data <- species_avg[tree_prune$tip.label,] #Sorts the rows of data so that they match the order of the tip labels. SUPER IMPORTANT


#### Ecology data ####
#load ecology data
eco_data <- read.csv(file = "data_ecology.csv")
rownames(eco_data) <- eco_data[,1] #make column 1 rownames
eco_data_prune <- data.frame(treedata(phy = tree_prune, data = eco_data, sort=FALSE)$data)

#match and combine pruned Procrustes and ecological data
data$preysize <- eco_data_prune$prey_size[match(data$species,eco_data_prune$species)]
data$locomotion3 <- eco_data_prune$locomotion3[match(data$species,eco_data_prune$species)]

preysize <- setNames(data$preysize, data$species)
locomotion3 <- setNames(data$locomotion3, data$species)

family <- setNames(data$family, data$species)


data %>% group_by(preysize) %>% tally()
data %>% group_by(locomotion) %>% tally()


geomean <- setNames(log(data$geomean), data$species)


#### cranial traits #### 
lnCBL <- setNames(log(data$CBL), data$species)
lnPOC <- setNames(log(data$POC), data$species)
lnZB <- setNames(log(data$ZB), data$species)
lnMB <- setNames(log(data$MB), data$species)
lnAFH <- setNames(log(data$AFH), data$species)
lnPFH <- setNames(log(data$PFH), data$species)
lnCVH <- setNames(log(data$CVH), data$species)
lnPW <- setNames(log(data$PW), data$species)
lnBCL <- setNames(log(data$BCL), data$species)

lnCBL_lsr <- setNames(log(data$CBL/data$geomean), data$species)
lnPOC_lsr <- setNames(log(data$POC/data$geomean), data$species)
lnZB_lsr <- setNames(log(data$ZB/data$geomean), data$species)
lnMB_lsr <- setNames(log(data$MB/data$geomean), data$species)
lnCVH_lsr <- setNames(log(data$CVH/data$geomean), data$species)
lnPW_lsr <- setNames(log(data$PW/data$geomean), data$species)
lnBCL_lsr <- setNames(log(data$BCL/data$geomean), data$species)

cran_data <- cbind(lnCBL,
                   lnPOC,
                   lnZB,
                   lnMB,
                   lnCVH,
                   lnPW,
                   lnBCL)

cran_lsr_data <- cbind(lnCBL_lsr,
                       lnPOC_lsr,
                       lnZB_lsr,
                       lnMB_lsr,
                       lnCVH_lsr,
                       lnPW_lsr,
                       lnBCL_lsr)

#### mandible traits #### 
lnMAT <- setNames(log(data$MAT), data$species)
lnMAM <- setNames(log(data$MAM), data$species)
lnMAM2 <- setNames(log(data$MAM2), data$species)
lnMOL <- setNames(log(data$MOL), data$species)
lnCOL <- setNames(log(data$COL), data$species)
lnMW <- setNames(log(data$MW), data$species)
lnAMW <- setNames(log(data$AMW), data$species)

lnMAT_lsr <- setNames(log(data$MAT/data$geomean), data$species)
lnMAM_lsr <- setNames(log(data$MAM/data$geomean), data$species)
lnMAM2_lsr <- setNames(log(data$MAM2/data$geomean), data$species)
lnMOL_lsr <- setNames(log(data$MOL/data$geomean), data$species)
lnCOL_lsr <- setNames(log(data$COL/data$geomean), data$species)
lnMW_lsr <- setNames(log(data$MW)/data$geomean, data$species)
lnAMW_lsr <- setNames(log(data$AMW/data$geomean), data$species)

mand_data <- cbind(lnMAT,
                   lnMAM, 
                   lnMAM2,
                   lnMOL, 
                   lnCOL, 
                   lnMW,
                   lnAMW)

mand_lsr_data <- cbind(lnMAT_lsr,
                       lnMAM_lsr, 
                       lnMAM2_lsr,
                       lnMOL_lsr, 
                       lnCOL_lsr, 
                       lnMW_lsr,
                       lnAMW_lsr)


#### C3 vertebrae #### 
lnC3_VW <- setNames(log(data$C3_VW), data$species)
lnC3_VH <- setNames(log(data$C3_VH), data$species)
lnC3_CL <- setNames(log(data$C3_CL), data$species)
lnC3_CW <- setNames(log(data$C3_CW), data$species)
lnC3_CH <- setNames(log(data$C3_CH), data$species)
lnC3_ZL <- setNames(log(data$C3_ZL), data$species)
lnC3_NSH <- setNames(log(data$C3_NSH), data$species)
lnC3_TPDV <- setNames(log(data$C3_TPDV), data$species)
lnC3_PZA <- setNames(log(data$C3_PZA), data$species)
lnC3_NSA <- setNames(log(data$C3_NSA), data$species)
lnC3_TPAP <- setNames(log(data$C3_TPAP), data$species)

lnC3_VW_lsr <- setNames(log(data$C3_VW/data$geomean), data$species)
lnC3_VH_lsr <- setNames(log(data$C3_VH/data$geomean), data$species)
lnC3_CL_lsr <- setNames(log(data$C3_CL/data$geomean), data$species)
lnC3_CW_lsr <- setNames(log(data$C3_CW/data$geomean), data$species)
lnC3_CH_lsr <- setNames(log(data$C3_CH/data$geomean), data$species)
lnC3_ZL_lsr <- setNames(log(data$C3_ZL/data$geomean), data$species)
lnC3_NSH_lsr <- setNames(log(data$C3_NSH/data$geomean), data$species)
lnC3_TPDV_lsr <- setNames(log(data$C3_TPDV), data$species)
lnC3_PZA_lsr <- setNames(log(data$C3_PZA), data$species)
lnC3_NSA_lsr <- setNames(log(data$C3_NSA), data$species)
lnC3_TPAP_lsr <- setNames(log(data$C3_TPAP), data$species)

C3_data <- cbind(
  lnC3_VW,
  lnC3_VH,
  lnC3_CL,
  lnC3_CW,
  lnC3_CH,
  lnC3_ZL,
  lnC3_NSH,
  lnC3_TPDV,
  lnC3_PZA,
  lnC3_NSA,
  lnC3_TPAP)

C3_lsr_data <- cbind(
  lnC3_VW_lsr,
  lnC3_VH_lsr,
  lnC3_CL_lsr,
  lnC3_CW_lsr,
  lnC3_CH_lsr,
  lnC3_ZL_lsr,
  lnC3_NSH_lsr,
  lnC3_TPDV_lsr,
  lnC3_PZA_lsr,
  lnC3_NSA_lsr,
  lnC3_TPAP_lsr)

#### C5 vertebrae #### 
lnC5_VW <- setNames(log(data$C5_VW), data$species)
lnC5_VH <- setNames(log(data$C5_VH), data$species)
lnC5_CL <- setNames(log(data$C5_CL), data$species)
lnC5_CW <- setNames(log(data$C5_CW), data$species)
lnC5_CH <- setNames(log(data$C5_CH), data$species)
lnC5_ZL <- setNames(log(data$C5_ZL), data$species)
lnC5_NSH <- setNames(log(data$C5_NSH), data$species)
lnC5_TPDV <- setNames(log(data$C5_TPDV), data$species)
lnC5_PZA <- setNames(log(data$C5_PZA), data$species)
lnC5_NSA <- setNames(log(data$C5_NSA), data$species)
lnC5_TPAP <- setNames(log(data$C5_TPAP), data$species)

lnC5_VW_lsr <- setNames(log(data$C5_VW/data$geomean), data$species)
lnC5_VH_lsr <- setNames(log(data$C5_VH/data$geomean), data$species)
lnC5_CL_lsr <- setNames(log(data$C5_CL/data$geomean), data$species)
lnC5_CW_lsr <- setNames(log(data$C5_CW/data$geomean), data$species)
lnC5_CH_lsr <- setNames(log(data$C5_CH/data$geomean), data$species)
lnC5_ZL_lsr <- setNames(log(data$C5_ZL/data$geomean), data$species)
lnC5_NSH_lsr <- setNames(log(data$C5_NSH/data$geomean), data$species)
lnC5_TPDV_lsr <- setNames(log(data$C5_TPDV), data$species)
lnC5_PZA_lsr <- setNames(log(data$C5_PZA), data$species)
lnC5_NSA_lsr <- setNames(log(data$C5_NSA), data$species)
lnC5_TPAP_lsr <- setNames(log(data$C5_TPAP), data$species)

C5_data <- cbind(
  lnC5_VW,
  lnC5_VH,
  lnC5_CL,
  lnC5_CW,
  lnC5_CH,
  lnC5_ZL,
  lnC5_NSH,
  lnC5_TPDV,
  lnC5_PZA,
  lnC5_NSA,
  lnC5_TPAP)

C5_lsr_data <- cbind(
  lnC5_VW_lsr,
  lnC5_VH_lsr,
  lnC5_CL_lsr,
  lnC5_CW_lsr,
  lnC5_CH_lsr,
  lnC5_ZL_lsr,
  lnC5_NSH_lsr,
  lnC5_TPDV_lsr,
  lnC5_PZA_lsr,
  lnC5_NSA_lsr,
  lnC5_TPAP_lsr)

#### T1 vertebrae #### 
lnT1_VW <- setNames(log(data$T1_VW), data$species)
lnT1_VH <- setNames(log(data$T1_VH), data$species)
lnT1_CL <- setNames(log(data$T1_CL), data$species)
lnT1_CW <- setNames(log(data$T1_CW), data$species)
lnT1_CH <- setNames(log(data$T1_CH), data$species)
lnT1_ZL <- setNames(log(data$T1_ZL), data$species)
lnT1_NSH <- setNames(log(data$T1_NSH), data$species)
lnT1_TPDV <- setNames(log(data$T1_TPDV), data$species)
lnT1_PZA <- setNames(log(data$T1_PZA), data$species)
lnT1_NSA <- setNames(log(data$T1_NSA), data$species)
lnT1_TPAP <- setNames(log(data$T1_TPAP), data$species)

lnT1_VW_lsr <- setNames(log(data$T1_VW/data$geomean), data$species)
lnT1_VH_lsr <- setNames(log(data$T1_VH/data$geomean), data$species)
lnT1_CL_lsr <- setNames(log(data$T1_CL/data$geomean), data$species)
lnT1_CW_lsr <- setNames(log(data$T1_CW/data$geomean), data$species)
lnT1_CH_lsr <- setNames(log(data$T1_CH/data$geomean), data$species)
lnT1_ZL_lsr <- setNames(log(data$T1_ZL/data$geomean), data$species)
lnT1_NSH_lsr <- setNames(log(data$T1_NSH/data$geomean), data$species)
lnT1_TPDV_lsr <- setNames(log(data$T1_TPDV), data$species)
lnT1_PZA_lsr <- setNames(log(data$T1_PZA), data$species)
lnT1_NSA_lsr <- setNames(log(data$T1_NSA), data$species)
lnT1_TPAP_lsr <- setNames(log(data$T1_TPAP), data$species)

T1_data <- cbind(
  lnT1_VW,
  lnT1_VH,
  lnT1_CL,
  lnT1_CW,
  lnT1_CH,
  lnT1_ZL,
  lnT1_NSH,
  lnT1_TPDV,
  lnT1_PZA,
  lnT1_NSA,
  lnT1_TPAP)

T1_lsr_data <- cbind(
  lnT1_VW_lsr,
  lnT1_VH_lsr,
  lnT1_CL_lsr,
  lnT1_CW_lsr,
  lnT1_CH_lsr,
  lnT1_ZL_lsr,
  lnT1_NSH_lsr,
  lnT1_TPDV_lsr,
  lnT1_PZA_lsr,
  lnT1_NSA_lsr,
  lnT1_TPAP_lsr)


#### midT vertebrae #### 
lnmidT_VW <- setNames(log(data$midT_VW), data$species)
lnmidT_VH <- setNames(log(data$midT_VH), data$species)
lnmidT_CL <- setNames(log(data$midT_CL), data$species)
lnmidT_CW <- setNames(log(data$midT_CW), data$species)
lnmidT_CH <- setNames(log(data$midT_CH), data$species)
lnmidT_ZL <- setNames(log(data$midT_ZL), data$species)
lnmidT_NSH <- setNames(log(data$midT_NSH), data$species)
lnmidT_TPDV <- setNames(log(data$midT_TPDV), data$species)
lnmidT_PZA <- setNames(log(data$midT_PZA), data$species)
lnmidT_NSA <- setNames(log(data$midT_NSA), data$species)
lnmidT_TPAP <- setNames(log(data$midT_TPAP), data$species)

lnmidT_VW_lsr <- setNames(log(data$midT_VW/data$geomean), data$species)
lnmidT_VH_lsr <- setNames(log(data$midT_VH/data$geomean), data$species)
lnmidT_CL_lsr <- setNames(log(data$midT_CL/data$geomean), data$species)
lnmidT_CW_lsr <- setNames(log(data$midT_CW/data$geomean), data$species)
lnmidT_CH_lsr <- setNames(log(data$midT_CH/data$geomean), data$species)
lnmidT_ZL_lsr <- setNames(log(data$midT_ZL/data$geomean), data$species)
lnmidT_NSH_lsr <- setNames(log(data$midT_NSH/data$geomean), data$species)
lnmidT_TPDV_lsr <- setNames(log(data$midT_TPDV), data$species)
lnmidT_PZA_lsr <- setNames(log(data$midT_PZA), data$species)
lnmidT_NSA_lsr <- setNames(log(data$midT_NSA), data$species)
lnmidT_TPAP_lsr <- setNames(log(data$midT_TPAP), data$species)

midT_data <- cbind(
  lnmidT_VW,
  lnmidT_VH,
  lnmidT_CL,
  lnmidT_CW,
  lnmidT_CH,
  lnmidT_ZL,
  lnmidT_NSH,
  lnmidT_TPDV,
  lnmidT_PZA,
  lnmidT_NSA,
  lnmidT_TPAP)

midT_lsr_data <- cbind(
  lnmidT_VW_lsr,
  lnmidT_VH_lsr,
  lnmidT_CL_lsr,
  lnmidT_CW_lsr,
  lnmidT_CH_lsr,
  lnmidT_ZL_lsr,
  lnmidT_NSH_lsr,
  lnmidT_TPDV_lsr,
  lnmidT_PZA_lsr,
  lnmidT_NSA_lsr,
  lnmidT_TPAP_lsr)


#### tranT vertebrae #### 
lntranT_VW <- setNames(log(data$tranT_VW), data$species)
lntranT_VH <- setNames(log(data$tranT_VH), data$species)
lntranT_CL <- setNames(log(data$tranT_CL), data$species)
lntranT_CW <- setNames(log(data$tranT_CW), data$species)
lntranT_CH <- setNames(log(data$tranT_CH), data$species)
lntranT_ZL <- setNames(log(data$tranT_ZL), data$species)
lntranT_NSH <- setNames(log(data$tranT_NSH), data$species)
lntranT_PZA <- setNames(log(data$tranT_PZA), data$species)
lntranT_NSA <- setNames(log(data$tranT_NSA), data$species)
lntranT_TPDV <- setNames(log(data$tranT_TPDV), data$species)

lntranT_VW_lsr <- setNames(log(data$tranT_VW/data$geomean), data$species)
lntranT_VH_lsr <- setNames(log(data$tranT_VH/data$geomean), data$species)
lntranT_CL_lsr <- setNames(log(data$tranT_CL/data$geomean), data$species)
lntranT_CW_lsr <- setNames(log(data$tranT_CW/data$geomean), data$species)
lntranT_CH_lsr <- setNames(log(data$tranT_CH/data$geomean), data$species)
lntranT_ZL_lsr <- setNames(log(data$tranT_ZL/data$geomean), data$species)
lntranT_NSH_lsr <- setNames(log(data$tranT_NSH/data$geomean), data$species)
lntranT_TPDV_lsr <- setNames(log(data$tranT_TPDV), data$species)
lntranT_PZA_lsr <- setNames(log(data$tranT_PZA), data$species)
lntranT_NSA_lsr <- setNames(log(data$tranT_NSA), data$species)

tranT_data <- cbind(
  lntranT_VW,
  lntranT_VH,
  lntranT_CL,
  lntranT_CW,
  lntranT_CH,
  lntranT_ZL,
  lntranT_NSH,
  lntranT_PZA,
  lntranT_NSA)

tranT_lsr_data <- cbind(
  lntranT_VW_lsr,
  lntranT_VH_lsr,
  lntranT_CL_lsr,
  lntranT_CW_lsr,
  lntranT_CH_lsr,
  lntranT_ZL_lsr,
  lntranT_NSH_lsr,
  lntranT_TPDV_lsr,
  lntranT_PZA_lsr,
  lntranT_NSA_lsr)

#### lastT vertebrae #### 
lnlastT_VW <- setNames(log(data$lastT_VW), data$species)
lnlastT_VH <- setNames(log(data$lastT_VH), data$species)
lnlastT_CL <- setNames(log(data$lastT_CL), data$species)
lnlastT_CW <- setNames(log(data$lastT_CW), data$species)
lnlastT_CH <- setNames(log(data$lastT_CH), data$species)
lnlastT_ZL <- setNames(log(data$lastT_ZL), data$species)
lnlastT_NSH <- setNames(log(data$lastT_NSH), data$species)
lnlastT_PZA <- setNames(log(data$lastT_PZA), data$species)
lnlastT_NSA <- setNames(log(data$lastT_NSA), data$species)

lnlastT_VW_lsr <- setNames(log(data$lastT_VW/data$geomean), data$species)
lnlastT_VH_lsr <- setNames(log(data$lastT_VH/data$geomean), data$species)
lnlastT_CL_lsr <- setNames(log(data$lastT_CL/data$geomean), data$species)
lnlastT_CW_lsr <- setNames(log(data$lastT_CW/data$geomean), data$species)
lnlastT_CH_lsr <- setNames(log(data$lastT_CH/data$geomean), data$species)
lnlastT_ZL_lsr <- setNames(log(data$lastT_ZL/data$geomean), data$species)
lnlastT_NSH_lsr <- setNames(log(data$lastT_NSH/data$geomean), data$species)
lnlastT_PZA_lsr <- setNames(log(data$lastT_PZA), data$species)
lnlastT_NSA_lsr <- setNames(log(data$lastT_NSA), data$species)

lastT_data <- cbind(
  lnlastT_VW,
  lnlastT_VH,
  lnlastT_CL,
  lnlastT_CW,
  lnlastT_CH,
  lnlastT_ZL,
  lnlastT_NSH,
  lnlastT_PZA,
  lnlastT_NSA)

lastT_lsr_data <- cbind(
  lnlastT_VW_lsr,
  lnlastT_VH_lsr,
  lnlastT_CL_lsr,
  lnlastT_CW_lsr,
  lnlastT_CH_lsr,
  lnlastT_ZL_lsr,
  lnlastT_NSH_lsr,
  lnlastT_PZA_lsr,
  lnlastT_NSA_lsr)

#### L1 vertebrae #### 
lnL1_VW <- setNames(log(data$L1_VW), data$species)
lnL1_VH <- setNames(log(data$L1_VH), data$species)
lnL1_CL <- setNames(log(data$L1_CL), data$species)
lnL1_CW <- setNames(log(data$L1_CW), data$species)
lnL1_CH <- setNames(log(data$L1_CH), data$species)
lnL1_ZL <- setNames(log(data$L1_ZL), data$species)
lnL1_NSH <- setNames(log(data$L1_NSH), data$species)
lnL1_TPDV <- setNames(log(data$L1_TPDV), data$species)
lnL1_PZA <- setNames(log(data$L1_PZA), data$species)
lnL1_NSA <- setNames(log(data$L1_NSA), data$species)
lnL1_TPAP <- setNames(log(data$L1_TPAP), data$species)

lnL1_VW_lsr <- setNames(log(data$L1_VW/data$geomean), data$species)
lnL1_VH_lsr <- setNames(log(data$L1_VH/data$geomean), data$species)
lnL1_CL_lsr <- setNames(log(data$L1_CL/data$geomean), data$species)
lnL1_CW_lsr <- setNames(log(data$L1_CW/data$geomean), data$species)
lnL1_CH_lsr <- setNames(log(data$L1_CH/data$geomean), data$species)
lnL1_ZL_lsr <- setNames(log(data$L1_ZL/data$geomean), data$species)
lnL1_NSH_lsr <- setNames(log(data$L1_NSH/data$geomean), data$species)
lnL1_TPDV_lsr <- setNames(log(data$L1_TPDV), data$species)
lnL1_PZA_lsr <- setNames(log(data$L1_PZA), data$species)
lnL1_NSA_lsr <- setNames(log(data$L1_NSA), data$species)
lnL1_TPAP_lsr <- setNames(log(data$L1_TPAP), data$species)


L1_data <- cbind(
  lnL1_VW,
  lnL1_VH,
  lnL1_CL,
  lnL1_CW,
  lnL1_CH,
  lnL1_ZL,
  lnL1_NSH,
  lnL1_TPDV,
  lnL1_PZA,
  lnL1_NSA,
  lnL1_TPAP)

L1_lsr_data <- cbind(
  lnL1_VW_lsr,
  lnL1_VH_lsr,
  lnL1_CL_lsr,
  lnL1_CW_lsr,
  lnL1_CH_lsr,
  lnL1_ZL_lsr,
  lnL1_NSH_lsr,
  lnL1_TPDV_lsr,
  lnL1_PZA_lsr,
  lnL1_NSA_lsr,
  lnL1_TPAP_lsr)

#### midL vertebrae #### 
lnmidL_VW <- setNames(log(data$midL_VW), data$species)
lnmidL_VH <- setNames(log(data$midL_VH), data$species)
lnmidL_CL <- setNames(log(data$midL_CL), data$species)
lnmidL_CW <- setNames(log(data$midL_CW), data$species)
lnmidL_CH <- setNames(log(data$midL_CH), data$species)
lnmidL_ZL <- setNames(log(data$midL_ZL), data$species)
lnmidL_NSH <- setNames(log(data$midL_NSH), data$species)
lnmidL_TPDV <- setNames(log(data$midL_TPDV), data$species)
lnmidL_PZA <- setNames(log(data$midL_PZA), data$species)
lnmidL_NSA <- setNames(log(data$midL_NSA), data$species)
lnmidL_TPAP <- setNames(log(data$midL_TPAP), data$species)

lnmidL_VW_lsr <- setNames(log(data$midL_VW/data$geomean), data$species)
lnmidL_VH_lsr <- setNames(log(data$midL_VH/data$geomean), data$species)
lnmidL_CL_lsr <- setNames(log(data$midL_CL/data$geomean), data$species)
lnmidL_CW_lsr <- setNames(log(data$midL_CW/data$geomean), data$species)
lnmidL_CH_lsr <- setNames(log(data$midL_CH/data$geomean), data$species)
lnmidL_ZL_lsr <- setNames(log(data$midL_ZL/data$geomean), data$species)
lnmidL_NSH_lsr <- setNames(log(data$midL_NSH/data$geomean), data$species)
lnmidL_TPDV_lsr <- setNames(log(data$midL_TPDV), data$species)
lnmidL_PZA_lsr <- setNames(log(data$midL_PZA), data$species)
lnmidL_NSA_lsr <- setNames(log(data$midL_NSA), data$species)
lnmidL_TPAP_lsr <- setNames(log(data$midL_TPAP), data$species)

midL_data <- cbind(
  lnmidL_VW,
  lnmidL_VH,
  lnmidL_CL,
  lnmidL_CW,
  lnmidL_CH,
  lnmidL_ZL,
  lnmidL_NSH,
  lnmidL_TPDV,
  lnmidL_PZA,
  lnmidL_NSA,
  lnmidL_TPAP)

midL_lsr_data <- cbind(
  lnmidL_VW_lsr,
  lnmidL_VH_lsr,
  lnmidL_CL_lsr,
  lnmidL_CW_lsr,
  lnmidL_CH_lsr,
  lnmidL_ZL_lsr,
  lnmidL_NSH_lsr,
  lnmidL_TPDV_lsr,
  lnmidL_PZA_lsr,
  lnmidL_NSA_lsr,
  lnmidL_TPAP_lsr)

#### lastL vertebrae #### 
lnlastL_VW <- setNames(log(data$lastL_VW), data$species)
lnlastL_VH <- setNames(log(data$lastL_VH), data$species)
lnlastL_CL <- setNames(log(data$lastL_CL), data$species)
lnlastL_CW <- setNames(log(data$lastL_CW), data$species)
lnlastL_CH <- setNames(log(data$lastL_CH), data$species)
lnlastL_ZL <- setNames(log(data$lastL_ZL), data$species)
lnlastL_NSH <- setNames(log(data$lastL_NSH), data$species)
lnlastL_TPDV <- setNames(log(data$lastL_TPDV), data$species)
lnlastL_PZA <- setNames(log(data$lastL_PZA), data$species)
lnlastL_NSA <- setNames(log(data$lastL_NSA), data$species)
lnlastL_TPAP <- setNames(log(data$lastL_TPAP), data$species)

lnlastL_VW_lsr <- setNames(log(data$lastL_VW/data$geomean), data$species)
lnlastL_VH_lsr <- setNames(log(data$lastL_VH/data$geomean), data$species)
lnlastL_CL_lsr <- setNames(log(data$lastL_CL/data$geomean), data$species)
lnlastL_CW_lsr <- setNames(log(data$lastL_CW/data$geomean), data$species)
lnlastL_CH_lsr <- setNames(log(data$lastL_CH/data$geomean), data$species)
lnlastL_ZL_lsr <- setNames(log(data$lastL_ZL/data$geomean), data$species)
lnlastL_NSH_lsr <- setNames(log(data$lastL_NSH/data$geomean), data$species)
lnlastL_TPDV_lsr <- setNames(log(data$lastL_TPDV), data$species)
lnlastL_PZA_lsr <- setNames(log(data$lastL_PZA), data$species)
lnlastL_NSA_lsr <- setNames(log(data$lastL_NSA), data$species)
lnlastL_TPAP_lsr <- setNames(log(data$lastL_TPAP), data$species)


lastL_data <- cbind(
  lnlastL_VW,
  lnlastL_VH,
  lnlastL_CL,
  lnlastL_CW,
  lnlastL_CH,
  lnlastL_ZL,
  lnlastL_NSH,
  lnlastL_TPDV,
  lnlastL_PZA,
  lnlastL_NSA,
  lnlastL_TPAP)

lastL_lsr_data <- cbind(
  lnlastL_VW_lsr,
  lnlastL_VH_lsr,
  lnlastL_CL_lsr,
  lnlastL_CW_lsr,
  lnlastL_CH_lsr,
  lnlastL_ZL_lsr,
  lnlastL_NSH_lsr,
  lnlastL_TPDV_lsr,
  lnlastL_PZA_lsr,
  lnlastL_NSA_lsr,
  lnlastL_TPAP_lsr)

#### sacral traits ####
lnsacral_L <- setNames(log(data$sacral_L), data$species)
lnsacral_H <- setNames(log(data$sacral_H), data$species)
lnsacral_W <- setNames(log(data$sacral_W), data$species)

lnsacral_L_lsr <- setNames(log(data$sacral_L/data$geomean), data$species)
lnsacral_H_lsr <- setNames(log(data$sacral_H/data$geomean), data$species)
lnsacral_W_lsr <- setNames(log(data$sacral_W/data$geomean), data$species)

sacral_lsr_data <- cbind(
  lnsacral_L_lsr,
  lnsacral_H_lsr,
  lnsacral_W_lsr
)



#### forelimb traits ####
lnscap_L <- setNames(log(data$scap_L), data$species)
lnscap_W <- setNames(log(data$scap_W), data$species)
lnhum_L <- setNames(log(data$hum_L), data$species)
lnhum_D <- setNames(log(data$hum_D), data$species)
lnhum_prox <- setNames(log(data$hum_prox), data$species)
lnhum_dist <- setNames(log(data$hum_dist), data$species)
lnul_L <- setNames(log(data$ul_L), data$species)
lnul_D <- setNames(log(data$ul_D), data$species)
lnul_OL <- setNames(log(data$ul_OL), data$species)
lnrad_L <- setNames(log(data$rad_L), data$species)
lnrad_D <- setNames(log(data$rad_D), data$species)
lnMC3L <- setNames(log(data$MC3L), data$species)
lnMC3W <- setNames(log(data$MC3W), data$species)

lnscap_L_lsr <- setNames(log(data$scap_L/data$geomean), data$species)
lnscap_W_lsr <- setNames(log(data$scap_W/data$geomean), data$species)
lnhum_L_lsr <- setNames(log(data$hum_L/data$geomean), data$species)
lnhum_D_lsr <- setNames(log(data$hum_D/data$geomean), data$species)
lnhum_prox_lsr <- setNames(log(data$hum_prox/data$geomean), data$species)
lnhum_dist_lsr <- setNames(log(data$hum_dist/data$geomean), data$species)
lnul_L_lsr <- setNames(log(data$ul_L/data$geomean), data$species)
lnul_D_lsr <- setNames(log(data$ul_D/data$geomean), data$species)
lnul_OL_lsr <- setNames(log(data$ul_OL/data$geomean), data$species)
lnrad_L_lsr <- setNames(log(data$rad_L/data$geomean), data$species)
lnrad_D_lsr <- setNames(log(data$rad_D/data$geomean), data$species)
lnMC3L_lsr <- setNames(log(data$MC3L/data$geomean), data$species)
lnMC3W_lsr <- setNames(log(data$MC3W/data$geomean), data$species)


fore_data <- cbind(
  lnscap_L,
  lnscap_W,
  lnhum_L,
  lnhum_D,
  lnhum_prox,
  lnhum_dist,
  lnul_L, 
  lnul_D,
  lnul_OL,
  lnrad_L,
  lnrad_D,
  lnMC3L, 
  lnMC3W 
)

fore_lsr_data <- cbind(
  lnscap_L_lsr,
  lnscap_W_lsr,
  lnhum_L_lsr,
  lnhum_D_lsr,
  lnhum_prox_lsr,
  lnhum_dist_lsr,
  lnul_L_lsr, 
  lnul_D_lsr,
  lnul_OL_lsr,
  lnrad_L_lsr,
  lnrad_D_lsr,
  lnMC3L_lsr, 
  lnMC3W_lsr 
)


#### hindlimb traits ####
lnpel_L <- setNames(log(data$pel_L), data$species)
lnfem_L <- setNames(log(data$fem_L), data$species)
lnfem_D <- setNames(log(data$fem_D), data$species)
lnfem_EB <- setNames(log(data$fem_EB), data$species)
lnfem_GT <- setNames(log(data$fem_GT), data$species)
lntib_L <- setNames(log(data$tib_L), data$species)
lntib_D <- setNames(log(data$tib_D), data$species)
lntib_prox <- setNames(log(data$tib_prox), data$species)
lntib_dist <- setNames(log(data$tib_dist), data$species)
lnfib_L <- setNames(log(data$fib_L), data$species)
lncal_L <- setNames(log(data$cal_L), data$species)
lnMT3L <- setNames(log(data$MT3L), data$species)
lnMT3D <- setNames(log(data$MT3D), data$species)

lnpel_L_lsr <- setNames(log(data$pel_L/data$geomean), data$species)
lnfem_L_lsr <- setNames(log(data$fem_L/data$geomean), data$species)
lnfem_D_lsr <- setNames(log(data$fem_D/data$geomean), data$species)
lnfem_EB_lsr <- setNames(log(data$fem_EB/data$geomean), data$species)
lnfem_GT_lsr <- setNames(log(data$fem_GT/data$geomean), data$species)
lntib_L_lsr <- setNames(log(data$tib_L/data$geomean), data$species)
lntib_D_lsr <- setNames(log(data$tib_D/data$geomean), data$species)
lntib_prox_lsr <- setNames(log(data$tib_prox/data$geomean), data$species)
lntib_dist_lsr <- setNames(log(data$tib_dist/data$geomean), data$species)
lnfib_L_lsr <- setNames(log(data$fib_L/data$geomean), data$species)
lncal_L_lsr <- setNames(log(data$cal_L/data$geomean), data$species)
lnMT3L_lsr <- setNames(log(data$MT3L/data$geomean), data$species)
lnMT3D_lsr <- setNames(log(data$MT3D/data$geomean), data$species)


hind_data <- cbind(
  lnpel_L,
  lnfem_L ,
  lnfem_D ,
  lnfem_EB,
  lnfem_GT,
  lntib_L ,
  lntib_D ,
  lntib_prox,
  lntib_dist,
  lnfib_L,
  lncal_L,
  lnMT3L,  
  lnMT3D 
)
hind_lsr_data <- cbind(
  lnpel_L_lsr,
  lnfem_L_lsr ,
  lnfem_D_lsr ,
  lnfem_EB_lsr,
  lnfem_GT_lsr,
  lntib_L_lsr ,
  lntib_D_lsr ,
  lntib_prox_lsr,
  lntib_dist_lsr,
  lnfib_L_lsr,
  lncal_L_lsr,
  lnMT3L_lsr,  
  lnMT3D_lsr 
)


all_data <- cbind(cran_data,
                  mand_data,
                  fore_data,
                  hind_data,
                  C3_data,
                  C5_data,
                  T1_data,
                  midT_data,
                  tranT_data,
                  lastT_data,
                  L1_data,
                  midL_data,
                  lastL_data)


all_lsr_data <- cbind(cran_lsr_data,
                      mand_lsr_data,
                      fore_lsr_data,
                      hind_lsr_data,
                      C3_lsr_data,
                      C5_lsr_data,
                      T1_lsr_data,
                      midT_lsr_data,
                      tranT_lsr_data,
                      lastT_lsr_data,
                      L1_lsr_data,
                      midL_lsr_data,
                      lastL_lsr_data)



#### Calculate Functional Proxies ####




temMA_can <- setNames((data$MAT/data$COL), data$species)
temMA_mol <- setNames((data$MAT/data$MOL), data$species)
masMA_can <- setNames((data$MAM2/data$COL), data$species)
masMA_mol <- setNames((data$MAM2/data$MOL), data$species)

mand_biomech <- cbind(temMA_can,
                      temMA_mol,
                      masMA_can,
                      masMA_mol)

C3_sagSMA <- setNames((pi*data$C3_CW*data$C3_CH^3)/4, data$species) #second moment of area = sagittal stiffness
C3_latSMA <- setNames((pi*data$C3_CH*data$C3_CW^3)/4, data$species) #second moment of area = latittal stiffness
C3_JTA <- setNames((360-abs(270-data$C3_PZA)), data$species) 
C3_JV <- setNames((abs(180-data$C3_PZA)), data$species) 

C3_biomech <- cbind(C3_sagSMA,
                    C3_latSMA, 
                    C3_JTA,
                    C3_JV)

C5_sagSMA <- setNames((pi*data$C5_CW*data$C5_CH^3)/4, data$species) #second moment of area = sagittal stiffness
C5_latSMA <- setNames((pi*data$C5_CH*data$C5_CW^3)/4, data$species) #second moment of area = latittal stiffness
C5_JTA <- setNames((360-abs(270-data$C5_PZA)), data$species) 
C5_JV <- setNames((abs(180-data$C5_PZA)), data$species) 

C5_biomech <- cbind(C5_sagSMA,
                    C5_latSMA, 
                    C5_JTA,
                    C5_JV)

T1_sagSMA <- setNames((pi*data$T1_CW*data$T1_CH^3)/4, data$species) #second moment of area = sagittal stiffness
T1_latSMA <- setNames((pi*data$T1_CH*data$T1_CW^3)/4, data$species) #second moment of area = latittal stiffness
T1_JTA <- setNames((360-abs(270-data$T1_PZA)), data$species) 
T1_JV <- setNames((abs(180-data$T1_PZA)), data$species) 

T1_biomech <- cbind(T1_sagSMA,
                    T1_latSMA, 
                    T1_JTA,
                    T1_JV)

midT_sagSMA <- setNames((pi*data$midT_CW*data$midT_CH^3)/4, data$species) #second moment of area = sagittal stiffness
midT_latSMA <- setNames((pi*data$midT_CH*data$midT_CW^3)/4, data$species) #second moment of area = latittal stiffness
midT_JTA <- setNames((360-abs(270-data$midT_PZA)), data$species) 
midT_JV <- setNames((abs(180-data$midT_PZA)), data$species) 

midT_biomech <- cbind(midT_sagSMA,
                      midT_latSMA, 
                      midT_JTA,
                      midT_JV)

tranT_sagSMA <- setNames((pi*data$tranT_CW*data$tranT_CH^3)/4, data$species) #second moment of area = sagittal stiffness
tranT_latSMA <- setNames((pi*data$tranT_CH*data$tranT_CW^3)/4, data$species) #second moment of area = latittal stiffness
tranT_JTA <- setNames((360-abs(270-data$tranT_PZA)), data$species) 
tranT_JV <- setNames((abs(180-data$tranT_PZA)), data$species) 

tranT_biomech <- cbind(tranT_sagSMA,
                       tranT_latSMA, 
                       tranT_JTA,
                       tranT_JV)

lastT_sagSMA <- setNames((pi*data$lastT_CW*data$lastT_CH^3)/4, data$species) #second moment of area = sagittal stiffness
lastT_latSMA <- setNames((pi*data$lastT_CH*data$lastT_CW^3)/4, data$species) #second moment of area = latittal stiffness
lastT_JTA <- setNames((360-abs(270-data$lastT_PZA)), data$species) 
lastT_JV <- setNames((abs(180-data$lastT_PZA)), data$species) 

lastT_biomech <- cbind(lastT_sagSMA,
                       lastT_latSMA, 
                       lastT_JTA,
                       lastT_JV)

L1_sagSMA <- setNames((pi*data$L1_CW*data$L1_CH^3)/4, data$species) #second moment of area = sagittal stiffness
L1_latSMA <- setNames((pi*data$L1_CH*data$L1_CW^3)/4, data$species) #second moment of area = latittal stiffness
L1_JTA <- setNames((360-abs(270-data$L1_PZA)), data$species) 
L1_JV <- setNames((abs(180-data$L1_PZA)), data$species) 

L1_biomech <- cbind(L1_sagSMA,
                    L1_latSMA, 
                    L1_JTA,
                    L1_JV)

midL_sagSMA <- setNames((pi*data$midL_CW*data$midL_CH^3)/4, data$species) #second moment of area = sagittal stiffness
midL_latSMA <- setNames((pi*data$midL_CH*data$midL_CW^3)/4, data$species) #second moment of area = latittal stiffness
midL_JTA <- setNames((360-abs(270-data$midL_PZA)), data$species) 
midL_JV <- setNames((abs(180-data$midL_PZA)), data$species) 

midL_biomech <- cbind(midL_sagSMA,
                      midL_latSMA, 
                      midL_JTA,
                      midL_JV)

lastL_sagSMA <- setNames((pi*data$lastL_CW*data$lastL_CH^3)/4, data$species) #second moment of area = sagittal stiffness
lastL_latSMA <- setNames((pi*data$lastL_CH*data$lastL_CW^3)/4, data$species) #second moment of area = latittal stiffness
lastL_JTA <- setNames((360-abs(270-data$lastL_PZA)), data$species) 
lastL_JV <- setNames((abs(180-data$lastL_PZA)), data$species) 

lastL_biomech <- cbind(lastL_sagSMA,
                       lastL_latSMA, 
                       lastL_JTA,
                       lastL_JV)

SI <- setNames(data$scap_L/data$scap_W, data$species)
BI <- setNames(data$rad_L/data$hum_L, data$species)
HRI <- setNames(data$hum_D/data$hum_L, data$species)
HEI <- setNames(data$hum_dist/data$hum_L, data$species)
OLI <- setNames(data$ul_OL/data$ul_L, data$species)
URI <- setNames(data$ul_D/data$ul_L, data$species)
MANUS <- setNames(data$MC3L/data$hum_L, data$species)

fore_biomech <- cbind(SI,
                      BI,
                      HRI,
                      HEI,
                      OLI,
                      URI,
                      MANUS)

CI <- setNames(data$tib_L/data$fem_L, data$species)
FRI <- setNames(data$fem_D/data$fem_L, data$species)
GI <- setNames(data$fem_D/data$fem_GT, data$species)
FEI <- setNames(data$fem_EB/data$fem_L, data$species)
TRI <- setNames(data$tib_D/data$tib_L, data$species)
PES <- setNames(data$MT3L/data$fem_L, data$species)
IM <- setNames((data$hum_L+data$rad_L)/(data$fem_L+data$tib_L), data$species)

hind_biomech <- cbind(CI,
                      FRI,
                      GI,
                      FEI,
                      TRI,
                      PES,
                      IM)


data_biomech <- cbind(mand_biomech,
                      C3_biomech,
                      C5_biomech,
                      T1_biomech,
                      midT_biomech,
                      tranT_biomech,
                      lastT_biomech,
                      L1_biomech,
                      midL_biomech,
                      lastL_biomech,
                      fore_biomech,
                      hind_biomech)



