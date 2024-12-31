library(tidyverse)
library(plyr)
library(plotrix)
library(phytools)
library(geiger)
library(phangorn)
library(geomorph)
library(mvMORPH)

cols_family = c('black',#01Nandiniidae
                # "#ff4500", #02Prionodontidae
                "#FFD700", #03Felidae
                'red', #04Viverridae
                "#D2691E", #05Hyaenidae
                "green", #06Eupleridae
                'purple', #07Herpestidae
                'blue', #08Canidae
                'white', #09Ursidae
                #'orange', 
                #'yellow', 
                #'#00008b', 
                'orange', 
                'red', 
                'lightblue', 
                'lightgreen')

shapes_family <- c("01Nandiniidae" = 22,
                   #"02Prionodontidae" = 24,
                   "03Felidae" = 23,
                   "04Viverridae" = 21,
                   "05Hyaenidae" = 22,
                   "06Eupleridae" = 24,
                   "07Herpestidae" = 21,
                   "08Canidae" = 22,
                   "09Ursidae" = 24,
                   "13Mephitidae" = 21,
                   "14Ailuridae" = 22,
                   "22Procyonidae" = 24,
                   "16Mustelidae" = 21)

cols_preysize = c("carn_large" = 'red',
                  "carn_medium" = "orange", 
                  "carn_small" = "yellow", 
                  "herbivore" = 'green', 
                  "insectivore" = "black", 
                  "omnivore" = "white", 
                  "piscivore" = 'blue')

shapes_preysize = c("carn_large" = 22,
                    "carn_medium" = 24, 
                    "carn_small" = 23, 
                    "herbivore" = 21, 
                    "insectivore" = 22, 
                    "omnivore" = 24, 
                    "piscivore" = 23)


cols_locomotion3 <- c("arboreal" = "green",
                      "cursorial" = "purple",
                      "scansorial" = "yellow",
                      "semi_aquatic" = "lightblue",
                      "semi_fossorial" = "orange",
                      "terrestrial_hunt" = "black",
                      "terrestrial_nohunt" = "white")

shapes_locomotion3 <- c("arboreal" = 21,
                        "cursorial" = 22,
                        "scansorial" = 23,
                        "semi_aquatic" = 24,
                        "semi_fossorial" = 25,
                        "terrestrial_hunt" = 21,
                        "terrestrial_nohunt" = 22)

#### all_lsr adaptive landscapes ####
##### all - create PCA based on emperical data ##### 
#all_lsr PCA
all_lsr_pca <- gm.prcomp(all_lsr_data)
all_lsr_phylomorphospace <- gm.prcomp(all_lsr_data, phy = tree_prune)

#pca plot - locomotion3
plot(all_lsr_pca$x[,2] ~ all_lsr_pca$x[,1], pch = shapes_locomotion3[as.factor(locomotion3)], bg = cols_locomotion3[as.factor(locomotion3)], cex = log(data$geomean), phylo.par = list(tip.labels = FALSE, node.labels = FALSE), main = "all_lsr PCA", las = 1, cex.names = 3)

plot(all_lsr_phylomorphospace, phylo = TRUE, pch = shapes_locomotion3[as.factor(locomotion3)], bg = cols_locomotion3[as.factor(locomotion3)], cex = log(data$geomean), phylo.par = list(tip.labels = FALSE, node.labels = FALSE), main = "all_lsr PCA", las = 1, cex.names = 3)


all_pca <- gm.prcomp(all_data)
all_phylomorphospace <- gm.prcomp(all_data, phy = tree_prune)

#pca plot - family
plot(all_pca$x[,2] ~ all_pca$x[,1], pch = shapes_family[as.factor(family)], bg = cols_family[as.factor(family)], cex = log(data$geomean), phylo.par = list(tip.labels = FALSE, node.labels = FALSE), main = "all PCA", las = 1, cex.names = 3)

phylomorphospace(tree_prune, cbind(all_pca$x[,1], all_pca$x[,2]), label = "off", control = list(pch = shapes_family[as.factor(family)], bg = cols_family[as.factor(family)], cex = log(data$geomean)))


#pca plot - locomotion3
plot(all_pca$x[,2] ~ all_pca$x[,1], pch = shapes_locomotion3[as.factor(locomotion3)], bg = cols_locomotion3[as.factor(locomotion3)], cex = log(data$geomean), phylo.par = list(tip.labels = FALSE, node.labels = FALSE), main = "all PCA", las = 1, cex.names = 3)

#pca plot - preysize
plot(all_pca$x[,2] ~ all_pca$x[,1], pch = shapes_preysize[as.factor(preysize)], bg = cols_preysize[as.factor(preysize)], cex = log(data$geomean), phylo.par = list(tip.labels = FALSE, node.labels = FALSE), main = "all PCA", las = 1, cex.names = 3)



#save PC 1 and 2 scores
all_lsr_PC1_2 <- all_lsr_pca$x[,1:2]

#### ANOVAs ####
#locomotion
all_lsr_PC1_2_locomotion_rrpp <- lm.rrpp(all_lsr_PC1_2 ~ locomotion3, iter = 999, SS.type = "II")
anova(all_lsr_PC1_2_locomotion_rrpp)
all_lsr_PC1_2_locomotion_rrpp_pw <- pairwise(all_lsr_PC1_2_locomotion_rrpp, groups = locomotion3 )
summary(all_lsr_PC1_2_locomotion_rrpp_pw, show.vectors = TRUE, stat.table = FALSE)

#preysize
all_lsr_PC1_2_preysize_rrpp <- lm.rrpp(all_lsr_PC1_2 ~ preysize, iter = 999,  SS.type = "II")
anova(all_lsr_PC1_2_preysize_rrpp)
all_lsr_PC1_2_preysize_rrpp_pw <- pairwise(all_lsr_PC1_2_preysize_rrpp, groups = preysize )
summary(all_lsr_PC1_2_preysize_rrpp_pw, show.vectors = TRUE, stat.table = FALSE)


#### mvPGLS ####
#locomotion
all_lsr_PC1_2_locomotion_mvpanova <- mvgls(all_lsr_PC1_2 ~ locomotion3, tree = tree_prune, method = "LOOCV", model = "lambda")
all_lsr_PC1_2_locomotion_mvpanova_test <- manova.gls(all_lsr_PC1_2_locomotion_mvpanova, test = "Pillai", nbcores = 3, verbose = TRUE)
all_lsr_PC1_2_locomotion_mvpanova_test

#preysize
all_lsr_PC1_2_preysize_mvpanova <- mvgls(all_lsr_PC1_2 ~ preysize, tree = tree_prune, method = "LOOCV", model = "lambda")
all_lsr_PC1_2_preysize_mvpanova_test <- manova.gls(all_lsr_PC1_2_preysize_mvpanova, test = "Pillai", nbcores = 3, verbose = TRUE)
all_lsr_PC1_2_preysize_mvpanova
all_lsr_PC1_2_preysize_mvpanova_test


