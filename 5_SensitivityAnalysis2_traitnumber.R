library(Morphoscape)
library(geomorph)
library(MASS)


#### limb LDA ####
#lda code from Slater, G. J. 2022. Topographically distinct adaptive landscapes for teeth, skeletons, and size explain the adaptive radiation of Carnivora (Mammalia). Evolution, doi: 10.1111/evo.14577.

limb_biomech <- cbind(fore_biomech, hind_biomech)
limb_biomech_lda <- lda(locomotion3 ~ limb_biomech)


(limb_biomech_lda$svd^2)/sum(limb_biomech_lda$svd^2)

limb_biomech_lda_df1 <- apply(limb_biomech, 1, function(x) sum(x*limb_biomech_lda$scaling[,1]))
limb_biomech_lda_df2 <- apply(limb_biomech, 1, function(x) sum(x*limb_biomech_lda$scaling[,2]))

cbind(limb_biomech_lda$scaling[,1], limb_biomech_lda$scaling[,2])
plot(limb_biomech_lda_df1, limb_biomech_lda_df2, pch = shapes_locomotion3[as.factor(locomotion3)], bg = cols_locomotion3[as.factor(locomotion3)], cex=1.5)

limb_pred <- MASS:::predict.lda(limb_biomech_lda, as.data.frame(limb_biomech))$class

limb_biomech_ldaCV <- lda(locomotion3 ~ limb_biomech, CV=TRUE)

limb_biomech_lda_table <- (cbind(locomotion3, as.character(limb_pred), round(limb_biomech_ldaCV$posterior,2)))

limb_biomech_lda_conftable <- table(locomotion3, limb_biomech_ldaCV$class, dnn = c('Actual Group','Predicted Group'))


#### vert LDA ####
#lda code from Slater, G. J. 2022. Topographically distinct adaptive landscapes for teeth, skeletons, and size explain the adaptive radiation of Carnivora (Mammalia). Evolution, doi: 10.1111/evo.14577.

vert_biomech <- cbind(C5_biomech,
                      midT_biomech,
                      tranT_biomech,
                      midL_biomech)
vert_biomech_lda <- lda(locomotion3 ~ vert_biomech)

(vert_biomech_lda$svd^2)/sum(vert_biomech_lda$svd^2)

vert_biomech_lda_df1 <- apply(vert_biomech, 1, function(x) sum(x*vert_biomech_lda$scaling[,1]))
vert_biomech_lda_df2 <- apply(vert_biomech, 1, function(x) sum(x*vert_biomech_lda$scaling[,2]))

cbind(vert_biomech_lda$scaling[,1], vert_biomech_lda$scaling[,2])
plot(vert_biomech_lda_df1, vert_biomech_lda_df2, pch = shapes_locomotion3[as.factor(locomotion3)], bg = cols_locomotion3[as.factor(locomotion3)], cex=1.5)

vert_pred <- MASS:::predict.lda(vert_biomech_lda, as.data.frame(vert_biomech))$class

vert_biomech_ldaCV <- lda(locomotion3 ~ vert_biomech, CV=TRUE)

vert_biomech_lda_table <- (cbind(locomotion3, as.character(vert_pred), round(vert_biomech_ldaCV$posterior,2)))

vert_biomech_lda_conftable <- table(locomotion3, vert_biomech_ldaCV$class, dnn = c('Actual Group','Predicted Group'))


#### load morphological data ####
#save PC 1 and 2 scores
all_lsr_PC1 <- all_lsr_pca$x[,1]
all_lsr_PC2 <- all_lsr_pca$x[,2]

all_dat <- data.frame(all_lsr_PC1, all_lsr_PC2, family, locomotion3, preysize)

#### create theoretical data ####
#save eigenvectors 1 and 2
all_lsr_pca_loadings <- all_lsr_pca$rotation
all_lsr_PC1_eig <- all_lsr_pca_loadings[,1]
all_lsr_PC2_eig <- all_lsr_pca_loadings[,2]

#load theoretical PC scores 
all_tPC <- read.csv("data/tdata_reduced/theoretical_PC_all_lsr_9x7.csv")

#calculate mean data of all_lsr 
all_lsr_meandata <- colMeans(all_lsr_data)

#### #9x7 create dataframe of theoretical trait/geomean data ####
all_sr_tdata <- data.frame(rbind(
  exp(-3.150*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + 01.560*all_lsr_PC2_eig + all_lsr_meandata),
  
  exp(-3.150*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + 00.930*all_lsr_PC2_eig + all_lsr_meandata),
  
  exp(-3.150*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + 00.310*all_lsr_PC2_eig + all_lsr_meandata),
  
  exp(-3.150*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + -0.320*all_lsr_PC2_eig + all_lsr_meandata),
  
  exp(-3.150*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + -0.950*all_lsr_PC2_eig + all_lsr_meandata),
  
  exp(-3.150*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + -1.570*all_lsr_PC2_eig + all_lsr_meandata),
  
  exp(-3.150*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-2.568*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.987*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-1.406*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.825*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(-0.244*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.338*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(00.919*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata),
  exp(01.500*all_lsr_PC1_eig + -2.200*all_lsr_PC2_eig + all_lsr_meandata)))


#### 14 traits ####

##### theoretical functional proxies #####
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

temMA_theo <- scale_values(all_sr_tdata$lnMAT_lsr/all_sr_tdata$lnMOL_lsr)
masMA_theo <- scale_values(all_sr_tdata$lnMAM2_lsr/all_sr_tdata$lnMOL_lsr)
HRI_theo <- scale_values(all_sr_tdata$lnhum_D_lsr/all_sr_tdata$lnhum_L_lsr)
HEI_theo <- scale_values(all_sr_tdata$lnhum_dist_lsr/all_sr_tdata$lnhum_L_lsr)
OLI_theo <- scale_values(all_sr_tdata$lnul_OL_lsr/all_sr_tdata$lnul_L_lsr)
FRI_theo <- scale_values(all_sr_tdata$lnfem_D_lsr/all_sr_tdata$lnfem_L_lsr)
FEI_theo <- scale_values(all_sr_tdata$lnfem_EB_lsr/all_sr_tdata$lnfem_D_lsr)
TRI_theo <- scale_values(all_sr_tdata$lntib_D_lsr/all_sr_tdata$lntib_L_lsr)
C5_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lnC5_PZA))
C5_JV_theo <- scale_values(abs(180-all_sr_tdata$lnC5_PZA))
midT_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lnmidT_PZA))
midT_JV_theo <- scale_values(abs(180-all_sr_tdata$lnmidT_PZA))
tranT_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lntranT_PZA))
tranT_JV_theo <- scale_values(abs(180-all_sr_tdata$lntranT_PZA))

#dataframe of PC positions of theoretical functional data 
all_wraps <- data.frame(cbind(all_tPC,
                              temMA_theo,
                              masMA_theo,
                              HRI_theo,
                              HEI_theo,
                              OLI_theo,
                              FRI_theo,
                              FEI_theo,
                              TRI_theo,
                              C5_JTA_theo,
                              C5_JV_theo,
                              midT_JTA_theo,
                              midT_JV_theo,
                              tranT_JTA_theo,
                              tranT_JV_theo
))

#convert data frame to fnc_df object
all_fnc <- as_fnc_df(all_wraps, func.names = c(
  "temMA_theo",
  "masMA_theo",
  "HRI_theo",
  "HEI_theo",
  "OLI_theo",
  "FRI_theo",
  "FEI_theo",
  "TRI_theo",
  "C5_JTA_theo",
  "C5_JV_theo",
  "midT_JTA_theo",
  "midT_JV_theo",
  "tranT_JTA_theo",
  "tranT_JV_theo"
))

###### create performance surface based on theoretical functional data ###### 
# create grid from PC positions of theoretical functional data 
all_grid <- resample_grid(all_wraps, hull = NULL, padding = 1.1)

# Create kriged surface of theoretical functional data 
all_kr_surf <- krige_surf(all_fnc, grid = all_grid)
all_kr_surf
plot(all_kr_surf)

# Create kriged surface of theoretical functional data with emperical trait data (PC 1 and 2) 
all_kr_surf <- krige_new_data(all_kr_surf, new_data = data.frame(cbind(all_lsr_PC1, all_lsr_PC2)))
all_kr_surf

all_kr_surf_plot <- plot(all_kr_surf, countour = FALSE) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 shape= locomotion3,
                 color = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)
all_kr_surf_plot

######create adaptive landscapes ###### 
# create matrix containing weight combinations 
all_weights <- generate_weights(step = .15, data = all_kr_surf)

# Create adaptive landscapes 
all_all_landscapes <- calc_all_lscps(all_kr_surf, grid_weights = all_weights, file = "all_all_landscapes_reduced.RData")

#Calculate with optimally weighted adaptive landscapes 
all_wprime_all <- calcGrpWprime(all_all_landscapes)
all_wprime_all

all_wprime_all_plot <- plot(all_wprime_all, countour = FALSE) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 shape= locomotion3,
                 color = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)
all_wprime_all_plot

#Calculate with optimally weighted adaptive landscapes by locomotion3
all_wprime_by_locomotion3 <- calcWprimeBy(all_all_landscapes, by = locomotion3)
all_wprime_by_locomotion3
summary(all_wprime_by_locomotion3)
write.csv(summary(all_wprime_by_locomotion3), "senProxies14_landscapes_locomotion3.csv")

all_wprime_by_locomotion3_plot <- plot(all_wprime_by_locomotion3, countour = FALSE, ncol = 3) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 color= locomotion3,
                 shape = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)+
  ggtitle("all locomotion3")
all_wprime_by_locomotion3_plot

#test for optimal adaptive landscape differences among locomotion3 groups  
all_wprime_by_locomotion3_test <- multi.lands.grp.test(all_wprime_by_locomotion3)
all_wprime_by_locomotion3_test



#### 6 traits ####

##### calculate theoretical functional traits #####
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

HRI_theo <- scale_values(all_sr_tdata$lnhum_D_lsr/all_sr_tdata$lnhum_L_lsr)
HEI_theo <- scale_values(all_sr_tdata$lnhum_dist_lsr/all_sr_tdata$lnhum_L_lsr)
FRI_theo <- scale_values(all_sr_tdata$lnfem_D_lsr/all_sr_tdata$lnfem_L_lsr)
C5_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lnC5_PZA))
C5_JV_theo <- scale_values(abs(180-all_sr_tdata$lnC5_PZA))
midT_JV_theo <- scale_values(abs(180-all_sr_tdata$lnmidT_PZA))


###### create performance surface based on theoretical functional data ###### 

#dataframe of PC positions of theoretical functional data 
all_wraps <- data.frame(cbind(all_tPC,
                              HRI_theo,
                              HEI_theo,
                              FRI_theo,
                              C5_JTA_theo,
                              C5_JV_theo,
                              midT_JV_theo
))



#convert data frame to fnc_df object
all_fnc <- as_fnc_df(all_wraps, func.names = c(
  "HRI_theo",
  "HEI_theo",
  "FRI_theo",
  "C5_JTA_theo",
  "C5_JV_theo",
  "midT_JV_theo"
))

# create grid from PC positions of theoretical functional data 
all_grid <- resample_grid(all_wraps, hull = NULL, padding = 1.1)

# Create kriged surface of theoretical functional data 
all_kr_surf <- krige_surf(all_fnc, grid = all_grid)
all_kr_surf
plot(all_kr_surf)

# Create kriged surface of theoretical functional data with emperical trait data (PC 1 and 2) 
all_kr_surf <- krige_new_data(all_kr_surf, new_data = data.frame(cbind(all_lsr_PC1, all_lsr_PC2)))
all_kr_surf

all_kr_surf_plot <- plot(all_kr_surf, countour = FALSE) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 shape= locomotion3,
                 color = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)
all_kr_surf_plot

###### create adaptive landscapes ###### 
# create matrix containing weight combinations 
all_weights <- generate_weights(step = .05, data = all_kr_surf)

# Create adaptive landscapes 
all_all_landscapes <- calc_all_lscps(all_kr_surf, grid_weights = all_weights, file = "all_all_landscapes_reduced.RData")

#Calculate with optimally weighted adaptive landscapes 
all_wprime_all <- calcGrpWprime(all_all_landscapes)
all_wprime_all

all_wprime_all_plot <- plot(all_wprime_all, countour = FALSE) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 shape= locomotion3,
                 color = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)
all_wprime_all_plot

#Calculate with optimally weighted adaptive landscapes by locomotion3
all_wprime_by_locomotion3 <- calcWprimeBy(all_all_landscapes, by = locomotion3)
all_wprime_by_locomotion3
summary(all_wprime_by_locomotion3)
write.csv(summary(all_wprime_by_locomotion3), "senProxies6_landscapes_locomotion3.csv")

all_wprime_by_locomotion3_plot <- plot(all_wprime_by_locomotion3, countour = FALSE, ncol = 3) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 color= locomotion3,
                 shape = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)+
  ggtitle("all locomotion3")
all_wprime_by_locomotion3_plot

#test for optimal adaptive landscape differences among locomotion3 groups  
all_wprime_by_locomotion3_test <- multi.lands.grp.test(all_wprime_by_locomotion3)
all_wprime_by_locomotion3_test

