library(Morphoscape)
library(geomorph)


##### load morphological data ##### 
#save PC 1 and 2 scores
all_lsr_PC1 <- all_lsr_pca$x[,1]
all_lsr_PC2 <- all_lsr_pca$x[,2]

all_dat <- data.frame(all_lsr_PC1, all_lsr_PC2, family, locomotion3, preysize)

#save eigenvectors 1 and 2
all_lsr_pca_loadings <- all_lsr_pca$rotation
all_lsr_PC1_eig <- all_lsr_pca_loadings[,1]
all_lsr_PC2_eig <- all_lsr_pca_loadings[,2]

#calculate mean data of all_lsr 
all_lsr_meandata <- colMeans(all_lsr_data)


##### create functional data ##### 
#calculate functional traits 
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

temMA <- scale_values(setNames((data$MAT/data$COL), data$species))
masMA <- scale_values(setNames((data$MAM2/data$COL), data$species))
SI <- scale_values(setNames(data$scap_L/data$scap_W, data$species))
BI <- scale_values(setNames(data$rad_L/data$hum_L, data$species))
HRI <- scale_values(setNames(data$hum_D/data$hum_L, data$species))
HEI <- scale_values(setNames(data$hum_dist/data$hum_L, data$species))
OLI <- scale_values(setNames(data$ul_OL/data$ul_L, data$species))
URI <- scale_values(setNames(data$ul_D/data$ul_L, data$species))
MANUS <- scale_values(setNames(data$MC3L/data$hum_L, data$species))
CI <- scale_values(setNames(data$tib_L/data$fem_L, data$species))
FRI <- scale_values(setNames(data$fem_D/data$fem_L, data$species))
GI <- scale_values(setNames(data$fem_D/data$fem_GT, data$species))
FEI <- scale_values(setNames(data$fem_EB/data$fem_L, data$species))
TRI <- scale_values(setNames(data$tib_D/data$tib_L, data$species))
PES <- scale_values(setNames(data$MT3L/data$fem_L, data$species))
C5_sagSMA <- scale_values(setNames((pi*data$C5_CW*data$C5_CH^3)/4, data$species)) 
C5_latSMA <- scale_values(setNames((pi*data$C5_CH*data$C5_CW^3)/4, data$species)) 
C5_SMA <- scale_values(rowMeans(cbind(C5_sagSMA, C5_latSMA)))
C5_JTA <- scale_values(setNames((360-abs(270-data$C5_PZA)), data$species))
C5_JV <- scale_values(setNames((abs(180-data$C5_PZA)), data$species))
midT_sagSMA <- scale_values(setNames((pi*data$midT_CW*data$midT_CH^3)/4, data$species)) 
midT_latSMA <- scale_values(setNames((pi*data$midT_CH*data$midT_CW^3)/4, data$species)) 
midT_SMA <- scale_values(rowMeans(cbind(midT_sagSMA, midT_latSMA)))
midT_JTA <- scale_values(setNames((360-abs(270-data$midT_PZA)), data$species) )
midT_JV <- scale_values(setNames((abs(180-data$midT_PZA)), data$species) )
tranT_sagSMA <- scale_values(setNames((pi*data$tranT_CW*data$tranT_CH^3)/4, data$species)) 
tranT_latSMA <- scale_values(setNames((pi*data$tranT_CH*data$tranT_CW^3)/4, data$species)) 
tranT_SMA <- scale_values(rowMeans(cbind(tranT_sagSMA, tranT_latSMA)))
tranT_JTA <- scale_values(setNames((360-abs(270-data$tranT_PZA)), data$species) )
tranT_JV <- scale_values(setNames((abs(180-data$tranT_PZA)), data$species))
midL_sagSMA <- scale_values(setNames((pi*data$midL_CW*data$midL_CH^3)/4, data$species)) 
midL_latSMA <- scale_values(setNames((pi*data$midL_CH*data$midL_CW^3)/4, data$species)) 
midL_SMA <- scale_values(rowMeans(cbind(midL_sagSMA, midL_latSMA)))
midL_JTA <- scale_values(setNames((360-abs(270-data$midL_PZA)), data$species) )
midL_JV <- scale_values(setNames((abs(180-data$midL_PZA)), data$species))


##### create performance surface based on functional data ##### 

#dataframe of PC positions of functional data 
all_wraps <- data.frame(cbind(cbind(all_lsr_PC1,all_lsr_PC2),
                              temMA,
                              masMA,
                              SI,
                              BI,
                              HRI,
                              HEI,
                              OLI,
                              MANUS,
                              CI,
                              FRI,
                              GI,
                              FEI,
                              TRI,
                              PES,
                              URI,
                              C5_SMA,
                              C5_JTA,
                              C5_JV,
                              midT_SMA,
                              midT_JTA,
                              midT_JV,
                              tranT_SMA,
                              tranT_JTA,
                              tranT_JV,
                              midL_SMA,
                              midL_JTA,
                              midL_JV))

#convert data frame to fnc_df object
all_fnc <- as_fnc_df(all_wraps, func.names = c(
  "temMA",
  "masMA",
  "SI",
  "BI",
  "HRI",
  "HEI",
  "OLI",
  "MANUS",
  "CI",
  "FRI",
  "GI",
  "FEI",
  "TRI",
  "PES",
  "URI",
  "C5_SMA",
  "C5_JTA",
  "C5_JV",
  "midT_SMA",
  "midT_JTA",
  "midT_JV",
  "tranT_SMA",
  "tranT_JTA",
  "tranT_JV",
  "midL_SMA",
  "midL_JTA",
  "midL_JV"))

# create grid from PC positions of functional data 
all_grid <- resample_grid(all_wraps, hull = NULL, padding = 1.1)

# Create kriged surface of functional data 
all_kr_surf <- krige_surf(all_fnc, grid = all_grid)
all_kr_surf
plot(all_kr_surf)


##### create performance surface based on functional data - remove bad traits ##### 

#dataframe of PC positions of functional data 
all_wraps <- data.frame(cbind(cbind(all_lsr_PC1,all_lsr_PC2),
                              temMA,
                              # masMA,
                              # SI,
                              BI,
                              HRI,
                              HEI,
                              OLI,
                              MANUS,
                              CI,
                              FRI,
                              # GI,
                              FEI,
                              #  TRI,
                              PES,
                              URI,
                              C5_SMA,
                              C5_JTA,
                              C5_JV,
                              midT_SMA,
                              midT_JTA,
                              midT_JV,
                              tranT_SMA,
                              tranT_JTA,
                              tranT_JV,
                              midL_SMA,
                              midL_JTA,
                              midL_JV))

#convert data frame to fnc_df object
all_fnc <- as_fnc_df(all_wraps, func.names = c(
  "temMA",
  # "masMA",
  #"SI",
  "BI",
  "HRI",
  "HEI",
  "OLI",
  "MANUS",
  "CI",
  "FRI",
  # "GI",
  "FEI",
  # "TRI",
  "PES",
  "URI",
  "C5_SMA",
  "C5_JTA",
  "C5_JV",
  "midT_SMA",
  "midT_JTA",
  "midT_JV",
  "tranT_SMA",
  "tranT_JTA",
  "tranT_JV",
  "midL_SMA",
  "midL_JTA",
  "midL_JV"))

# recreate grid from PC positions of functional data 
all_grid <- resample_grid(all_wraps, hull = NULL, padding = 1.1)

# recreate kriged surface of functional data 
all_kr_surf <- krige_surf(all_fnc, grid = all_grid)
all_kr_surf
plot(all_kr_surf)

# Create kriged surface of functional data (PC 1 and 2) 
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


##### create adaptive landscapes ##### 
# create matrix containing weight combinations 
all_weights <- generate_weights(n = 4 , data = all_kr_surf)

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


###### locomotion3 landscapes ##### 
#Calculate with optimally weighted adaptive landscapes by locomotion3
all_wprime_by_locomotion3 <- calcWprimeBy(all_all_landscapes, by = locomotion3)
all_wprime_by_locomotion3
summary(all_wprime_by_locomotion3)
write.csv(summary(all_wprime_by_locomotion3), "landscapes_locomotion3_empiricaldata.csv")

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


###### preysize landscapes ##### 
#Calculate with optimally weighted adaptive landscapes by preysize
all_wprime_by_preysize <- calcWprimeBy(all_all_landscapes, by = preysize)
all_wprime_by_preysize
summary(all_wprime_by_preysize)
write.csv(summary(all_wprime_by_preysize), "landscapes_preysize_empiricaldata.csv")

all_wprime_by_preysize_plot <- plot(all_wprime_by_preysize, countour = FALSE, ncol = 4) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 color= preysize,
                 shape = preysize)) +
  scale_shape_manual(values=shapes_preysize)+
  scale_color_manual(values = cols_preysize)+
  ggtitle("all preysize")
all_wprime_by_preysize_plot

#test for optimal adaptive landscape differences among preysize groups  
all_wprime_by_preysize_test <- multi.lands.grp.test(all_wprime_by_preysize)
all_wprime_by_preysize_test


###### family landscapes ##### 
#Calculate with optimally weighted adaptive landscapes by family
all_wprime_by_family <- calcWprimeBy(all_all_landscapes, by = family)
all_wprime_by_family
summary(all_wprime_by_family)
write.csv(summary(all_wprime_by_family), "landscapes_family_empiricaldata.csv")

all_wprime_by_family_plot <- plot(all_wprime_by_family, countour = FALSE, ncol=3) +
  geom_point(data = all_dat,
             aes(x = all_lsr_PC1, y = all_lsr_PC2,
                 color= family,
                 shape = family)) +
  scale_shape_manual(values=shapes_family)+
  scale_color_manual(values = cols_family)+
  ggtitle("all family")
all_wprime_by_family_plot

#test for optimal adaptive landscape differences among families   
all_wprime_by_family_test <- multi.lands.grp.test(all_wprime_by_family)
all_wprime_by_family_test



