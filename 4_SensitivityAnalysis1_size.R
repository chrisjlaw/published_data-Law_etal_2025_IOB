library(Morphoscape)
library(geomorph)


##### load morphological data ##### 
#save PC 1 and 2 scores
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

all_PC1 <- scale_values(all_pca$x[,1])
all_PC2 <- scale_values(all_pca$x[,2])

all_dat <- data.frame(all_PC1, all_PC2, family, locomotion3, preysize)

#####  create theoretical data ##### 
#save eigenvectors 1 and 2
all_pca_loadings <- all_pca$rotation
all_PC1_eig <- all_pca_loadings[,1]
all_PC2_eig <- all_pca_loadings[,2]

#load theoretical PC scores 
all_tPC <- read.csv("data/tdata_reduced/theoretical_PC_all_9x7_scaled.csv")

#calculate mean data of all 
all_meandata <- colMeans(all_data)

#create dataframe of theoretical trait/geomean data 
all_sr_tdata <- data.frame(rbind(
  exp(-14.28*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + 02.550*all_PC2_eig + all_meandata),
  
  exp(-14.28*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + 01.860*all_PC2_eig + all_meandata),
  
  exp(-14.28*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + 01.170*all_PC2_eig + all_meandata),
  
  exp(-14.28*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + 00.480*all_PC2_eig + all_meandata),
  
  exp(-14.28*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + -0.220*all_PC2_eig + all_meandata),
  
  exp(-14.28*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + -0.910*all_PC2_eig + all_meandata),
  
  exp(-14.28*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(-10.65*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(-07.03*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(-03.40*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(000.23*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(003.86*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(007.49*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(011.11*all_PC1_eig + -1.600*all_PC2_eig + all_meandata),
  exp(014.74*all_PC1_eig + -1.600*all_PC2_eig + all_meandata)))

#calculate theoretical functional traits 
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

temMA_theo <- scale_values(all_sr_tdata$lnMAT/all_sr_tdata$lnMOL)
masMA_theo <- scale_values(all_sr_tdata$lnMAM2/all_sr_tdata$lnMOL)
SI_theo <- scale_values(all_sr_tdata$lnscap_L/all_sr_tdata$lnscap_W)
BI_theo <- scale_values(all_sr_tdata$lnrad_L/all_sr_tdata$lnhum_L)
HRI_theo <- scale_values(all_sr_tdata$lnhum_D/all_sr_tdata$lnhum_L)
HEI_theo <- scale_values(all_sr_tdata$lnhum_dist/all_sr_tdata$lnhum_L)
OLI_theo <- scale_values(all_sr_tdata$lnul_OL/all_sr_tdata$lnul_L)
URI_theo <- scale_values(all_sr_tdata$lnul_D/all_sr_tdata$lnul_L)
MANUS_theo <- scale_values(all_sr_tdata$lnMC3L/all_sr_tdata$lnhum_L)
CI_theo <- scale_values(all_sr_tdata$lntib_L/all_sr_tdata$lnfem_L)
FRI_theo <- scale_values(all_sr_tdata$lnfem_D/all_sr_tdata$lnfem_L)
GI_theo <- scale_values(all_sr_tdata$lnfem_D/all_sr_tdata$lnfem_GT)
FEI_theo <- scale_values(all_sr_tdata$lnfem_EB/all_sr_tdata$lnfem_D)
TRI_theo <- scale_values(all_sr_tdata$lntib_D/all_sr_tdata$lntib_L)
PES_theo <- scale_values(all_sr_tdata$lnMT3L/all_sr_tdata$lnfem_L)
C5_sagSMA_theo <- scale_values((pi*all_sr_tdata$lnC5_CW*all_sr_tdata$lnC5_CH^3)/4) 
C5_latSMA_theo <- scale_values((pi*all_sr_tdata$lnC5_CH*all_sr_tdata$lnC5_CW^3)/4) 
C5_SMA_theo <- scale_values(rowMeans(cbind(C5_sagSMA_theo, C5_latSMA_theo)))
C5_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lnC5_PZA))
C5_JV_theo <- scale_values(abs(180-all_sr_tdata$lnC5_PZA))
midT_sagSMA_theo <- scale_values((pi*all_sr_tdata$lnmidT_CW*all_sr_tdata$lnmidT_CH^3)/4) 
midT_latSMA_theo <- scale_values((pi*all_sr_tdata$lnmidT_CH*all_sr_tdata$lnmidT_CW^3)/4) 
midT_SMA_theo <- scale_values(rowMeans(cbind(midT_sagSMA_theo, midT_latSMA_theo)))
midT_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lnmidT_PZA))
midT_JV_theo <- scale_values(abs(180-all_sr_tdata$lnmidT_PZA))
tranT_sagSMA_theo <- scale_values((pi*all_sr_tdata$lntranT_CW*all_sr_tdata$lntranT_CH^3)/4) 
tranT_latSMA_theo <- scale_values((pi*all_sr_tdata$lntranT_CH*all_sr_tdata$lntranT_CW^3)/4) 
tranT_SMA_theo <- scale_values(rowMeans(cbind(tranT_sagSMA_theo, tranT_latSMA_theo)))
tranT_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lntranT_PZA))
tranT_JV_theo <- scale_values(abs(180-all_sr_tdata$lntranT_PZA))
midL_sagSMA_theo <- scale_values((pi*all_sr_tdata$lnmidL_CW*all_sr_tdata$lnmidL_CH^3)/4) 
midL_latSMA_theo <- scale_values((pi*all_sr_tdata$lnmidL_CH*all_sr_tdata$lnmidL_CW^3)/4) 
midL_SMA_theo <- scale_values(rowMeans(cbind(midL_sagSMA_theo, midL_latSMA_theo)))
midL_JTA_theo <- scale_values(360-abs(270-all_sr_tdata$lnmidL_PZA))
midL_JV_theo <- scale_values(abs(180-all_sr_tdata$lnmidL_PZA))

#dataframe of PC positions of theoretical functional data 
all_wraps <- data.frame(cbind(all_tPC,
                              temMA_theo,
                              masMA_theo,
                              SI_theo,
                              BI_theo,
                              HRI_theo,
                              HEI_theo,
                              OLI_theo,
                              MANUS_theo,
                              CI_theo,
                              FRI_theo,
                              GI_theo,
                              FEI_theo,
                              TRI_theo,
                              PES_theo,
                              URI_theo,
                              C5_SMA_theo,
                              C5_JTA_theo,
                              C5_JV_theo,
                              midT_SMA_theo,
                              midT_JTA_theo,
                              midT_JV_theo,
                              tranT_SMA_theo,
                              tranT_JTA_theo,
                              tranT_JV_theo,
                              midL_SMA_theo,
                              midL_JTA_theo,
                              midL_JV_theo))



#convert data frame to fnc_df object
all_fnc <- as_fnc_df(all_wraps, func.names = c(
  "temMA_theo",
  "masMA_theo",
  "SI_theo",
  "BI_theo",
  "HRI_theo",
  "HEI_theo",
  "OLI_theo",
  "MANUS_theo",
  "CI_theo",
  "FRI_theo",
  "GI_theo",
  "FEI_theo",
  "TRI_theo",
  "PES_theo",
  "URI_theo",
  "C5_SMA_theo",
  "C5_JTA_theo",
  "C5_JV_theo",
  "midT_SMA_theo",
  "midT_JTA_theo",
  "midT_JV_theo",
  "tranT_SMA_theo",
  "tranT_JTA_theo",
  "tranT_JV_theo",
  "midL_SMA_theo",
  "midL_JTA_theo",
  "midL_JV_theo"))

##### all - create performance surface based on theoretical functional data ##### 
# create grid from PC positions of theoretical functional data 
all_grid <- resample_grid(all_wraps, hull = NULL, padding = 1.1)

# Create kriged surface of theoretical functional data 
all_kr_surf <- krige_surf(all_fnc, grid = all_grid)
all_kr_surf
plot(all_kr_surf)

# Create kriged surface of theoretical functional data with emperical trait data (PC 1 and 2) 
all_kr_surf <- krige_new_data(all_kr_surf, new_data = data.frame(cbind(all_PC1, all_PC2)))
all_kr_surf

all_kr_surf_plot <- plot(all_kr_surf, countour = FALSE) +
  geom_point(data = all_dat,
             aes(x = all_PC1, y = all_PC2,
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
             aes(x = all_PC1, y = all_PC2,
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
write.csv(summary(all_wprime_by_locomotion3), "senSize_landscapes_locomotion3.csv")

all_wprime_by_locomotion3_plot <- plot(all_wprime_by_locomotion3, countour = FALSE, ncol = 3) +
  geom_point(data = all_dat,
             aes(x = all_PC1, y = all_PC2,
                 color= locomotion3,
                 shape = locomotion3)) +
  scale_shape_manual(values=shapes_locomotion3)+
  scale_color_manual(values = cols_locomotion3)+
  ggtitle("all locomotion3")
all_wprime_by_locomotion3_plot

#test for optimal adaptive landscape differences among locomotion3 groups  
all_wprime_by_locomotion3_test <- multi.lands.grp.test(all_wprime_by_locomotion3, quantile = .05)
all_wprime_by_locomotion3_test


###### preysize landscapes ##### 
#Calculate with optimally weighted adaptive landscapes by preysize
all_wprime_by_preysize <- calcWprimeBy(all_all_landscapes, by = preysize)
all_wprime_by_preysize
summary(all_wprime_by_preysize)
write.csv(summary(all_wprime_by_preysize), "senSize_landscapes_preysize.csv")

all_wprime_by_preysize_plot <- plot(all_wprime_by_preysize, countour = FALSE, ncol = 3) +
  geom_point(data = all_dat,
             aes(x = all_PC1, y = all_PC2,
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
write.csv(summary(all_wprime_by_family), "senSize_landscapes_family.csv")

all_wprime_by_family_plot <- plot(all_wprime_by_family, countour = FALSE, ncol=3) +
  geom_point(data = all_dat,
             aes(x = all_PC1, y = all_PC2,
                 color= family,
                 shape = family)) +
  scale_shape_manual(values=shapes_family)+
  scale_color_manual(values = cols_family)+
  ggtitle("all family")
all_wprime_by_family_plot

#test for optimal adaptive landscape differences among families   
all_wprime_by_family_test <- multi.lands.grp.test(all_wprime_by_family)
all_wprime_by_family_test

