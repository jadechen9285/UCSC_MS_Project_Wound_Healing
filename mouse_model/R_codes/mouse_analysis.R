# II.
# mouse_wound Affy pre-processed data analysis
# skin vs tongue in 8 time point 
# author: Jianhong Chen
# date: 08/27/2019

#######################################################################
# intial set up
ls() # list all objects
gc() # force R to release memory of deleted variables
rm(list = ls()) # delete all variables in the environment
# ---------------------------------------------
library('ggplot2')
library('tidyverse') # inlcude "tidyr", "readr", "dplyr", "purrr"
library('reshape2') 
library('plotly') # 3D plotting
library("kernlab") # kernel-based machine learning techniques 
library("factoextra") # multivariable analysis in R [hkmean]
library("dendextend") # create dendrogram object for hclustering plot 
library("cluster") # clustering functions
library("fpc") # Density-Based Spatial Clustering and application with Noise
library("dbscan")# Density-Based Spatial Clustering and application with Noise

library(RColorBrewer)

# load the data set
dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/cleaned_data"
setwd(dir)
mw_df = read_csv("mw_DE_fc1_v1.csv")
# change time to factor value first
time_lab = c("0hr", "6hrs", "12hrs", "24hrs", "3days", "5days", "7days", "10days")
mw_df$time = factor(mw_df$time,levels = c(0,6, 12, 24,72,120,168,240))

mw_wide = spread(mw_df, key = time, value = intensity)
colnames(mw_wide)[-(1:4)] = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")

# compute the difference in intensity bewteen each wounding time vs unwound time points
mw_diff = matrix(0, dim(mw_wide)[1], length(time_lab)-1) # create the 0 matrix with the right size
for (i in 1:(length(time_lab)-1)){
  mw_temp = mw_wide[, -(1:4)] # create dummy stroage df
  mw_diff[, i] = unlist(abs(mw_temp[, i+1] - mw_temp[, 1]))
  rm(mw_temp) # remove dummy df
}
colnames(mw_diff) = c("d1", "d2", "d3", "d4", "d5","d6", "d7")
mw_diff = cbind(mw_wide[,1:4], data.frame(mw_diff))

# Skin tissue with all samples
# perform cluastering using the difference in intensity as the similarity measurement
# NOTICE: using the original intensity data set without any additional computation
skin_diff = filter(mw_diff, tissue =="skin" )
skin_original = filter(mw_wide, tissue == "skin")
skin_cl = kmeans(skin_diff[, -(1:4)], centers = 4)
skin_PCA = prcomp(skin_diff[, -(1:4)], scale = T) # prcomp uses svd for PCA
#calculate the eigenalues percentage
PCA_eigen_per = percentVar = round(100*skin_PCA$sdev^2 / sum(skin_PCA$sdev^2), 1)
#
skin_GG = cbind(skin_diff[,c(1:2,4)], data.frame(skin_PCA$x)) %>%
  mutate("cluster" = factor(skin_cl$cluster))
#
ggplot(skin_GG, aes(PC1, PC2)) +
  geom_point(alpha = 0.6,aes(colour = cluster)) +
  ggtitle("Visualization in PCA axis for clustering of Skin(All Samples)") +
  xlab(paste0("PC1, VarExp: ", PCA_eigen_per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", PCA_eigen_per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = brewer.pal(4, "Set2"))

# Tongue tissue with all samples
# perform cluastering using the difference in intensity as the similarity measurement
# NOTICE: using the original intensity data set without any additional computation
tongue_diff = filter(mw_diff, tissue == "tongue")
tongue_cl = kmeans(tongue_diff[, -(1:4)], centers = 4)
tongue_PCA = prcomp(tongue_diff[, -(1:4)], scale = T) 
PCA_eigen_per = percentVar = round(100*tongue_PCA$sdev^2 / sum(tongue_PCA$sdev^2), 1)
#
tongue_GG = cbind(tongue_diff[,c(1:2,4)], data.frame(tongue_PCA$x)) %>%
  mutate("cluster" = factor(tongue_cl$cluster))
#
ggplot(tongue_GG, aes(PC1, PC2)) +
  geom_point(alpha = 0.6,aes(colour = cluster)) +
  ggtitle("Visualization in PCA axis for clustering of Tongue(All Samples)") +
  xlab(paste0("PC1, VarExp: ", PCA_eigen_per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", PCA_eigen_per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = brewer.pal(4, "Set2"))

# -----------------------------
# after clustering: analysis the grouped geens and try to annote the function of the these 
#       genes to check if they belong to one wound healing stage

# skin sample1:
skin_gene_g1 = filter(skin_s1_GG, cluster ==1)[, 1:2]














