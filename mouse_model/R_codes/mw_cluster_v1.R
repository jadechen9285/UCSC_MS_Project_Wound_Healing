# Author: Jianhong Chen
# date: August 7 2019
# perform quick kmeans cluster and PCA on the mouse wound data set

# load the necessary libraries 
library(tidyverse)
library(cluster)
library(kernlab)
# for plotting 
library(RColorBrewer)
library(gplots)
library(ggplot2)


# -----------------------------------------------------------------------------------------------------
# intial set up
ls() # list all objects 
gc() # force R to release memory of deleted variables 
rm(list = ls()) # delete all variables in the environment
# -----------------------------------------------------------------------------------------------------

# load the data set into dataframe 
dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/cleaned_data"
setwd(dir) 
group_df = read_csv("groupV2.csv")
mw_df = read_csv("mw_final.csv")
skin_df = filter(mw_df, tissue == "skin")
time_lab = unique(skin_df$time)
skin_df_t1 = filter(mw_df, tissue == "skin" & time == time_lab[1])

# skin_df_t1 = filter(group_df, tissue == "skin" )
# 
# skin_df_t1_cl = select(skin_df_t1, intensity) %>%
#   kmeans(centers = 3)
# 
# skin_df_t1 = mutate(skin_df_t1, "cluster" = skin_df_t1_cl$cluster)
# 
# ggplot(skin_df_t1, aes(x = factor(group), y = intensity, color = factor(cluster))) +
#   geom_jitter()

# find the best k numer
# optimization k-value of the clustering algorithm 
total_withinSS_km = map_dbl(1:20, function(k){
  model_km = kmeans(skin_df_t1$intensity, centers = k)
  model_km$tot.withinss
})

km_elbow = data.frame( 'k' = 1:20, 'tot_withinSS' = total_withinSS_km)
# the plot implies k = 2 as the optimized value
ggplot(km_elbow, aes(x = k, y = tot_withinSS)) +
  ggtitle("elbow plot for Spearman Correlation Distance") +
  geom_line() +
  scale_x_continuous(breaks = 1:20) + 
  theme(text = element_text(size = 25))

time_list = list()
for (t in 1:8) {
  time_list[[t]] = filter(mw_df, time == time_lab[[t]])
}

mw_wide = time_list[[1]] %>%
  select(c(symbol, gene_des, tissue, sample))
for (t in 1:8){
  mw_wide = cbind(mw_wide, time_list[[t]]$intensity)
}
colnames(mw_wide)[5:12] = time_lab

skin_wide = filter(mw_wide, tissue == "skin")

skin_wide_cl = skin_wide[, 5:12] %>%
  kmeans(centers = 5)
skin_wide = mutate(skin_wide, "cluster" = skin_wide_cl$cluster)

skin_long = skin_wide %>%
  gather(time, intensity, 5:12) 
skin_long$time = factor(skin_long$time, levels = rev(time_lab))

ggplot(skin_long, aes(x = cluster, y = intensity, group = cluster )) + 
  geom_boxplot() + 
  facet_wrap(~time) + 
  ggtitle("Clustering Results over time for Skin Tissue") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))
  

# Prinple component analysis of the dataset 
#exp_raw = log2(Biobase::exprs(raw_data))
skin_wide_list = tongue_wide_list = list()
for (s in 1:3){
  skin_wide_list[[s]] = filter(mw_wide, tissue == "skin" & sample == s)
  tongue_wide_list[[s]] = filter(mw_wide, tissue == "tongue" & sample == s)
}


PCA_data = skin_wide_list[[1]] %>%
  select(c(symbol, gene_des))
for (s in 1:3){
  PCA_data = cbind(PCA_data, skin_wide_list[[s]][,5:12], tongue_wide_list[[s]][, 5:12])
}

tissue_PCA = rep(c(rep("skin", 8), rep("tongue", 8)),3)
sample_PCA = factor(c(rep(1,16), rep(2, 16), rep(3, 16)), levels = c(1,2,3))
time_PCA = factor(rep(c("10 days", "7 days",  "5 days", "3 days",  
                    "24 hrs", "12 hrs",  "6 hrs",   "0 hr"), 6), 
              levels = c("10 days", "7 days",  "5 days", "3 days",  
                         "24 hrs", "12 hrs",  "6 hrs",   "0 hr"))

PCA_raw = prcomp(t(PCA_data[, 3:50]), scale. = F) #note: t() = transpose 
# PCA, R use SVD to compute PCs, output singular value in "sdev"
percentVar = round(100*PCA_raw$sdev^2 / sum(PCA_raw$sdev^2), 1) #calculate the eigenalues percentage
sd_ratio = sqrt(percentVar[2]/percentVar[1]) # ratio between the top eigenvalues scores PCs

dataGG = data.frame(PC1 = PCA_raw$x[, 1], PC2 = PCA_raw$x[, 2],
                    Time = time_PCA, 
                    Tissue = tissue_PCA,
                    Sample = sample_PCA)
# visualize the result of PCA of the top 2 PCs
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Tissue, colour = Time)) +
  ggtitle("PCA plot of the log-transformed expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = brewer.pal(8, "RdGy"),
                     name = "Time",
                     labels = c("10 days", "7 days",  "5 days", "3 days",  
                                "24 hrs", "12 hrs",  "6 hrs",   "0 hr"))



# output useful dataframes into csv files
write_csv(x = mw_wide, "mw_wide.csv")
write_csv(x = skin_long, "skin_long.csv")
  
  
  
  
  
  


