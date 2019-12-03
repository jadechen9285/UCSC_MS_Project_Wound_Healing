# ----------------------------------------------------------
# complie CEL files into a csv or txt file
# mouse skin vs tongue tissue sample in 8 time point(in hr): 0, 6, 12, 24, 72, 120, 168, 240
# 
# summer 2019 wound healing research program
# author: Jianhong Chen
# Date : 07_23_2019


#library(oligo)

#library(simpleaffy)
library(affy)
library('ggplot2')
library('tidyverse') # inlcude "tidyr", "readr", "dplyr", "purrr"
library("readxl")
library(tidyverse)

rm(list=ls())
gc()

dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/cleaned_data"
setwd(dir) 
#  using oligo & simpleaffy library method
# cel_files = list.celfiles()
# affyraw = read.celfiles(cel_files)
# 
# library('pd.mogene.2.0.st') 
# #load the specific annotation for the target microarray chip
# exp_rma = oligo::rma(affyraw)

# use affy library:
# mydata = ReadAffy()
# eset = rma(mydata)
# eset2 = mas5(mydata)
# eset2_PMA = mas5calls(mydata)
# test = data.frame(exprs(eset2))
#  Exporting data from ExpressionSet objects
# x = data.frame(exprs(eset2), exprs(eset2_PMA), assayDataElement(eset2_PMA, "se.exprs"))
# x <- x[,sort(names(x))]
# eset2_p = data.frame(assayDataElement(eset2_PMA, "se.exprs"))
# output the data with gene expression lvl and p-value in excel files
mw_clean1 = read_csv("mw_clean1.csv")
mouse_wound = read_csv("mouse_wound.csv")
# -------------------------------------------------
# data preprocessing before analysis
# clean up label of p-value data frame
# create the new description label: based on (tissue type, time(hr) , sample#)
desc_label = c("tongue 0 1", "tongue 0 2", "tongue 0 3", "skin 0 1", "skin 0 2", "skin 0 3",
               "tongue 6 1", "tongue 6 2", "tongue 6 3", "skin 6 1", "skin 6 2", "skin 6 3",
               "tongue 12 1", "tongue 12 2", "tongue 12 3", "skin 12 1", "skin 12 2", "skin 12 3",
               "tongue 24 1", "tongue 24 2", "tongue 24 3", "skin 24 1", "skin 24 2", "skin 24 3",
               "tongue 72 1", "tongue 72 2", "tongue 72 3", "skin 72 1", "skin 72 2", "skin 72 3",
               "tongue 120 1", "tongue 120 2", "tongue 120 3", "skin 120 1", "skin 120 2", "skin 120 3",
               "tongue 168 1", "tongue 168 2", "tongue 168 3", "skin 168 1", "skin 168 2", "skin 168 3",
               "tongue 240 1", "tongue 240 2", "tongue 240 3", "skin 240 1", "skin 240 2", "skin 240 3")
col_label = rev(desc_label)
colnames(mw_clean1)[1:48] = col_label

# filter skin and tongue tissue to apply fold change >2 filtering criteria.
skin_clean1 = mw_clean1[ ,grepl("skin", colnames(mw_clean1))]
tongue_clean1 = mw_clean1[, grepl("tongue", colnames(mw_clean1))]

# skin_s1 = skin_clean1[, c(1, 4, 7, 10, 13, 16, 19, 22)]
# skin_s1 = filter(skin_s1,
#                  abs(skin_s1[,1] - skin_s1[, 8]) >= 1 | abs(skin_s1[,2] - skin_s1[, 8]) >= 1 |
#                   abs(skin_s1[,3] - skin_s1[, 8]) >= 1 | abs(skin_s1[,4] - skin_s1[, 8]) >= 1 |
#                   abs(skin_s1[,5] - skin_s1[, 8]) >= 1 | abs(skin_s1[,6] - skin_s1[, 8]) >= 1 |
#                   abs(skin_s1[,7] - skin_s1[, 8]) >= 1 )# fulfill condition 1 and filter significance differentiable genes
# 
# skin_s2 = skin_clean1[, c(2, 5, 8, 11, 14, 17, 20, 23)]
# skin_s2 = filter(skin_s2,
#                  abs(skin_s2[,1] - skin_s2[, 8]) >= 1 | abs(skin_s2[,2] - skin_s2[, 8]) >= 1 |
#                    abs(skin_s2[,3] - skin_s2[, 8]) >= 1 | abs(skin_s2[,4] - skin_s2[, 8]) >= 1 |
#                    abs(skin_s2[,5] - skin_s2[, 8]) >= 1 | abs(skin_s2[,6] - skin_s2[, 8]) >= 1 |
#                    abs(skin_s2[,7] - skin_s2[, 8]) >= 1 )# fulfill condition 1 and filter significance differentiable genes
# 



mw_clean1_long = gather(mw_clean1, key = "desc", value = "intensity", 1:48)
# create the label for the description of the dataset
tissue = map_chr(strsplit(mw_clean1_long$desc, split = " "), function(x) x[1])
time = factor(map_chr(strsplit(mw_clean1_long$desc, split = " "), function(x) x[2]), 
              levels = c(0,6, 12, 24,72,120,168,240))
sample = map_chr(strsplit(mw_clean1_long$desc, split = " "), function(x) x[3])



mw_final = mw_clean1_long[, -3] %>%
  mutate(tissue, time , sample)

group1_keyword = c("wound", "inflammatory", "cytokine", "toll-like receptor",
                 "jak-stat", "chemotaxis")  # early up-regulated genes (table#3,4)
group2_keyword = c("transcription regulation", "DNA binding") # early down-regulated genes (table #5,6)
group3_keyword = c("extracellular matrix", "collagen", "structural protein", "ECM-receptor",
                   "cell communication", "peptidase activity") # late up-regulated genes (table#7,8)

group_keywords = c("inflammatory cytokines", "neutrophils", "macrophages",
                   "angiogenesis")

g1V2 = filter(mouse_wound, grepl("macrophage", mouse_wound$`Gene Title`)) %>%
  select(c(`Gene Symbol`, `Gene Title`, tissue, time, sample, intensity)) %>%
  mutate("group" = rep(1, dim(g1V2)[1]))
g2V2 = filter(mouse_wound, grepl("angiogenesis", mouse_wound$`Gene Title`)) %>%
  select(c(`Gene Symbol`, `Gene Title`, tissue, time, sample, intensity)) %>%
  mutate("group" = rep(2, dim(g2V2)[1]))
g3V2 = filter(mouse_wound, grepl("cytokine", mouse_wound$`Gene Title`)) %>%
  select(c(`Gene Symbol`, `Gene Title`, tissue, time, sample, intensity)) %>%
  mutate("group" = rep(3, dim(g3V2)[1]))
g4V2 = filter(mouse_wound, grepl("neutrophils", mouse_wound$`Gene Title`)) # no genes associated with this keywords

groupV2 = rbind(g1V2, g2V2, g3V2)


group1_df = filter(mw_final, grepl(paste(group1_keyword, collapse = "|"), 
                                   mw_final$gene_des))
group1_df = mutate(group1_df, "group" = rep(1, dim(group1_df)[1]))

group2_df = filter(mw_final, grepl(paste(group2_keyword, collapse = "|"), 
                                   mw_final$gene_des))
group2_df = mutate(group2_df, "group" = rep(2, dim(group2_df)[1]))

group3_df = filter(mw_final, grepl(paste(group3_keyword, collapse = "|"), 
                                   mw_final$gene_des))
group3_df = mutate(group3_df, "group" = rep(3, dim(group3_df)[1]))

group_final_df = rbind(group1_df, group2_df, group3_df)


# -----------------------------------------------------------------------
skin_final = filter(mw_final, tissue == "skin")
skin_list = skin_cl_list = list()
time_lab = unique(skin_final$time)
for (t in 1:8) {
  skin_list[[t]] = filter(skin_final, time == time_lab[t])
  skin_cl_list[[t]] = kmeans(skin_list[[t]]$intensity, centers = 5)
  skin_list[[t]] = mutate(skin_list[[t]], "cluster" = factor(skin_cl_list[[t]]$cluster))
}

# create a dataframe that contains the centroild for the clustered results
skin_center = data.frame("center" = skin_cl_list[[1]]$centers)

for (t in 2:8){
  skin_center = cbind(skin_center, skin_cl_list[[t]]$centers)
}
colnames(skin_center) = time_lab

# normalize the centroids of the clustered results
col_mean = colMeans(skin_center)
col_sd = sapply(skin_center, sd)

for (t in 1:8){
  skin_center[, t] = (skin_center[, t] - col_mean[t]) / col_sd[t]
}

skin_center_long = mutate(skin_center, "center" = factor(1:5)) %>%
  gather("time", "value", -center) 

ggplot(skin_center_long, aes(x = time, y = value, group = center, color = center)) +
  geom_point() + 
  geom_line()

# mw_df_long = gather(mw_df, key = "desc", value = "intensity", names(mw_df)[5]:names(mw_df)[52])
# 
# # construct the clean-up data based on the two condtions by the authors
# # 1) at least 1 time point is 2 fold change (log2 scale) of intensity valeu compare to normal tissue (t0)
# # 2) p_adjust (FDR corrected) should be < 1E-5 or 0.05 ??????
# 
# BH_test = p.adjust(eset2_p_long$p_value, method= "BH")
# mw_long= mw_df_long[, -c(5,6)] %>%
#   mutate(tissue, time, sample, 
#          "intensity" = mw_df_long$intensity, "p_value" = eset2_p_long$p_value,
#          "p_adjust" = BH_test)
# 
# # need to make the dataframe wide for time to satisfy conditon 1
# wide_inten_lab = c("intensity_t1","p_adjust_t1", "intensity_t2","p_adjust_t2",
#                    "intensity_t3","p_adjust_t3", "intensity_t4","p_adjust_t4",
#                    "intensity_t5","p_adjust_t5", "intensity_t6","p_adjust_t6",
#                    "intensity_t7","p_adjust_t7", "intensity_t8","p_adjust_t8")
# t_label = unique(mw_clean_df$time)
# mw_wide = filter(mw_long, time == 0)[, c(1:5,7)]
# wide_index = grep(paste(c("intensity", "p_adjust"), collapse = "|"), names(filter(mw_long, time ==0)))
# for (t in 1:length(t_label)){ 
#   mw_wide = cbind(mw_wide, filter(mw_long, time == t_label[t])[,wide_index]) # construct wide data structure
# }
# names(mw_wide)[7:22] = wide_inten_lab
# 
# 
# mw_clean1 = filter(mw_wide, 
#       abs(intensity_t2 - intensity_t1) >= 1 | abs(intensity_t3 - intensity_t1) >= 1 |
#       abs(intensity_t4 - intensity_t1) >= 1 | abs(intensity_t5 - intensity_t1) >= 1 |
#       abs(intensity_t6 - intensity_t1) >= 1 | abs(intensity_t7 - intensity_t1) >= 1 |
#       abs(intensity_t8 - intensity_t1) >= 1 ) # fulfill condition 1 and filter significance differentiable genes
# 
# # need some tricks to compile the data into long format in order to filter based on p_value
# time2 = c(rep(0, dim(mw_clean1)[1]),rep(6, dim(mw_clean1)[1]),rep(12, dim(mw_clean1)[1]),
#           rep(24, dim(mw_clean1)[1]),rep(72, dim(mw_clean1)[1]),rep(120, dim(mw_clean1)[1]),
#           rep(168, dim(mw_clean1)[1]),rep(240, dim(mw_clean1)[1]))
# p_adjust_long = mw_clean1[, sort(names(mw_clean1))] %>% # sort out the name label of the dataframe
#   gather(p_time, p_adjust, 12:19) # only compile the adjusted p-values from different time pt
#                      
# intensity_long = mw_clean1[, sort(names(mw_clean1))] %>% # sort out the name label of the dataframe
#   gather(int_time, intensity, 4:11) # only compile the adjusted p-values from different time pt
# 
# mw_clean2 = data.frame("gene" = p_adjust_long$`Gene Title`, 
#                        "tissue" = p_adjust_long$tissue, 
#                        "sample" = p_adjust_long$sample) %>%
#   mutate("time" = time2,"intensity" = intensity_long$intensity, 
#          "p_adjust" = p_adjust_long$p_adjust) %>%
#   filter(p_adjust < 0.05)
# 
# gene_count = length(unique(mw_clean2$`Gene Symbol`))
# 
# group1_keyword = c("wound", "inflammatory", "cytokine", "toll-like receptor",
#                  "jak-stat", "chemotaxis") # early up-regulated genes (table#3,4)
# group2_keyword = "DNA binding" # early down-regulated genes (table #5,6)
# group3_keyword = c("extracellular matrix", "collagen", "structural protein", "ECM-receptor", 
#                    "cell communication", "peptidase activity") # late up-regulated genes (table#7,8)
# 
# group1_df = filter(mw_clean2, grepl(paste(group1_keyword, collapse = "|"), mw_clean2$gene))
# group2_df = filter(mw_clean2, grepl(paste(group2_keyword, collapse = "|"), mw_clean2$gene))
# group3_df = filter(mw_clean2, grepl(paste(group3_keyword, collapse = "|"), mw_clean2$gene))
# 




# output all datafiles 
dir2 = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/cleaned_data"
setwd(dir2)
write_csv(x = mw_final, "mw_final.csv")
write_csv(x = group_final_df, "groupV1.csv")
write_csv(x = groupV2, "groupV2.csv" )
