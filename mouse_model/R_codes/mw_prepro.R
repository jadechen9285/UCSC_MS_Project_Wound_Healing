# I.
# mouse_wound Affy data pre-processing
# skin vs tongue in 8 time point 
# author: Jianhong Chen
# date: 07/31/2019

#######################################################################
# intial set up
ls() # list all objects
gc() # force R to release memory of deleted variables
rm(list = ls()) # delete all variables in the environment

# load the required packages
#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(pd.mouse430.2 )
#library(mogene20sttranscriptcluster.db)
library(mouse4302.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
library(splines)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

# data formating 
library(tidyverse)

#######################################################################
# download & load the data set 

# create a temporary directory to store the raw data
raw_data_dir <- tempdir()
if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}

anno_AE = getAE("E-GEOD-23006", path = raw_data_dir, type ="raw") # anno_AE download data from fttp link


# extra information from the raw data to construct a dataframe
# I. obtain info from SDRF file 
sdrf_location = file.path(raw_data_dir, "E-GEOD-23006.sdrf.txt") #similar to "paste" (but faster in R)
SDRF = read.delim(sdrf_location)

rownames(SDRF) = SDRF$Array.Data.File # replace generic machine label of the files
SDRF = AnnotatedDataFrame(SDRF) # create a biobased dataframe object

# II. create Expression Set object of array data
raw_data = oligo::read.celfiles(filenames = file.path(raw_data_dir, SDRF$Array.Data.File),
                                verbose = F, phenoData = SDRF)
# pre-selecte the useful corresponding columns for the datasets
Biobase::pData(raw_data) = Biobase::pData(raw_data)[, c("Source.Name", 
                                                        "Comment..Sample_source_name.",
                                                        "Comment..Sample_description.",
                                                        "FactorValue..TIME.POST.WOUND.",
                                                        "FactorValue..TISSUE.")]

#----------------
# quality control of the raw data

# I. perform log2 intensity value and PCA 
exp_raw = log2(Biobase::exprs(raw_data))
PCA_raw = prcomp(t(exp_raw), scale. = F) #note: t() = transpose 
# PCA, R use SVD to compute PCs, output singular value in "sdev"
percentVar = round(100*PCA_raw$sdev^2 / sum(PCA_raw$sdev^2), 1) #calculate the eigenalues percentage
sd_ratio = sqrt(percentVar[2]/percentVar[1]) # ratio between the top eigenvalues scores PCs

dataGG = data.frame(PC1 = PCA_raw$x[, 1], PC2 = PCA_raw$x[, 2],
                    Time = pData(raw_data)$FactorValue..TIME.POST.WOUND., 
                    Tissue = pData(raw_data)$FactorValue..TISSUE.,
                    Sample = pData(raw_data)$Source.Name)

# visualize the result of PCA of the top 2 PCs
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Tissue, colour = Time)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = brewer.pal(8, "RdGy"),
                     name = "Time",
                     labels = c("10 days", "7 days",  "5 days", "3 days",  
                                "24 hrs", "12 hrs",  "6 hrs",   "0 hr"))

# boxplot of each microarray chip
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensities for the raw data")

# create a html report for the quality control of the raw data
# testdir = "/Users/jianhongchen/Desktop"
# arrayQualityMetrics(expressionset = raw_data,
#                     outdir = testdir,
#                     force = TRUE, do.logtransform = TRUE,
#                     intgroup = c("FactorValue..TIME.POST.WOUND.", 
#                                  "FactorValue..TISSUE."))
#----------------
# background adjustment, summarization and annotation
mw_eset = oligo::rma(raw_data, normalize = F) # no normalization
# using robust multichip average(RMA) algorithm for background adjustment,
# normalization and summarization

# relative log expression (RLE) analysis
#  calculate the row medians of expression value for each microarray
row_medians_assayData = 
  Biobase::rowMedians(as.matrix(Biobase::exprs(mw_eset)))
RLE_data = sweep(Biobase::exprs(mw_eset), 1, row_medians_assayData)

RLE_data_gathered = as.data.frame(RLE_data) %>%
  gather(sample_array, log2_expression_deviation)

ggplot(RLE_data_gathered, aes(x = sample_array, 
                              y =log2_expression_deviation )) +
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2.5,2.5)) +
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"),
        text = element_text(size = 15))

# RMA with background-correct, normalize and summarize:
mw_eset_norm = oligo::rma(raw_data)

# PCA for the corrected dataset
exp_mw = Biobase::exprs(mw_eset_norm)
PCA = prcomp(t(exp_mw), scale. = F)

# PCA, R use SVD to compute PCs, output singular value in "sdev"
percentVar = round(100*PCA$sdev^2 / sum(PCA$sdev^2), 1) #calculate the eigenalues percentage
sd_ratio = sqrt(percentVar[2]/percentVar[1]) # ratio between the top eigenvalues scores PCs
dataGG = data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2],
                    Time = pData(raw_data)$FactorValue..TIME.POST.WOUND., 
                    Tissue = pData(raw_data)$FactorValue..TISSUE.,
                    Sample = pData(raw_data)$Comment..Sample_description.)

# visualize the result of PCA of the top 2 PCs
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Tissue, colour = Time)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = brewer.pal(8, "RdGy"),
                     name = "Time",
                     labels = c("10 days", "7 days",  "5 days", "3 days",  
                                "24 hrs", "12 hrs",  "6 hrs",   "0 hr"))

# Heatmap clustering analysis
annotation_for_heatmap = 
  data.frame(Tissue = pData(mw_eset_norm)$FactorValue..TISSUE., 
             Time = pData(mw_eset_norm)$FactorValue..TIME.POST.WOUND.)
row.names(annotation_for_heatmap) = row.names(pData(mw_eset_norm))
# compute the distance among sample to sample
dists = as.matrix(dist(t(exp_mw), method = "manhattan"))
rownames(dists) <- row.names(pData(mw_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) = NULL
diag(dists) = NA

ann_time = levels(annotation_for_heatmap$Time)
ann_colors <- list(
  Tissue = c(skin = "chartreuse4", "mucosa (tongue)" = "burlywood3"))

  # Time = c("10 days" = brewer.pal(8, "YlOrRd")[1], "7 days" = brewer.pal(8, "YlOrRd")[2],
  #             "5 days" = brewer.pal(8, "YlOrRd")[3], "3 days" = brewer.pal(8, "YlOrRd")[4],
  #             "24 hrs" = brewer.pal(8, "YlOrRd")[5], "12 hrs" = brewer.pal(8, "YlOrRd")[6],
  #             "6 hrs" = brewer.pal(8, "YlOrRd")[7], "0 hrs" = brewer.pal(8, "YlOrRd")[8]
  #             ))


# visualize the results in kmeans cluster fashion
pheatmap(dists, col = (hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = T,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = T), max(dists, na.rm = T)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")

# ------------------------- !![skip manual filter on low expressed gens]]!!
# filtering based on intensity: filter out lowly expressed genes 
# mw_medians = rowMedians(Biobase::exprs(mw_eset_norm))
# man_threshold = 4 # mannual cutoff median intensity for filtering
# hist_res = hist(mw_medians, 100, col = "cornsilk1", freq = F,
#                 main = "Histogram of the median intensities",
#                 border = "antiquewhite4",
#                 xlab = "Median intensities")
# abline(v = man_threshold, col = "coral4", lwd = 2)
# 
# # now perform the filtering after the visualization
# no_of_samples = table(paste0(pData(mw_eset_norm)$FactorValue..TISSUE., "_",
#                              pData(mw_eset_norm)$FactorValue..TIME.POST.WOUND.))
# samples_cutoff = min(no_of_samples) # threshold for smallest experimental group
# # create a table to check which gene pass the mannual threshold of the median intensity
# idx_man_threshold = apply(Biobase::exprs(mw_eset_norm), MARGIN = 1,
#                           function(x){
#                             sum(x > man_threshold) >= samples_cutoff
#                           })
# 
# mw_manfiltered = subset(mw_eset_norm, idx_man_threshold)

# -------------------
# add annotation to the data:
# create df for annotation label of gene symbols and short description
anno_mw = AnnotationDbi::select(mouse4302.db,
                                      keys = (featureNames(mw_eset_norm)),
                                      columns = c("SYMBOL", "GENENAME"),
                                      keytype = "PROBEID") 
anno_mw = subset(anno_mw, !is.na(SYMBOL)) # remove NA label 
# 1. get rid of one prob ID maps to multiple genes
anno_grouped = group_by(anno_mw, PROBEID) %>%
  dplyr::summarize(no_of_matches = n_distinct(SYMBOL)) # summary multiple mappings of each prob
probe_stats= filter(anno_grouped, no_of_matches > 1) # extra multiple mapping genes
# for probes that mapped to multiple genes,it's difficult to decide which mapping is "correct", 
# therefore, we excludes these probes!!!
ids_to_exclude = (featureNames(mw_eset_norm) %in% probe_stats$PROBEID)
mw_filter1 = subset(mw_eset_norm, !ids_to_exclude) # remove multiple mappings
fData(mw_filter1)$PROBEID = rownames(fData(mw_filter1)) # fData() = featureData()
fData(mw_filter1) = left_join(fData(mw_filter1), anno_mw, by = "PROBEID")
# left_join keeps the rows and columns of the firt argument then add the corresponding

# 2. get rid of multiple same genes that map to different probe IDS
gene_unique = unique(fData(mw_filter1)$SYMBOL)
filter1_grouped = group_by(fData(mw_filter1), SYMBOL) %>%
  dplyr::summarise(no_of_matches = n_distinct(PROBEID))
gene_stats = filter(filter1_grouped, no_of_matches > 1)
genes_to_exclude = (fData(mw_filter1)$SYMBOL %in% gene_stats$SYMBOL)
# the final data set contain 1-to-1 probe:gene maping
mw_final = subset(mw_filter1, !genes_to_exclude) 
# column entries of the second argument.
rownames(fData(mw_final)) = fData(mw_final)$PROBEID

# ----------------------------
# linear models for microarrays
#     
# using time course fiting (limma user guide 9.6.1)
lvl = levels(pData(mw_final)$Comment..Sample_source_name.)
sample_lab = factor(pData(mw_final)$Comment..Sample_source_name., levels = lvl)
design = model.matrix(~0 + sample_lab)
lvl = str_remove_all(lvl, "post-wounding") %>%
  str_remove_all("[[:punct:]]")%>%
  str_remove_all(" ")
colnames(design) = lvl
fit = lmFit(mw_final, design)

# create contract matrix
cont_tongue = makeContrasts("tongue10days-unwoundedtongue", "tongue7days-unwoundedtongue",
                            "tongue5days-unwoundedtongue", "tongue3days-unwoundedtongue",
                            "tongue24hrs-unwoundedtongue", "tongue12hrs-unwoundedtongue",
                            "tongue6hrs-unwoundedtongue",
                            levels = design)

cont_skin = makeContrasts("skin10days-unwoundedskin", "skin7days-unwoundedskin",
                          "skin5days-unwoundedskin", "skin3days-unwoundedskin",
                          "skin24hrs-unwoundedskin", "skin12hrs-unwoundedskin",
                          "skin6hrs-unwoundedtongue",
                          levels = design)

fit_skin = eBayes(contrasts.fit(fit, cont_skin))
fit_tongue = eBayes(contrasts.fit(fit, cont_tongue))
skin_table = topTable(fit_skin, number = Inf, adjust="BH") # applied "BH" adjusted p-value
tongue_table = topTable(fit_tongue, number = Inf, adjust = "BH")

# --the results sort by absolute value of t-value for each row (largest go first)
hist(skin_table$adj.P.Val, col = brewer.pal(3, name = "Set2")[1],
     main = "Histogram of All Data in Adjusted P-value -- Skin Tissue", 
     xlab = "Adjusted p-value", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(tongue_table$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "Histogram of All Data in Adjusted P-value -- Tongue Tissue", 
     xlab = "Adjusted p-value", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# filtered skin dataframe (adj.p-value <1E-5)
skin_table_filtered = subset(skin_table, adj.P.Val < 0.00001) %>%
  filter(!is.na(SYMBOL))
tongue_table_filtered= subset(tongue_table, adj.P.Val < 0.00001) %>%
  filter(!is.na(SYMBOL))
# -- visualize the filtered data distribution via histogram
hist(skin_table_filtered$adj.P.Val, col = brewer.pal(3, name = "Set2")[1],
     main = "Histogram of Filtered Data -- Skin Tissue", 
     xlab = "Adjusted p-value", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(tongue_table_filtered$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "Histogram of Filtered Data -- Tongue Tissue", 
     xlab = "Adjusted p-value", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

# visualization of differentially expressed genes (volcano plot)
skin_fit_genes = ifelse(abs(fit_skin$coefficients)>=2,
                        fit_skin$genes$SYMBOL, NA)
volcanoplot(fit_skin[,7], coef = 1L, style = "p-value", highlight = 100,
            names = skin_fit_genes[, 7], xlab = "Log2 Fold Change",
            main = "Differentally Expressed Genes at T7(FC > 2) -- Skin",
            ylab = NULL, pch = 16, cex = 0.35,
            xlim = c(-5,5), 
            cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(v = -2, col = "coral4", lwd = 2)
abline(v = 2, col = "coral4", lwd = 2)

tongue_fit_genes = ifelse(abs(fit_tongue$coefficients)>=2,
                        fit_tongue$genes$SYMBOL, NA)
volcanoplot(fit_tongue[,7], coef = 1L, style = "p-value", highlight = 1000,
            names = tongue_fit_genes[,7], xlab = "Log2 Fold Change",
            main = "Differentally Expressed Genes at T7(FC > 2) -- Tongue",
            ylab = NULL, pch = 16, cex = 0.35,
            xlim = c(-5,5),
            cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(v = -2, col = "coral4", lwd = 2)
abline(v = 2, col = "coral4", lwd = 2)


# ----------------------------
# final extraction of the genes for both skin and tongue tissues
skin_probes = unique(skin_table_filtered$PROBEID)
tongue_probes = unique(tongue_table_filtered$PROBEID)
all_unique_probes = unique(c(skin_probes, tongue_probes))

skin_clean = subset(mw_final, pData(mw_final)$FactorValue..TISSUE. == "skin")
tongue_clean = subset(mw_final, pData(mw_final)$FactorValue..TISSUE. == "mucosa (tongue)")

mw_clean = subset(mw_final, featureNames(mw_final) %in% all_unique_probes)

# mw_clean_df = data.frame(exprs(mw_clean)) %>%
#   mutate("symbol" = fData(mw_clean)$SYMBOL,
#          "gene_des" = fData(mw_clean)$GENENAME)

# extracting high differentially expressed genes based on the foldchange threshold
skin_DE_gene = unique(as.vector(skin_fit_genes))[-1] # [-1] remove NA term
tongue_DE_gene = unique(as.vector(tongue_fit_genes))[-1] #[-1] remove NA term
all_DE_gene = unique(c(skin_DE_gene, tongue_DE_gene))

mw_DE = subset(mw_final, fData(mw_final)$SYMBOL %in% all_DE_gene)
mw_DE_df = data.frame(exprs(mw_DE)) %>%
  mutate("probeID" = fData(mw_DE)$PROBEID,
        "gene_symbol" = fData(mw_DE)$SYMBOL, 
         "gene_des" = fData(mw_DE)$GENENAME)

# convert final dataframe into human readable annotation format
desc_label = c("tongue 0 1", "tongue 0 2", "tongue 0 3", "skin 0 1", "skin 0 2", "skin 0 3",
               "tongue 6 1", "tongue 6 2", "tongue 6 3", "skin 6 1", "skin 6 2", "skin 6 3",
               "tongue 12 1", "tongue 12 2", "tongue 12 3", "skin 12 1", "skin 12 2", "skin 12 3",
               "tongue 24 1", "tongue 24 2", "tongue 24 3", "skin 24 1", "skin 24 2", "skin 24 3",
               "tongue 72 1", "tongue 72 2", "tongue 72 3", "skin 72 1", "skin 72 2", "skin 72 3",
               "tongue 120 1", "tongue 120 2", "tongue 120 3", "skin 120 1", "skin 120 2", "skin 120 3",
               "tongue 168 1", "tongue 168 2", "tongue 168 3", "skin 168 1", "skin 168 2", "skin 168 3",
               "tongue 240 1", "tongue 240 2", "tongue 240 3", "skin 240 1", "skin 240 2", "skin 240 3")
col_label = rev(desc_label)
colnames(mw_DE_df)[1:48] = col_label
mw_DE_long = gather(mw_DE_df, key = "desc", value = "intensity", 1:48)
# create the label for the description of the dataset
tissue = map_chr(strsplit(mw_DE_long$desc, split = " "), function(x) x[1])
time = factor(map_chr(strsplit(mw_DE_long$desc, split = " "), function(x) x[2]), 
              levels = c(0,6, 12, 24,72,120,168,240))
sample = map_chr(strsplit(mw_DE_long$desc, split = " "), function(x) x[3])

mw_DE_final = mw_DE_long[, -4] %>%
  mutate(tissue, time , sample)

# # try to select genes that are part of the macrophage family
# skin_macro = filter(skin_DE_df, grepl("macrophage", skin_DE_df$gene_des))
# tongue_macro = filter(tongue_DE_df, grepl("macrophage", tongue_DE_df$gene_des))
# 
# ------------------------------
# output data into csv file
dir2 = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/cleaned_data"
setwd(dir2)
write_csv(x = mw_DE_final, "mw_DE_fc1_v1.csv")

