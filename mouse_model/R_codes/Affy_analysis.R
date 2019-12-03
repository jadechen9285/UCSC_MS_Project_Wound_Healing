# an end-to-end Affymetrix microarray differential expression workflow 
# using Bioconductor packages.


# author: Jianhong Chen
# date: 07/30/2019

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
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

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

# dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/test_micro_array"
# setwd(dir)

# create a temporary directory to store the raw data
raw_data_dir <- tempdir()
if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}

anno_AE = getAE("E-MTAB-2967", path = raw_data_dir, type ="raw")

# extra information from the raw data to construct a dataframe
# I. obtain info from SDRF file 
sdrf_location = file.path(raw_data_dir, "E-MTAB-2967.sdrf.txt") #similar to "paste" (but faster in R)
SDRF = read.delim(sdrf_location)

rownames(SDRF) = SDRF$Array.Data.File # replace generic machine label of the files
SDRF = AnnotatedDataFrame(SDRF) # create a biobased dataframe object

# II. create Expression Set object of array data
raw_data = oligo::read.celfiles(filenames = file.path(raw_data_dir, SDRF$Array.Data.File),
                                verbose = F, phenoData = SDRF)
# pre-selecte the useful corresponding columns for the datasets
Biobase::pData(raw_data) = Biobase::pData(raw_data)[, c("Source.Name", 
                                                        "Characteristics.individual.",
                                                        "Factor.Value.disease.",
                                                        "Factor.Value.phenotype.")]
#----------------
# quality control of the raw data

# I. perform log2 intensity value and PCA 
exp_raw = log2(Biobase::exprs(raw_data))
PCA_raw = prcomp(t(exp_raw), scale. = F) #note: t() = transpose 
# PCA, R use SVD to compute PCs, output singular value in "sdev"
percentVar = round(100*PCA_raw$sdev^2 / sum(PCA_raw$sdev^2), 1) #calculate the eigenalues percentage
sd_ratio = sqrt(percentVar[2]/percentVar[1]) # ratio between the top eigenvalues scores PCs

dataGG = data.frame(PC1 = PCA_raw$x[, 1], PC2 = PCA_raw$x[, 2],
                    Disease = pData(raw_data)$Factor.Value.disease., 
                    Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                    Individual = pData(raw_data)$Characteristics.individual.)
# visualize the result of PCA of the top 2 PCs
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

# boxplot of each microarray chip
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensities for the raw data")

# create a more detailted quality control report of the data set
# arrayQualityMetrics(expressionset = raw_data,
#                     outdir = tempdir(),
#                     force = TRUE, do.logtransform = TRUE,
#                     intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))


#----------------
# background adjustment, summarization and annotation
palmieri_eset = oligo::rma(raw_data, target = "core", normalize = F) # no normalization
# using robust multichip average(RMA) algorithm for background adjustment and summarization

# relative log expression (RLE) analysis
#  calculate the row medians of expression value for each microarray
row_medians_assayData = 
  Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))
RLE_data = sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

RLE_data_gathered = as.data.frame(RLE_data) %>%
  gather(patient_array, log2_expression_deviation)

ggplot(RLE_data_gathered, aes(x = patient_array, 
                              y =log2_expression_deviation )) +
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2,2)) +
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))

# RMA with background-correct, normalize and summarize:
palmieri_eset_norm = oligo::rma(raw_data, target = "core")

# PCA for the corrected dataset
exp_palmieri = Biobase::exprs(palmieri_eset_norm)
PCA = prcomp(t(exp_palmieri), scale = F)

# PCA, R use SVD to compute PCs, output singular value in "sdev"
percentVar = round(100*PCA$sdev^2 / sum(PCA$sdev^2), 1) #calculate the eigenalues percentage
sd_ratio = sqrt(percentVar[2]/percentVar[1]) # ratio between the top eigenvalues scores PCs

dataGG = data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2],
                    Disease = pData(raw_data)$Factor.Value.disease., 
                    Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                    Individual = pData(raw_data)$Characteristics.individual.)
# visualize the result of PCA of the top 2 PCs
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

# Heatmap clustering analysis
phenotype_names = ifelse(str_detect(pData(palmieri_eset_norm)$Factor.Value.phenotype.,
                                    "non"), "non_infl.", "infl.")
disease_names = ifelse(str_detect(pData(palmieri_eset_norm)$Factor.Value.disease.,
                                  "Crohn"), "CD", "UC")
annotation_for_heatmap = 
  data.frame(Phenotype = phenotype_names, Disease = disease_names)
row.names(annotation_for_heatmap) = row.names(pData(palmieri_eset_norm))

# compute the distance among sample to sample
dists = as.matrix(dist(t(exp_palmieri), method = "manhattan"))
rownames(dists) <- row.names(pData(palmieri_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) = NULL
diag(dists) = NA

ann_colors <- list(
  Phenotype = c(non_infl. = "chartreuse4", infl. = "burlywood3"),
  Disease = c(CD = "blue4", UC = "cadetblue2")
)

# visualize the results in kmeans cluster fashion
pheatmap(dists, col = (hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = T,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = T), max(dists, na.rm = T)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")


# -------------------------
# filtering based on intensity: filter out lowly expressed genes 
palmieri_medians = rowMedians(Biobase::exprs(palmieri_eset_norm))
man_threshold = 4 # mannual cutoff median intensity for filtering
hist_res = hist(palmieri_medians, 100, col = "cornsilk1", freq = F,
                main = "Histogram of the median intensities",
                border = "antiquewhite4",
                xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)

# now perform the filtering after the visualization
no_of_samples = table(paste0(pData(palmieri_eset_norm)$Factor.Value.disease., "_",
                            pData(palmieri_eset_norm)$Factor.Value.phenotype.))
samples_cutoff = min(no_of_samples) # threshold for smallest experimental group
# create a table to check which gene pass the mannual threshold of the median intensity
idx_man_threshold = apply(Biobase::exprs(palmieri_eset_norm), MARGIN = 1,
                          function(x){
                            sum(x > man_threshold) >= samples_cutoff
                          }) 

palmieri_manfiltered = subset(palmieri_eset_norm, idx_man_threshold)

# -------------------
# add annotation to the data:
# create df for annotation label of gene symbols and short description
anno_palmieri = AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                      keys = (featureNames(palmieri_manfiltered)),
                                      columns = c("SYMBOL", "GENENAME"),
                                      keytype = "PROBEID") 
anno_palmieri = subset(anno_palmieri, !is.na(SYMBOL)) # remove NA label 
anno_grouped = group_by(anno_palmieri, PROBEID) %>%
  dplyr::summarize(no_of_matches = n_distinct(SYMBOL)) # summary multiple mappings of each prob
probe_stats= filter(anno_grouped, no_of_matches > 1)
# for probes that mapped to multiple genes,
# it's difficult to decide which mapping is "correct", 
# therefore, we excludes these probes!!!
ids_to_exclude = (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)
palmieri_final = subset(palmieri_manfiltered, !ids_to_exclude) # remove multiple mappings
fData(palmieri_final)$PROBEID = rownames(fData(palmieri_final)) # fData() = featureData()
fData(palmieri_final) = left_join(fData(palmieri_final), anno_palmieri, by = "PROBEID")
# left_join keeps the rows and columns of the firt argument then add the corresponding
# column entries of the second argument.
rownames(fData(palmieri_final)) = fData(palmieri_final)$PROBEID

# ----------------------------
# linear models for microarrays

# first create abbreviations for labeling
individual = as.character(Biobase::pData(palmieri_final)$Characteristics.individual.)
tissue = str_replace_all(Biobase::pData(palmieri_final)$Factor.Value.phenotype.,
                         " ", "_") # conncet empty space with "_"
tissue = ifelse(tissue == "non-inflamed_colonic_mucosa", "nI", "I")
disease = str_replace_all(Biobase::pData(palmieri_final)$Factor.Value.disease.,
                          " ", "_")# conncet empty space with "_"
disease = ifelse(str_detect(Biobase::pData(palmieri_final)$Factor.Value.disease.,
                            "Crohn"), "CD", "UC")
  
# design metric for the model
# for "Crohn's disease
i_CD = individual[disease == "CD"]
design_palmieri_CD = model.matrix(~ 0 + tissue[disease == "CD"] + i_CD)
colnames(design_palmieri_CD)[1:2] = c("I", "nI")
rownames(design_palmieri_CD) = i_CD
# for Ulcerative colitis
i_UC = individual[disease == "UC"]
design_palmieri_UC = model.matrix(~0 + tissue[disease == "UC"] + i_UC)
colnames(design_palmieri_UC)[1:2] = c("I", "nI")
rownames(design_palmieri_UC) = i_UC

# analysis of differential expression based on a single gene: CRAT gene (PROBEID = 8164535)
tissue_CD = tissue[disease == "CD"]
crat_expr = Biobase::exprs(palmieri_final)["8164535", disease == "CD"]
crat_data = as.data.frame(crat_expr)
colnames(crat_data)[1] = "org_value"
crat_data = mutate(crat_data, individual = i_CD, 
                   tissue_CD = factor(tissue_CD, levels = c("nI", "I")))
ggplot(crat_data, aes(x = tissue_CD,y = org_value,
                      group = individual, color = individual)) + 
  geom_line() + 
  ggtitle("Expression changes for the CRAT gene")

# now fit a linear model
crat_coef = lmFit(palmieri_final[, disease == "CD"],
                  design = design_palmieri_CD)$coefficients["8164535", ]
crat_fitted = design_palmieri_CD %*% crat_coef # calculated the fitted value
rownames(crat_fitted) = names(crat_expr)
colnames(crat_fitted) = "fitted_value"
crat_data$fitted_value = crat_fitted
ggplot(crat_data, aes(x = tissue_CD, y = fitted_value,
                      group = individual, color = individual))+
  geom_line() +
  ggtitle("Fitted expression changes for the CRAT gene")

# perform paired t-test to check the difference between non-inflamed and 
#   inflamed tissued differs significantly from 0
crat_noninflamed = na.exclude(crat_data$org_value[tissue == "nI"])
crat_inflamed = na.exclude(crat_data$org_value[tissue == "I"])
res_t = t.test(crat_noninflamed, crat_inflamed, paired = T)
# p-value is this case is ~0, thus crat gene is differentially expressed between
#  non-inflamed and inflamed tissue.

# next: perform the same test for all genes
contrast_matrix_CD = makeContrasts(I-nI, levels = design_palmieri_CD)
palmieri_fit_CD = eBayes(contrasts.fit(lmFit(palmieri_final[, disease == "CD"],
                                             design = design_palmieri_CD),
                                       contrast_matrix_CD))
contrast_matrix_UC = makeContrasts(I-nI, levels = design_palmieri_UC)
palmieri_fit_UC = eBayes(contrasts.fit(lmFit(palmieri_final[, disease == "UC"],
                                             design = design_palmieri_UC),
                                       contrast_matrix_UC))
# eBayes = "Empirical Bayes" vairance moderation method to the model

# finally: extract the number of differentially expressed genes
table_CD = topTable(palmieri_fit_CD, number = Inf)
# the results sort by absolute value of t-value for each row (largest go first)
hist(table_CD$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflamed vs non-inflamed - Crohn's disease", xlab = "p-value")
table_UC = topTable(palmieri_fit_UC, number = Inf)
hist(table_UC$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "inflamed vs non-inflamed - Ulcerative colitis", xlab = "p-values")

# ---------------------------------
# Multiple testing FDR
# notice: in the orginal paper of these data set, a p-value of 0.001 was used as cut-off
total_genenumber_CD = length(subset(table_CD, P.Value < 0.001)$SYMBOL)
total_genenumber_UC = length(subset(table_UC, P.Value < 0.001)$SYMBOL)

# visualization of differentially expressed genes (volcano plot)
volcano_names = ifelse(abs(palmieri_fit_CD$coefficients)>=1,
                       palmieri_fit_CD$genes$SYMBOL, NA)
volcanoplot(palmieri_fit_CD, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names, xlab = "Log2 Fold Change",
            ylab = NULL, pch = 16, cex = 0.35)

# ------------------------------------
# 14.Gene ontology(GO) based enrichment analysis 
DE_genes_CD = subset(table_CD, adj.P.Val < 0.1)$PROBEID # this time we adjusted p-value for FDR cutoff
back_genes_idx = genefilter::genefinder(palmieri_final, as.character(DE_genes_CD),
                                        method = "manhattan", scale = "none")
# find background genes with similiar in expression of the DE genes
back_genes_idx = sapply(back_genes_idx, function(x) x$indices) #obtain the index of the background genes
back_genes = featureNames(palmieri_final)[back_genes_idx] %>%
  setdiff(DE_genes_CD) # separate DE genes from the background genes

intersect(back_genes, DE_genes_CD)# check if the separation is done correctly

# visualized if background genes show similiar results to DE genes
multidensity(list(all = table_CD[, "AveExpr"],
                  fore = table_CD[DE_genes_CD, "AveExpr"],
                  back = table_CD[rownames(table_CD) %in% back_genes, "AveExpr"]),
             col = c("#e46981", "#ae7ee2", "#a7ad4a"),
             xlab = "mean expression",
             main = "DE genes for CD-background-matching")

# actual GO genes tesing by topGO packages:
# the GO has three top ontologies: 1) Cellular component(CC),
#         2) biological processes(BP), and 3) molecular function(MF).

# I. create a named vector
gene_IDs = rownames(table_CD)
# find the genes index  that are universal for both background and DE enrichment analysis
in_universe = gene_IDs %in% c(DE_genes_CD, back_genes)
#  find the genes index only in the DE genes 
in_selection = gene_IDs %in% DE_genes_CD
# first select all the elements from "in_selection" that are in "in_universe
all_genes = in_selection[in_universe]
all_genes = factor(as.integer(in_selection[in_universe]))
names(all_genes) = gene_IDs[in_universe]

# initialize the togGo data set which contained in the annotation 
#   data base for the chip we are using here
top_GO_data = new("topGOdata", ontology = "BP", allGenes = all_genes,
                  nodeSize = 10, annot = annFUN.db,
                  affyLib = "hugene10sttranscriptcluster.db")

result_top_GO_elim = runTest(top_GO_data, algorithm = "elim",
                             statistic = "Fisher") # perform elim algorithm for statical test
result_top_GO_classic = runTest(top_GO_data, algorithm = "classic",
                                statistic = "Fisher")
res_top_GO = GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                      Fisher.classic = result_top_GO_classic,
                      orderBy = "Fisher.elim", topNodes = 100)
genes_top_GO = printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                          chip = "hugene10sttranscriptcluster.db", geneCutOff = 1000)
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
        collapse = "")
}) # notice: raw p-value = 2 means the genes are significant; 
# while raw p-value = 1 means the genes are insignificant.

# visualization fo the GO-analysis results
showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
               useInfo = 'def')




