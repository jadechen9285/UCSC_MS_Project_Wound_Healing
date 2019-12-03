########################################################################
#
# mouse_wound Affy data quick pre-processing without any plots
# skin vs tongue in 8 time point 
# author: Jianhong Chen
# date: 09/25/2019
# this script focus on the tongue sample
#
#######################################################################
# intial set up
ls() # list all objects
gc() # force R to release memory of deleted variables
rm(list = ls()) # delete all variables in the environment
set.seed(78)

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
library('plotly') # 3D plotting

# data formating 
library(tidyverse)

# clustering
library("kernlab") # kernel-based
library(clValid) # clustering evaulation

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

# RMA with background-correct, normalize and summarize:
mw_eset_norm = oligo::rma(raw_data)

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
# linear models for microarrays: extracting differentally expressed genes
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

# --------------------------------------------------------
# load the preprossed data set
dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/summer_2019/data/cleaned_data"
setwd(dir)
mw_df = read_csv("mw_DE_fc1_v1.csv")
# change time to factor value first
time_lab = c("0hr", "6hrs", "12hrs", "24hrs", "3days", "5days", "7days", "10days")
mw_df$time = factor(mw_df$time,levels = c(0,6, 12, 24,72,120,168,240))

mw_wide = spread(mw_df, key = time, value = intensity)
colnames(mw_wide)[-(1:5)] = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")

# compute the difference in intensity bewteen each wounding time vs unwound time points
mw_diff = matrix(0, dim(mw_wide)[1], length(time_lab)-1) # create the 0 matrix with the right size
for (i in 1:(length(time_lab)-1)){
  mw_temp = mw_wide[, -(1:5)] # create dummy stroage df
  mw_diff[, i] = unlist(log2(abs(mw_temp[, i+1] - mw_temp[, 1])+1)) #log2 transformation
  rm(mw_temp) # remove dummy df
}
colnames(mw_diff) = c("d1", "d2", "d3", "d4", "d5","d6", "d7")
mw_diff = cbind(mw_wide[,1:5], data.frame(mw_diff))

# ------------------------------------------------------------
# perform clustering: focus on the tongue tissue
tongue_diff = filter(mw_diff, tissue =="tongue")

# ----elbow method--------
# optimization k-value of the clustering algorithm 
total_withinSS_km = map_dbl(1:20, function(k){
  model_km = kmeans(tongue_diff[, -(1:5)], centers = k)
  model_km$tot.withinss
})

km_elbow = data.frame( 'k' = 1:20, 'tot_withinSS' = total_withinSS_km)
# the plot implies k = 2 as the optimized value
ggplot(km_elbow, aes(x = k, y = tot_withinSS)) +
  ggtitle("elbow plot for k-means") +
  geom_line() +
  scale_x_continuous(breaks = 1:20) + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20))

# set K = number of clusters
#------ input parameters
K = 5
cl_input = scale(tongue_diff[, -c(1:5)])
cl_dist = dist(cl_input, method = "euclidean")
#cl_input = filter(mw_wide, tissue == "skin")[, -c(1:5)]
# -----
tongue_hcl = hclust(cl_dist, method = "average")
#plot(as.dendrogram(tongue_hcl))
tongue_hcl_cl = cutree(tongue_hcl, k = K)
tongue_kkm = kkmeans(as.matrix(cl_input), centers = K, kernel = "rbfdot")
tongue_km = kmeans(cl_input, centers = K)

# -- compute Dunn Index
hcl_dunn = dunn(distance = cl_dist, clusters = tongue_hcl_cl)
km_dunn = dunn(Data = cl_input, clusters = tongue_km$cluster)
kkm_dunn = dunn(Data = cl_input, clusters = tongue_kkm)

tongue_PCA = prcomp(tongue_diff[, -(1:5)], scale = T) # prcomp uses svd for PCA
#calculate the eigenalues percentage
PCA_eigen_per = percentVar = round(100*tongue_PCA$sdev^2 / sum(tongue_PCA$sdev^2), 1)
# cluster results dataframe
cl_res = factor(tongue_km$cluster)

tongue_norm = filter(mw_wide, tissue == "tongue") %>%
  mutate("cluster" = cl_res)

tongue_GG = cbind(tongue_diff[,c(1:3,5)], data.frame(tongue_PCA$x)) %>%
  mutate("cluster" = cl_res)

tongue_df = filter(mw_df, tissue == "tongue" & sample == 1) %>%
  mutate("cluster" = rep(filter(tongue_GG, sample ==1)$cluster, 8))

# plot1: cluster plot
ggplot(tongue_norm, aes(t1, t2)) +
  geom_point(alpha = 0.6, shape =20, aes(colour = cluster)) +
  ggtitle("Visualization in norm axis for Km clustering(log2) of Tongue(All Samples)") +
  xlab("t1") +
  ylab("t2") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20))+
  scale_color_manual(values = brewer.pal(8, "Set2"))

ggplot(tongue_GG, aes(PC1, PC2)) +
  geom_point(alpha = 0.6, shape =20, aes(colour = cluster)) +
  ggtitle("Visualization in PCA axis for Km clustering(log2) of Tongue(All Samples)") +
  xlab(paste0("PC1, VarExp: ", PCA_eigen_per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", PCA_eigen_per[2], "%")) +
  ylim(-6.5, 6.5)+
  xlim(-4,12) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20))+
  scale_color_manual(values = brewer.pal(8, "Set2"))

# plot2: boxplot for case study genes in time 
ggplot(tongue_df, aes(x = time, y = intensity)) +
  geom_boxplot() +
  facet_wrap(~cluster) + 
  ggtitle("Boxplot for Tongue(sample 1) in Km Clustering") + 
  xlab("time(hrs)") + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20))

# --------------------------------------------------------------------------------
# perform Gene ontology (GO) based enrichment analysis for 1 sample at a time

# -----------function
TopGO_Analysis = function(DE_genes_id, gene_table, clean_expr){
  # -----------------------
  # function perform TopGO analysis!!!
  #     require differentially expressed genes and background genes for contract
  #     inputs: DE_genes_id = probe ID of the differentially expressed genes
  #             gene_table = a table that contractins all the genes and the significan p-value
  #             clean_expr = cleaned expression set data structure 
  
  #     output: return topGO object that contains the results
  # ----------------------
  # matching DE genes wrt the background set of genes
  back_genes_idx = genefilter::genefinder(clean_expr,
                                          as.character(DE_genes_id),
                                          method = "manhattan",
                                          scale = "none")
  
  # extract PROBE IDs for the selected backgroun genes indicies
  back_genes_idx = sapply(back_genes_idx, function(x)x$indices)
  back_genes = featureNames(clean_expr)[back_genes_idx]
  # setdiff helps eliminate foreground genes 
  back_genes = setdiff(back_genes, DE_genes_id)
  # intersect check if the separation was done successfully
  # return 0 == correct; >0 == incorrect
  intersect(back_genes, DE_genes_id)
  
  # visualized the result: multidensity plot with mean expression 
  multidensity(list(
    all = gene_table[, "AveExpr"],
    fore = gene_table[DE_genes_id, "AveExpr"],
    back = gene_table[rownames(gene_table) %in% back_genes, "AveExpr"]),
    col = c("#e46981", "#ae7ee2", "#a7ad4a"),
    xlab = "mean expression",
    main = "DE genes for tongue cluster1-background-matching")
  
  # performing topGO
  gene_IDs= rownames(gene_table) # IDs for all skin genes
  in_universe = gene_IDs %in% c(DE_genes_id, back_genes)
  in_selection = gene_IDs %in% DE_genes_id
  
  all_genes = in_selection[in_universe]
  all_genes = factor(as.integer(in_selection[in_universe]))
  names(all_genes) = gene_IDs[in_universe]
  
  # create topGO object
  top_GO_data = new("topGOdata", ontology = "BP", 
                    allGenes = all_genes,
                    nodeSize = 10,
                    annot = annFUN.db,
                    affyLib = "mouse4302.db")
  
  result_top_GO_elim = runTest(top_GO_data, algorithm = "elim",
                               statistic = "Fisher")
  result_top_GO_classic = runTest(top_GO_data, algorithm = "classic",
                                  statistic = "Fisher")
  
  res_top_GO = GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                        Fisher.classic = result_top_GO_classic,
                        orderBy = "Fisher.elim", topNodes = 100)
  genes_top_GO = printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                            chip = "mouse4302.db", geneCutOff = 1000)
  res_top_GO$sig_genes = sapply(genes_top_GO, function(x){
    str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"], ":"), 
          collapse = "")
  })
  
  return(res_top_GO)
}
Genes_Functional_Analysis = function(topGO_res, GO_term, K){
  # A funciton that extract all the selected funtional genes set
  # inputs:   topGO_res (list) = resulted data file from top GO 
  #           GO_term (string) = selected gene funciton
  #           K (numeric) = number of cluster
  # outputs:  genes_fun_all (string vector) = genes from all cluster witin 
  #                                             the selected functional group
  
  # allocate variables 
  genes_fun = list()
  genes_fun_all = c()
  cl_label = c()
  for (k in 1:K){
    genes_fun_temp = filter(topGO_res[[k]], Term == GO_term)$sig_genes %>%
      str_split(":")
    # create cluster labels for the list 
    if (length(genes_fun_temp) > 0){
      cl_label = c(cl_label, k)
    }
    # unroll the list into a vector & remove last empty element
    genes_fun[[k]] = unlist(genes_fun_temp)[-length(unlist(genes_fun_temp))]
    rm(genes_fun_temp) # release dummy variable
  }
  names(genes_fun) = cl_label
  genes_fun_all = unlist2(genes_fun)
  
  return(genes_fun_all)
}
# ----------------
# --input parameters
GO_input = filter(tongue_GG, sample == 1)

tongue_DE = list()
topGO_res = list()
for (k in 1:K){ # loop through each cluster
  tongue_DE[[k]] = filter(GO_input, cluster == k)$probeID
  topGO_res[[k]] = TopGO_Analysis(tongue_DE[[k]], tongue_table, mw_final)
  
}

allGO_terms = c()
for (k in 1:K){
  allGO_terms = c(allGO_terms, topGO_res[[k]]$Term)
}
allGO_terms = unique(allGO_terms)



# ----- only selecte the top 5 GO terms frome each cluster
top5_GO = unique(unlist(map(topGO_res, function(x){x$Term[1:5]})))
# ------------------------------------------------------------------------
# after TopGO analysis, we need evaulate the reults 
#   by performing selected functional genes set analysis
Selected_GO_Plot = function(selected_GO, topGO_result, sample_GG_s1, mw_df, K){
  # -----------------------
  # function plot the selected GO terms in cluster and boxplots
  #     require selected GO terms, DE genes and origanl genes sample
  #
  #     inputs: selected_GO = extracted GO terms based on the users standard
  #             topGO_result = resulted file of topGO analysis
  #             sample_GG_s1 = dataframe of any individual sample with intensity values
  #             mw_df = dataframe of mouse wound model with all tissue sample intensity
  #             K = number of centers of clustering
  #
  #     output: return topGO object that contains the results
  # ----------------------
  # save plots as pdf into selected diretory 
  dir_out = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/mouse_model/temp_plots"
  setwd(dir_out)
  file.remove(list.files()) # delete all files in the directory
  
  for (term in 1:length(selected_GO)){
    funct_name = selected_GO[term]
    temp_genes = Genes_Functional_Analysis(topGO_result, selected_GO[term], K)
    temp_PCA = filter(sample_GG_s1, gene_symbol %in% temp_genes)
    temp_df = filter(mw_df, gene_symbol %in% temp_genes) %>%
      filter(tissue == "tongue" & sample == 1) %>%
      mutate("cluster" = rep(filter(sample_GG_s1, gene_symbol %in% temp_genes)$cluster, 8))
    
    p_cluster = 
      ggplot(temp_PCA, aes(PC1, PC2)) +
      geom_point(alpha = 0.6,aes(colour = cluster)) +
      ggtitle(funct_name) +
      xlab(paste0("PC1, VarExp: ", PCA_eigen_per[1], "%")) +
      ylab(paste0("PC2, VarExp: ", PCA_eigen_per[2], "%")) +
      ylim(-6.5,6.5) +
      xlim(-4, 12) +
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            legend.text = element_text(size = 20))+
      scale_color_manual(values = brewer.pal(5, "Set2"))
    
    p_box = 
      ggplot(temp_df, aes(x = time, y = intensity)) +
      geom_boxplot() +
      facet_wrap(~cluster) +
      ggtitle(funct_name) +
      xlab("time(hrs)") +
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            legend.text = element_text(size = 20))
    
    print(selected_GO[term])
    print(p_cluster)
    ggsave(str_c(funct_name, "--cluster.pdf"))
    
    print(p_box)
    ggsave(str_c(funct_name, "--box.pdf"))
  }
}

Selected_GO_Plot(top5_GO, topGO_res, GO_input, mw_df, K)





