########################################################################
# mouse_wound Affy data quick pre-processing without any plots
# skin vs tongue in 8 time point 
# author: Jianhong Chen
# date: 09/25/2019

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
library(GEOquery)
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

gse <- getGEO("GSE460", GSEMatrix = TRUE)

filePaths = getGEOSuppFiles("GSE460") # No raw data from this experiment
anno_AE = getAE("E-GEOD-460", path = raw_data_dir, type ="raw") # anno_AE download data from fttp link


# extra information from the raw data to construct a dataframe
# I. obtain info from SDRF file 
sdrf_location = file.path(raw_data_dir, "E-GEOD-460.sdrf.txt") #similar to "paste" (but faster in R)
SDRF = read.delim(sdrf_location)

# create a new column names the CEL files label
data_file_name = SDRF$Source.Name %>%
  sapply(function(x){paste0(x, ".CEL")})
SDRF = mutate(SDRF, "Array.Data.File" = data_file_name)
rownames(SDRF) = SDRF$Array.Data.File # replace generic machine label of the files

SDRF = AnnotatedDataFrame(SDRF) # create a biobased dataframe object
# II. create Expression Set object of array data
raw_data = oligo::read.celfiles(filenames = file.path(raw_data_dir, SDRF$Array.Data.File),
                                verbose = F, phenoData = SDRF)









