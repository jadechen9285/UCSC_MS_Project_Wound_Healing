# -----------------Mouse Wound Healing Data ----------------
# Modeling time series of mouse wound healing data with transitional probability model.
#
# Jianhong Chen
# 10-22-2019
# --------------------------------------------------------

# Data wrangling
library('tidyverse') # inlcude "tidyr", "readr", "dplyr", "purrr"
library('reshape2') 

# clustering
library("dendextend") # create dendrogram object for hclustering plot 
library("cluster") # clustering functions
library("fpc") # Density-Based Spatial Clustering and application with Noise
library("dbscan")# Density-Based Spatial Clustering and application with Noise

# classifier packages
library(class) # kNN
library("factoextra") # multivariable analysis in R [hkmean]
library(e1071) # naiveBayes 
library(rpart) # decision tree
library(ROCR) # better evalulize the result of the classifier models
library("kernlab") # kernel-based machine learning techniques 

# plotting packages
library("RColorBrewer")
library('ggplot2')
library('plotly') # 3D plotting

# -----------------------------------------------------------------------------------------------------
# intial set up
ls() # list all objects 
gc() # force R to release memory of deleted variables 
rm(list = ls()) # delete all variables in the environment
# -----------------------------------------------------------------------------------------------------


# --load data set 

input_dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Wound_healing/mouse_model/data"
setwd(input_dir)
model_df = read_csv("mw_model.csv")
model_df = model_df[order(model_df$functional),] # sort the functional group in order

time_lab = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")
# compute the difference in intensity bewteen each wounding time vs unwound time points
model_diff = matrix(0, dim(model_df)[1], length(time_lab)-1) # create the 0 matrix with the right size
for (i in 1:(length(time_lab)-1)){
  model_temp = model_df[, -(1:5)] # create dummy stroage df
  model_diff[, i] = unlist(log2(abs(model_temp[, i+1] - model_temp[, 1])+0.05)) #log2 transformation
  rm(model_temp) # remove dummy df
}
colnames(model_diff) = c("d1", "d2", "d3", "d4", "d5", "d6", "d7")
model_diff = cbind(model_df[,1:5], data.frame(model_diff)) %>%
  mutate("functional" = model_df$functional)

model_PCA = prcomp(select(model_diff, d1:d7), scale = T)
#calculate the eigenalues percentage
PCA_eigen_per = percentVar = round(100*model_PCA$sdev^2 / sum(model_PCA$sdev^2), 1)

# --split into training and testing set
SplitTrainTest = function(x_df, split_pcent){
  smp_size = floor(split_pcent*nrow(x_df))
  train_ind <- sample(seq_len(nrow(x_df)), size = smp_size)
  x_train = x_df[train_ind,]
  x_test = x_df[-train_ind, ]
  return(list(x_train, x_test)) 
}
# 1. normal intensity value
# class_input = filter(model_df, functional %in% c("angiogenesis", 
#                                                  "inflammatory response")) %>%
#   select(c(functional, t1:t8)) %>%
#   transform(functional = factor(functional)) 

# 2. time diff. intensity value (use fold change on the dataset)
# class_input = filter(model_diff, functional %in% c("angiogenesis", 
#                                                  "inflammatory response")) %>%
#   select(c(functional, d1:d7)) %>%
#   transform(functional = factor(functional))

# 3. PCA top2 component values
class_input = data.frame(model_PCA$x,
                         "functional" = model_diff$functional) %>%
  filter(functional %in% c("angiogenesis", "inflammatory response")) %>%
  transform(functional = factor(functional))

# class_label = class_input$functional
# class_input = as.data.frame(scale(select(class_input, -functional))) %>%
#   mutate("functional" = class_label)

model_train = SplitTrainTest(class_input, 0.75)[[1]] 
model_test = SplitTrainTest(class_input, 0.75)[[2]]

# -------train svm classifier ---------
ComputeFMeasure = function(prediction, ground_T){
  # compute the F-measure of a binary classifier
  # precision, recall, accuracy, f1-score
  #
  # inputs: prediction = the predicted class label by the model
  #         ground_T = ground true label from the actual file
  #
  # outputs: A list contains the four F-measure metrics. 
  tab = table(prediction, ground_T)
  
  TP = tab[1,1]
  FN = tab[2,1]
  FP = tab[1,2]
  TN = tab[2,2]
  # precision = true pos. / (total predict pos.)
  prec = TP / (TP + FP)
  # recall = true pos./(total actual pos.)
  recall = TP / (TP+FN) 
  
  accuracy = (TP + TN) /(TP + FN + FP + TN)
  # f1 = 2*(prec. * recall) / (prec. + recall)
  f1 = 2*prec*recall/(prec + recall)
  
  result = list(prec, recall, accuracy, f1)
  names(result) = c("precision", "recall", "accuracy", "f1_score")
  return(result)
}

# 1. normal intensity value 
# model_svm = svm(functional ~ t1+t2+t3+t4+t5+t6+t7+t8, model_train,
#                 method="C-classification", kernal="rbfdot",
#                 gamma=0.1, cost=10)

# 2. time diff. intensity value
# model_svm = svm(functional ~ d1+d2+d3+d4+d5+d6+d7, model_train,
#                 method="C-classification", kernal="rbfdot",
#                 gamma=0.1, cost=10)

# 3. PCA top2 component values
model_svm = svm(functional ~. , model_train,
                method="C-classification", kernal="rbfdot",
                gamma=0.1, cost=10)

model_logreg = glm(functional~. , model_train, 
                   family = "binomial" )

logreg_prob = predict(model_logreg, model_test, type = "response")
logreg_pred = ifelse(logreg_prob > 0.5, "angiogenesis", "inflammatory response")
logreg_f1 = ComputeFMeasure(logreg_pred, model_test$functional)

plot(model_svm, model_train, PC1~PC2)

train_pred = predict(model_svm, model_train)
train_f1 = ComputeFMeasure(train_pred, model_train$functional)

test_pred = predict(model_svm, model_test)
test_f1 = ComputeFMeasure(test_pred, model_test$functional)


# --useful parameters
GO_terms = unique(model_df$functional)
time_label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")

J = length(GO_terms) # number of functional group
TT = length(time_label) # time point

# ---compute the global mean
global_mean = c()
for(j in 1:J){
  temp_df = filter(model_df, functional == GO_terms[j])
  global_mean = c(global_mean, mean(unlist(select(temp_df, t1:t8))))
}
global_mean = data.frame(GO_terms, global_mean)

model_s1 = filter(model_df, sample == 1)

# --- convert into Boolean state
ConvertBoolean = function(x_df, threshold_df, time_label){
  
  J = length(threshold_df$GO_terms) # number of functional group
  TT = length(time_label) # time point
  state_res = c()
  # convert into Boolean Values
  for (t in 1:TT){
    for (j in 1:J){
      temp_df = filter(model_s1, functional == global_mean$GO_terms[j]) # only one sample
      local_mean = mean(temp_df[[time_label[t]]])
      if (local_mean > global_mean$global_mean[j]){
        state_res = c(state_res, 1)
      } else {
        state_res = c(state_res, 0)
        
      }
    }
  }
  state_res = matrix(state_res, nrow = 6)
  rownames(state_res) = global_mean$GO_terms
  colnames(state_res) = time_label
  
  return(state_res)
}

s1_res = ConvertBoolean(model_s1, global_mean, time_label)
s2_res = ConvertBoolean(filter(model_df, sample ==2), global_mean, time_label)
s3_res = ConvertBoolean(filter(model_df, sample ==3), global_mean, time_label)






