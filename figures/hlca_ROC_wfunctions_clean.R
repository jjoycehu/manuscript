rm(list=ls())
# setwd("/Users/joycehu/Dev/hlca") #<-----


###########################################################################################

library(data.table)
library(tidyverse)
library(ROCR)

###########################################################################################

all_clusters <- fread("data/cluster_order_HLCA.csv", sep=",")$cluster_order #<-----
all_clusters_unassigned <- c(all_clusters, "unassigned")
n <- length(all_clusters)

df_temp <- data.frame(query = rep(all_clusters, each=n+1), prediction = rep(all_clusters_unassigned, times=n))
## check
nrow(df_temp) == n*(n+1) #TRUE



#####

get_df_summary<- function(df){
  # double groupby
  df_summary <- df %>% group_by(hlca_true, hlca_pred) %>%
    summarise(mean.score=mean(score)) %>% ungroup() %>% 
    mutate(truth = ifelse(hlca_true==hlca_pred,1,0))

  # df_summary <- df_temp %>%
  #   left_join(df_summary, by = c("query"="hlca_true", "prediction"="hlca_pred")) %>% #<-----
  # mutate(mean.score = replace_na(mean.score, 0)) %>% mutate(truth=ifelse(query==prediction,1,0))
  
  return(df_summary)
}

get_pred <- function(df_summary){
  return(prediction(df_summary$mean.score, df_summary$truth))
}

auc_func <- function(pred){
  return(performance(pred, measure = "auc")@y.values[[1]])
}

plot_roc <- function(df, title){
  
  df_summary <- get_df_summary(df)
  
  pred <- get_pred(df_summary)
  
  perf <- performance(pred, "tpr", "fpr")
  
  
  plot(perf,
       colorize = TRUE,
       print.cutoffs.at = c(0.8),
       text.adj=c(-0.5,0.5), # move position of cutoff label
       lwd = 3,
       main = title) #<-----
  abline(a = 0, b = 1, col="grey", lwd =3, lty = 3)
  
  ## AUC
  auc <- auc_func(pred)
  legend("bottomright", legend = paste("AUC =", round(auc,3)), bty = "n") 
  return(auc)
}
#### 

##############
## FR-Match ##
##############
fr <- fread("data/frmatch_hlca.csv")
auc_frmatch <- plot_roc(fr, "FR-Match ROC Curve")


#############
## Azimuth ##
#############
az <- fread("data/azimuth_hlca.csv")
auc_azimuth <- plot_roc(az, "Azimuth ROC Curve")

################
## CellTypist ##
################
ct <- fread("data/celltypist_hlca.csv")
auc_celltypist <- plot_roc(ct, "Celltypist ROC Curve")


##############
## scArches ##
##############
sc <- fread("data/scarches_hlca.csv")
auc_scarches <- plot_roc(sc, 'scArches ROC Curve')


##############
## singleR ##
##############
sr <- fread("data/singleR_hlca.csv") %>% replace_na(list(hlca_pred="unassigned"))
auc_singleR <- plot_roc(sr, 'singleR ROC Curve')

##############
## singleR ##
##############
sp <- fread("data/scPred_hlca.csv") %>% replace_na(list(hlca_pred="unassigned"))
auc_scPred <- plot_roc(sp, 'scPred ROC Curve')

###########################################################################################
## ALL methods
###########################################################################################

pdf("ROC plot.pdf", width = 6, height = 6)
par(mar = c(5.1, 4.1, 2.1, 2.1), mgp = c(2, .7, 0), oma = c(1,1,1,0))
# RColorBrewer::display.brewer.pal(n = 4, name = 'Dark2') # <------
# mycolors <- RColorBrewer::brewer.pal(n = 4, name = "Dark2")[c(4,2,1,3)]
mycolors <- c('#A4CE95', '#F4EDCC','#FF7070', '#965F8A', '#4AC6B7', '#F3AE4B')

## Azimuth
az_summary <- get_df_summary(az)
az_pred <- get_pred(az_summary) 

## ROC curve
az_perf <- performance(az_pred, "tpr", "fpr")
plot(az_perf,
     col = mycolors[3],
     lwd = 3,
     main = "ROC Curves")
abline(a = 0, b = 1, col="grey", lwd =3, lty = 3)


## CellTypist
ct_summary <- get_df_summary(ct)
ct_pred <- get_pred(ct_summary)
## ROC curve
ct_perf <- performance(ct_pred, "tpr", "fpr")
plot(ct_perf,
     col = mycolors[4],
     lwd = 3, 
     add=TRUE)


## FR-Match
fr_summary <- get_df_summary(fr)
fr_pred <- get_pred(fr_summary)

## ROC curve
fr_perf <- performance(fr_pred, "tpr", "fpr")
plot(fr_perf,
     col = mycolors[6],
     lwd = 3,
     add = TRUE)


## scArches
sc_summary <- get_df_summary(sc)
sc_pred <- get_pred(sc_summary)

## ROC curve
sc_perf <- performance(sc_pred, "tpr", "fpr")
plot(sc_perf,
     col = mycolors[5],
     lwd = 3,
     add = TRUE)


## singleR
sr_summary <- get_df_summary(sr)
sr_pred <- get_pred(sr_summary)

## ROC curve
sr_perf <- performance(sr_pred, "tpr", "fpr")
plot(sr_perf,
     col = mycolors[2],
     lwd = 3,
     add = TRUE)

## scPred
sp_summary <- get_df_summary(sp)
sp_pred <- get_pred(sp_summary)

## ROC curve
sp_perf <- performance(sp_pred, "tpr", "fpr")
plot(sp_perf,
     col = mycolors[1],
     lwd = 3,
     add = TRUE)


legend("bottomright", 
       legend = paste(c("scPred", "singleR", "Azimuth", "CellTypist", "scArches", "FR-Match"), "=", 
                      round(c(auc_scPred, auc_singleR, auc_azimuth, auc_celltypist, auc_scarches, auc_frmatch),3)), #<-----
       title = "AUC statistic", col = mycolors[1:6], lwd =3, bty = "n")
dev.off()
