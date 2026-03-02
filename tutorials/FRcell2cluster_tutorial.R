rm(list=ls())
setwd("~/Documents/Lung-matching-2024/Analysis/FR-Match_CellRef_HLCA_binarygenes/FRcell2cluster_CellRef_HLCA_binarygenes/") #<-----

###########################################################################################

library(dplyr)
library(tibble)
library(magrittr)
library(FRmatch)

###########################################################################################
## CellRef ---cell2cluster--> HLCA using binary genes
###########################################################################################

load("../../sce_CellRef_query_for_HLCA.rda") #<-----

## my ten-fold split
df_ind <- read.csv("cellref_tenfold_idx.csv") #<-----

## take care of small clusters, except Chondrocyte
ind_small <- which(sce_CellRef_query_for_HLCA@colData$cluster_membership %in% c("Megakaryocyte/Platelet", "Innate lymphoid cell"))

myind_lst <- vector(length = 10, mode = "list")
for (i in 1:10){
  ind <- df_ind[,i]
  myind <- setdiff(ind, ind_small)
  myind_lst[[i]] <- myind
}
myind_lst[[1]] <- c(myind_lst[[1]], ind_small) #<=====USE THIS INDEX LIST!!!

## check
sum(sapply(myind_lst, FUN = function(z) length(z))) #347970 = total number of CellRef cells


###############################################################################################################
## run 10 folds
###############################################################################################################

## load data
load("../../sce_CellRef_query_for_HLCA.rda") #<-----
load("../../sce_HLCA_binarygenes.rda") #<-----

## set up data
sce.query <- sce_CellRef_query_for_HLCA #<-----
sce.ref <- sce_HLCA_binarygenes #<-----

## normalization on the whole dataset
norm.by = "mean" #<=== FINAL!!!
sce.query.norm <- normalization(sce.query, 
                                norm.by = norm.by
) 
sce.ref.norm <- normalization(sce.ref, 
                              norm.by = norm.by
)

## set up ten-fold
for(i in 1:10){
  ind <- myind_lst[[i]] #<-----fold i <===USE MY INDEX LIST!!!
  sce.query.norm.i <- sce.query.norm[,ind] #<-----

  size = 10
  iter = 2000
  k = 5
  rst_cell2cluster_CellRef2HLCA <- FRmatch_cell2cluster(sce.query.norm.i, sce.ref.norm,
                                                        prefix = c("",""), #prefix = c("CellRef.","HLCA."),
                                                        subsamp.size=size, subsamp.iter=iter,
                                                        subsamp.iter.custom=TRUE, subsamp.iter.custom.k = k)
  
  # ## load pre-generated results for plots
  # load(paste0("rst_FRcell2cluster_CellRef2HLCA_size",size,"_iter",iter,"custom",k,"_",p.adj.method,"_fold",i,".rda"))
  
  ## order query clusters for plots
  oo <- read.csv("cellref_cluster_order_to_plot.csv")$cluster_order
  rst_cell2cluster_CellRef2HLCA$cell2cluster %<>% mutate(query.cluster=factor(query.cluster, levels = oo))
  
  ## plots
  sig.level = .05
  p.adj.method = "BY"
  outputs <- plot_FRmatch_cell2cluster(rst_cell2cluster_CellRef2HLCA, sig.level = sig.level, p.adj.method = p.adj.method,
                                       return.value = T, reorder = F, width = 15, height=15,
                                       filename = paste0("plot_FRcell2cluster_size",size,"_iter",iter,"custom",k,"_",p.adj.method,"_fold",i,".pdf"))
  df_cell2cluster_adj <- outputs$cell2cluster.adj
  
  plot_FRmatch_cell2cluster(rst_cell2cluster_CellRef2HLCA, sig.level = sig.level, type = "padj", p.adj.method = p.adj.method, reorder = F, 
                            filename = paste0("plot_FRcell2cluster_padj_size",size,"_iter",iter,"custom",k,"_",p.adj.method,"_fold",i,".pdf"))
  
  plot_FRmatch_cell2cluster(rst_cell2cluster_CellRef2HLCA, sig.level = sig.level, type = "match.all", 
                            reorder = F, width = 15, height=15,
                            filename = paste0("plot_FRcell2cluster_matchall_size",size,"_iter",iter,"custom",k,"_fold",i,".pdf"))
  ## save
  save(rst_cell2cluster_CellRef2HLCA, df_cell2cluster_adj,
       file = paste0("rst_FRcell2cluster_CellRef2HLCA_size",size,"_iter",iter,"custom",k,"_",p.adj.method,"_fold",i,".rda"))
}


###############################################################################################################
## combine ten folds and save matching results for each cell
###############################################################################################################

df_all <- as.data.frame(c())

for(i in 1:10){
  load(paste0("rst_FRcell2cluster_CellRef2HLCA_size",size,"_iter",iter,"custom",k,"_",p.adj.method,"_fold",i,".rda"))
  df_i <- df_cell2cluster_adj %>% mutate(fold = i)
  cat(nrow(df_i), "\n")
  df_all <- rbind(df_all, df_i)
}
## save
write.csv(df_all, file="FR-Match_cell2cluster_CellRef.csv", row.names = F)



