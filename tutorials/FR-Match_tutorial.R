rm(list=ls())
setwd("~/Documents/Lung-matching-2024/Analysis/FR-Match_CellRef_HLCA_binarygenes/") #<-----

###########################################################################################

library(dplyr)
library(tibble)
library(FRmatch)

###########################################################################################
## CellRef --> HLCA using binary genes
###########################################################################################

##########
## data ##
##########

load("../sce_CellRef_query_for_HLCA.rda") #<-----
load("../sce_HLCA_binarygenes.rda") #<-----

sce.E1 <- sce.query <- sce_CellRef_query_for_HLCA #<-----
sce.E2 <- sce.ref <- sce_HLCA_binarygenes #<-----
name.E1 = name.query = "CellRef." #<-----
name.E2 = name.ref = "HLCA." #<-----

## check binary genes
df_binarygenes <- sce.ref@metadata$cluster_marker_info
setdiff(df_binarygenes$markerGene, rownames(sce.query))

binarygenes <- intersect(df_binarygenes$markerGene, rownames(sce.query))
length(binarygenes) #490

#################
## basic plots ##
#################

## cluster size
plot_clusterSize(sce.query, sce.ref, name.E1=name.query, name.E2=name.ref, width = 20, height = 12,
                 filename="plot_clusterSize.pdf")

# ## original distribution
# plot_exprDist(sce.query, sce.ref, markers = binarygenes,
#               name.E1=name.query, name.E2=name.ref,
#               filename="plot_exprDist_orig.pdf")

###################
## barcode plots ##
###################

## create barcode plot folder
folder_name <- "barcode_plots_HLCA"
if(!dir.exists(folder_name)){dir.create(folder_name)}

## cross barcode
query.clusters <- unique(colData(sce.query)$cluster_membership)
for(cl in query.clusters){
  cl_filename <- gsub("/","_", cl)
  plot_cluster_by_markers(sce.ref, sce.query, cluster.name = cl, name.cross = name.query,
                          filename = paste0(folder_name,"/",name.query,cl_filename,".pdf"))
}
## self barcode
ref.clusters <- unique(colData(sce.ref)$cluster_membership)
for(cl in ref.clusters){
  plot_cluster_by_markers(sce.ref, cluster.name = cl, name.self = name.ref,
                          filename = paste0(folder_name,"/",name.ref,cl,".pdf"))
}

###################
## normalization ## 
###################

norm.by = "mean" #<=== FINAL!!!
sce.query.norm <- normalization(sce.query, 
                                norm.by = norm.by
                                ) 
sce.ref.norm <- normalization(sce.ref, 
                              norm.by = norm.by
                              )

# ## post-normalization distribution
# plot_exprDist(sce.query.norm, sce.ref.norm, markers=binarygenes,
#               name.E1=name.query, name.E2=name.ref, xlim=c(0,1), ylim=c(0,20),
#               filename="plot_exprDist_norm.pdf")

##---------------##
## barcode plots ##
##---------------##

## create barcode plot folder
folder_name <- "barcode_plots_HLCA_norm"
if(!dir.exists(folder_name)){dir.create(folder_name)}

## cross barcode
query.clusters <- unique(colData(sce.query)$cluster_membership)
for(cl in query.clusters){
  cl_filename <- gsub("/","_", cl)
  plot_cluster_by_markers(sce.ref.norm, sce.query.norm, cluster.name = cl, name.cross = name.query,
                          filename = paste0(folder_name,"/",name.query,cl_filename,".pdf"))
}
## self barcode
ref.clusters <- unique(colData(sce.ref)$cluster_membership)
for(cl in ref.clusters){
  plot_cluster_by_markers(sce.ref.norm, cluster.name = cl, name.self = name.ref,
                          filename = paste0(folder_name,"/",name.ref,cl,".pdf"))
}

#############
## FRmatch ##
#############

size = 10
iter = 1000

rst_FRmatch_CellRef2HLCA <- FRmatch(sce.query.norm, sce.ref.norm, #prefix = c(name.query, name.ref), 
                            use.cosine=TRUE, filter.size = 10, #<=== filter out Chondrocyte
                            subsamp.size=size, subsamp.iter=iter) ### FINAL: subsamp.size=10, subsamp.iter=1000

##------------------------------------------------------------------------------
## save
save(rst_FRmatch_CellRef2HLCA, file = paste0("rst_FRmatch_CellRef2HLCA_binarygenes_size",size,"_iter",iter,"_norm_",norm.by,".rda"))
##------------------------------------------------------------------------------

sig.level = .05
p.adj.method = "bonferroni"

outputs <- plot_FRmatch(rst_FRmatch_CellRef2HLCA, sig.level = sig.level, p.adj.method = p.adj.method, 
                        main = "CellRef2HLCA", return.value = T, 
                        filename = paste0("plot_FRmatch_CellRef2HLCA_binarygenes_size",size,"_iter",iter,"_norm_",norm.by,"_",p.adj.method,".pdf"))

pmat <- plot_FRmatch(rst_FRmatch_CellRef2HLCA, sig.level = sig.level, type = "padj", p.adj.method = p.adj.method,
                     main = "CellRef2HLCA", return.value = T, 
                     filename = paste0("plot_FRmatch_CellRef2HLCA_binarygenes_padj_size",size,"_iter",iter,"_norm_",norm.by,"_",p.adj.method,".pdf"))

pheatmap::pheatmap(pmat, cluster_rows = F, cluster_cols = F, cellwidth = 10, cellheight=10,  
                   main = "CellRef2HLCA adjusted p-values",
                   filename = paste0("plot_FRmatch_CellRef2HLCA_binarygenes_padjmat_size",size,"_iter",iter,"_norm_",norm.by,"_",p.adj.method,".pdf"))



###########################################################################################
## HLCA --> CellRef using binary genes
###########################################################################################

##########
## data ##
##########

load("../sce_CellRef_binarygenes.rda") #<-----
load("../sce_HLCA_query_for_CellRef.rda") #<-----

sce.E1 <- sce.query <- sce_HLCA_query_for_CellRef #<-----
sce.E2 <- sce.ref <- sce_CellRef_binarygenes #<-----
name.E1 = name.query = "HLCA." #<-----
name.E2 = name.ref = "CellRef." #<-----

## binary genes
df_binarygenes <- sce.ref@metadata$cluster_marker_info
setdiff(df_binarygenes$markerGene, rownames(sce.query))

binarygenes <- intersect(df_binarygenes$markerGene, rownames(sce.query))
length(binarygenes) #432

###################
## barcode plots ##
###################

## create barcode plot folder
folder_name <- "barcode_plots_CellRef"
if(!dir.exists(folder_name)){dir.create(folder_name)}

## cross barcode
query.clusters <- unique(colData(sce.query)$cluster_membership)
for(cl in query.clusters){
  plot_cluster_by_markers(sce.ref, sce.query, cluster.name = cl, name.cross = name.query,
                          filename = paste0(folder_name,"/",name.query,cl,".pdf"))
}
## self barcode
ref.clusters <- unique(colData(sce.ref)$cluster_membership)
for(cl in ref.clusters){
  cl_filename <- gsub("/","_", cl)
  plot_cluster_by_markers(sce.ref, cluster.name = cl, name.self = name.ref,
                          filename = paste0(folder_name,"/",name.ref,cl_filename,".pdf"))
}

###################
## normalization ## 
###################

norm.by = "mean" #<=== FINAL!!!
sce.query.norm <- normalization(sce.query, 
                                norm.by = norm.by
) 
sce.ref.norm <- normalization(sce.ref, 
                              norm.by = norm.by
)

##---------------##
## barcode plots ##
##---------------##

## create barcode plot folder
folder_name <- "barcode_plots_CellRef_norm"
if(!dir.exists(folder_name)){dir.create(folder_name)}

## cross barcode
query.clusters <- unique(colData(sce.query)$cluster_membership)
for(cl in query.clusters){
  plot_cluster_by_markers(sce.ref.norm, sce.query.norm, cluster.name = cl, name.cross = name.query,
                          filename = paste0(folder_name,"/",name.query,cl,".pdf"))
}
## self barcode
ref.clusters <- unique(colData(sce.ref)$cluster_membership)
for(cl in ref.clusters){
  cl_filename <- gsub("/","_", cl)
  plot_cluster_by_markers(sce.ref.norm, cluster.name = cl, name.self = name.ref,
                          filename = paste0(folder_name,"/",name.ref,cl_filename,".pdf"))
}

#############
## FRmatch ##
#############

size = 10
iter = 1000

rst_FRmatch_HLCA2CellRef <- FRmatch(sce.query.norm, sce.ref.norm, #prefix = c(name.query, name.ref), 
                            use.cosine=TRUE, filter.size = 10, #<=== filter out Chondrocyte
                            subsamp.size=size, subsamp.iter=iter) ### FINAL: subsamp.size=10, subsamp.iter=1000

##------------------------------------------------------------------------------
## save
save(rst_FRmatch_HLCA2CellRef, file = paste0("rst_FRmatch_HLCA2CellRef_binarygenes_size",size,"_iter",iter,"_norm_",norm.by,".rda"))
##------------------------------------------------------------------------------

sig.level = .05
p.adj.method = "bonferroni"

outputs <- plot_FRmatch(rst_FRmatch_HLCA2CellRef, sig.level = sig.level, p.adj.method = p.adj.method,
                        main = "HLCA2CellRef", return.value = T, 
                        filename = paste0("plot_FRmatch_HLCA2CellRef_binarygenes_size",size,"_iter",iter,"_norm_",norm.by,"_",p.adj.method,".pdf"))

pmat <- plot_FRmatch(rst_FRmatch_HLCA2CellRef, sig.level = sig.level, type = "padj", p.adj.method = p.adj.method,
                     main = "HLCA2CellRef", return.value = T, 
                     filename = paste0("plot_FRmatch_HLCA2CellRef_binarygenes_padj_size",size,"_iter",iter,"_norm_",norm.by,"_",p.adj.method,".pdf"))

# hist(as.vector(pmat), breaks = 50)
# abline(v=sig.level, col=2)

pheatmap::pheatmap(pmat, cluster_rows = F, cluster_cols = F, cellwidth = 10, cellheight=10, 
                   main = "HLCA2CellRef adjusted p-values",
                   filename = paste0("plot_FRmatch_HLCA2CellRef_binarygenes_padjmat_size",size,"_iter",iter,"_norm_",norm.by,"_",p.adj.method,".pdf"))



###########################################################################################
## FR-Match two-way matching
###########################################################################################

load("rst_FRmatch_CellRef2HLCA_binarygenes_size10_iter1000_norm_mean.rda")
load("rst_FRmatch_HLCA2CellRef_binarygenes_size10_iter1000_norm_mean.rda")

sig.level = .05
p.adj.method = "bonferroni"

plot_bi_FRmatch(rst_FRmatch_CellRef2HLCA, rst_FRmatch_HLCA2CellRef, sig.level = sig.level, p.adj.method = p.adj.method,
                name.E1 = "CellRef.", name.E2 = "HLCA.",
                filename = "plot_FRmatch_CellRef_HLCA_binarygenes_iter1000_twoway.pdf")

plot_bi_FRmatch(rst_FRmatch_CellRef2HLCA, rst_FRmatch_HLCA2CellRef, sig.level = sig.level, p.adj.method = p.adj.method,
                name.E1 = "CellRef.", name.E2 = "HLCA.", two.way.only = TRUE,
                filename = "plot_FRmatch_CellRef_HLCA_binarygenes_iter1000_twoway_only.pdf")

###########################################################################################


