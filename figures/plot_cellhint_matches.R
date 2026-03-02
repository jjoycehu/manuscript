rm(list=ls())
setwd("/gpfs/gsfs12/users/zhangy71/Kidney-2025/CellHint/")

library(dplyr)
library(pheatmap)
library(RColorBrewer)

df <- data.table::fread("cellhint_matches.csv") %>% as.data.frame()



load("/gpfs/gsfs12/users/zhangy71/Kidney-2025/FR-Match/prep_data/sce_ref_Lake_NSForest.rda")
load("/gpfs/gsfs12/users/zhangy71/Kidney-2025/FR-Match/prep_data/sce_ref_mBDRC_NSForest.rda")

clusters_lake <- sce_ref_Lake_NSForest@metadata$cluster_order
clusters_mBDRC <- sce_ref_mBDRC_NSForest@metadata$cluster_order

mat <- matrix(0, nrow = length(clusters_lake), ncol = length(clusters_mBDRC))
colnames(mat) <- clusters_mBDRC
rownames(mat) <- clusters_lake

for (i in 1:nrow(df)){
  mat[df[i,2],df[i,1]] <- 2
}

colnames(mat) <- paste0("mBDRC.", colnames(mat))
rownames(mat) <- paste0("Lake.", rownames(mat))

unassigned <- 2*as.numeric(colSums(mat)==0)
mat <- rbind(mat, unassigned)

mat <- FRmatch:::reorder(mat)
pheatmap(mat,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(1,1,7)]))(3),
         breaks=seq(0,2,length.out=3),
         legend_breaks=c(0,2),
         legend_labels=c("No match", "Match"),
         cluster_rows=F, cluster_cols=F,
         gaps_row=nrow(mat)-1,
         cellwidth=10, cellheight=10,
         main="CellHint harmonization",
         filename="plot_cellhint_matches.pdf")
