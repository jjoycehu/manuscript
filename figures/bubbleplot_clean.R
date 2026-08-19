library(ggplot2)
library(dplyr)
library(viridis)
rm(list = ls())

setwd("/Users/joycehu/Dev/hlca")


## hlca order
order_df <- read.csv('bubble_plots/cluster_order_HLCA.csv')
order_hlca <- order_df$cluster_order

## cellref order
cellref_order_df <- read.csv('materials/cellref_cluster_order_to_plot.csv')
order_cellref <- cellref_order_df$cluster_order

# hlca azimuth 
{
  hlca <- read.csv('azimuth/hlca/azimuth_hlca.csv')
  bubble_prop(hlca, 'hlca_pred', 'hlca_true', order_hlca, order_hlca, 'HLCA Annotations by Proportion Using Azimuth')
}

# hlca celltypist
{hlca <- read.csv('celltypist/hlca/celltypist_hlca.csv')
  bubble_prop(hlca, 'hlca_pred', 'hlca_true', order_hlca, order_hlca, 'HLCA Annotations by Proportion Using Celltypist')
}

# hlca scarches
{
  hlca <- read.csv('scArches/hlca/scarches_hlca.csv')
  bubble_prop(hlca, 'hlca_pred', 'hlca_true', order_hlca, order_hlca, 'HLCA Annotations by Proportion Using scArches')
}

# hlca fr match 
{
  hlca <- read.csv('miscellaneous/frmatch_hlca.csv')
  bubble_prop(hlca, 'hlca_pred', 'hlca_true', order_hlca, order_hlca, 'HLCA Annotations by Proportion Using FR-Match')
}



bubble_prop <- function(df, pred_label, true_label, pred_order, true_order, title, x_name='true_label', y_name='predicted_label'){
  # table of matches
  tab.match <- table(df[[pred_label]], df[[true_label]]) # true labels in columns, prediction in rows
  tab.match.prop <- sweep(tab.match,2,colSums(tab.match),"/") #column sums should be 1
  
  ## order predicted type in rows
  oo <- match(c(pred_order,"unassigned"), rownames(tab.match.prop))
  oo.names <- rownames(tab.match.prop)[oo] %>% na.omit() #some ref clusters may not have matched query cells
  tab.match.prop <- tab.match.prop[oo.names,]
  
  ## order columns 
  cc <- match(c(true_order,"unassigned"), colnames(tab.match.prop))
  cc.names <- colnames(tab.match.prop)[cc] %>% na.omit() 
  tab.match.prop <- tab.match.prop[,cc.names]
  
  # select cols 
  gt <- tab.match.prop %>% as.data.frame() %>%
    select(predicted_label=Var1, true_label=Var2, Prop=Freq)%>%
    mutate(true_label=factor(true_label, levels = c(true_order, "unassigned"))) %>%
    mutate(predicted_label=factor(predicted_label, levels=rev(c(pred_order, "unassigned"))))
  
  ggplot(gt, aes(x=true_label, y=predicted_label, size=Prop, fill=Prop)) +
    geom_point(alpha=0.7, shape=21, color="black") +
    scale_size_continuous(range = c(0, 10), limits = c(0, 1)) +
    scale_fill_viridis(option="D", guide = "legend", limits = c(0, 1)) +
    scale_y_discrete(label=rev(pred_order), limits=rev(pred_order), drop=FALSE) + #show all ref clusters even if no match
    scale_x_discrete(label=true_order, limits=true_order, drop=FALSE) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = title, x = x_name, y = y_name)
  
}



