library(dbplyr)
library(dplyr)
library(data.table)
# DEA results, low vs. high ITH, from DEAith.R are in variable exp_data
# collect number of genes with adjP > 0.05 in each cancer

signi_exp_data <- list()
cancers <- cancers[!cancers %in% c("LAML", "PCPG")]  # LAML, PCPG removed due to low samples 
for (cancer in cancers){
  if (sum(exp_data[[cancer]][["results"]]["adj.P.Val"]<0.05)>10){     # filter cancers < 10 significant DE genes
    signi_exp_data[[cancer]] <- list()
    signi_exp_data[[cancer]] <- exp_data[[cancer]][["results"]][exp_data[[cancer]][["results"]]["adj.P.Val"]<0.05,]
  }
}
# 55 common DE genes from top 8 cancers with most abundant samples:
top8 <- Reduce(intersect, list(rownames(signi_exp_data$LUAD), 
                               rownames(signi_exp_data$STAD), 
                               rownames(signi_exp_data$LUSC), 
                               rownames(signi_exp_data$BRCA), 
                               rownames(signi_exp_data$HNSC), 
                               rownames(signi_exp_data$COAD), 
                               rownames(signi_exp_data$BLCA), 
                               rownames(signi_exp_data$PRAD)))

