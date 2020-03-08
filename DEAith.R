library(edgeR)
library(dplyr)
library(data.table)

exp_data <- list()    #DEA results
expr_model <- list()  #expression for model fit
#min_samples <- 4
#SPs <- as.data.frame(fread("EXPANDS_results_nozero.tab"))       # requires this results table, SPs per case
#colnames(SPs) <- c("cancer", "patient", "SPs")
#cancers <- unique(SPs$cancer) 
          # NB! LAML and PCPG have few samples

source("ITHclass.R") # run my ITH classifier
cancers <- cancers[!cancers %in% c("LAML", "PCPG")]  # LAML, PCPG removed due to low samples 

for (cancer in cancers){          
  print(cancer)
  
  expr <- read.table(paste("Pancancer_GDC/Transcriptome/TCGA-",cancer,"_htseq_counts.tab",
                           sep=""),
                     sep="\t", header=T, quote = "", comment.char = '', stringsAsFactors = F)          # RNA-seq data, pats as colnames
  expr <- expr[-c(1:5),]                                                                               # remove first 5 rows of quality data
  rownames(expr) <- expr[,1]      # set genes, 1st col, as rownames
  expr <- expr[-c(1)]             # remove 1st col
  expr_model[[cancer]] <- list()
  expr_model[[cancer]] <- expr
 
  ITH_set <- SPlist[[cancer]]
  ITH_set$patient <- gsub("-", ".", ITH_set$patient)            # sample names between RNA and ITH sets have different sep.
  RNA_set <- data.frame("samples" = colnames(expr)[grep(".01A", colnames(expr))],
                        "patient" = substr(colnames(expr)[grep(".01A", colnames(expr))], 1, 12), stringsAsFactors = F)  # only tumor samples "01A"
  
  ITH_low <- ITH_set[ITH_set$ITHclass == "low",]
  ITH_high<- ITH_set[ITH_set$ITHclass == "high",]
  
  RNA_low <- RNA_set[RNA_set$patient %in% ITH_low$patient ,]              # RNA samples with low ITH
  RNA_low$ITHclass <- "low"
  RNA_high<- RNA_set[RNA_set$patient %in% ITH_high$patient ,]             # RNA samples with high ITH
  RNA_high$ITHclass <- "high"
    
#  if(length(RNA_high) > min_samples |
#     length(RNA_low ) > min_samples){                  
    
    exp_data[[cancer]] <- list()
    
    expr_matched <- as.data.frame(expr[,c(RNA_low$samples,                # here, the above, stringsAsFactors = F, is crucial
                                          RNA_high$samples)])             # read counts: low ITH vs high ITH
    colData <- data.frame(rbind(RNA_low, RNA_high), stringsAsFactors = F)
    rownames(colData) <- colData$samples
    colData$samples <- NULL
    colData$ITHclass <- factor(colData$ITHclass, levels = c("low","high"))
    
    ## edgeR
    y <- DGEList(counts=expr_matched)
    #Maybe choose another critera for filtering low expressed genes... since we have so many samples
    keep <- filterByExpr(y)                         # Genes that have sufficiently large counts to be retained in a statistical analysis
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)                         # Calculate normalization factors to scale the raw library sizes.
    cpm <- cpm(y, normalized.lib.sizes = TRUE)
    # With so many samples, we can use other methods to make the statistical test...
    # Voom, or a simple wilcoxon!

    ## limma    
    design <- model.matrix(~ ITHclass, data=colData)    # ??????????????????????????
    y <- voom(y, design, plot = T)
    fit <- lmFit(y, design)                         # fit where variance of low exp genes is shrunk towards fit with increasing penalty
    tmp <- eBayes(fit)                              # DEA
    top.table <- topTable(tmp, sort.by = "P", n = Inf, coef="ITHclasshigh")   # see limma guide 
    
    ## summary
    exp_data[[cancer]][["cpm"]] <- cpm
    exp_data[[cancer]][["results"]] <- top.table
    exp_data[[cancer]][["metadata"]] <- colData
    
    gc()
#  }
  
}
