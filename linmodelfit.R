library(ISLR)
library(glmnet)
library(pheatmap)
## build the matrix to perform linear model fit

# per cancer: cols: patient; genes; ITH

# use SPlist from ITHclass.R

#source("DEAith.R")        # cancers ITH list: SPlist
#source("DE_genes.R")      # selected DE genes: top8

# with genes from top8 as identifiers, select and normalize expression for these

head(SPlist$LUSC)         # SPlist already has most of the structure needed. only genes expression missing
DEGs <- list()
normalized_counts <- list()
cancers <- cancers[!cancers %in% c("LAML", "PCPG")]  # LAML, PCPG removed due to low samples 
for (cancer in cancers) {
  DEGs[[cancer]] <- as.data.frame(expr_model[[cancer]][top8,])  # expr_model is reads of all genes and all samples; top8 is the ones that I select
  DEGs[[cancer]] <- DEGs[[cancer]][,grep(".01A", colnames(DEGs[[cancer]]))]
  z <- DGEList(counts=DEGs[[cancer]])
  keep <- filterByExpr(z)                         # Genes that have sufficiently large counts to be retained in a statistical analysis
  z <- z[keep, , keep.lib.sizes=FALSE]
  z <- calcNormFactors(z)                         # Calculate normalization factors to scale the raw library sizes.

  normalized_counts[[cancer]] <- cpm(z, normalized.lib.sizes = TRUE)      # counts per million
  colnames(normalized_counts[[cancer]]) <- substr(colnames(normalized_counts[[cancer]]), 1, 12)
  normalized_counts[[cancer]] <- as.data.frame(t(normalized_counts[[cancer]]))
  normalized_counts[[cancer]] <- cbind(row.names(normalized_counts[[cancer]]), normalized_counts[[cancer]])
  colnames(normalized_counts[[cancer]])[1] <- "patient"
}
model_matrix <- list()  # prepare matrix for model
lm_matrix <- list()
for (cancer in cancers){
  SPlist[[cancer]][["patient"]] <- gsub("-", ".", SPlist[[cancer]][["patient"]])
  model_matrix[[cancer]] <- merge(SPlist[[cancer]], normalized_counts[[cancer]], 
                                  by = "patient")
  lm_matrix[[cancer]] <- model_matrix[[cancer]][, -c(1,2,4)]       # remove excess columns
  lm_matrix[[cancer]][,-1] <- scale(lm_matrix[[cancer]][,-1])      # scales with mean at 0
}

## linear model fit

fitcols <- lm_matrix$LUAD
fittest <- lm(SPs~ENSG00000093009.8, data=fitcols)

## apply lasso to find signals in the data and come up wth interpretable subsets among the features presentad to us

x <- model.matrix(SPs~.-1, data=fitcols)
y <- fitcols$SPs

fit.lasso <- glmnet(x,y, alpha=1)          # glmnet does a shrinkage on the variables and performs subset selection by shrinking the coefficients towards zero
plot(fit.lasso, xvar = "lambda", label = TRUE)
print(fit.lasso)
#plot(fit.lasso, xvar = "dev", label = TRUE)
cv.lasso <- cv.glmnet(x,y,alpha=1)         # cv.glmnet does cross validation to pick a appropriate lambda
plot(cv.lasso)
coefs <- coef(cv.lasso)                    # coef() gives coeficients of "best" model (model with least variables within 1 std of model with minimum MSE)

## pick the genes selected by model fit with lasso for each cancer and compare


model_res <- list()
model_coefs <- list()
for (cancer in names(lm_matrix)){
  data <- lm_matrix[[cancer]]
  x <- model.matrix(SPs~.-1, data=data)
  y <- data$SPs
                                           # glmnet does a shrinkage on the variables and performs subset selection by shrinking the coefficients towards zero
  cv.lasso <- cv.glmnet(x,y)               # cv.glmnet does cross validation to pick an appropriate lambda (shrinkage factor)
  model_res[[cancer]] <- glmnet(x,y, lambda = cv.lasso$lambda.1se)    # lasso summary
  tmp_coefs <- coef(cv.lasso)
  coefs <- data.frame(name = tmp_coefs@Dimnames[[1]][tmp_coefs@i + 1], coefficient = tmp_coefs@x)
  if (nrow(coefs) > 1){
  model_coefs[[cancer]] <- tmp_coefs      # coeficients for genes from shirkage, from this build COEF MATRIX
  }
  coef_matrix <- matrix(
    nrow = length(model_coefs),
    ncol = length(gene_group),
    dimnames = list(names(model_matrix), top8)
    
}

## heatmap for plot?



