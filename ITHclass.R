library(data.table)

## Classifying patients' ITH according to quartiles PER CANCER TYPE

SPlist <- list()
SPs <- as.data.frame(fread("EXPANDS_results_nozero.tab"))
colnames(SPs) <- c("cancer", "patient", "SPs")
cancers <- unique(SPs$cancer) 

# chol <- SPs[SPs$cancer == "CHOL",]  # doesn't work... WHY?!??!?!?!??!

# classifier loop
for (cancer in cancers){
  SPlist[[cancer]] <- list()
  SPlist[[cancer]] <- SPs[SPs$cancer == cancer,]
  ITHclass <- c()
    for (SP in SPlist[[cancer]][["SPs"]]){                                                                          # write as function, and apply on SPlist
      if (SP < summary(SPlist[[cancer]][["SPs"]])["1st Qu."] | SP == summary(SPlist[[cancer]][["SPs"]])["Min."]){   # "low" if < 1st Quartile
        ITHclass <- append(ITHclass, "low")
        } else if (SP > summary(SPlist[[cancer]][["SPs"]])["3rd Qu."]){                                             # "high" if > 3rd Quartile
          ITHclass <- append(ITHclass, "high")
        } else {                                                                                                    # "moderate" else
          ITHclass <- append(ITHclass, "moderate")
      }
    }
  SPlist[[cancer]] <- cbind(SPlist[[cancer]], ITHclass)
}

# Classifying patients' ITH with classes determined universally for ALL CANCERS