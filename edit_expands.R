#This may be necessary in some systems
#Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk1.8.0_231/")

library(data.table)
library(expands)
#Now try to apply this to TCGA data (eg. ACC)

args = commandArgs(trailingOnly=TRUE)

index = 1
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
  index = args[1]
}

data_links <- fread("sorted_data_links.txt") # see how to avoid using this table and use list of files directly


type <- data_links$Cancer_type[index]
print(type)

# data loop
#for (type in data_links$Cancer_type[index]){
  time_start <- Sys.time()
  cnv_data <- fread(paste("PanCancer_GDC/CNV/", data_links$CNV_file[grep(type, data_links$CNV_file)], sep="")) # see how to avoid using this table and use list of files directly
  # CNV - chr; startpos; endpos; CN_Estimate
  cnv_data$CN_Estimate <- 2^(cnv_data$Segment_Mean)
  colnames(cnv_data)[c(2,3,4)] <- c("chr","startpos","endpos") 
  #ignore sex chromosomes
  cnv_data <- cnv_data[!cnv_data$chr %in% c("X","Y"),]
  cnv_data$Sample <- substr(cnv_data$Sample, 1, 12) # only use case ID
  
  snv_data <- fread(paste("PanCancer_GDC/SNV/", data_links$SNV_file[grep(type, data_links$SNV_file)], sep=""))
  snv_data$Chromosome <- gsub("chr", "", snv_data$Chromosome)
  # Ignore sex chromosomes
  snv_data <- snv_data[!snv_data$Chromosome %in% c("X","Y"),]
  snv_data<-snv_data[,c("Chromosome","Start_Position","End_Position","t_depth","t_ref_count", "t_alt_count","n_depth", "n_ref_count","n_alt_count","case_id","Tumor_Sample_Barcode")]
  # SNV - chr; startpos; endpos; AF_Tumour; PN_B
  snv_data$AF_Tumor <- snv_data$t_alt_count / snv_data$t_depth
  snv_data$PN_B <- 0
  colnames(snv_data)[c(1,2,3)] <- c("chr","startpos","endpos")
  snv_data <- snv_data[,c("chr","startpos","endpos","case_id","Tumor_Sample_Barcode","AF_Tumor","PN_B")]
  snv_data$Tumor_Sample_Barcode <- substr(snv_data$Tumor_Sample_Barcode, 1, 12) # only use case ID
  
  #print(length(unique(snv_data$Tumor_Sample_Barcode)))
  #print(length(unique(cnv_data$Sample)))
  dir.create(paste("results/", type, "_results", sep=""))
  # analysis loop
  matches <- c()
  nr_aliquots <- c()
  for (case in unique(snv_data$Tumor_Sample_Barcode)[1:5]){ 
    if (case %in% unique(cnv_data$Sample)){
      case_snv <- subset(snv_data, Tumor_Sample_Barcode == case)
      case_snv <- subset(case_snv, case_id == unique(case_snv$case_id)[1]) #add line to analyse combination of every case of snv's and cnv's
      #nr_aliquots$snv <- append(nr_aliquots$snv, length(unique(case_snv$case_id)))
      dmcase_snv <- data.matrix(case_snv[,c("chr","startpos","endpos","AF_Tumor","PN_B")])
      case_cnv <- subset(cnv_data, Sample == case)
      case_cnv <- subset(case_cnv, GDC_Aliquot == unique(case_cnv$GDC_Aliquot)[1]) #add line to analyse combination of every case of snv's and cnv's
      dmcase_cnv <- data.matrix(case_cnv[,c("chr","startpos","endpos", "CN_Estimate")])
      #nr_aliquots$cnv <- append(nr_aliquots$cnv, length(unique(case_cnv$GDC_Aliquot)))
      matches <- append(matches, 1)
      #print(head(case_snv))
      dm <- assignQuantityToMutation(dmcase_snv,dmcase_cnv,"CN_Estimate")
      # Let's assume these parameters as they were(?) used for TCGA data
      max_PM=6; maxS=0.7; precision=0.018; plotF=1;
      cfd <- computeCellFrequencyDistributions(dm, max_PM=max_PM, p=precision, nc=6)
      toUseIdx <- which(apply(is.finite(cfd$densities),1,all) )
      SPs <- clusterCellFrequencies(cfd$densities[toUseIdx,], p=precision)

      SPs <- SPs[SPs[,"score"] <= maxS,]; 
      #print(SPs)
        #add line: if file missing
      write.table(SPs,paste("results/",type,"_results/",case,"_SPs.txt", sep=""))
      #aM= assignMutations( dm, SPs, verbose = F)
      #o=plotSPs(aM$dm, snvF,cex=1)
    } else {
      matches <- append(matches, 0)
    }
  }
  cases_run <- print(paste("For cancer type ", type, ", out of ", length(matches), " cases, ", sum(matches), " were matches!", sep=""))
  #print(nr_aliquots)
  time_end <- Sys.time()
#}


# useful functions: append(), list(), write.table(), paste(), rbind(), cbind(), grep(), intersect() to match patients in both lists

