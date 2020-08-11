#creates a manhattan and qq for MetaXcan results
"%&%" = function(a,b) paste(a,b,sep="")
library(readr)
library(data.table)
library(stringr)
library(dplyr)
library(qqman)
#source("/home/angela/forJournals/GWAS/qqman.r")
options(warn=-1)

load_bp_chrom <- function(pop){
	BP_Chrome <- read.table(paste("/home/egeoffroy/Summary_Statistics/BP_Chrome_files/", pop, "_attempt_newBP.txt", sep = ''), header = T)
	colnames(BP_Chrome) <- c('PROBE_ID', 'gene_name', 'chr', 'BP')
  #BP_Chrome <- transform(BP_Chrome, chr=as.numeric(chr))
  #BP_Chrome %>% mutate(CHR = str_extract(chr, "[0-22]+")
	  BP_Chrome <- transform(BP_Chrome, BP=as.numeric(BP))
	  BP_Chrome$GENE <- gsub("\\..*","",BP_Chrome$PROBE_ID)
	  print(BP_Chrome)
  #BP_Chrome$gene_name <- NULL
}

make_plots <- function(S_Pred_file, directory, title, name){
  print(S_Pred_file)
  pheno1 <- str_split(S_Pred_file, '/')[[1]][2]
  pheno <- str_replace_all(pheno1, '_', ' ')
  pheno <- str_replace_all(pheno, '53', '')
  if(pheno == 'Waist-hip50'){
	pheno <- 'Waist-hip ratio female'
  }
  if(pheno == 'Waist-hip51'){
        pheno <- 'Waist-hip ratio male'
  }
  if(pheno == 'Waist-hip52'){
        pheno <- 'Waist-hip ratio'
  }
  if(pheno == 'White blood'){
	pheno <- 'WBC Count'
  }
  if(pheno == 'PR'){
	pheno <- 'PR interval'
  }
  if(pheno == 'Platelet'){
	pheno <- 'Platelet count'
  }
  if(pheno == 'End renal'){
	pheno <- 'End renal kidney failure'
  }
  print(pheno) 
  if(!str_detect(S_Pred_file, "GWAS")){
  if(!str_detect(S_Pred_file, "sig_genes")){
  if(str_detect(S_Pred_file, "ALL")){
        model <- "ALL"
	BP_Chrome <- load_bp_chrom('ALL')
  }
  if(str_detect(S_Pred_file, "SMulti")){
	model <- "S-MultiXcan"
	BP_Chrome <- load_bp_chrom('ALL')   	
  }
  if(str_detect(S_Pred_file, "AFA")){
        model <- "AFA"
	BP_Chrome <- load_bp_chrom('AFA')   
  }
  if(str_detect(S_Pred_file, "HIS")){
        model <- "HIS"
	BP_Chrome <- load_bp_chrom('HIS')   
  }
  if(str_detect(S_Pred_file, "AFHI")){
        model <- "AFHI"
	BP_Chrome <- load_bp_chrom('AFHI')   
  }
  if(str_detect(S_Pred_file, "CAU")){
	model <- "EUR"
	BP_Chrome <- load_bp_chrom('CAU')   
  } 
   if(str_detect(S_Pred_file, "GTEX")){
	model <- str_split(S_Pred_file, '/')[[1]][4]
	model <- str_remove(model, '.csv')
	BP_Chrome <- load_bp_chrom('ALL')   
   }
   if(str_detect(S_Pred_file, 'SMulti')){
	model <- 'GTEX_v7'
	BP_Chrome <- load_bp_chrom('ALL')
   }
   if(str_detect(S_Pred_file, "HapMap_ALL")){
	model <- "HapMap_ALL"
	BP_Chrome <- load_bp_chrom('ALL')   
   }
   if(str_detect(S_Pred_file, "HapMap_ASN")){
	model <- "HapMap_ASN"
	BP_Chrome <- load_bp_chrom('ALL')   
   }
   if(str_detect(S_Pred_file, "HapMap_AFR")){
	model <- "HapMap_AFR"
	BP_Chrome <- load_bp_chrom('ALL')   
   } 
   print(model)
   S_Pred_file <- read.table(S_Pred_file, header = T,  sep = ',')
 
   S_Pred_file$GENE <- gsub("\\..*","",S_Pred_file$gene) 
  
   S_Pred_file <- S_Pred_file[-c( 9)]
   names(S_Pred_file)[names(S_Pred_file) == 'pvalue'] <- 'P'
   
   GWAS <- merge(S_Pred_file, BP_Chrome, by = c('GENE', 'gene_name'))
   GWAS <- na.omit(GWAS)
   colnames(GWAS)[14] <- "CHR"
   GWAS <- GWAS %>%  #added by Jenny
    transform(CHR = str_replace(CHR, "chr", "")) #added by Jenny
   GWAS<- transform(GWAS, CHR=as.numeric(CHR)) #added by Jenny
   
   threshold <- -log10(0.05/dim(GWAS)[1])
  #threshold <- 1.7e-07

   thresh <- 0.05/dim(GWAS)[1]
   a_threshold <- data.frame(thresh, model, pheno)
   return(a_threshold)
   print(dim(GWAS)[1])
  #signif <- filter(GWAS,-log10(P) > -log10(1.7e-07))
   signif <- filter(GWAS, -log10(P) > -log10(0.05/dim(GWAS)[1]))
   signif <- unique(signif)  
  #if(nrow(signif) > 0 ){

   png(file = paste(directory, '/', model, 'man_qq.png', sep=''))
	par(mfrow=c(2,1))
#  	png(file = paste( directory, "/", model,  "_man_final.png", sep=""))
  	manhattan(GWAS, main = paste('MESA', model, pheno, sep = ' '), suggestiveline = 0, genomewideline = threshold, col = c("blue4", "orange3"))
#  	dev.off()
  
#	png(file = paste(directory, "/", pheno1, model, "_qq_final.png", sep=""))
  	qq(GWAS$P, main = paste('MESA', model, pheno, sep = ' '))
    dev.off()
   #}
    #write.csv(signif, paste(directory, "/", model, "_sig_genes_1e7.csv", sep=""))
}
}
}
