#!/usr/bin/Rscript

#make histogram to compare all of the significant genes found by each model
library(stringr)
library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(data.table)
#source("/home/egeoffroy/Summary_Statistics/scripts/Man_QQ.R")
#all_dir <- list.dirs("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik/", full.names=FALSE)
setwd("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik")
all_files = list.files(pattern = "sig_genes.csv", include.dirs=TRUE, full.names = TRUE, recursive = TRUE)
print(all_files)

list_group <- list()

for(files in all_files){
	if(!str_detect(files, "old")){
		if(str_detect(files, "HIS")){
			#print(files)
			list_group <- append(list_group, files)
		}
	}
}


print(list_group)
new_table <- data.frame(phenotype=character(), significant_hits=integer(), stringsAsFactors=FALSE)

for(i in 1:length(list_group)){
  file <- paste("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik/", str_remove(list_group[i], "./"), sep = "")
  print(file)
  S_Pred_file <- fread(file, header = T,  sep = ',')
  S_Pred_file <- subset(S_Pred_file, select=c("gene_name"))
  pheno <- strsplit(file, "/")
  pheno <- pheno[[1]][7]
  pheno <- gsub("_", ' ', pheno)
  if(pheno == "Waist-hip50"){
  	pheno = "Waist-hip ratio Female"
  }
  if(pheno == "Waist-hip51"){
 	pheno = "Waist-hip ratio Male"
  }
  if(pheno == "Waist-hip52"){
	pheno = "Waist-hip ratio"
  }
  if(pheno == "White blood"){
	pheno = "White Blood"
  }
  if(pheno == 'PR'){
	pheno = 'PR interval'
  }
  pheno <- gsub('[[:digit:]]+', '', pheno)
  print(pheno)
  if(nrow(S_Pred_file) == 0){
	S_Pred_file <- data.frame(phenotype=pheno, significant_hits=0)
  }
  else{
  	S_Pred_file$phenotype <- rep(pheno, nrow(S_Pred_file))
  
  	S_Pred_file$significant_hits <- rep(nrow(S_Pred_file), nrow(S_Pred_file))
  #num_signif <- nrow(S_Pred_file)
  	S_Pred_file <- subset(S_Pred_file, select=c("phenotype", "significant_hits"))
  }
  new_table <- rbind(new_table, S_Pred_file)
  print(new_table)
}
new_table <- new_table[order(significant_hits), ]
#print(new_table)
new_table <- new_table %>% filter(phenotype != "QRS interval")
new_table <- unique(new_table)
print(new_table)
positions <- c("Smoking", "End renal", "Fasting insulin", "Fasting glucose", "Glomerular", "Coffee", "Body mass index", "QT interval", "Hypertension", "Systolic blood", "Diastolic blood", "QRS duration", "Waist-hip ratio", "Waist-hip ratio Female", "Waist-hip ratio Male", "PR interval", "Mean corpuscular hemoglobin", "Hemoglobin", "Diabetes", "Chronic kidney", "Platelet", "C-reactive", "Triglyceride", "Total cholesterol", "LDL cholesterol", "HDL cholesterol", "Height", "White Blood")
print(positions)
ggplot(new_table, aes(x = phenotype, y = significant_hits, fill=phenotype)) + coord_flip()+ geom_bar(stat = "identity") + labs(title="Wojcik/HIS MESA Model", 
                                                                                                                  x="Phenotype", y = "Number of Hits")+theme(axis.text.y = element_text(color = "grey20", size = 10, hjust = .5, vjust = .5, face = "plain")) + scale_x_discrete(limits = positions)#+ theme(legend.position = "none")#+ scale_y_continuous(name="", limits=c(0, 45))
ggsave("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik/HIS_MESA_Wojcik.png")

