#Run through all of the GTEx files
#library(lubridate)
library(dplyr)
source("/home/egeoffroy/Summary_Statistics/scripts/Man_QQ.R")
all_dir <- list.dirs("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik/", full.names=FALSE)
#print(all_dir)
#for(dirs in all_dir){
	#if(str_detect(dirs, "QT_interval")){
	#print(paste("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik", dirs, sep='/'))
	setwd(paste("/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik"))
	all_files = list.files(pattern = ".csv", full.names = TRUE, recursive = TRUE) #recursive is usually true except for plotting SMultiXcan results
#	print(all_files)	
thresholds <- data.frame()	
Sum_Stat_Name <- "Wojcik" #change for each GWAS Summary Statistic
for (file in all_files) {
 #print(file)
 if(str_detect(file, 'MESA') && !str_detect(file, 'old') && !str_detect(file, 'sig_genes') && !str_detect(file, 'GWAS') && (str_detect(file, 'AFHI2') || str_detect(file, 'CAU'))){  
	print(file)
	directory <- paste(normalizePath(dirname(file)), '/', sep = '')
 	print(directory) 
 # read in the csv
#	thresholds <- data.frame()
        thresh <-  make_plots(file, directory, sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file)), sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file)))
  	make_plots(file, directory, sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file)), sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file)))
	print(thresh)
	thresholds <- rbind(thresholds, thresh)
	write.csv(thresholds, 'Wojcik_thresholds_man.csv', quote = F, row.names=F)
  }
}

