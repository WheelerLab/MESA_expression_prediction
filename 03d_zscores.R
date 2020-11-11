library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(wesanderson)

all_files <- list.files('/home/egeoffroy/Summary_Statistics/Summary_Stats/Wojcik', pattern = 'MESA', recursive = T, full.names=T)
#print(all_files)
output <- data.frame()
for(file in all_files){
        if(!str_detect(file, 'old') && !str_detect(file, 'GWAS') && !str_detect(file, 'sig_genes') && !str_detect(file, 'HapMap')){
        if(str_detect(file, '.csv')){
        if(str_detect(file, 'CAU') || str_detect(file, 'AFHI2') || str_detect(file, 'ALL2')){
        print(file)
        file1 <- fread(file, sep = ',', header = T, stringsAsFactors=F)
        model <- str_split(file, '/')[[1]][8]
        model <- str_remove(model, '2')
        model <- str_remove(model, 'MESA_')
        model <- str_remove(model, '.csv')
        phenotype <- str_split(file, '/')[[1]][7]
        phenotype <- str_replace_all(phenotype, '_', ' ')
        file1$phenotype <- rep(phenotype, nrow(file1))
        file1$model <- rep(model, nrow(file1))
        output <- rbind(output, file1)
}
}
}
}

print(output)
output <- output %>% select(gene, gene_name, model, phenotype, zscore, effect_size)
output$phenotype <- str_replace_all(output$phenotype, 'C-reactive', 'C-reactive protein levels')
output$phenotype <- str_replace_all(output$phenotype, 'White blood', 'WBC count')
ALL <- output %>% filter(model == 'ALL')
AFHI <- output %>% filter(model == 'AFHI')
CAU <- output %>% filter(model == 'CAU')

total <- merge(CAU, AFHI, by =c('gene', 'phenotype'))
#total <- merge(total, CAU, by = c('gene', 'phenotype'))
print(head(total))

pearson1 <- cor.test(total$zscore.x,total$zscore.y, method="pearson")
print(pearson1)
p <- ggplot(total, aes(x=zscore.x, y=zscore.y)) + geom_point()+ labs( x="EUR Gene Z-scores", y = "AFHI Gene Z-scores")  + annotate(geom="text", size = 6, x=-10, y =15, label=paste("R = ", round(pearson1$estimate, digits = 3), sep = ''))+geom_smooth(method = 'lm', fullrange=TRUE, aes(group=1)) +geom_density_2d() +theme_minimal()+ theme_bw(20)+ylim(-20,20) + xlim(-20, 20) + ggsave('AFHI_EUR_zscore.tiff', device = 'tiff')
pearson <- cor.test(total$effect_size.x,total$effect_size.y, method="pearson")
print(pearson)
p <- ggplot(total, aes(x=effect_size.x, y=effect_size.y)) + geom_point()+ labs( x="EUR Gene Effect-size", y = "AFHI Gene Effect-size")  + annotate(geom="text", x=-10, y =15, label=paste("R = ", round(pearson$estimate, digits = 3), sep = ''))+theme_bw(20)+geom_smooth(method = 'lm', fullrange=TRUE, aes(group=1)) + theme_minimal() + geom_density_2d() + ylim(-20,20) + xlim(-20, 20) + ggsave('AFHI_EUR_betas.png')
write.csv(total, 'zscores_total.csv', quote = F, row.names=F)

SPred <- fread('/home/egeoffroy/coloc/SPred_total_all_pheno.csv', header = T, sep = ',', stringsAsFactors = F)
print(head(SPred))
SPred$Phenotype <- str_replace_all(SPred$Phenotype, 'White_blood', 'WBC count')
SPred$Phenotype <- str_replace_all(SPred$Phenotype, '53', '')
SPred$Phenotype <- str_replace_all(SPred$Phenotype, 'C-reactive', 'C-reactive protein levels')
SPred$Phenotype <- str_replace_all(SPred$Phenotype, '_', ' ')

ALL2 <- SPred %>% filter(Model == 'ALL')
AFHI2 <- SPred %>% filter(Model == 'AFHI')
CAU2 <- SPred %>% filter(Model == 'CAU')
total <- merge(CAU2, AFHI2, by =c('gene', 'Phenotype'))
#total <- merge(total, CAU2, by = c('gene', 'Phenotype'))
print(head(total))
write.csv(total, 'zscores_sig_total.csv', quote = F, row.names=F)

pearson <- cor.test(total$zscore.x,total$zscore.y, method="pearson")
print(pearson)
#This is figure S1 in the paper
p <- ggplot(total, aes(x=zscore.x, y=zscore.y, color=Phenotype)) + geom_point() + theme_minimal()+
  scale_color_manual(values = wes_palette("GrandBudapest1", 11, type = "continuous"))+ geom_smooth(method='lm', fullrange=TRUE, aes(group=1)) + annotate(geom="text", x=-10, y =15, label=paste("r = ", round(pearson$estimate, digits = 3), sep = ''))+ labs( x="EUR Significant Gene Z-scores", y = "AFHI Significant Gene Z-scores") +theme(legend.position="bottom")+ ylim(-20,20) + xlim(-20, 20) +ggsave('AFHI_EUR_sig_gene_color_zscore.png', width = 8, height = 6)
p <- ggplot(total, aes(x=zscore.x, y=zscore.y, label = gene_name.x)) + geom_point() + geom_label() + geom_smooth(method = 'lm') + labs( x="EUR Significant Gene Z-scores", y = "AFHI Significant Gene Z-scores") +ylim(-20,20) + xlim(-20, 20) +theme_minimal()+ggsave('AFHI_EUR_sig_gene_label_zscore.png')
p <- ggplot(total, aes(x=zscore.x, y=zscore.y)) + geom_point() + geom_smooth(method = 'lm') + labs( x="EUR Significant Gene Z-scores", y = "AFHI Significant Gene Z-scores") +ylim(-20,20) + xlim(-20, 20) +ggsave('AFHI_EUR_sig_gene_zscore.tiff', device = 'tiff')

#Trying to do new thing when I only take all EUR-sig genes in AFHI whether they are sig or not in AFHI
CAU2$pairs <- paste(CAU2$gene_name, CAU2$Phenotype)
AFHI2$pairs <- paste(AFHI2$gene_name, AFHI2$Phenotype)
CAU_genes <- as.list(CAU2$pairs)
AFHI_genes <- as.list(AFHI2$pairs)
gene_list <- c(CAU_genes, AFHI_genes)
length(gene_list)
print(nrow(AFHI))
print(unique(AFHI$phenotype))
print(unique(AFHI2$Phenotype))
AFHI$pairs <- paste(AFHI$gene_name, AFHI$phenotype)
print(head(AFHI))
AFHI3 <- AFHI[AFHI$pairs %in% gene_list, ] #subset(AFHI, pairs %in% AFHI_genes)
print(AFHI3$pairs)
CAU$pairs <- paste(CAU$gene_name, CAU$phenotype)
CAU3 <- subset(CAU, pairs %in% gene_list)
print(CAU3$pairs)

#This is the figure 1 in the paper
total1 <- merge(CAU3, AFHI3, by =c('gene', 'phenotype'))
pearson1 <- cor.test(total1$zscore.x,total1$zscore.y, method="pearson")
print(pearson1)
p <- ggplot(total1, aes(x=zscore.x, y=zscore.y)) + geom_point()+ labs( x="EUR Gene Z-scores", y = "AFHI Gene Z-scores")  + annotate(geom="text",size=6, x=-10, y =15, label=paste("R = ", round(pearson1$estimate, digits = 3), sep = ''))+geom_smooth(method = 'lm', fullrange=TRUE, aes(group=1)) +geom_density_2d() + theme_minimal()+ylim(-20,20) + theme_bw(20)+ xlim(-20, 20) + ggsave('AFHI_EUR_sig_genes_2_zscore.tiff', device = 'tiff')
