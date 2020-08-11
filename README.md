## Transcriptome Models: 
#### Multi-Ethic Study of Atherosclerosis (MESA)
1. ALL model: European, African American, and Hispanic populations combined
2. AFA model: African American population, n = 233
3. AFHI model: African American population and Hispanic population combined
4. HIS model: Hispanic population, n = 352
5. EUR (CAU) model: European population, n = 578


## [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/31217584)
### GWAS Summary Statistics Tested:
- All Wojcik et al 2019 PAGE GWAS summary statistics (28 clinical and behavioral phenotypes)
- Body mass index, chronic kidney disease, c-reactive protein levels, type-II diabetes, diastolic blood pressure, end-renal kidney failure, fasting insulin, fasting glucose, glomerular filtration estimation rate, hemoglobin levels, hypertension, HDL cholesterol, LDL cholesterol, height, QT interval, QRS duration, PR interval, platelet count, coffee consumption, smoking behavior, total cholesterol, mean corpuscular hemoglobin concentration, triglyceride levels, white blood cell (WBC) count, waist-hip ratio in males, waist-hip ratio in females, and waist-hip ratio. 


## Scripts: 
### 01: Run S-PrediXcan with the GWAS Summary Statistics from the study.
*Note that these files were accessed before the availability of the harmonised form.*

### 02: Generate Manhattan and QQ plots for each of the S-PrediXcan tests. 

### 03: Generate barplots to display S-PrediXcan significant genes results.

### 04: Combine all of the significant genes files together

### 05: Check if the new harmonised data released by GWAS Catalog matches the old GWAS Summary Statistics data and output by running S-PrediXcan on the White Blood Cell Count phenotype.

### 06-09: Preliminary script to properly format the input for LOCUS Comparer. 

### 10: Run Locus Comparer with eQTL and GWAS data to depict LD.

### 11: Merge the Locus Compare Plots for each Population into a pdf.

### 12: Merge PhenomeXcan files and create PhenomeXcan Table for AFHI and EUR results.
