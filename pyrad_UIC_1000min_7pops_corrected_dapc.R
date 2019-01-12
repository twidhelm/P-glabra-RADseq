setwd("~/Dropbox/P.glabra/Rad-seq/dapc/pyrad_UIC_min1000/")
getwd()

### Copying data files
#system("cp /home/fgrewe/Usnea/Usnea_RAD1234/6-pyrad7/outfiles/c90d6m4p3.vcf .")
#system("cp /home/fgrewe/Usnea/Usnea_RAD1234/9-DAPC/Usnea.pop .")
system("vcftools --vcf glabra_UIC_1000min.vcf --max-missing 0.5 --maf 0.05 --out glabra_UIC_min1000_corrected_7pop.vcftools.vcf --recode 2>&1", intern = T)

# [1] ""                                                        
# [2] "VCFtools - 0.1.15"                                       
# [3] "(C) Adam Auton and Anthony Marcketta 2009"               
# [4] ""                                                        
# [5] "Parameters as interpreted:"                              
# [6] "\t--vcf glabra_UIC_1000min.vcf"                          
# [7] "\t--maf 0.05"                                            
# [8] "\t--max-missing 0.5"                                     
# [9] "\t--out glabra_data2-noOG.vcftools.vcf"                  
# [10] "\t--recode"                                              
# [11] ""                                                        
# [12] "Eighth Header entry should be INFO: INFO    "            
# [13] "After filtering, kept 216 out of 216 Individuals"        
# [14] "Outputting VCF file..."                                  
# [15] "After filtering, kept 2440 out of a possible 88090 Sites"
# [16] "Run Time = 3.00 seconds"   

### Required packages
library(vcfR)
library(adegenet)
library(hierfstat)
library(qqman)
library(mmod)

### Loading vcf file into R genind object  
vcf <- read.vcfR("glabra_data2-noOG.vcftools.vcf.recode.vcf")

# Scanning file to determine attributes.
# File attributes:
#   meta lines: 11
# header_line: 12
# variant count: 2440
# column count: 225
# Meta line 11 read in.
# All meta lines processed.
# gt matrix initialized.
# Character matrix gt created.
# Character matrix gt rows: 2440
# Character matrix gt cols: 225
# skip: 0
# nrows: 2440
# row_num: 0
# Processed variant: 2440
# All variants processed

data.genlight <- vcfR2genlight(vcf, n.cores = 3)

# Warning message:
#   In vcfR2genlight(vcf, n.cores = 3) :
#   Found 29 loci with more than two alleles.
# Objects of class genlight only support loci with two alleles.
# 29 loci will be omitted from the genlight object.

pop.file <- read.table('pop_file_1000min_7pops.txt', header=F)
pop(data.genlight) <- pop.file[,2]
## To get this to work, I had to remove the "." from the loci names 
data.genind <- df2genind(as.data.frame(data.genlight), pop=pop(data.genlight), ploidy=1, ind.names=indNames(data.genlight), loc.names=locNames(data.genlight))
data.genind

# /// GENIND OBJECT /////////
  
#   // 216 individuals; 2,411 loci; 4,822 alleles; size: 5.1 Mb

# // Basic content
# @tab:  216 x 4822 matrix of allele counts
# @loc.n.all: number of alleles per locus (range: 2-2)
# @loc.fac: locus factor for the 4822 columns of @tab
# @all.names: list of allele names for each locus
# @ploidy: ploidy of each individual  (range: 1-1)
# @type:  codom
# @call: df2genind(X = as.data.frame(data.genlight), ind.names = indNames(data.genlight), 
#                  loc.names = locNames(data.genlight), pop = pop(data.genlight), 
#                  ploidy = 1)

# // Optional content
# @pop: population of each individual (group size range: 1-66)

### Pairwise Fst
pwD <- pairwise_D(data.genind)
pwGstN <- pairwise_Gst_Nei(data.genind)
pwGstH <- pairwise_Gst_Hedrick(data.genind)
pwGstN
pwGstH
pwD

# > pwD <- pairwise_D(data.genind)
# There were 50 or more warnings (use warnings() to see the first 50)
# > pwGstN <- pairwise_Gst_Nei(data.genind)
# There were 50 or more warnings (use warnings() to see the first 50)
# > pwGstH <- pairwise_Gst_Hedrick(data.genind)
# There were 50 or more warnings (use warnings() to see the first 50)
# > pwGstN
# CAM        CHI        NZN        NZS        PHO        TAS
# CHI 0.31474889                                                       
# NZN 0.24295927 0.14727745                                            
# NZS 0.21299886 0.11662000 0.06802304                                 
# PHO 0.66314342 0.41690387 0.29779390 0.33469248                      
# TAS 0.46419314 0.27857919 0.22120872 0.19968760 0.49791119           
# VIC 0.36443624 0.23400576 0.17938065 0.15724973 0.42247216 0.05967977
# > pwGstH
# CAM       CHI       NZN       NZS       PHO       TAS
# CHI 0.5463886                                                  
# NZN 0.4671132 0.3262028                                        
# NZS 0.4112208 0.2618454 0.1660587                              
# PHO 0.8509863 0.6847103 0.5537737 0.5973598                    
# TAS 0.6928019 0.5170144 0.4465989 0.4049405 0.7497377          
# VIC 0.6014184 0.4603033 0.3835010 0.3380814 0.6859595 0.1326597
# > pwD
# CAM        CHI        NZN        NZS        PHO        TAS
# CHI 0.12968395                                                       
# NZN 0.12507140 0.09345387                                            
# NZS 0.09251908 0.06694940 0.04432351                                 
# PHO 0.26428292 0.23385672 0.17529944 0.19225187                      
# TAS 0.16052337 0.14400116 0.13222154 0.10799143 0.25337761           
# VIC 0.14432007 0.13055627 0.11397770 0.09106509 0.22650667 0.02256386

obj_seplocus <- seploc(data.genind);
#Calculation of Gst Nei
obj_pwGst <- lapply(obj_seplocus, pairwise_Gst_Nei);
obj_pairwiseGstnum <- sapply(obj_pwGst, as.numeric);
obj_pairwiseGstnum[obj_pairwiseGstnum<0] <- 0;
#Calculation of Gst Hedrick
obj_pwGstH <- lapply(obj_seplocus, pairwise_Gst_Hedrick);
obj_pairwiseGstHnum <- sapply(obj_pwGstH, as.numeric);
obj_pairwiseGstHnum[obj_pairwiseGstHnum<0] <- 0;
#Calculation of D
obj_pwD <- lapply(obj_seplocus, pairwise_D);
obj_pairwiseDnum <- sapply(obj_pwD, as.numeric);
obj_pairwiseDnum[obj_pairwiseDnum<0] <- 0;

#Priting 3 figures arranged in 3 rows and 1 column
par(mfrow=c(3,1))
hist(obj_pairwiseGstnum, breaks=100, main="Nei's Gst", xlab="Pairwise Values of Nei's Gst", col='lightgreen');
hist(obj_pairwiseGstHnum, breaks=100, main="Hedrick's G'st", xlab="Pairwise Values of Hedrick's G'st", col='lightblue');
hist(obj_pairwiseDnum, breaks=100, main="Jost's D", xlab="Pairwise Values of Jost's D", col='lightyellow');

### DAPC
dapc1 <- dapc(data.genind, n.pca = 60, n.da = 2)
mycol <- c('darkorange', 'darkgreen', 'red1', 'red4', 'black', 'royalblue4', 'royalblue2')

#print DAPC
par(mfrow=c(2,1))
scatter(dapc1, col = mycol, scree.pca = TRUE, posi.pca = "topright", posi.da = "bottomright")
compoplot(dapc1, posi='topleft', ncol = 1, col.pal = mycol, cleg = 0.8, cex.names = 0.4)

### save all in pdf
pdf("glabra_UIC_pyrad_1000min_7pop_corrected.pdf")
par(mfrow=c(3,1))
hist(obj_pairwiseGstnum, breaks=100, main="Nei's Gst", xlab="Pairwise Values of Nei's Gst", col='lightgreen');
hist(obj_pairwiseGstHnum, breaks=100, main="Hedrick's G'st", xlab="Pairwise Values of Hedrick's G'st", col='lightblue');
hist(obj_pairwiseDnum, breaks=100, main="Jost's D", xlab="Pairwise Values of Jost's D", col='lightyellow');
par(mfrow=c(2,1))
scatter(dapc1, col = mycol, scree.pca = TRUE, posi.pca = "topright", posi.da = "bottomright")
compoplot(dapc1, posi='bottomleft', ncol = 1, col = mycol, cleg = 0.8, cex.names = 0.4)
dev.off()

