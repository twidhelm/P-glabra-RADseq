setwd("~/Dropbox/P.glabra/Rad-seq/dapc/pyrad_UIC_migrate_reduced//")
getwd()

### Copying data files
#system("cp /home/fgrewe/Usnea/Usnea_RAD1234/6-pyrad7/outfiles/c90d6m4p3.vcf .")
#system("cp /home/fgrewe/Usnea/Usnea_RAD1234/9-DAPC/Usnea.pop .")
system("vcftools --vcf glabra_migrate_reduced.vcf --max-missing 0.5 --maf 0.05 --out glabra_migrate_reduced.vcftools.vcf --recode 2>&1", intern = T)

# [1] ""                                                         
# [2] "VCFtools - 0.1.15"                                        
# [3] "(C) Adam Auton and Anthony Marcketta 2009"                
# [4] ""                                                         
# [5] "Parameters as interpreted:"                               
# [6] "\t--vcf glabra_migrate_reduced.vcf"                       
# [7] "\t--maf 0.05"                                             
# [8] "\t--max-missing 0.5"                                      
# [9] "\t--out glabra_migrate_reduced.vcftools.vcf"              
# [10] "\t--recode"                                               
# [11] ""                                                         
# [12] "Eighth Header entry should be INFO: INFO    "             
# [13] "After filtering, kept 30 out of 30 Individuals"           
# [14] "Outputting VCF file..."                                   
# [15] "rval"
# [16] "Run Time = 1.00 seconds"

### Required packages
library(vcfR)
library(adegenet)
library(hierfstat)
library(qqman)
library(mmod)

### Loading vcf file into R genind object  
vcf <- read.vcfR("glabra_migrate_reduced.vcftools.vcf.recode.vcf")

# Scanning file to determine attributes.
# File attributes:
#   meta lines: 11
# header_line: 12
# variant count: 14229
# column count: 39
# Meta line 11 read in.
# All meta lines processed.
# gt matrix initialized.
# Character matrix gt created.
# Character matrix gt rows: 14229
# Character matrix gt cols: 39
# skip: 0
# nrows: 14229
# row_num: 0
# Processed variant: 14229
# All variants processed

data.genlight <- vcfR2genlight(vcf, n.cores = 3)

# Warning message:
# In vcfR2genlight(vcf, n.cores = 3) :
#   Found 200 loci with more than two alleles.
# Objects of class genlight only support loci with two alleles.
# 200 loci will be omitted from the genlight object.

pop.file <- read.table('pop_migrate_30indiv.txt', header=F)
pop(data.genlight) <- pop.file[,2]
## To get this to work, I had to remove the "." from the loci names 
data.genind <- df2genind(as.data.frame(data.genlight), pop=pop(data.genlight), ploidy=1, ind.names=indNames(data.genlight), loc.names=locNames(data.genlight))
data.genind

# /// GENIND OBJECT /////////

# // 30 individuals; 14,029 loci; 28,058 alleles; size: 9.7 Mb

# // Basic content
# @tab:  30 x 28058 matrix of allele counts
# @loc.n.all: number of alleles per locus (range: 2-2)
# @loc.fac: locus factor for the 28058 columns of @tab
# @all.names: list of allele names for each locus
# @ploidy: ploidy of each individual  (range: 1-1)
# @type:  codom
# @call: df2genind(X = as.data.frame(data.genlight), ind.names = indNames(data.genlight), 
#                  loc.names = locNames(data.genlight), pop = pop(data.genlight), 
#                  ploidy = 1)

# // Optional content
# @pop: population of each individual (group size range: 10-10)

### Pairwise Fst
pwD <- pairwise_D(data.genind)
# There were 50 or more warnings (use warnings() to see the first 50)
pwGstN <- pairwise_Gst_Nei(data.genind)
# There were 46 warnings (use warnings() to see them)
pwGstH <- pairwise_Gst_Hedrick(data.genind)
# There were 46 warnings (use warnings() to see them)
pwGstN
pwGstH
pwD

# > pwGstN
# AUS       CHI
# CHI 0.2809691          
# NZ  0.2119006 0.1620554
# > pwGstH
# AUS       CHI
# CHI 0.5261499          
# NZ  0.4397738 0.3568245
# > pwD
# AUS       CHI
# CHI 0.1558257          
# NZ  0.1385118 0.1080490

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
dapc1 <- dapc(data.genind, n.pca = 5, n.da = 1)
dapc2 <- dapc(data.genind, n.pca = 5, n.da = 2)
mycol <- c('royalblue2', 'darkgreen', 'red1')

#print DAPC
par(mfrow=c(2,1))
scatter(dapc1, col = mycol, scree.pca = TRUE, posi.pca = "topright", posi.da = "bottomright")
compoplot(dapc1, posi='topleft', ncol = 1, col.pal = mycol, cleg = 0.8, cex.names = 0.4)

par(mfrow=c(2,1))
scatter(dapc2, col = mycol, scree.pca = TRUE, posi.pca = "topright", posi.da = "bottomright")
compoplot(dapc2, posi='topleft', ncol = 1, col.pal = mycol, cleg = 0.8, cex.names = 0.4)


### save all in pdf
pdf("glabra_UIC_1000min_5pops_corrected.pdf")
par(mfrow=c(3,1))
hist(obj_pairwiseGstnum, breaks=100, main="Nei's Gst", xlab="Pairwise Values of Nei's Gst", col='lightgreen');
hist(obj_pairwiseGstHnum, breaks=100, main="Hedrick's G'st", xlab="Pairwise Values of Hedrick's G'st", col='lightblue');
hist(obj_pairwiseDnum, breaks=100, main="Jost's D", xlab="Pairwise Values of Jost's D", col='lightyellow');
par(mfrow=c(2,1))
scatter(dapc1, col = mycol, scree.pca = TRUE, posi.pca = "topright", posi.da = "bottomright")
compoplot(dapc1, posi='topleft', ncol = 1, col.pal = mycol, cleg = 0.8, cex.names = 0.4)
dev.off()

