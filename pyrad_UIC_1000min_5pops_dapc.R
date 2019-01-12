setwd("~/Dropbox/P.glabra/Rad-seq/dapc/pyrad_UIC_min1000/")
getwd()

### Copying data files
#system("cp /home/fgrewe/Usnea/Usnea_RAD1234/6-pyrad7/outfiles/c90d6m4p3.vcf .")
#system("cp /home/fgrewe/Usnea/Usnea_RAD1234/9-DAPC/Usnea.pop .")
system("vcftools --vcf glabra_UIC_1000min_5pop.vcf --max-missing 0.5 --maf 0.05 --out glabra_UIC_1000min_5pops.vcftools.vcf --recode 2>&1", intern = T)

# [1] ""                                                         "VCFtools - 0.1.15"                                       
# [3] "(C) Adam Auton and Anthony Marcketta 2009"                ""                                                        
# [5] "Parameters as interpreted:"                               "\t--vcf glabra_UIC_1000min_5pop.vcf"                     
# [7] "\t--maf 0.05"                                             "\t--max-missing 0.5"                                     
# [9] "\t--out glabra_UIC_1000min_5pops.vcftools.vcf"            "\t--recode"                                              
# [11] ""                                                         "Eighth Header entry should be INFO: INFO    "            
# [13] "After filtering, kept 202 out of 202 Individuals"         "Outputting VCF file..."                                  
# [15] "After filtering, kept 2255 out of a possible 81083 Sites" "Run Time = 2.00 seconds"  

### Required packages
library(vcfR)
library(adegenet)
library(hierfstat)
library(qqman)
library(mmod)

### Loading vcf file into R genind object  
vcf <- read.vcfR("glabra_UIC_1000min_5pops.vcftools.vcf.recode.vcf")

# Scanning file to determine attributes.
# File attributes:
#   meta lines: 11
# header_line: 12
# variant count: 2255
# column count: 211
# Meta line 11 read in.
# All meta lines processed.
# gt matrix initialized.
# Character matrix gt created.
# Character matrix gt rows: 2255
# Character matrix gt cols: 211
# skip: 0
# nrows: 2255
# row_num: 0
# Processed variant: 2255
# All variants processed

data.genlight <- vcfR2genlight(vcf, n.cores = 3)

# Warning message:
#   In vcfR2genlight(vcf, n.cores = 3) :
#   Found 31 loci with more than two alleles.
# Objects of class genlight only support loci with two alleles.
# 31 loci will be omitted from the genlight object.

pop.file <- read.table('pop_file_1000min_5pops.txt', header=F)
pop(data.genlight) <- pop.file[,2]
## To get this to work, I had to remove the "." from the loci names 
data.genind <- df2genind(as.data.frame(data.genlight), pop=pop(data.genlight), ploidy=1, ind.names=indNames(data.genlight), loc.names=locNames(data.genlight))
data.genind

# /// GENIND OBJECT /////////
  
#   // 202 individuals; 2,224 loci; 4,448 alleles; size: 4.5 Mb

# // Basic content
# @tab:  202 x 4448 matrix of allele counts
# @loc.n.all: number of alleles per locus (range: 2-2)
# @loc.fac: locus factor for the 4448 columns of @tab
# @all.names: list of allele names for each locus
# @ploidy: ploidy of each individual  (range: 1-1)
# @type:  codom
# @call: df2genind(X = as.data.frame(data.genlight), ind.names = indNames(data.genlight), 
#                  loc.names = locNames(data.genlight), pop = pop(data.genlight), 
#                  ploidy = 1)

# // Optional content
# @pop: population of each individual (group size range: 11-60)

### Pairwise Fst
pwD <- pairwise_D(data.genind)
pwGstN <- pairwise_Gst_Nei(data.genind)
pwGstH <- pairwise_Gst_Hedrick(data.genind)
pwGstN
pwGstH
pwD

# > pwGstN
# CHI        NZN        NZS        TAS
# NZN 0.14616298                                 
# NZS 0.11546406 0.06814561                      
# TAS 0.27533185 0.21892748 0.19753898           
# VIC 0.23198384 0.17710142 0.15449780 0.06126183
# > pwGstH
# CHI       NZN       NZS       TAS
# NZN 0.3309072                              
# NZS 0.2644481 0.1701334                    
# TAS 0.5198048 0.4508762 0.4077721          
# VIC 0.4638908 0.3863554 0.3384012 0.1376280
# > pwD
# CHI        NZN        NZS        TAS
# NZN 0.10183161                                 
# NZS 0.07241562 0.04875867                      
# TAS 0.15491221 0.14304747 0.11619887           
# VIC 0.14002090 0.12222243 0.09661457 0.02507160

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
mycol <- c('darkgreen', 'red1', 'red4', 'royalblue4', 'royalblue2')

#print DAPC
par(mfrow=c(2,1))
scatter(dapc1, col = mycol, scree.pca = TRUE, posi.pca = "topright", posi.da = "bottomright")
compoplot(dapc1, posi='topleft', ncol = 1, col.pal = mycol, cleg = 0.8, cex.names = 0.4)

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

