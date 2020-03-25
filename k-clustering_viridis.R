### Clustering plots
library(vcfR)
glabra.VCF <- read.vcfR("glabra_no_CAM_PFR_PHO_min1000.recode.vcf")

## Scanning file to determine attributes.
## File attributes:
##   meta lines: 11
## header_line: 12
## variant count: 3756
## column count: 310
## Meta line 11 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
## Character matrix gt rows: 3756
## Character matrix gt cols: 310
## skip: 0
## nrows: 3756
## row_num: 0
##  variant: 3756
## All variants processed

# Read in pop file
pop.data <- read.table("pop_file_301samples.txt", sep = "\t", header = TRUE)

# We can now check that all the samples in the VCF and the population data frame are included:
all(colnames(glabra.VCF@gt)[-1] == pop.data$AccessID)

## [1] TRUE

# Convert to a genlight object
gl.glabra <- vcfR2genlight(glabra.VCF)

## Loading required namespace: adegenet
## Warning message:
##   In vcfR2genlight(glabra.VCF) : Found 39 loci with more than two alleles.
## Objects of class genlight only support loci with two alleles.
## 39 loci will be omitted from the genlight object.

# Specify ploidy
#ploidy(gl.glabra) <- 1

#pop(gl.glabra) <- pop.data$Country

#gl.glabra

## /// GENLIGHT OBJECT /////////
##   
##   // 301 genotypes,  3,717 binary SNPs, size: 2.7 Mb
## 418932 (37.44 %) missing data
## 
## // Basic content
## @gen: list of 301 SNPbin
## @ploidy: ploidy of each individual  (range: 1-1)
## 
## // Optional content
## @ind.names:  301 individual labels
## @loc.names:  3717 locus labels
## @chromosome: factor storing chromosomes of the SNPs
## @position: integer storing positions of the SNPs
## @pop: population of each individual (group size range: 21-77)
## @other: a list containing: elements without names 

# K-means clustering
library(adegenet)
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.glabra, n.pca = 400, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

# Visualizing K- means clustering
library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1

# DAPC
my_k <- 3:6

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.glabra, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.glabra, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl.glabra, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

#scatter plot
my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

library(viridis)
my_pal <- cividis(6)
  
#RColorBrewer::brewer.pal(n=8, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c("#00204DFF", "#31446BFF", "#666970FF", "#958F78FF", "#CBBA69FF", "#FFEA46FF"))
p2 <- p2 + scale_fill_manual(values=c("#00204DFF", "#31446BFF", "#666970FF", "#958F78FF", "#CBBA69FF", "#FFEA46FF"))
p2

#barplot
tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data$Country
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data$Country
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Dark2")
p3 <- p3 + scale_fill_manual(values=c("#00204DFF", "#31446BFF", "#666970FF", "#958F78FF", "#CBBA69FF", "#FFEA46FF"))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3

# Bring all figures together for publication
library(ggpubr)
#tiff('dapc__k3_5_dapc.tiff', width=6.5, height=6.5, units='in', compression='lzw', res=300)
ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)
