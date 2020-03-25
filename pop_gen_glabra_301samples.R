library(vcfR)
library(poppr)
library(ape)
library(wesanderson)
library(adegenet)

glabra.VCF <- read.vcfR("glabra_no_CAM_PFR_PHO_min1000.vcf")

## Scanning file to determine attributes.
## File attributes:
##   meta lines: 11
## header_line: 12
## variant count: 119298
## column count: 310
## Meta line 11 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
## Character matrix gt rows: 119298
## Character matrix gt cols: 310
## skip: 0
## nrows: 119298
## row_num: 0
## Processed variant: 119298
## All variants processed

glabra.VCF

## ***** Object of Class vcfR *****
## 301 samples
## 28253 CHROMs
## 119,298 variants
## Object size: 283.7 Mb
## 74.13 percent missing data
## *****        *****         *****

# Modify VCF file
system("vcftools --vcf glabra_no_CAM_PFR_PHO_min1000.vcf --max-missing 0.5 --maf 0.05 --out glabra_no_CAM_PFR_PHO_min1000 --recode 2>&1", intern = T)

# Read in modified VCF file
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
## Processed variant: 3756
## All variants processed

glabra.VCF

## ***** Object of Class vcfR *****
## 301 samples
## 1977 CHROMs
## 3,756 variants
## Object size: 9 Mb
## 37.44 percent missing data
## *****        *****         *****

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
ploidy(gl.glabra) <- 1

pop(gl.glabra) <- pop.data$Country

gl.glabra

##  /// GENLIGHT OBJECT /////////
## 
## // 301 genotypes,  3,717 binary SNPs, size: 2.7 Mb
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



### Population genetic analyses for GBS data
# Distance Matrix

glabra.dist <- poppr::bitwise.dist(gl.glabra)

tree <- aboot(gl.glabra, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

tree <- root(tree, node = 374)
tree <- ladderize(tree, right = FALSE)

cols <- wes_palette("Cavalcanti1", n = nPop(gl.glabra))
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.glabra)], type = "phylogram")
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("CHI", "NZN", "NZS", "TAS", "VIC"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

### Minimum spanning networks

library(igraph)

glabra.dist <- bitwise.dist(gl.glabra)
glabra.msn <- poppr.msn(gl.glabra, glabra.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(gl.glabra))
names(node.size) <- indNames(gl.glabra)
vertex.attributes(glabra.msn$graph)$size <- node.size

set.seed(9)
plot_poppr_msn(gl.glabra, glabra.msn , palette = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15"), 
               gadj = 70, inds = "cool", nodelab = 1000)

### Principal components analysis

glabra.pca <- glPca(gl.glabra, nf = 4)
barplot(100*glabra.pca$eig/sum(glabra.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

glabra.pca.scores <- as.data.frame(glabra.pca$scores)
glabra.pca.scores$pop <- pop(gl.glabra)

library(ggplot2)
set.seed(9)
p <- ggplot(glabra.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=3)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()

p

### DAPC

glabra.dapc <- dapc(gl.glabra, n.pca = 4, n.da = 4)

scatter(glabra.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "topright", scree.pca = TRUE,
        posi.pca = "bottomleft", cleg = 0.75)

compoplot(glabra.dapc, col = cols, posi = 'top')

dapc.results <- as.data.frame(glabra.dapc$posterior)
dapc.results$pop <- pop(gl.glabra)
dapc.results$indNames <- rownames(dapc.results)

library(reshape2)
dapc.results <- melt(dapc.results)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

