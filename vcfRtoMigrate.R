setwd("~/Dropbox/P.glabra/Rad-seq/migrate-n_glabra/")

library(vcfR)

vcf <- read.vcfR(file = "glabra_migrate_reduced.vcf")

# Scanning file to determine attributes.
# File attributes:
#   meta lines: 11
# header_line: 12
# variant count: 40591
# column count: 39
# Meta line 11 read in.
# All meta lines processed.
# gt matrix initialized.
# Character matrix gt created.
# Character matrix gt rows: 40591
# Character matrix gt cols: 39
# skip: 0
# nrows: 40591
# row_num: 0
# Processed variant: 40591
# All variants processed

my_pop <- as.factor(c("AUS","AUS","AUS","AUS","AUS","AUS","AUS","AUS","AUS","AUS",
                      "CHI","CHI","CHI","CHI","CHI","CHI","CHI","CHI","CHI","CHI",
                      "NZ","NZ","NZ","NZ","NZ","NZ","NZ","NZ","NZ","NZ"))

vcfR2migrate(vcf = vcf, pop = my_pop, in_pop = c("AUS", "CHI", "NZ"), out_file = "glabra_migrateN_infile_H", method = 'H')

vcfR2migrate(vcf = vcf, pop = my_pop, in_pop = c("AUS", "CHI", "NZ"), out_file = "glabra_migrateN_infile_N", method = 'N')
