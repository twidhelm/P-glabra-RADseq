library(phrynomics)

setwd("~/Dropbox/P.glabra/BFD*_SNAPP/data_formatting/")

snpdata50 <- ReadSNP(file = "glabra_UIC_1000min_SP-delim_OG_50_unlinked_snp.phy", extralinestoskip = 1, fileFormat = "phy")
snpdata90 <- ReadSNP(file = "glabra_UIC_1000min_SP-delim_OG_90_unlinked_snp.phy", extralinestoskip = 1, fileFormat = "phy")

# Prepare data 50% dataset for SNPAPP
snpdata50a <- RemoveNonBinary(snpdata50)
snpdata50a <- TranslateBases(snpdata50a, translateMissingChar="?", ordered=TRUE)
WriteSNP(snpdata50a, file="snpdata50.nex", format="nexus", missing="?")
# Removed 18 nonbinary sites

# Prepare data 90% dataset for SNPAPP
snpdata90a <- RemoveNonBinary(snpdata90)
snpdata90a <- TranslateBases(snpdata90a, translateMissingChar="?", ordered=TRUE)
WriteSNP(snpdata90a, file="snpdata90.nex", format="nexus", missing="?")
# Removed 3 nonbinary sites