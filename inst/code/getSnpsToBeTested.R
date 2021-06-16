library(VariantAnnotation)

# From stack overflow: 
bind.pad <- function(l,
                     side="r",
                     len=max(vapply(l,length, FUN.VALUE=1))
){
    if (side %in% c("b", "r")) {
        out <- sapply(l, 'length<-', value=len)
    } else {
        out <- sapply(sapply(sapply(l, rev), 'length<-', value=len, simplify=F), rev)}
    if (side %in% c("r", "l")){
        out <- t(out)
    }
    out
}

inputfile <- "00-common_all_snps_only_maf1perc.vcf.gz"
svp <- ScanVcfParam(info=c("CAF"))
vcf <- readVcf(inputfile, param=svp)

set.seed(10)
snps <- sample(names(vcf),1000)
vcf_subset <- vcf[snps]
gr <- rowRanges(vcf_subset)
snps_df <- data.frame(chr=seqnames(gr),
                      pos=start(gr),
                      nuc_ref=as.character(gr$REF))



# Getting alternative alleles:
alt <- rowRanges(vcf_subset)$ALT
alt <- CharacterList(alt)
alt <- bind.pad(alt)

# Constructing CAF matrix:
caf <- info(vcf_subset)$CAF
caf <- bind.pad(caf)
caf <- as.data.frame(caf, stringsAsFactors=FALSE)
colnames(caf) <- c("FREQREF", "FREQALT1","FREQALT2","FREQALT3")
for (i in 1:4){
    caf[,i] <- as.numeric(caf[,i])
}
caf_alt <- caf[,2:4] #alt allele only


# Only keeping one alt allele (maximum CAF per row):
x <- as.matrix(caf_alt)
x[is.na(x)] <- 0 #rowRanks doesn't accept NAs
os <- matrixStats::rowRanks(-x)
indices <- which(t(os)==1)

# Selected allele:
alt <- t(alt)[indices]
caf_alt_selected <- t(as.matrix(caf_alt))[indices]

#Creating bi-allelic CAF matrix:
caf <- cbind(caf[,1], caf_alt_selected)
colnames(caf) <- c("FREQREF", "FREQALT")
caf[,1] <- as.numeric(caf[,1])
caf[,2] <- as.numeric(caf[,2])

# Adding minor allele info:
snps_df$nuc_alt <- alt
snps_df$maf_ref <- caf[,1]
snps_df$maf_alt <- caf[,2] 
snpsToBeTested <- snps_df
save(snpsToBeTested, file="snpsToBeTested.rda")

