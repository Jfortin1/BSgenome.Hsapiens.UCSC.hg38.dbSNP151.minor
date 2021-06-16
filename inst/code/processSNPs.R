library(VariantAnnotation)
library(matrixStats)
#assuming HTSlib in installed on system
snvFilter <- function(x){
	isSNV(x, singleAltOnly=FALSE)
}

mafFilter <- function(x){
	maf_cutoff=0.01
	#Extracting allele frequencies:
	caf <- info(x)$CAF
	nalleles <- sapply(caf, length)
	maf <- sapply(caf, function(y){
		y <- as.numeric(y[2:length(y)])
		max(y, na.rm=TRUE)
	})
	wh  <- maf>=maf_cutoff
	return(wh)
}

# From stack overflow: 
bind.pad <- function(l, side="r", len=max(sapply(l,length))){
  if (side %in% c("b", "r")) {
    out <- sapply(l, 'length<-', value=len)
  } else {
    out <- sapply(sapply(sapply(l, rev), 'length<-', value=len, simplify=F), rev)}
  if (side %in% c("r", "l")) out <- t(out)
  out
}


########### First step: only keeping SNPs ##############
inputdir  <- "./"
inputfile  <- file.path(inputdir,"00-common_all.vcf.gz")
outputfile <- file.path(inputdir,"00-common_all_snps_only.vcf")
filters <- FilterRules(list(isSNV=isSNV))
filters <- FilterRules(list(isSNV=snvFilter))
tabix.file <- TabixFile(inputfile, yieldSize=1000000)
filterVcf(tabix.file,
	destination=outputfile,
	filters=filters,
	verbose=TRUE
)
system("module load HTSlib;bgzip 00-common_all_snps_only.vcf")
system("module load HTSlib;tabix -p vcf 00-common_all_snps_only.vcf.gz")
#########################################################


######### Second step: only keeping SNPs with at ########
######### least one ALT allele wth AF>0.01   ############
inputfile  <- "00-common_all_snps_only.vcf.gz"
outputfile <- "00-common_all_snps_only_maf1perc.vcf"
filters <- FilterRules(list(maf=mafFilter))
tabix.file <- TabixFile(inputfile, yieldSize=1000000)
filterVcf(tabix.file,
	destination=outputfile,
	filters=filters,
	verbose=TRUE
)
system("module load HTSlib;bgzip 00-common_all_snps_only_maf1perc.vcf")
system("module load HTSlib;tabix -p vcf 00-common_all_snps_only_maf1perc.vcf.gz")
##########################################################


#####################################
# Cleaning to keep two top alleles
#####################################

# Reading in data:
inputfile <- "00-common_all_snps_only_maf1perc.vcf.gz"
svp <- ScanVcfParam(info=c("CAF"))
vcf <- readVcf(inputfile, param=svp)

# Keeping at most 4 alleles:
len <- sapply(info(vcf)$CAF, length)
vcf <- vcf[len<=4]
caf <- info(vcf)$CAF

# Constructing CAF matrix:
caf <- bind.pad(caf)
caf <- as.data.frame(caf, stringsAsFactors=FALSE)
colnames(caf) <- c("FREQREF", "FREQALT1","FREQALT2","FREQALT3")
for (i in 1:4) caf[,i] <- as.numeric(caf[,i])
caf_alt <- caf[,2:4] #alt allele only

# Constructing nucleotide matrix (alt allele)
alt <- rowRanges(vcf)$ALT
alt <- CharacterList(alt)
alt <- bind.pad(alt)

# Only keeping one alt allele (maximum CAF per row):
x <- as.matrix(caf_alt)
x[is.na(x)] <- 0 #rowRanks doesn't accept NAs
os <- matrixStats::rowRanks(-x)
indices <- which(t(os)==1)

# Selected allele:
alt_selected     <- t(alt)[indices]
caf_alt_selected <- t(as.matrix(caf_alt))[indices]

#Creating bi-allelic CAF matrix:
caf <- cbind(caf[,1], caf_alt_selected)
colnames(caf) <- c("FREQREF", "FREQALT")
caf[,1] <- as.numeric(caf[,1])
caf[,2] <- as.numeric(caf[,2])


#Replacing CAF in vcf file by bi-allelic CAF:
CAF <- caf
CAF[,1] <- as.character(CAF[,1])
CAF[,2] <- as.character(CAF[,2])
colnames(CAF) <- NULL
dfs <- CharacterList(asplit(CAF,1))
info(vcf)$CAF <- dfs


#Replacing GR in vcf file by bi-allelic SNP:
gr  <- rowRanges(vcf)
alt <- CharacterList(strsplit(alt_selected,""))
alt <- DNAStringSetList(alt)
VariantAnnotation::fixed(vcf)$ALT <- alt


#####################################
# Creating major/minor alleles VCFs
#####################################

# Genotype encoding for major allele genome:
major <- rep("0|0", nrow(vcf))
major[as.data.frame(caf)$FREQREF < as.data.frame(caf)$FREQALT] <- "1|1"
major <- as.matrix(major)
colnames(major) <- "major"
rownames(major) <- names(vcf)
#     0|0      1|1 
#10736543  1837807 

# Genotype encoding for minor allele genome:
minor <- rep("1|1", nrow(vcf))
minor[as.data.frame(caf)$FREQREF < as.data.frame(caf)$FREQALT] <- "0|0"
minor <- as.matrix(minor)
colnames(minor) <- "minor"
rownames(minor) <- names(vcf)


pheno.major <- DataFrame(name="major")
pheno.minor <- DataFrame(name="minor")
rownames(pheno.major) <- "major"
rownames(pheno.minor) <- "minor"

# From VariantAnnotation package, there is an example data in hg38.
# We can use that to get a proper header and geno info.
hd <- header(vcf)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
hd2 <- header(readVcf(fl))
geno(hd) <- geno(hd2)[1,,drop=FALSE]
info(hd) <- info(hd)["CAF",,drop=FALSE]


#Changing chromosom names:
gr <- rowRanges(vcf)
seqlevels(gr) <- paste0("chr", seqlevels(gr))

vcf_major <- VCF(rowRanges=gr,
 		         colData = pheno.major,
 		         exptData = list(header=hd),
                 fixed = VariantAnnotation::fixed(vcf),
                 info = info(vcf),
                 geno = SimpleList(GT=major))
vcf_minor <- VCF(rowRanges=gr,
 		         colData = pheno.minor,
 		         exptData = list(header=hd),
                 fixed = VariantAnnotation::fixed(vcf),
                 info = info(vcf),
                 geno = SimpleList(GT=minor))

writeVcf(vcf_major, file="hg38.dbsnp151.major.allele.vcf")
writeVcf(vcf_minor, file="hg38.dbsnp151.minor.allele.vcf")

system("module load HTSlib;bgzip hg38.dbsnp151.major.allele.vcf")
system("module load HTSlib;bgzip hg38.dbsnp151.minor.allele.vcf")
system("module load HTSlib;tabix -p vcf hg38.dbsnp151.major.allele.vcf.gz")
system("module load HTSlib;tabix -p vcf hg38.dbsnp151.minor.allele.vcf.gz")

