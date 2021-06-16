#samtools, BCFtools needed
#module load BCFtools
#module load samtools
#module load HTSlib

#Downloading hg39 fasta file:
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

#Downloading common variants:
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz.tbi

#Processing major/minor SNPs:
Rscript processSNPs.R

#To output minor allele genome:
vcf="hg38.dbsnp151.minor.allele.vcf.gz"
outputgenome="hg38_minor.fa"
refgenome="hg38.fa"
bcftools consensus -f $refgenome -s minor $vcf > $outputgenome

#To split fastas for BSgenome:
pip install pyfaidx
faidx -x hg38_minor.fa