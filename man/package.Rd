\name{BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor}
\docType{package}

\alias{BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor-package}
\alias{BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor}
\alias{Hsapiens}

\title{Full genome sequences for Homo sapiens (UCSC version hg38, based on GRCh38.p12) with injected minor alleles (dbSNP151)}

\description{
  Full genome sequences for Homo sapiens (Human) as provided by UCSC (hg38, based on GRCh38.p12) with minor alleles injected from dbSNP151, and stored in Biostrings objects. Only common
  single nucleotide variants (SNVs) with at least one alternate allele with frequency greater than 0.01 were considered. For SNVs with more than 1 alternate allele, the most frequent allele was chosen as the minor allele to be injected into the reference genome.
}


\author{Jean-Philippe Fortin}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor
genome_min <- BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor
head(seqlengths(genome_min))

# Getting nucleotide at SNP rs12813551 (C/T, MAF>0.5)
# Minor allele genome has a C:
chr <- "chr12"
pos <- 25241845L
getSeq(genome_min, chr, start=pos, end=pos)

# Reference genome has the minor allele, C:
if (require(BSgenome.Hsapiens.UCSC.hg38)){
    genome_ref <- BSgenome.Hsapiens.UCSC.hg38
    getSeq(genome_ref, chr, start=pos, end=pos)
}

# Major allele genome has a T:
if (require(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)){
    genome_maj <- BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major
    getSeq(genome_maj, chr, start=pos, end=pos)
}
}


\keyword{package}
\keyword{data}
