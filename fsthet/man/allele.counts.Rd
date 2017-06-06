\name{allele.counts}
\alias{allele.counts}
\title{This counts the number of alleles at a locus. }
\description{
  This counts the number of times each allele occurs at a locus from a list of genotypes (the sum of all the counts is 2*number of individuals).
}
\usage{
allele.counts(genotypes)
}
\arguments{
  \item{genotypes}{A list of genotypes.}
}
\value{
	\item{AlleleCounts}{The number of times each allele is recorded at the locus.}
}
\examples{
  #create a random sample of genotypes
  genotypes<-sample(c("0101","0102","0202"),50,replace=TRUE)
  counts<-allele.counts(genotypes)
}