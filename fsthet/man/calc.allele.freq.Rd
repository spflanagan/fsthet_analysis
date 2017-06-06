\name{calc.allele.freq}
\alias{calc.allele.freq}
\title{This calculates allele frequencies. }
\description{
  This calculates allele frequencies from a list of genotypes.
}
\usage{
calc.allele.freq(genotypes)
}
\arguments{
  \item{genotypes}{A list of genotypes.}
}
\value{
	\item{obs.af}{A list of observed allele frequencies in the genotypes list.}
}
\examples{
  #create a random sample of genotypes
  genotypes<-sample(c("0101","0102","0202"),50,replace=TRUE)
  af<-calc.allele.freq(genotypes)
}