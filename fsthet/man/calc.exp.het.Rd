\name{calc.exp.het}
\alias{calc.exp.het}
\title{This calculates expected heterozygosities. }
\description{
  This calculates expected heterozygosities from a list of allele frequencies.
}
\usage{
calc.exp.het(af)
}
\arguments{
  \item{af}{is a list of allele frequencies.}
}
\value{
	\item{ht}{The expected heterozygosity under Hardy-Weinberg expectations. This is a single numerical value.}
}
\examples{
  #create a random sample of genotypes
  genotypes<-sample(c("0101","0102","0202"),50,replace=TRUE)
  af<-calc.allele.freq(genotypes)
  hs<-calc.exp.het(af)
}