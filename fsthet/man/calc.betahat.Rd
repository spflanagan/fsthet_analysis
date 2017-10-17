\name{calc.betahat}
\alias{calc.betahat}
\title{This calculates beta-hat, the Fst value used in Lositan. }
\description{
	This calculates Weir & Cockerham (1993)'s beta-hat.
	Beaumont & Nichols (1996) used this formulation in FDIST2 (and is implemented in Lositan)
	See the vignette for details on the calculation of beta.
}
\usage{
calc.betahat(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{HB (or 1-F1). This is a single numerical value.}
	\item{fst}{The calculated betahat value ((F0-F1)/(1-F1))for this locus.}
}
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20),
    loc0=sample(c("0101","0102","0202"),40,replace=TRUE))
  bh<-calc.betahat(gpop, 3) #calculate betahat for the SNP
  \donttest{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.gpop(gfile)
    beta1<-calc.betahat(gpop,3) #calculate betahat for the first SNP
  }
}