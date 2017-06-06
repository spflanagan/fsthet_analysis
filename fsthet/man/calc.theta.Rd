\name{calc.theta}
\alias{calc.theta}
\title{This calculates theta.}
\description{
	This calculates Weir (1990)'s theta.
	See the vignette for details on the calculation of beta.
}
\usage{
calc.theta(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{T2. This is a single numerical value.}
	\item{fst}{The calculated theta value (T1/T2) for this locus.}
}
\examples{
\dontrun{
  gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
  gpop<-my.read.gpop(gfile)
  theta1<-calc.theta(gpop,3) #calculate theta for the first SNP
}
}