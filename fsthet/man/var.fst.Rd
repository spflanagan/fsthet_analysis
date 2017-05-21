\name{var.fst}
\alias{var.fst}
\title{This calculates Cockerham & Weir's Beta. }
\description{
	This calculates Weir & Cockerham (1993)'s Fst. 
	The caluclation is based on variance in allele frequencies.
	See the vignette for details on the calculation of beta.
}
\usage{
var.fst(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{2pbar(1-pbar). This is a single numerical value.}
	\item{fst}{The calculated beta value for this locus.}
}