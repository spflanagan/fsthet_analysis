\name{wc.fst}
\alias{wc.fst}
\title{This calculates Weir & Cockerham's Fst. }
\description{
	This calculates Weir & Cockerham (1993)'s Fst. 
	The caluclation is done as beta=(F0-F1)/(1-F1), 
		where:
		F0=(M*X-N)/((M-1)*N),
		F1=(Y-X)/(N*(N-1)),
		N=number of populations
		M=average number of individuals per population
		X=sum(sum(pij^2)
		Y=sum((sum(pij)^2)
			where pij is the frequency of allele i in group j
	This calculation is used in bootstrapping functions.
}
\usage{
wc.fst(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{1-F1. This is a single numerical value.}
	\item{fst}{The calculated Fst value ((F0-F1)/(1-F1))for this locus.}
}