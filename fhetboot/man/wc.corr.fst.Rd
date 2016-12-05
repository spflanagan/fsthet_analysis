\name{wc.corr.fst}
\alias{wc.corr.fst}
\title{This calculates Beaumont & Nichols's Fst. }
\description{
	This calculates Beaumont & Nichols (1996)'s Fst. 
	This is just a sample-size corrected version of the Weir & Cockerham (1993)'s beta.
	The caluclation is done as beta=(q2-q3)/(1-q3), 
		where:
		q3=2Y/(N(N-1)),
		q2=x0/N,
		N=number of populations
		x0=sum(sum((cij*cij)-nj)/(nj(nj-1)))
			where cij is the number of allele i (allele count i) in group j
			and nj is the sample size of population j
		Y=sum((sum(cij*cik)/(nj*nk)))
			where cij is allele count i in population j and cik is allele count i in population k
			and nj is the sample size in pop j and nk is the sample size in pop k.
			
	This calculation is used in bootstrapping functions.
}
\usage{
wc.corr.fst(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{1-q3. This is a single numerical value.}
	\item{fst}{The calculated Fst value ((q2-q3)/(1-q3))for this locus.}
}