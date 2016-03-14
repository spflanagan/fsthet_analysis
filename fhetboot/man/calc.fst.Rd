\name{calc.fst}
\alias{calc.fst}
\title{This calculates Fst. }
\description{
	This calculates Fst. 
	The caluclation is done as (Ht-Hs)/Ht, 
		where Ht is the expected heterozygosity for all populations 
		and Hs is the expected heterozygosity for each population.
	This calculation is used in bootstrapping functions.
}
\usage{
calc.fst(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{The expected heterozygosity under Hardy-Weinberg expectations. This is a single numerical value.}
	\item{fst}{The calculated Fst value for this locus.}
}
\examples{
	gpop<-read.my.genepop(example)
	fst<-calc.fst(gpop,3)
}