\name{fst.boot.onecol}
\alias{fst.boot.onecol}
\title{This bootstraps across all individuals to calculate a bootstrapped Fst for a randomly-sampled locus. }
\description{
	This calculates Fst using calc.fst.
	It randomly selects a column containing genotype information for all individuals.
	It then calculates Fst and Ht for that locus.
}
\usage{
fst.boot.onecol(df)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
}
\value{
	\item{ht.fst}{A vector containin Ht and Fst}
}