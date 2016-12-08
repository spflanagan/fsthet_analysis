\name{fst.boot.means}
\alias{fst.boot.means}
\title{Calculates mean values within the bins. }
\description{
	This calculates mean heterozygosity and Fst values for each bin used in bootstrapping.
}
\usage{
fst.boot.means(boot.out)
}
\arguments{
  \item{boot.out}{The first item in the output lists from fst.boot (aka boot.out[[1]]. }
}
\value{
	\item{bmu}{A dataframe containing four columns:
		heterozygosity
		Fst
		the number of loci in the bin
		the lower Ht value for the bin
		and the upper Ht value for the bin.}
}