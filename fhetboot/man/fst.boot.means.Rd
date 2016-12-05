\name{fst.boot.means}
\alias{fst.boot.means}
\title{Calculates mean values within the bins. }
\description{
	This calculates mean heterozygosity and Fst values for each bin used in bootstrapping.
}
\usage{
fst.boot.means(boot.out[[1]])
}
\arguments{
  \item{boot.out.means}{The first item in the output lists from fst.boot. }
}
\value{
	\item{bmu}{A dataframe containing four columns:
		heterozygosity
		Fst
		the number of loci in the bin
		the lower Ht value for the bin
		and the upper Ht value for the bin.}
}