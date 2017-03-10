\name{make.bins}
\alias{make.bins}
\title{This sorts Fst values into a designated number of overlapping heterozygosity bins. }
\description{
	This breaks up Fst values into a designated number of overlapping heterozygosity bins.
	It returns a list containing a data.frame called bins a list called bin.fst with the Fst values for each of the Het categories.
}
\usage{
make.bins(fsts,num.breaks=25, Ht.name="Ht", Fst.name="Fst")
}
\arguments{
  \item{fsts}{A dataframe containing at least the columns with heterozygosity and Fst values.}
  \item{num.breaks}{ The number of breaks used to create bins (default is 25)}
  \item{Ht.name}{Provide the name of the column with the heterozygosity values, unless the column is named "Ht".}
  \item{Fst.name}{Provide the name of the column with the Fst values, unless the column is named "Fst".}

}
\value{
	\item{list(bins,bin.fst)}{A list with a data.frame called bins with the upper and lower Fst and Ht values
	  and a list called bin.fst with the Fst values for each of the Het categories.}
}