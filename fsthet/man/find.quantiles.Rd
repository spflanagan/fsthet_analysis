\name{find.quantiles}
\alias{find.quantiles}
\title{Generates quantiles from binned Fst values}
\description{
	This takes the output from make.bins and calculates the smoothed quantiles.
}
\usage{
find.quantiles(bins,bin.fst,ci=0.05)
}
\arguments{
  \item{bins}{A dataframe containing with upper and lower het and Fst values for each bin (output from make.bins).}
  \item{bin.fst}{ A list with the Fst values for each bin (output from make.bins). }
  \item{ci}{A value for the confidence intervals alpha (default is 0.05).}
}
\value{
	\item{fst.CI}{A list of data.frames, one for each ci value with the upper and lower Fst quantiles for each Heterozygosity bin.}
}
\examples{
\dontrun{
load(bins)
fst.CI<-find.quantiles(bins$bins,bins$bin.fst)
}
}