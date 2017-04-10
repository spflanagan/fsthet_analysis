\name{bins}
\docType{data}
\alias{bins}
\title{Example heterozygosity bins from fsthet.}
\description{
   This is a list with a data.frame of bins (the lower and upper bounds for each heterozygosity bin) and a list of fsts that fall into each bin, with the name of each set of Fst values being the upper heterozygosity bound from the data.frame of bins.
}
\usage{bins}
\format{list}
\source{
	bins<-make.bins(fsts)
}
\references{
   See Flanagan & Jones
}
