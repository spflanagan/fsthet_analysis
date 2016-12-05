\name{fst.boot}
\alias{fst.boot}
\title{This is the major bootstrapping function to calculate confidence intervals. }
\description{
	This randomly samples all of the loci, with replacement 
		(so if you have 200 loci, it will choose 200 loci to calculate Fst for, but some may be sampled multiply)
	It makes use of fst.boot.onerow.
	To calculate the confidence intervals, this function bins the Fst values based on heterozygosity values.
	The bins are overlapping and each bin is the width of smooth.rate. 
	The Fst value which separates the top 100*(ci/2) and bottom 100*(ci/2) percent in each bin are the upper and lower CIs.
	This function can be slow. We recommend running it 10 times to generate confidence intervals for analysis.
}
\usage{
fst.boot(df)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{fst.choice}{
	A character defining which fst calculation is to be used. 
	The three options are:
		Nei's Fst (nei,Nei,NEI,N)
		Weir and Cockerham 1993's beta (WeirCockerham,weircockerham,wc,WC)
		Corrected Weir and Cockerham 1993's beta from Beaumont and Nichols 1996 (WeirCockerhamCorrected, weircockerhamcorrected,corrected,wcc,WCC)
	Default is Nei's.}
  \item{ci}{A value for the confidence intervals alpha (default is 0.05).}
  \item{smooth.rate}{The distance between the start of bins (default is 0.1)}
}
\value{
	\item{Fsts}{The bootstrapped Fst and Ht values}
	\item{Bins}{A dataframe containing the bins start and stop Ht values.}
	\item{fst.CI}{A list of dataframes containing the lower and upper confidence intervals' Ht values.}
}