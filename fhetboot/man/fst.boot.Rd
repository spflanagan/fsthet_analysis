\name{fst.boot}
\alias{fst.boot}
\title{This is the major bootstrapping function to calculate confidence intervals. }
\description{
	This randomly samples all of the loci, with replacement 
		(so if you have 200 loci, it will choose 200 loci to calculate Fst for, but some may be sampled multiply)
	It makes use of fst.boot.onerow.
	To calculate the confidence intervals, this function bins the Fst values based on heterozygosity values.
	The bins are overlapping and each bin is 0.1 Ht wide. 
	The Fst value which separates the top 2.5 and bottom 2.5 percent in each bin are the upper and lower 95 percent CIs.
	The Fst value which separates the top 0.5 and bottom 0.5 percent in each bin are the upper and lower 99 percent CIs.
	This function can be slow. We recommend running it 10 times to generate confidence intervals for analysis.
}
\usage{
fst.boot(df)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
}
\value{
	\item{Fsts}{The bootstrapped Fst and Ht values}
	\item{CI95}{A matrix with the lower 95 CI in column 1 and the upper 95 percent CI in column 2.}
	\item{CI99}{A matrix with the lower 99 CI in column 1 and the upper 99 percent CI in column 2.}
}