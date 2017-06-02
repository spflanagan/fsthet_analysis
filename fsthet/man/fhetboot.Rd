\name{fhetboot}
\alias{fhetboot}
\title{This is a wrapper to run the bootstrapping and plot the confidence intervals and significant loci.}
\description{
	This calcualtes global Fsts from a genepop dataframe and then does:
		p-value calculations
		plots the Heterozygosity-Fst relationship with smoothed CIs
		outputs the loci lying outside the confidence intervals.
	Returns a data frame containing Locus ID, Ht, Fst, P-value, a Benjamini-Hochberg-corrected P-value, and a true/false value of whether it's an outlier.
}
\usage{
	fhetboot(gpop, fst.choice="fst", alpha=0.05,nreps=10)
}
\arguments{
  \item{gpop}{Provide the genepop dataframe (from my.read.genepop).}
  \item{fst.choice}{Specify which type of fst calculation should be used. See fst.options.print for the choices.}
  \item{alpha}{The alpha value for the confidence intervals and the p-value adjustment calculations (default is 0.05).}
  \item{nreps}{The number of bootstrap replicates to use. The default is 10.}
}
\value{
	\item{fsts}{This returns a dataframe with Locus, Ht, Fst, P-value, correcte P-value, and True/False of whether it's an outlier.}
}