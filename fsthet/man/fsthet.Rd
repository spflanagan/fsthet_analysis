\name{fsthet}
\alias{fsthet}
\title{This is a wrapper to generate and plot the smoothed quantiles and identify outliers.}
\description{
	This calcualtes global Fsts from a genepop dataframe and then does:
		calculates smoothed quantiles
		plots the Heterozygosity-Fst relationship with smoothed quantiles
		outputs the loci lying outside the quantiles.
	Returns a data frame containing Locus ID, Ht, Fst, and a true/false value of whether it's an outlier.
}
\usage{
	fsthet(gpop, fst.choice="wright", alpha)
}
\arguments{
  \item{gpop}{Provide the genepop dataframe (from my.read.genepop).}
  \item{fst.choice}{Specify which type of fst calculation should be used. See fst.options.print for the choices.}
  \item{alpha}{The alpha value for the quantiles (default is 0.05 to generate 95 percent quantiles).}
}
\value{
	\item{fsts}{This returns a dataframe with Locus, Ht, Fst, and True/False of whether it's an outlier.}
}