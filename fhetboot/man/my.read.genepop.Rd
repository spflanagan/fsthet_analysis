\name{my.read.genepop}
\alias{read.genepop}
\title{This reads a genepop file into R}
\description{
  This reads a genepop file into R. 
  It was adapted from a similar functionin adegenet.
}
\usage{
my.read.genepop(file, ncode = 2L, quiet = FALSE)
}
\arguments{
  \item{file}{is the filename of the genpop file.}
  \item{quiet}{If quiet = FALSE updates will be printed. 
  If quiet = T status updates will not be printed.
  }
  \item{ncode}{Do not change this argument.}
}
\value{
	\item{res}{A dataframe with the Population ID in column 1, the Individual ID in column 2, 
	and the genotypes in columns following that. There is one row per individual.}
}
\references{
	\url{ http://adegenet.r-forge.r-project.org/ }
}