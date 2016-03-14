\name{mean.cis}
\alias{mean.cis}
\title{This calculates the average confidence intervals from multiple bootstrap outputs. }
\description{
	This calculates the mean upper and lower confidence intervals from a list of bootstrap CI matrices.
}
\usage{
mean.cis(boot.out.list)
}
\arguments{
  \item{boot.out.list}{A list of matrices. Each matrix is either the 95% or 99% CIs from fst.boot.}
}
\value{
	\item{avg.cil}{A list of the average lower CI values}
	\item{avg.ciu}{A list of the average upper CI values}
}