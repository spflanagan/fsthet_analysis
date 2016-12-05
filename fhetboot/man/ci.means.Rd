\name{ci.means}
\alias{ci.means}
\title{This calculates the average confidence intervals from multiple bootstrap outputs. }
\description{
	This calculates the mean upper and lower confidence intervals from a list of bootstrap CI matrices.
}
\usage{
ci.means(boot.out.list)
}
\arguments{
  \item{boot.out.list}{A list of matrices. Each matrix is the CIs from fst.boot (boot.out[[3]]).}
}
\value{
	\item{avg.cil}{A list of the average lower CI values}
	\item{avg.ciu}{A list of the average upper CI values}
}