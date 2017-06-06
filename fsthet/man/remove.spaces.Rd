\name{remove.spaces}
\alias{remove.spaces}
\title{This removes spaces from a character vector}
\description{
  This removes spaces from a before and after words in a character vector. 
  It was adapted from a similar function in adegenet.
}
\usage{
remove.spaces(charvec)
}
\arguments{
  \item{charvec}{is a vector of characters containing spaces to be removed.}
}
\value{
	\item{charvec}{A vector of characters without spaces}
}
\references{
	\url{ http://adegenet.r-forge.r-project.org/ }
}
\examples{
charvec<-c("this ", " is"," a"," test")
remove.spaces(charvec)
}