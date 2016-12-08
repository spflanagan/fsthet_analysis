\name{ci.df}
\docType{data}
\alias{ci.df}
\title{Example dataframe of confidence intervals from bootstrapping.}
\description{
   Example dataframe of confidence intervals from bootstrapping.
   This file contains a dataframe with 22 columns and 2 rows. The rownames are the Ht x-value for plotting. 
   The columns are the upper and lower confidence intervals.
}
\usage{ci.df}
\format{data.frame}
\source{
	ci.df<-data.frame(do.call(cbind(boot.out[[3]])))
}
\references{
   See Flanagan & Jones
}
