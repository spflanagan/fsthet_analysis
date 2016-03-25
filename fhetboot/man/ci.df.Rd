\name{ci.df}
\docType{data}
\alias{ci.df}
\title{Example dataframe of confidence intervals from bootstrapping.}
\description{
   Example dataframe of confidence intervals from bootstrapping.
   This file contains a dataframe with 22 columns and 4 rows. The rownames are the Ht x-value for plotting. 
   The columns are the upper and lower 95 and 99 percent confidence intervals in the following order: 
   low95, upp95, low99, upp99.
}
\usage{ci.df}
\format{data.frame}
\source{
	ci.df<-data.frame(low95=ci.means(boot.out[[2]])[[1]],upp95=ci.means(boot.out[[2]])[[2]],
		low99=ci.means(boot.out[[3]])[[1]],upp99=ci.means(boot.out[[3]])[[1]])
}
\references{
   See Flanagan & Jones
}
