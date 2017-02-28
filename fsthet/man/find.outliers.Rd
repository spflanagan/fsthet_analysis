\name{find.outliers}
\alias{find.outliers}
\title{This identifies all of the SNPs outside of the confidence intervals in the dataset.}
\description{
	This identifies all of the SNPs outside of the confidence intervals in the dataset.
}
\usage{
find.outliers(df, boot.out, ci.df = NULL, file.name = NULL)
}
\arguments{
  \item{df}{Provide the dataframe with Ht and Fst values.}
  \item{boot.out}{Bootstrap output. You must provide this.}
  \item{ci.df}{List of confidence intervals. You may provide this in addition to bootstrap output to save a small amount of time.}
  \item{file.name}{You may provide a file name to output the outliers to a csv file. Otherwise, the function will only return the outliers.}
}
\value{
	\item{out}{A list of the outlier loci}
}