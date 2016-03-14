\name{find.outliers}
\alias{find.outliers}
\title{This identifies all of the SNPs outside of the confidence intervals in the dataset.}
\description{
	This identifies all of the SNPs outside of the confidence intervals in the dataset.
}
\usage{
find.outliers(df, boot.out, file.name="Test")
}
\arguments{
  \item{df}{Provide the dataframe with Ht and Fst values.}
  \item{boot.out}{Bootstrap output. You must either provide this or a list of confidence interval values.}
  \item{ci.list}{List of confidence intervals. You must either provide this or bootstrap output.}
  \item{file.name}{You may provide a file name to output the outliers to a csv file. Otherwise, the function will only return the 95% and 99% outliers.}
}
\value{
	\item{out95}{A list of the 95% outlier loci}
	\item{out99}{A list of the 99% outlier loci}
}