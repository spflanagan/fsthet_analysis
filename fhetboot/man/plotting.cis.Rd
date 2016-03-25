\name{plotting.cis}
\alias{plotting.cis}
\title{This plots a dataframe of fsts with bootstrapped confidence intervals.}
\description{
	This plots a dataframe of fsts with bootstrapped confidence intervals.
}
\usage{
plotting.cis(df, boot.out, ci.list, Ht.name, Fst.name, ci.col, pt.pch,file.name,make.file)
}
\arguments{
  \item{df}{A dataframe of Fst and Ht values. It must have at least two columns, one named "Ht" and one named "Fst". Or you must pass the column names to the function}
  \item{boot.out}{Bootstrap output. You must either provide this or a list of confidence interval values.}
  \item{ci.list}{List of confidence intervals. You must either provide this or bootstrap output.}
  \item{Ht.name}{Provide the name of the column with the heterozygosity values, unless the column is named "Ht".}
  \item{Fst.name}{Provide the name of the column with the Fst values, unless the column is named "Fst".}
  \item{ci.col}{You can input the colors of the confidence intervals to be plotted. First is the 95 percent CI, second is the 99 percent CI. Defaults are "red" and "gold".}
  \item{pt.pch}{You can change the point shape here. Default is 1 (open circles)}
  \item{file.name}{You can provide the filename. If not provided, default is "OutlierLoci" in the current directory.}
  \item{make.file}{A boolean value (TRUE or FALSE). If TRUE, a file will be created with the plot. If FALSE, the plot will be made in R only (and can be further annotated).}
}
