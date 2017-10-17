\name{plotting.cis}
\alias{plotting.cis}
\title{This plots a dataframe of fsts with bootstrapped confidence intervals.}
\description{
	This plots a dataframe of fsts with bootstrapped confidence intervals.
}
\usage{
plotting.cis(df,boot.out,ci.df=NULL,sig.list=NULL,Ht.name="Ht",Fst.name="Fst",
	ci.col="red", pt.pch=1,file.name=NULL,sig.col=ci.col,make.file=TRUE)
}
\arguments{
  \item{df}{A dataframe of Fst and Ht values. It must have at least two columns, one named "Ht" and one named "Fst". Or you must pass the column names to the function}
  \item{boot.out}{Bootstrap output. You must either provide this or a list of confidence interval values.}
  \item{ci.df}{Data frame of confidence intervals. You must either provide this or bootstrap output.}
  \item{sig.list}{List of significant locus names (this acts as a way to highlight particular loci). This is optional and colors some of the points using the same shape as pt.pch and the color of sig.col (default sig.color is same as ci.col).}
  \item{Ht.name}{Provide the name of the column with the heterozygosity values, unless the column is named "Ht".}
  \item{Fst.name}{Provide the name of the column with the Fst values, unless the column is named "Fst".}
  \item{ci.col}{You can input the colors of the confidence intervals to be plotted. First is the 95 percent CI, second is the 99 percent CI. Defaults are "red" and "gold".}
  \item{pt.pch}{You can change the point shape here. Default is 1 (open circles)}
  \item{sig.col}{The color of the significant loci, if that option is taken. The default is the same color as the confidence interval.}
  \item{file.name}{You can provide the filename. If not provided, default is "OutlierLoci" in the current directory.}
  \item{make.file}{A boolean value (TRUE or FALSE). If TRUE, a file will be created with the plot. If FALSE, the plot will be made in R only (and can be further annotated).}
}
\examples{
 gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20),
                  loc0=sample(c("0101","0102","0202"),40,replace=TRUE),
                  loc1=sample(c("0101","0102","0202"),40,replace=TRUE))
  fsts<-calc.actual.fst(gpop)
  bins<-make.bins(fsts)
  cis<-find.quantiles(bins = bins$bins,bin.fst = bins$bin.fst)
  quant.list<-cis$CI0.95
  plotting.cis(df=fsts,ci.df=quant.list,make.file=F)
  \dontrun{
  load(fsts)
  bins<-make.bins(fsts)
  cis<-find.quantiles(bins = bins$bins,bin.fst = bins$bin.fst)
  quant.list<-cis$CI0.95
  plotting.cis(df=fsts,ci.df=quant.list,make.file=F)
  }
}