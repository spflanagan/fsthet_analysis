\name{find.outliers}
\alias{find.outliers}
\title{This identifies all of the SNPs outside of the smoothed quantiles in the dataset.}
\description{
	This identifies all of the SNPs outside of the smoothed quantiles in the dataset.
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
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20),
    loc0=sample(c("0101","0102","0202"),40,replace=TRUE),
    loc1=sample(c("0101","0102","0202"),40,replace=TRUE))
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.gpop(gfile)
    fsts<-calc.actual.fst(gpop)
    boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
    outliers<-find.outliers(fsts,boot.out)
  }
}
