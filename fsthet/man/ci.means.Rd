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
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20),
    loc0=sample(c("0101","0102","0202"),40,replace=TRUE),
    loc1=sample(c("0101","0102","0202"),40,replace=TRUE))
  quant.out<-fst.boot(gpop, bootstrap = FALSE)
  quant.list<-ci.means(quant.out[[3]])
\dontrun{
  gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
  gpop<-my.read.gpop(gfile)
  quant.out<-fst.boot(gpop, bootstrap = FALSE)
  quant.list<-ci.means(quant.out[[3]])
}
}