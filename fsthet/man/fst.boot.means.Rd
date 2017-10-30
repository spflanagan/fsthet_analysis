\name{fst.boot.means}
\alias{fst.boot.means}
\title{Calculates mean values within the bins. }
\description{
	This calculates mean heterozygosity and Fst values for each bin used in bootstrapping.
}
\usage{
fst.boot.means(boot.out)
}
\arguments{
  \item{boot.out}{The first item in the output lists from fst.boot (aka boot.out[[1]]. }
}
\value{
	\item{bmu}{A dataframe containing four columns:
		heterozygosity
		Fst
		the number of loci in the bin
		the lower Ht value for the bin
		and the upper Ht value for the bin.}
}
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20))
   for(i in 1:40){
    gpop[1:20,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
    gpop[21:40,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
   }
    fsts<-calc.actual.fst(gpop)
  boot.out<-as.data.frame(t(replicate(1, fst.boot(gpop))))
  outliers<-find.outliers(fsts,boot.out)
\dontrun{
  gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
  gpop<-my.read.genepop(gfile)
  fsts<-calc.actual.fst(gpop)
  boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
  outliers<-find.outliers(fsts,boot.out)
}
}
