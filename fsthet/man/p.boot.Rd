\name{p.boot}
\alias{p.boot}
\title{Calculates mean values within the bins. }
\description{
	This calculates mean heterozygosity and Fst values for each bin used in bootstrapping.
}
\usage{
p.boot(actual.fsts, boot.out,boot.means=NULL)
}
\arguments{
  \item{actual.fsts}{The first item in the output lists from fst.boot. }
  \item{boot.out}{The output from a bootstrapping run. Either supply this or boot.means. }
  \item{boot.means}{The output from fst.boot.means. Either supply this or bootstrapping output. }
}
\value{
	\item{pvals}{A numeric containing uncorrected p-values for each locus. The names attribute are the locus names.}
}
\examples{
  \dontshow{
    gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20))
     for(i in 1:40){
      gpop[1:20,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
      gpop[21:40,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
     }
    fsts<-calc.actual.fst(gpop)
    boot.out<-as.data.frame(t(replicate(1, fst.boot(gpop, num.breaks=20))))
    boot.pvals<-p.boot(fsts,boot.out=boot.out)
  }
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.genepop(gfile)
    fsts<-calc.actual.fst(gpop)
    boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
    boot.pvals<-p.boot(fsts,boot.out=boot.out)
  }
}