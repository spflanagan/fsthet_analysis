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
\dontrun{
  gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
  gpop<-my.read.gpop(gfile)
  fsts<-calc.actual.fst(gpop)
  boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
  boot.pvals<-p.boot(fsts,boot.out=boot.out)
}
}