\name{make.bins}
\alias{make.bins}
\title{This sorts Fst values into a designated number of overlapping heterozygosity bins. }
\description{
	This breaks up Fst values into a designated number of overlapping heterozygosity bins.
	It returns a list containing a data.frame called bins a list called bin.fst with the Fst values for each of the Het categories.
}
\usage{
make.bins(fsts,num.breaks=25, Ht.name="Ht", Fst.name="Fst")
}
\arguments{
  \item{fsts}{A dataframe containing at least the columns with heterozygosity and Fst values.}
  \item{num.breaks}{ The number of breaks used to create bins (default is 25)}
  \item{Ht.name}{Provide the name of the column with the heterozygosity values, unless the column is named "Ht".}
  \item{Fst.name}{Provide the name of the column with the Fst values, unless the column is named "Fst".}

}
\value{
	\item{list(bins,bin.fst)}{A list with a data.frame called bins with the upper and lower Fst and Ht values
	  and a list called bin.fst with the Fst values for each of the Het categories.}
}
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20))
     for(i in 1:40){
      gpop[1:20,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
      gpop[21:40,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
     }
  fsts<-calc.actual.fst(gpop)
  nloci<-(ncol(gpop)-2)
  boot.out<-as.data.frame(t(replicate(nloci, fst.boot.onecol(gpop,"fst"))))
  bins<-make.bins(boot.out,25,Ht.name="V1",Fst.name="V2")
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.gpop(gfile)
    fsts<-calc.actual.fst(gpop)
    nloci<-(ncol(gpop)-2)
    boot.out<-as.data.frame(t(replicate(nloci, fst.boot.onecol(gpop))))
    make.bins(boot.out,25)
  }
}