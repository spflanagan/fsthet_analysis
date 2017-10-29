\name{find.quantiles}
\alias{find.quantiles}
\title{Generates quantiles from binned Fst values}
\description{
	This takes the output from make.bins and calculates the smoothed quantiles.
}
\usage{
find.quantiles(bins,bin.fst,ci=0.05)
}
\arguments{
  \item{bins}{A dataframe containing with upper and lower het and Fst values for each bin (output from make.bins).}
  \item{bin.fst}{ A list with the Fst values for each bin (output from make.bins). }
  \item{ci}{A value for the confidence intervals alpha (default is 0.05).}
}
\value{
	\item{fst.CI}{A list of data.frames, one for each ci value with the upper and lower Fst quantiles for each Heterozygosity bin.}
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
  fst.CI<-find.quantiles(bins$bins,bins$bin.fst)

}