\name{fst.boot}
\alias{fst.boot}
\title{This is the major bootstrapping function to calculate confidence intervals. }
\description{
	This randomly samples all of the loci, with replacement 
		(so if you have 200 loci, it will choose 200 loci to calculate Fst for, but some may be sampled multiply)
	It makes use of fst.boot.onerow.
	To calculate the confidence intervals, this function bins the Fst values based on heterozygosity values.
	The bins are overlapping and each bin is the width of smooth.rate. 
	The Fst value which separates the top 100*(ci/2) and bottom 100*(ci/2) percent in each bin are the upper and lower CIs.
	This function can be slow. We recommend running it 10 times to generate confidence intervals for analysis.
}
\usage{
fst.boot(df,fst.choice="fst",ci=0.05,num.breaks=25, bootstrap = TRUE,min.per.bin=20)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{fst.choice}{
	A character defining which fst calculation is to be used. See fst.options.print() for the choices.}
  \item{ci}{A value for the confidence intervals alpha (default is 0.05).}
  \item{num.breaks}{The number of breaks used to create bins (default is 25)}
  \item{bootstrap}{A TRUE/FALSE statement telling the program whether to bootstrap and then determine the bins or to calculate bins and confidence intervals from the empirical dataset without bootstrapping. The default is TRUE, which means bootstrapping occurs.}
  \item{min.per.bin}{The minimum number of loci that are required for a bin to be retained. Default is 20.}
}
\value{
	\item{Fsts}{The bootstrapped Fst and Ht values}
	\item{Bins}{A dataframe containing the bins start and stop Ht values.}
	\item{fst.CI}{A list of dataframes containing the lower and upper confidence intervals' Ht values.}
}
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20))
   for(i in 1:40){
    gpop[1:20,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
    gpop[21:40,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
   }
  fsts<-calc.actual.fst(gpop)
  quant.out<-as.data.frame(t(replicate(1, fst.boot(gpop,bootstrap=FALSE))))
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.genepop(gfile)
    fsts<-calc.actual.fst(gpop)
    quant.out<-as.data.frame(t(replicate(1, fst.boot(gpop,bootstrap=FALSE))))
  }
}