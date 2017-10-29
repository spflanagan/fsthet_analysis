\name{fst.boot.onecol}
\alias{fst.boot.onecol}
\title{This bootstraps across all individuals to calculate a bootstrapped Fst for a randomly-sampled locus. }
\description{
	This calculates Fst using calc.fst.
	It randomly selects a column containing genotype information for all individuals.
	It then calculates Fst and Ht for that locus.
}
\usage{
fst.boot.onecol(df, fst.choice)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{fst.choice}{
	A character defining which fst calculation is to be used. 
	The three options are:
		Wright's Fst (Wright, wright, WRIGHT, W, w)
		Weir and Cockerham 1993's beta (WeirCockerham,weircockerham,wc,WC)
		Corrected Weir and Cockerham 1993's beta from Beaumont and Nichols 1996 (WeirCockerhamCorrected, weircockerhamcorrected,corrected,wcc,WCC)}
}
\value{
	\item{ht.fst}{A vector containin Ht and Fst}
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
\dontrun{
  gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
  gpop<-my.read.genepop(gfile)
  fsts<-calc.actual.fst(gpop)
  nloci<-(ncol(gpop)-2)
  boot.out<-as.data.frame(t(replicate(nloci, fst.boot.onecol(gpop,"fst"))))
  outliers<-find.outliers(fsts,boot.out)
}
}
