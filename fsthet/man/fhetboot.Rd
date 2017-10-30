\name{fhetboot}
\alias{fhetboot}
\title{This is a wrapper to run the bootstrapping and plot the confidence intervals and significant loci.}
\description{
	This calcualtes global Fsts from a genepop dataframe and then does:
		p-value calculations
		plots the Heterozygosity-Fst relationship with smoothed CIs
		outputs the loci lying outside the confidence intervals.
	Returns a data frame containing Locus ID, Ht, Fst, P-value, a Benjamini-Hochberg-corrected P-value, and a true/false value of whether it's an outlier.
}
\usage{
	fhetboot(gpop, fst.choice="fst", alpha=0.05,nreps=10)
}
\arguments{
  \item{gpop}{Provide the genepop dataframe (from my.read.genepop).}
  \item{fst.choice}{Specify which type of fst calculation should be used. See fst.options.print for the choices.}
  \item{alpha}{The alpha value for the confidence intervals and the p-value adjustment calculations (default is 0.05).}
  \item{nreps}{The number of bootstrap replicates to use. The default is 10.}
}
\value{
	\item{fsts}{This returns a dataframe with Locus, Ht, Fst, P-value, correcte P-value, and True/False of whether it's an outlier.}
}
\examples{
  \dontshow{
   gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20))
   for(i in 1:40){
    gpop[1:20,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
    gpop[21:40,(i+2)]<-sample(c("0101","0102","0202"),20,replace=TRUE)
   }
   out.dat<-fhetboot(gpop, fst.choice="fst", alpha=0.05,nreps=1)
  }
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.genepop(gfile)
    out.dat<-fhetboot(gpop, fst.choice="fst", alpha=0.05,nreps=10)
  }
}