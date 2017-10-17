\name{var.fst}
\alias{var.fst}
\title{This calculates Cockerham & Weir's Beta. }
\description{
	This calculates Weir & Cockerham (1993)'s Fst. 
	The caluclation is based on variance in allele frequencies.
	See the vignette for details on the calculation of beta.
}
\usage{
var.fst(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{2pbar(1-pbar). This is a single numerical value.}
	\item{fst}{The calculated beta value for this locus.}
}
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20),
                  loc0=sample(c("0101","0102","0202"),40,replace=TRUE),
                  loc1=sample(c("0101","0102","0202"),40,replace=TRUE))
  var1<-var.fst(gpop,3) 
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.gpop(gfile)
    var1<-var.fst(gpop,3) #calculate variance-based for the first SNP
  }
}