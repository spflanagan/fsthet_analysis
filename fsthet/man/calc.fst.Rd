\name{calc.fst}
\alias{calc.fst}
\title{This calculates Fst. }
\description{
	This calculates Fst. 
	The caluclation is done as (Ht-Hs)/Ht, 
		where Ht is the expected heterozygosity for all populations 
		and Hs is the expected heterozygosity for each population.
	This calculation is used in bootstrapping functions.
}
\usage{
calc.fst(df, i)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
  \item{i}{Column number containing genotype information.}
}
\value{
	\item{ht}{The expected heterozygosity under Hardy-Weinberg expectations. This is a single numerical value.}
	\item{fst}{The calculated Fst value for this locus.}
}
\examples{
  gpop<-data.frame(popinfo=c(rep("POP 1", 20),rep("POP 2", 20)),ind.names=c(1:20,1:20),
    loc0=sample(c("0101","0102","0202"),40,replace=TRUE))
  fst1<-calc.fst(gpop,3)
  \dontrun{
    gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
    gpop<-my.read.gpop(gfile)
    fst1<-calc.fst(gpop,3) #calculate fst for the first SNP
  }
}