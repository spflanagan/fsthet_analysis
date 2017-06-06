\name{calc.actual.fst}
\alias{calc.actual.fst}
\title{This calcualtes global Fsts from a genepop dataframe.}
\description{
	This calcualtes global Fsts from a genepop dataframe.
	This does not include bootstrapping.
}
\usage{
calc.actual.fst(df, fst.choice="fst")
}
\arguments{
  \item{df}{Provide the genepop dataframe (from my.read.genepop).}
  \item{fst.choice}{Specify which type of fst calculation should be used. See fst.options.print for the choices.}
}
\value{
	\item{fsts}{This returns a dataframe with Locus, Ht, and Fst characters.}
}
\examples{
\dontrun{gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
gpop<-my.read.genepop(gfile)
fsts<-calc.actual.fst(gpop)}
}