\name{remove.spaces}
\alias{remove.spaces}
\title{This removes spaces from a character vector}
\description{
  This removes spaces from a character vector. 
  It was adapted from a similar functionin adegenet.
}
\usage{
remove.spaces(charvec)
}
\arguments{
  \item{charvec}{is a vector of characters containing spaces to be removed.}
}
\value{
	\item{charvec}{A vector of characters without spaces}
}
\references{
	\url{ http://adegenet.r-forge.r-project.org/ }
}


\name{my.read.genepop}
\alias{read.genepop}
\title{This reads a genepop file into R}
\description{
  This reads a genepop file into R. 
  It was adapted from a similar functionin adegenet.
}
\usage{
my.read.genepop(file, ncode = 2L, quiet = FALSE)
}
\arguments{
  \item{file}{is the filename of the genpop file.}
  \item{quiet}{If quiet = FALSE updates will be printed. 
  If quiet = T status updates will not be printed.
  }
  \item{ncode}{Do not change this argument.}
}
\value{
	\item{res}{A dataframe with the Population ID in column 1, the Individual ID in column 2, 
	and the genotypes in columns following that. There is one row per individual.}
}
\references{
	\url{ http://adegenet.r-forge.r-project.org/ }
}

\name{calc.allele.freq}
\alias{calc.allele.freq}
\title{This calculates allele frequencies. }
\description{
  This calculates allele frequencies from a list of genotypes.
}
\usage{
calc.allele.freq(genotypes)
}
\arguments{
  \item{genotypes}{A list of genotypes.}
}
\value{
	\item{obs.af}{A list of observed allele frequencies in the genotypes list.}
}

\name{calc.exp.het}
\alias{calc.exp.het}
\title{This calculates expected heterozygosities. }
\description{
  This calculates expected heterozygosities from a list of allele frequencies.
}
\usage{
calc.exp.het(af)
}
\arguments{
  \item{af}{is a list of allele frequencies.}
}
\value{
	\item{ht}{The expected heterozygosity under Hardy-Weinberg expectations. This is a single numerical value.}
}


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
	gpop<-read.my.genepop(example)
	fst<-calc.fst(gpop,3)
}

\name{fst.boot.onecol}
\alias{fst.boot.onecol}
\title{This bootstraps across all individuals to calculate a bootstrapped Fst for a randomly-sampled locus. }
\description{
	This calculates Fst using calc.fst.
	It randomly selects a column containing genotype information for all individuals.
	It then calculates Fst and Ht for that locus.
}
\usage{
fst.boot.onecol(df)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
}
\value{
	\item{ht.fst}{A vector containin Ht and Fst}
}

\name{fst.boot}
\alias{fst.boot}
\title{This is the major bootstrapping function to calculate confidence intervals. }
\description{
	This randomly samples all of the loci, with replacement 
		(so if you have 200 loci, it will choose 200 loci to calculate Fst for, but some may be sampled multiply)
	It makes use of fst.boot.onerow.
	To calculate the confidence intervals, this function bins the Fst values based on heterozygosity values.
	The bins are overlapping and each bin is 0.1 Ht wide. 
	The Fst value which separates the top 2.5% and bottom 2.5% in each bin are the upper and lower 95% CIs.
	The Fst value which separates the top 0.5% and bottom 0.5% in each bin are the upper and lower 99% CIs.
	This function can be slow. We recommend running it 10 times to generate confidence intervals for analysis.
}
\usage{
fst.boot(df)
}
\arguments{
  \item{df}{A dataframe containing the genepop information, where the first column is the population ID.}
}
\value{
	\item{Fsts}{The bootstrapped Fst and Ht values}
	\item{CI95}{A matrix with the lower 95% CI in column 1 and the upper 95% CI in column 2.}
	\item{CI99}{A matrix with the lower 99% CI in column 1 and the upper 99% CI in column 2.}
}

\name{mean.cis}
\alias{mean.cis}
\title{This calculates the average confidence intervals from multiple bootstrap outputs. }
\description{
	This calculates the mean upper and lower confidence intervals from a list of bootstrap CI matrices.
}
\usage{
mean.cis(boot.out.list)
}
\arguments{
  \item{boot.out.list}{A list of matrices. Each matrix is either the 95% or 99% CIs from fst.boot.}
}
\value{
	\item{avg.cil}{A list of the average lower CI values}
	\item{avg.ciu}{A list of the average upper CI values}
}


\name{plot.cis}
\alias{plot.cis}
\title{This plots }
\description{
	This calculates the mean upper and lower confidence intervals from a list of bootstrap CI matrices.
}
\usage{
mean.cis(boot.out.list)
}
\arguments{
  \item{boot.out.list}{A list of matrices. Each matrix is either the 95% or 99% CIs from fst.boot.}
}
\value{
	\item{avg.cil}{A list of the average lower CI values}
	\item{avg.ciu}{A list of the average upper CI values}
}

Here is the R manual page for one the datasets:

\name{bcdat.rad}
\docType{data}
\alias{bcdat.rad}
\title{Title string for bcdat.rad}
\description{
   This is the description of the bcdat.rad dataset.
}
\usage{bcdat.rad}
\format{This describes the format of bcdat.rad.}
\source{This gives the source of bcdat.rad}
\references{
   Any references for bcdat.rad go here
}