
# fsthet

The goal of fsthet is to calculate smoothed quantiles from the existing dataset to identify loci with extreme Fst values relative to their heterozygosity. This allows the user to identify loci that have extreme Fst values without relying on a specific population genetics model. Please see Flanagan and Jones (2017) for details on this issue.

Flanagan SP and Jones AG. 2017. Constraints on the Fst-heterozygosity outlier approach. Journal of Heredity esx048. doi: 10.1093/jhered/esx048

## Installation instructions

To install the latest version, use `devtools::install_github("spflanagan/fsthet_analysis/fsthet")` to download this repository.

## Example

This is a basic example of how to use fsthet. Please see the vignette for more detail.

```{r example}
library(fsthet)
gfile<-system.file("extdata", "example.genepop.txt",package = 'fsthet')
gpop<-my.read.genepop(gfile)
out.dat<-fsthet(gpop)

```
## Issues, bugs, and suggestions

If you run into trouble using fsthet, first ask yourself the following questions:

1. Are my data biallelic SNP data?
2. Do I have many loci in the dataset (say, at least 50 loci)?
3. Do I have many independent populations (at least 10)?

If you answer No to any of those questions, then fsthet is likely not the best program for you to be using. If you answer yes to all of those questiosn and you've still have an issue running your analysis, think you have found a bug in the program, or have a suggestion for improvement, the best thing to do is open up an issue on this github page. Alternatively, you can email me directly (spflanagan.phd@gmail.com). If relevant, please include your dataset (or a subset of it) and reproducible code. 
