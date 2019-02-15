---
output:
  md_document:
    variant: markdown_github
---


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
