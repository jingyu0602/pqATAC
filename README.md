
# bulkATACquality

<!-- badges: start -->
<!-- badges: end -->

The goal of bulkATACquality is to generate diagnostic plots and related statistics for human brain bulk ATAC-seq sample auality, and use a logistic model to predict the quality of sample. An HTML file can be generated to enable users to obtain prediction result, diagnostic plots and other quality-related information of bulk ATAC-seq samples.


## Installation

You can install the released version of bulkATACquality from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bulkATACquality")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bulkATACquality)
bamfile<-system.file("extdata","ex0001.bam",package = "bulkATACquality", mustWork = TRUE)
outdir<-substr(basename(bamfile),1,6)
dir.create(outdir)
bulkATAC(bamfile,outdir)
## basic example code
```

