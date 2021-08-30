
# bulkATACquality

<!-- badges: start -->
<!-- badges: end -->

The goal of bulkATACquality is to generate diagnostic plots and related statistics for human brain bulk ATAC-seq sample auality, and use a logistic model to predict the quality of sample. An HTML file can be generated to enable users to obtain prediction result, diagnostic plots and other quality-related information of bulk ATAC-seq samples.


## Installation

You can install the released version of bulkATACquality with:

``` r
devtools::install_github("jingyu0602/pqATAC")
```

## Example

This is the main function to generate dignostic plots and statistics

``` r
library(bulkATACquality)
bamfile<-system.file("extdata","ex0001.bam",package = "bulkATACquality", mustWork = TRUE) # read internal data
outdir<-substr(basename(bamfile),1,6) # set the output directory
dir.create(outdir)
bulkATAC(bamfile,outdir) #generate dignostic plots and statistics
```

Downsampling bamfile and generate related plots can be achieved by:
``` r
down(bamfile,outdir)
```

To generate HTML file with prediction result of sample quality, all dignostic plots and statistics:
``` r
generateHTML(outdir) # the directory should be same as above output directory, an HTML file can be generated in this directory
```
