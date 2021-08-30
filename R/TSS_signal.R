#' Generate TSS signal plot of nucleosome-bound region and nucleosome-free region
#' Calculate signal difference between the highest peak and the second highest peak
#' Calculate cov of signal
#'
#' @param gal1 GAl object
#' @param bamfile path of filtered bamfile
#' @param dir output directory
#' @return signal difference between the highest peak and the second highest peak and cov of signal
#' @import ChIPpeakAnno
#' @import GenomeInfoDb
#' @import GenomicFeatures
#' @import graphics
#' @import stats
#' @export
#'
#' @examples
#' \dontrun{
#' TSSsignal_result<-TSS_signal(bamfile,gAl,dir0)
#' TSSsignal_cov<-TSSsignal_result$cov_tsssignal
#' TSSsignal_difference<-TSSsignal_result$TSS_diff
#' }


# split reads
## run program for chromosome 1 only
TSS_signal<-function(bamfile,gal1,dir){
  genome <- Hsapiens
  txs <- bulkATACquality:::txs_chr1
  outPath <- paste0("splited_",strsplit(basename(bamfile),"_")[[1]][1])
  ## split the reads into NucleosomeFree, mononucleosome,
  ## dinucleosome and trinucleosome.
  ## and save the binned alignments into bam files.
  objs <- ATACseqQC::splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)

  bamfiles <- file.path(outPath,
                        c("NucleosomeFree.bam",
                          "mononucleosome.bam",
                          "dinucleosome.bam",
                          "trinucleosome.bam"))

  TSS <- promoters(txs, upstream=0, downstream=1)
  TSS <- unique(TSS)
  (librarySize <- ChIPpeakAnno::estLibSize(bamfiles))

  ## calculate the signals around TSSs.
  NTILE <- 101
  dws <- ups <- 1010
  seqlev <- "chr1" ## subsample data for quick run
  sigs <- ATACseqQC::enrichedFragments(gal=objs[c("NucleosomeFree",
                                       "mononucleosome",
                                       "dinucleosome",
                                       "trinucleosome")],
                            TSS=TSS,
                            librarySize=librarySize,
                            seqlev=seqlev,
                            TSS.filter=0.5,
                            n.tile = NTILE,
                            upstream = ups,
                            downstream = dws)

  ## get signals normalized for nucleosome-free and nucleosome-bound regions.
  out <- ChIPpeakAnno::featureAlignedDistribution(sigs,
                                    reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt=.5, n.tile=NTILE, type="l",
                                    ylab="Averaged coverage")
  ## rescale the nucleosome-free and nucleosome signals to 0~1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  out <- apply(out, 2, range01)

  a<-out[,1]
  b<-c(a[2:101],a[101])
  d<-c((b-a)[2:101],0)
  data_TSSsignal<-data.frame(a,b,b-a,d)
  colnames(data_TSSsignal)<-c("a","b","diff","p_n")
  data_TSSsignal_peak<-data_TSSsignal[(data_TSSsignal$p_n)*(data_TSSsignal$diff)<0,]
  maxTSS<-max(data_TSSsignal_peak$b)
  sec_maxTSS<-max(data_TSSsignal_peak$b[data_TSSsignal_peak$b!=max(data_TSSsignal_peak$b)])
  TSS_diff<-maxTSS-sec_maxTSS
  print(paste0("7./The standardized difference of highest and second highest peak is ", TSS_diff))
  #TSS_Diff<-c(TSS_Diff,TSS_diff)

  png(paste0(dir,"_TSS_signal_plot.png"),res=150, width = 2400,height = 1200)
  matplot(out, type="l", xaxt="n",
          xlab="Position (bp)",
          ylab="Fraction of signal")
  axis(1, at=seq(0, 100, by=10)+1,
       labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
  text(500,0.8,paste0("The standadized difference between 2 peaks is ", TSS_diff))
  abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
  abline(h=maxTSS,col="red",lwd=2)
  abline(h=sec_maxTSS,col="orange",lwd=2)
  dev.off()

  # calculate the coefficient of variation
  cov_tsssignal<-sd(out[,1])/mean(out[,1])
  #cov_TSSsignal<-c(cov_TSSsignal,cov_tsssignal)
  print(paste0("8./The cov of TSSsignal is ", cov_tsssignal))

  return(list(TSS_diff=TSS_diff,cov_tsssignal=cov_tsssignal))
}
