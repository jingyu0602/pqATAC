#' shifted all reads in input bamfile for downstream analysis, such as peak-calling and footprinting

#' The function shiftGAlignmentsList can be used to shift the reads. By default, all reads aligning to the positive strand are offset by +4bp,
#' and all reads aligning to the negative strand are offset by -5bp. The adjusted reads will be written into a new bamfile for peak calling
#' or footprinting.


#' @param bamfile of sample to be shifted and predicted
#'
#' @return shifted alignment gal1
#' @import lattice
#' @import ATACseqQC
#' @import Rsamtools
#' @export
#'
#' @examples
#' library(lattice)
#' library(ATACseqQC)
#' library(Rsamtools)
#' \dontrun{
#' bamfile<-system.file("extdata","exp.bam",package = "bulkATACquality")
#' gAl<-GAl(bamfile)
#' }

GAl<-function(bamfile){
  ## bamfile tags to be read in
  possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2",
                                  "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                  "TC", "UQ"),
                      "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                    "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                    "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                    "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                    "U2"))

  bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                       param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0]


  ## files will be output into outPath
  outPath <- paste0("splited_",strsplit(basename(bamfile),"_")[[1]][1])
  dir.create(outPath)
  ## shift the coordinates of 5'ends of alignments in the bamfile
  seqlev <- "chr1" ## subsample data for quick run
  which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
  gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
  shiftedBamfile <- file.path(outPath, "shifted.bam")
  gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
  return(gal1)
}
