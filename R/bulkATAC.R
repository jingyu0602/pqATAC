#' Filtering BAM files and call peaks; Generate metrics plot and statistics used to predict sample quality
#'
#' @param original_bamfile The BAM file of sample to be predicted
#' @param outdir output directory of all results
#' @param path_to_samtools the intalled path of samtools e.g. "/tools/bin/samtools"
#' @param path_to_sambamba the intalled path of sambamba e.g. "/tools/bin/sambamba"
#' @param path_to_macs2 the intalled path of macs2 e.g. "/tools/bin/macs2"
#' @param path_to_bedtools the intalled path of bedtools e.g. "/tools/bin/bedtools"
#'
#' @return null
#' @import utils
#' @import ATACseqQC
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#'
#' @examples
#' \dontrun{
#' library(bulkATACquality)
#' library(ATACseqQC)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(phastCons100way.UCSC.hg19)
#' library(ChIPpeakAnno)
#' library(dplyr)
#' bamfile<-system.file("extdata","ex0001.bam",package = "bulkATACquality")
#' outdir<-substr(basename(bamfile),1,6)
#' dir.create(outdir)
#' bulkATAC(bamfile,outdir)
#' }
#' @export

### QC, Filtering and peak calling
bulkATAC<-function(original_bamfile,outdir,path_to_samtools=NULL,path_to_sambamba=NULL,path_to_macs2=NULL,path_to_bedtools=NULL){
  # check the path of softwares
  if(!is.null(path_to_samtools)){
    samtools<-path_to_samtools
    print(samtools)
  } else{
    if(.checkPath("samtools", throwError = FALSE) == TRUE){
      samtools<-"samtools"
      print(samtools)
    } else {
      write("ERROR:Please install samtools", stderr())
    }
  }

  if(!is.null(path_to_sambamba)){
    sambamba<-path_to_sambamba
    print(sambamba)
  } else{
    if(.checkPath("sambamba", throwError = FALSE) == TRUE){
      sambamba<-"sambamba"
      print(sambamba)
    } else {
      write("ERROR:Please install sambamba", stderr())
    }
  }

  if(!is.null(path_to_bedtools)){
    bedtools<-path_to_bedtools
    intersectBed<-paste0(gsub("bedtools","",path_to_bedtools),"intersectBed")
    mergeBed<-paste0(gsub("bedtools","",path_to_bedtools),"mergeBed")
    print(bedtools)
    print(intersectBed)
    print(mergeBed)
  } else{
    if(.checkPath("bedtools", throwError = FALSE) == TRUE){
      bedtools<-"bedtools"
      intersectBed<-"intersectBed"
      mergeBed<-"mergeBed"
      print(bedtools)
      print(intersectBed)
      print(mergeBed)
      } else {
      write("ERROR:Please install bedtools", stderr())
    }
  }

  if(!is.null(path_to_macs2)){
    macs2<-path_to_macs2
    print(macs2)
  } else{
    if(.checkPath("macs2", throwError = FALSE) == TRUE){
      macs2<-"macs2"
      print(macs2)
    } else {
      write("ERROR:Please install macs2", stderr())
    }
  }


  RL<-substr(basename(original_bamfile),start=1,stop=6)

  # load required BED files
  promoters<-system.file("bed","promoterhg19_2kb.bed", package = "bulkATACquality")
  black_M<-system.file("bed","M_size.bed", package = "bulkATACquality")
  blacklist<-system.file("bed","ENCFF000KJP.bed", package = "bulkATACquality")
  chromSizes<-system.file("bed","hg19.chromSizes", package = "bulkATACquality")
  chrBed<-system.file("bed","chromS.bed", package = "bulkATACquality")

  dir4=paste0(outdir, "/4.Alignment")
  dir5=paste0(outdir, "/5.Results")

  dir.create(dir4)
  dir.create(dir5)
  system(paste0("cp ",original_bamfile, " ", dir4))

  # 1. STATS
  system('echo -e "\t1 Flag Stats"')
  system(paste0(samtools," flagstat ", dir4, "/", RL, ".bam  > ", dir4, "/", RL, "_flagStats"))

  # 2. Count MT reads
  system('echo -e "\t2. Count MT reads"')
  system(paste0(samtools," index ", dir4, "/", RL, ".bam"))
  system(paste0(samtools, " view -c " , dir4, "/", RL, ".bam chrM > ", dir4, "/", RL, "_mitoReads"))

  # 3. Proper Pairs
  system('echo -e "\t3. Data is PE ---> $pairs\n\tFilter alignment, proper Pairs and stats"')
  system(paste0(sambamba," view -t 4 -f bam -F ", noquote("\"proper_pair\" \ "), dir4, "/", RL, ".bam -o ", dir4, "/", RL, "_pairs.bam"))

  # 4. PCR dups
  system('echo -e "\t4. Remove PCR dups and stats"')
  system(paste0(sambamba, " markdup -p -t 4 -r --hash-table-size=1000000  ",dir4, "/", RL, "_pairs.bam ", dir4, "/", RL, "_pairs_dedup.bam"))

  # 5. Remove blacklist
  system('echo -e "\t5. Intercept BlackListed regions"')
  system(paste0(samtools, " view -q 30 -L ", blacklist, " -U ", dir4, "/", RL, "_pairs_dedup_filt.bam -b ", dir4, "/", RL, "_pairs_dedup.bam > ", dir4, "/", RL, "_pairs_dedup_crap"))

  # 6. Remove MT reads
  system(paste0(samtools, " view -L ", black_M, " -U ", dir4, "/", RL, "_pairs_dedup_filt_noMT.bam -b ", dir4, "/", RL, "_pairs_dedup_filt.bam > ", dir4, "/", RL, "_pairs_dedup_2crap"))
  system(paste0(samtools, " flagstat ", dir4, "/", RL, "_pairs_dedup_filt_noMT.bam > ", dir4, "/", RL, "_pairs_dedup_filt_noMT_flagStats"))
  system(paste0(samtools, " index ", dir4, "/", RL, "_pairs_dedup_filt_noMT.bam"))

  # 7. RMs
  system('echo -e "\t7. Remove Intermediate Files, some bam also"')
  system(paste0("rm ", dir4, "/", RL, "_*crap"))
  system(paste0("rm ", dir4, "/", RL, "_pairs.bam*"))
  system(paste0("rm ", dir4, "/", RL, "_pairs_dedup.bam*"))
  system(paste0("rm ", dir4, "/", RL, "_pairs_dedup_filt.bam"))

  # Just Keep original alignment + final filtered bam + index final + flagStats all + coverage.bw final
  fBam=paste0(dir4, "/", RL, "_pairs_dedup_filt_noMT.bam")

  # 8. Reads in Promotors
  system('echo -e "\t8. Counting reads in promoters"')
  system(paste0(samtools, " view -c -L ", promoters, " -b ", fBam ," > ", dir4 , "/", RL, "_ReadsInProm"))

  # 9. MACS2 NFR peaks
  system('echo -e "\t9. MACS2:::params for Tn5 insert sites, params ENCODE pipe for NFRs\n\tCalcule Reads in NFR peaks\n\tCalcule Reads in NFR peaks overlapping promotors"')
  system(paste0(macs2, " callpeak --nomodel --extsize 150 --shift -75 -t ", fBam, " -f BAM -n ", dir4, "/", RL,"\"_NFR\" --keep-dup all --gsize hs"))
  system(paste0("cut -f1-3 ", dir4, "/", RL, "_NFR_peaks.narrowPeak | ", samtools, " view -L - -c -b ", fBam, " > ", dir4, "/", RL, "_ReadsInNFRpeaks"))
  system(paste0("cut -f1-3 ", dir4, "/", RL, "_NFR_peaks.narrowPeak | ", intersectBed, " -wa -a - -b ", promoters, " | ", samtools, " view -L - -c -b ", fBam, " > ", dir4, "/", RL, "_ReadsInNFRpeaks_overProm"))

  # 10. MACS BAMPE peaks
  system('echo -e "\t10. MACS2:::params for PE fragments TLEN"')
  system(paste0(macs2, " callpeak -t ", fBam, " -f BAMPE -n ", dir4, "/", RL, "_BAMPE  --keep-dup all --gsize hs"))

  system(paste0("cut -f1-3 ", dir4, "/", RL, "_BAMPE_peaks.narrowPeak | ", samtools,  " view -L - -c -b ", fBam, " > ", dir4, "/", RL, "_ReadsInBAMPEpeaks"))
  system(paste0("cut -f1-3 ", dir4, "/", RL, "_BAMPE_peaks.narrowPeak | ", intersectBed, " -wa -a - -b ", promoters, " | ", samtools, " view -L - -c -b ", fBam, " > ", dir4, "/", RL, "_ReadsInBAMPEfrag_overProm"))

  # 11. Reads and Peaks In Background
  system('echo -e "\t11. Reads in Background"')
  system(paste0("cut -f1-3 ", dir4, "/", RL, "_BAMPE_peaks.narrowPeak ", dir4, "/", RL, "_NFR_peaks.narrowPeak | sort -k1,1 -k2,2n  | ", mergeBed, " -i - | grep -v chrM | ", bedtools, " subtract -wb -a ", chrBed, " -b - | wc -l > ", dir4, "/", RL, "_PeaksInBackground"))
  system(paste0("cut -f1-3 ", dir4, "/", RL, "_BAMPE_peaks.narrowPeak ", dir4, "/", RL, "_NFR_peaks.narrowPeak | sort -k1,1 -k2,2n  | ", mergeBed, " -i - | grep -v chrM | ", bedtools, " subtract -a ",chrBed, " -b - | ", samtools, " view -L - -c -b ", fBam, " > ", dir4, "/", RL, "_ReadsInBackground"))

  system(paste0("echo \"num_reads\", \"mito_reads\", \"readsAfterFilter\", \"readsInPro\", \"readsInNFR\", \"readsInNFRoverPro\", \"readsInBAMPE\", \"readsInBAMPEoverPro\", \"readsInBack\" >  ", dir4, "/", RL, "\"_colnames.csv\""))
  system(paste0("echo $(grep QC ",dir4, "/", RL, "_flagStats | sed 's/+.*//g' ), $(cat ", dir4, "/", RL, "_mitoReads), $(grep QC ", dir4, "/", RL, "_pairs_dedup_filt_noMT_flagStats | sed 's/+.*//g' ), $(cat ", dir4, "/", RL, "\"_ReadsInProm\"), $(cat ", dir4, "/", RL, "_ReadsInNFRpeaks), $(cat ", dir4, "/", RL, "_ReadsInNFRpeaks_overProm), $(cat ", dir4, "/", RL,"_ReadsInBAMPEpeaks), $(cat ", dir4, "/", RL, "_ReadsInBAMPEfrag_overProm), $(cat ", dir4, "/", RL, "_ReadsInBackground) > ", dir4, "/", RL, "_numreads.csv"))
  system(paste0("cat ", dir4, "/", RL, "_colnames.csv ", dir4, "/", RL, "_numreads.csv > ", dir5, "/", RL, "_reads.csv"))

  bamfile<-paste0(outdir,"/4.Alignment/",RL,"_pairs_dedup_filt_noMT.bam")
  gAl<-GAl(bamfile)
  dir0<-paste0(outdir,"/5.Results/",RL)
  txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

  ### PTscore plot
  pt <- PTscore(gAl, txs)

  p_PT<-pt_result(pt, outliers=TRUE, dir0)
  p_PT_excout<-pt_result(pt, outliers=FALSE, dir0)


  ### NFR score plot
  nfr <- NFRscore(gAl, txs)

  NFR_T<-nfr_result(nfr,T,dir0)
  p_NFR<-NFR_T$p_nfr
  num_passingpoint<-NFR_T$passpoints_NFR
  count_hexbins<- NFR_T$num_hexbins
  var_countinhexbins<-NFR_T$var_hexbins
  cov_countinhexbins<-NFR_T$cov_hexbins

  p_NFR_excout<-nfr_result(nfr,F,dir0)


  ### TSSEscore plot
  tsse <- TSSEscore(gAl, txs)
  TSSE_result<-TSSE(dir0, tsse)
  TSSE_cov<-TSSE_result$cov_TSSscore
  TSSE_wid<-TSSE_result$wid


  ### TSS signal plot

  TSSsignal_result<-TSS_signal(bamfile,gAl, dir0)
  TSSsignal_cov<-TSSsignal_result$cov_tsssignal
  TSSsignal_difference<-TSSsignal_result$TSS_diff

  result_table<-data.frame(p_PT,p_PT_excout,p_NFR,p_NFR_excout,num_passingpoint,count_hexbins,var_countinhexbins,cov_countinhexbins,TSSE_cov,TSSE_wid,TSSsignal_cov,TSSsignal_difference)
  rownames(result_table)<-substr(basename(bamfile),start=1,stop=10)
  write.csv(result_table,file = paste0(dir0,"_metrics.csv"))
  system("rm -r splited*")
}
