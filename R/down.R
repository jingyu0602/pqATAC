#' downsampling bam files and generate metrics when downsampling 6000000, 9000000, 12000000,
#' 15000000, 16000000, 17000000, 18000000, 19000000, 20000000 reads, see how metrics change
#'
#' @param bam bamfile to be downsampled
#' @param outdir directory of output files
#' @param path_to_samtools the intalled path of samtools e.g. "/tools/bin/samtools"
#' @param path_to_sambamba the intalled path of sambamba e.g. "/tools/bin/sambamba"
#' @param path_to_macs2 the intalled path of macs2 e.g. "/tools/bin/macs2"
#' @param path_to_bedtools the intalled path of bedtools e.g. "/tools/bin/bedtools"
#'
#' @return null
#'
#' @import ATACseqQC
#' @import Rsamtools
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import phastCons100way.UCSC.hg19
#' @import ChIPpeakAnno
#' @import ggplot2
#' @import utils
#' @export
#'
#' @examples
#' \dontrun{
#' down(bam,"/project/exa_data/")
#' }

down<-function(bam,outdir,path_to_samtools=NULL,path_to_sambamba=NULL,path_to_macs2=NULL,path_to_bedtools=NULL){

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


  # load required BED files
  promoters<-system.file("bed","promoterhg19_2kb.bed", package = "bulkATACquality")
  black_M<-system.file("bed","M_size.bed", package = "bulkATACquality")
  blacklist<-system.file("bed","ENCFF000KJP.bed", package = "bulkATACquality")
  chromSizes<-system.file("bed","hg19.chromSizes", package = "bulkATACquality")
  chrBed<-system.file("bed","chromS.bed", package = "bulkATACquality")

RL=basename(substr(bam,start=1,stop=6))

reads<-read.csv(paste0(RL, "/5.Results/", RL, "_reads.csv"))
totalReads<-reads$num_reads


num_downsampling<-c(6000000,9000000,12000000,15000000,16000000,17000000,18000000,19000000,20000000)
proportion<-num_downsampling/totalReads

dir6<-paste0(outdir, "/6.Downsampling/")
dir.create(dir6)

p_PT<-c()
p_PT_excout<-c()
p_NFR<-c()
p_NFR_excout<-c()
num_passingpoint<-c()
count_hexbins<-c()
var_countinhexbins<-c()
cov_countinhexbins<-c()
TSSE_cov<-c()
TSSE_wid<-c()
TSSsignal_cov<-c()
TSSsignal_difference<-c()
rn<-c()

for (i in 1:length(proportion)){
  m=paste0(num_downsampling[i]/1000000,"m")
  print(m)
  print(paste0("*** start downsampling ", RL))
  print(paste0("The percentage of downsampling is ",proportion[i]))

  if (proportion[i] >=0.1){
    system(paste0("samtools view -bs 42.",round(10e6*proportion[i]), " ", bam, "  > ", dir6, RL, "_", m,".bam"))
  } else if (proportion[i] >= 0.01 & proportion[i] < 0.1){
    system(paste0("samtools view -bs 42.0",round(10e6*proportion[i]), " ", bam, "  > ", dir6, RL, "_", m,".bam"))
  } else if (proportion[i] < 0.01) {
    system(paste0("samtools view -bs 42.00",round(10e6*proportion[i]), " ", bam, "  > ", dir6, RL, "_", m,".bam"))
  }
  print(paste0("*** finish downsampling ", RL))


  ### 2. Filtering and peak calling

  sub_bam<-paste0(dir6, RL, "_", m,".bam")
  print(paste0("Start filting BAM file ", basename(sub_bam)))

  dir7<-paste0(dir6, RL,"_",m,"/")
  system(paste0("mkdir ", dir7))
#  system(paste0("./downsampling.sh ",sub_bam," ",num_downsampling[i]," ",dir7))

  system('echo -e "\t1. Flag Stats"')
  system(paste0(samtools," flagstat ", sub_bam, "  > ", dir7, "/", RL, "_flagStats"))

  system('echo -e "\t2. Count MT reads"')
  system(paste0(samtools, " index ", sub_bam))
  system(paste0(samtools, " view -c " , sub_bam," chrM > ", dir7, "/", RL, "_mitoReads"))

  system('echo -e "\t3. Data is PE ---> $pairs\n\tFilter alignment, proper Pairs and stats"')
  system(paste0(sambamba," view -t 4 -f bam -F ", noquote("\"proper_pair\" \ "), sub_bam," -o ", dir7, "/", RL, "_pairs.bam"))

  system('echo -e "\t4. Remove PCR dups and stats"')
  system(paste0(sambamba, " markdup -p -t 4 -r --hash-table-size=1000000  ",dir7, "/", RL, "_pairs.bam ", dir7, "/", RL, "_pairs_dedup.bam"))

  system('echo -e "\t5. Intercept BlackListed regions"')
  system(paste0(samtools, " view -q 30 -L ", blacklist, " -U ", dir7, "/", RL, "_pairs_dedup_filt.bam -b ", dir7, "/", RL, "_pairs_dedup.bam > ", dir7, "/", RL, "_pairs_dedup_crap"))

  system(paste0(samtools, " view -L ", black_M, " -U ", dir7, "/", RL, "_pairs_dedup_filt_noMT.bam -b ", dir7, "/", RL, "_pairs_dedup_filt.bam > ", dir7, "/", RL, "_pairs_dedup_2crap"))
  system(paste0(samtools, " flagstat ", dir7, "/", RL, "_pairs_dedup_filt_noMT.bam > ", dir7, "/", RL, "_pairs_dedup_filt_noMT_flagStats"))
  system(paste0(samtools, " index ", dir7, "/", RL, "_pairs_dedup_filt_noMT.bam"))

  system('echo -e "\t7. Remove Intermediate Files, some bam also"')
  system(paste0("rm ", dir7, "/", RL, "_*crap"))
  system(paste0("rm ", dir7, "/", RL, "_pairs.bam*"))
  system(paste0("rm ", dir7, "/", RL, "_pairs_dedup.bam*"))
  system(paste0("rm ", dir7, "/", RL, "_pairs_dedup_filt.bam"))

  fBam=paste0(dir7, "/", RL, "_pairs_dedup_filt_noMT.bam")

  system('echo -e "\t8. Counting reads in promoters"')
  system(paste0(samtools, " view -c -L ", promoters, " -b ", fBam ," > ", dir7 , "/", RL, "_ReadsInProm"))

  system('echo -e "\t9. MACS2:::params for Tn5 insert sites, params ENCODE pipe for NFRs\n\tCalcule Reads in NFR peaks\n\tCalcule Reads in NFR peaks overlapping promotors"')
  system(paste0(macs2, " callpeak --nomodel --extsize 150 --shift -75 -t ", fBam, " -f BAM -n ", dir7, "/", RL,"\"_NFR\" --keep-dup all --gsize hs"))
  system(paste0("cut -f1-3 ", dir7, "/", RL, "_NFR_peaks.narrowPeak | ", samtools, " view -L - -c -b ", fBam, " > ", dir7, "/", RL, "_ReadsInNFRpeaks"))
  system(paste0("cut -f1-3 ", dir7, "/", RL, "_NFR_peaks.narrowPeak | ", intersectBed, " -wa -a - -b ", promoters, " | ", samtools, " view -L - -c -b ", fBam, " > ", dir7, "/", RL, "_ReadsInNFRpeaks_overProm"))

  system('echo -e "\t10. MACS2:::params for PE fragments TLEN"')
  system(paste0(macs2, " callpeak -t ", fBam, " -f BAMPE -n ", dir7, "/", RL, "_BAMPE  --keep-dup all --gsize hs"))

  system(paste0("cut -f1-3 ", dir7, "/", RL, "_BAMPE_peaks.narrowPeak | ", samtools,  " view -L - -c -b ", fBam, " > ", dir7, "/", RL, "_ReadsInBAMPEpeaks"))
  system(paste0("cut -f1-3 ", dir7, "/", RL, "_BAMPE_peaks.narrowPeak | ", intersectBed, " -wa -a - -b ", promoters, " | ", samtools, " view -L - -c -b ", fBam, " > ", dir7, "/", RL, "_ReadsInBAMPEfrag_overProm"))

  system('echo -e "\t11. Reads in Background"')
  system(paste0("cut -f1-3 ", dir7, "/", RL, "_BAMPE_peaks.narrowPeak ", dir7, "/", RL, "_NFR_peaks.narrowPeak | sort -k1,1 -k2,2n  | ", mergeBed, " -i - | grep -v chrM | ", bedtools, " subtract -wb -a ", chrBed, " -b - | wc -l > ", dir7, "/", RL, "_PeaksInBackground"))
  system(paste0("cut -f1-3 ", dir7, "/", RL, "_BAMPE_peaks.narrowPeak ", dir7, "/", RL, "_NFR_peaks.narrowPeak | sort -k1,1 -k2,2n  | ", mergeBed, " -i - | grep -v chrM | ", bedtools, " subtract -a ",chrBed, " -b - | ", samtools, " view -L - -c -b ", fBam, " > ", dir7, "/", RL, "_ReadsInBackground"))

  system(paste0("echo \"num_reads\", \"mito_reads\", \"readsAfterFilter\", \"readsInPro\", \"readsInNFR\", \"readsInNFRoverPro\", \"readsInBAMPE\", \"readsInBAMPEoverPro\", \"readsInBack\" >  ", dir7, "/", RL, "\"_colnames.csv\""))
  system(paste0("echo $(grep QC ",dir7, "/", RL, "_flagStats | sed 's/+.*//g' ), $(cat ", dir7, "/", RL, "_mitoReads), $(grep QC ", dir7, "/", RL, "_pairs_dedup_filt_noMT_flagStats | sed 's/+.*//g' ), $(cat ", dir7, "/", RL, "_ReadsInProm), $(cat ", dir7, "/", RL, "_ReadsInNFRpeaks), $(cat ", dir7, "/", RL, "_ReadsInNFRpeaks_overProm), $(cat ", dir7, "/", RL,"_ReadsInBAMPEpeaks), $(cat ", dir7, "/", RL, "_ReadsInBAMPEfrag_overProm), $(cat ", dir7, "/", RL, "_ReadsInBackground) > ", dir7, "/", RL, "_numreads.csv"))

  system(paste0("cat ", dir7, RL,"_colnames.csv ", dir7 ,"/", RL, "_numreads.csv > ",  dir7, "/", RL, "_reads.csv"))

  ### 3. Genreting Metrics Plots

  bamfile<-paste0(dir7, RL, "_pairs_dedup_filt_noMT.bam")
  print(bamfile)

  print(paste0("Start running metrics of sample ",RL))


  ### call functions above

  gAl<-GAl(bamfile)

  dir8<-paste0(dir7,"Plots/")
  dir.create(dir8)
  dir0<-paste0(dir8,RL)

  txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

  ### PTscore plot
  pt <- PTscore(gAl, txs)
  p_PT0<-pt_result(pt, outliers=TRUE, dir0)
  p_PT_excout0<-pt_result(pt, outliers=FALSE, dir0)

  ### NFR score plot
  nfr <- NFRscore(gAl, txs)

  NFR_T<-nfr_result(nfr,T,dir0)
  p_NFR0<-NFR_T$p_nfr
  num_passingpoint0<-NFR_T$passpoints_NFR
  count_hexbins0<- NFR_T$num_hexbins
  var_countinhexbins0<-NFR_T$var_hexbins
  cov_countinhexbins0<-NFR_T$cov_hexbins

  p_NFR_excout0<-nfr_result(nfr,F,dir0)

  ### TSS signal plot

  TSSsignal_result<-TSS_signal(bamfile,gAl,dir0)

  TSSsignal_cov0<-TSSsignal_result$cov_tsssignal
  TSSsignal_difference0<-TSSsignal_result$TSS_diff

  ### TSSEscore plot
  tsse <- TSSEscore(gAl, txs)
  TSSE_result<-TSSE(dir0, tsse)

  TSSE_cov0<-TSSE_result$cov_TSSscore
  TSSE_wid0<-TSSE_result$wid


  p_PT[i]<-p_PT0
  p_PT_excout[i]<-p_PT_excout0
  p_NFR[i]<-p_NFR0
  p_NFR_excout[i]<-p_NFR_excout0
  num_passingpoint[i]<-num_passingpoint0
  count_hexbins[i]<-count_hexbins0
  var_countinhexbins[i]<-var_countinhexbins0
  cov_countinhexbins[i]<-cov_countinhexbins0
  TSSE_cov[i]<-TSSE_cov0
  TSSE_wid[i]<-TSSE_wid0
  TSSsignal_cov[i]<-TSSsignal_cov0
  TSSsignal_difference[i]<-TSSsignal_difference0
  rn[i]<-m
  #result_table<-data.frame(p_PT,p_PT_excout,p_NFR,p_NFR_excout,num_passingpoint,count_hexbins,var_countinhexbins,cov_countinhexbins,TSSE_cov,TSSE_wid,TSSsignal_cov,TSSsignal_difference)
  #rownames(result_table)<-substr(basename(bamfile),start=1,stop=9)
  #write.csv(result_table,file = paste0(dir0,"_result_table.csv"))

  system(paste0("rm ",dir6,"*bam*"))
  system(paste0("rm -r ","splited_",RL))
  system(paste0("rm ",dir7,"*bam*"))
  system(paste0("echo \"num_reads\", \"mito_reads\", \"readsAfterFilter\", \"readsInPro\", \"readsInNFR\", \"readsInNFRoverPro\", \"readsInBAMPE\", \"readsInBAMPEoverPro\", \"readsInBack\" >  ", dir6, "/", RL, "\"_colnames.csv\""))
  system(paste0("cat ", dir7, RL, "_numreads.csv >> ", dir6, "Reads.csv"))
}

system(paste0("cat ",dir6,"/", RL, "_colnames.csv ", dir6, "Reads.csv > ", dir6, "Reads_downsampling.csv"))
system(paste0("rm ",dir6,"Reads.csv", " ",dir6, "/", RL, "_colnames.csv "))

result_table<-data.frame(p_PT,p_PT_excout,p_NFR,p_NFR_excout,num_passingpoint,count_hexbins,var_countinhexbins,cov_countinhexbins,TSSE_cov,TSSE_wid,TSSsignal_cov,TSSsignal_difference)
rownames(result_table)<-rn
write.csv(result_table,paste0(dir6,"Metrics_down.csv"))
res_table<-data.frame(result_table, "num_downsampling"=1:9)

reads_table<-read.csv(paste0(dir6,"Reads_downsampling.csv"),header=T)
reads_table$mito_reads<-reads_table$mito_reads/reads_table$num_reads
reads_table[,4:9]<-reads_table[,4:9]/reads_table$readsAfterFilter
reads_table<-data.frame(reads_table, "num_downsampling"=1:9)

### Downsampling plots:
p1<-ggplot(data=res_table,aes(x=num_downsampling,y=num_passingpoint))+
  geom_line() +   geom_point() +
  scale_x_continuous(name="Number of downsampling",labels=c("6000000","9000000","12000000","15000000","16000000","17000000","18000000","19000000","20000000"),limits=c(0,10), breaks=seq(1,9,1))+theme_bw() +
  scale_y_continuous(name="Number of passing points") +
  ggtitle("Number of passing points for downsampling") +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5,  face="bold", size=14))
ggsave(paste0(dir6,"Number_passingpoints_downsampling.png"),p1)

p2<-ggplot(data=reads_table,aes(x=num_downsampling,y=readsInNFR))+
  geom_line() +   geom_point() +
  scale_x_continuous(name="Number of downsampling",labels=c("6000000","9000000","12000000","15000000","16000000","17000000","18000000","19000000","20000000"),limits=c(0,10), breaks=seq(1,9,1))+theme_bw() +
  scale_y_continuous(name="Percentage of reads in NFR peaks") +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5)) +
  ggtitle("Percentage of reads in NFR peaks for downsampling") +
  theme(plot.title = element_text(hjust = 0.5,  face="bold", size=14))
ggsave(paste0(dir6,"Perc_readsInNFR_downsampling.png"),p2)

p3<-ggplot(data=reads_table,aes(x=num_downsampling,y=readsInNFRoverPro))+
  geom_line() +   geom_point() +
  scale_x_continuous(name="Number of downsampling",labels=c("6000000","9000000","12000000","15000000","16000000","17000000","18000000","19000000","20000000"),limits=c(0,10), breaks=seq(1,9,1))+theme_bw() +
  scale_y_continuous(name="Percentage of reads in NFR peaks overlapping promoters") +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5)) +
  ggtitle("Percentage of reads in NFR peaks overlapping promoters for downsampling") +
  theme(plot.title = element_text(hjust = 0.5,  face="bold", size=14))
ggsave(paste0(dir6,"Perc_readsInNFRoverPro_downsampling.png"),p3)

p4<-ggplot(data=reads_table,aes(x=num_downsampling,y=readsInBAMPE))+
  geom_line() +   geom_point() +
  scale_x_continuous(name="Number of downsampling",labels=c("6000000","9000000","12000000","15000000","16000000","17000000","18000000","19000000","20000000"),limits=c(0,10), breaks=seq(1,9,1))+theme_bw() +
  scale_y_continuous(name="Percentage of reads in BAMPE peaks") +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5)) +
  ggtitle("Percentage of reads in BAMPE peaks for downsampling") +
  theme(plot.title = element_text(hjust = 0.5,  face="bold", size=14))
ggsave(paste0(dir6,"Perc_readsInBAMPE_downsampling.png"),p4)

p5<-ggplot(data=reads_table,aes(x=num_downsampling,y=readsInBAMPEoverPro))+
  geom_line() +   geom_point() +
  scale_x_continuous(name="Number of downsampling",labels=c("6000000","9000000","12000000","15000000","16000000","17000000","18000000","19000000","20000000"),limits=c(0,10), breaks=seq(1,9,1))+theme_bw() +
  scale_y_continuous(name="Percentage of reads in BAMPE peaks overlapping promoters") +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5)) +
  ggtitle("Percentage of reads in BAMPE peaks overlapping promoters for downsampling") +
  theme(plot.title = element_text(hjust = 0.5,  face="bold", size=14))
ggsave(paste0(dir6,"Perc_readsInBAMPEoverPro_downsampling.png"),p5)}
