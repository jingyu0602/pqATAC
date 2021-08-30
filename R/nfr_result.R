#' Generate PT score plot and fit a line of all points
#' Generate density plot (hexbins) of points to visualize the distribution of points
#' Count the number of points with high coverage and high NFR score
#' Calculate variance and cov of hexbins
#'
#' @param nfr NFR object
#' @param outliers True —— including all points ; False —— excluding outliers
#' @param dir output directory
#'
#' @return p value of fitted line (whether the slope is significantly bigger than 0);
#'         counts of points with high coverage and high NFR score;
#'
#' @import hexbin
#' @import grDevices
#' @export
#'
#' @examples
#' \dontrun{
#' NFR_T<-nfr_result(nfr,T)
#' p_NFR<-NFR_T$p_nfr
#' num_passingpoint<-NFR_T$passpoints_NFR
#' count_hexbins<- NFR_T$num_hexbins
#' var_countinhexbins<-NFR_T$var_hexbins
#' cov_countinhexbins<-NFR_T$cov_hexbins
#'
#' p_NFR_excout<-nfr_result(nfr,F,dir0)
#' }

nfr_result<-function(nfr, outliers,dir){

  if (outliers==TRUE){
    lm_NFR<-lm(nfr$NFR_score ~ nfr$log2meanCoverage)
    p_nfr<-pt(summary(lm_NFR)$coefficients[2,3], lm_NFR$df, lower.tail = F) #已知t值和自由度的话,可以使用pt()计算p value; 如果设为 FALSE 则是计算 P[X > x] 的概率.
    print(paste0("4./The P value of NFR score is ", p_nfr))

    NFR_data<-data.frame(nfr$log2meanCoverage,nfr$NFR_score)
    colnames(NFR_data)<-c("coverage","nfr_score")

    png(paste0(dir,"_NFR_score.png"),res=150, width = 2400,height = 1200)
    plot(nfr$log2meanCoverage, nfr$NFR_score,
         xlab="log2 mean coverage",
         ylab="Nucleosome Free Regions score",
         main="NFRscore for 200bp flanking TSSs"
    )
    text(-5, 2, paste0("P value is ", p_nfr))
    abline(lm_NFR,col="red",lwd=2)
    dev.off()

    png(paste0(dir,"_NFR_score_hexagonalbins.png"))
    bin<-hexbin(nfr$log2meanCoverage,nfr$NFR_score, xbins=15)
    plot(bin, main="NFRscore_Hexagonalbins")
    dev.off()

    bin<-hexbin(nfr$log2meanCoverage,nfr$NFR_score, xbins=15)
    plot(bin, main="NFRscore_Hexagonalbins")
    num_hexbins<-length(bin@count)
    print(paste0("The number of hexbins is ", num_hexbins))
    var_hexbins<-var(bin@count)
    print(paste0("The variance of hexbins is ", var_hexbins))
    cov_hexbins<-var(bin@count)/mean((bin@count))
    print(paste0("The cov of hexbins is ", cov_hexbins))

    NFR_pass<-NFR_data[NFR_data$coverage>-5 & NFR_data$nfr_score>0,]
    passpoints_NFR<-dim(NFR_pass)[1]
    print(paste0("The number of passing points is ", passpoints_NFR))

    return(list(p_nfr=p_nfr, passpoints_NFR=passpoints_NFR, num_hexbins=num_hexbins, var_hexbins=var_hexbins,cov_hexbins=cov_hexbins))

  }

  if(outliers==FALSE){
    out_NFRscore<-boxplot.stats(nfr$log2meanCoverage)$out
    NFR_score_rmout<-nfr$NFR_score[which(! nfr$log2meanCoverage %in% out_NFRscore)]
    NFR_log2meanCoverage_rmout<-nfr$log2meanCoverage[which(! nfr$log2meanCoverage %in% out_NFRscore)]

    NFR_data<-data.frame(NFR_log2meanCoverage_rmout,NFR_score_rmout)
    colnames(NFR_data)<-c("coverage","nfr_score")
    lm_NFR<-lm(nfr_score ~ coverage,data=NFR_data)
    p_nfr<-pt(summary(lm_NFR)$coefficients[2,3], lm_NFR$df, lower.tail = F) #已知t值和自由度的话,可以使用pt()计算p value; 如果设为 FALSE 则是计算 P[X > x] 的概率.
    print(paste0("3./The P value of NFR score excluding outliers is ", p_nfr))

    png(paste0(dir,"_NFR_score_excout.png"),res=150, width = 2400,height = 1200)
    plot(NFR_log2meanCoverage_rmout, NFR_score_rmout,
         xlab="log2 mean coverage",
         ylab="Nucleosome Free Regions score",
         main="NFRscore for 200bp flanking TSSs")
    text(-5, 2, paste0("P value is ", p_nfr))
    abline(lm_NFR,col="red",lwd=2)
    dev.off()

    png(paste0(dir,"_NFR_score_excout_hexagonalbins.png"))
    bin<-hexbin(NFR_log2meanCoverage_rmout,NFR_score_rmout, xbins=15)
    plot(bin, main="NFRscore_excout_Hexagonalbins")
    dev.off()

    return(p_nfr)
  }

}
