#' Generate TSSE score plot and calculate width of curve and cov of score of all points
#'
#' @param tsse TSSE object
#' @param dir output directory
#'
#' @return width of curve at mid-height and cov of score of all points
#' @importFrom stats lm
#' @importFrom stats sd
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics abline
#' @export
#'
#' @examples
#' \dontrun{
#' TSSE_result<-TSSE(dir0,tsse)
#' TSSE_cov<-TSSE_result$cov_TSSscore
#' TSSE_wid<-TSSE_result$wid
#' }

TSSE<-function(dir, tsse){

  # calculate width of two points on the height of mid points

  TSSEvalue<-max(tsse$value[10],tsse$value[11])
  mid<-TSSEvalue-((TSSEvalue-min(tsse$value))/2)
  t<-c(tsse$values[2:20],tsse$values[20])
  diff_tsse<-t-tsse$values
  TSS<-data.frame(position=100*(-9:10-.5),tsse_values=tsse$values,difference=tsse$values-mid,tag=rep(0,20),diff_tsse=diff_tsse)
  TSS_neg<-TSS[1:10,]
  for (i in 1:9){
    TSS_neg$tag[i]<-TSS_neg$difference[i]*TSS_neg$difference[i+1]}

  TSS_pos<-TSS[11:20,]
  for (i in 1:9){
    TSS_pos$tag[i]<-TSS_pos$difference[i]*TSS_pos$difference[i+1]}

  ind_neg<-which(TSS_neg$tag<0 & TSS_neg$diff_tsse>0)
  ind_neg<-ind_neg[length(ind_neg)]
  neg<-TSS_neg[(ind_neg):(ind_neg+1),]
  fit_neg<-lm(tsse_values ~ position, data=neg)

  ind_pos<-which(TSS_pos$tag<0 & TSS_pos$diff_tsse<0)
  ind_pos<-ind_pos[1]
  pos<-TSS_pos[(ind_pos):(ind_pos+1),]
  fit_pos<-lm(tsse_values ~ position, data=pos)

  png(paste0(dir,"_TSSEscore.png"),res=150, width = 2400,height = 1200)
  plot(100*(-9:10-.5), tsse$values, type="b",
       xlab="distance to TSS",
       ylab="aggregate TSS score")
  abline(fit_neg,col="grey")
  abline(fit_pos,col="grey")
  abline(h=mid,col="red",lwd=2)
  dev.off()
  x1<-(mid-summary(fit_neg)$coefficients[1,1])/summary(fit_neg)$coefficients[2,1]
  x2<-(mid-summary(fit_pos)$coefficients[1,1])/summary(fit_pos)$coefficients[2,1]
  wid<-x2-x1
  print(paste0("5./The width in median  is ", wid))

  #calculate coefficient of variation

  cov_TSSscore<-sd(tsse$values)/mean(tsse$values)
  print(paste0("6./The cov of TSSE score is ", cov_TSSscore))

  return(list(wid=wid,cov_TSSscore=cov_TSSscore))
}

