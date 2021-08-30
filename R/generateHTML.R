#' Generate HTML of sample quality prediction
#'
#' @param dir folder name of all results (should be same as outdir of previous function)
#' @return null
#' @import rmarkdown
#'
#'
#' @examples
#' \dontrun{
#' generateHTML("/project/RL0001/")
#' }
#' @export


generateHTML<-function(dir){
rmd <- system.file("rmd", "bulkATACqualityPrediction.Rmd", package="bulkATACquality", mustWork=TRUE)
rmarkdown::render(rmd, params = list (folder= dir), output_file = paste0(dir, "/",basename(dir), "_HTML"))
}
