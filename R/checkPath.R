#' Check whether tools exist
#'
#' @param u name of tools
#' @param path path of tools
#' @param throwError whether throw errors when cannot find tools
#' @return out True: find this tools in your machine; False: cannot find this tools in your machine;
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' .checkPath("samtools", throwError = FALSE)

.checkPath <- function(u = NULL, path = NULL, throwError = TRUE){
  if(is.null(u)){
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE){
    if (Sys.which(x) == "") {
      if(!is.null(path) && file.exists(file.path(path, x))){
        o <- TRUE
      }else{
        if(throwError){
          stop(x, " not found in path, please add ", x, " to path!")
        }else{
          o <- FALSE
        }
      }
    }else{
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all #library(dplyr)
  return(out)
}
