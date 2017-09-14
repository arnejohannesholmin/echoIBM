#*********************************************
#*********************************************
#' Implemented ind strings in echoIBM and realted functions
#' 
#' @param ind see ind.expand()
#'
#' @importFrom sonR sonR_implemented
#'
#' @export
#' @rdname sonR_implemented_ind
#'
sonR_implemented_ind <- function(ind){
	if(is.character(ind)){
		if(sonR_implemented(ind, c("sbe", "mbe", "mbs"))){
			list(-(1:100),NULL)
			}
		else if(sonR_implemented(ind, "ofs")){
			list(-(1:300),-(32:34))
			}
		}
	else{
		ind
		}
	}
