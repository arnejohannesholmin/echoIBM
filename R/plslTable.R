#*********************************************
#*********************************************
#' Gets pulse duration from the max range of the sonar.
#'
#' @param range	The max range of the sonar.
#' @param mode	The mode of the sonar. Use plslTable() to get avaialbe modes.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.setup
#'
plslTable <- function(range=NULL, mode="CWNormal"){
	ranges <- c(    150, 300, 450, 600, 900, 1200, 1500, 2000, 2500, 3000, 3500, 4500, 6000, 8000)
	plsl <- list(
		CWShort   = c(1,   2,   3,   4,   6,    7,    8,   12,   14,   17,   20,   25,   34,   43),
		CWNormal  = c(2,   4,   6,   8,  12,   14,   17,   23,   29,   34,   40,   51,   58,   69),
		CWLong    = c(6,  12,  17,  24,  34,   41,   51,   68,   85,   85,   85,   85,   85,   85),
		FMAuto    = c(6,  12,  18,  24,  36,   48,   60,   80,   85,   85,   85,   85,   85,   85),
		FMShort   = c(2,   4,   6,   8,  12,   16,   20,   27,   33,   40,   47,   60,   80,   85),
		FMNormal  = c(4,   8,  12,  16,  24,   32,   40,   54,   67,   80,   85,   85,   85,   85),
		FMLong    = c(6,  12,  18,  24,  36,   48,   60,   80,   85,   85,   85,   85,   85,   85)
	)
	if(length(range)==0){
		return(list(ranges=ranges, plsl=plsl))
	}
	
	matchRange <- match(range, ranges)
	if(length(matchRange)==0){
		warning(paste0("range not matching any of the following ranges: ", paste(ranges, collapse=", ")))
	}
	if(! mode %in% names(plsl)){
		warning(paste0("mode not matching any of the following modes: ", paste(names(plsl), collapse=", ")))
	}
	plsl[[mode]][matchRange] * 1e-3
}