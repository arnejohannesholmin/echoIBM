#*********************************************
#*********************************************
#' Discards fish outside of the sampling region of the sonar/echosounder.
#'
#' @param dynschool  is a list of the dynamic fish information to be subsetted.
#' @param data  is a list holding the variables 'esnm', and vessel spesifications.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom sonR is.sonar rotate3D
#' @importFrom TSD car2sph
#'
#' @export
#' @rdname echoIBM.fishInside
#'
echoIBM.extraxtRandSel <- function(dynschool, data, rand.sel, dynschoolnames, staticschoolnames){
	# Extract a random selection of the targets using 'rand.sel':
	Nl <- max(length(dynschool$psxf), length(dynschool$psyf), length(dynschool$pszf))
	if(0<rand.sel[1] && rand.sel[1]<1){
		if(length(rand.sel)>1){
			set.seed(rand.sel[2])
		}
		affected.variables <- c(dynschoolnames, staticschoolnames)
		#selection <- sample(c(TRUE, FALSE), Nl, TRUE, c(rand.sel[1], 1-rand.sel[1]))
		selection <- sample(seq_len(Nl), round(Nl * rand.sel[1]))
		for(j in seq_along(affected.variables)){
			thisvar <- affected.variables[j]
			if(length(data[[thisvar]])==Nl && !is.function(data[[thisvar]])){
				data[[thisvar]] <- data[[thisvar]][selection]
			}
			if(length(dynschool[[thisvar]])==Nl && !is.function(dynschool[[thisvar]])){
				dynschool[[thisvar]] <- dynschool[[thisvar]][selection]
			}
		}
	}
	
	list(dynschool=dynschool, data=data)
}