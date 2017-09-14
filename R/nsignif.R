#*********************************************
#*********************************************
#' Returns an estimate of the number of significant values in 'x', where the significance of the values are measured as the fractional magnitude relative to the maximum value.
#'
#' @param x  is a numeric vector.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname nsignif
#'
nsignif=function(x,na.rm=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-28 - Clean version.
	########### DESCRIPTION: ###########
	# Returns an estimate of the number of significant values in 'x', where the significance of the values are measured as the fractional magnitude relative to the maximum value.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is a numeric vector.
		

	##################################################
	##################################################
	sum(x/max(x,na.rm=na.rm),na.rm=na.rm)
	##################################################
	##################################################
	}
