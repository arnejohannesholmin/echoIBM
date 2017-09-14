#*********************************************
#*********************************************
#' Histogram without plot, faster.
#'
#' @param x  is a vector.
#' @param breaks  is the number of breaks (in which case pretty() is used), or a vector of breaks.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname hist_simple
#'
hist_simple<-function(x,breaks=10){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-10-24 - Clean version.
	########### DESCRIPTION: ###########
	# Histogram without plot, faster.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is a vector.
	# ---breaks--- is the number of breaks (in which case pretty() is used), or a vector of breaks.
		

	##################################################
	##################################################
	if(length(breaks)==1){
		breaks=pretty(c(min(x),max(x)),n=breaks)
		}
	x=findInterval(x,breaks,rightmost.closed=TRUE)
	list(breaks=breaks,counts=tabulate(x))
	##################################################
	##################################################
	}
