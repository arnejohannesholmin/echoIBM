#*********************************************
#*********************************************
#' Returns omnidirectional directional factor of dimension equal to the first input.
#'
#' @param data  is a list containing the elements 'dire' (elevarion angle in radians [0,pi/2]) and 'wnsz' (product of wave number and size), and optionally 'dira' (azimuth angle in radians [0,pi/2]), 'indi' (index of the fish / transducer element) and 'obln' (oblongness of the ).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname pointSource.TSD
#'
pointSource.TSD<-function(data){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-06-11 - Clean version.
	########### DESCRIPTION: ###########
	# Returns omnidirectional directional factor of dimension equal to the first input.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the elements 'dire' (elevarion angle in radians [0,pi/2]) and 'wnsz' (product of wave number and size), and optionally 'dira' (azimuth angle in radians [0,pi/2]), 'indi' (index of the fish / transducer element) and 'obln' (oblongness of the ).
	
	
	##################################################
	##################################################
	##### Preparation #####
	l=length(data$wnsz)
	
	
	##### Execution and output #####
	double(l)+1
	##################################################
	##################################################
	}
