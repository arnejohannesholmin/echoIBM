#*********************************************
#*********************************************
#' Returns the theoretical far field approximation of the acoustic intensity from a line source of length 'L', as a function of incidence data$direle 'data$dire'.
#'
#' @param data  is a list containing the elements 'dire' (elevarion angle in radians [0,pi/2]) and 'wnsz' (product of wave number and size), and optionally 'dira' (azimuth angle in radians [0,pi/2]), 'indi' (index of the fish / transducer element)-
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname lineSource.TSD
#'
lineSource.TSD<-function(data){
	
	############### LOG: ###############
	# Start: 2010-01-25 - Clean version.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# pi/2 is subtracted from 'data$dire' to specify maximum to be perpendicular to the line source:
	data$dire=sin(data$dire-pi/2)/2 * data$wnsz
	
	
	##### Execution and output #####
	out=(sin(data$dire)/data$dire)^2
	out[data$dire==0]=1
	out
	##################################################
	##################################################
	}
