#*********************************************
#*********************************************
#' Gets the beam pattern values from the transducer, either at emision or reception, depending on the value of direction.
#'
#' @param data  	is the list of data for the simulation, containing all of the following elements:
#' \itemize{
#'   \item fishdirT
#'   \item fishposV
#'   \item dira
#'   \item dire
#'   \item wavenumber
#'   \item lthesel
#'   \item sllf
#'   \item indi
#'   \item bpt*
#'   \item rad*
#'   \item pbp*
#' }
#' @param direction	Integer: 1 is emision and 2 is reception.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.oneping.j2
#'
echoIBM_setB_T <- function(data, direction=1, beamnr=1){

	# Set the names of the beam pattern function, the radial dimension and the parametric beam pattern:
	bpt_name <- paste0("bpt", direction)
	rad_name <- paste0("rad", direction)
	pbp_name <- paste0("pbp", direction)
	
	# Fish positions in (T) for the current beam, giving the fish directions 'data$fishdirT' in the current beam:
	data$fishdirT <- rotate3D(
		data$fishposV, 
		"zx", 
		cbind(data$dira[beamnr] - pi/2, -data$dire[beamnr]), 
		drop.out=FALSE, 
		sph.out=TRUE
	)

	# Define the beam pattern input data. For circularPiston_ellipticRadius() sllf is not used, and for circularPiston() neither ssil nor dire are used:
	beamPatternInput <- list(
		# The azimuth and elevation angles [number of fish * number of beams]:
		dira = c(data$fishdirT[,2,]), 
		dire = c(data$fishdirT[,3,]), 
		# The product of wave number and size, and Side Lobe Level Factor [number of fish * number of beams, 2]:
		#wnsz = data$uniquek[i] * data$[[rad_name]][beamnr,,drop=FALSE], 
		wnsz = repm(data$wavenumber[beamnr] * data[[rad_name]][beamnr,,drop=FALSE], data$lthesel, byrow=TRUE), 
		sllf = repm(data$sllf[beamnr,,drop=FALSE], data$lthesel, byrow=TRUE), 
		indi = rep(data$indi[beamnr]) # Rarely used
	)
	
	# Set the approximation funciton for the beam pattern:
	if(data[[pbp_name]] %in% c("circularPiston_ellipticRadius_sidelobefit", "circularPiston_ellipticRadius", "circularPiston")){
		appr = circularPiston_appr
	}
	# If circular piston is not used, the beam pattern needs to be linked to beam index number 'indi':
	else{
		appr = NULL
	}

# Run the beam pattern function and return:
 data[[bpt_name]](beamPatternInput, appr=appr)	
}

#echoIBM_setB_T <- function(data, direction=1, beamnr=1){
#	
#	# Set the names of the beam pattern function, the radial dimension and the parametric beam pattern:
#	bpt_name <- paste0("bpt", direction)
#	rad_name <- paste0("rad", direction)
#	pbp_name <- paste0("pbp", direction)
#	
#	if(data[[pbp_name]] == "circularPiston_ellipticRadius_sidelobefit"){
#		# Use as input the 'dira' [lthesel], 'dire' [lthesel], 'wnsz' [1,2] and 'sllf' [1,2]:
#		data[[bpt_name]](list(
#			dira = data$fishdirT[,2], 
#			dire = data$fishdirT[,3], 
#			wnsz = data$uniquek[i] * data$[[rad_name]][beamnr,,drop=FALSE], 
#			sllf = data$sllf[beamnr,,drop=FALSE]), 
#			appr = circularPiston_appr
#		)
#	}
#	else if(data$[[pbp_name]] == "circularPiston_ellipticRadius"){
#	# Use as input the 'dira' [lthesel], 'dire' [lthesel] and 'wnsz' [1,2]:
#	 data[[bpt_name]](list(
#		 dira = data$fishdirT[,2], 
#		 dire = data$fishdirT[,3], 
#		 wnsz = data$uniquek[i] * data$[[rad_name]][beamnr,,drop=FALSE]), 
#		 appr = circularPiston_appr
#		 )
#	 }
#	else if(data$[[pbp_name]] == "circularPiston"){
#		# Use as input the 'dire' [lthesel] and 'wnsz' [1,2]:
#		 data[[bpt_name]](list(
#			 dira = data$fishdirT[,2], 
#			 wnsz = data$uniquek[i] * data$[[rad_name]][beamnr,,drop=FALSE]), 
#			 appr = circularPiston_appr
#		 )
#	}
#	else{
#		# Use as input the 'dira' [lthesel], 'dire' [lthesel], 'indi' [1]:
#		 data[[bpt_name]](list(
#			 dira = data$fishdirT[,2], 
#			 dire = data$fishdirT[,3], 
#			 indi = data$indi[beamnr])
#		 )
#	}
#}