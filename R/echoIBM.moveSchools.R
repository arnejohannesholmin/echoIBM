#*********************************************
#*********************************************
#' Calculates new positions of school specified compactly by the dynamic variables of the school as a whole, given time elapsed. 
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD sph2car
#' @importFrom stats approx
#'
#' @export
#' @rdname echoIBM.moveSchools
#' 
echoIBM.moveSchools <- function(data, utim){
	##########################################################################################
	########## Modify this to accept multiple time steps specified for each school. ##########
	##########################################################################################
	approxOne <- function(x, utmS, utim){
		ndim <- length(dim(x))
		if(ndim==2){
			apply(x, 1, function(xx) approx(utmS, xx, utim)$y)
		}
		else{
			approx(utmS, x, utim)$y
		}
	}
	# If only one time step is given in the compactly stored schools, return the data from this time step if the speed 'ispS' is not given, and move the school otherwise:
	if(length(data$utmS)==1){
		if(length(data$ispS)==0 || any(is.na(data$ispS))){
			timespan <- utim - data$utmS
			xyz <- cbind(data$psxS, data$psyS, data$pszS) + sph2car(cbind(timespan*data$aspS, data$thtS, data$phiS))
			data$psxS <- xyz[,1]
			data$psyS <- xyz[,2]
			data$pszS <- xyz[,3]
			data
		}
	}
	else{
		notUtim <- names(data) != "utmS"
		data[notUtim] <- lapply(data[notUtim], approxOne, data$utmS, utim)
	}
	return(data)
	### timespan = utim - data$utmS
	### xyz = cbind(data$psxS, data$psyS, data$pszS) + sph2car(cbind(timespan*data$aspS, data$thtS, data$phiS))
	### data$psxS = xyz[,1]
	### data$psyS = xyz[,2]
	### data$pszS = xyz[,3]
	### data
	}
