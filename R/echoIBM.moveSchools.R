#*********************************************
#*********************************************
#' Moves compactly stored school data by time.
#'
#' If only one time step is given in the compactly stored school data, calculate new positions of school given time elapsed. If multiple time steps are given, interpolate all variables at the requested time 'utim'. 
#'
#' @param data	The compatly stored school data.
#' @param utim	The time in UNIX time.
#'
#' @return Modified 'data'.
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
	# Funciton that interpolates one variable in 'data' by the given time 'utim'. Time is assumed organized along the second dimension:
	approxOne <- function(var, data, utim){
		out <- double(nrow(data[[var]]))
		for(i in seq_along(out)){
			out[i] <- approx(data$utmS[i,], data[[var]][i,], utim)$y
		}
		out
	}
	
	# Modify the school data only if the UNIX time is present:
	if(length(data$utmS)){
		if(NCOL(data$psxS)==1){
			timespan <- utim - data$utmS
			xyz <- cbind(data$psxS, data$psyS, data$pszS) + sph2car(cbind(timespan * data$aspS, data$thtS, data$phiS))
			data$psxS <- xyz[,1]
			data$psyS <- xyz[,2]
			data$pszS <- xyz[,3]
		}
		else{
			notUNIXnames <- setdiff(names(data), "utmS")
			data[notUNIXnames] <- lapply(notUNIXnames, approxOne, data=data, utim=utim)
		}
	}
	
	
	return(data)
	
	### ##########################################################################################
	### ########## Modify this to accept multiple time steps specified for each school. ##########
	### ##########################################################################################
	### approxOne <- function(y, utmS, utim){
	### 	ndim <- length(dim(y))
	### 	if(ndim==2){
	### 		# Interpolate the information in a matrix with schools along the rows and time steps along the columns [#schools, #timesteps]:
	### 		apply(y, 1, function(yy) approx(utmS, y=yy, xout=utim)$y)
	### 	}
	### 	else{
	### 		approx(x=utmS, y=y, xout=utim)$y
	### 	}
	### }
	### # If only one time step is given in the compactly stored schools, return the data from this time step if the speed 'ispS' is not given, and move the school otherwise:
	### if(length(data$utmS)==1){
	### 	#if(length(data$ispS)==0 || any(is.na(data$ispS))){
	### 		timespan <- utim - data$utmS
	### 		xyz <- cbind(data$psxS, data$psyS, data$pszS) + sph2car(cbind(timespan * data$aspS, data$thtS, data$phiS))
	### 		data$psxS <- xyz[,1]
	### 		data$psyS <- xyz[,2]
	### 		data$pszS <- xyz[,3]
	### 		data
	### 	#}
	### }
	### else{
	### 	notUtim <- names(data) != "utmS"
	### 	data[notUtim] <- lapply(data[notUtim], approxOne, data$utmS, utim)
	### }
	### return(data)
	### ### timespan = utim - data$utmS
	### ### xyz = cbind(data$psxS, data$psyS, data$pszS) + sph2car(cbind(timespan*data$aspS, data$thtS, data$phiS))
	### ### data$psxS = xyz[,1]
	### ### data$psyS = xyz[,2]
	### ### data$pszS = xyz[,3]
	### ### data
	}
