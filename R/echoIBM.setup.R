#*********************************************
#*********************************************
#' Prepare for simulation of underwater acoustic data using \code{\link{echoIBM}}.
#'
#' @param event	The name of the event, or the path to all sub events event/esnm/"tsd", in which case \code{dir} is ignored.
#' @param dir		Directory in which to place the event.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR read.event
#' @importFrom TSD car2global catchWarnings ftim2utim NAs ones write.TSD zeros
#' @importFrom stats runif
#'
#' @export
#' @rdname echoIBM_SX90_BlindZone
#'
echoIBM.setup <- function(
	event,
	dir = NULL,
	esnm = "EK60", 
	t = "all",
	pre = TRUE, 
	vessel = list(), 
	school = list(), 
	ctd = FALSE, 
	beams = list(), 
	noise = list(), 
	calibration = list()
	){

	# 0. Create the events
	if(!is.list(event) && !c("path", "esnm", "name") %in% names(event)){
		event <- createEvent(event, dir, esnm, ow=FALSE)
	}
	files <- list()
	
	# Set the global variables:
	globalVar <- list(event=event)
	
	
	
	#######################
	##### C1. Vessel: #####
	#######################
	if(!identical(vessel, FALSE)){
		files$vessel <- do.call(echoIBM.setVessel, c(globalVar, vessel))
	}
	#######################
	
	
	######################
	##### S1. Beams: #####
	######################
	if(!identical(beams, FALSE)){
		files$beams <- do.call(echoIBM.setBeams, c(globalVar, beams))
	}
	######################
	
	
	#######################
	##### C2. School: #####
	#######################
	if(!identical(school, FALSE)){
		files$school <- do.call(echoIBM.setSchool, c(globalVar, school))
	}
	#######################
	
	
	
	####################
	##### C3. CTD: #####
	####################
	if(!identical(ctd, FALSE)){
		files$ctd <- do.call(echoIBM.setCTD, c(globalVar, ctd))
	}
	####################
	
	
	
	############################################################
	########## Specific for the acoustic instruments: ##########
	############################################################
	
	##### 5. Noise: #####
	if(!identical(noise, FALSE)){
		files$noise <- do.call(echoIBM.setNoise, c(globalVar, noise))
	}
			
	##### 6. Calibration: #####
	if(!identical(calibration, FALSE)){
		files$calibration <- do.call(echoIBM.setCalibration, c(globalVar, calibration))
	}
	############################################################
	############################################################
	
	
	############################################################
	########### Common for all acoustic instruments: ###########
	############################################################
	
	############################################################
	############################################################
	
	
	
	return(list(event=event, files=files))



	
	### ##################################
	### ######### (5) Simulation: ########
	### ##################################
	### # Create dumpfile directories:
	### dumpfiledir = file.path(events, "dumpfiles")
	### for(i in seq_along(dumpfiledir)){
	### 	suppressWarnings(dir.create(dumpfiledir[i]))
	### 	}
	### dumpfilename = file.path(dumpfiledir,"dump.txt")
    ### 
	### 
	### # Add noise if specified:
	### noise <- if(addnoise) c("ms","ca","bg") else "ms"
	### 
	### 
	### ###########
	### ## SU90: ##
	### ###########
	### if("SU90" %in% esnm){
	### 	# Simulatte in parallel:
	### 	event = c(events["SU90"==esnm], if(length(events)>1 && "SU90"!=esnm[1]) events[1])
	### 	system.time(warn <- catchWarnings(echoIBM(event=event, t=t, noise=noise, max.memory=2e9, dumpfile=dumpfilename["SU90"==esnm], parlist=list(pre=TRUE, seed="ind"), rand.sel=10, cores=cores)))
	### 	#    user   system  elapsed 
	### 	# 216.750   76.885 4722.775 
	### 	# Only warnings related to not writing variables with not extactly four characters in the variable names:
	### 	echoIBM.addwarnings(dumpfilename["SU90"==esnm], warn)
	### 	}
	### ###########
    ### 
	### ###########
	### ## EK60: ##
	### ###########
	### if("EK60" %in% esnm){
	### 	# Simulatte in parallel:
	### 	event = c(events["EK60"==esnm], if(length(events)>1 && "EK60"!=esnm[1]) events[1])
	### 	system.time(warn <- catchWarnings(echoIBM(event, t=t, noise=noise, max.memory=2e9, dumpfile=dumpfilename["EK60"==esnm], parlist=list(pre=TRUE, seed="ind"), max.radius=0.4, cores=cores)))
	### 	#   user  system elapsed 
	### 	# 44.517   7.256 996.010 
	### 	# Only warnings related to not writing variables with not extactly four characters in the variable names:
	### 	echoIBM.addwarnings(dumpfilename["EK60"==esnm], warn)
	### 	}
	### ###########
	### events
	}
