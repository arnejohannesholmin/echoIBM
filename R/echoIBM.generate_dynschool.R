#*********************************************
#*********************************************
#' Generates the dynamic variables of schools positioned close enough to enter the sampling volume of the sonar.
#'
#' @param data  is a list containing the necessary data.
#' @param t  is the current time step.
#' @param vesselutim  is a vector of the vessel utim points (all time steps in the event)
#' @param margin  is an extra margin to add, so that the school is detected outside of (positive value) or or inside of (negative value) of the sampling volume.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom sonR sonR_implemented
#' @importFrom TSD labl.TSD NAs prettyIntegers
#' @importFrom utils write.table
#'
#' @export
#' @rdname echoIBM.generate_dynschool
#'
echoIBM.generate_dynschool <- function(data, t=1, vesselutim=NULL, adds=list(), margin=0, dumpfile=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-02-05 - First version.
	########### DESCRIPTION: ###########
	# Generates the dynamic variables of schools positioned close enough to enter the sampling volume of the sonar.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the necessary data.
	# ---t--- is the current time step.
	# ---vesselutim--- is a vector of the vessel utim points (all time steps in the event)
	# ---margin--- is an extra margin to add, so that the school is detected outside of (positive value) or or inside of (negative value) of the sampling volume.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Check for the required variables:
	requiredVariables <- c(
	"utmS", "aspS"
	,"psxS", "psyS", "pszS"
	,"thtS", "phiS"
	,"szxS", "szyS", "szzS"
	#,"rtxS","rtyS","rtzS"
	#,"shpS"
	)
	densityVariables <- c("rhoS", "nbfS")
	sizeVariables <- c("MEsz", "SDsz", "seed")
	polarizationVariables <- c("plHS", "SDxf", "SDyf", "SDzf", "sEdx", "sEdy", "sEdz")
	if(!all(requiredVariables %in% names(data))){
		warning(paste0("All of the school specification variables (", paste0(requiredVariables,collapse=", "), ") must be present in the data. No schools generated."))
		}
	if(!any(densityVariables %in% names(data))){
		warning(paste0("At least one of ", paste(densityVariables,collapse=", "), " must be present in the data. No schools generated."))
		}
	if(!all(sizeVariables %in% names(data))){
		warning(paste0("All of the fish size specification variables (", paste(sizeVariables,collapse=", "), ") must be present in the data. No schools generated."))
		}
	if(!polarizationVariables[1] %in% names(data) && !all(polarizationVariables[-1] %in% names(data))){
		warning(paste0("Either ", polarizationVariables[1], " or all of ", paste(densityVariables,collapse=", "), " must be present in the data. No schools generated."))
		}
	#range <- data$sint*data$asps/2*max(data$lenb)
	range <- soundbeam_range(data, pos="max")
	### # Get the horizontal range of the sonar/echosounder:
	###if(sonR_implemented(data, "SBE")){
	###	# Assume down to 500 m with a 7 degree opening angle (sin(3.5*pi/180) * 500), and add some to cover further than the -3 dB-angle:
	###	range <- 50
	###	}
	###if(sonR_implemented(data, "MBE")){
	###	# Assume 45 degrees angle of the upper beams:
	###	range <- range * sin(45 / 180*pi)
	###	}
		
	# Special case if seeds were given using sEds:
	if(!length(data$seed) && length(data$sEds)){
		data$seed <- data$sEds
	}	
	
		
	########## Execution ##########
	# Check whether any of the schools are within range of the sonar:
	# Get the current positions of the schools:
	if(t>1){
		timediff <- vesselutim[t] - vesselutim[t-1]
		}
	else if(length(vesselutim)>1){
		timediff <- vesselutim[t+1] - vesselutim[t]
		}
	else{
		timediff <- 1
		}
	
	# Move the schools to match the time of the current ping:
	data <- echoIBM.moveSchools(data, vesselutim[t])
	# This is crude, since the case that a very large school is present will result in a large value of 'maxsize'. However, this is just to get the schools that are considered in each ping, and including more schools than are inside the sampling volume causes no harm:
	maxsize <- apply(cbind(data$szxS, data$szyS, data$szzS), 1, max)
	
	# Get the distance to the schools along the sea surface:
	distToSchools <- sqrt((data$psxS - data$psxv)^2 + (data$psyS - data$psyv)^2)
	inside <- distToSchools < range + margin + maxsize
	### # Get the schools that are valid for the current time step:
	### insideTime <- vesselutim[t] >= data$utmS
	### if(length(data$ut9S)>0){
	### 	insideTime <- insideTime & vesselutim[t] <= data$ut9S
	### 	}
	### # Discard the schools for which the current time is outside of the time range of the schools given by 'utmS' (start time) and 'ut9S' (end time):
	### inside <- which(inside & insideTime)
	
	# For the schools inside the volume relevant to sampling, draw fish positions and sizes:
	#tempdynschool <- list()
	#if(length(data$rhoS)==0){
	#	data$rhoS <- NAs(length(data$rhoS))
	#	}
	#for(i in seq_along(inside)){
	#	thisdata <- echoIBM.generate_oneschool(data,school=inside[i],dt=timediff,seed=data$seed+t,dumpfile=dumpfile)
	#	tempdynschool[paste(names(thisdata),inside[i],sep="")]=thisdata
	#	# Update 'rhoS' (in the case that it was not given), and 'nbfS':
	#	data$rhoS[inside[i]] <- thisdata$rhoS[inside[i]]
	#	data$nbfS[inside[i]] <- thisdata$nbfS[inside[i]]
	#	data$plHS <- data$plHS*180/pi
	#	}
	
	# For some reason the unlist method above caused time lag, so we instead try appending:
	# For the schools inside the volume relevant to sampling, draw fish positions and sizes:
	echoIBM.generate_oneschool_labl <- labl.TSD("echoIBM.generate_oneschool_labl")
	dynschool <- vector("list", length(echoIBM.generate_oneschool_labl))
	names(dynschool) <- echoIBM.generate_oneschool_labl
	if(length(data$rhoS)==0){
		data$rhoS <- NAs(length(inside))
		}
	if(length(data$nbfS)==0){
		data$nbfS <- NAs(length(inside))
		}
	for(i in seq_along(inside)){
		thisdata <- echoIBM.generate_oneschool(data, school=inside[i], dt=timediff, dumpfile=dumpfile)
		for(j in seq_along(echoIBM.generate_oneschool_labl)){
			dynschool[[echoIBM.generate_oneschool_labl[j]]] <- c(dynschool[[echoIBM.generate_oneschool_labl[j]]], thisdata[[echoIBM.generate_oneschool_labl[j]]])
			}
		# Update 'rhoS' (in the case that it was not given), and 'nbfS':
		data$rhoS[inside[i]] <- thisdata$rhoS[inside[i]]
		data$nbfS[inside[i]] <- thisdata$nbfS[inside[i]]
		}
	
	# Convert to degrees:
	data$plHS <- data$plHS * 180/pi
	
	# Write school information to the dumpfile:
	cat("The following schools were simulated at the current time step (",t,"): ", prettyIntegers(inside), "\n", sep="")
	if(length(dumpfile)>0){
		cat("\nThe following schools were simulated at the current time step (", t, "): \n\n", sep="", file=dumpfile, append=TRUE)
		if(any(inside)){
			write.table(t(round(as.data.frame(c(list(indS=inside), lapply(data[c("utmS", "aspS", "psxS", "psyS", "pszS", "thtS", "phiS", "szxS", "szyS", "szzS", "volS", "rhoS", "nbfS", "MEsz", "SDsz", "seed", "plHS")], "[", inside))), digits=1)), dumpfile, append=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
			}
		else{
			cat("none", file=dumpfile, append=TRUE)
			}
		cat("\n\n", file=dumpfile, append=TRUE)
		}
	
	# Apply the school information in 'adds':
	dynschoolLengths <- unlist(lapply(dynschool, length))
	if(length(adds)>0 && length(dynschoolLengths)>0 && !any(dynschoolLengths==0)){
		# Add the data in 'adds' (overrides the data read in the school files):
		dynschoolinadds <- intersect(names(adds), names(dynschool))
		Nl <- max(dynschoolLengths)
		dynschool[dynschoolinadds] <- lapply(adds[dynschoolinadds], rep, length.out=Nl)
		}
	
	
	########## Output ##########
	dynschool
	##################################################
	##################################################
	}