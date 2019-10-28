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
	
	############### LOG: ###############
	# Start: 2014-02-05 - First version.
	
	##################################################
	##################################################
	########## Preparation ##########
	# Check for the required variables:
	requiredVariables <- c(
	"utmS", "aspS",
	"psxS", "psyS", "pszS",
	#"thtS", "phiS",
	"szxS", "szyS", "szzS"
	#,"rtxS","rtyS","rtzS"
	#,"shpS"
	)
	densityVariables <- c("rhoS", "nbfS")
	sizeVariables <- c("MEsz", "SDsz")
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
	
	# Repeat varables that should be given for each time step (such as rhoS):
	Nschools <- length(data$psxS)
	echoIBM.generate_oneschool_labl <- labl.TSD("echoIBM.generate_oneschool_labl")
	data[echoIBM.generate_oneschool_labl] <- lapply(data[echoIBM.generate_oneschool_labl], function(x) if(length(x)==1) rep(x, Nschools) else x)
		
	# Special case if seeds were given using sEds:
	if(!length(data$seed) && length(data$sEds)){
		data$seed <- data$sEds
	}	
	
		
	########## Execution ##########
	# Move the schools to match the time of the current ping:
	data <- echoIBM.moveSchools(data, vesselutim[t])
	
	
	
	# This is crude, since the case that a very large school is present will result in a large value of 'maxsize'. However, this is just to get the schools that are considered in each ping, and including more schools than are inside the sampling volume causes no harm:
	maxsize <- apply(cbind(data$szxS, data$szyS, data$szzS), 1, max)
	
	# Get the distance to the schools along the sea surface:
	distToSchools <- sqrt((data$psxS - c(data$psxv))^2 + (data$psyS - c(data$psyv))^2)
	inside <- distToSchools < range + margin + maxsize
	
	# Get the schools that are valid for the current time step:
	insideTime <- vesselutim[t] >= data$utmS
	if(length(data$ut9S)>0){
		insideTime <- insideTime & vesselutim[t] <= data$ut9S
		}
	# Discard the schools for which the current time is outside of the time range of the schools given by 'utmS' (start time) and 'ut9S' (end time):
	inside <- which(inside & insideTime)
	
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
	
	# For some reason the unlist method above caused time lag, so we instead try appending. Not so many schools inside the sonar range at each time anyhow:
	# For the schools inside the volume relevant to sampling, draw fish positions and sizes:
	dynVars <- labl.TSD("echoIBM_dynVars")
	dynschool <- vector("list", length(dynVars))
	names(dynschool) <- dynVars
	if(length(data$rhoS)==0){
		data$rhoS <- NAs(length(inside))
	}
	if(length(data$nbfS)==0){
		data$nbfS <- NAs(length(inside))
	}
	for(i in seq_along(inside)){
		thisdata <- echoIBM.generate_oneschool(data, school=inside[i], dumpfile=dumpfile)
		if(length(thisdata)==0){
			inside[i] <- NA
			next
		}
		for(j in seq_along(dynVars)){
			dynschool[[dynVars[j]]] <- c(dynschool[[dynVars[j]]], thisdata[[dynVars[j]]])
		}
		
		# Update 'rhoS' (in the case that it was not given), and 'nbfS':
		data$rhoS[inside[i]] <- thisdata$rhoS[inside[i]]
		data$nbfS[inside[i]] <- thisdata$nbfS[inside[i]]
	}
	
	# Subset 'data' to the school inside the observation range:
	inside <- inside[!is.na(inside)]
	#data <- lapply(data[intersect(names(data), labl.TSD("cs"))], "[", inside)
	
	data <- lapply(data[intersect(names(data), labl.TSD("cs"))], function(x) if(length(x) > 1) x[inside] else x)
	
	
	
	
	
	
	# Add the data stored in 'data' to the individual fish dynamics:
	dynschool <- c(dynschool, data)
	
	# Apply the school information in 'adds':
	dynschoolLengths <- unlist(lapply(dynschool, length))
	if(length(adds)>0 && length(dynschoolLengths)>0 && !any(dynschoolLengths==0)){
		# Add the data in 'adds' (overrides the data read in the school files):
		dynschoolInAdds <- intersect(names(adds), names(dynschool))
		Nl <- max(dynschoolLengths)
		dynschool[dynschoolInAdds] <- lapply(adds[dynschoolInAdds], rep, length.out=Nl)
	}
	
	# Write school information to the dumpfile. Here we do not dump the generated fish positions, but only the compact school variables, where number of fish and packing density has been updated above:
	if(length(inside)){
		cat("The following schools were simulated at the current time step (", t, "): ", prettyIntegers(inside), "\n", sep="")
	}
	if(length(dumpfile)>0){
		# Convert to degrees in the dumped output:
		data$plHS <- data$plHS * 180/pi
	
		cat("\nThe following schools were simulated at the current time step (", t, "): \n\n", sep="", file=dumpfile, append=TRUE)
		if(any(inside)){
			write.table(t(as.data.frame(c(list(indS=inside), data))), dumpfile, append=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
		}
		else{
			cat("none", file=dumpfile, append=TRUE)
		}
		cat("\n\n", file=dumpfile, append=TRUE)
	}
	
	
	########## Output ##########
	dynschool
	##################################################
	##################################################
}
