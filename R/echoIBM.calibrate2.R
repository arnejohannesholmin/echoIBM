#*********************************************
#*********************************************
#' Calibrates (and writes the calibration data to file) the simulation of observations from the MS70 sonar performed by echoIBM(). Not used in any of the other echoIBM code files.
#'
#' @param directory  the path to the directory in which to put the calibration file (use TRUE for the default directory file.path(echoIBM_datasets_directory(), "Resources", "Calibration")).
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param type  is "sv" for calibrating the simulation model with respect to volume scattering 'sv' and "ts" for calibrating the simulation model with respect to backscattering cross section / target strength.
#' @param runs  is an integer giving the number of simulated pings of which the mean is taken as the final calibration factors.
#' @param add  is the angle in degrees to add on either end of the insonifying volume of the acoustical instrument, defining the volume in which point targets are distributed.
#' @param N  is the number of point targets used in the calibration.
#' @param j  is the sampling interval to use for the calibration (should be larger than 2).
#' @param sgbs  is the backscattering cross section of the point targets.
#' @param max.memory  is used in echoIBM to set the maximum memory for echoIBM to occupy.
#' @param filetag  is a string to add to the file name.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom sonR get.specs.esnm read.event runif.sph
#' @importFrom TSD NAs sph2car strff write.TSD zeros
#'
#' @export
#' @rdname echoIBM.calibrate
#'
echoIBM.calibrate2 <- function(directory=NULL, event=1, type=c("sv","ts"), dire=NULL, runs=1, add=30, N=1e5, j=100, sgbs=1e-6, max.memory=1e9, method=c("closest","linear"), max.radius=0.2, seed=0, cores=1, usemean=FALSE, nameadd=""){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-03-20 - Clean version.
	# Update: 2011-06-25 - Major changes: The option 'type' added, for choosing between calibrating for sv or for TS. A factor 1/(4*pi) added to the true 'sv'. File written to the directory 'directory', while the previous version only had the option to write to a specific file. The file name is as follows: paste("CalibrationTable_",toupper(esnm),"_Cruise_",cruise,"_sigma_bs.beams",sep="").
	# Update: 2011-06-26 - Added the option 'runs', for running multiple simulations of fewer calibration targets (saving memory).
	# Update: 2011-08-24 - Changed from using 'epsl', 'acca' and 'ssil' to using 'epss' and 'sgbs'. Tested using the script "Test echoIBM.calibrate.R".
	# Update: 2012-02-13 - Changed the method radically, by using echoIBM exactly in the way an ordinary simulation would be performed, in order to keep the same control over changes in the simulation model without making the same changes in this function (the previous version used echoIBM.oneping.oneschool). The data are then read back in to the function, and the calibration factors calculated.
	# Update: 2013-02-25 - Saved some minor changes several months ago.
	# Last: 2015-02-22 - Added support for multiple tilt angles through 'dire'.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Function for generating school data for all runs (pings):
	generateSchool <- function(data, caldir, calname, thisdire=list(), add=30, sgbs=1e-6, N=1e5){
		region <- get.specs.esnm(c(thisdire, data), add=add)$region
		#dr <- data$asps*data$sint/2
		r0 <- (j - 2) * data$rres
		r1 <- j * data$rres
		# Total volume in which the target positions are drawn. Also see appendix 2 of the documentation of echoIBM:
		V <- vol.sph(r=c(r0,r1), theta=region[,1], phi=region[,2])$volx
		# The true volume backscatter, equal to N*sigma_bs/V = N*A_0/V for spherical targets:
		truesv <- sgbs * N / V
	
		# Generate the target positions and write to the school file at each run:
		cat("Write the dynamic school data...\n")
		schooldynamicFile <- file.path(caldir, paste0(calname, "_pointTargets_dynamic.school"))
	
		for(i in seq_len(runs)){
			cat("Run ", i, ":\t", sep="")
			### Generate positions corresponding to uniformly distributed objects in space, from which to simulate observations from the MS70 sonar. See appendix 2 of the documentation of echoIBM: ###
			xyz <- runif.sph(N, r=c(r0, r1), theta=region[,1], phi=region[,2], sph.out=FALSE)
			#plot3d(xyz,size=0.3)
		
			# Add transducer position and cut the positions above the sea level:
			xyz[,3] <- xyz[,3] + data$psze[1]
			xyz <- xyz[xyz[,3] < 0,]
			# Update 'N':
			thisN <- nrow(xyz)
			cat(thisN, "point sources drawn\n")
		
			# Add the target information (all rotations left out because of points targets):
			schooldynamic <- list(psxf=xyz[,1], psyf=xyz[,2], pszf=xyz[,3], utim=i - 1)
		
			# Write the dynamic school data to file:
			if(i==1){
				write.TSD(schooldynamic, schooldynamicFile, reserve=runs-1, numt=1)
			}
			else{
				write.TSD(schooldynamic, schooldynamicFile, numt=1, append=TRUE)
			}
		}
		return(truesv)
	}
	
	if(isTRUE(directory)){
		directory <- file.path(echoIBM_datasets_directory(), "Resources", "Calibration")
	}
	
	# 'beams' variables must be given:
	event <- event.path(event)
	data <- read.event(var="beams", event=event$event)
	data$rres <- soundbeam_range(data, "rres")
	if(j - 2 > min(data$lenb)){
		stop(paste0("The sampling interval index 'j' chosen too large for the minimum beam length ", min(data$lenb)))
	}
	
	# Get the number of beams:
	Ni <- length(data$dira)
		
		
	########## Execution and output ##########
	# The calibration procedure of echoIBM for volume backscattering coefficient 'sv' is as follows: Draw 'N' uniformly distributed spherical targets of backscattering cross section 'sgbs' covering the volume affecting all voxels corresponding to sampling interval 'j' (see appendix 2 of the documentation of echoIBM). Define the true 'sv' to be N*sgbs/V, where 'V' is the size of the volume covered by sampling interval 'j' in the system of voxels. Then run the simulations and compare to the true 'sv' by defining calibration factors equal to true 'sv' / colSums(simulated 'sv').
	
	### (1) Create the directory of the calibration simulation: ###
	calname <- paste("Calibration", data$esnm, "Event", event$eventname, "sv", sep="_")
	
	caldir <- file.path(tempdir(), "Calibration", "Events", calname, data$esnm, "tsd")
	unlink(caldir, recursive=TRUE, force=TRUE)
	dir.create(caldir,recursive=TRUE)
		
	
	### (2) Write the beam-configuration data: ###
	cat("Write the beam-configuration data...\n")
	beamsfile <- file.path(caldir, paste(calname, "beams", sep="."))
	write.TSD(data, beamsfile, numt=1)
	
	
	### (3) Write the vessel data: ###
	cat("Write the vessel data...\n")
	vessel <- list(psxv=zeros(runs), psyv=zeros(runs), pszv=zeros(runs), rtzv=zeros(runs), utim=sequence(runs) - 1)
	vesselfile <- file.path(caldir, paste(calname, "vessel", sep="."))
	write.TSD(vessel, vesselfile, numt=runs)
	
	
	### (4) Write the static school data: ###
	schoolstatic <- list(pbpf="ps", sgbs=sgbs, obln=1, gamw=0, gaml=0)
	cat("Write the static school data...\n")	
	schoolstaticFile <- file.path(caldir, paste0(calname, "_pointTargets_static.school"))
	write.TSD(schoolstatic, schoolstaticFile, numt=1)
	
	# If 'dire' is given, run through its elements and calibrate for each value, equal for all beams:
	ndire <- max(1, length(dire))
	cali <- NAs(Ni, ndire)
	
	for(k in seq_len(ndire)){
		
		if(ndire==1){
			cat("Calibration ...\n")
			# Use the input elevation angle:
			thisdire <- list()
		}
		else{
			cat("Calibration for tilt angle ", k, " (of", length(dire), ") ...\n")
			# Use the input elevation angle:
			thisdire <- list(dire=rep(dire[k], length(data$lenb)))
		}
		
		### (5) Generate the target positions and write the school file: ###
		truesv <- generateSchool(data=data, caldir=caldir, calname=calname, thisdire=thisdire, add=add, sgbs=sgbs, N=N)
		
		### (6) Run the simulation: ###
		echoIBM(caldir, adds=thisdire, t="all", TVG.exp=2, calibrate=FALSE, noise="", max.memory=max.memory, max.radius=max.radius, method=method, parlist=list(seed=seed), keep.temp=FALSE)
		
		### (7) Read the simulated data: ###
		sv <- zeros(max(data$lenb), length(data$lenb))
		for(i in seq_len(runs)){
			thissv <- read.event(event=caldir, t=i, var="vbsc")$vbsc
			sv <- sv + thissv
		}	
				
		### (8) Store the calibration factors of the current tilt: ###
		cali[,k] <- truesv / sv[j,] * runs
	}
	
	### (9) Write and return the calibration data: ###
	# Apply the mean of all beams if usemean==TRUE (omnidirectional sonars such as Simrad SX90):
	if(usemean){
		out <- c(list(cali=apply(cali, 2, function(xx) rep(mean(xx), length(xx))), cal0=cali, grde=dire, ctim=unclass(Sys.time())))
	}
	else{
		out <- c(list(cali=cali, grde=dire, ctim=unclass(Sys.time())))
	}
	# Write the calibration data:
	if(length(directory)>0 && !identical(directory,FALSE)){
		con <- file.path(directory, paste("CalibrationTable_", toupper(data$esnm), "_Event_", event$eventname, if(nchar(nameadd)>0) "_", nameadd, "_sv.beams", sep=""))
		write.TSD(out, con=con, numt=1)
	}
	
	# Return the calibration data:
	c(out, list(truesv=truesv, caldir=caldir))
	##################################################
	##################################################
}
