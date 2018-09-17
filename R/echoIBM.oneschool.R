#*********************************************
#*********************************************
#' Simulates echo sounder observations based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school.
#'
#' @param files  is a vector of the paths to the files on which the simulation is based.
#' @param t  is a vector of the numbers of the pings to be simulated, as listed from 1 to the number of time steps of the simulated school. If t=="all", all time steps are simulated.
#' @param pingsSchool  is a list of four-column matrices:
#' @param adds  is an optional list of variables overriding the variables located in the files given by 'files'.
#' @param pingsdir  is a character string holding the path to the directory in which the .pings files should be put.
#' @param pingsname  is a character string holding the name of the .pings files, after which the time point of the first ping of each file will be added.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#' @param compensated  specifies which rotations are compensated for by the echo sounder:
#' @param filesize  is the maximum size of the files, used both for temporal and final files.
#' @param calibrate  is FALSE if calibration data are to be discarded if present.
#' @param noise  See echoIBM.add.noise().
#' @param max.memory  is the maximum amount of memory to allow for the function to occupy. The level of for loops is chosen according to this value, i.e. if the method using only two loops, one for the radial distances and one for the unique frequencies, demands more memory than 'max.memory', the method usint three for loops is chosen.
#' @param ow  is TRUE if the user wish to overwrite existing file.
#' @param origin  is a vector of two elements representing the origin of the global coordinate system (G), or the numbering index of the ping regarded as the origin of (G) (ignoring heave so that the x-y-plane of (G) is on the surface of the sea).
#' @param dumpfile  is the name of the file to which information and warnings about the simulation is written.
#' @param rand.sel  is a numeric specifying a random selection of the school to use in the simulations. If 'rand.sel' is not >0 and <1 no random selection is done.	If given as a vector of length 2, the second element is regarded as the seed for the random selection.
#' @param method  is "closest" if the beam pattern value of the closest grid point is to be selected, and "linear" if linear interpolation should be used to extract the beam pattern value (time demanding and only available for 2D grids).
#' @param ask  is TRUE if the used should be asked to for approval if the memory of the least memory demanding calculation method of the individual radial sampling intervals exceed the memory limit 'max.memory'.
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#' @param bptfile  is the name of the file to which 'sllf', 'rad1', 'rad2', 'pbp1' and 'pbp2' is written (use NULL for no writing).
#' @param cores  is an integer specifying the number of cores to run the simulations over in parallel (should be lower than the number of cores in the computer).
#' @param saveBrkt  is TRUE to save the variables 'Brkt' and 'nsig', which contain data only when the number of significant targets is low enough (<=2) for the exact Barakat distribution to kick in instead of the approximation of the exponential distribution.
#' @param msg				Logical: If TRUE, show messages to the console and to dump files.
#' @param discardOutside	A vector of three elements specifying the size of the margin around the -3 dB main beam, specifide in units of the range/angular resolution, inside which targets are considered. Using e.g. c(2, 3, 3) will discard targest outside of a volume extending three beam widths around the beam width of the system (and to 2 times the radial range beyond the maximum range). This may save process time, but is not recommended without testing, as higher order side lobes encounter a larger volume than lower order side lobes.
#' @param fishReaction		A list of variables specifying the fish reaction to the vessel. Must contain the following variables: "magR", "mxrR", "mxdR", "facR". See info.TSD(c("magR", "mxrR", "mxdR", "facR")) for definition of these variables. 
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom sonR echoIBM.getSchoolfileType read.event
#' @importFrom TSD ftim.TSD global2car labl.TSD prettyIntegers read.TSD read.TSDs splitSeqIntoBlocks strff write.TSD zeropad mergeTextFiles
#' @importFrom utils head
#'
#' @export
#' @rdname echoIBM.oneschool
#'
echoIBM.oneschool <- function(files, event, t=1, tvessel, vesselutim, pingsSchool, areCompact, adds=NULL, pingsdir, pingsname, esnm=NULL, TVG.exp=2, compensated=c("pitch", "roll"), filesize=3e8, calibrate=TRUE, noise=c("nr", "bg", "ex"), mode=c("active", "passive"), max.memory=1e9, ow=TRUE, origin=1, dumpfile="dump.txt", dumpsize=10e6, timedumpfile="timedump.txt", rand.sel=1, scls=1, volsc=c(1, 1, 1), method=c("closest", "linear"), ask=FALSE, parlist=list(), bptfile=NULL, max.radius=0.2, cores=1, saveBrkt=FALSE, msg=FALSE, discardOutside=c(r=Inf,az=Inf,el=Inf), fishReaction=list()){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-09-30 - First version.
	# Update: 2010-03-16 - Added exponential noise implementation.
	# Update: 2010-03-16 - Separated from the old function named echoSim(), and simulates the echo from one school based on the files given by 'files', which need to include ".beams"-files of the acoustic instrument, calibration and noise, in addition to ".ctd"-files and ".vessel"-files.
	# Update: 2011-01-11 - Changed the method of selecting the apropriate time steps for the relevant files. Earlier a two column matrix named 'filet' was created, holding the file numbers and time step numbers to read. Now the method is similar to the method used in read.event() where a list of time step indexes is created for the files, and for each time step in 't' all files and time steps in the files are located using which(). As a consequence the parameters 'vesselutim' and 'schoolindt' were added.
	# Update: 2011-03-04 - Fixed a bug when using 'tvessel' to specify the output time points. The t value included in the name of the output files did not change to correspond to 'tvessel'. In this version 'vesselutim' is unchanged by 'tvessel' (as oposed to the previous version) and tvessel is used to specify the index when writing vesselutim to file.
	# Update: 2011-05-19 - Fixed a bug when static variables are present in the dynamic school files. The former variable 'schooltype' is replaced by the two variables 'schooltypeD' and 'schooltypeS', so that a file can be read as both dynamic and static school file.
	# Update: 2011-05-20 - Added the option 'rand.sel' for selecting randomly only a fraction of the school specified by the numeric 'rand.sel' (in the range (0, 1)).
	# Update: 2011-08-24 - Changed the handeling of the dumpfile from the old method of writing the dumpfile and deleting at the end of the function, to not writing, in the case that dumpfile has length 0 or 0 characters. This was forced through by the fact that echoIBM.calibrate() uses echoIBM.oneschool.oneping(), which wrote to the dumpfile unintendedly.
	# Update: 2012-02-08 - Added the option bptfile, which specifies the file name of the file to which the values 'sllf', 'rad1', 'rad2', 'pbp1' and 'pbp2' are to be written (NULL for not writing the file).
	# Update: 2012-03-28 - Removed the parameters 'lenkL', 'max.cells', 'max.radius', 'pres'.	
	# Update: 2013-09-01 - Added the parameter 'scls' used to reduce CPU time.
	# Update: 2013-11-01 - Added parallel processing using the option 'cores'.
	# Update: 2013-11-07 - Changed 'recycle' to possilby specify several time steps to recycle, a vector or a function taking the index numbers of the time steps as the first input.
	# Update: 2014-01-31 - Replaced 'recycle' by 'pingsSchool'.
	# Update: 2014-02-26 - Added the parameter 'areCompact', indicating compactly specified schools. Also changed according to the changes in echoIBM.generate_dynschool(), which now take the 'vesselutim' as imput.
	# Update: 2014-05-05 - Added 'vxIX' to compress the output.
	# Update: 2014-05-28 - Fixed bugs.
	# Update: 2014-05-28 - Added 'saveBrkt' with default FALSE, to avoid writing the mostly empty Barakat values, which only occur when there is wery few targets.
	# Last: 2015-03-17 - Removed conversion from dirl to dire.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	##### Functions: #####
	# Function for extracting a time step from a vector, and leave it untouched if of length 1:
	extract_timestep <- function(x, t){
		if(length(x)>1){
			x[t]
		}
		else{
			x
		}
	}
	
	# Function used to simulate and write to file one ping:
	simulateWriteOne <- function(thist, thist_vec, data, pingsSchool, pingfile, event, vesselfiles, origin, mode, dynschoolfiles, adds, dynschoolnames, staticschoolnames, dumpfile, tempdump, nchart, esnm, TVG.exp, compensated, calibrate, noise, parlist, rand.sel, max.memory, lt, tvessel, scls, discardOutside, fishReaction){
		
		# Store the process time to save the time used for the simulation:
		ptm <- proc.time()[3]
		thiswarning <- NULL
		thisdumpfile <- file.path(tempdump, paste(substr(basename(dumpfile), 1, nchar(basename(dumpfile))-4), "_", zeropad(thist, nchart), ".txt", sep=""))
	
		# Get the dynamic school file number 'thisfilenr' and the time step number of the file 'thisindtoffile' for the current time step:
		thisfilenr <- pingsSchool[[thist]][, 2]
		thisindtoffile <- pingsSchool[[thist]][, 3]
	
		# Simulate echo sounder observation for ping thist:
		vessel <- read.TSDs(vesselfiles, t=thist, var="all", header=FALSE, info=FALSE)
		if(length(vessel)==0){
			warning(paste("No vessel data present in the following vessel files for time step ", thist, ":\n", paste(vesselfiles, collapse="\n"), sep=""))
		}
		else if(!all(sapply(vessel, length)==1)){
			stop("The vessel file(s) may not be organized by time steps, causing more than one value to be read for the vessel information at each time step")
		}
		if(is.null(vessel$psxv) && is.null(vessel$psyv)){
			vessel0 <- read.TSDs(vesselfiles, t=origin, var=c("lonv", "latv"), header=FALSE)
			posv <- global2car(cbind(vessel$lonv, vessel$latv), vessel0)
			vessel$psxv <- posv[, 1]
			vessel$psyv <- posv[, 2]
		}
		# Add the vessel information to the data:
		data[names(vessel)] <- vessel
	
		### # 2017-10-16: Reading in beams data at each time step was a bad idea, since the beams data are modified using echoIBM.default.oneschool() in echoIBM.oneschool(). Thus we only accept one beam setting, i.e., one ping of beams data, and all the work on implemening an indp file will not be used for now, but may be used in the future? ###
		### # Read beams data for the current time step:
		### beams <- read.event(event, t=thist, var="beams")
		### data[names(beams)] <- beams
	
		# Write info to the current dumpfile:
		write(paste("\n\n\n\n\n########## TIME STEP ", thist, " (", paste(ftim.TSD(vessel, format="yyyy-mm-dd HH:MM:SS.FFF"), " GMT", sep=""), "): ##########", sep=""), thisdumpfile, append=FALSE)

		# Get the dynamic school variables, if mode is "active":
		# Add the seed of the randomness added to the school in echoIBM.generate_oneschool():
		data$seed <- parlist$seed[thist]
		data$discardOutside <- discardOutside
		data$rand.sel <- rand.sel
		
		dynschool <- list()
		if(strff("a", mode[1])){
			if(areCompact){
				# Apply the 'utim' information used as input to echoIBM.moveSchools() if 'recycle' was given as a numeric in echoIBM(), in which case 'pingsSchool' contains only the numerics 'recycle':
				dynschool <- echoIBM.generate_dynschool(data, t=thisindtoffile[1], vesselutim=vesselutim, adds=adds, dumpfile=thisdumpfile)
			}
			else if(length(thisfilenr)>0){
				if(length(thisfilenr)==1){
					dynschool <- read.TSD(dynschoolfiles[thisfilenr[1]], t=thisindtoffile[1], var="all", header=FALSE, dimension=FALSE)
				}
				else{
					tempdynschool <- list()
					for(i in seq_along(thisfilenr)){
						thisdata <- read.TSD(dynschoolfiles[thisfilenr[i]], t=thisindtoffile[i], var="all", header=FALSE, dimension=FALSE)
						tempdynschool[paste(names(thisdata), i, sep="")] <- thisdata
					}
					# Merge the data read in separate steps in the for loop:
					names(tempdynschool) <- substr(names(tempdynschool), 1, 4)
					namestempdynschool <- names(tempdynschool)
					for(i in seq_along(namestempdynschool)){
						dynschool[[namestempdynschool[i]]] <- unlist(tempdynschool[namestempdynschool==namestempdynschool[i]])
					}
				}
			}
		}
		gentime <- proc.time()[3]-ptm
		cat("Time used for generating schools:", gentime, "sec \n")
		
		# Add the data in 'adds' (overrides the data read in the school files):
		#dynschoolinadds <- intersect(names(adds), dynschoolnames)
		#dynschool[dynschoolinadds] <- adds[dynschoolinadds]
		
		# Removed since the info is written to the individual dump files:
		### Write info to the dumpfile:
		###if(length(dumpfile)>0 && nchar(dumpfile)>0){
		###	write(paste("\n\n\n\n\n########## TIME STEP ", thist, " (", paste(ftim.TSD(vessel, format="yyyy-mm-dd HH:MM:SS.FFF"), " GMT", sep=""), "): ##########", sep=""), dumpfile, append=TRUE)
		###	}
		
					## Discard fish outside of the radial and angular range of the sonar/echosounder:
					##discardOutside <- c(2, 3, 3)
					#temp <- echoIBM.fishInside(dynschool=dynschool, data=data, dumpfile=thisdumpfile, discardOutside=discardOutside, rand.sel)
					#dynschool <- temp$dynschool
					#data <- temp$data
		            #
					## Add randomness in fish positions and orientations:
					#dynschool$seed <- parlist$seed[thist]
					#dynschool <- echoIBM.addRandomness(dynschool, grsd_plHS=NULL)

		# Apply reaction of the fish to the vessel:
		thisFishReaction <- lapply(fishReaction, extract_timestep, t=thist)
		dynschool <- echoIBM.fishReaction(c(dynschool, vessel, thisFishReaction, data["psze"]))

		# Keep only the dynamic school variables:
		# Strip 'dynschool' of any static variables:
		dynschool <- dynschool[intersect(dynschoolnames, names(dynschool))]
		
		# If mode is "passive" the funciton echoIBM.oneping.oneschool.passive() is used:
		if(sum(unlist(lapply(dynschool, length)))==0){
			if(strff("a", mode[1])){
				thiswarning <- c(thiswarning, paste("School dynamics missing for time step ", thist, sep=""))
			}
			sv <- echoIBM.oneping.oneschool.passive(data, esnm=esnm, TVG.exp=TVG.exp, compensated=compensated, calibrate=calibrate, noise=noise, dumpfile=thisdumpfile, parlist=parlist)
		}
		else{
			# Simulate the i'th time step of the current school (the order of the input list c(dynschool, data, vessel) is important, because the dynamic school data may be present in files containing static variables, and in the order used the specifically read dynamic data has precedence):
			
			### # Extract a random selection of the targets using 'rand.sel':
			### Nl <- max(length(dynschool$psxf), length(dynschool$psyf), length(dynschool$pszf))
			### if(0<rand.sel[1] && rand.sel[1]<1){
			### 	if(length(rand.sel)>1){
			### 		set.seed(rand.sel[2])
			### 	}
			### 	affected.variables <- c(dynschoolnames, staticschoolnames)
			### 	#selection <- sample(c(TRUE, FALSE), Nl, TRUE, c(rand.sel[1], 1-rand.sel[1]))
			### 	selection <- sample(seq_len(Nl), round(Nl * rand.sel[1]))
			### 	for(j in seq_along(affected.variables)){
			### 		thisvar <- affected.variables[j]
			### 		if(length(data[[thisvar]])==Nl && !is.function(data[[thisvar]])){
			### 			data[[thisvar]] <- data[[thisvar]][selection]
			### 		}
			### 		if(length(dynschool[[thisvar]])==Nl && !is.function(dynschool[[thisvar]])){
			### 			dynschool[[thisvar]] <- dynschool[[thisvar]][selection]
			### 		}
			### 	}
			### }
			
			# Repeat data$epss and data$epsl if any of them are not NULL or a function:
			Nl <- max(length(dynschool$psxf), length(dynschool$psyf), length(dynschool$pszf))
			if(length(data$epss)>0 && !is.function(data$epss)){
				data$epss <- rep(data$epss, length.out=Nl)
			}
			if(length(data$epsl)>0 && !is.function(data$epsl)){
				data$epsl <- rep(data$epsl, length.out=Nl)
			}
	
			# Simulate the current time step (adding dynamic vessel data from 'adds'):
			#data <- c(list(scls=scls), dynschool, lapply(adds[intersect(names(adds), dynvesselnames)], function(x) x[thist]), data[setdiff(names(data), c(names(dynschool), names(vessel)))], vessel)
			thisdata <- c(dynschool, lapply(adds[intersect(names(adds), dynvesselnames)], function(x) x[thist]), data[setdiff(names(data), c(names(dynschool), names(vessel)))], vessel)
			
			# If scls is given in the school data and in the input parameter 'scls', multiply these:
			if(scls!=1 && length(thisdata$scls)){
				thisdata$scls <- thisdata$scls * scls
			}
			
			#---# # If the 'indp' is given, which links each time step to presets in the beams data, select the appropriate preset. Otherwise extract the current time step, if any:
			#---# beamsnames <- labl.TSD("rb")
			#---# if(length(dim(data$freq))>1){
			#---# 	areBeamsVar <- intersect(names(data), beamsnames)
			#---# 	data[areBeamsVar] <- extractTimeStep(data[areBeamsVar], if(length(data$indp)) data$indp[thist] else thist)
			#---# }
		
			sv <- echoIBM.oneping.oneschool(thisdata, esnm=esnm, TVG.exp=TVG.exp, compensated=compensated, calibrate=calibrate, noise=noise, max.memory=max.memory, dumpfile=thisdumpfile, ask=FALSE, parlist=parlist, msg=msg)
		}
		simtime <- proc.time()[3]-ptm - gentime
		cat("Time used for simulation:", simtime, "sec \n")
		
		# Compress the acoustic data:
		vbsc <- c(sv$sv)
		vxIX <- which(vbsc>0)
		if(length(vxIX)>=length(vbsc)/2){
			vxIX <- NULL
		}
		else{
			vbsc <- vbsc[vxIX]
		}
		# Write the ping:
		###if(lt==length(tvessel)){
		###	x_time <- list(indt=as.double(tvessel[t==thist]), utim=vesselutim[tvessel[t==thist]])
		###}
		###else{
		x_time <- list(indt=as.double(thist), utim=vesselutim[thist])
		###}
		if(saveBrkt){
			x <- c( x_time, list(numb=data$numb, lenb=data$lenb, vbsc=vbsc, vxIX=vxIX, Brkt=sv$Brkt, nsig=sv$nsig), vessel[c("psxv", "psyv", "pszv", "rtxv", "rtyv", "rtzv", "ispv")] )
		}
		else{
			x <- c( x_time, list(numb=data$numb, lenb=data$lenb, vbsc=vbsc, vxIX=vxIX), vessel[c("psxv", "psyv", "pszv", "rtxv", "rtyv", "rtzv", "ispv")] )
		}
		# Assure that no duplicated elements names of 'x' are present when writing to file (Favouring the first of duplicatly named elements):
		x <- x[!duplicated(names(x))]
		
		# Write the data to file:
		write.TSD(con=pingfile, x=x, numt=1, ow=ow, keep.float=TRUE, append=thist!=thist_vec[1], reserve=length(thist_vec))
		
		# Print and write the time used for the current simulation:
		simtime <- proc.time()[3]-ptm
		cat("Time used in total:", simtime, "sec \n------------------------------\n")
		
		if(length(thisdumpfile)>0 && nchar(thisdumpfile)>0){
			write(c(thist, simtime), timedumpfile, append=TRUE, sep="\t")
		}
		# Return the warnings:
		return(thiswarning)
	}
	
	
	
	simulateWrite <- function(thist_vec, data, pingsSchool, event, vesselfiles, origin, mode, dynschoolfiles, adds, dynschoolnames, staticschoolnames, dumpfile, tempdump, nchart, esnm, TVG.exp, compensated, calibrate, noise, parlist, rand.sel, max.memory, lt, tvessel, scls, discardOutside, fishReaction){
		# Define the name of the file to which the simulated data are written (containing the first time step index):
		pingfile <- file.path(pingsdir, paste(pingsname, "_", esnm, "_T", zeropad(thist_vec[1], nchar(lt)), ".pings", sep=""))
		# Run through the time steps:
		outputWarn <- lapply(thist_vec, simulateWriteOne, thist_vec=thist_vec, data=data, pingsSchool=pingsSchool, pingfile=pingfile, event=event, vesselfiles=vesselfiles, origin=origin, mode=mode, dynschoolfiles=dynschoolfiles, adds=adds, dynschoolnames=dynschoolnames, staticschoolnames=staticschoolnames, dumpfile=dumpfile, tempdump=tempdump, nchart=nchart, esnm=esnm, TVG.exp=TVG.exp, compensated=compensated, calibrate=calibrate, noise=noise, parlist=parlist, rand.sel=rand.sel, max.memory=max.memory, lt=lt, tvessel=tvessel, scls=scls, discardOutside=discardOutside, fishReaction=fishReaction)
		# Return the warnings:
		return(pingfile)		
	}
	
	
	# If this function is the top level function, create the warnings object:
	echoIBM.warnings_warninglist <- NULL
	
	# 'rand.sel' may be given higher than 1, implying that the targets remaining after a random selection are scaled correspondingly:
	if(rand.sel[1]>1){
		scls <- rand.sel[1]
		rand.sel[1] <- 1/rand.sel[1]
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist, "'scls' set to 'rand.sel' and 'rand.sel' set to 1/rand.sel")
		#warning("'scls' set to 'rand.sel' and 'rand.sel' set to 1/rand.sel")
		}
	
	# The dynamic variable names of the vessel (elements in 'adds' named by one of these names will be subsetted for each ping):
	dynvesselnames <- labl.TSD(c("v", "t"), list.out=FALSE)
	# The dynamic variable names of the school, and legal time variable names:
	dynschoolnames <- labl.TSD("ds")
	# The compactly specified dynamic variables of the school:
	compactschoolnames <- labl.TSD("cs")
	# The static variable names of the school:
	staticschoolnames <- labl.TSD("ss")
	
	
	filegroups <- echoIBM.getFileTypes(files, ext=c("beams", "ctd", "school", "vessel"), labl=list(beams=c("cali", "rad1"), school=c(dynschoolnames, staticschoolnames)))
	
	### # Separate 'files' into the apropriate file types:
	### ext = tools::file_ext(files)
	### #ext = lapply(strsplit(files, ".", fixed=TRUE), function(x) tail(x, 1))
	### beamsfiles = files[which(ext=="beams")]
	### ctdfiles = files[which(ext=="ctd")]
	### vesselfiles = files[which(ext=="vessel")]
	### #schoolfiles = files[which(ext=="school")]
	### schoolfiles = files[echoIBM.is.schoolfiles(files, c(dynschoolnames, staticschoolnames))]
	
	# Read the .beams files:
	#beams <- read.TSDs(filegroups$beams, t="all", var="all", drop=FALSE)
	### # 2017-10-16: Reading in beams data at each time step was a bad idea, since the beams data are modified using echoIBM.default.oneschool() in echoIBM.oneschool(). Thus we only accept one beam setting, i.e., one ping of beams data, and all the work on implemening an indp file will not be used for now, but may be used in the future? ###
	### beams0 <- read.event(event, t="all", var="beams")
	
	beams0 <- read.TSDs(filegroups$beams, t=1, var="all")
	
	if(is.null(esnm)){
		esnm <- head(beams0$esnm, 1)
	}
	
	# Read the .ctd files:
	ctd <- read.TSDs(filegroups$ctd, t="all", var="all")
	# Transform airpressure to Pascal (is airpressure is given and is less than 100 (implying decibar)):
	if(!is.null(ctd$hpr0) && ctd$hpr0<100){
		ctd$hpr0 <- ctd$hpr0*10000
		}
		
	# Detect dynamic and static school files:
	schooltype <- echoIBM.getSchoolfileType(filegroups$school, dynschoolnames, c(compactschoolnames,staticschoolnames))
	dynschoolfiles <- filegroups$school[schooltype$schooltypeD==1]
	staticschoolfiles <- filegroups$school[schooltype$schooltypeS==1]
	if(length(dynschoolfiles)==0 && length(staticschoolfiles)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist, "No school files found in the simulation files")
		}
	
	# Read static school variables. Here the compactly saved dynamic school data are also read:
	staticschool <- read.TSDs(staticschoolfiles, var="all")
		
	# Strip 'staticschool' of any dynamic variables:
	staticschool <- staticschool[intersect(c(compactschoolnames,staticschoolnames), names(staticschool))]
	
	# Merge the empirical beam patterns, and add these to the start of the list while removing alreaddy existing empirical beam pattern information:
	names_ebpf <- c("graf", "gref", "grsf", "ebpf")
	staticschool <- c(merge_ebpf(staticschoolfiles), staticschool[setdiff(names(staticschool), names_ebpf)])
		
	# Check whether any vessel dynamics supplied in 'adds' have the required length:
	dynvesselinadds <- intersect(names(adds), dynvesselnames)
	addstooshort <- sapply(adds[dynvesselinadds], function(x) length(x)>0 & length(x)<max(t))
	if(any(addstooshort)){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist, paste("Additional vessel dynamics", paste(dynvesselinadds[addstooshort], collapse=""), "supplied in 'adds' are to short and are recycled to max(t)", sep=""))
		adds[dynvesselinadds[addstooshort]] <- lapply(adds[dynvesselinadds[addstooshort]], function(x) rep(x, length.out=max(t)))
		}
	
		
	########## Execution and output ##########
	# Simulate and write:
	if(length(t)>0){
		totallength <- 0
		newfile <- TRUE
		
		# Set defaults and calculate stuff that need only be caculated once (sonar specific and static school parameters). Also report errors and warnings:
		cat("Treating beam patterns and related parameters...\n")
		data <- c(adds[setdiff(names(adds), dynvesselnames)], beams0, ctd, staticschool)
		data <- echoIBM.default.oneschool(data, esnm=esnm, compensated=compensated, method=method, dumpfile=dumpfile, bptfile=bptfile, max.radius=max.radius)
		
		# Dump information before going through the time steps:
		if(length(dumpfile)>0 && nchar(dumpfile)>0){
			# Information about the acoustical instrument:
			echoIBM.dump_summary(data[labl.TSD(c("b", "rb"), list.out=FALSE)], dumpfile, type="beams", append=TRUE)
			# CTD-information:		
			echoIBM.dump_summary(data[labl.TSD("ctd")], dumpfile, type="ctd", append=TRUE)
			# Information about static variables of the targets:
			echoIBM.dump_summary(data[labl.TSD("ss")], dumpfile, type="staticschool", append=TRUE)
			
			# Print the warnings to the dumpfile:
			write("\n\n\n##### WARNINGS in \"echoIBM.oneschool\" PRIOR TO SIMULATION: #####", dumpfile, append=TRUE)
			if(length(echoIBM.warnings_warninglist)>0){
				for(i in seq_along(echoIBM.warnings_warninglist)){
					write(paste(i, ": ", echoIBM.warnings_warninglist[i], sep=""), dumpfile, append=TRUE)
				}
			}
			else{
				write("none", dumpfile, append=TRUE)
			}
			echoIBM.warnings_warninglist <- NULL
		}
		
		# Create a directory of temporary dump files:
		tempdump <- file.path(dirname(dumpfile), "tempDumps")
		if(!file.exists(tempdump)){
			dir.create(tempdump)
		}
			
		
		# Parallel processing using the pblapply() function in the pbapply package:
		parallel <- cores > 1
		maxt <- max(t)
		# Parallel processing using the pblapply() function in the pbapply package:
		
		# Detect the number of cores and use the minimum of this and the number of requested cores:	
		cores <- min(cores, length(t), detectCores())
		# Split 't' into a list of blocks of time steps:
		if(length(filesize)>0){
			nBytes <- 4
			fact_total_vs_vbsc <- 1.1
			t <- splitSeqIntoBlocks(t=t, size1=beams0$numb[1] * max(beams0$lenb) * nBytes * fact_total_vs_vbsc, size=filesize, blocks=cores)
			cat("Time steps grouped as follows:\n")
			cat(paste(seq_along(t), sapply(t, prettyIntegers), sep=": "), sep="\n")
		}
		
		if(parallel){
			cat("Parallel simulations on", cores, "cores:\n")
			cores <- makeCluster(cores)
		}
		# Progress bar parallel processing (if cores>1):
		else{
			cat("Simulating:\n")
		}
		out <- pblapply(t, simulateWrite, data=data, pingsSchool=pingsSchool, event=event, vesselfiles=filegroups$vessel, origin=origin, mode=mode, dynschoolfiles=dynschoolfiles, adds=adds, dynschoolnames=dynschoolnames, staticschoolnames=staticschoolnames, dumpfile=dumpfile, tempdump=tempdump, nchart=max(nchar(unlist(t))), esnm=esnm, TVG.exp=TVG.exp, compensated=compensated, calibrate=calibrate, noise=noise, parlist=parlist, rand.sel=rand.sel, max.memory=max.memory, lt=maxt, tvessel=tvessel, scls=scls, discardOutside=discardOutside, fishReaction=fishReaction, cl=cores)	
		if(parallel){
			stopCluster(cores)
		}
		
		# Merge the dump files:
		# Use a maximum file size:
		dumpfiles <- c(dumpfile, list.files(tempdump, full.names=TRUE))
		mergeTextFiles(dumpfiles, con=dumpfile, maxsize=dumpsize)
		unlink(tempdump, recursive=TRUE)
					
		# Print the warnings to the dumpfile:
		if(length(dumpfile)>0 && nchar(dumpfile)>0){
			write("\n\n\n##### WARNINGS in \"echoIBM.oneschool\" AFTER SIMULATION: #####", dumpfile, append=TRUE)
			if(length(echoIBM.warnings_warninglist)>0){
				for(i in seq_along(echoIBM.warnings_warninglist)){
					write(paste(i, ": ", echoIBM.warnings_warninglist[i], sep=""), dumpfile, append=TRUE)
				}
			}
			else{
				write("none", dumpfile, append=TRUE)
			}
		}
		return(out)
	}
	else{
		return(NULL)
	}
}
