#*********************************************
#*********************************************
#' Simulates echo sounder observations based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school. The variables of the school, vessel, acoustic instrument and sea are given as files in the TSD format located in the directory specified by 'event'. More than one school may be given located in separate directories in 'event'. The simulated echogram will be stored as a TSD file in 'event'. The function includes an option to plot the vessel path relative to the centers of mass of the school. Schools can be given compactly in order to save space, in which case the individual fish information is generated on the fly (see "Example of compact schools in echoIBM.R").
#'
#' @param event  is a vector of strings, where the first is the path to the directory of the simulation files, and any additional stings are directories or files containing for example school files located elsewhere than the main directory given by the first string in 'event'. Additional schools must be given as directories in 'event', not as paths to individual files, as this implies that these files "belong" to any existing school located in event[1]. Files containing beam configuration, noise values, calibration values, ctd-values and vessel dynamics need to be given in the main directory! Also the directories need to be exact, and not directories above the actual directory, like in read.event(). School information may be given also ...???
#' @param t  is a vector of the numbers of the pings to be simulated, as listed from 1 to the number of time steps of the simulated school. If t=="all", all time steps are simulated.
#' @param adds  is an optional list of variables overriding the variables located in the 'event' directory.
#' @param rph  is a matrix of two columns of length 3 representing the mean (column 1) and the standard deviation (column 2) of the roll values (rtxv), pitch values (rtyv) and heave values (przv) of the vessel, used in gaussian simulation of the roll, pitch and heave values.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#' @param compensated  specifies which rotations are compensated for by the echo sounder:
#' @param filesize  is the maximum size of the merged files.
#' @param calibrate  is FALSE if calibration data are to be discarded if present.
#' @param noise  See echoIBM.add.noise().
#' @param mode  is one of "active" and "passive".
#' @param tvessel  has two different interpretations depending on the option 'path': In the case that path==TRUE 'tvessel' is the UNIX time values to be written to the .vessel-file. When path==FALSE 'tvessel' is a vector of the index numbers of the time steps to assign to the simulated pings. If tvessel==NULL the school time points are used.
#' @param scan.volume  is TRUE if the volume selected when path=TRUE is to be scanned to see if the school is inside the volume at each ping.
#' @param margin  is a vector of the margins on either side of the span of 'x' (recycled if not of length 4).
#' @param smooth  has one possible value "spline", smooths the x-values and the y-values of the generated vessel path separately using the spline function.
#' @param max.memory  is the maximum amount of memory to allow for the function to occupy. The level of for loops is chosen according to this value, i.e. if the method using only two loops, one for the radial distances and one for the unique frequencies, demands more memory than 'max.memory', the method usint three for loops is chosen.
#' @param ow  is TRUE if the user wish to overwrite existing file.
#' @param origin  is a vector of two elements representing the origin of the global coordinate system (G), or the numbering index of the ping regarded as the origin of (G) (ignoring heave so that the x-y-plane of (G) is on the surface of the sea).
#' @param recycle  can be given in a number of ways. If more than one school is to be simulated, 'recycle' must either be a list, or is converted to a list, of length equal to the number of schools / folders including school files. Each list element can be (1) a function selecting the appriopriate time steps, as indexed for each school (say if school nr. 2 has the time step indices 4, 5, 7 with respect to the time steps of the vessel, and the function is to alternate between the second and third time step, the result is 5, 7, 5, 7, and so on). A second possibility (2) is that 'recycle' is given as TRUE, to indicate recycling the first time step. Also (3) a vector of time step indices is accepted, which for the example above would be 2:3. If schools are given compactly (see "Example of compact schools in echoIBM.R"), 'recycle' can be given as a single numeric to freeze the schools at a specific time step.
#' @param keep.temp  has two elements, where the first is TRUE if the temporary directory holding the noiseless data directory is to be kept, and where the second is TRUE if the noise-added data directory is to be kept.
#' @param cores  is an integer specifying the number of cores to run the simulations over in parallel (should be lower than the number of cores in the computer).
#' @param rand.sel  is a numeric specifying a random selection of the school to use in the simulations. If rand.sel>1, 1/rand.sel of the targets are selected, and the selected targets are scaled by 'rand.sel' to keep the original backscatter from the school. If given as a vector of length 2, the second element is regarded as the seed for the random selection.
#' @param scls  is a factor by which the backscattering cross sectional area of the targets are scaled.
#' @param method  is "closest" if the beam pattern value of the closest grid point is to be selected, and "linear" if linear interpolation should be used to extract the beam pattern value (time demanding and only available for 2D grids).
#' @param ask  is TRUE if the used should be asked to for approval if the memory of the least memory demanding calculation method of the individual radial sampling intervals exceed the memory limit 'max.memory'.
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#' @param bptfile  is the name of the file to which 'sllf', 'rad1', 'rad2', 'pbp1' and 'pbp2' is written (use NULL for no writing).
#' @param path  is TRUE to allow the user to select vessel positions interactively.
#' @param pathnr  is the number of the vessel path drawn, to be used in the name of the ".vessel" file. If a new ".vessel" file is to be written, 'pathnr' must be different than the pathnr of existing files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR echoIBM.getSchoolfileType
#' @importFrom TSD labl.TSD prettyIntegers read.TSD read.TSDs s2dhms strff utim.TSD
#'
#' @export
#' @rdname echoIBM
#'
echoIBM <- function(event, t=1, adds=NULL, rph=NULL, esnm=NULL, TVG.exp=2, compensated=c("pitch","roll"), filesize=3e8, calibrate=TRUE, noise=c("bg","cex"), mode=c("active","passive"), tvessel=NULL, scan.volume=TRUE, margin=500, smooth="spline", max.memory=1e9, ow=TRUE, origin=1, recycle=FALSE, keep.temp=c(TRUE,FALSE), dumpsize=10e6, cores=1, rand.sel=1, scls=1, method=c("closest","linear"), ask=FALSE, parlist=list(), bptfile=TRUE, max.radius=0.2, path=FALSE, pathnr=1, onlyMerge=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-09-30 - First version.
	# Update: 2010-03-16 - Added exponential noise implementation.
	# Update: 2010-10-19 - Generalized from the old funciton echoSim() to having more than one school as input, and to recycle the first  positions of schools at each time step (useful when simulating bottom echo).
	# Update: 2011-01-07 - Cleaned up bugs, and tested for multiple school input. Also separated the case of single school input, to avoid merging at the end of the simulation. Introduced set.path.vessel.cm(), separating the method of selecting vessel positions from the function body.
	# Update: 2011-01-11 - Sorted out bug when path==TRUE: 'numt' and 't' (if t="all") defined from school-files if path==TRUE and from vessel-files if path==FALSE. Also no longer need for the distinction between path==TRUE and path==FALSE when utilizing echoIBM.oneschool(). Also the new version of echoIBM.oneschool() (2011-01-11) is used.
	# Update: 2011-02-24 - Fixed some bugs (occuring when reading files of different 'indt', which were generated by default in read.TSD() when var="time"), and added dumping the histograms of the spherical transduces postitions for the fish.
	# Update: 2011-03-04 - Fixed a bug when using 'tvessel' to specify the output time points. The t value included in the name of the output files did not change to correspond to 'tvessel'. In this version 'vesselutim' is unchanged by 'tvessel' (as oposed to the previous version) and tvessel is used to specify the index when writing vesselutim to file.
	# Update: 2011-05-19 - Fixed a bug when static variables are present in the dynamic school files. The former variable 'schooltype' is replaced by the two variables 'schooltypeD' and 'schooltypeS', so that a file can be read as both dynamic and static school file.
	# Update: 2011-05-20 - Added the option 'rand.sel' for selecting randomly only a fraction of the school specified by the numeric 'rand.sel' (in the range (0,1)).
	# Update: 2011-08-05 - Added the option of inputing the maximum backscattering cross section sigma_0='sgbs' (according to the notation in Holmin et al 2012) which is equal to A_0/chi = 'acca'/'ssil'. This is advantageous because it discards the need for estimating the spherical surface integral ssil of the targets at a number of frequencies. The variables 'acca' and 'ssil' should only be used if 'acca' is required as an interpretable variable, such as the physical size of a bubble of air. The introduction of 'sgbs' also involves the possibly frequency dependent variable epsilon_sigma = sigma_bs/S^m = 'epss' which links the optimal backscattering cross section to the target size to the power 'm'. The power 'm' is called 'spow', and has default 'm'=2 (squared target size). In conclusion the old method of giving 'acca', 'ssil' and 'epsl' is complemented by the method of giving 'sgbs' and 'epss' and not caring about the spherical surface integral 'ssil'. Thus the new method should be faster and simpler with regards to the input.
	# Update: 2011-08-24 - Changed the handeling of the dumpfile from the old method of writing the dumpfile and deleting at the end of the function, to not writing, in the case that dumpfile has length 0 or 0 characters. This was forced through by the fact that echoIBM.calibrate() uses echoIBM.oneschool.oneping(), which wrote to the dumpfile unintendedly.
	# Update: 2011-11-28 - Changed to always write the data not added any noise to the tempfile, and as default keep these data.
	# Update: 2012-02-08 - Added the option bptfile, which specifies the file name of the file to which the values 'sllf', 'rad1', 'rad2', 'pbp1' and 'pbp2' are to be written (NULL for not writing the file).
	# Update: 2012-02-11 - Added the function getSchoolfileType() for retreiving the correct school-file type.	
	# Update: 2012-02-20 - Changed to use the function echoIBM.getSchoolfileType().	
	# Update: 2012-03-28 - Removed the parameters 'lenkL', 'max.cells', 'max.radius', 'pres'.	
	# Update: 2012-09-28 - Added the option of applying variable correlation throughout the volume, consistant with the correlation of the periodic noise in the MS70 sonar.
	# Update: 2013-10-01 - Added the parameter 'scls' used to reduce CPU time.
	# Update: 2013-11-01 - Added parallel processing using the option 'cores'.
	# Update: 2014-01-31 - Changed 'recycle' to possilby specify several time steps to recycle, given as a list of vectors if more than one school are to be simulated.
	# Update: 2014-01-31 - Added the actual school pings to the dumpfile.
	# Update: 2014-02-05 - Added the option of specifying compressed school information, causing echoIBM to generate fish positions along the way. This is done by specidying a set of required variables in a separate school file.
	# Update: 2014-02-05 - Changed to support recycling compactly specified schools.
	# Update: 2014-02-05 - Added report of time usage to the file 'timedumpfile'.
	# Last: 2015-02-22 - Changed to only read calibration data from calibration files.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Hard code the dumpfile names:
	dumpfile <- "dump.txt"
	timedumpfile <- "timedump"
	
	# Process time variables to be reported to the dumpfile:
	starttime <- Sys.time()
	startptime <- proc.time()[3]
	
	if(length(keep.temp) == 1){
		keep.temp <- c(keep.temp,FALSE)
		}
	
	# Define the warnings (because the usual system of R for specifying warnings only returns the warnings at the return of the top level function):
	echoIBM.warnings_warninglist <- NULL
	
	# Get the name of the event to be simulated ('event' may have more than one element, in which case the first string is considered to be the main directory of the simulation experimen and the other strings are files or directories located elswhere containing information needed in the simulation experiment. The name of the simulation experiment as exracted from the first string, as the third to last subdirectory (the name of the acoustic instrument simulated is the second to last and "tsd" is the last)):
	eventname <- rev(strsplit(event[1], "/", fixed=TRUE)[[1]])
	if(tolower(eventname[1]) != "tsd"){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,"The path to the event specified by 'event' does not end with \"tsd\". Naming of the .pings files may be unexpected")
		}
	eventname <- eventname[3]
	
	# Define the name of the file to which the radii of the trasducer elements, the beam pattern types, and the side lobe level factors are written:
	if(isTRUE(bptfile)){
		#bptfile <- paste(event[1],"/",eventname,"_beampattern.beams",sep="")
		bptfile <- file.path(event[1], paste0(eventname, "_beampattern.beams"))
		}
	else if(length(bptfile)>0 && nchar(bptfile)>0){
		#bptfile <- paste(event[1],"/",bptfile,".beams",sep="")
		bptfile <- file.path(event[1], paste0(bptfile, ".beams"))
		}
	
	# The dynamic variable names of the vessel (elements in 'adds' named by one of these names will be subsetted for each ping):
	dynvesselnames <- labl.TSD(c("v","t"), list.out=FALSE)
	# The dynamic variable names of the school, and legal time variable names:
	dynschoolnames <- labl.TSD("ds")
	# The compactly specified dynamic variables of the school:
	compactschoolnames <- labl.TSD("cs")
	# The static variable names of the school:
	staticschoolnames <- labl.TSD("ss")
	
	# Get the names of the files to read (all other string elements of 'event' than the first is additional directories or files, and thus added to the list extracted by list.files()):
	filesanddir <- c(list.files(event[1], full.names=TRUE), event[-1])
	
	# Register schools as separate directories:
	isdir <- file.info(filesanddir)$isdir
	schoolfiles <- list()
	schooldirs <- NULL
	underevents <- NULL
	for(i in which(isdir)){
		filesi <- list.files(filesanddir[i], full.names=TRUE, recursive=TRUE)
		#valid_schoolfiles <- echoIBM.is.schoolfiles(filesi, c(dynschoolnames,staticschoolnames))
		valid_schoolfiles <- echoIBM.is.schoolfiles(filesi, c(dynschoolnames, compactschoolnames))
		if(any(valid_schoolfiles)){
			schooldirs <- c(schooldirs,filesanddir[i])
			schoolfiles <- c(schoolfiles,list(filesi[valid_schoolfiles]))
			underevents <- c(underevents,basename(filesanddir[i]))
			}
		}
	# Get file extensions:
	#ext <- tools::file_ext(files)
	filegroups <- echoIBM.getFileTypes(filesanddir, ext=c("beams", "ctd", "vessel"), key=list(beams=c("bgns", "nrnp", "cali")))
	#ext <- lapply(strsplit(filesanddir[!isdir],".",fixed=TRUE),function(x) tail(x,1))
	
	# Add ".school"-files that are not located in separate directories (see at the beginning of 'schoolfiles'), but only if both static and dynamic school files are present:
	valid_schoolfiles <- echoIBM.is.schoolfiles(filesanddir[!isdir], c(dynschoolnames,staticschoolnames))
	if(any(valid_schoolfiles)){
		schoolfiles_main <- filesanddir[!isdir][valid_schoolfiles]
		# 'schooltype_mainD' and 'schooltype_mainS' denote dynamic school files and static school files in the main directory:
		schooltype <- echoIBM.getSchoolfileType(schoolfiles_main,dynschoolnames,staticschoolnames)
		schooltype_mainD <- schooltype$schooltypeD
		schooltype_mainS <- schooltype$schooltypeS
		schooltype_mainB <- schooltype$schooltypeB
		# The school files in the main directory are considered to be of a school to be simulated only if both static and dynamic school files are present, or if only one school file is present in the entire list of files.
		#if(sum(valid_schoolfiles)==1 || (any(c(schooltype_mainD==1, schooltype_mainB==1)) && any(schooltype_mainS==1))){
		if(any(c(schooltype_mainD == 1, schooltype_mainB == 1)) && any(schooltype_mainS == 1)){
			schoolfiles <- c(list(schoolfiles_main),schoolfiles)
			}
		# Else, add to the other schools:
		else if(any(schooltype_mainD == 1) || any(schooltype_mainS == 1)){
			schoolfiles <- lapply(schoolfiles,function(x) c(x,schoolfiles_main))
			}
		}
	
	
	# 'schooltype' is 0 for dynamic school files and 1 for static school files:
	schooltype <- schoolfiles
	# For loop through the school files:
	for(i in seq_along(schoolfiles)){
		schooltype[[i]] <- echoIBM.getSchoolfileType(schoolfiles[[i]],dynschoolnames,staticschoolnames)$schooltypeS
		}
	
	# Assure that 'recycle' is a list, and repeat to the apropriate length (the number of separate directories holding school files):
	if(length(schoolfiles)>1){
		# Transform to list:
		if(is.function(recycle)){
			recycle <- list(recycle)
			}
		else if(!is.list(recycle)){
			recycle <- as.list(recycle)
			}
		# Expand the list to the appropriate length:
		if(length(recycle)<length(schoolfiles)){
			recycle <- rep(recycle, length.out=max(length(recycle),length(schoolfiles)))
			}
		}
	# If only one school is given and 'recycle' is not a list, simply transform to a list:
	else if(!is.list(recycle)){
		recycle <- list(recycle)
		}
	
	# Read the .beams files (used when merging simulations of different school, and used for extracting the name 'esnm' of the acoustic instrument):
	# Read one by one, ignoring calibration data and storing the beam configuration of the calibration data to cross-reference with the other beam configuration information (no longer used):
	### beams <- list()
	### #cali <- list()
	### for(i in seq_along(filegroups$beams)){
	### 	thisbeams <- read.TSD(filegroups$beams[i],t="all",var="all",dimension=TRUE,header=FALSE)
	### 	if(!"cali" %in% names(thisbeams)){
	### 		beams <- c(beams,thisbeams)
	### 		}
	### 	#if("cali" %in% names(thisbeams)){
	### 	#	cali <- c(cali,thisbeams[c("cali","grde")])
	### 	#	}
	### 	#else{
	### 	#	beams <- c(beams,thisbeams)
	### 	#	}
	### 	}
	
	beams <- read.TSDs(filegroups$beams, t="all", var="all", clean=FALSE)
	# Add the beams information in 'adds':
	beamsinadds <- intersect(names(beams),names(adds))
	beams[beamsinadds] <- adds[beamsinadds]
	if(is.null(esnm)){
		esnm <- beams$esnm
		}
	
	# Read the .ctd files (used when merging simulations of different school):
	ctd <- read.TSDs(filegroups$ctd, t="all", var="all")
	ctdinadds <- intersect(names(ctd),names(adds))
	ctd[ctdinadds] <- adds[ctdinadds]
	
	# Read the time points of the vessel, and add these to 'uniqeutim' later:
	vesselutim <- utim.TSD(c(adds[intersect(names(adds),dynvesselnames)], read.TSDs(filegroups$vessel, t="all", var="time", header=FALSE)))
	if(is.list(vesselutim)){
		vesselutim <- vesselutim[[1]]
		}
	# The vessel files contains all time steps, and 'numt' is thus extracted from the vessel time points:
	indt <- seq_along(vesselutim)
	
	# Get the time points of the school files:
	areCompact <- logical(length(schoolfiles))
	schoolutim <- vector("list", length(schoolfiles))
	for(i in seq_along(schoolfiles)){
		thisutim <- suppressWarnings(read.TSDs(schoolfiles[[i]], var=c("time","psxS"), t="all", clean=FALSE))
		# Apply the vessel UNIX time points if the schools are given in the compact form:
		if(length(thisutim$psxS)>0){
			thisutim <- list(sort(vesselutim))
			areCompact[i] <- TRUE
			}
		else{
			thisutim <- utim.TSD(thisutim, keep.list=TRUE)
			}
		if(length(thisutim)>0){
			schoolutim[[i]] <- thisutim
			}
		}
	
	# If the simulations are done in passive mode, simply override any dynamic school variables by psxf=Inf, psyf=Inf, pszf=Inf, rtzf=0, sgbs=0 (WHY??????):
	if(strff("p",mode[1])){
		adds <- c(list(psxf=Inf, psyf=Inf, pszf=Inf, rtzf=0, sgbs=0),adds)
		schoolutim <- list(list(sort(vesselutim)))
		}
	suppressWarnings(uniqeschoolutim <- sort(unique(unlist(schoolutim))))
	
	
	# If path==TRUE, draw and write the vessel positions to file. Use school-files as the basis for extracting 'numt' when path==TRUE:
	if(path){
		# The number of time steps:
		numt <- length(uniqeschoolutim)
		# If t=="all", all time points are read:
		if(identical(t,"all")){
			t <- 1:numt
			}
		t <- t[t>=1 & t<=numt]
		
		# Get the dynamic school files to be used in when plotting the vessel positions. Else the dynamic school files are treated in echoIBM.oneschool:
		dynschoolfiles <- schoolfiles
		for(i in seq_along(schoolfiles)){
			dynschoolfiles[[i]] <- schoolfiles[[i]][schooltype[[i]] == 0]
		}
		# Set the path of the vessel:
		output <- set.path.vessel.cm(schooldirs=schooldirs, dynschoolfiles=dynschoolfiles, recycle=recycle, t=round(t), tvessel=tvessel, numt=numt, smooth=smooth, margin=margin, scan.volume=scan.volume, beams=beams, ctd=ctd, rph=rph, event=event, eventname=eventname, pathnr=pathnr)
		return(output)
	}
	# Use vessel-files as the basis for extracting 'numt' when path==FALSE:
	else{
		# Test whether all time points of the school files are present in the vessel files (ok for empty 'uniqeschoolutim'):
		if(!all(uniqeschoolutim %in% vesselutim)){
			stop(paste0("Time steps found in school files that are not present in the vessel files: ",paste(setdiff(uniqeschoolutim, vesselutim), collapse=", ")))
		}
		# The number of time steps:
		numt <- length(vesselutim)
		# If t=="all", all time points are read:
		if(identical(t,"all")){
			t <- 1:numt
		}
		t <- round(t[t>=1 & t<=numt])
		
		
		# Create a list for each school holding the school file numbers and the time step indices of the school files correponding to each time step in the vessel data:
		pingsSchool <- vector("list",length(schoolutim))
		for(i in seq_along(pingsSchool)){
			# Create a matrix of 4 columns:
			# (1) the utim information of the school files of the current school (after unlisting)
			# (2) the file number holding the time steps
			# (3) the time step indices in the school files
			# (4) the time step indices of the entire school
			filelengths <- sapply(schoolutim[[i]],length)
			utimFilenrIndt <- cbind(unlist(schoolutim[[i]]), rep(seq_along(schoolutim[[i]]),filelengths), sequence(filelengths))
			utimFilenrIndt <- cbind(utimFilenrIndt,match(utimFilenrIndt[,1],unique(utimFilenrIndt[,1])))
			# Split into a list for each unique school time step:
			utimFilenrIndt <- lapply(split(utimFilenrIndt,utimFilenrIndt[,1]), matrix,ncol=4)
			uniqueSchoolUtim <- sapply(utimFilenrIndt,"[",1,1)
			uniqueSchoolIndt <- sapply(utimFilenrIndt,"[",1,4)
			
			# Match with the vessel time points (creates as list of length equal to the number of vessel time steps):
			pingsSchool[[i]] <- utimFilenrIndt[match(vesselutim,uniqueSchoolUtim)]
			# At the requected time steps, aply the recycling:
			if(is.function(recycle[[i]])){
				cat("'recycle' given as a function for school nr. ", i, "\n", sep="")
				pingsSchool[[i]][t] <- utimFilenrIndt[recycle[[i]](seq_along(t))]
			}
			# If 'recycle' is given as a numeric, repeate the given time steps:
			else if(sum(recycle[[i]])>0){
				pingsSchool[[i]][t] <- utimFilenrIndt[rep(as.numeric(recycle[[i]]),length.out=length(t))]
				cat("The following time steps simulated for school nr. ",i,":\n",prettyIntegers(uniqueSchoolIndt, sep="-", collapse=", ", force=TRUE), "\n", sep="")
			}
		}
	}
	# Warning if 't' is empty:
	if(length(t) == 0){
		warning("'t' is empty. Possibly the range of the vessel time points is not covering the requested 't'")
		}
		
	
	########## Execution and output ##########
	# Get the path to the dumpfile, adding the time step numbers: 
	#tToAddToFiles <- prettyIntegers(if(length(tvessel)>0) tvessel else t,sep="-",collapse=", ",force=TRUE)
	tToAddToFiles <- paste(range(if(length(tvessel)>0) tvessel else t), sep="-", collapse=", ")
	dumpfile <- file.path(event[1], "dumpfiles", paste(dumpfile, tToAddToFiles, ".txt", sep="_"))
	# Create the directory of the dumpfiles:
	suppressWarnings(dir.create(dirname(dumpfile)))
	
	# Get the path to the timedumpfile, adding the time step numbers: 
	timedumpfile <- file.path(event[1],paste(timedumpfile, tToAddToFiles, ".txt", sep="_"))
	
	# Write to the dumpfile (1):
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		if(!file.exists(dirname(dumpfile))){
			dir.create(dirname(dumpfile))
			}
		write(paste("########## DUMP FROM THE SIMULATION EXPERIMENT BUILT FROM THE FILES LOCATED IN THE DIRECTORY (WITH SUBDIRECTORIES): ##########\n\n",event[1],
			ngettext(length(t),"\n\n##### SIMULATION DONE FOR TIME STEP: #####\n","\n\n##### SIMULATIONS DONE FOR TIME STEPS: #####\n"), 
			prettyIntegers(t, sep="-", collapse=",", force=TRUE),
			sep=""), 
			dumpfile, append=FALSE)
		
		for(i in seq_along(schoolfiles)){
			if(sum(recycle[[i]])>0){
				write(paste0(
					ngettext(length(t), paste0("\nUsing the following time step for school ", i, ":\n"), paste("\nUsing the following time steps for school ", i, ":\n")), 
					prettyIntegers(unlist(lapply(pingsSchool[[i]][t],function(xx) xx[,4])), sep="-", collapse=",", force=TRUE)), 
				dumpfile, append=TRUE)
				}
			}
		write(paste0("\n\n##### SIMULATION STARTED: #####\n", format(as.POSIXlt(starttime,tz="GMT"), usetz=TRUE)), dumpfile, append=TRUE)
		write("\n\n\n##### FILES USED IN THE SIMULATION: #####", dumpfile, append=TRUE)
		if(is.list(schoolfiles) && length(schoolfiles)>1){
			for(i in seq_along(schoolfiles)){
				write(paste0("SCHOOL NR. ", i, ":"), dumpfile, append=TRUE)
				write(unlist(schoolfiles[[i]]), dumpfile, append=TRUE)
				}
			}
		else{
			write(unlist(schoolfiles), dumpfile, append=TRUE)
			}
		write(filegroups$beams, dumpfile, append=TRUE)
		write(filegroups$vessel, dumpfile, append=TRUE)
		write(filegroups$ctd, dumpfile, append=TRUE)
		}
	
	### # 'rand.sel' may be given higher than 1, implying that the targets remaining after a random selection are scaled correspondingly:
	### if(rand.sel[1]>1){
	### 	scls <- rand.sel[1]
	### 	rand.sel[1] <- 1/rand.sel[1]
	### 	warning("'scls' set to 'rand.sel' and 'rand.sel' set to 1/rand.sel")
	### 	}
	# Write to the dumpfile (2):
	if(length(dumpfile)>0 && nchar(dumpfile)>0 && any(rand.sel[1] != 1, scls != 1)){
		# Report that a random selection of the targets is simulated, and that the backscatter from the targets is scaled:
		echoIBM.dump_summary(list(rand.sel=rand.sel, scls=scls), dumpfile, type="rand.sel_scale", append=TRUE)
		}
	
	# Default values for 'parlist':
	parlist <- echoIBM_rexp_defaults(noise=noise, indt=indt, data=beams, parlist=parlist)
	
	# Write to the dumpfile (3):
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# Report the noise applied to the simulations:
		implemented_noise <- intersect(noise,c("ms","acex","aex","cex","ex","nr","bg","pn","bk","phase"))
		if(length(implemented_noise) == 0){
			implemented_noise <- "none"
			}
		write(paste0("\n\n\n\n\n##### NOISE APPLIED TO THE SIMULATIONS: #####\n", paste(implemented_noise, collapse="\n")), dumpfile, append=TRUE)
		}	
	# Dump information before going through the time steps:
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# Information about the acoustical instrument:
		echoIBM.dump_summary(parlist, dumpfile, type="noise", append=TRUE)
		}
			
	# Create a temporary directory holding the separately generated simulations, to be merged at the end of the function. This directory is deleted if keep.temp[1]==FALSE:
	tempevent <- file.path(dirname(event[1]),"temp")
	dir.create(tempevent)
	
	# Write to the dumpfile:
	if(length(schoolfiles) == 0 || strff("p",mode[1])){
		# Set the 'schoolfiles' to a list of one empty element to enable the for loop below:
		schoolfiles <- vector("list",1)
		if(length(dumpfile)>0 && nchar(dumpfile)>0){
			write(paste0("\n\n\n##### SIMULATIONS DONE IN PASSIVE MODE #####\n\n\n"), dumpfile, append=TRUE)
			}
		}
	else{
		if(length(dumpfile)>0 && nchar(dumpfile)>0){
			write(paste0("\n\n\n##### SCHOOL NR. ", i, ": #####\n\n\n"), dumpfile, append=TRUE)
			}	
		}
	
	# Move through the schools:
	if(!onlyMerge){
		for(i in seq_along(schoolfiles)){
			# Get the list of files for the current school (merging general files for the simulation with the school specific files located in the individual school directories):
			files <- unique(c(filesanddir[!isdir],unlist(schoolfiles[i])))
			# No noise, add afterwards:
			echoIBM.oneschool(files=files, t=t, tvessel=tvessel, vesselutim=vesselutim, pingsSchool=pingsSchool[[i]], areCompact=areCompact[i], adds=adds, pingsdir=tempevent, pingsname=if(length(schooldirs)>1) paste(eventname,underevents[i],sep="_") else eventname, esnm=esnm, TVG.exp=TVG.exp, compensated=compensated, filesize=filesize, calibrate=calibrate, noise="", mode=mode, max.memory=max.memory, ow=ow, origin=origin, dumpfile=dumpfile, dumpsize=dumpsize, timedumpfile=timedumpfile, rand.sel=rand.sel, scls=scls, method=method, ask=ask, parlist=parlist, bptfile=bptfile, max.radius=max.radius, cores=cores)
		}
	}

	# Sum up the simulations and add noise:
	cat("\n")
	TVGinout <- sum(TVG.exp)>0
	
	# Merge the pings and add noise and randomness:
	#echoIBM.merge(event=event[1], t=if(length(tvessel)==0) t else tvessel, beams=beams, ctd=ctd, TVG.in=TVGinout, TVG.out=TVGinout, TVG.exp=TVG.exp, filesize=filesize, noise=noise, ow=ow, parlist=parlist, cores=cores, keep.temp=keep.temp[2])
	echoIBM.add.noise.TSD(event=event[1], fileind=NULL, beams=beams, ctd=ctd, TVG.in=TVGinout, TVG.out=TVGinout, TVG.exp=TVG.exp, noise=noise, scls=scls, ow=ow, parlist=parlist, cores=cores)
		
		
	# Delete the temporary directory if required
	if(!keep.temp[1]){
		unlink(tempevent,TRUE)
		}
	
	# Write to the dumpfile (6):
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# Print any echoIBM.warnings:
		write("\n\n\n##### WARNINGS IN \"echoIBM\": #####", dumpfile, append=TRUE)
		if(length(echoIBM.warnings_warninglist)>0){
			for(i in seq_along(echoIBM.warnings_warninglist)){
				write(paste0(i, ": ", echoIBM.warnings_warninglist[i]), dumpfile, append=TRUE)
				}
			}
		else{
			write("none", dumpfile, append=TRUE)
			}
		
		# Print the end time and time used by the simulation:
		endtime <- Sys.time()
		write(paste0("\n\n\n##### SIMULATION ENDED: #####\n",format(as.POSIXlt(endtime,tz="GMT"), usetz=TRUE)), dumpfile, append=TRUE)
		ptime=s2dhms(proc.time()[3]-startptime,names=c("day","hrs","min","sec"),strip=TRUE)
		write("\n\n\n##### TIME USED: #####", dumpfile, append=TRUE)
		write(colnames(ptime), ncolumns=5, file=dumpfile, sep="\t", append=TRUE)
		write(ptime, file=dumpfile, sep="\t", append=TRUE)
		}
	##################################################
	##################################################
	}
