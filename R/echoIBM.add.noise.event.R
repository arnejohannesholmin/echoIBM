#*********************************************
#*********************************************
#' Merges time steps of acoustic data located in the neighbor directory "temp" of the "tsd"-directory given by 'event'.
#'
#' @param event  is the path to the directory of the simulated files. If the length of 'event' is 2, the second string is the path to the temp-directory holding the files to merge, while the first is the directory in which to put the merged files.
#' @param beams0  is the list of beams inputs as returned from read.TSD. For the treatment of noise at least one of the following variables are required: bgns, nrn0.
#' @param ctd  is the list of ctd inputs as returned from read.TSD.
#' @param TVG.in  is TRUE if Time Varied Gain compensation is already applied to the input.
#' @param TVG.out  is TRUE if Time Varied Gain compensation should be applied to the output.
#' @param TVG.exp  is the exponent in the TVG amplification (usually 2 for Sv and 4 for TS).
#' @param noisetypes  See echoIBM.add.noise().
#' @param ow  is TRUE if the user wish to overwrite existing file(s).
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom SimradRaw apply.TVG
#' @importFrom sonR UNIX_time read.event_unzip_vbsc read.event
#' @importFrom TSD read.TSD read.TSDs s2dhms splitSeqIntoBlocks strff write.TSD zeropad
#' @importFrom utils tail
#'
#' @export
#' @rdname echoIBM
#'
echoIBM.add.noise.event <- function(event, pingsfiles=NULL, fileind=NULL, beams0=NULL, ctd=NULL, TVG.in=TRUE, TVG.out=TRUE, TVG.exp=2, noisetypes=c("nr","bg","cex"), ow=TRUE, parlist=list(), cores=1){
	
	########## Preparation ##########
	# Hard code the dumpfile name of the merging:
	dumpfile <- "dump_merge.txt"
	
	### Functions ###
	# (1) Function for reading simulated noiseless and deterministic school data and summing them up:
	readAddNoiseWrite <- function(fileind, pingsfiles, indtlist, beams0, ctd, noise, noisetypes, TVG.exp, parlist){
		#thisdata <- read.TSD(pingsfiles[fileind], t="all", drop.out=FALSE)
		#thisdata <- read.event_unzip_vbsc(thisdata, pad=TRUE, split=TRUE, filename=pingsfiles[fileind], t="all", fill=0)
		#thisdata[["vxIX"]] <- NULL
		thisPingsFiles <- sapply(pingsfiles, "[[", fileind)
		#thisdata <- read.TSDs(pingsfiles[fileind], t="all", drop.out=FALSE)
		thisdata <- read.TSDs(thisPingsFiles, t="all", drop.out=FALSE, clean = NULL)
		
		# Expand the pings:
		#thisdata <- read.event_unzip_vbsc(thisdata, pad=TRUE, split=TRUE, filename=pingsfiles[fileind], t="all", fill=0)
		thisdata <- lapply(thisdata, read.event_unzip_vbsc, pad=TRUE, split=TRUE, filename=pingsfiles[[1]][fileind], t="all", fill=0)
		# Sum the data:
		vbscList <- lapply(thisdata, "[[", "vbsc")
		vbsc <- do.call("+", vbscList)
		
		thisdata <- thisdata[[1]]
		thisdata$vbsc <- vbsc
		thisdata$vxIX <- NULL
		
		# Remove TVG if TVG compensation is applied to the input data:
		if(TVG.in){
			thisdata$vbsc <- SimradRaw::apply.TVG(thisdata$vbsc, beams=c(beams0, ctd), rm=TRUE, TVG.exp=TVG.exp)
		}
		
		### Add noise to the summed data: ###
		if(length(dim(thisdata$vbsc))==2){
			dim(thisdata$vbsc) <- c(dim(thisdata$vbsc), 1)
		}
		# Run through the time steps:
		for(i in seq_len(tail(dim(thisdata$vbsc), 1))){
			### # 2017-10-16: Reading in beams data at each time step was a bad idea, since the beams data are modified using echoIBM.default.oneschool() in echoIBM.oneschool(). Thus we only accept one beam setting, i.e., one ping of beams data, and all the work on implemening an indp file will not be used for now, but may be used in the future? ###
			### beams <- read.event(event, t=indtlist[[fileind]][i], var="beams")
			### thisdata$vbsc[,,i] <- echoIBM.add.noise(sv=thisdata$vbsc[,,i], noisetypes=noisetypes, data=c(parlist, noise, beams), indt=indtlist[[fileind]][i], parlist=parlist)
			thisdata$vbsc[,,i] <- echoIBM.add.noise(sv=thisdata$vbsc[,,i], noisetypes=noisetypes, data=c(parlist, noise, beams0), indt=indtlist[[fileind]][i], parlist=parlist)
		}
		
		# Add TVG if required:
		if(TVG.out){
			thisdata$vbsc <- SimradRaw::apply.TVG(thisdata$vbsc, beams=c(beams0, ctd), TVG.exp=TVG.exp)
		}
		
		# Write the data to file:
		pingfileWithNoise <- file.path(dirname(dirname(dirname(pingsfiles[[1]][fileind]))), "tsd", basename(pingsfiles[[1]][fileind]))
		numt <- tail(dim(thisdata$vbsc), 1)
		bytes <- TSD::write.TSD(con=pingfileWithNoise, x=thisdata, numt=numt, ow=ow, keep.null=TRUE)
	}	
	### End of functions ###
	
	
	# Process time variables to be reported to the dumpfile:
	starttime <- Sys.time()
	startptime <- proc.time()[3]
	cat("MERGING PINGS-FILES:\n")
	
	# Get the names of the .beams and .ctd files to read if 'beams' and/or 'ctd' is missing:
	cat("Listing files of the event...\n")
	filesanddir <- list.files(event[1],full.names=TRUE)
	isdir <- file.info(filesanddir)$isdir
	filesanddir <- filesanddir[!isdir]
	filegroups <- echoIBM.getFileTypes(filesanddir)
	
	# Read the .beams and noise files:
	cat("Reading beams and noise files...\n")
	if(length(beams0)==0){
		beams0 <- read.event(event, t="all", var="beams")
	}
	
	# Read all the noise data in one go. other = TRUE is needed to read the noise files:
	noise <- read.event(event, t="all", var=labl.TSD("noise"), other=TRUE)
		
	# Error if pulselength 'sint' and/or absorption factor 'absr' is missing:
	if(any(is.null(beams0$sint), is.null(beams0$absr))){
		stop(paste0("Elements \"sint\" and/or \"absr\" missing in the event \"", event[1], "\""))
	}	
	# Read the .ctd files:
	cat("Reading ctd files...\n")
	if(length(ctd)==0){
		ctd <- read.TSDs(filegroups$ctd, t="all", var="all")
	}
	# Add default speed of sound if no ctd-data are present:
	if(length(ctd)==0){
		ctd <- list(asps=1500)
	}	
	
	
	########## Execution and output ##########
	### Sum up the simulations and add noise: ###
	# Get the names of the .pings files located in the "temp"-directoy:
	tempevent <- file.path(dirname(event[1]), "temp")
	if(length(pingsfiles)==0){
		cat("Listing temporary files...\n")
		files <- list.files(tempevent, full.names=TRUE, recursive=TRUE)
		pingsfiles <- echoIBM.getFileTypes(files, ext=c("pings"))$pings
	}
	if(length(pingsfiles)==0){
		warning(paste0("There are no pings-files in the directory of temp-files ", tempevent))
	}
	#if(length(fileind)){
	#	pingsfiles <- pingsfiles[fileind]
	#}
	
	
	# First a quick test to see if there are an equal number of files as time steps, and file names correspond to a continuous sequence of time steps:
	#utim <- read.TSDs(pingsfiles, var="utim", t="all", clean=FALSE)
	utim <- read.TSDs(pingsfiles[[1]], var="utim", t="all", clean=FALSE)
	utim <- utim[names(utim)=="utim"]
	uniqeutim <- unique(unlist(utim))
	# Get the time step indexes of the files located in the temp-directory:
	indtlist <- lapply(utim, match, uniqeutim)
	indt <- unlist(indtlist, use.names=FALSE)
	
	# The default dumpfile name is paste0("dump_merge_",date(),".txt"): 
	if(!file.exists(dirname(dumpfile))){
		dir.create(dirname(dumpfile))
	}
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		write(paste0(
			"########## DUMP FROM MERGING AND ADDING NOISE TO SIMULATION FILES LOCATED IN THE DIRECTORY (WITH SUBDIRECTORIES): ##########\n\n", 
			tempevent, 
			"\n\n##### MERGING STARTED: #####\n",
			format(as.POSIXlt(starttime,tz="GMT"),usetz=TRUE)
			), 
			dumpfile)
		write("\n\n\n##### FILES MERGED: #####", dumpfile, append=TRUE)
		write(unlist(pingsfiles), dumpfile, append=TRUE)
	}
	
	# Write to the dumpfile (3):
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		
		# Report the noise applied to the simulations:
		implemented_noise <- intersect(gsub("_fun", "", noisetypes),c("ms","acex","aex","cex","ex","nr","bg","pn","bk","phase"))
		if(length(implemented_noise)==0){
			implemented_noise <- "none"
		}
		write(paste0("\n\n\n\n\n##### NOISE APPLIED TO THE SIMULATIONS: #####\n", paste(implemented_noise, collapse="\n")), dumpfile, append=TRUE)
	}	
	# Dump information before going through the time step:
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# Information about the noise:
		echoIBM.dump_summary(c(parlist, noise), dumpfile, type="noise", append=TRUE)
	}
	
	# Parallel processing using the pblapply() function in the pbapply package:
	parallel <- cores > 1
	# Parallel processing using the pblapply() function in the pbapply package:
	if(parallel){
		# Detect the number of cores and use the minimum of this and the number of requested cores:	
		cores = min(cores, length(pingsfiles[[1]]), detectCores())
		cores <- makeCluster(cores)
	}
	# Progress bar parallel processing (if cores>1):
	cat("Adding noise to simulated data:\n")
	pingsFileSequence <- seq_along(pingsfiles[[1]])
	out <- pblapply(pingsFileSequence, readAddNoiseWrite, pingsfiles=pingsfiles, indtlist=indtlist, beams0=beams0, ctd=ctd, noise=noise, noisetypes=noisetypes, TVG.exp=TVG.exp, parlist=parlist, cl=cores)
	
	if(parallel){
		stopCluster(cores)
	}
	
	# Print the end time and time used by the simulation:
	endtime <- Sys.time()
	write(paste0("\n\n\n##### MERGING ENDED: #####\n", format(as.POSIXlt(endtime, tz="GMT"), usetz=TRUE)), dumpfile, append=TRUE)
	ptime <- s2dhms(proc.time()[3]-startptime, names=c("day","hrs","min","sec"), strip=TRUE)
	write("\n\n\n##### TIME USED: #####", dumpfile, append=TRUE)
	write(colnames(ptime), ncolumns=5, file=dumpfile, sep="\t", append=TRUE)
	write(ptime, file=dumpfile, sep="\t", append=TRUE)
	##################################################
	##################################################
	}
