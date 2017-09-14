#*********************************************
#*********************************************
#' Merges time steps of acoustic data located in the neighbor directory "temp" of the "tsd"-directory given by 'event'.
#'
#' @param event  is the path to the directory of the simulated files. If the length of 'event' is 2, the second string is the path to the temp-directory holding the files to merge, while the first is the directory in which to put the merged files.
#' @param t  is a vector of the numbers of the pings to be merged, as given in the simulation event. If t=="all", all time steps are merged.
#' @param beams  is the list of beams inputs as returned from read.TSD. For the treatment of noise the following variables are required:
#' @param ctd  is the list of ctd inputs as returned from read.TSD.
#' @param TVG.in  is TRUE if Time Varied Gain compensation is already applied to the input.
#' @param TVG.out  is TRUE if Time Varied Gain compensation should be applied to the output.
#' @param TVG.exp  is the exponent in the TVG amplification (usually 2 for Sv and 4 for TS).
#' @param filesize  is the maximum size of the merged files.
#' @param calibrate.in  is TRUE if input data is calibrated.
#' @param calibrate.out  is TRUE if output data should be calibrated.
#' @param noise  See echoIBM.add.noise().
#' @param ow  is TRUE if the user wish to overwrite existing file(s).
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#' @param reserve  is FALSE if time steps should not be reserved, requiring that the lengths and dimensions of the variables in the files to merge match exactly. If variable length of the data for each time step, set reserve=TRUE!
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom SimradRaw apply.TVG
#' @importFrom sonR UNIX_time
#' @importFrom TSD read.TSD read.TSDs s2dhms splitSeqIntoBlocks strff write.TSD zeropad
#'
#' @export
#' @rdname echoIBM.merge
#'
echoIBM.merge <- function(event, t=1, beams=NULL, ctd=NULL, TVG.in=TRUE, TVG.out=TRUE, TVG.exp=2, filesize=3e8, noise=c("nr","bg","cex"), scls=1, ow=TRUE, parlist=list(), cores=1, reserve=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-09-30 - First version.
	# Update: 2010-03-16 - Added exponential noise implementation.
	# Update: 2011-11-29 - Replaced the method of reading variables that are not bbackscatter data, by a single line using read.TSDs().
	# Update: 2011-12-09 - Changed the treatment of time steps to consider 't' as the time step indexes of the total event, and not of the temp-directory only.
	# Update: 2011-12-10 - Added support for including the temp-directory in 'event'.
	# Update: 2012-06-28 - Added 'indt' as input to echoIBM.add.noise().
	# Update: 2013-11-01 - Added parallel processing using the option 'cores'.
	# Update: 2013-11-01 - Added 'scls'.
	# Update: 2014-05-05 - Added reading 'vxIX' to decompress the acoustic data.
	# Last: 2014-10-20 - Speeded up the function.
	########### DESCRIPTION: ###########
	# Merges time steps of acoustic data located in the neighbor directory "temp" of the "tsd"-directory given by 'event'.
	########## DEPENDENCIES: ###########
	# read.TSD(), read.TSDs(), zeros(), utim.TSD(), set.path(), echoIBM.oneping(), write.TSD()
	############ VARIABLES: ############
	# ---event--- is the path to the directory of the simulated files. If the length of 'event' is 2, the second string is the path to the temp-directory holding the files to merge, while the first is the directory in which to put the merged files.
	# ---t--- is a vector of the numbers of the pings to be merged, as given in the simulation event. If t=="all", all time steps are merged.
	# ---beams--- is the list of beams inputs as returned from read.TSD. For the treatment of noise the following variables are required:
	#		 For noise=="bg": 'lenb', 'numb', and 'bgns'
	#		 For noise=="pn": 'sint', 'lenb', 'pns1', 'pns2', 'pns3', 'acfq', 'harm', and 'bgns' or 'numb'
	#		 For noise=="hi": 'lenb', 'numb', 'freq', 'hins' and 'hini'
	#		 For noise=="nr","np","cp": 'nr0p' or 'cnpM', 'lenb', 'numb'
	#		 For noise=="na","ca": 'nr0a' or 'cnaM', 'lenb', 'numb'
	#		 For noise==""acex","aex","cex": 'rate','prob','l','buffer','esnm','rho','w','wC','luqf','numb','nuqf','Cind'	
	#		 For noise=="ms": 'L','N','P','esnm','olpn','w','scale'
	# ---ctd--- is the list of ctd inputs as returned from read.TSD.
	# ---TVG.in--- is TRUE if Time Varied Gain compensation is already applied to the input.
	# ---TVG.out--- is TRUE if Time Varied Gain compensation should be applied to the output.
	# ---TVG.exp--- is the exponent in the TVG amplification (usually 2 for Sv and 4 for TS).
	# ---filesize--- is the maximum size of the merged files.
	# ---calibrate.in--- is TRUE if input data is calibrated.
	# ---calibrate.out--- is TRUE if output data should be calibrated.
	# ---noise--- See echoIBM.add.noise().
	# ---ow--- is TRUE if the user wish to overwrite existing file(s).
	# ---parlist--- is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
	# ---reserve--- is FALSE if time steps should not be reserved, requiring that the lengths and dimensions of the variables in the files to merge match exactly. If variable length of the data for each time step, set reserve=TRUE!
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Hard code the dumpfile name of the merging:
	dumpfile="dump_merge.txt"
	
	### Functions ###
	# (1) Function for reading simulated noiseless and deterministic school data and summing them up:
	readAddNoise <- function(i, requested_indt, indt, index, pingsfiles, beams, ctd, noise, TVG.exp, parlist){
		# Read all the simulated noiseless and deterministic school data and sum them up:
		nschools <- 0
		thissv <- 0
		# Define the output data in case there is no data read in the for loop below:
		thisdata <- list()
		thisindex <- index[[requested_indt[i]]]
		if(length(thisindex)>0){
			for(j in seq_len(nrow(thisindex))){
				thisdata <- read.TSD(pingsfiles[thisindex[j,1]], t=thisindex[j,2])
				# Treat compressed data
				if(length(thisdata$vbsc)==0 && length(thisdata$vxIX)==0){
					thisdata$vbsc <- double(max(beams$lenb) * beams$numb[1])
					}
				else if(length(thisdata$vbsc)>0 && length(thisdata$vxIX)>0){
					temp <- double(max(beams$lenb) * beams$numb[1])
					temp[thisdata$vxIX] <- thisdata$vbsc
					thisdata$vbsc <- temp
					}
				# Remove TVG if TVG compensation is applied to the input data:
				if(TVG.in){
					thissv <- thissv + SimradRaw::apply.TVG(thisdata$vbsc, beams=c(beams,ctd), rm=TRUE, TVG.exp=TVG.exp)
					}
				else{
					thissv <- thissv + thisdata$vbsc
					}
				}
			}
		# Update the acoustic data by the sum over all files for the current time step:
		thisdata$vbsc <- thissv
		# Scale the signal:
		if(scls!=1){
			thisdata$vbsc <- thisdata$vbsc * scls
			}
		
		### Add noise to the summed data: ###
		# Set seed for the current time step:
		if(length(parlist$seed)>0){
			parlist$seed <- seed[i]
			}
		# Apply the Barakat PDF only if exactly one school is simulated:
		if(nschools==1){
			parlist$Brkt <- thisdata$Brkt
			parlist$nsig <- thisdata$nsig
			}
		
		# Add noise:
		thisdata$vbsc <- echoIBM.add.noise(sv=thisdata$vbsc, noise=noise, data=c(parlist,beams), indt=requested_indt[i], parlist=parlist)
		# Add TVG if required:
		if(TVG.out){
			thisdata$vbsc <- SimradRaw::apply.TVG(thisdata$vbsc, beams=c(beams,ctd), TVG.exp=TVG.exp)
			}
		
		# Return the data:
		thisdata[names(thisdata)!="vxIX"]
		}	
	
	# (2) Define a function called by parLapply(), which merges simulated data for each time step (if more than one school are included in the simulations), adds noise, and writes the resulting simulated ping to file. After applying parLapply() all the files are merged to larger files:
	readAddNoiseWrite <- function(i_vec, requested_indt, indt, index, pingsfiles, beams, ctd, noise, TVG.exp, parlist, tempevent, ow){
		for(i in i_vec){
			# Sum up data for the current time step and add noise:
			thisdata <- readAddNoise(i,requested_indt, indt, index, pingsfiles, beams, ctd, noise, TVG.exp, parlist)
			# Write ping i:
			if(i==i_vec[1]){
				pingfile <- file.path(tempevent, paste0(eventname, "_", beams$esnm, "_T", TSD::zeropad(requested_indt[i], nchar(max(requested_indt))), ".pings"))
				}
			# Write the data to file:
			bytes <- TSD::write.TSD(con=pingfile, x=thisdata, numt=1, ow=ow, keep.null=TRUE, append=i > i_vec[1])
			}
		}
	### End of functions ###
	
	
	# Process time variables to be reported to the dumpfile:
	starttime <- Sys.time()
	startptime <- proc.time()[3]
	cat("MERGING PINGS-FILES:\n")
	
	# The tempevent can be given as the second string of 'event':
	if(length(event)>1){
		eventname <- rev(strsplit(event[1], "/", fixed=TRUE)[[1]])[3]
		tempevent <- event[2]
		if(basename(tempevent)=="tsd"){
			tempevent <- file.path(dirname(tempevent),"temp")
			}
		}
	else{
		# Get the name of the simulated event and the path to the "temp".directory holding the simulated files to be merged:
		eventname <- rev(strsplit(event[1],"/",fixed=TRUE)[[1]])[3]
		tempevent <- file.path(dirname(event[1]),"temp")
		}
	
	# Get the names of the .beams and .ctd files to read if 'beams' and/or 'ctd' is missing:
	cat("Listing files of the event...\n")
	filesanddir <- list.files(event[1],full.names=TRUE)
	isdir <- file.info(filesanddir)$isdir
	filesanddir <- filesanddir[!isdir]
	#ext = tools::file_ext(filesanddir)
	#ext <- lapply(strsplit(filesanddir, ".", fixed=TRUE), function(x) tail(x,1))
	
	
	filegroups <- echoIBM.getFileTypes(filesanddir, ext=c("beams", "ctd"), key=list(beams=c("bgns", "nrnp", "cali")))
	
	# Read the .beams files:
	cat("Reading beams-files...\n")
	if(length(beams)==0){
		#beamsfiles <- filesanddir[which(ext=="beams")]
		beams <- read.TSDs(filegroups$beams, t="all", var="all", dimension=TRUE)
		}
		
	# Error if pulselength 'sint' and/or absorption factor 'absr' is missing:
	if(any(is.null(beams$sint), is.null(beams$absr))){
		stop(paste0("Elements \"sint\" and/or \"absr\" missing in the event \"", event[1], "\""))
		}	
	# Read the .ctd files:
	cat("Reading ctd-files...\n")
	if(length(ctd)==0){
		#ctdfiles <- filesanddir[which(ext=="ctd")]
		ctd <- read.TSDs(filegroups$ctd, t="all", var="all", dimension=TRUE)
		}
	# Add default speed of sound if no ctd-data are present:
	if(length(ctd)==0){
		ctd <- list(asps=1500)
		}	
	# Read the time points of the event:
	cat("Reading time points...\n")
	utim_all <- UNIX_time(event[1])$U000
	
	
	########## Execution and output ##########
	### Sum up the simulations and add noise: ###
	# Get the names of the .pings files located in the "temp"-directoy:
	cat("Listing temporary files...\n")
	files <- list.files(tempevent, full.names=TRUE, recursive=TRUE)
	filegroups <- echoIBM.getFileTypes(files, ext=c("pings"))
	#ext = tools::file_ext(files)
	#ext <- lapply(strsplit(files, ".", fixed=TRUE), function(x) tail(x,1))
	pingsfiles <- filegroups$pings
	if(length(pingsfiles)==0){
		warning(paste0("There are no pings-files in the directory of temp-files ", tempevent))
		}
	
	# Check whether there is one pingsfile per time step:
	# Get the unix time points of the .pings files:
	cat("Reading time points of the files...\n")
	
	# First a quick test to see if there are an equal number of files as time steps, and file names correspond to a continuous sequence of time steps:
	utim1 <- read.TSD(pingsfiles[1], var="indt")
	utim999 <- read.TSD(pingsfiles[length(pingsfiles)], var="indt", t="end")
	if(utim999$indt - utim1$indt + 1 == length(pingsfiles)){
		indt <- seq(utim1$indt, utim1$indt + length(pingsfiles) - 1)
		utim <- as.list(utim_all)
		}
	else{
		utim <- read.TSDs(pingsfiles, var="utim", t="all", clean=FALSE)
		utim <- utim[names(utim)=="utim"]
		# uniqeutim=nearly.unique(unlist(utim),tolerance=sqrt(.Machine$double.eps))
		uniqeutim <- unique(unlist(utim))
		# Get the time step indexes of the files located in the temp-directory:
		# The old version used indt=which(utim_all %in% uniqeutim), which only worked correctly if the files were ordered in the order of time steps.
		indt <- match(uniqeutim, utim_all)
		}
	
	# 'numt' is the number of time steps, and if t=="all", 't' is set to 1:numt:
	numt <- length(utim_all)
	# If t=="all", all time points are read:
	if(identical(t,"all")){
		t <- 1:numt
		}
	t=t[t>=1 & t<=numt]
	requested_indt <- intersect(t, indt)
	requested_utim <- utim_all[requested_indt]
	lrequested_utim <- length(requested_utim)
	if(length(setdiff(t, indt))>0){
		warning(paste0("The following time steps are not present in the temp-directory:", paste(setdiff(t, indt), collapse=", ")))
		}
	
	# Generate the list of file indices, where there are one list element for each time step, and each element contains a matrix of two columns, the first being the file indices and the second being the time step index in the corresponding file. 
	cat("Generating file and time indices...\n")
	index <- vector("list",numt)
	for(i in seq_along(pingsfiles)){
		for(j in seq_along(utim[[i]])){
			thisindt <- which(utim_all==utim[[i]][j])
			index[[thisindt]] <- rbind(index[[thisindt]], c(i,j))
			}
		}
			
	# Default values for 'parlist':
	parlist <- echoIBM_rexp_defaults(noise=noise, data=beams, parlist=parlist)
	# The seed is not included in the defaults:
	if(!is.null(parlist$seed)){
		if(strff("ind", parlist$seed)){
			seed <- requested_indt
			}
		else if(length(parlist$seed)==length(requested_indt)){
			seed <- parlist$seed
			}
		else{
			warning("Random seed 'seed' in the parameter list 'parlist' did not have the same length as the number of time steps, and was repeated to this length")
			seed <- rep(parlist$seed, length.out=length(requested_indt))
			}
		}
	
	# The default dumpfile name is paste0("dump_merge_",date(),".txt"): 
	if(!file.exists(dirname(dumpfile))){
		dir.create(dirname(dumpfile))
		}
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		write(paste0(
			"########## DUMP FROM MERGING AND ADDING NOISE TO SIMULATION FILES LOCATED IN THE DIRECTORY (WITH SUBDIRECTORIES): ##########\n\n", 
			tempevent, 
			ngettext(length(t), 
				"\n\n##### SIMULATION FILES MERGED FOR TIME STEP: #####\n", 
				"\n\n##### SIMULATION FILES MERGED FOR TIME STEPS: #####\n"),
			paste(t,collapse=", "),
			"\n\n##### MERGING STARTED: #####\n",
			format(as.POSIXlt(starttime,tz="GMT"),usetz=TRUE)
			), 
			dumpfile)
		write("\n\n\n##### FILES MERGED: #####", dumpfile, append=TRUE)
		write(unlist(pingsfiles), dumpfile, append=TRUE)
		}
	
	# Write to the dumpfile (3):
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# If the signal is scaled:
		if(scls!=1){
			write(paste0("\n\n# Signal scaled by the factor ", scls, "\n\n"), dumpfile, append=TRUE)
			}
		
		# Report the noise applied to the simulations:
		implemented_noise <- intersect(gsub("_fun","",noise),c("ms","acex","aex","cex","ex","nr","bg","pn","bk","phase"))
		if(length(implemented_noise)==0){
			implemented_noise <- "none"
			}
		write(paste0("\n\n\n\n\n##### NOISE APPLIED TO THE SIMULATIONS: #####\n", paste(implemented_noise, collapse="\n")), dumpfile, append=TRUE)
		}	
	# Dump information before going through the time step:
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# Information about the noise:
		echoIBM.dump_summary(c(parlist,beams), dumpfile, type="noise", append=TRUE)
		}
	
	
	# Parallel processing using the pblapply() function in the pbapply package:
	if(cores>1){
		# Detect the number of cores and use the minimum of this and the number of requested cores:	
		cores = min(cores, length(t), detectCores())
		# Send blocks of time steps to readAddNoiseWrite():
		nBytes <- 4
		fact_total_vs_vbsc <- 1.1
		i_vec <- splitSeqIntoBlocks(t=seq_along(requested_indt), size1=beams$numb[1] * beams$lenb[1] * nBytes * fact_total_vs_vbsc, size=filesize, blocks=cores)
	}
	else{
		i_vec <- seq_along(requested_indt)
	}
	# Progress bar parallel processing (if cores>1):
	cat("Merging simulated data:\n")
	out <- pblapply(i_vec, readAddNoiseWrite, requested_indt=requested_indt, indt=indt, index=index, pingsfiles=pingsfiles, beams=beams, ctd=ctd, noise=noise, TVG.exp=TVG.exp, parlist=parlist, tempevent=event, ow=ow, cl=cores)
	
	
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
