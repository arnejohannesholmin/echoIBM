#*********************************************
#*********************************************
#' Returns a list of the following strings: (1) the path to the event, (2) the event name, (3) the event number, (4) the path to the cruise, and (5) the cruise name.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param dir.type  is the name of the directory holding the data files (usually one of "tsd" and "raw")
#' @param ...  is used in agrep() for locating events based on approximate string matching.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD pathparts rm.na
#'
#' @export
#' @rdname echoIBM.setup
#'
echoIBM.setBeams <- function(
	eventDir, 
	eventName = basename(eventDir), 
	utim = NULL, 
	esnm = "EK60", 
	data = list(), 
	bmmd = NULL){
	
	############### LOG: ###############
	# Start: 2017-03-29 - Clean version.
	copyForOneEsnm <- function(esnm, files, data, bmmd, numt, eventDir, eventName){
		# Match 'esnm' against the pre-defined systems:
		files <- files[grep(esnm, basename(files), ignore.case=TRUE)]
	
		# Read the resource beams file if there was a match. These files must have one time step per beam mode 'bmmd'. Use drop.out=FALSE to allow for selecting time steps using 'bmmd' below:
		if(length(files)){
			beams <- read.TSD(files[1], t="all", drop.out=FALSE)
		}
		else{
			warning(paste0("No beam configuration files matching the specified system (", esnm, "). Available files are the following:", paste(basename(files), collapse="\n")))
		}
		# Add data:
		beams <- c(data, beams)
		if(length(beams)==0){
			stop("'data' must be given as a list of beam configuration data. One list per system (if not only one system is used) specifying variables such as those named by labl.TSD(\"rb\").")
		}
		# Repeat the first ping if bmmd is not given:
		if(length(bmmd)==0){
			bmmd <- rep(1, numt)
		}
		# Here we utilize the drop.out=FALSE used in read.TSD() above:
		beams <- lapply(beams, function(xx) xx[,bmmd])
	
		subevent <- file.path(eventDir, esnm, "tsd")
		file <- file.path(subevent, paste0(eventName, "_", esnm, ".beams"))
		write.TSD(beams, file, numt=numt)
		return(file=file, subevent=subevent)
	}
	
	# 'utim' must be given:
	if(length(utim)==0){
		stop("'utim' must be given and must contain all the time steps of the simulation")
	}
	else{
		numt <- length(utim)
	}
	
	# Get available beams files:
	files <- list.files(system.file("extdata", "beams", package="echoIBM"), full.names=TRUE)
	
	outfiles <- rep(NA, length(esnm))
	outdirs <- rep(NA, length(esnm))
	for(i in seq_along(esnm)){
		temp <- copyForOneEsnm(esnm[i], files=files, data=if(length(esnm)==1) data else data[[esnm[i]]], bmmd=bmmd, numt=numt, eventDir=eventDir, eventName=eventName)
		outfiles[i] <- temp$file
		outdirs[i] <- temp$subevent
	}
	names(outfiles) <- esnm
	names(outfiles) <- esnm
	return(list(outfiles=outfiles, outdirs=outdirs))
	### # Match 'esnm' against the pre-defined systems:
	### availableBeamsFiles <- grep(esnm, availableBeamsFiles, ignore.case=TRUE)
	### 
	### # Read the resource beams file if there was a match:
	### if(length(availableBeamsFiles)){
	### 	beams <- read.TSD(availableBeamsFiles, t="all", drop.out=FALSE)
	### }
	### # Add data:
	### beams <- c(data, beams)
	### # Repeat the first ping if bmmd is not given:
	### if(length(bmmd)==0){
	### 	bmmd <- rep(1, numt)
	### }
	### beams <- lapply(beams, function(xx) xx[,bmmd])
	### 
	### beamsfile <- file.path(eventDir, paste(esnm, "beams", sep="."))
	### write.TSD(beams, beamsfile, numt=numt)
}












### 
### 
### 
### 
### ############################################################################################
### ########## Create the beams resource files for the different implemented systems: ##########
### ############################################################################################
### #rep.col<-function(x,n){
### #   matrix(rep(x, each=n), ncol=n, byrow=TRUE)
### #}
### 
### 
### getMissingNames <- function(x){
### 	beamsnames <- c("esnm", "utim", "asps", "numb", "indi", "freq", "absr", "sint", "rres", "plsl", "psze", "lenb", "dirx", "diry", "dirz", "dira", "dire", "bwtl", "bwtt", "bwth", "bwtv", "bwtx", "bwty", "eqba", "sacr", "tpow", "gain", "gai1", "gai2", "bmmd")
### 	n <- names(x)
### 	n <- n[sapply(x, length)>0]
### 	setdiff(beamsnames, n)
### }
### 
### beamsdir <- "/Applications/echoIBM/Frameworks/R_packages/echoIBM/data/beams"
### 
### 
### 
### get.rres.TSD <- function(x){
### 	if(is.list(x)){
### 		x$asps * x$sint / 2
### 	}
### 	else{
### 		return(NULL)
### 	}
### }
### 
### 	# EK60:
### 	beamsEK60 <- read.TSD("~/Data/echoIBM/Resources/Beams/S2014119_D200141030_E0021_EK60.beams", t=1)
### 	beamsEK60$bmmd <- 1
### 	beamsEK60$utim <- utim.TSD(beamsEK60)
### 	beamsEK60$rres <- get.rres.TSD(beamsEK60)
### 	getMissingNames(beamsEK60) # "bwth" "bwtv" "gai1" "gai2"
### 	f <- file.path(beamsdir, "EK60.TSD")
### 	write.TSD(beamsEK60, f, numt=1, keep.null=FALSE)
### 	dim_all(read.TSD(f, drop.out=FALSE))
### 
### # ME70:
### beamsEK60 <- read.TSD("~/Data/echoIBM/Resources/Beams/S2014119_D200141030_E0021_EK60.beams", t=1)
### beamsEK60$bmmd <- 1
### beamsEK60$utim <- utim.TSD(beamsEK60)
### beamsEK60$rres <- get.rres.TSD(beamsEK60)
### getMissingNames(beamsEK60)
### f <- file.path(beamsdir, "EK60.TSD")
### write.TSD(beamsEK60, f, numt=1, keep.null=FALSE)
### dim_all(read.TSD(f, drop.out=FALSE))
### 
### 	# MS70:
### 	beamsMS70 <- read.TSD("~/Data/echoIBM/Resources/Beams/S2009116_D20091113_E0001_MS70.beams", t=1)
### 	beamsMS70$bmmd <- 3
### 	beamsMS70$utim <- utim.TSD(beamsMS70)
### 	beamsMS70$rres <- get.rres.TSD(beamsMS70)
### getMissingNames(beamsMS70)
### f <- file.path(beamsdir, "MS70.TSD")
### write.TSD(beamsMS70, f, numt=1, keep.null=FALSE)
### dim_all(read.TSD(f, drop.out=FALSE))
### 
### # SX90:
### beamsEK60 <- read.TSD("~/Data/echoIBM/Resources/Beams/S2014119_D200141030_E0021_EK60.beams", t=1)
### beamsEK60$bmmd <- 1
### beamsEK60$utim <- utim.TSD(beamsEK60)
### beamsEK60$rres <- get.rres.TSD(beamsEK60)
### getMissingNames(beamsEK60)
### f <- file.path(beamsdir, "EK60.TSD")
### write.TSD(beamsEK60, f, numt=1, keep.null=FALSE)
### dim_all(read.TSD(f, drop.out=FALSE))
### 
### # SU90:
### 	beamsSU90 <- read.TSD("/Volumes/Acoustics/S2015116_PG.O.Sars[4174]/Events/E0022/SU90/tsd/E0022_SU90.beams", t=1:2)
### 	beamsSU90$utim <- utim.TSD(beamsSU90)
### 	beamsSU90$rres <- get.rres.TSD(beamsSU90)
### 	getMissingNames(beamsSU90)
### 	f <- file.path(beamsdir, "SU90.TSD")
### write.TSD(beamsSU90, f, numt=1, keep.null=FALSE)
### dim_all(read.TSD(f, drop.out=FALSE))
### 
### 
### 
### 
### 
### 
### 
### 
### 
### 
### 
### 
### # Write new TSD file from EK60, ME70, MS70, SX90 and SU90, and use the beams-files:
### rawfiles <- list(
### 	EK60 <- "/Volumes/Acoustics/Generate_beams_files_for_echoIBM/Events/S2009116_D20091113_E0001/EK60/raw",
### 	ME70 <- "/Volumes/Acoustics/Generate_beams_files_for_echoIBM/Events/X2010001_TobisTest_E0005_20100425/ME70/raw",
### 	MS70 <- "/Volumes/Acoustics/Generate_beams_files_for_echoIBM/Events/S2009116_D20091113_E0001/MS70/raw", 
### 	SX90 <- "/Volumes/Acoustics/Generate_beams_files_for_echoIBM/Events/D20130316_D20130316_E0001/SX90/raw", 
### 	SU90 <- "/Volumes/Acoustics/Generate_beams_files_for_echoIBM/Events/S2015116_PG.O.Sars[4174]_E0001/SU90/raw")
### 
### 
### 	lapply(rawfiles, EKRaw2TSD)
### 
### 
### 
### 
### 
### 
### 
### c("esnm", "utim", "asps", "numb", "indi", "freq", "absr", "sint", "rres", "plsl", "psze", "lenb", "dirx", "diry", "dirz", "dira", "dire", "bwtl", "bwtt", "bwth", "bwtv", "bwtx", "bwty", "eqba", "sacr", "tpow", "gain", "gai1", "gai2", "bmmd")
### 
### 
### 
### 
### 
### 
### 
### esnm
### utim
### asps
### numb
### indi
### freq
### absr
### sint
### rres
### plsl
### psze
### lenb
### dirx
### diry
### dirz
### dira
### dire
### bwtl
### bwtt
### bwth
### bwtv
### bwtx
### bwty
### eqba
### sacr
### tpow
### gain
### gai1
### gai2
### bmmd
### 
### 
### 
### 
### list(utim=utim, ispv=ispv, psxv=psxyv[,1], psyv=psxyv[,2], pszv=heave, rtxv=rtxv, rtyv=rtyv, rtzv=rtzv, lonv=lonlat[,1], latv=lonlat[,2], lon0=origin[1], lat0=origin[2])
### 
