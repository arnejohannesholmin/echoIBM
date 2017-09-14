echoIBM.setCTD <- function(events, 
	eventName, 
	data = list(),
	...){

	copyForOneEsnm <- function(esnm, files, data, eventDir, eventName){
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


	
	
	files <- list.files(system.file("extdata", "CTD", package="echoIBM"), full.names=TRUE)

	CTD<- read.TSD(files[1], t="all", drop.out=FALSE)


	CTD <- c(data, CTD)
	if(length(CTD)==0){
		stop("'data' must be given as a list of beam CTD data specifying variables such as those named by labl.TSD(\"ctd\").")
	}
	
	subevent <- file.path(eventDir, esnm, "tsd")
	file <- file.path(subevent, paste0(eventName, "_", esnm, ".beams"))
	write.TSD(beams, file, numt=1)
	return(beamsfile=beamsfile, subevent=subevent)
	
	
	
	
	
	
	
	
	
	
	
	
		
	vessel <- echoIBM.getPath(starttime=starttime, origin=origin, heading=heading, distance=distance, speed=speed, duration=duration, nodesLocal=nodesLocal, nodesEarth=nodesEarth, utim=utim, pingduration=pingduration, heave=heave, rtx=rtx, rty=rty, rtz=rtz, type="v")
	vesselfiles <- file.path(event, paste0(eventName, ".vessel"))

	# Run through the events and write the vessel file:
	lapply(events, function(event) write.TSD(vessel, con=vesselfiles, numt=vessel$numt))

	# Add the files to the output:
	vessel$files <- vesselfiles
	return(vessel)
	
		if(!identical(CTD, FALSE)){
			CTDFile = "/Volumes/Acoustics/S2014119_PG.O.Sars[4174]/Events/S2014119_D200141030_E0021/SX90/tsd/S2014119_D200141030_E0021_SX90.ctd"
		for(i in seq_along(events)){
				file.copy(CTDFile, events[i], overwrite=TRUE)
				}
			}
		}