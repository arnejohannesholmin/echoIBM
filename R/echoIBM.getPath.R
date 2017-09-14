#*********************************************
#*********************************************
#' Returns a list of the following strings: (1) the path to the event, (2) the event name, (3) the event number, (4) the path to the cruise, and (5) the cruise name.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param dir.type  is the name of the directory holding the data files (usually one of "tsd" and "raw")
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
echoIBM.getPath <- function(
	starttime="2015-01-01 00:00:00", 
	origin=NULL, 
	heading=NULL, 
	distance=NULL, 
	speed=NULL, 
	duration=NULL, 
	nodesLocal=NULL, 
	nodesEarth=NULL, 
	utim=NULL, 
	pingduration=1,
	heave=0, 
	rtx=0, 
	rty=0, 
	rtz=NULL, 
	type="v"){
	
	############### LOG: ###############
	# Start: 2017-03-29 - Clean version.

	# Function for applying randmoness to vessel dynamics:
	setVesselDynamics <- function(x, x0, nums, segmentind){
		# Set the default values, usually 0 but for rtz it should be the heading:
		x0 <- rep(x0, length.out=nums)[segmentind]
		# Apply a function of the number of time steps, or simply repeat the numeric value, which is assumed to correspond to the vessel track segments:
		if(is.numeric(x)){
			x <- rep(x, length.out=nums)[segmentind]
		}
		else if(is.function(x)){
			x <- x(length(segmentind))
		}
		# Add the default and the altered values:
		x + x0
	}


	# Set the speed if missing:
	if(length(speed)==0){
		warning("Speed not given and was set to 1 knot = 1 * 1852 / 3600 = 5.14444 m/s")
		speed <- 1 * 1852 / 3600
	}
	
	# If given as nodes in the coordinate system of the Earth, convert to local nodes:
	if(length(nodesEarth)){
		origin <- nodesEarth[1, 1:2]
		nodesLocal <- global2car(nodesEarth, origin=origin)
	}
	# Add the origin to the nodesLocal if not present:
	if(length(nodesLocal) && !all(nodesLocal[1, 1:2] == c(0,0))){
		nodesLocal <- rbind(c(0,0), nodesLocal)
	}
	# If local nodes and origin are given or deduced, extract heading and distance
	if(length(nodesLocal) && length(origin)){
		segments <- diff(nodesLocal[,1:2])
		heading <- atan2(segments[,2], segments[,1])
		distance <- sqrt(rowSums(segments^2))
		duration <- distance / speed
	}
	# If heading, duration and origin is not given, issue an error:
	if(!any(length(duration), length(distance)) && !all(length(origin), length(heading))){
		stop("Dynamics need to be given either as (1) 'nodesEarth', (2) 'nodesLocal' and 'origin', or (3) 'origin', 'heading', and one of 'duration' and 'distance'")
	}
	# If these are present, extract the local nodes:
	else{
		if(length(distance)==0){
			distance <- duration * speed
		}
		nodesLocal <- cbind(c(0, cumsum(distance * cos(heading))), c(0, cumsum(distance * sin(heading))))
	}
	
	# If origin, heading, speed and distance is given, extract duration instead of distance:
	if(length(duration)==0 && length(origin) && length(heading) && length(speed) && length(distance)){
		duration <- distance / speed
	}
	totalduration <- sum(duration)
	
	
	nums <- length(duration)
	# Repeat the speed
	
	# With origin, heading and distance, get the positions from the times (add utim0 later):
	if(length(utim)==0){
		utim0 <- ftim2utim(starttime)
		utim <- seq(0, totalduration, pingduration)
		numt <- length(utim)
		# ingduration <- rep(pingduration, numt)
		# utim0 <- ftim2utim(starttime)
		# utim <- c(0, cumsum(pingduration[-1]))
	}
	else{
		numt <- length(utim)
		utim0 <- utim[1]
		utim <- utim - utim0
	}
	
	# Get cummulative duration:
	cduration <- cumsum(duration)
	# And intervals of duration:
	cduration0 <- c(0, cduration)
	
	# Get the indices of the segments of each time step:
	segmentind <- findInterval(utim, cduration0)
	# And the duration in the present segment for each time step:
	durationInSegment <- utim - cduration0[segmentind]
	# Add the first time to the output time:
	utim <- utim0 + utim
	
	# repeat speed, heave, rtx and rty by segments:
	speed <- rep(speed, nums)[segmentind]
	isp <- speed[segmentind]

	# Set vessel dynamics not derived elsewhere in the function:
	heave <- setVesselDynamics(heave, 0, nums, segmentind)
	rtx <- setVesselDynamics(rtx, 0, nums, segmentind)
	rty <- setVesselDynamics(rty, 0, nums, segmentind)
	rtz <- setVesselDynamics(rtz, heading, nums, segmentind)

	
	# Get the positions
	psxy <- nodesLocal[segmentind,1:2] + durationInSegment * speed[segmentind] * cbind(cos(heading[segmentind]), sin(heading[segmentind]))
	lonlat <- car2global(psxy, origin=origin)
	
	out <- list(isp=isp, psx=psxy[,1], psy=psxy[,2], psz=heave, rtx=rtx, rty=rty, rtz=rtz, lon=lonlat[,1], lat=lonlat[,2])
	names(out) <- paste0(names(out), type)
	out <- c(list(utim=utim, numt=numt), out, list(lon0=origin[1], lat0=origin[2]))
	return(out)
}
