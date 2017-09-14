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
echoIBM.setVessel_old <- function(events, 
	eventName, 
	starttime="2015-01-01 00:00:00", 
	origin=NULL, 
	heading=NULL, 
	distance=NULL, 
	speed=NULL, 
	duration=NULL, 
	nodesLocal=NULL, 
	nodesEarth=NULL, 
	utim=NULL, 
	numt=1,
	pingduration=1,
	heave=0, 
	rtxv=0, 
	rtyv=0, 
	rtzv=NULL, 
	...){
	
	############### LOG: ###############
	# Start: 2017-03-29 - Clean version.

	# Set the speed if missing:
	if(length(speed)==0){
		warning("Speed not given and was set to 10 knots = 10 * 1852 / 3600 = 5.14444 m/s")
		speed <- 10 * 1852 / 3600
	}
	
	# If given as nodes in the coordinate system of the Earth, convert to local nodes:
	if(length(nodesEarth)){
		origin <- nodesEarth[1, 1:2]
		nodesLocal <- global2car(nodesEarth, origin=origin)
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
		stop("Vessel dynamics need to be given either as (1) 'nodesEarth', (2) 'nodesLocal' and 'origin', or (3) 'origin', 'heading', and one of 'duration' and 'distance'")
	}
	# If these are present, extract the local nodes:
	else{
		if(length(distance)==0){
			distance <- duration * speed
		}
		#nodesLocal <- rbind(c(0,0), distance * cbind(cos(heading), sin(heading)))
		nodesLocal <- cbind(c(0, cumsum(distance * cos(heading))), c(0, cumsum(distance * sin(heading))))
	}
	
	# If origin, heading, speed and distance is given, extract duration instead of distance:
	if(length(duration)==0 && length(origin) && length(heading) && length(speed) && length(distance)){
		duration <- distance / speed
	}
	
	
	nsegments <- length(duration)
	# Repeat the speed
	
	# With origin, heading and distance, get the positions from the times:
	if(length(utim)==0){
		pingduration <- rep(pingduration, numt)
		utim <- c(0, cumsum(pingduration[-1]))
	}
	
	# Get cummulative duration:
	cduration <- cumsum(duration)
	# And intervals of duration:
	cduration0 <- c(0, cduration)
	
	# Get the indices of the segments of each time step:
	segmentind <- findInterval(utim, cduration0)
	# And the duration in the present segment for each time step:
	durationInSegment <- utim - cduration0[segmentind]
	
	if(length(rtzv)==0){
		rtzv <- heading[segmentind]
	}
	else{
		rtzv  <- rep(rtzv, nsegments)[segmentind]
	}
	
	# repeat speed, heave, rtx and rty by segments:
	speed <- rep(speed, nsegments)[segmentind]
	ispv <- speed[segmentind]
	heave <- rep(heave, nsegments)[segmentind]
	rtxv <- rep(rtxv, nsegments)[segmentind]
	rtyv <- rep(rtyv, nsegments)[segmentind]
	utim <- ftim2utim(starttime) + utim
	
	
	# Get the positions
	psxyv <- nodesLocal[segmentind,1:2] + durationInSegment * speed[segmentind] * cbind(cos(heading[segmentind]), sin(heading[segmentind]))
	lonlat <- car2global(psxyv, origin=origin)
	
	vessel <- list(utim=utim, ispv=ispv, psxv=psxyv[,1], psyv=psxyv[,2], pszv=heave, rtxv=rtxv, rtyv=rtyv, rtzv=rtzv, lonv=lonlat[,1], latv=lonlat[,2], lon0=origin[1], lat0=origin[2])
	lapply(events, function(event) write.TSD(vessel, file.path(event, paste0(eventName, ".vessel")), numt=numt))
}
