#*********************************************
#*********************************************
#' Generates and writes vessel files for all acoustic instruments of a simulation event. The function \code{echoIBM.getPath} generates the vessel information for a specific esnm.
#'
#' The vessel information can be given in 4 ways: 
#' \itemize{
#'	\item nodesEarth + speed
#'	\item origin + nodesLocal + speed
#'	\item origin + heading + distance + speed
#'	\item origin + heading + duration + speed
#'	}
#'
#' @param eventName		The name of the simulation event, used in the file name of the vessel file.
#' @param starttime		The start time of the event, i.e., the time of the first ping/vessel position, either given as a time object or a string such as "2015-01-01 00:00:00". 
#' @param utim			The UNIX time points of the event. If not given, the starttime and pingduration will define the UNIX time points.
#' @param origin		The origin of the event, i.e., a vector of two elements giving the longitude and latitude of the first ping.
#' @param heading		The heading of the vessel in radians counter clockwise from North, either given as a single numeric, or as a vector of headings associated with vessel track segments.
#' @param distance		The distance associated with each vessel track segment.
#' @param duration		The duration of each vessel track segment.
#' @param speed			The speed of the vessel in knots associated with each vessel track segment.
#' @param nodesLocal	The nodes defining the vessel track segments, given as local nodes relative to the origin (can include the origin c(0, 0))
#' @param nodesEarth	The nodes defining the vessel track segments, given as global nodes in the coordinate system of the earth relative to the origin (can include the origin c(0, 0)).
#' @param pingduration	The duration of the pings in seconds. Currently only one single fixed value is allowed.
#' @param heave			The heave of the vessel in meters, either given as a single value, a vector of length equal to the number of vessel track segments, or alternatively a funciton of the number of time steps such as rnorm() or the more appropriate function(x) {set.seed(x); runif(x)}, which sets the seed as the number of time steps.
#' @param rtxv			The pitch of the vessel in radians positive for uppwards pitch (bow lifting). See info.TSD("rtxv")
#' @param rtyv			The roll of the vessel in radians positive for starboard side tilted down. See info.TSD("rtyv")
#' @param rtzv			The orientation of the vessel in radians counter clockwise, IN ADDITION to the heading of the vessel specified in \code{heading}. See info.TSD("rtzv")
#' @param ...			Further parameters passed to or from other methods.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.setup
#'
echoIBM.setVessel <- function(
	events, 
	eventName, 
	starttime = "2015-01-01 00:00:00", 
	utim = NULL, 
	origin = NULL, 
	heading = NULL, 
	distance = NULL, 
	duration = NULL, 
	speed = NULL, 
	nodesLocal = NULL, 
	nodesEarth = NULL, 
	pingduration = 1,
	heave = 0, 
	rtxv = 0, 
	rtyv = 0, 
	rtzv = 0, 
	...){
	
	############### LOG: ###############
	# Start: 2017-03-29 - Clean version.
	
	vessel <- echoIBM.getPath(starttime=starttime, origin=origin, heading=heading, distance=distance, speed=speed, duration=duration, nodesLocal=nodesLocal, nodesEarth=nodesEarth, utim=utim, pingduration=pingduration, heave=heave, rtx=rtx, rty=rty, rtz=rtz, type="v")
	vesselfiles <- file.path(event, paste0(eventName, ".vessel"))
	
	# Run through the events and write the vessel file:
	lapply(events, function(event) write.TSD(vessel, con=vesselfiles, numt=vessel$numt))
	
	# Add the files to the output:
	vessel$files <- vesselfiles
	return(vessel)
}
#'
#' @importFrom TSD car2global ftim2utim
#'
#' @export
#' @rdname echoIBM.setup
#'
echoIBM.getPath <- function(
	starttime = "2015-01-01 00:00:00", 
	origin = NULL, 
	heading = NULL, 
	distance = NULL, 
	speed = NULL, 
	duration = NULL, 
	nodesLocal = NULL, 
	nodesEarth = NULL, 
	utim = NULL, 
	pingduration = 1,
	heave = 0, 
	rtx = 0, 
	rty = 0, 
	rtz = 0, 
	type = "v"){
	
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
		warning("Speed not given and was set to 10 knots = 10 * 1852 / 3600 = 5.14 m/s")
		speed <- 10 * 1852 / 3600
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
