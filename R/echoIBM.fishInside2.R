#*********************************************
#*********************************************
#' Discards fish outside of the sampling region of the sonar/echosounder FOR EACH SCHOOL.
#'
#' @param data  is a list holding the variables 'esnm', and vessel spesifications.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom sonR is.sonar rotate3D
#' @importFrom TSD car2sph
#'
#' @export
#' @rdname echoIBM.fishInside
#'
echoIBM.fishInside2 <- function(data, dumpfile, discardOutside=c(r=Inf,az=Inf,el=Inf), rand.sel=1){
	
	############### LOG: ###############
	# Start: 2014-02-27 - Clean version.
	
	# Extract a random selection of the targets using 'rand.sel':
	Nl <- max(length(data$psxf), length(data$psyf), length(data$pszf))
	selection <- seq_len(Nl)
	
	if(0<rand.sel[1] && rand.sel[1]<1){
		if(length(rand.sel)>1){
			set.seed(rand.sel[2])
		}
		selection <- sample(seq_len(Nl), round(Nl * rand.sel[1]))
	}
	

	thisdata <- data[c("dira", "dire")]
	
	
	if(!any(is.infinite(discardOutside[2:3]))){
		
		# Expand to the edges of the observation volume, expanded by the 'discardOutside':
		max_range <- (max(data$lenb) + discardOutside[1]) * soundbeam_range(data, pos="rres")
		max_beamwidth_azimuth <- max(data$bwtl) * pi/180
		max_beamwidth_elevation <- max(data$bwtt) * pi/180
		mean_dira <- mean(thisdata$dira)
		mean_dire <- mean(thisdata$dire)
	
	
		
		dira_full <- c(outer(thisdata$dira, c(-1, 0, 1) * max_beamwidth_azimuth * discardOutside[2], "+"))
		dire_full <- c(outer(thisdata$dire, c(-1, 0, 1) * max_beamwidth_elevation * discardOutside[3], "+"))
	
		# Assure that the angles are on the unit sphere:
		dira_full[dira_full < 0] <- 0
		dira_full[dira_full > 2*pi] <- 2*pi
		dire_full[dire_full < 0] <- 0
		dire_full[dire_full > pi] <- pi
	
		# Add the axes between the beams:
		addAxes <- function(ang){
			axis1 <- ceiling(min(ang) / (pi/2))
			axis2 <- floor(max(ang)   / (pi/2))
			if(axis2 >= axis1){
				ang <- c(ang, seq(axis1, axis2) * (pi/2))
			}
			ang
		}
		#temp <- addAxes(dira=dira_full, dire=dire_full, mean_dira=mean_dira, mean_dire=mean_dire)
		dira_full <- addAxes(dira_full)
		dire_full <- addAxes(dire_full)
		#dira_full <- temp$dira
		#dire_full <- temp$dire
		edges <- cbind(max_range, expand.grid(dira_full, dire_full))
	
	
		# Convert to Cartesian positions in the global coordinate system, and use the resulting box to discard fish:
		edgesXY <- sph2car(edges)
		# Add position of the sonar/echosounder:
		#edgesXY <- rbind(c(0, 0, data$psze), edgesXY)
		edgesXY <- rbind(0, edgesXY)
		edgesXY[,3] <- edgesXY[,3] + data$psze
		# Add the position of the vessel:
		edgesXY <- edgesXY + matrix(c(data$psxv, data$psyv, data$pszv), ncol=3, nrow=nrow(edgesXY), byrow=TRUE)
	
		# Define a box surrounding the edges of the sonar outside which fish are discarded:
		box <- apply(edgesXY, 2, range)
	
		inside <- 
			data$psxf[selection] >= box[1,1] & 
			data$psxf[selection] <= box[2,1] & 
			data$psyf[selection] >= box[1,2] & 
			data$psyf[selection] <= box[2,2] & 
			data$pszf[selection] >= box[1,3] & 
			data$pszf[selection] <= box[2,3]
	
		if(length(dumpfile)>0 && nchar(dumpfile)>0){
			write(paste0("\n\n# Proportion of fish discarded outside of the observation volume: ", mean(!inside)),dumpfile,append=TRUE)
			}
				
		selection <- selection[inside]
	}
	
	
	# The dynamic variable names of the school, and legal time variable names:
	dynschoolnames <- labl.TSD("ds")
	# The static variable names of the school:
	staticschoolnames <- labl.TSD("ss")
	
	# Return the dynamic school data and the other data:
	affected.variables <- c(dynschoolnames, staticschoolnames)
	for(j in seq_along(affected.variables)){
		thisvar <- affected.variables[j]
		if(length(data[[thisvar]])==Nl && !is.function(data[[thisvar]])){
			data[[thisvar]] <- data[[thisvar]][selection]
		}
	}
	
	return(data)
}
