#*********************************************
#*********************************************
#' The aspect factor (cross sectional area relative to the maximum cross sectional area) of a cylindrical object of rounded ends, a prolate spheroid or a point source (circular). See the documentation of echoIBM for details (Section 3.1.4).
#'
#' @param data  is a list containing the elements "obln", specifying the oblongness of the object (ratio of length and thickness) and "transducerposL", giving the position of the transducet in the coordinate system of the target in spherical coordinates. Also "pbpf" may be given in 'data'.
#' @param pbpf  is a string naming the theoretical acoustical model of the target (may be included in 'data').
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ones
#'
#' @export
#' @rdname etaOmega
#'
etaOmega <- function(data, pbpf=c("ls","pr","ps")){
	
	############### LOG: ###############
	# Start: 2010-12-17 - Clean version.
	# Update: 2011-01-18 - Added depth dependence on the aspect factor.
	# Last: 2017-05-03 - Cleaned ut the code.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(!is.list(data)){
		stop("'data' must be a list")
		}
	if(length(grep("pbpf",names(data)))>0){
		pbpf <- data$pbpf
		}
	# Depth compression factor of the oblongness:
	data$obln <- data$obln * etaCompression(data, type="obln")
	pbpf <- tolower(pbpf)
	
	
	########## Execution and output ##########
	##### 1. Line source: #####
	if(data$pbpf == "ls"){
		# Multiply the oblongness by 2 to get the ratio L/a, where 'L' is the length and 'a' is the radius of the line source:
		# Subtract the hemispheres at the ends:
		xi2 <- 2 * data$obln - 2
		out <- (pi + 2 * xi2 * sin(data$transducerposL[,3])) / (pi + 2 * xi2)
	}
	##### 2. Prolate spheriod: #####
	else if(data$pbpf == "pr"){
		# Set the angle of the tangent point on the prolate spheroid of the widest axis of the ellipse projected on the plane normal to the viewpoint:
		obln2 <- data$obln^2
		theta <- atan(1/tan(data$transducerposL[,3])/obln2)
		# See the documentation of echoIBM for details (Section 3.1.4):
		out <- sqrt(1 / (obln2-(obln2-1) * cos(theta)^2)) * sin(theta + data$transducerposL[,3])
	}
	##### 3. Point source: #####
	else if(data$pbpf == "ps"){
		out <- ones(length(data$transducerposL[,3]))
	}
	return(out)
	
	
	### if(tolower(data$pbp) %in% c("l", "ls", "linesource", "linesource.tsd")){
	### 	# Multiply the oblongness by 2 to get the ratio L/a, where 'L' is the length and 'a' is the radius of the line source:
	### 	# Subtract the hemispheres at the ends:
	### 	xi2 <- 2 * data$obln - 2
	### 	(pi + 2 * xi2 * sin(data$transducerposL[,3])) / (pi + 2 * xi2)
	### 	}
	### else if(tolower(data$pbp) %in% c("l", "ls", "linesource", "linesource.tsd")){
	### 	ones(length(data$transducerposL[,3]))
	### 	}
	### else if(length(grep("prolatespheroid",pbpf[1]))>0 || length(grep("psph",pbpf[1]))>0){
	### 	# Set the angle of the tangent point on the prolate spheroid of the widest axis of the ellipse projected on the plane normal to the viewpoint:
	### 	obln2 <- data$obln^2
	### 	theta <- atan(1/tan(data$transducerposL[,3])/obln2)
	### 	# See the documentation of echoIBM for details (Section 3.1.4):
	### 	sqrt(1 / (obln2-(obln2-1) * cos(theta)^2)) * sin(theta + data$transducerposL[,3])
	### 	}
	### else if(length(agrep("linesource",pbpf[1]))>0){
	### 	# Multiply the oblongness by 2 to get the ratio L/a, where 'L' is the length and 'a' is the radius of the line source:
	### 	# Subtract the hemispheres at the ends:
	### 	xi2 <- 2 * data$obln - 2
	### 	(pi + 2 * xi2 * sin(data$transducerposL[,3])) / (pi + 2 * xi2)
	### 	}
	### else if(length(agrep("prolatespheroid",pbpf[1]))>0){
	### 	# Set the angle of the tangent point on the prolate spheroid of the widest axis of the ellipse projected on the plane normal to the viewpoint:
	### 	obln2 <- data$obln^2
	### 	theta <- atan(1/tan(data$transducerposL[,3])/obln2)
	### 	# See the documentation of echoIBM for details (Section 3.1.4):
	### 	sqrt(1/(obln2-(obln2-1)*cos(theta)^2)) * sin(theta+data$transducerposL[,3])
	### 	}
	### else if(length(agrep("pointsource",pbpf[1]))>0){
	### 	ones(length(data$transducerposL[,3]))
	### 	}
	### else{
	### 	warning("'type' not recognized using fuzzy matching. Defaulted to line source")
	### 	# Multiply the oblongness by 2 to get the ratio L/a, where 'L' is the length and 'a' is the radius of the line source:
	### 	# Subtract the hemispheres at the ends:
	### 	xi2 <- 2 * data$obln - 2
	### 	(pi + 2 * xi2 * sin(data$transducerposL[,3])) / (pi + 2 * xi2)
	### 	}
	##################################################
	##################################################
	}
