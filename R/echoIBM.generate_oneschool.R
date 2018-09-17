#*********************************************
#*********************************************
#' Generates the position variables of one school. The fish are first put on an iso grid defined by the packing density. Then noise is added to these positions by the Gaussian distribution with mean 0 and standard deviation specified the inputs 'SDxf', 'SDyf', 'SDzf', (repeated to length 3 and defaulted to 1 if missing), all in units of the grid size. Finally, the orientations of the fish is set independently of the positions, by generating two sets of positions which are gridded positions separated by one grid unit in the x-direction and added noise by a standard deviation producing the desired polarization (by a pre-generated table).
#'
#' @param data  is a list containing the necessary data: 
#' @param school  is the index number of the school.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR vl2rt.TSD
#' @importFrom TSD labl.TSD read.TSD sph2car strff
#' @importFrom stats approx rnorm
#'
#' @export
#' @rdname echoIBM.generate_oneschool
#'
echoIBM.generate_oneschool <- function(data, school=1, seed=NULL, dumpfile=NULL, aboveSurface=c("cut","force")){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-02-05 - First version.
	# Last: 2014-07-01 - Changed to generate noise on positions and orientation of fish separately.
	########### DESCRIPTION: ###########
	# Generates the position variables of one school. The fish are first put on an iso grid defined by the packing density. Then noise is added to these positions by the Gaussian distribution with mean 0 and standard deviation specified the inputs 'SDxf', 'SDyf', 'SDzf', (repeated to length 3 and defaulted to 1 if missing), all in units of the grid size. Finally, the orientations of the fish is set independently of the positions, by generating two sets of positions which are gridded positions separated by one grid unit in the x-direction and added noise by a standard deviation producing the desired polarization (by a pre-generated table).
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the necessary data: 
	#	One of 'nbfS' or 'rhoS', 
	#	all of 'szxS', 'szyS', 'szzS', 
	#	'psxS', 'psyS', 'pszS', 
	#	'nbfS', 'SDxf', 'SDyf', 'SDzf', 
	#	'thtS', 'phiS', 
	#	and optionally 'plHS'.
	# ---school--- is the index number of the school.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Select the given school:
	cs <- labl.TSD("cs")
	data[cs] <- lapply(data[cs], function(x) if(length(x)==1) x else x[school])
	# Set the seed of the school:
	#set.seed(data$seed)
	
		
	########## Execution ##########
	# If the packing density is missing, calculate it from the number of fish:
	if(length(data$rhoS)==0 && length(data$nbfS)>0 && !is.na(data$nbfS)){
		# The volume of the school:
		if(strff("e",data$shpS)){
			# Volume of an ellipsoid:
			semiaxisX <- data$szxS/2
			semiaxisY <- data$szyS/2
			semiaxisZ <- data$szzS/2
			data$volS <- 4/3*pi * semiaxisX * semiaxisY * semiaxisZ
		}
		else if(strff("r",data$shpS)){
			# Volume of a box:
			data$volS <- data$szxS * data$szyS * data$szzS
		}
		data$rhoS <- data$nbfS/data$volS
	}
	spacing <- data$rhoS^(-1/3)
	# First generate gridded positions:
	getGridOneDimensionRandom <- function(size, spacing){
		# Pick a random position of the grid
		start <- runif(1, -size/2, -size/2 + spacing)
		
		numGrid <- ceiling(size / spacing)
		
		grid <- seq(0, numGrid) * spacing + start
		grid <- grid[-size/2 <= grid  &  grid <= size/2]
		
		grid
	}
	# First generate gridded positions:
	getGridOneDimension <- function(size, spacing){
		nhalf <- floor(size/2 / spacing)
		
		grid <- seq(-nhalf, nhalf) * spacing
		
		grid
	}
	
	
	
	
	grid <- as.matrix(expand.grid(
		getGridOneDimension(size=data$szxS, spacing=spacing), 
		getGridOneDimension(size=data$szyS, spacing=spacing), 
		getGridOneDimension(size=data$szzS, spacing=spacing)
		))
		
		#grid <- as.matrix(expand.grid(seq(-data$szxS/2, data$szxS/2, spacing), seq(-data$szyS/2, data$szyS/2, spacing), seq(-data$szzS/2, data$szzS/2, spacing)))
	
	if(strff("e",data$shpS)){
		# Crop off to the ellipsoid:
		d <- matrix(c(data$szxS, data$szyS, data$szzS) / 2, byrow=TRUE, nrow=nrow(grid), ncol=ncol(grid))
		inside <- rowSums((grid / d)^2) <= 1
		#inside <- (grid[,1]/(data$szxS/2))^2 + (grid[,2]/(data$szyS/2))^2 + (grid[,3]/(data$szzS/2))^2 <= 1
		grid <- grid[inside,]
	}
	
	# Rotate the school:
	if(length(data$rtzS) && !all(data$rtzS==0)){
		grid <- rotate3D(grid, by="z", ang=data$rtzS)
	}		
	if(length(data$rtxS) && !all(data$rtxS==0)){
		grid <- rotate3D(grid, by="x", ang=data$rtxS)
	}		
	if(length(data$rtyS) && !all(data$rtyS==0)){
		grid <- rotate3D(grid, by="y", ang=data$rtyS)
	}		
	
	# Add the position of the school:
	data$psxf <- grid[,1] + data$psxS
	data$psyf <- grid[,2] + data$psyS
	data$pszf <- grid[,3] + data$pszS
	if(length(grid)==0){
		warning(paste0("School nr ", school, " was too small to be generated with the current packing density."))
		return()
	}
	#data$nbfS <- nrow(grid)
	
	
	# Fish above the sea surface:
	if(tolower(substr(aboveSurface[1],1,1))=="f"){
		data$pszf[data$pszf>0] <- 0
	}
	else{
		valid <- data$pszf<=0
		long <- sapply(data,length)==length(valid)
		data[long] <- lapply(data[long], "[", valid)
		#data$nbfS <- sum(valid)
	}
	
	
	# Discard fish outside of the radial and angular range of the sonar/echosounder:
	data <- echoIBM.fishInside2(data=data, dumpfile=dumpfile, discardOutside=data$discardOutside, rand.sel=data$rand.sel)
	
	# Add randomness in fish positions and orientations:
	data <- echoIBM.addRandomness(data, grsd_plHS=NULL)
	
	
	########## Output ##########
	data[labl.TSD("echoIBM.generate_oneschool_labl")]
	##################################################
	##################################################
	}
