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
echoIBM.generate_oneschool<-function(data, school=1, seed=NULL, dumpfile=NULL, aboveSurface=c("cut","force")){
	
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
	data[cs] <- lapply(data[cs], "[", school)
	# Set the seed of the school:
	set.seed(data$seed)
	
		
	########## Execution ##########
	# Apply the default standard deviation value:
	if(all(sapply(data[c("SDxf","SDyf","SDzf")],length)==0)){
		data$SDxf <- 1
	}
	# Repeat the standard deviation value:
	u <- unlist(data[c("SDxf","SDyf","SDzf")])
	if(length(u)<3){
		u <- rep(u,length.out=3)
		data$SDxf <- u[1]
		data$SDyf <- u[2]
		data$SDzf <- u[3]
	}
	
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
	grid <- as.matrix(expand.grid(seq(-data$szxS/2, data$szxS/2, spacing), seq(-data$szyS/2, data$szyS/2, spacing), seq(-data$szzS/2, data$szzS/2, spacing)))
	
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
	grid[,1] <- grid[,1] + data$psxS
	grid[,2] <- grid[,2] + data$psyS
	grid[,3] <- grid[,3] + data$pszS
	if(length(grid)==0){
		warning(paste0("School nr ", school, " was too small to be generated with the current packing density."))
		return()
	}
	data$nbfS <- nrow(grid)
	
	# Get the positions of the fish, added noise:
	data$psxf <- grid[,1] + rnorm(data$nbfS, 0, data$SDxf)
	data$psyf <- grid[,2] + rnorm(data$nbfS, 0, data$SDyf)
	data$pszf <- grid[,3] + rnorm(data$nbfS, 0, data$SDzf)
	
	# Fish above the sea surface:
	if(tolower(substr(aboveSurface[1],1,1))=="f"){
		data$pszf[data$pszf>0] <- 0
	}
	else{
		valid <- data$pszf<=0
		long <- sapply(data,length)==length(valid)
		data[long] <- lapply(data[long], "[", valid)
	}
	
	# Add noise to the fish positions. If the polarization is given, calculate the standard deviation on the positions:
	if(length(data$plHS)>0){
		#write(c(school, data$plHS, data$plHS), "test", append=TRUE, ncolumns=20)
		#write("\n\n", "test", append=TRUE)
		file <- system.file("extdata", "grsd_plHS.TSD", package="echoIBM")
		#echoIBM_datadir_ <- file.path(echoIBM_frameworks, "R", "Functions", "echoIBM Main", "Utilities")
		#filebasename <- "grsd_plHS.TSD"
		#file <- file.path(echoIBM_datadir_, filebasename)
		grsd_plHS <- read.TSD(file)
		# Get the standard deviation of the displacement of the fish:
		### if(min(grsd_plHS$plHS)>data$plHS || max(grsd_plHS$plHS)<data$plHS){
		### 	stop(paste("Polarization values outside of the valid range (0,pi/2) radians. Range:", paste(round(range(data$plHS),digits=3),collapse=", ")))
		### }
		# Temporarily redefine standard deviations to use when generating the polarization. This does not affect the fish positions, but is only used to derive the velocities of the fish:
		data$SDxf <- data$SDyf <- data$SDzf <- approx(grsd_plHS$plHS, grsd_plHS$grsd, data$plHS, rule=2)$y
		
		# Get the second positions of the fish, added noise and the movement of the school, used to calculate the velocities of the fish:
		# Displacement equal to 1 is used in "grsd_plHS.TSD":
		displacement_car <- sph2car(c(1, data$thtS, data$phiS))
		data$psxf1 <- grid[,1] + rnorm(data$nbfS, 0, data$SDxf)
		data$psyf1 <- grid[,2] + rnorm(data$nbfS, 0, data$SDyf)
		data$pszf1 <- grid[,3] + rnorm(data$nbfS, 0, data$SDzf)
		data$psxf2 <- grid[,1] + rnorm(data$nbfS, 0, data$SDxf) + displacement_car[1]
		data$psyf2 <- grid[,2] + rnorm(data$nbfS, 0, data$SDyf) + displacement_car[2]
		data$pszf2 <- grid[,3] + rnorm(data$nbfS, 0, data$SDzf) + displacement_car[3]
		data$vlxf <- data$psxf2 - data$psxf1
		data$vlyf <- data$psyf2 - data$psyf1
		data$vlzf <- data$pszf2 - data$pszf1
	}
	# Set the fish to simulated omnidirectional if polarization is missing:
	else{
		data$vlxf <- rnorm(data$nbfS)
		data$vlyf <- rnorm(data$nbfS)
		data$vlzf <- rnorm(data$nbfS)
	}
	
	# Add rotation data:
	data[c("rtzf","rtxf")] <- vl2rt.TSD(data[c("rtzf","rtxf","vlxf","vlyf","vlzf")])[c("rtzf","rtxf")]
	
	# Add sizes of the fish:
	if(length(data$PDsz)==0){
		thisPDsz <- "rnorm"
	}
	else{
		thisPDsz <- data$PDsz
	}
	data$size <- do.call(thisPDsz, list(data$nbfS, data$MEsz, data$SDsz))
	#data$size <- rnorm(data$nbfS,data$MEsz,data$SDsz)
	
	
	########## Output ##########
	data[labl.TSD("echoIBM.generate_oneschool_labl")]
	##################################################
	##################################################
	}
