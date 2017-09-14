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
#' @importFrom utils tail
#' @importFrom stats runif
#'
#' @export
#' @rdname echoIBM.setup
#'

# 1. school size distribution + school density + packing density + school position distribution



echoIBM.setSchool <- function(
	events, 
	eventName, 
	vessel, 
	surveyRegion = NULL, 
	target = "herring",
	schools = "layer",
	schoolSize = list(x = 100, y = 100, z = 20, type = "Weibull", shape = 5), 
	data = list(static = list(), dynamic = list()),
	seed = 0, 
	onlyFirstEvent = FALSE){
	
	############### LOG: ###############
	# Start: 2017-03-29 - Clean version.
	
	
	############################################################
	############### (1) Treat static school data: ##############
	############################################################
	# Apply default settings specified by 'target':
	# Herring:
	if(tolower(target) == "herring"){
		schoolStatic <- list(
			pbpf = "pr",
			obln = 5,
			tilt = 0, 
			# Define compression of the swim bladder in case the simulate targets have one (defaults taken from Ona 2003, with compression only radially in the swim bladder):
			gamw = -0.23,
			gaml = 0, 
			# Ratio of swim bladder versus total length (Gorska and Ona 2003, Modelling the effect of swimbladder compression on the acoustic backscattering from herring at normal or near-normal dorsal incidences).
			zeta = 0.26,
			# The backscattering coefficient of each individual target can be given directly through the sgbs variable:
			sgbs = NULL,
			# Otherwise use the link between fish size and sgbs:
			epss = function(f){
				10^(-6.54)*100^2 * (f/38000)^(-0.4)
			}
			)
		schoolDynamic <- list(
			MEsz = 0.32,
			SDsz = 0.02,
			PDsz = "rnorm"
			)
	}
	else if(tolower(target) == "mackerel"){
		schoolStatic <- list(
			pbpf = "pr",
			obln = 5,
			tilt = 0, 
			# Define compression of the swim bladder in case the simulate targets have one (defaults taken from Ona 2003, with compression only radially in the swim bladder):
			gamw = 0,
			gaml = 0, 
			# Ratio of swim bladder versus total length, set to 1 when no swimbladder:
			zeta = 1,
			# The backscattering coefficient of each individual target can be given directly through the sgbs variable:
			sgbs = NULL,
			# Otherwise use the link between fish size and sgbs:
			epss = function(f){
				NA
			}
			)
		schoolDynamic <- list(
			MEsz = 0.32,
			SDsz = 0.02,
			PDsz = "rnorm"
			)
	}
	if(tolower(target) == "point"){
		schoolStatic <- list(
			pbpf = "ps",
			obln = 1,
			tilt = 0, 
			# Define compression of the swim bladder in case the simulate targets have one (defaults taken from Ona 2003, with compression only radially in the swim bladder):
			gamw = 0,
			gaml = 0, 
			# Gorska and Ona 2003, Modelling the effect of swimbladder compression on the acoustic backscattering from herring at normal or near-normal dorsal incidences.
			zeta = 1,
			# The backscattering coefficient of each individual target can be given directly through the sgbs variable:
			sgbs = NULL,
			# Otherwise use the link between fish size and sgbs:
			epss = function(f){
				10^(-6.54)*100^2 * (f/38000)^(-0.4)
			}
			)
		schoolDynamic <- list(
			MEsz = 0.32,
			SDsz = 0,
			PDsz = "rnorm"
			)
	}
	
	# If prolate spheroid is given for the parametric beam pattern of the targets, include the beam pattern file of the given aspect ratio:
	if(schoolStatic$pbpf == "pr"){
		targetBeamPatternFiles <- list.files(system.file("extdata", "TargetBeamPattern", package="echoIBM"), full.names=TRUE)
		# Split the basenames by "_": 
		targetBeamPatternbasenameSplitted <- strsplit(basename(targetBeamPatternFiles), "_")
		obln <- as.numeric(sapply(targetBeamPatternbasenameSplitted, function(x) x[which(tolower(x)=="obln")[1] + 1]))
		selected <- which.min(abs(obln - schoolStatic$obln))
		if(obln[selected] != schoolStatic$obln){
			warning(paste0("The target beam pattern file of the closest oblongness chosen (", targetBeamPatternFiles[selected], ")"))
		}
		lapply(events, function(event) file.copy(targetBeamPatternFiles[selected], event))
	}
		
	# Add custom static data:
	schoolStatic <- c(data$static, schoolStatic)
	# Write the static school file to the events:
	lapply(if(onlyFirstEvent) events[1] else events, function(event) write.TSD(schoolStatic, file.path(event, paste0(eventName, "_SchoolStatic", ".school")), numt=1))
	############################################################
	############################################################
	
	
	############################################################
	############## (2) Treat dynamic school data: ##############
	############################################################
	# Get total rectangular survey region:
	if(length(surveyRegion)==0){
		surveyRegion <- data.frame(x=range(vessel$psxv), y=range(vessel$psyv))
	}
	
	# Special case if a layer of fish is requsted. In that case rectangular schools are defined, so that the entire observation region is filled 
	if(schools=="layer"){
		# Divide the survey region into rectangles with sizes given by 'schoolSize':
		gridx <- seq(surveyRegion$x[1], surveyRegion$x[2], schoolSize$x[1])
		if(tail(gridx, 1) < surveyRegion$x[2]){
			gridx <- c(gridx, surveyRegion$x[2])
		}
		gridy <- seq(surveyRegion$y[1], surveyRegion$y[2], schoolSize$y[1])
		if(tail(gridy, 1) < surveyRegion$y[2]){
			gridy <- c(gridy, surveyRegion$y[2])
		}
		# Get mid points and sizes from the grid:
		midgridx <- (gridx[-1] + gridx[-length(gridx)]) / 2
		midgridy <- (gridy[-1] + gridy[-length(gridy)]) / 2
		midgridz <- mean(schoolSize$z[1:2])
		midgrid <- expand.grid(midgridx, midgridy, midgridz)
		Nschool <- nrow(midgrid)
		sizex <- diff(gridx)
		sizey <- diff(gridy)
		sizez <- diff(schoolSize$z[1:2])
		size <- expand.grid(sizex, sizey, sizez)
		
		schoolDynamic <- c(schoolDynamic, list(
			# Layer block positions:
			psxS = midgrid[,1],
			psyS = midgrid[,2],
			pszS = midgrid[,3],
			# Layer block sizes:
			szxS = size[,1],
			szyS = size[,2],
			szzS = size[,3],
			# Layer block rotations (no rotation):
			rtxS = 0,
			rtyS = 0, 
			rtzS = 0, 
			# Use rectangular shaped "schools" for the layer:	 
			shpS = "r",
			# SD of mean 0-Gaussian distribution of positions:
			SDxf = 1,
			SDxf = 1,
			SDxf = 1, 
			# Speed, heading, density, polarization:
			ispS = 0,
			thtS = runif(N, 0, 2*pi),
			phiS = pi/2,
			rhoS = 1,
			plHS = Inf,
			# Add volS also as a global variable for use in nbfS. Volume 'volS' is only used to derive the approximate memory used:
			volS = volS <- schools$szxS * schools$szyS * schools$szzS,
			nbfS = volS * volS
			))
		
	}
	else if(schools=="uniform"){
		
		schoolDynamic <- c(schoolDynamic, list(
			# Layer block positions:
			psxS = midgrid[,1],
			psyS = midgrid[,2],
			pszS = midgrid[,3],
			# Layer block sizes:
			szxS = size[,1],
			szyS = size[,2],
			szzS = size[,3],
			# Layer block rotations (no rotation):
			rtxS = 0,
			rtyS = 0, 
			rtzS = 0, 
			# Use rectangular shaped "schools" for the layer:	 
			shpS = "r",
			# SD of mean 0-Gaussian distribution of positions:
			SDxf = 1,
			SDxf = 1,
			SDxf = 1, 
			# Speed, heading, density, polarization:
			ispS = 0,
			thtS = runif(N, 0, 2*pi),
			phiS = pi/2,
			rhoS = 1,
			plHS = Inf,
			# Add volS also as a global variable for use in nbfS. Volume 'volS' is only used to derive the approximate memory used:
			volS = volS <- schools$szxS * schools$szyS * schools$szzS,
			nbfS = volS * volS
			))
		
		
		
		
		
		
		
		
		
		
		
		
		
		# Model the schools as elipsoids:
		volS = 4/3*pi * schools$szxS/2 * schools$szyS/2 * schools$szzS/2
		
		
		
		
		
		
		# Ellipsoidal schools (http://keisan.casio.com/has10/SpecExec.cgi?path=05000000.Mathematics%2F01000300.Volume%20and%20surface%20area%2F13000700.Volume%20of%20an%20ellipsoidal%20cap%2Fdefault.xml&charset=utf-8):
		axA <- schools$szxS/2
		axB <- schools$szyS/2
		axC <- schools$szzS/2
		above <- axC - schools$pszS
		volS <- pi/3 * axA*axB * above^2/axC^2 * (3*axC - above)
		volE <- 4*pi/3 * axA * axB * axC
		#### <- schools$rhoS * volS
		
		
		
	}
	else if(schools=="flatUniform"){
		
	}
	
	
	# Define the seeds to use in when generating the schools and individual fish:
	set.seed(seed)
	seed <- as.list(round(runif(length(Nschool), 0, 1e6)))
	
	
	
	
	
	# Add the maximum school size on all sides of the surevy region:
	max
	
	
	Nschool <- 120
	
	
	
	# Set the sizes of the schools:
	size <- as.matrix(sapply(schoolSize[c("x", "y", "z")], rep, length.out=Nschool))
	if(strff("Wei", schoolSize$type)){
		# w=1 implies no correlation between schools:
		size <- rexp_MultSines(J=Nschool, I=3, L=3, N=40, P=10, w=1, olpn=c(0.5,0.5,0.5), shape=schoolSize$shape[1], mean=size, seed=seed)
	}
	szxS <- size[,1]
	szyS <- size[,2]
	szzS <- size[,3]
	
	
	#pp(1)
	#plot(NULL, xlim=c(0,1300), ylim=c(0,1300))
	#for(i in seq_len(nrow(rr))){
	#	lines(ellipse(runif(2, 100,1200), rr[i,1], rr[i,2]))
	#}
	#
	#hist(rr[,1]/rr[,2], breaks=1000, xlim=c(0,5))
	#		summary(rr[,1]/rr[,2], breaks=40)
	#		sd(rr[,1]/rr[,2])
	#		cvar(rr[,1]/rr[,2])
	
	
	
	
	
	
	
	
	##############################
	#### (3e) SCHOOL DYNAMIC: ####
	if(length(schools)==0){
		# The schools are here specified in compact form, with one set of specifications for each school, which are used at each time in the simulations to generate the school on the fly, and not pre-generating it with ARsim_school() or some other method.

		# The following variables should be set for the compatly specified schools:
		# (I) Unix time point of the following data (utmS)
		# (II) Size (szxS, szyS, szzS)
		# (III) Position (psxS, psyS, pszS)
		# (IV) Orientation (thtS, phiS)
		# For now all schools have constant speed, and changing this sould involve assigning one speed value and one direction to each of a set of time intervals.
		# (V) Speed (aspS)
		# (VI) Density (rhoS)
		# (VII) Polarization (Huth and Wissel 1992 definition, also used by Parrish et al. 2002 and Holmin et al. 2012) (plHS)
		# (VIIIa) Mean fish size (MEsz)
		# (VIIIb) Standard deviation of fish sizes (SDsz)
		# (VIIIc) Seed for the fish sizes (seed)
		# (IX) (optional) Volume (volS)
		# (X) (optional) Number of fish (nbfS)
		trackLength = diff(range(vessel$psxv))
		N = ceiling(trackLength / schooDist)

		# The y positions should be uniformly distributed, and also schools that are only partially inside the sonar volume should be included. This will depend on the sizes of the schools. So if a school that has center 50 m outside of the maximum range of the sonar has size <100 m, it will not intersect with the sonar volume, whereas if it has size >100 m, it will. Therefore, some schools will be too far away, and will result in empty simulations. To avoid this, we generate too many schools, and replace those that are to far away by closer ones.
		# We need to define the maximum range of all acoustic instruments to be simulated, and the maximum size of the schools, and use these to find the preliminary number of schools:

		# Similarly, the schools that intersect with the surface should be treated with care. If these schools are removed, or relocated randomly so that the do not intersect with the surface, there will be an upper layer with lower density of fish. To acheive uniformly distributed fish positions, schools will be positioned also above the surface, and all fish above the surface will be excluded when echoIBM.generate_oneschool() is run.

		# Thus we need to make the school centers span maxRange + maxSize/2 to each side of the vessel track, and the following preliminary number of schools should be sufficient. However, we will apply the 
		preN = 3 * N


		###############
		## (I) utmS: ##
		###############
		# Pick the start time of the simulation project as the time of all of the schools:
		schools$utmS=rep(starttime, N)

		############################
		## (II) szxS, szyS, szzS: ##
		############################
		# To include a correlation between dimensions of the schools we use a Gaussian copula with beta distributed marginals. The marginals specify the lower (l) and upper (u) limit of the size in each dimension, and the position of the mode (m).

		# Here we use lower correlation between z and (x,y), to allow for some tall but small schools:
		schools$szxS = rep(maxSize, N)
		schools$szyS = rep(maxSize, N)
		schools$szzS = rep(maxHeight, N)

		#############################
		## (III) psxS, psyS, pszS: ##
		#############################
		set.seed(seed)
		# Place schools along the vessel track regularly by the value of 'schooDist':
		schools$psxS = schooDist/2 + schooDist * seq(0, N-1)
		# Draw uniformly distributed y and z positions:
		preRange = maxRange + maxSize/2
		schools$psyS = runif(preN, -preRange, preRange)
		schools$pszS = sample(seq(minDepth, maxHeight/2, l=preN))

		# Remove the schools that never intersect with the sonar volume:
		valid = which(abs(schools$psyS) - schools$szyS/2 < maxRange)[seq_len(N)]
		N = length(valid)
		#bufferdepth = -1
		#valid = which(abs(schools$psyS) - schools$szyS/2 < maxRange  &  schools$pszS - schools$szzS/2 < bufferdepth)[seq_len(N)]

		# Select sizes for only the valid schools:
		schools$szxS = schools$szxS[valid]
		schools$szyS = schools$szyS[valid]
		schools$szzS = schools$szzS[valid]
		# And for the positions:
		schools$psyS = schools$psyS[valid]
		schools$pszS = schools$pszS[valid]
		schools$pszS = sample(seq(minDepth, maxHeight/2, l=length(valid)))


		######################
		## (IV) thtS, phiS: ##
		######################
		schools$thtS = runif(N, 0, 2*pi)
		schools$phiS = zeros(N) + pi/2

		###############
		## (V) aspS: ##
		###############
		schools$aspS = zeros(N) + 1e-9

		################
		## (VI) rhoS: ##
		################
		# Gives approcimately mean density of one fish per cubic meters:
		schools$rhoS = ones(N)

		#################
		## (VII) plHS: ##
		#################
		# Set the polarization values by the beta distribution between 10 and 50 with mode at 30, in degrees. This seems simple, and does not exaggerate the polarization as this would increase variance in the results:
		# Must be given in radians:
		schools$plHS = zeros(N) + 30 * pi/180

		###################
		## (VIIIa) MEsz: ##
		###################
		schools$MEsz = zeros(N) + 0.32

		###################
		## (VIIIb) SDsz: ##
		###################
		schools$SDsz = rep(0.02,N)

		###################
		## (VIIIc) seed: ##
		###################
		schools$seed = seq_len(N)

		################
		## (IX) volS: ##
		################
		schools$volS = 4/3*pi * schools$szxS/2 * schools$szyS/2 * schools$szzS/2

		###############
		## (X) nbfS: ##
		###############
		schools$nbfS=NAs(N)
		for(i in seq_len(N)){
			cat(i,", ",sep="")
			spacing = schools$rhoS[i]^(-1/3)
			# First generate gridded positions:
			grid = as.matrix(expand.grid(seq(-schools$szxS[i]/2,schools$szxS[i]/2,spacing), seq(-schools$szyS[i]/2,schools$szyS[i]/2,spacing), seq(-schools$szzS[i]/2,schools$szzS[i]/2,spacing)))
			# Crop off to the ellipsoid, and remove fish above the surface:
			inside = (grid[,1]/(schools$szxS[i]/2))^2 + (grid[,2]/(schools$szyS[i]/2))^2 + (grid[,3]/(schools$szzS[i]/2))^2 <= 1
			below = grid[,3] + schools$pszS[i] <= 0
			schools$nbfS[i] = sum(inside & below)
			}
		cat("\n")
		}
	
	# Write the school file:
	schoolFiles = file.path(events, paste0(project, "_", eventspec, "_schools.school"))
	write.TSD(schools, schoolFiles[1], numt=1)
	
	
}
