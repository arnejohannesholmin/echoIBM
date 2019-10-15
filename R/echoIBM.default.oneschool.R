#*********************************************
#*********************************************
#' Performs default handeling of the echoIBM acoustic simulator.
#'
#' @param data  is a list of data as read by echoIBM.oneschool().
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param compensated  specifies which rotations are compensated for by the echo sounder. Legal values are "pitch" and "roll", which are both compensated for in the MS70 sonar.
#' @param max.radius  is the estimated maxumum radius of the transducer elements, used to limit the range over which the radius is calculated if only beam width is given.
#' @param pres  is the presicion used in the function integrateonsphere() finding the surface integral of the beam pattern of the fish at r=1.
#' @param max.cells  is the maximum number of cells in the grid used to calculate 'chi'.
#' @param method  is used in beamPattern.TSD().
#' @param dumpfile  is the name of the file to which information and warnings about the simulation is written.
#' @param bptfile  is the name of the file to which 'sllf', 'rad1', 'rad2', 'pbp1' and 'pbp2' is written (use NULL for no writing).
#' @param  Static variables  NA
#' @param  Beam pattern  NA
#' @param  Configuring variables  NA
#' @param  Beam pattern at emission  NA
#' @param  Beam pattern at reception  NA
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR sonR_implemented Simrad_bwt Simrad_dir
#' @importFrom TSD strff write.TSD
#'
#' @export
#' @rdname echoIBM.default.oneschool
#'
echoIBM.default.oneschool <- function(data, esnm="MS70", compensated=c("pitch","roll"), max.radius=0.2, pres=1e-6, max.cells=1e7, method=c("closest","linear"), dumpfile="dump.txt", bptfile=NULL){

	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-03-16 - First version.
	# Update: 2011-08-24 - Changed the handeling of the dumpfile from the old method of writing the dumpfile and deleting at the end of the function, to not writing, in the case that dumpfile has length 0 or 0 characters. This was forced through by the fact that echoIBM.calibrate() uses echoIBM.oneschool.oneping(), which wrote to the dumpfile unintendedly.
	# Update: 2012-02-07 - Changes to the treatment of beam patterns and transducer radii, and added the side lobe level factors 'sllf' used for adjusting side lobe levels to the desired level.
	# Update: 2012-02-08 - Added the option bptfile, which specifies the file name of the file to which the values 'sllf', 'rad1', 'rad2', 'pbp1' and 'pbp2' are to be written (NULL for not writing the file).
	# Last: 2014-07-02 - Added the function Simrad_dir().
	
	
	##################################################
	##################################################
	##### Defaults: #####
	
	# School static:
	default.spow <- 2 # Backscattering cross section proportional to the square of fish size.
	default.tilt <- 0
	# Assuming no significant compression of the length and based on Ona 2003 (herring):
	default.gamw <- -0.23
	default.gaml <- 0
	default.obln <- 5 # Herring?
	default.pbpf <- "ls"
	
	# Beams:
	default.psze <- -7 # G. O. Sars
	default.pbp1 <- "circularPiston"
	default.pbp2 <- "circularPiston"
	sideLobeLevelCircularPiston <- -10*log10(0.01749786)
	default.sllf_MS70 <- c(25/sideLobeLevelCircularPiston,1)
	default.sllf <- c(1,1)
	
	# CTD:
	default.rho0 <- 1026
	default.gacc <- 9.82 # Norewegian sea
	default.hpr0 <- 101325
	default.asps <- 1500
	
	
	##### Error and default handling: #####
	# Information about missing variales without defaults will be collected and printed as an error:
	errors <- NULL
	echoIBM.warnings_warninglist <- NULL
	
	
	### INPUTS REPRESENTING THE SCHOOL: ###
	# --- Static variables ---
	# Extract the possibly frequency dependent function for 'epsl'. If data$epsl is a single numeric value, it is assumed to be frequency independent:
	data$epsl_table <- data$epsl
	data$epsl <- fun_eps(data[c("epsl","grff")],method="linear")
	
	# Extract the possibly frequency dependent function for 'epss'. If data$epss is a single numeric value, it is assumed to be frequency independent:
	# 'epss_table' is only used to print to the dumpfile:
	if(!is.function(data$epss)){
		data$epss_table <- data$epss
		}
	data$epss <- fun_eps(data[c("epss","grff")],method="linear")
	
	# If the optimal backscattering cross section 'sgbs' (sigma_bs) is present, no relation to target size is specified:
	if(!is.null(data$sgbs)){
		data$sigma0mode <- 1
		}
	# If the coefficient 'epss' linking 'sgbs' to fish size is present, it may be a function of frequency, or simply a numeric vector:
	else if(!is.null(data$epss)){
		data$sigma0mode <- 2
		}
	# If the maximum acoustical cross sectional area 'acca' (A_0) is present, no relation to target size is specified:
	else if(!is.null(data$acca)){
		data$sigma0mode <- 3
		}
	# If the coefficient 'epsl' linking 'acca' to fish size is present, it may be a function of frequency, or simply a numeric vector:
	else if(!is.null(data$epsl)){
		data$sigma0mode <- 4
		}
	# Else an error is issued:
	else{
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,"One of \"sgbs\", \"epss\", \"acca\" or \"epsl\" should be present in the input data. \"epss\" set to the function(f) {return(1)}")
		data$epss <- function(f) 1
		data$sigma0mode <- 2
		#errors <- c(errors,"One of \"sgbs\", \"epss\", \"acca\" or \"epsl\" must be present in the input data")
		}
	
	# If 'data$spow' (the power of target size related to the backscatter) is missing, it is defaulted to default.spow=2 for all fish:
	if(length(data$spow)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$spow' was defaulted to ",default.spow,sep=""))
		data$spow <- default.spow
		}
	# If 'data$tilt' (usually tilt of swim bladder) is missing, it is defaulted to default.tilt=0 for all fish:
	if(length(data$tilt)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$tilt' was defaulted to ",default.tilt,sep=""))
		data$tilt <- default.tilt
		}
	# If 'data$gamw' (radial hydrostatic compression factor) is missing, it is defaulted to default.gamw=-0.23 for all fish:
	if(length(data$gamw)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$gamw' was defaulted to ",default.gamw,sep=""))
		data$gamw <- default.gamw
		}
	# If 'data$gaml' (hydrostatic compression factor of the length) is missing, it is defaulted to default.gaml=0 for all fish:
	if(length(data$gaml)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$gaml' was defaulted to ",default.gaml,sep=""))
		data$gaml <- default.gaml
		}
	# If 'data$obln' (oblongness) is missing, it is defaulted to default.obln=5 for all fish:
	if(length(data$obln)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$obln' was defaulted to ",default.obln,sep=""))
		data$obln <- default.obln
		}
	# Error if data$obln<1:
	if(any(data$obln<1)){
		stop("The oblongness of the line source representing the fish (Length/radius) cannot be less than 1")
		}
	# --- Beam pattern ---
	# If 'data$pbpf' (parametric beam pattern) is missing, it is defaulted to default.pbpf="lineSource":
	if(length(data$pbpf)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$pbpf' was defaulted to ",default.pbpf,sep=""))
		data$pbpf <- default.pbpf
		}
	# Interpret the pbpf:
	else{
		# 1. Line source:
		if(tolower(data$pbpf) %in% c("ls", "linesource")){
			data$pbpf <- "ls"
		}
		# 2. Prolate spheriod:
		else if(tolower(data$pbpf) %in% c("pr", "psph", "prsph", "prolatespheroid")){
			data$pbpf <- "pr"
		}
		# 3. Point source:
		else if(tolower(data$pbpf) %in% c("ps", "pointsource")){
			data$pbpf <- "ps"
		}
		else{
			warning("'type' not recognized (must be one of \"ls\"/\"linesource\", \"ps\"/\"pointsource\", or \"pr\"/\"psph\"/\"prsph\"/\"prolatespheroid\"). Defaulted to line source")
			data$pbpf <- "ls"
		}
	}
	# If the beam pattern of the fish is based on the index numbers of the fish, and these are missing, an error is collected:
	if(!length(data$grif)==0 && length(data$indl)==0){
		errors <- c(errors,"'data$indl' missing with no default")
		}
	
	
	### INPUTS REPRESENTING THE ECHO SOUNDER: ###
	# --- Configuring variables ---
	# If the alternative direction angles 'data$dirl' and 'data$dirt' are given, 'data$dira' and 'data$dire' are extracted from these angles. Only the MS70 multibeam sonar and the EK60 mutifrequency echosounder are implemented in this version:
	if(length(esnm)>0){
		data$esnm[1] <- esnm
		}
	
	# Get the azimuth and elevation angles of the beams if missing:
	if(any(length(data$dira)==0, length(data$dire)==0)){
		data <- Simrad_dir(data[c("dirx","diry","dirl","dirt","esnm")])
		}
	# If all direction variables are present, the number of beams 'data$Ni' is obtained as the length of the longest variable, and the shorter are recycled to match the longer. Else an error is collected and the number of beams 'data$Ni' is defaulted to 1 to complete the default actions and error handeling of missing input variables:
	if(!any(length(data$dira)==0,length(data$dire)==0)){
		data$Ni <- max(length(data$dira),length(data$dire))
		data$dira <- rep(data$dira,length.out=data$Ni)
		data$dire <- rep(data$dire,length.out=data$Ni)
		}
	else{
		errors <- c(errors,"One or both of 'data$dira' or 'data$dirl' (azimuth angle) and 'data$dire' or 'data$dirt' (elevation angle) missing with no default")
		data$Ni <- 1
		}
	# If 'data$freq' (frequency) is missing an error is collected:
	if(length(data$freq)==0){
		errors <- c(errors,"'data$freq' (frequency) missing with no default")
		}
	# If 'data$absr' (absorption coefficient) is missing an error is collected:
	if(length(data$absr)==0){
		errors <- c(errors,"'data$absr' (absorption coefficient) missing with no default")
		}
	# If 'data$psze' (z-coordinate of echo sounder) is missing, it is defaulted to default.psze=-7 for all beams:
	if(length(data$psze)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$psze' was defaulted to ",default.psze,sep=""))
		data$psze <- default.psze
		}
	# If 'data$sint' (pulslength) is missing an error is collected:
	if(length(data$sint)==0){
		errors <- c(errors,"'data$sint' (pulslength) missing with no default")
		}
	# If 'data$lenb' (length of beams) is missing an error is collected:
	if(length(data$lenb)==0){
		errors <- c(errors,"'data$lenb' (length of beams) missing with no default")
		}
	# If none of 'data$rad1' (transducer element radii at emission), 'data$bwtx' (transducer element beam widths at emission), or the SIMRAD given 'data$bwtl' and 'data$bwtt' are present, an error is collected. Conversion from beam width to radius is done by the function get.transducerRadius() in the "execution" section:
	if(all(length(data$rad1)==0,length(data$bwtx)==0,length(data$bwtl)==0,length(data$bwtt)==0)){
		errors <- c(errors,"none of 'data$rad1' (after beamforming \"transducer element\" radii at emission), 'data$bwtx' (beam widths at emission), or the 'data$bwtl' and 'data$bwtt' given by SIMRAD are present")
		}
	# If none of 'data$rad2' (transducer element radii at reception) or 'data$bwty' (transducer element beam widths at reception), or the SIMRAD given beam widhts alongship and athwarthship are present, an error is collected. Conversion from beam width to radius is done by the function get.transducerRadius() in the "execution" section:
	if(all(length(data$rad2)==0,length(data$bwty)==0,length(data$bwtl)==0,length(data$bwtt)==0)){
		errors <- c(errors,"none of 'data$rad2' (after beamforming \"transducer element\" radii at reception), 'data$bwty' (beam widths at reception), or the 'data$bwtl' and 'data$bwtt' given by SIMRAD are present")
		}
	
	# Extract the beam widths along the x axis (bwtx) and y axis (bwty) of the coordinate system of the beam, if not already present in the data:
	if(any(length(data$bwtx)==0, length(data$bwty)==0)){
		data = Simrad_bwt(data)
		}
		
	# Default the side lobe level factor for adjusting the side lobe level to the desired values:
	if(length(data$sllf)==0){
		#bptfile=c(bptfile,"new")
		if(sonR_implemented(esnm, "MBS")[1]){
			echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$sllf' was defaulted by",default.sllf_MS70))
			data$sllf <- matrix(default.sllf_MS70,nrow=data$Ni,ncol=2,byrow=TRUE)
			}
		else{
			echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$sllf' was defaulted by",default.sllf))
			data$sllf <- matrix(default.sllf,nrow=data$Ni,ncol=2,byrow=TRUE)
			}
		}
		
	# --- Beam pattern at emission ---
	# If 'data$pbp1' (parametric beam pattern) is missing, it is defaulted according to the type of acoustical instrument:
	if(length(data$pbp1)==0){
		#bptfile <- c(bptfile,"new")
		if(sonR_implemented(esnm)[1]){
			data$pbp1 <- "circularPiston_ellipticRadius_sidelobefit"
			}
		else if(strff("ms70_highsidelobe",data$esnm[1])){
			data$pbp1 <- "circularPiston_ellipticRadius"
			}
		else if(any(strff(c("ek60_circular","me70_circular","ms70_circular"),data$esnm[1]))){
			data$pbp1 <- "circularPiston"
			}
		else{
			echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$pbp1' was defaulted by",default.pbp1))
			data$pbp1 <- default.pbp1
			}
		}
	# --- Beam pattern at reception ---
	# If 'data$pbp1' (parametric beam pattern) is missing, it is defaulted according to the type of acoustical instrument:
	if(length(data$pbp2)==0){
		#bptfile <- c(bptfile,"new")
		if(sonR_implemented(esnm)[1]){
			data$pbp2 <- "circularPiston_ellipticRadius_sidelobefit"
			}
		else if(strff("ms70_highsidelobe",data$esnm[1])){
			data$pbp2 <- "circularPiston_ellipticRadius"
			}
		else if(any(strff(c("ek60_circular","me70_circular","ms70_circular"),data$esnm[1]))){
			data$pbp2 <- "circularPiston"
			}
		else{
			echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$pbp2' was defaulted by",default.pbp2))
			data$pbp2 <- default.pbp2
			}
		}
	# If the beam pattern of the echosounder/sonar is based on the index numbers of the beams, and these are missing, they are defaulted to 1:Ni:
	if(!length(data$gri1)==0 || !length(data$gri2)==0 && length(data$indi)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,"'data$indi' was defaulted")
		data$indi <- 1:data$Ni
		}
	
	
	### INPUTS REPRESENTING THE SEA: ###
	# If 'data$rho0' (equilibrium mass density of the sea) is missing, it is defaulted to default.rho0=1026:
	if(length(data$rho0)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$rho0' was defaulted to ",default.rho0,sep=""))
		data$rho0 <- default.rho0
		}
	# If 'data$gacc' (gravitational accelleration) is missing, it is defaulted to default.gacc=9.82:
	if(length(data$gacc)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$gacc' was defaulted to ",default.gacc,sep=""))
		data$gacc <- default.gacc
		}
	# If 'data$hpr0' (hydrostatic pressure at sealevel) is missing, it is defaulted to default.hpr0=101325:
	if(length(data$hpr0)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$hpr0' was defaulted to ",default.hpr0,sep=""))
		data$hpr0 <- default.hpr0
		}
	# If 'data$asps' (average speed of sound) is missing, it is defaulted to default.asps=1500:
	if(length(data$asps)==0){
		echoIBM.warnings_warninglist <- c(echoIBM.warnings_warninglist,paste("'data$asps' was defaulted to ",default.asps,sep=""))
		data$asps <- default.asps
		}
	
	
	### Report errors and terminate function if needed: ###
	if(!length(errors)==0){
		cat("The following errors occured:\n")
		for(i in 1:length(errors)){
			cat(paste(i,"-",errors[i]),"\n")
			}
		return()
		}
	
	
	### Beam patterns: ###
	# Generating the beam pattern function, which is either parametric or empirical, and has one of four modes numbered by 1, 2, 3 and 4: 
	#	(1) phi,ka - function of one angle of direction 'phi' and the product 'ka' of wave number and size of the source
	#	(2) phi,ind - function of one angle of direction 'phi' and the numbering index 'ind' of the source
	#	(3) theta,phi,ka - function of two angles of direction 'theta' and 'phi', and the product 'ka' of wave number and size of the source
	#	(4) theta,phi,ind - function of two angles of direction 'theta' and 'phi', and the numbering index 'ind' of the source
	if(length(data$ebpf)>0){
		data$ebpf <- data$ebpf/max(data$ebpf)
		}
	data$bptf <- beamPattern.TSD(data[c("pbpf","graf","gref","grsf","grif","dbpf","ebpf")],method=method)
	
	# Generating the beam pattern function, which is either parametric or empirical, and has one of four modes numbered by 1, 2, 3 and 4: 
	#	(1) phi,ka - function of one angle of direction 'phi' and the product 'ka' of wave number and size of the source
	#	(2) phi,ind - function of one angle of direction 'phi' and the numbering index 'ind' of the source
	#	(3) theta,phi,ka - function of two angles of direction 'theta' and 'phi', and the product 'ka' of wave number and size of the source
	#	(4) theta,phi,ind - function of two angles of direction 'theta' and 'phi', and the numbering index 'ind' of the source
	data$bpt1 <- beamPattern.TSD(data[c("pbp1","gra1","gre1","grs1","gri1","dbp1","ebp1")],method=method)
	
	# Generating the beam pattern function, which is either parametric or empirical, and has one of four modes numbered by 1, 2, 3 and 4: 
	#	(1) phi,ka - function of one angle of direction 'phi' and the product 'ka' of wave number and size of the source
	#	(2) phi,ind - function of one angle of direction 'phi' and the numbering index 'ind' of the source
	#	(3) theta,phi,ka - function of two angles of direction 'theta' and 'phi', and the product 'ka' of wave number and size of the source
	#	(4) theta,phi,ind - function of two angles of direction 'theta' and 'phi', and the numbering index 'ind' of the source
	data$bpt2 <- beamPattern.TSD(data[c("pbp2","gra2","gre2","grs2","gri2","dbp2","ebp2")],method=method)
	
	
	##### Transformations and rotations: #####
	# Wave number for each beam:
	data$wavenumber <- rep(2 * pi * data$freq / matrix(data$asps, nrow=NROW(data$freq), ncol=NCOL(data$freq), byrow=TRUE), l=data$Ni)
	# Identifying beams of equal frequency:
	data$uniquek <- unique(data$wavenumber)
	data$luniquek <- length(data$uniquek)
	# The frequencies will be arranged in a matrix having equal frequencies along the first or the second dimension. If equal frequencies are along the first dimension the array of beam patterns calculated at a later stage in the simulation must be permuted using aperm().
	# Assign indexes for the beams relating them to the frequencies in data$uniquek:
	data$matchk <- match(data$wavenumber,data$uniquek)
	# Frequency table of the frequencues. If not all equal, there is a varying number of beams for each frequency, which is not supported (because it requires using lists instead of arrays for the calculations, which is a lot slower):
	data$kfreq <- table(data$matchk)
	if(!all(data$kfreq==data$kfreq[1])){
		stop("Unequal number of beams for each frequency is not supported in this version")
		}
	# Only keep the first value of 'data$kfreq' for convenience:
	data$kfreq <- data$kfreq[1]
		
	# Calculate the radii of the transducer elements, as modeled by the circular piston with radius varying elliptically as a function of azimuth angle for the MS70 sonar. These are given as a two column matrix of beam widths in the x-direction and the y-direction in the coordinate system of the transducer elements:
	# The default is to assume identical beam pattern at emission and reception. This code must be changed if exeptions to this assumption are implemented in the future:
	# The beam angles called "beam width along ship" ('bwtl' in the TSD convension) and "beam width atwarth ship" ('bwtt' in the TSD convension) are defined for the MS70 sonar in a different way than for echo sounders and the ME70. The MS70 is mounted so that 'bwtl' is the beam angle vertically, while 'bwtt' is the beam angle horizontally, that is oposite the intuitive definition. Also the beam widths stored in the raw-files are one-way, and thus the beam angles needed in the simulation are known directly from the raw-files (as oposed to what was beleived to be the case before 2011-03-16)
	if(length(data$rad1)==0){
		#bptfile <- c(bptfile,"new")
		data$rad1 <- get.transducerRadius2.TSD(c(beampattern=data$bpt1, max.radius=max.radius, twoway=FALSE, data[c("bwtx","bwty","sllf","wavenumber")]))
		}
	if(length(data$rad2)==0){
		#bptfile <- c(bptfile,"new")
		# Here it is assumed that the beam patterns at emission and reception are equal. Change this code if the oposite is implemented in the future:
		# data$rad2 <- get.transducerRadius2.TSD(c(beampattern=data$bpt2, max.radius=max.radius, twoway=FALSE, data[c("bwtx","bwty","sllf","wavenumber")]))
		data$rad2 <- data$rad1
		}
	
	# If data$rad1 and data$rad2 are equal a simplification is done in the caluculations:
	if(identical(data$rad1,data$rad2)){
		data$equalbp_em_re <- TRUE
		}
	else{
		data$equalbp_em_re <- FALSE
		}
	# Calculating the surface integrals 'ssif':
	if(length(data$ssil)==0 && length(data$ssif)==0 && any(data$sigma0mode %in% 3:4)){
		data$ssif <- write.ssif(x=data, outfile=FALSE, pres=pres, max.cells=max.cells, onlyssif.out=TRUE)
		}
		
	# Print the warnings to the dumpfile:
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		write("\n\n\n##### WARNINGS IN \"echoIBM.default.oneschool\": #####",dumpfile,append=TRUE)
		if(length(echoIBM.warnings_warninglist)>0){
			for(i in seq_along(echoIBM.warnings_warninglist)){
				write(paste(i,": ",echoIBM.warnings_warninglist[i],sep=""),dumpfile,append=TRUE)
				}
			}
		else{
			write("none",dumpfile,append=TRUE)
			}
		}
	
	# Write the beam pattern specific variables calculated in this function to a file for copying to other simulation cases which use identical beam configurations. The calculation of these variables is time consuming, and time is saved if this file is alreaddy present at the start of the simulation:
	if(length(bptfile)>0 && !file.exists(bptfile)){
		write.TSD(data[c("sllf","rad1","rad2","pbp1","pbp2")], bptfile[1], numt=1)
	}
	
	# Output:
	data
	##################################################
	##################################################
	}
