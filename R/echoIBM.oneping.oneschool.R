#*********************************************
#*********************************************
#' Simulates one echo sounder observation of based on positions, orientations, sizes and other specifics of each fish in one known (simulated) school.
#'
#' @param data  see below.
#' @param esnm  is the name of the echo sounder. Currently implemented are "MS70", "ME70", "EK60", "MS70_circular" (may be given in lover case)
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#' @param compensated  specifies which rotations are compensated for by the echo sounder. Legal values are "pitch" and "roll", which are both compensated for in the MS70 sonar.
#' @param calibrate  is FALSE if calibration data are to be discarded if present.
#' @param noise  is a vector of character strings of length 2, specifying which types of noise to apply to the data:
#' @param max.memory  is the maximum amount of memory to allow for the function to occupy. The level of for loops is chosen according to this value, i.e. if the method using only two loops, one for the radial distances and one for the unique frequencies, demands more memory than 'max.memory', the method usint three for loops is chosen.
#' @param dumpfile  is the name of the file to which information and warnings about the simulation is written.
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#' @param data  is the list of input parameters giving the information about the school consisting of 'Nl' fish.
#' @param ask  is TRUE if the used should be asked to for approval if the memory of the least memory demanding calculation method of the individual radial sampling intervals exceed the memory limit 'max.memory'.
#' @param  Dynamic variables  NA
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom gdata object.size
#' @importFrom SimradRaw apply.TVG soundbeam_range
#' @importFrom sonR rotate3D vl2rt.TSD
#' @importFrom TSD car2sph global2car labl.TSD ones zeros
#'
#' @export
#' @rdname echoIBM.oneping.oneschool
#'
echoIBM.oneping.oneschool <- function(data, esnm=NULL, TVG.exp=2, compensated=c("pitch", "roll"), calibrate=TRUE, noise=c("nr", "bg", "ex"), max.memory=1e9, dumpfile="dump.txt", ask=FALSE, parlist=list()){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-09-16 - First version.
	# Update: 2010-03-22 - Changed method to reduce memory occupied by the function.
	# Update: 2011-02-20 - Changed to use beamPattern.TSD() which takes as input a list of inputs named by the convencional TSD-names. Thus the parameter "mod" as a list element of the beam pattern function is no longer needed.
	# Update: 2011-03-17 - Corrected the calculation of ciston radii, based on the new information about the beam widths in the raw files for the MS70. EK60 and ME70 not yet altered.
	# Update: 2011-08-24 - Changed the handeling of the dumpfile from the old method of writing the dumpfile and deleting at the end of the function, to not writing, in the case that dumpfile has length 0 or 0 characters. This was forced through by the fact that echoIBM.calibrate() uses echoIBM.oneschool.oneping(), which wrote to the dumpfile unintendedly.
	# Update: 2011-10-14 - Fixed bug when 'size' not given. 'size' and 'lenl' are now repeated to the number of fish.
	# Update: 2011-12-01 - Added the parameter 'scls' in 'data' to reduce CPU time.
	# Update: 2012-12-03 - Added speed by subsetting the data used in echoIBM.oneping.j1() and echoIBM.oneping.j2().
	# Update: 2014-03-14 - Changed to comply with the new echoIBM.oneping.j1() using svnext, which calculates the echo only once for the two voxels that are affected by each target.
	# Update: 2015-02-22 - Changed to support grid calibration values.
	# Last: 2015-02-24 - Fixed bug in the radialindex, which had one too many null elements at the beginning of the list.
	########### DESCRIPTION: ###########
	# Simulates one echo sounder observation of based on positions, orientations, sizes and other specifics of each fish in one known (simulated) school.
	########## DEPENDENCIES: ###########
	# global2car(), rotate3D(), car2sph(), integrateonsphere(), zeros()
	############ VARIABLES: ############
	# 	'Nl' is the number of fish in the school.
	#	'Ni' is the number of beams of the sonar.
	#
	# ---data--- see below.
	# ---esnm--- is the name of the echo sounder. Currently implemented are "MS70", "ME70", "EK60", "MS70_circular" (may be given in lover case)
	# ---TVG.exp--- is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
	# ---compensated--- specifies which rotations are compensated for by the echo sounder. Legal values are "pitch" and "roll", which are both compensated for in the MS70 sonar.
	# ---calibrate--- is FALSE if calibration data are to be discarded if present.
	# ---noise--- is a vector of character strings of length 2, specifying which types of noise to apply to the data:
	#		"nr" - Near-range noise.
	#		"bg" - Background noise.
	#		"ex" - Exponential noise due to acoustic interference (applies to near-range noise ("nr"), background noise ("bg") and the simulated echo).
	# ---max.memory--- is the maximum amount of memory to allow for the function to occupy. The level of for loops is chosen according to this value, i.e. if the method using only two loops, one for the radial distances and one for the unique frequencies, demands more memory than 'max.memory', the method usint three for loops is chosen.
	# ---dumpfile--- is the name of the file to which information and warnings about the simulation is written.
	# ---parlist--- is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
	# ---data--- is the list of input parameters giving the information about the school consisting of 'Nl' fish.
	#
	#		INPUTS REPRESENTING THE SCHOOL:
	#		
	#		--- Dynamic variables ---
	#		[['psxf']] - A vector of length 'Nl' representing the x-positions of the fish in the school, given in the global coordinate system (G).
	#		[['psyf']] - A vector of length 'Nl' representing the y-positions of the fish in the school, given in the global coordinate system (G).
	#		[['pszf']] - A vector of length 'Nl' representing the z-positions of the fish in the school, given in the global coordinate system (G).
	#		[['vlxf']] - A vector of length 'Nl' representing the velocity in the x-direction of the fish in the school, given in the global coordinate system (G) (optional).
	#		[['vlyf']] - A vector of length 'Nl' representing the velocity in the y-direction of the fish in the school, given in the global coordinate system (G) (optional).
	#		[['vlzf']] - A vector of length 'Nl' representing the velocity in the z-direction of the fish in the school, given in the global coordinate system (G) (optional).
	#		[['rtxf']] - A vector of length 'Nl' representing the x-rotation angles 'b' of the fish in the school, given in the global coordinate system (G). The order of rotation is set to be "z-x-y" by the angles a,b,c.
	#		[['rtyf']] - A vector of length 'Nl' representing the y-rotation angles 'c' of the fish in the school, given in the global coordinate system (G). The order of rotation is set to be "z-x-y" by the angles a,b,c.
	#		[['rtzf']] - A vector of length 'Nl' representing the z-rotation angles 'a' of the fish in the school, given in the global coordinate system (G). The order of rotation is set to be "z-x-y" by the angles a,b,c.
	#		
	#		--- Static variables ---
	#		[['indl']] - A vector of length 'Nl' representing a numbering index of the fish.
	#		[['acca']] - A vector of optimal acoustic cross sectional area A_0 of the fish in m^2.
	#		[['size']] - A vector of length 1 or 'Nl' representing the lengths of the fish in the school.
	#		[['tilt']] - A vector of length 1 or 'Nl' representing the tilt angle of the swim bladder of the fish in the school, needed if the fish have swim bladder and the beam pattern from a line source or some other parametric model is used for the beam pattern of the fish.
	#		[['epsl']] - A vector of length 1 or 'Nl' representing the factor linking optimal scattering area A0 and size squared: A0 = epsilon * S^2.
	#		[['gamw']] - A vector of length 1 or 'Nl' representing the compression factor for the hydrostatic compression of the width of the swim bladder.
	#		[['gaml']] - A vector of length 1 or 'Nl' representing the  compression factor for the hydrostatic compression of the length of the swim bladder.
	#		[['obln']] - A vector of length 1 or 'Nl' representing the oblongness of the fish, i. e. the ratio of length and radius of the rounded cylinder representing the acoustic properties of the fish.
	#		[['zeta']] - A vector of length 1 or 'Nl' representing the ratio of length of the line source representing the fish and the length of the fish.
	#		
	#		--- Beam pattern ---
	#		[['pbpf']] - The beam pattern of the fish as a function of one or two angle variables and a variable representing the relative size of the fish. Given either as a function or as the name of a predefined function (one of "lineSource" or "circularPiston").
	#		[['graf']] - A vector of arbitrary length representing the grid azimuth angle to the line sources representing the fish. Only required if empirical beam pattern as a function of two angles of direction is given.
	#		[['gref']] - A vector of arbitrary length representing the grid elevation (incidence) angle to the line sources representing the fish. Must be given along with the empirical beam pattern.
	#		[['grsf']] - A vector of arbitrary length representing the grid size of the source relative to the wavelength, represented by the product of size of the source and wave number.
	#		[['grif']] - A vector of arbitrary length representing the grid numbering index of the fish.
	#		[['dbpf']] - A vector representing the dimension of the empirical beam pattern.
	#		[['ebpf']] - An array of dimension no less than c(length(data$graf), length(data$gref), length(data$rszf)) representing the empirical beam pattern values of the fish.
	#		[['scls']] - A single numeric used to scale the backscatter from the targets. If scls>1, the CPU time can be reduced by selecting less fish. Too large value of 'scls' is not recommended.
	#
	#
	#		INPUTS REPRESENTING THE ECHO SOUNDER:
	#		
	#		--- Configuration variables ---
	#		[['indi']] - A vector of length 'Ni' representing a numbering index of the beams.
	#		[['dira']] - A vector of length 'Ni' representing the azimuth angle of the direction of the beam in the spherical coordinate system of the vessel. 
	#		[['dire']] - A vector of length 'Ni' representing the elevation angle of the direction of the beam in the spherical coordinate system of the vessel. 
	#		[['freq']] - A vector of length 1 or 'Ni' representing the frequency of the beam.
	#		[['absr']] - A vector of length 1 or 'Ni' representing the absorption coefficient of the sea for each beam.
	#		[['psze']] - A vector of length 1 representing the z-coordinate of the echo sounder in the coordinate system of the vessel.
	#		[['sint']] - A vector of length 1 or 'Ni' representing the duration of the sound pulse for each beam.
	#		[['lenb']] - A vector of length 1 or 'Ni' representing the lengths of the beams.
	#		[['rad1']] - A vector of length 1 or 'Ni' representing the radius of the transducer elements at emission.
	#		[['rad2']] - A vector of length 1 or 'Ni' representing the radius of the transducer elements at reception.
	#		[['bwtx']] - A vector of length 1 or 'Ni' representing the beam widths (-3 dB beam angles) at emission.
	#		[['bwty']] - A vector of length 1 or 'Ni' representing the beam widths (-3 dB beam angles) at reception.
	#		
	#		--- Beam pattern ar emission ---
	#		[['pbp1']] - The beam patterns of the echo sounder at emission as a function of one or two angle variables and a variable representing the relative size of the transducer element. Given either as a function or as the name of a predefined function (one of "lineSource" or "circularPiston").
	#		[['gra1']] - A vector of arbitrary length representing the grid azimuth angle to the transducer elements at emission. Only required if empirical beam pattern as a function of two angles of direction is given.
	#		[['gre1']] - A vector of arbitrary length representing the grid elevation (incidence) angle to the transducer elements at emission. Must be given along with the empirical beam pattern.
	#		[['grs1']] - A vector of arbitrary length representing the grid size of the source relative to the wavelength at emission, represented by the product of size of the source and wave number.
	#		[['gri1']] - A vector of arbitrary length representing the grid numbering index of the beams at emission.
	#		[['dbp1']] - A vector representing the dimension of the empirical beam pattern.
	#		[['ebp1']] - An array of dimension no less than c(length(data$gra1), length(data$gre1), length(data$rsz1)) representing the empirical beam pattern values of the transducer elements at emission.
	#		
	#		--- Beam pattern ar reception ---
	#		[['pbp2']] - The beam patterns of the echo sounder at reception as a function of one or two angle variables and a variable representing the relative size of the transducer element. Given either as a function or as the name of a predefined function (one of "lineSource" or "circularPiston").
	#		[['gra2']] - A vector of arbitrary length representing the grid azimuth angle to the transducer elements at reception. Only required if empirical beam pattern as a function of two angles of direction is given.
	#		[['gre2']] - A vector of arbitrary length representing the grid elevation (incidence) angle to the transducer elements at reception. Must be given along with the empirical beam pattern.
	#		[['grs2']] - A vector of arbitrary length representing the grid size of the source relative to the wavelength at reception, represented by the product of size of the source and wave number.
	#		[['gri2']] - A vector of arbitrary length representing the grid numbering index of the beams at reception.
	#		[['dbp2']] - A vector representing the dimension of the empirical beam pattern.
	#		[['ebp2']] - An array of dimension no less than c(length(data$gra2), length(data$gre2), length(data$rsz2)) representing the empirical beam pattern values of the transducer elements at reception.
	#
	#		--- Beam pattern ar reception ---
	#		[[bgns]] - Estimated background noise.
	#		[[nrns]] - Estimated near-range noise.
	#
	#
	#		INPUTS REPRESENTING THE VESSEL:
	#		
	#		[['psxv']] - A vector of length 1 representing the x-positions of the vessel, given in the global coordinate system (G).
	#		[['psyv']] - A vector of length 1 representing the y-positions of the vessel, given in the global coordinate system (G).
	#		[['pszv']] - A vector of length 1 representing the z-positions of the vessel, given in the global coordinate system (G).
	#		[['rtxv']] - A vector of length 1 representing the x-rotation angles 'b' of the vessel, given in the global coordinate system (G). The order of rotation is set to be "z-x-y" by the angles a,b,c.
	#		[['rtyv']] - A vector of length 1 representing the y-rotation angles 'c' of the vessel, given in the global coordinate system (G). The order of rotation is set to be "z-x-y" by the angles a,b,c.
	#		[['rtzv']] - A vector of length 1 representing the z-rotation angles 'a' of the vessel, given in the global coordinate system (G). The order of rotation is set to be "z-x-y" by the angles a,b,c.
	#		[['lonv']] - A vector of length 1 representing the longitude positions of the vessel given in decimal degrees east (positive or negative)s.
	#		[['latv']] - A vector of length 1 representing the latitude positions of the vessel given in decimal degrees north (positive or negative).
	#		[['orgn']] - A vector of 2 elements latitude and longitude representing the origin of the global coordinate system (G)
	#
	#
	#		INPUTS REPRESENTING THE SEA:
	#		
	#		[['rho0']] - A vector of length 1 representing the equilibrium mas density of the sea in kg/m^3.
	#		[['gacc']] - A vector of length 1 representing the gravitational constant (defaulted to 9.82, representing the North Sea).
	#		[['hpr0']] - A vector of length 1 representing the hydrostatic pressure at sea level.
	#		[['asps']] - A vector of length 1 representing the speed of sound of the sea, assumed (as a simplification) to be constant through the water column.
	# ---ask--- is TRUE if the used should be asked to for approval if the memory of the least memory demanding calculation method of the individual radial sampling intervals exceed the memory limit 'max.memory'.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Define dump data:
	dumpdata = NAs(7)
	names(dumpdata) = c("etaj", "B_L", "etaa", "epss", "epsl", "chi", "sigma0mode")
	
	##### Defaults: #####
	# School dynamic:
	default.rtxf = 0
	default.rtyf = 0
	default.size = 0.3 # Herring?
	default.zeta = 0.26 # Gorska and Ona 2003, Modelling the effect of swimbladder compression on the acoustic backscattering from herring at normal or near-normal dorsal incidences.
	# Vessel:
	default.pszv = 0
	default.rtxv = 0
	default.rtyv = 0
	
	
	##### Error and default handling: #####
	echoIBM.warnings_warninglist = NULL
	
	# Error is 'data' is not a list:
	if(!is.list(data)){
		stop("'data' not a list")
		}
	
	# Information about missing variales without defaults will be collected and printed as an error:
	errors = NULL
	
	
	### INPUTS REPRESENTING THE SCHOOL: ###
	# --- Dynamic variables ---
	# If all position variables are present, the number of fish 'Nl' is obtained as the length of the longest variable, and the shorter are recycled to match the longer. Else the number of fish 'Nl' is defaulted to 1 to complete the default actions and error handeling of missing input variables:
	if(!any(is.null(data$psxf),is.null(data$psyf),is.null(data$pszf))){
		Nl = max(length(data$psxf),length(data$psyf),length(data$pszf))
		data$psxf = rep(data$psxf,length.out=Nl)
		data$psyf = rep(data$psyf,length.out=Nl)
		data$pszf = rep(data$pszf,length.out=Nl)
		}
	else{
		errors = c(errors,"Some of 'data$posx' or 'data$psxf' (fish x-position), 'data$posy' or 'data$psyf' (fish y-position) and 'data$posz' or 'data$pszf' (fish z-position) missing with no default")
		Nl = 1
		}
	# If any fish are located above the sea surface, en error is produced:
	if(any(data$pszf>0)){
		warning(paste("Some fish are located above sea level. Range of 'pszf': (",min(data$pszf),", ",max(data$pszf),")",sep=""))
		valid = data$pszf <= 0
		correctLength = sapply(data, length) == length(valid)
		data[correctLength] = lapply(data[correctLength], "[", valid)
		#stop(paste("Some fish are located above sea level. Range of 'pszf': (",data$pszf,")",sep=""))
		#stop(paste("Some fish are located above sea level. Range of 'pszf': (",min(data$pszf),", ",max(data$pszf),")",sep=""))
		}
	# If all three velocity vectors are given, the z- and x-rotations are obtained from the velocities. Otherwise, if 'data$rtzf' (heading) is not present, an error is collected: 
	if(any(is.null(data$rtxf),is.null(data$rtzf)) && !any(is.null(data$vlxf),is.null(data$vlyf),is.null(data$vlzf))){
		# Get fish angle data if not present:
		data[c("rtzf", "rtxf")] = vl2rt.TSD(data[c("rtzf", "rtxf", "vlxf", "vlyf", "vlzf")])[c("rtzf", "rtxf")]
		# Remove velocities from memory:
		data$vlxf = NULL
		data$vlyf = NULL
		data$vlzf = NULL
		}
	else if(is.null(data$rtzf) && tolower(data$pbpf)!="ps"){
		errors = c(errors,"'data$rtzf' (fish heading) missing with no default")
		}
	# If 'data$rtxf' (pitch) is missing, it is defaulted to default.rtxf = 0 for all fish:
	if(is.null(data$rtxf)){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"'data$rtxf' was defaulted")
		data$rtxf = rep(default.rtxf,Nl)
		}
	# If the number of targets is positive and 'epss' or 'epsl' are to be utilized or 'data$lenl' (length of the line source representing the fish) are missing, 'data$size' must be present or defaulted:
	if(any(data$sigma0mode %in% c(2,4)) || length(data$lenl)==0){
		# If 'data$size' (length of fish) is missing, it is defaulted to default.size = 0.3 for all fish:
		if(length(data$size)==0){
			echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,paste("'data$size' was defaulted to ",default.size," for all fish",sep=""))
			data$size = rep(default.size,Nl)
			}
		}
	# If the size vector of the fish does not have the correct length, it is repeated to length 'Nl' with a warning:
	if(length(data$size)!=0 && length(data$size)!=Nl){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,paste("'data$size' was repeated to the number of fish ",Nl,sep=""))
		data$size = rep(data$size,length.out=Nl)
		}
	
	# If 'data$lenl' (length of the line source representing the fish) is missing, 'data$zeta' must be present or defaulted:
	if(length(data$lenl)==0){
		# If 'data$zeta' (line source length as a fraction of fish length) is missing, it is defaulted to default.zeta=1/5 for all fish:
		if(length(data$zeta)==0){
			echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,paste("'data$zeta' was defaulted to ",default.zeta,sep=""))
			data$zeta = default.zeta
			}
		}
	# Length of the line sources representing the fish:
	if(length(data$lenl)==0){
		data$lenl = data$zeta * data$size
		data$zeta = NULL
		}
	# If the length vector of the line sources does not have the correct length, it is repeated to length 'Nl' with a warning:
	if(length(data$lenl)!=0 && length(data$lenl)!=Nl){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,paste("'data$lenl' was repeated to the number of fish ",Nl,sep=""))
		data$lenl = rep(data$lenl,length.out=Nl)
		}
	
	
	### INPUTS REPRESENTING THE VESSEL: ###
	# If longitude and latitude positions and origin are all given, cartesian x- and y-positions relative to the origin are calculated:
	if(!any(is.null(data$lonv),is.null(data$latv),is.null(data$orgn))){
		data$psyv = global2car(cbind(data$lonv,data$latv),data$orgn) # This is what was former names 'pos', but changed to avoid memory leak.
		data$psxv = data$psyv$x
		data$psyv = data$psyv$y
		}
	# If one or both of data$posx or data$psxv (x-position) and data$posy or data$psyv  (y-position) are missing, an error is collected:
	if(any(is.null(data$psxv),is.null(data$psyv))){
		errors = c(errors,"One or both of 'data$posx'/'data$psxv' (x-position) and 'data$posy'/'data$psyv' (y-position) missing with no default")
		}
	# If 'data$pszv' (z-position of the vessel) is missing, it is defaulted to default.pszv=0:
	if(is.null(data$pszv)){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"'data$pszv' was defaulted")
		data$pszv = default.pszv
		}
	# If 'data$rtzv' (z-rotation of the vessel) is missing an error is collected:
	if(is.null(data$rtzv)){
		errors = c(errors,"'data$rtzv' (z-rotation of the vessel) missing with no default")
		}
	# The user may choose to ignore pitch and roll of the vessel, simulating the compensation of pitch and roll performed in Simrad systems:
	if(sum(grep("pitch",compensated))>0){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"'data$rtxv' (compensated for and set to 0)")
		data$rtxv = default.rtxv
		}
	if(sum(grep("roll",compensated))>0){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"'data$rtyv' (compensated for and set to 0)")
		data$rtyv = default.rtyv
		}
	# If 'data$rtxv' (x-rotation of the vessel) is missing, it is defaulted to default.rtxv=0:
	if(is.null(data$rtxv)){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"'data$rtxv' was defaulted")
		data$rtxv = default.rtxv
		}
	# If 'data$rtyv' (t-rotation of the vessel) is missing, it is defaulted to default.rtyv=0:
	if(is.null(data$rtyv)){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"'data$rtyv' was defaulted")
		data$rtyv = default.rtyv
		}
	########## Execution and output ##########
	##### Transformations and rotations: #####
	
	### Position of the transducer in the coordinate systems of the fish: ###
	cat(paste("Transformations and rotations (number of targets: ",length(data$psxf),") ... \n",sep=""))
	# Add the position of the echo sounder:
	transducerposG = c(data$psxv,data$psyv,data$pszv)+c(0,0,data$psze[1])
	
	# Subtract the positions of the fish:
	transducerposL = matrix(transducerposG,ncol=3,nrow=Nl,byrow=TRUE)-cbind(data$psxf,data$psyf,data$pszf)
	
	# Discard fish that are outside of the range of the sonar:
	maxlenb = max(data$lenb)
	# The radial resolution 'data$rres':
	data$rres = soundbeam_range(data, "rres")
	#data$rres = data$sint * data$asps/2
	sonarRange = maxlenb * data$rres
	indicesRange = rowSums(transducerposL^2)<sonarRange^2
	Nl = sum(indicesRange)
	
	if(Nl>0){
		transducerposL = transducerposL[indicesRange,,drop=FALSE]
		# Rotate as given in the documentation:
		if(is.null(data$rtyf)){
			data$transducerposL = rotate3D(transducerposL, by="zx", ang=cbind(data$rtzf,data$rtxf+data$tilt-pi/2), paired=TRUE, drop.out=FALSE)
			}
		else{
			data$transducerposL = rotate3D(transducerposL, by="zxyx", ang=cbind(data$rtzf,data$rtxf,data$rtyf,data$tilt-pi/2), paired=TRUE, drop.out=FALSE)
			data$rtyf = NULL
			}
		# Remove objects no longer in use:
		rm(transducerposG)
		rm(transducerposL)
		data$rtxf = NULL
		data$rtzf = NULL
		data$tilt = NULL
		#gc()
		
		# Transform to spherical cordinates:
		data$transducerposL = car2sph(data$transducerposL)
		
		
		##### Outside the loop: #####
		cat("Frequency independent calculations ... \n")
		
		data$fish = zeros(Nl)
		# Orientation factor:
		#data$etaomega = etaOmega(data)
		data$etaomega = etaOmega(data[c("pbpf", "obln", "transducerposL", "rho0", "gacc", "pszf", "hpr0", "gaml", "gamw")])
		
		# Depth compression factor 'etaC':
		data$etaC = etaCompression(data[c("rho0", "gacc", "pszf", "hpr0", "gaml", "gamw")],type="area")
		#data$etaC = etaCompression(data,type="area")
		#etaC = (1-data$rho0 * data$gacc * data$pszf/data$hpr0)^(data$gamw+data$gaml)
		# Attenuation:
		data$etar4 = data$transducerposL[,1]^(-4)
		
		# If the optimal backscattering cross section 'sgbs' (sigma_bs) is present, no relation to target size is specified. In this case 'sgbs' is applied at this stage and the echo ability 'epss' is NOT USED in the function echoIBM.oneping.ji() where i = 0, 1, 2:
		if(data$sigma0mode==1){
			# data$sgbs already present:	
			}
		# If the coefficient 'epss' linking 'sgbs' to fish size is present, it may be a function of frequency, or simply a numeric vector. In any case the target size to the power 'spow' is applied at this stage, and 'epss' is applied in the function echoIBM.oneping.ji() where i = 0, 1, 2:
		else if(data$sigma0mode==2){
			data$sgbs = data$size^data$spow
			}
		# If the optimal acoustic cross sectional area 'acca' (A_0) is present, no relation to target size is specified: In this case 'acca' is applied at this stage and the echo ability 'epsl' is NOT USED in the function echoIBM.oneping.ji() where i = 0, 1, 2:
		else if(data$sigma0mode==3){
			data$sgbs = data$acca
			}
		# If the coefficient 'epsl' linking 'acca' to fish size is present, it may be a function of frequency, or simply a numeric vector: In any case the target size to the power 'spow' is applied at this stage, and 'epsl' is applied in the function echoIBM.oneping.ji() where i = 0, 1, 2:
		else if(data$sigma0mode==4){
			data$sgbs = data$size^data$spow
			}
			
		# Merge non-frequency independent values:	
		data$fish = data$etar4 * data$etaomega * data$etaC * data$sgbs
			
		# Scale the backscatter of the targets (scls>1 reduces CPU time):
		if(length(data$scls)==1 && data$scls!=1){
			data$fish = data$fish * data$scls
			}
		
		# Dump:
		if(length(dumpfile)>0 && nchar(dumpfile)>0){
			# Dynamic information about the school:
			echoIBM.dump_summary(data[labl.TSD("v")], dumpfile, type="vessel", append=TRUE)
			# Dynamic information about the vessel:
			echoIBM.dump_summary(data[labl.TSD("dse")], dumpfile, type="dynschool", append=TRUE)
			}
			
		# Clear from memory:
		data$size = NULL
		data$etaC = NULL
		data$etar4 = NULL
		data$etaomega = NULL
		# gc()
		
		### Defining indexes for which radial paritions the fish are located in:
		cat("Radial paritioning of targets ... \n")
		# Indexes locating the fish at the correct radial layer, assuming constant speed of sound:
		##### HERE IT WAS DISCOVERED A CHANGE IN THE EFFECT OF THE OCDE SINCE THE CHANGE AT 2014-03-14, AFTER WHICH THE ECHO FROM TARGETS WERE NO LONGER CALCULATED TWICE, BUT ONLY ONCE, AND THE RADIAL WHEITHING ETA_J WAS APPLIED TWICE INSTEAD, ONCE USING ETA_J AND ONCE USING 1-ETA_J. HOWEVER, THIS REVEALED THE FACT THAT data$radialindex DOES NOT REPRESENT THE VOXELS BUT RATHER THE BOUNDARIES OF THE RADIAL POSITIONS COVERED BY THE VOXELS. THUS THE CODE IN THE CALIBRATION ROUTINES SHOULD BE CHANGED ACCORDINGLY
		data$radialindex = findInterval(data$transducerposL[,1]/data$rres,0:maxlenb)
		# An object contributes to the time interval located above and the next (if r/data$rres=0.8, it contributes to time interval 1 and 2). 'data$validr' is a vector of the time intervals that are affected by the fish:
		data$validr = min(data$radialindex):(max(data$radialindex)+1)
		# Adjust 'data$validr' and 'data$radialindex' to the length of the beams:
		data$validr = data$validr[data$validr<=maxlenb]
		# Divide 'data$radialindex' into a list of indexes for the radial layers/shells. 'NULLlist' assists in filling in NULL for layers with no fish:
		data$radialindex = split(order(data$radialindex),sort(data$radialindex))
		# Adjust 'data$validr' and 'data$radialindex' to the length of the beams:
		data$radialindex = data$radialindex[as.numeric(names(data$radialindex))<=maxlenb]
		
		# Insert null elements in the list of radial indices, so that voxels with targets only in the first half of the valid radial region of the voxel are counted as well:
		namesradialindex = as.numeric(names(data$radialindex))
		namesNULLlist = seq(min(namesradialindex),max(namesradialindex)+1)
		NULLlist = vector("list",length(namesNULLlist))
		names(NULLlist) = namesNULLlist
		NULLlist[as.character(namesradialindex)] = data$radialindex
		data$radialindex = NULLlist
		
		# update 'Nl':
		Nl = sum(unlist(sapply(data$radialindex,length)))
		}
	else{
		data$validr = NULL
		}
	
	# Output voxel system:
	sv = zeros(maxlenb,data$Ni)
	nsig = zeros(maxlenb,data$Ni)
	
	#memory_basic = 8 * maxlenb * data$numb
	memory_basic = object.size(data)+object.size(sv)
	# A vector of the echos from the last radial layer:
	fromlast = 0
	
	##### The loop through time intervals/radial distances #####
	cat("r[j]:\n")
	write("\nr[n]:\n",dumpfile,append=TRUE)
	for(j in seq_along(data$validr)){
		
		# Which of the fish that are to be treated. According to expression (3.9) in the documentation, the fish of radial distance in the range (j-2) * delta_r to j * delta_r will contribute to the j'th time interval, which is realized by data$radialindex[j+0:1]:
		data$thesel = data$radialindex[[j]]
		data$lthesel = length(data$thesel)
		if(data$lthesel>0 || fromlast>0){
		
			# Input only a subset of the data to the functions echoIBM.oneping.j1() and echoIBM.oneping.j2():
			thesedata = data[c("psxf", "psyf", "pszf", "epss", "numb", "indi", "dira", "dire", "absr", "psze", "asps", "grsf", "ebpf", "pbpf", "ssif", "sigma0mode", "sllf", "pbp1", "pbp2", "bptf", "bpt1", "bpt2", "epsl", "wavenumber", "uniquek", "kfreq", "rad1", "rad2", "equalbp_em_re", "psxv", "psyv", "pszv", "rtzv", "lenl", "rtxv", "rtyv", "transducerposL", "fish", "rres", "validr", "thesel", "lthesel")]
			
			memory_1 = 8 * data$lthesel * (5+length(data$wavenumber)/data$luniquek*7) + memory_basic
			memory_2 = 8 * data$lthesel * 12 + memory_basic
			
			# Empirically found adjustment factors for the different architectures:
			if(.Platform$r_arch=="x86_64"){
				memory_1 = memory_1*4
				memory_2 = memory_2*4
				}
			else if(.Platform$r_arch=="i386"){
				memory_1 = memory_1*2
				memory_2 = memory_2*2
				}
			
			if(memory_1<max.memory){
				# Print the currently processed radial layer:
				cat(data$validr[j],"[",data$lthesel,"]\t",sep="")
				cat(data$validr[j],"[",data$lthesel,"]\t",sep="",file=dumpfile,append=TRUE)
				
				# Get the svs from this radial layer, and those extending into the next:
				thissv = echoIBM.oneping.j1(j,thesedata)
				fromthis = thissv$sv
				sv[data$validr[j],] = fromlast + fromthis
				fromlast = thissv$svnext
				
				dumpdata = rowSums(cbind(dumpdata, thissv$dumpdata[names(dumpdata)] * data$lthesel), na.rm=TRUE)
				
				nsig[data$validr[j],] = thissv$nsig
				}
			else if(memory_2<max.memory){
				# Print the currently processed radial layer:
				cat(data$validr[j],"*[",data$lthesel,"]\t",sep="")
				cat(data$validr[j],"*[",data$lthesel,"]\t",sep="",file=dumpfile,append=TRUE)
				
				# Get the svs from this radial layer, and those extending into the next:
				thissv = echoIBM.oneping.j2(j,thesedata)
				fromthis = thissv$sv
				sv[data$validr[j],] = fromlast + fromthis
				fromlast = thissv$svnext
				
				dumpdata = rowSums(cbind(dumpdata, thissv$dumpdata[names(dumpdata)] * data$lthesel), na.rm=TRUE)
				nsig[data$validr[j],] = thissv$nsig
				}
			else if(ask){
				ans = readline(paste("Memory (",max.memory,") limit exceeded (",memory_2,") for the least memory demanding method. Continue? (y/n)",sep=""))
				if(ans=="n"){
					stop("Function terminated")
					}
				# Print the currently processed radial layer:
				cat(data$validr[j],"*[",data$lthesel,"]\t",sep="")
				cat(data$validr[j],"*[",data$lthesel,"]\t",sep="",file=dumpfile,append=TRUE)
				
				# Get the svs from this radial layer, and those extending into the next:
				thissv = echoIBM.oneping.j2(j,thesedata)
				fromthis = thissv$sv
				sv[data$validr[j],] = fromlast + fromthis
				fromlast = thissv$svnext
				
				dumpdata = rowSums(cbind(dumpdata, thissv$dumpdata[names(dumpdata)] * data$lthesel), na.rm=TRUE)
				nsig[data$validr[j],] = thissv$nsig
				}
			}
		} # End of for j.
	cat("\n")
	cat("\n",file=dumpfile,append=TRUE)
	
	# Apply the threshold for the number of significant scatterers:
	if(length(parlist$nsth)>0){
		parlist$Brkt = which(1<=nsig & nsig<parlist$nsth)
		parlist$nsig = nsig[parlist$Brkt]
		}
	rm(nsig)
	# Transform 'dumpdata' into a list after dividing by the number of fish to get the mean:
	echoIBM.dump_summary(c(as.list(dumpdata/(Nl)),list(Nl=Nl)), dumpfile, type="freqschool", append=TRUE)
	
	# Calibrate the beams:
	if(!is.null(data$cali) && calibrate){
		# Extract for the correct beam mode, which is only applied when the calibration data are stored in a list:
		if(is.list(data$cali)){
			data$cali = data$cali[[data$bmmd[1]+1]]
			}
		# Extract calibration values at the elevation angles closest to elevation angles in the sonar
		if(length(data$grde)>0){
			atgrde = apply(abs(outer(data$dire, data$grde, "-")), 1, which.min)
			cali = data$cali[cbind(seq_along(data$dire),atgrde)]
		}
		else if(ncol(sv)!=length(data$cali)){
			stop("Length of calibration data differ from the number of beams")
			}
		else{
			cali = data$cali
			}
		write("\nCalibration factors:\n",dumpfile,append=TRUE)
		write(cali,dumpfile,append=TRUE)
		sv = sv * outer(ones(maxlenb),cali)
		cat("Simulation calibrated", "\n")
		}
	else if(is.null(data$cali) && calibrate){
		echoIBM.warnings_warninglist = c(echoIBM.warnings_warninglist,"Calibration data missing")
		}
		
	# Add noise to the sv-values (No time step index applied, To apply time step dependetnt noise, use echoIBM.merge()):
	sv = echoIBM.add.noise(sv=sv, noise=noise, data=data, parlist=parlist)
		
	# Add TVG if required:
	if(sum(TVG.exp)>0){
		sv = apply.TVG(x=sv, beams=data[c("lenb", "sint", "absr", "asps")], TVG.exp=TVG.exp)
		}
	
	
	# Print the warnings to the dumpfile:
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		write("\n\n\n##### WARNINGS IN \"echoIBM.oneping.oneschool\" FOR THE CURRENT TIME STEP: #####",dumpfile,append=TRUE)
		if(length(echoIBM.warnings_warninglist)>0){
			for(i in seq_along(echoIBM.warnings_warninglist)){
				write(paste(i,": ",echoIBM.warnings_warninglist[i],sep=""),dumpfile,append=TRUE)
				}
			}
		else{
			write("none",dumpfile,append=TRUE)
			}
		}
	
			
	# Return:
	list(sv=sv,Brkt=parlist$Brkt,nsig=parlist$nsig)
	##################################################
	##################################################
	}
