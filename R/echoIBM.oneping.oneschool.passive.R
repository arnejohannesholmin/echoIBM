#*********************************************
#*********************************************
#' Simulates one echo sounder observation of based on positions, orientations, sizes and other specifics of each fish in one known (simulated) school.
#'
#' @param data  see below.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#' @param compensated  specifies which rotations are compensated for by the echo sounder. Legal values are "pitch" and "roll", which are both compensated for in the MS70 sonar.
#' @param calibrate  is FALSE if calibration data are to be discarded if present.
#' @param noise  is a vector of character strings of length 2, specifying which types of noise to apply to the data:
#' @param dumpfile  is the name of the file to which information and warnings about the simulation is written.
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#' @param data  is the list of input parameters giving the information about the school consisting of 'Nl' fish.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw apply.TVG
#' @importFrom TSD global2car zeros
#'
#' @export
#' @rdname echoIBM.oneping.oneschool.passive
#'
echoIBM.oneping.oneschool.passive<-function(data, esnm, TVG.exp=2, compensated=c("pitch","roll"), calibrate=TRUE, noise=c("nr","bg","ex"), dumpfile="dump.txt", parlist=list()){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-09-16 - First version.
	# Update: 2010-03-22 - Changed method to reduce memory occupied by the function.
	# Update: 2011-02-20 - Changed to use beamPattern.TSD() which takes as input a list of inputs named by the convencional TSD-names. Thus the parameter "mod" as a list element of the beam pattern function is no longer needed.
	# Update: 2011-03-17 - Corrected the calculation of ciston radii, based on the new information about the beam widths in the raw files for the MS70. EK60 and ME70 not yet altered.
	# Last: 2011-08-24 - Changed the handeling of the dumpfile from the old method of writing the dumpfile and deleting at the end of the function, to not writing, in the case that dumpfile has length 0 or 0 characters. This was forced through by the fact that echoIBM.calibrate() uses echoIBM.oneschool.oneping(), which wrote to the dumpfile unintendedly.
	########### DESCRIPTION: ###########
	# Simulates one echo sounder observation of based on positions, orientations, sizes and other specifics of each fish in one known (simulated) school.
	########## DEPENDENCIES: ###########
	# global2car(), car2sph(), integrateonsphere(), zeros()
	############ VARIABLES: ############
	# 	'Nl' is the number of fish in the school.
	#	'Ni' is the number of beams of the sonar.
	#
	# ---data--- see below.
	# ---TVG.exp--- is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
	# ---compensated--- specifies which rotations are compensated for by the echo sounder. Legal values are "pitch" and "roll", which are both compensated for in the MS70 sonar.
	# ---calibrate--- is FALSE if calibration data are to be discarded if present.
	# ---noise--- is a vector of character strings of length 2, specifying which types of noise to apply to the data:
	#		"nr" - Near-range noise.
	#		"bg" - Background noise.
	#		"ex" - Exponential noise due to acoustic interference (applies to near-range noise ("nr"), background noise ("bg") and the simulated echo).
	# ---dumpfile--- is the name of the file to which information and warnings about the simulation is written.
	# ---parlist--- is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
	# ---data--- is the list of input parameters giving the information about the school consisting of 'Nl' fish.
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
	
	
	##################################################
	##################################################
	########## Preparation ##########
	##### Error and default handling: #####
	echoIBM.warnings_warninglist=NULL
	
	# Error is 'data' is not a list:
	if(!is.list(data)){
		stop("'data' not a list")
		}
	
	# Information about missing variales without defaults will be collected and printed as an error:
	errors=NULL
	
	
	##### Defaults: #####
	# Vessel:
	default.pszv=0
	default.rtxv=0
	default.rtyv=0
	
	
	### INPUTS REPRESENTING THE VESSEL: ###
	# If longitude and latitude positions and origin are all given, cartesian x- and y-positions relative to the origin are calculated:
	if(!any(is.null(data$lonv),is.null(data$latv),is.null(data$orgn))){
		data$psyv=global2car(cbind(data$lonv,data$latv),data$orgn) # This is what was former names 'pos', but changed to avoid memory leak.
		data$psxv=data$psyv$x
		data$psyv=data$psyv$y
		}
	# If one or both of data$posx or data$psxv (x-position) and data$posy or data$psyv  (y-position) are missing, an error is collected:
	if(any(is.null(data$psxv),is.null(data$psyv))){
		errors=c(errors,"One or both of 'data$posx'/'data$psxv' (x-position) and 'data$posy'/'data$psyv' (y-position) missing with no default")
		}
	# If 'data$pszv' (z-position of the vessel) is missing, it is defaulted to default.pszv=0:
	if(is.null(data$pszv)){
		echoIBM.warnings_warninglist=c(echoIBM.warnings_warninglist,"'data$pszv' was defaulted")
		data$pszv=default.pszv
		}
	# If 'data$rtzv' (z-rotation of the vessel) is missing an error is collected:
	if(is.null(data$rtzv)){
		errors=c(errors,"'data$rtzv' (z-rotation of the vessel) missing with no default")
		}
	# The user may choose to ignore pitch and roll of the vessel, simulating the compensation of pitch and roll performed in Simrad systems:
	if(sum(grep("pitch",compensated))>0){
		echoIBM.warnings_warninglist=c(echoIBM.warnings_warninglist,"'data$rtxv' (compensated for and set to 0)")
		data$rtxv=default.rtxv
		}
	if(sum(grep("roll",compensated))>0){
		echoIBM.warnings_warninglist=c(echoIBM.warnings_warninglist,"'data$rtyv' (compensated for and set to 0)")
		data$rtyv=default.rtyv
		}
	# If 'data$rtxv' (x-rotation of the vessel) is missing, it is defaulted to default.rtxv=0:
	if(is.null(data$rtxv)){
		echoIBM.warnings_warninglist=c(echoIBM.warnings_warninglist,"'data$rtxv' was defaulted")
		data$rtxv=default.rtxv
		}
	# If 'data$rtyv' (t-rotation of the vessel) is missing, it is defaulted to default.rtyv=0:
	if(is.null(data$rtyv)){
		echoIBM.warnings_warninglist=c(echoIBM.warnings_warninglist,"'data$rtyv' was defaulted")
		data$rtyv=default.rtyv
		}
		
	########## Execution and output ##########
	# Dump:
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		# Dynamic information about the school:
		echoIBM.dump_summary(data, dumpfile, type="vessel", append=TRUE)
		}
		
	# Output voxel system:
	sv=zeros(data$lenb[1],data$Ni)
	
	# Add noise to the sv-values:
	sv=echoIBM.add.noise(sv=sv, noise=noise, data=data, parlist=parlist)
		
	# Add TVG if required:
	if(sum(TVG.exp)>0){
		sv=apply.TVG(x=sv, beams=data[c("lenb","sint","absr","asps")], TVG.exp=TVG.exp)
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
	list(sv=sv)
	##################################################
	##################################################
	}
