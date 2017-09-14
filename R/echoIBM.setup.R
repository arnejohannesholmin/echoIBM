#*********************************************
#*********************************************
#' ***DECRIPTION MISSING***
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR read.event
#' @importFrom TSD car2global catchWarnings ftim2utim NAs ones write.TSD zeros
#' @importFrom stats runif
#'
#' @export
#' @rdname echoIBM_SX90_BlindZone
#'
echoIBM.setup <- function(

	eventName,
	dir, 
	t = "all",
	pre = TRUE, 



	vessel = list(), 
	schoolStatic = list(), 
	schoolDynamic = list(), 
	ctd = list(), 
	beams = list(), 
	noise = list(), 
	calibration = list(), 
	
	
	addnoise = TRUE, 
	passive = FALSE, 
	cores = 1
	


   ### 
   ### 
   ### eventName,
   ### dir, 
   ### esnm, 
   ### starttime = ftim2utim(20150101),  # Start time of the event
   ### seed = 1, 
   ### lon0 = 5,  # Longitude of the start of the event
   ### lat0 = 70,  # Latitude of the start of the event
   ### beamPatternFish = FALSE, 
   ### pingduration = 5, # Seconds
   ### speedOfVessel = 10, # Knots
   ### duration = 12 * 3600, 
   ### schooDist = 1000, 
   ### maxRange = 450, 
   ### maxSize = 100, 
   ### maxHeight = 10, 
   ### minDepth = -150, 
   ### sonar_tilt = 3, 
   ### sonar_plsl = 0.006, 
   ### cores = 1, 
   ### t = "all",
   ### pre = TRUE, 
   ### schools = list(), 
   ### CTD = FALSE, 
   ### addnoise = TRUE
	){






		# 0. Create the events
		events <- createEvent(eventName, dir, esnm)
		files <- list()
		
		
		# 1. Vessel
		if(!identical(vessel, FALSE)){
			globalVar <- list(events=events, eventName=eventName)
			vessel <- do.call(echoIBM.setVessel, c(globalVar, vessel))
			files$vessel <- vessel$files # This was a nice palindrome...
		}
		
		# 2. School
		if(!identical(school, FALSE)){
			globalVar <- list(events=events, eventName=eventName, vessel=vessel)
			school <- do.call(echoIBM.setVessel, c(globalVar, school))
			files$school <- school$files
		}
		
		# 3. CTD
		
		
		# 4. Beams

		# 5. Noise

		# 6. Calibration





	
	##############################
	if(length(events)==0){
		events = getEvents(eventspec, dir, project, esnm, new=TRUE)
		}
	##############################
	############################################################

	if(pre==TRUE || pre=="only"){
		############################################################
		########### (3) COMMON FOR ALL ACOUSTIC SYSTEMS: ###########
		############################################################

		
		##############################
		########## (3a) CTD: #########
		##############################
		# Copy from an existig CTD file:
		if(CTD){
			CTDFile = "/Volumes/Acoustics/S2014119_PG.O.Sars[4174]/Events/S2014119_D200141030_E0021/SX90/tsd/S2014119_D200141030_E0021_SX90.ctd"
		for(i in seq_along(events)){
				file.copy(CTDFile, events[i], overwrite=TRUE)
				}
			}
		
		### Or generate the CTD file here:
		#CTD = list()
		#CTD$utim = starttime
		#CTD$lonc = lon0
		#CTD$latc = lat0
		#CTD$pszc = # A vector of depths (negative values)
		#CTD$ihpr = # A vector of hyrdostatic pressure, given in decibar (10 000 Pa, approximately equivalent to meters in depth)
		#CTD$temp = # A vector of temperatures in celcius
		#CTD$slty = # A vector of salinity in parts pr million
		#CTD$isps = # A vector of sound speed (if already calculated)
		#CTD$rho0 = # A vector of mass density of the sea
		#CTD$gacc = # The gravitational acceleration constant
		#CTD$hpr0 = # Hydrostatic pressure at the surface in Pascal

		# Write the CTD file:
		#CTDFiles = file.path(events, paste0(project, "_", eventspec, ".ctd"))
		#for(i in seq_along(events)){
		#	write.TSD(CTD, CTDFiles[i], numt=1)
		#	}
		##############################

		##############################
		## (3b) TARGET BEAMPATTERN: ##
		##############################
		if(beamPatternFish){
			targetBeamPatternFile = "/Users/arnejh/Data/echoIBM/Resources/TargetBeamPattern/obln5/psms_gasfilled_obln_5_phi_0-0.5-90_kL_0.2-1-47.2__kirk_gasfilled_obln_5_phi_0-0.5-90_kL_48.2-1-300.2_ssif.tsd"
			for(i in seq_along(events)){
				file.copy(targetBeamPatternFile, events[i], overwrite=TRUE)
				}
		}
		##############################

		# !!! NOTE THAT THE VESSEL FILE IS USUALLY DIFFERENT BETWEEN ACOUSTIC SYSTEMS IN REAL DATA, DUE TO DIFFERENT PING RATE AND THAT ONE DOES NOT PING AT THE SAME TIME. IN SIMULATIONS THERE IS HOWEVER NO APPARENT NEED FOR USING DIFFERENT VESSEL FILES, UNLESS DIFFERENT PING RATE IS NEEDED !!!

		##############################
		######## (3c) VESSEL: ########
		##############################
		# Basic variables:
		speedInMetersPerSecond = speedOfVessel * 0.514444 # Convert to meters per second

		# Set the vessel path to move from west to east during 12 hours:
		vessel = list()
		vessel$utim = seq(0, duration, pingduration) + starttime
		numt = length(vessel$utim)

		vessel$psxv = speedInMetersPerSecond * (vessel$utim - vessel$utim[1])
		vessel$psyv = zeros(numt)
		vessel$pszv = zeros(numt)
		vessel$ispv = rep(speedOfVessel, numt)
		vessel$rtxv = zeros(numt)
		vessel$rtyv = zeros(numt)
		vessel$rtzv = rep(-pi/2, numt)
		vessel$lon0 = rep(lon0, numt)
		vessel$lat0 = rep(lat0, numt)
		# Convert to latitude and longitude values, and remove the x and y:
		longlat = car2global(cbind(vessel$psxv, vessel$psyv), origin=c(vessel$lon0[1], vessel$lat0[1], 0))
		vessel$lonv = longlat[,1]
		vessel$latv = longlat[,2]

		# Write the vessel data:
		vesselFiles = file.path(events, paste0(project, "_", eventspec, ".vessel"))
		for(i in seq_along(events)){
			write.TSD(vessel, vesselFiles[i], numt=numt)
			}
		# vv = read.TSD(vesselFiles[1], t="all")
		##############################

		##############################
		#### (3d) SCHOOLS STATIC: ####
		##############################
		#schoolStaticFile = 
		#for(i in seq_along(events)){
		#	file.copy(schoolStaticFile, events[i], overwrite=TRUE)
		#	}
	
		### Or generate the static school file for herring here:
		schoolStatic = list()
		# Tilt angle of the swim bladders:
		schoolStatic$tilt=0
		# Assuming no significant compression of the length and based on Ona 2003:
			#schoolStatic$gamw=-0.23
		schoolStatic$gamw=0
		schoolStatic$gaml=0
		schoolStatic$obln=5 # Oblongness of the swim bladder model (length / width)
		schoolStatic$zeta=0.26 # Gorska and Ona 2003, Modelling the effect of swimbladder compression on the acoustic backscattering from herring at normal or near-normal dorsal incidences.
		if(beamPatternFish){
			#schoolStatic$pbpf="ps"
		}
		else{
			schoolStatic$pbpf="ps"
		}
		schoolStatic$epss=function(f){
			10^(-6.54)*100^2 * (f/38000)^(-0.4) # (Ona 2003, and the data collected in the years 1996 - 2010, given here in units of m, and not cm as in Ona 2003 (thus the 100^2))
			}

		# Write the static school file to the first event (the first acoustic system):
		schoolStaticFiles = file.path(events, paste0(project, "_", eventspec, "_herring_static.school"))
		
		# Add any static school variables given in 'schools':
		schoolStaticNames <- labl.TSD("ss")
		schoolStaticNamesPresent <- intersect(schoolStaticNames, names(schools))
		
		schoolStatic[schoolStaticNamesPresent] <- schools[schoolStaticNamesPresent]
		#schoolStaticNames <- intersect(names(schoolStatic), names(schools))
		#schoolStatic[commonNames] <- schools[commonNames]
		
		write.TSD(schoolStatic, schoolStaticFiles[1], numt=1)
		##############################


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
			# (VIIIc) Seed for the fish sizes (sEds)
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
			## (VIIIc) sEds: ##
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
		##############################
		#Ttha = writeTtha(events[1])

		############################################################
		############################################################
		############################################################


		############################################################
		########## (4) SPECIFIC FOR EACH ACOUSTIC SYSTEM: ##########
		############################################################

		##############################
		## (4c) BEAM CONFIGURATION: ##
		##############################
		# Critical settings, must be defeined?
		# 1. Set esnm
		# 2. Set utim
		# 3. Set beam mode
		# 4. Set tilt for omni
		# 5. Set vertical fan orientation for omni
		# 




		beamsFiles = file.path(events, paste0(project, "_", eventspec, "_", esnm, ".beams"))
		resourcesBeams = "~/Data/echoIBM/Resources/Beams"
		# !!! Use numt=1 since we only simulated the horizontal mode. This saves time when reading the beams-file, a problem that must be fixed: !!!
	
		###########
		## SU90: ##
		###########
		if("SU90" %in% esnm){
			manipulateBeams(
				reference_beamsFile=file.path(resourcesBeams, "S2015116_PG.O.Sars_E0001_SU90.beams"), 
				outfile=beamsFiles["SU90"==esnm], 
				numt=1, 
				tilt=sonar_tilt, 
				utim=vessel$utim, 
				plsl=sonar_plsl, 
				maxRange=maxRange)
			}
		###########
						
		###########
		## EK60: ##
		###########
		if("EK60" %in% esnm){
			manipulateBeams(reference_beamsFile=file.path(resourcesBeams, "S2014119_D200141030_E0021_EK60.beams"), 
			outfile=beamsFiles["EK60"==esnm], 
			numt=1, 
			utim=vessel$utim[1], 
			maxRange=maxRange*2)
			}
		###########


		##############################
		###### (4d) CALIBRATION: #####
		##############################
		calidir="/Users/arnejh/Data/echoIBM/Resources/Calibration"
		dire=(seq(90,135,5))*pi/180

		###########
		## SU90: ##
		###########
		#set.seed(90)
		#system.time(warn<-catchWarnings(cal<-echoIBM.calibrate(directory=calidir, event=events[1], cruise="SU90_BlindZone_freq_26000", N=1e5, runs=50, esnm="SU90", type="sv", dire=dire, usemean=TRUE)))
		# Imediately add the system warnings to the dumpfile:
		#echoIBM.addwarnings(cal$caldir,warn)
		#     user   system  elapsed 
		# 2848.324  699.516 3524.614 
		# The calibration factor:
		#cal$cali
		# Copy the calibration file to the event:
		if("SU90" %in% esnm){
			file.copy(file.path(calidir, "CalibrationTable_SU90_Cruise_SU90_BlindZone_freq_26000_sv.tsd"), events["SU90"==esnm], overwrite=TRUE)
			}
		###########

		###########
		## EK60: ##
		###########
		#set.seed(60)
		#system.time(warn<-catchWarnings(cal<-echoIBM.calibrate(directory=calidir, event=events[2], cruise="EK60_BlindZone", N=1e5, runs=50, esnm="EK60", type="sv", usemean=TRUE, max.radius= 0.4)))
		# Imediately add the system warnings to the dumpfile:
		#echoIBM.addwarnings(cal$caldir,warn)
		#   user  system elapsed 
		# 46.788  11.517  60.613
		# The calibration factor:
		#cal$cali
		# Copy the calibration file to the event:
		if("EK60" %in% esnm){
			file.copy(file.path(calidir, "CalibrationTable_EK60_Cruise_EK60_BlindZone_sv.tsd"), events["EK60"==esnm], overwrite=TRUE)
			}
		###########


		##############################
		######### (4e) NOISE: ########
		##############################
		noiseFiles = file.path(events, paste0(project, "_", eventspec, "_", esnm, "_noise.tsd"))
		resourcesNoise = "~/Data/echoIBM/Resources/Noise"
		
		###########
		## SU90: ##
		###########
		if(addnoise && "SU90" %in% esnm){
			beams = read.event(event=events["SU90"==esnm], var="beams")
			noise = list(nr0a = SX90_nr0a(beams, SNref=-65, rref=600, SN0=-25, rjoin=20, m=0, TVG.exp=2))
			write.TSD(noise, noiseFiles["SU90"==esnm])
			}
		###########

		###########
		## EK60: ##
		###########
		if(addnoise && "EK60" %in% esnm){
			file.copy(file.path(resourcesNoise, "S2009116_Noise_EK60_bgns45_pdns45_quantile02.tsd"), events["EK60"==esnm], overwrite=TRUE)
			}
		###########
	
	if(pre=="only"){
		return(TRUE)
		}
	}

	
	##################################
	######### (5) Simulation: ########
	##################################
	# Create dumpfile directories:
	dumpfiledir = file.path(events, "dumpfiles")
	for(i in seq_along(dumpfiledir)){
		suppressWarnings(dir.create(dumpfiledir[i]))
		}
	dumpfilename = file.path(dumpfiledir,"dump.txt")

	
	# Add noise if specified:
	noise <- if(addnoise) c("ms","ca","bg") else "ms"
	
	
	###########
	## SU90: ##
	###########
	if("SU90" %in% esnm){
		# Simulatte in parallel:
		event = c(events["SU90"==esnm], if(length(events)>1 && "SU90"!=esnm[1]) events[1])
		system.time(warn <- catchWarnings(echoIBM(event=event, t=t, noise=noise, max.memory=2e9, dumpfile=dumpfilename["SU90"==esnm], parlist=list(pre=TRUE, seed="ind"), rand.sel=10, cores=cores)))
		#    user   system  elapsed 
		# 216.750   76.885 4722.775 
		# Only warnings related to not writing variables with not extactly four characters in the variable names:
		echoIBM.addwarnings(dumpfilename["SU90"==esnm], warn)
		}
	###########

	###########
	## EK60: ##
	###########
	if("EK60" %in% esnm){
		# Simulatte in parallel:
		event = c(events["EK60"==esnm], if(length(events)>1 && "EK60"!=esnm[1]) events[1])
		system.time(warn <- catchWarnings(echoIBM(event, t=t, noise=noise, max.memory=2e9, dumpfile=dumpfilename["EK60"==esnm], parlist=list(pre=TRUE, seed="ind"), max.radius=0.4, cores=cores)))
		#   user  system elapsed 
		# 44.517   7.256 996.010 
		# Only warnings related to not writing variables with not extactly four characters in the variable names:
		echoIBM.addwarnings(dumpfilename["EK60"==esnm], warn)
		}
	###########
	events
	}









createEvent <- function(eventName, dir, esnm){
	eventPath <- file.path(dir, eventName, esnm, "tsd")
	dir.create(eventPath)
}