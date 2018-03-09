#*********************************************
#*********************************************
#' Generates and writes to file, fish positions which are AR(1) but not interacting, by the model in Holmin et al. 2012, "Simulations of multi-beam sonar echos from schooling individual fish in a quiet environment".
#'
#' @param con  is a the path to the directory in which to put the files.
#' @param size  is a three element vector giving the x, y, and z size of the school (the axes, not the semi-axes).
#' @param shape  is the shape of the school, either given as one of the two pre-defined shape strings "e" (ellipsoid) and "cc" (cubic compound school structure), as a three column matrix of initial fish positions relative to the school center, or as a function applied to the grid of fish positions inside the size of the school.
#' @param pd  is the packing density of the school in fish per cubic meter.
#' @param pol  is the polarization of the school, given either by one value representing the mean angle deviation between individual and group heading (Huth and Wissel 1982), or as two values representing the mean angle deviation between individual and group heading in the azimuth angle and the elevation angle separately. If the polarization is to change between time steps, provide 'pol' as a one or two column matrix with time steps in the rows.
#' @param vel  is the velocity of the school center (x, y, and z) given either as one vector for constant velocity, or as a three column matrix for varying velocity.
#' @param T  is the number of time steps.
#' @param dt  is the time step difference (constant).
#' @param t0  is the start time in UNIX time.
#' @param x0  is the initial position of the school center.
#' @param gamma  is the three element vector of AR(1)-coefficients.
#' @param spacing  is the spacing between the fish, deduced from 'pd' if not given. 'spacing' can be used to define different spacing vertically and horizontally.
#' @param seed  is the seed if the simulation.
#' @param form  is a funciton used for transforming the positions of the school at time step 'i'.
#' @param maxfilesize  is the maximum size of the files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD car2sph NAs write.TSD zeros
#' @importFrom tools file_ext
#' @importFrom stats rnorm
#'
#' @export
#' @rdname ARsim_school
#'
ARsim_school<-function(con=NULL, size=c(30,20,10), shape=c("e","cc"), pd=1, pol=16.9, sd=NULL, vel=c(0.5,0,0), current=c(0,0,0), T=100, dt=1, t0=unclass(Sys.time()), pre=10, x0=c(0,0,0), gamma=c(0.8,0.8,0.8), spacing=NULL, seed=NULL, form=NULL, maxfilesize=3e8){
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-03-11 - Finished.
	########### DESCRIPTION: ###########
	# Generates and writes to file, fish positions which are AR(1) but not interacting, by the model in Holmin et al. 2012, "Simulations of multi-beam sonar echos from schooling individual fish in a quiet environment".
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---con--- is a the path to the directory in which to put the files.
	# ---size--- is a three element vector giving the x, y, and z size of the school (the axes, not the semi-axes).
	# ---shape--- is the shape of the school, either given as one of the two pre-defined shape strings "e" (ellipsoid) and "cc" (cubic compound school structure), as a three column matrix of initial fish positions relative to the school center, or as a function applied to the grid of fish positions inside the size of the school.
	# ---pd--- is the packing density of the school in fish per cubic meter.
	# ---pol--- is the polarization of the school, given either by one value representing the mean angle deviation between individual and group heading (Huth and Wissel 1982), or as two values representing the mean angle deviation between individual and group heading in the azimuth angle and the elevation angle separately. If the polarization is to change between time steps, provide 'pol' as a one or two column matrix with time steps in the rows.
	# ---vel--- is the velocity of the school center (x, y, and z) given either as one vector for constant velocity, or as a three column matrix for varying velocity.
	# ---T--- is the number of time steps.
	# ---dt--- is the time step difference (constant).
	# ---t0--- is the start time in UNIX time.
	# ---x0--- is the initial position of the school center.
	# ---gamma--- is the three element vector of AR(1)-coefficients.
	# ---spacing--- is the spacing between the fish, deduced from 'pd' if not given. 'spacing' can be used to define different spacing vertically and horizontally.
	# ---seed--- is the seed if the simulation.
	# ---form--- is a funciton used for transforming the positions of the school at time step 'i'.
	# ---maxfilesize--- is the maximum size of the files.
	##################################################
	##################################################
	##### Preparation #####
	# Split the connection into a part without file extension and the extension.
	con_ext <- file_ext(con)
	con_rest <- substr(con,1,nchar(con) - nchar(con_ext)-1)
	
	# Set the seed:
	set.seed(seed)
	# Use the semi axes:
	semisize <- size/2
	# Repeate 'gamma' to length 3:
	gamma <- rep(gamma,3)
	# Get the speed of the school:
	speed <- sqrt(sum(vel^2))
	
	# Get the spacing of the fish positions:
	if(length(spacing)==0){
		spacing <- rep(pd^(-1/3),3)
		}
	
	# Get the initial fish positions:
	if(is.numeric(shape)){
		fishgrid <- shape
		}
	else{
		# Define the grid of fish positions:
		fishgrid <- as.matrix(expand.grid(seq(-semisize[1],semisize[1],spacing[1]), seq(-semisize[2],semisize[2],spacing[2]), seq(-semisize[3],semisize[3],spacing[3])))
		# Get the indices of the fish inside the object enclosing the school:
		if(shape[1]=="e"){
			inside <- rowSums(fishgrid^2/matrix(semisize^2, nrow=nrow(fishgrid), ncol=3, byrow=TRUE))<1
			}
		else if(shape[1]=="cc"){
			not1 <- (fishgrid[,1] < -3/5*semisize[1]) & (abs(fishgrid[,2]) < 1/3*semisize[2])
			not2 <- (fishgrid[,1] > 3/5*semisize[1]) & (abs(fishgrid[,2]) < 1/3*semisize[2])
			not3 <- (abs(fishgrid[,1]) < 1/5*semisize[1]) & (fishgrid[,2] < -1/3*semisize[2])
			not4 <- (abs(fishgrid[,1]) < 1/5*semisize[1]) & (fishgrid[,2] > 1/3*semisize[2])
			inside <- !(not1 | not2 | not3 | not4)
			}
		else if(is.function(shape)){
			inside <- shape(fishgrid)
			}
		# Select only fish inside the the object enclosing the school:
		fishgrid <- fishgrid[inside,]
		}
	N <- nrow(fishgrid)
	
	### Expand to all valid time steps: ###
	if(length(dim(current))==0){
		current <- matrix(current, nrow=T, ncol=3, byrow=TRUE)
		}
	if(length(dim(vel))==0){
		vel <- matrix(vel, nrow=T, ncol=3, byrow=TRUE)
		}
	if(length(pol)==1){
		pol <- matrix(pol, nrow=T, ncol=1, byrow=TRUE)
		}
	else if(length(pol)==2){
		pol <- matrix(pol, nrow=T, ncol=2, byrow=TRUE)
		}
	if(length(sd)==0){
		sd <- pol2sd(pol, speed=speed, dt=dt, gamma=mean(gamma), type="h")$sigma
		}
	if(length(dim(sd))==0){
		sd <- matrix(sd,nrow=T,ncol=3,byrow=TRUE)
		}
	###
	
	# Expand 'vel', the time, 'sd' to all time steps, also the preliminary time steps:
	# First the current, since 'x0' is altered in the next line:
	current <- rbind(matrix(current[1,], nrow=pre, ncol=ncol(current), byrow=TRUE), current)
	vel <- rbind(matrix(vel[1,], nrow=pre, ncol=ncol(vel), byrow=TRUE), vel)
	pol <- rbind(matrix(pol[1,], nrow=pre, ncol=ncol(pol), byrow=TRUE), pol)
	sd <- rbind(matrix(sd[1,], nrow=pre, ncol=ncol(sd), byrow=TRUE), sd)
	# Calculate the displacement of the school center:
	x_current <- -matrix((pre+1)*current[1,]*dt, nrow=pre+T,ncol=ncol(current),byrow=TRUE) + apply(current*dt,2,cumsum)
	x         <- -matrix((pre+1)*vel[1,]*dt,     nrow=pre+T,ncol=ncol(vel),    byrow=TRUE) + apply(vel*dt,2,cumsum)
	# Define the time steps and the time step indices:
	tseq <- seq_len(T+pre)
	t <- seq(-pre,T-1)*dt
	
	
	##### Execution and output #####
	psxf1 <- NAs(N)
	psyf1 <- NAs(N)
	pszf1 <- NAs(N)
	ARxf1 <- zeros(N)
	ARyf1 <- zeros(N)
	ARzf1 <- zeros(N)
		
	# Define parameters used when plotting the time bar for the generation of bottom points:
	infostring <- "Processing time steps:"
	cat(infostring,"\n",sep="")
	totalsteps <- length(t)
	stepfact <- nchar(infostring)/totalsteps
	oldvalue <- 0
	
	filenr <- 0
	# Move through the time steps:
	for(i in tseq){
		# Print a dot if the floor of the new value exceeds the old value in:
		thisvalue <- floor(i*stepfact)
		if(thisvalue > oldvalue){
			cat(rep(".", thisvalue-oldvalue), if(i==totalsteps) "\n", sep="")
			oldvalue <- thisvalue
			}
		# Transform the within school positions:
		if(is.function(form)){
			fishgrid <- form(fishgrid,i)
			}
		# Generate the noise in the perturbations:
		GNxf <- rnorm(N, sd=sd[i,1])
		GNyf <- rnorm(N, sd=sd[i,2])
		GNzf <- rnorm(N, sd=sd[i,3])
		# Autoregress the perturbations:
		ARxf2 <- ARxf1*gamma[1] + GNxf
		ARyf2 <- ARyf1*gamma[2] + GNyf
		ARzf2 <- ARzf1*gamma[3] + GNzf
		# Get the new positions of the fish:
		psxf2 <- x0[1] + x[i,1] + ARxf2 + fishgrid[,1]
		psyf2 <- x0[2] + x[i,2] + ARyf2 + fishgrid[,2]
		pszf2 <- x0[3] + x[i,3] + ARzf2 + fishgrid[,3]
		psxf2_withCurrent <- psxf2 + x_current[i,1]
		psyf2_withCurrent <- psyf2 + x_current[i,2]
		pszf2_withCurrent <- pszf2 + x_current[i,3]
		# If not in the preliminary time steps, store the data to file:
		totalbytes <- 0
		if(i>pre){
			# Get the rotation angles of the fish from the current and previous fish positions:
			rthetaphi <- car2sph(psxf2-psxf1, psyf2-psyf1, pszf2-pszf1)
			#towrite <- list(utim=t0+t[i], psxf=psxf2, psyf=psyf2, pszf=pszf2, rtzf=rthetaphi[,2], rtxf=rthetaphi[,3])
			# Subtract the current:
			towrite <- list(utim=t0+t[i], psxf=psxf2_withCurrent, psyf=psyf2_withCurrent, pszf=pszf2_withCurrent, vlxf=(psxf2-psxf1)/dt, vlyf=(psyf2-psyf1)/dt, vlzf=(pszf2-pszf1)/dt, rtzf=rthetaphi[,2]-pi/2, rtxf=pi/2-rthetaphi[,3], crrt=current[i,], velS=vel[i,], polS=pol[i,], sdAR=sd[i,])
			# Issue av warning if there are fish above the surface:
			if(any(pszf2_withCurrent>0)){
				warnings("Fish generated that are above the sea surface")
				}
			if(length(con)>0){
				if(totalbytes<maxfilesize){
					totalbytes <- totalbytes + write.TSD(towrite, con, append=i>pre+1)
					}
				else{
					filenr <- filenr + 1
					con <- paste0(con_rest, "_", filenr, ".", con_rest)
					totalbytes <- write.TSD(towrite, con)
					}
				
				}
			}
		# Store for the nest time step:
		psxf1 <- psxf2
		psyf1 <- psyf2
		pszf1 <- pszf2
		# Reset the previous AR-series values:
		ARxf1 <- ARxf2
		ARyf1 <- ARyf2
		ARzf1 <- ARzf2
		}
	
	# Return the last written data:
	invisible(towrite)
	##################################################
	##################################################
	}
	
