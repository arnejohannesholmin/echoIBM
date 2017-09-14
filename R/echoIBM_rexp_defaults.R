#*********************************************
#*********************************************
#' Creates a list of default settings for noise generation methods used in echoIBM, specifically the Multiple Sines - method ("ms"), the rearrangement method of uniform independent variables resulting in correlated and autocorrelated beams, and the simple independent exponential distribution (in which case only 'ssed' is defaulted).
#'
#' @param noise  is a vector of character strings of length 2, specifying which types of noise to apply to the data:
#' @param data  is a list of the required beam configuration information, including length 'lenb' of the beams, number 'numb' of beams and the frequency 'freq' of the beams, as well as the sonar/echosounder type 'esnm'.
#' @param parlist  is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR sonR_implemented noise.path
#' @importFrom TSD read.TSDs strff zeros
#'
#' @export
#' @rdname echoIBM_rexp_defaults
#'
echoIBM_rexp_defaults<-function(noise="ms", indt=NULL, data=list(), parlist=list()){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-02-28 - Clean version.
	# Update: 2012-11-12 - Changed to read default overlap values only if missing in the data.
	# Update: 2013-09-26 - Fundamental restructuring and implementation of sub-functions.
	# Last: 2013-11-05 - Changed to also return the input characters representing the type of overlap between the voxels ('input_olpn' and 'input_olps').
	########### DESCRIPTION: ###########
	# Creates a list of default settings for noise generation methods used in echoIBM, specifically the Multiple Sines - method ("ms"), the rearrangement method of uniform independent variables resulting in correlated and autocorrelated beams, and the simple independent exponential distribution (in which case only 'ssed' is defaulted).
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---noise--- is a vector of character strings of length 2, specifying which types of noise to apply to the data:
	##		"bg" - Background noise
	##		"pn" - Periodic noise
	##		"hi" - High intensity noise
	##		"nr", "np", "cp" - Near-range noise generated in passive mode
	##		"na", "ca" - Near-range noise generated in passive mode, opossibly originating from echo in the vessel
	##		"ex" - Independent exponential noise due to acoustic interference (applied on the near range noise / reverberation noise ("nr"), background noise ("bg") and the simulated echo)
	##		"acex", "aex", "cex" - Correlated exponential noise due to acoustic interference (applied on the near range noise / reverberation noise ("nr"), background noise ("bg") and the simulated echo). Generated using the rearrangement method developed by Holmin and Tjøstheim, in which independent uniformly distributed variables are drawn, and rearranged so that among the next 'l' voxels of a vector, the voxel following position 'k' is chosen to resemble the k'th voxel. This method has the ability to include both autocorrelation along vectors and correlation between vectors, because the voxel following the k'th voxel can be chosen according to a proximity criterion to the k'th voxel in several vector simulatneously, given specific weights that result in correlation between vectors. Of cource this method of rearrangement is ad-hoc and limited in which autororrelation and correlation structures that are achievable
	##		"ms" - Multiple sine waves: In this noise generation method the process of superimposed sine waves that occur when multiple targets scatter sound back to the receiver, is mimiced in order to produce the randomness present in real acoustic data. This method is more time consuming and depend on the number of targets contributing to each voxel (constant in this function), the  length of the sound pulse (affecting the autocorrelation along vectors) and the overlap between neighboring vectors, whereas the number of sample points and the number of preiods of the sine waves for each voxel affect the accuracy of the noise simulation. An advantage of this method is that the non-rayleigh distribution occuring when the number of targets is low can be acheived, but generally, in the simulations the number of targets is chosen high to produce exponentially distributed acoustical intensity
	# ---data--- is a list of the required beam configuration information, including length 'lenb' of the beams, number 'numb' of beams and the frequency 'freq' of the beams, as well as the sonar/echosounder type 'esnm'.
	# ---parlist--- is a list of input parameters to the function echoIBM.add.noise(), which generates noise and randomness. See echoIBM.add.noise() for explaination of the possible variables.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	### Functions: ###
	# A function for extracting overlap values for the signal or noise for the MS70 sonar (type="n" refers to noise, and type="s" refers to signal):
	get.olpn <- function(parlist, noise, esnm="ms70", type="n", input_out=FALSE){
		# Function for reading the overlap values of 'type' "n" (noise) or "s" (signal):
		read.olp<-function(data, olptype="olpb"){
			if(length(data[[olptype]])==0){
				suppressWarnings(read.TSDs(noise.path(esnm=esnm,utim=parlist$utim),var=olptype,dimension=TRUE)[[olptype]])
			}
			else{
				data[[olptype]]
			}
		}
			
		# Define the legal character values of 'olpn':
		# "c" - constant correlation between neighboring beams
		# "b" - variable correlation between neighboring beams  
		# "v" - variable correlation between neighboring beams for each voxel
		# "p" - variable correlation between neighboring beams for each voxel and each ping (used with periodic noise)
		if(sonR_implemented(esnm, type=c("MBE"))){
			if(strff("n",type)){
				legal <- c("c","b","v")
			}
			else if(strff("s",type)){
				legal <- c("c","b")
			}
			else{
				stop("Wrong value of 'type' in get.olpn()")
			}
		}
		# Other sonars:
		else if(sonR_implemented(esnm, type=c("MBE","SBE","OFS"))){
			legal <- c("c")
		}
		else{
			legal <- c("c")
			# Set parlist$olpn=1 later:
			warning("The acoustical system given by 'esnm' has not been fully implemented. Zero correlation between beams was applied for the noise")
		}
			
		# Defaults if not given, or not numeric and not in c(legal,"p"):
		if(length(parlist$olpn)==0 || (!is.numeric(parlist$olpn) && !any(strff(c(legal,"p"),parlist$olpn)))){
			if(sonR_implemented(esnm, type=c("MBS"))){
				parlist$olpn <- "v"
			}
			else if(sonR_implemented(esnm, type=c("OFS"))){
				#parlist$olpn <- "b" # Should be used, but requires estimation of the coorelation between beams from real data with little or no signal (preferably passive data), which is currently unavailable:
				parlist$olpn <- "c"
			}
			else if(sonR_implemented(esnm, type=c("MBE","SBE"))){
				parlist$olpn <- "c"
			}
			else{
				parlist$olpn <- "c"
			}
		}
		# If input_out==TRUE, return the character naming the type of overlap, to be used in echoIBM_rexp_MultSines() for locating the appropriate overlap file:
		if(input_out){
			if(length(parlist$olpn)==1 && is.character(parlist$olpn)){
				return(input_olpn <- parlist$olpn)
			}
			else{
				warning("Something did not proceed as expected. The variable olpn has allready been converted to numeric. olpn=\"v\" assumed")
				return(input_olpn="v")
			}
		}
		
		# Constant correlation for all voxels:
		if(strff("c",parlist$olpn) && "c" %in% legal){
			# See "Test_of_rexp_MultSines.R in the "extdata" directory of the echoIBM package.
			if(sonR_implemented(esnm, type="MBS")){
				parlist$olpn <- c(0.46,1,0.46)
			}
			else if(sonR_implemented(esnm, type="OFS")){
				parlist$olpn <- c(0.3,1,1,1,0.3)
			}
			else if(sonR_implemented(esnm, type=c("MBE","SBE"))){
				parlist$olpn <- 1
			}
			else{
				parlist$olpn <- 1
			}
		}
		# Variable correlation between beams:
		else if(strff("b",parlist$olpn) && "b" %in% legal){
			# Read the overlap values of the background noise 'olpb' (see "S2009116_PG.O.Sars[4174] - olp and acf.R"):
			parlist$olpn <- read.olp(data, olptype="olpb")
			# Set the dimension to [#correlated beams, #beams in each fan, #fans]:
			dim(parlist$olpn) <- c(dim(parlist$olpn)[1],parlist$luqf,parlist$nuqf)
		}
		# Variable correlation for all voxels:
		else if(strff("v",parlist$olpn) && "v" %in% legal){
			if(!"pn" %in% noise){
				warning("parlist$olpn=\"v\" usually implies that the periodic noise is applied to the simulations, while this is not specified in the parameter 'noise'")
			}
			# Read the overlap values of the total noise including periodic noise 'olpt' (see "S2009116_PG.O.Sars[4174] - olp and acf.R"):
			parlist$olpn <- read.olp(data, olptype="olpt")
			# If constant phase angels are applied (saves time), the overlap array stored in noise.path() is used. If variable phases are needed, this must be indicated by olpn=="pings" or "p", and is addressed in echoIBM.add.noise():
			########### Set the dimension to [#correlated beams, #beams in each fan, #fans]:
			dim(parlist$olpn) <- c(dim(parlist$olpn)[1:2],parlist$luqf,parlist$nuqf)
		}
		# Return the overlap values:
		parlist$olpn
	}
	### End of functions: ###
	
	# Extract length and number of beams:
	if(!all(c("lenb","numb","freq") %in% names(data))){
		stop("'data' must contain all of the variables \"lenb\", \"numb\" and \"freq\"")
		}
	
	# Treat the seed:
	if(length(parlist$seed)==0 || length(indt)==0){
		stop("'indt' and 'seed' must be given, either as a single integer or a vector of integers, or possibly a code word such as 'ind', implying to use the time step indices as seeds")
	}
	if(strff("ind", parlist$seed)){
		parlist$seed <- indt
	}
	else if(length(parlist$seed)!=length(indt)){
		warning("Random seed 'seed' in the parameter list 'parlist' did not have the same length as the number of time steps, and was repeated to this length")
		parlist$seed <- rep(parlist$seed, length.out=length(indt))
	}
		
	
	########## Execution ##########
	##### (1) Get the length of the beams: #####
	if(is.null(parlist$J)){
		parlist$J <- data$lenb[1]
		}	
	
	##### (2) Get the structure of the beams: #####
	if(is.null(parlist$nuqf) || is.null(parlist$luqf)){
		parlist$nuqf <- length(unique(data$freq))
		parlist$luqf <- data$numb/parlist$nuqf
		if(parlist$luqf %% 1 != 0){
			warning("The data indicate a non-rectangular grid of beams")
			}
		}	
	
	##### (3) Get the parameters for the correlated exponential distribution: #####
	if(any(c("ms") %in% noise)){
		### (3a, 3b, 3c) Standard defaults regardless of sonar/echosounder. The combination of these parameters provides a good exponential distribution with acceptable cpu-time. See "Test_of_rexp_MultSines.R": ###
		if(is.null(parlist$L)){
			parlist$L <- 3
			}
		if(is.null(parlist$N)){
			parlist$N <- 40
			}
		if(is.null(parlist$P)){
			parlist$P <- 10
			}
		
		### (3d) Get the duration of the sine waves: ###
		if(is.null(parlist$w)){
			# EK60 multifrequency echosounder:
			if(sonR_implemented(data$esnm, type=c("SBE"))){
				# Get lengths of the sine waves in the multisine-method for EK60 (correlations found in Holmin et al. (2013)):
				#corEK60 <- c(0.31,0.18,0.12,0.02,-0.06,0.15) # Old 2013-09-18
				corEK60 <- c(0.31,0.18,0.12,0,0,0)
				w <- zeros(6)
				for(i in seq_along(corEK60)){
					w[i] <- acftable[,1][which.min(abs(acftable[,3]-corEK60[i]))]
					}
				# [1] 2.48 1.84 1.62 1.23 1.23 1.23
				parlist$w <- w
				}
			# ME70 multibeam echosounder:
			else if(sonR_implemented(data$esnm, type=c("MBE"))){
				# Assume pulselength ≈ 2 ms and sampleinterval duration 0.128 ms, giving w = 2 /0.125 * 3/4 = 12:
				parlist$w <- rep(12,parlist$nuqf)
				}
			# MS70 multibeam sonar:
			else if(sonR_implemented(data$esnm, type=c("MBS"))){
				# Set the length of the pulses in units of sample intervals. In reality this is set to be 4, but the autocorrelation estimates of the MS70 sonar suggest 3:
				parlist$w <- rep(3,parlist$nuqf)
				}
			# SX90 multibeam sonar:
			else if(sonR_implemented(data$esnm, type=c("OFS"))){
				# From the experience with simulations of the MS70 sonar, which normally has 2 ms pulselength and 1/2 ms sampling interval, which should imply w=4, the value w=3 was more appropriate. Thus we propose to estimate 'w' by plsl/sint * 3/4:
				if(all(c("plsl","sint") %in% names(data))){
					parlist$w <- rep(data$plsl/data$sint * 3/4, parlist$nuqf)
					}
				# Otherwise, use the values for 600 meter maximum range of the SX90:
				else{
					plsl_600 <- 8
					sint_600 <- 0.25
					parlist$w <- rep(plsl_600/sint_600 *3/4, parlist$nuqf)
					}
				}
			# Other acoustical systems:
			else{
				warning("The acoustical system given by 'esnm' has not been fully implemented. Standard value w=3 used")
				parlist$w <- rep(3,parlist$nuqf)
				}
			}
		
		### (3e) Get the between beam overlap for the signal: ###
		parlist$input_olps <- get.olpn(parlist,noise,esnm=data$esnm[1],type="s",input_out=TRUE)
		parlist$olpn <- parlist$input_olps
		parlist$olps <- get.olpn(parlist,noise,esnm=data$esnm[1],type="s")
		
		### (3f) Get the between beam overlap for the noise: ###
		# Get the character value of the olp variable, and also set the variable parlist$olpn:
		parlist$input_olpn <- get.olpn(parlist,noise,esnm=data$esnm[1],input_out=TRUE)
		parlist$olpn <- parlist$input_olpn
		parlist$olpn <- get.olpn(parlist,noise,esnm=data$esnm[1])
		}
	##### End of "Default settings for the Multiple Sines method" #####
	
	
	##### (4) Default settings for the rearrangement method (no longer maintained to the extent of the "ms" method): #####
	else if(any(c("acex","aex","cex") %in% noise)){
		# Standard defaults regardless of sonar/echosounder. See "Test of echoIBM_rexp_Rearr":
		if(is.null(parlist$rate)){
			parlist$rate <- 1
			}
		if(is.null(parlist$prob)){
			parlist$prob <- 1
			}
		if(is.null(parlist$l)){
			parlist$l <- 100
			}
		# Add to values to avoid the strange error when (almost exactly) 1324 voxels and 500 beams are drawn:
		parlist$buffer <- 100
		
		# Sonar/echosounder specific parameters. See "Test_of_rexp_MultSines.R":
		if(sonR_implemented(data$esnm, type=c("MBS"))){
			if(is.null(parlist$rho)){
				parlist$rho <- c(1,1,1)
				}
			if(is.null(parlist$w)){
				parlist$w <- 3
				}
			# Get the number of significant beams:
			parlist$wC <- length(parlist$rho)
			# Define the "correlation" matrix:
			totheside <- (parlist$wC-1)/2
			S <- zeros(parlist$luqf, 2*totheside + parlist$luqf)
			for(i in seq_len(parlist$luqf)){
				S[i,seq_along(parlist$rho)+i-1] <- parlist$rho
				}
			S <- S[,totheside+seq_len(parlist$luqf)]
			parlist$C=zeros(data$numb,data$numb)
			for(i in seq_len(parlist$nuqf)){
				parlist$C[seq_len(parlist$luqf)+(i-1)*parlist$luqf, seq_len(parlist$luqf)+(i-1)*parlist$luqf] <- S
				}
			}
		else{
			if(is.null(parlist$rho)){
				parlist$rho <- 1
				}
			if(is.null(parlist$w)){
				parlist$w <- 3
				}
			parlist$wC <- 1
			# Define identity "correlation" matrix:
			parlist$C <- diag(data$numb)
			}
		
		# Create the array of indexes for which parlist$C>0 and the width of this array (defaulted to 3 for echoIBM()):
		if(length(parlist$Cind)==0){
			parlist$Cind <- matrix(0:(parlist$wC-1) - (parlist$wC-1)/2,nrow=data$numb,ncol=parlist$wC,byrow=TRUE) + 0:(data$numb-1)
			}
		}
	else if(any("bk" %in% noise)){
		# N=40 iterations seems to be sufficient for establishing the Barakat PDF to some degree of approximation, at lest for nsig>=3. Also the threshold for nsig, 'nsth' is set to 2 based on the plot in "Test Barakat.R":
		parlist$nsth <- 2
		parlist$N <- 40
		parlist$max.memory <- 1e9
		}
		
	
	########## Output ##########
	parlist
	##################################################
	##################################################
	}
