#*********************************************
#*********************************************
#' Simulates one radial sphere of one echo sounder observation based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school. One for loop used (for each frequency).
#'
#' @param j		is the sampling interval index.
#' @param data	is the list of data for the simulation, containing all of the following elements:
#' \itemize{
#'   \item numb
#'   \item lthesel
#'   \item transducerposL
#'   \item thesel
#'   \item rres
#'   \item validr
#'   \item psxf
#'   \item psyf
#'   \item pszf
#'   \item psxv
#'   \item psyv
#'   \item pszv
#'   \item rtzv
#'   \item rtxv
#'   \item rtyv
#'   \item psze
#'   \item sigma0mode
#'   \item uniquek
#'   \item wavenumber
#'   \item dira
#'   \item dire
#'   \item rad1
#'   \item bpt1
#'   \item indi
#'   \item equalbp_em_re
#'   \item rad2
#'   \item bpt2
#'   \item bptf
#'   \item lenl
#'   \item ssif
#'   \item grsf
#'   \item absr
#'   \item fish
#'   \item asps
#'   \item epss
#'   \item epsl
#' }
#' @param split	Logical: If TRUE, do the calculation of the transducer beam patterns for each beam, and not for all beams of the same frequency. This saves memory, but increases process time.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR rotate3D
#' @importFrom TSD NAs zeros
#' @importFrom stats approx
#'
#' @export
#'
echoIBM.oneping.oneschool.onesample <- function(j, data, split=FALSE){
	
	############### LOG: ###############
	# Start: 2008-03-10 - Clean version.
	# Update: 2011-02-20 - Changed to use beamPattern.TSD() which takes as input a list of inputs named by the convencional TSD-names. Thus the parameter "mod" as a list element of the beam pattern function is no longer neede.
	# Update: 2012-02-07 - Changes to the treatment of sonar beam patterns, and added the side lobe level factors 'sllf' used for adjusting side lobe levels to the desired level.
	# Update: 2012-11-29 - Expanded the dump.
	# Last: 2018-08-21 - Cleaned up the code, merging echoIBM.oneping.j1 and echoIBM.oneping.j2 to one function.
	
	##################################################
	##################################################
	# Function for appending the mean and max of the specified variable to the dumpdata:
	addMeanMaxToDumpdata <- function(name, dumpdata){
		x <- get(name)
		dumpdata[[name]] <- rbind(dumpdata[[name]], c(Mean=mean(x, na.rm=TRUE), Max=max(x, na.rm=TRUE)))
		dumpdata
	}
	
	# Get number of fish and beams:
	Nf <- data$lthesel
	Nb <- length(data$wavenumber)
	
	# Define outputs:
	out <- zeros(data$numb)
	outnext <- zeros(data$numb)
	nsig <- zeros(data$numb)
	
	# Declare the dump data to return:
	dumpnames <- c("etaj", "B_T1", "B_T2", "B_L", "etaa", "epss", "sigma0mode", "withoutetaj")
	
	
	if(Nf == 0){
		# If no fish are simulated return NAs:
		dumpdata <- as.data.frame(zeros(2, length(dumpnames)))
		names(dumpdata) <- dumpnames
	}
	# If more than 0 fish are to be treated:
	if(Nf > 0){
		
		# Declare the dumpdata to be a list of empty elements:
		dumpdata <- vector("list", length(dumpnames))
		names(dumpdata) <- dumpnames
	
		# The radial weighting due to the time interval of the incoming sound (to the receiver):
		etaj <- 1 - abs(data$transducerposL[data$thesel,1] / data$rres - data$validr[j] + 1)
		# Dump 'etaj':
		#dumpdata$etaj <- mean(etaj, na.rm=TRUE)
		dumpdata <- addMeanMaxToDumpdata("etaj", dumpdata)
		
		
		# Extracting the positions of the fish in the coordinate system of the vessel (V):
		data$fishposV <- rotate3D(
			cbind(data$psxf[data$thesel], data$psyf[data$thesel], data$pszf[data$thesel]) - matrix(c(data$psxv,data$psyv,data$pszv), nrow=Nf, ncol=3, byrow=TRUE), 
			by="zxy", 
			ang=c(data$rtzv, data$rtxv, data$rtyv),
			drop.out=TRUE)
		if(Nf==1){
			dim(data$fishposV) <- c(1,3)
			}
		data$fishposV[,3] <- data$fishposV[,3] - data$psze
		
		
		# For loop through the frequencies:
		for(i in seq_along(data$uniquek)){
			
			# Which beams are to be treated:
			thesebeams <- which(data$wavenumber==data$uniquek[i])
			Nub <- length(thesebeams)
			
			
			if(split){
				# Reserving memory for the beam patterns of the echo sounder as matrices of dimension c(data$lthesel,data$numb):
				B_T1 <- B_T2 <- zeros(Nf, Nub)
		
				# Extracting the beam pattern values at the fish:
				# 'iii' is the position in the matrices 'B_T1_i' and 'B_T2_i' in the for loop:
				for(i in seq_len(Nub)){
					# Beam pattern on emission from the transducer:
					B_T1[,i] <- echoIBM_setB_T(data, direction=1, beamnr=thesebeams[i])
					# Beam pattern on reception at the trasducer:
					if(data$equalbp_em_re){
						B_T2[,i] <- B_T1[,i]
					}
					else{
						B_T2[,i] <- echoIBM_setB_T(data, direction=2, beamnr=thesebeams[i])
					}
				} # End of for(ii in thesebeams).
			}
			else{
				# Beam pattern on emission from the transducer:
				B_T1 <- echoIBM_setB_T(data, direction=1, beamnr=thesebeams)
				# Beam pattern on reception at the trasducer:
				if(data$equalbp_em_re){
					B_T2 <- B_T1
				}
				else{
					B_T2 <- echoIBM_setB_T(data, direction=2, beamnr=thesebeams)
				}
			}
			
			
			# Take the sum of the beam pattern of all beams of the current frequency for each fish:
			dim(B_T1) <- c(Nf, Nub)
			B_T1 <- matrix(rowSums(B_T1), nrow=Nf, ncol=Nub)
			
			
			# Beam pattern of the fish:
			# If each fish is exposed to different frequencies, there are different beam patterns from the same fish for the different frequencies:
			B_L <- data$bptf(
				list(
					dira = data$transducerposL[data$thesel,2], 
					dire = data$transducerposL[data$thesel,3], 
					wnsz = data$uniquek[i] * data$lenl[data$thesel]
				)
			)
			
			# Absorption, assuming that the transducer and receiver is located on the same place (the same physical object). The expression differs from the one in the documentation due to the definition of the absorption coefficient used by Simrad, which is absr=alpha/10, where 'alpha' is the absorption coefficient used in the documentation:
			etaa <- 10^(outer(-0.2*data$transducerposL[data$thesel,1] , data$absr[thesebeams],"*"))
			
			
			### Add to the output:
			# data$sigma0mode==1: If the optimal backscattering cross section 'sgbs' (sigma_bs) is present, no relation to target size is specified:
			# data$sigma0mode==2: If the coefficient 'epss' linking 'sgbs' to fish size is present, it may be a function of frequency, or simply a numeric vector:
			# See also echoIBM.oneping.oneschool() and echoIBM.default.oneschool() for explaination of data$sigma0mode:
			#if(Nf>100) browser()
			fish <- data$fish[data$thesel]
			if(data$sigma0mode==1){
				# Simply multiply the varialbes 'fish', 'etaj', 'B_T1', 'B_L', 'B_T2' and 'etaa':
				withoutetaj <- fish * B_T1 * B_L * B_T2 * etaa
			}
			else if(data$sigma0mode==2){
			    frequency <- data$uniquek[i] * data$asps / (2*pi)
			    # Include the echo ability 'epss' as a function of frequency:
				if(is.function(data$epss)){
				    epss <- data$epss(frequency)
					withoutetaj <- fish * B_T1 * B_L * B_T2 * etaa * epss
				}
			    # Include the echo ability 'epss' as from a table of frequency and response:
			    else if(length(dim(data$epss)) == 2){
			        epss <- approx(data$epss, xout = frequency)$y
			        withoutetaj <- fish * B_T1 * B_L * B_T2 * etaa * epss
			    }
			    # Include the echo ability 'epss' as a vector of the same length as the number of fish:
			    else{
			        epss <- data$epss[data$thesel]
			        withoutetaj <- fish * B_T1 * B_L * B_T2 * etaa * epss
			    }
			    # Dump 'esps':
				#dumpdata$epss <- c(dumpdata$epss, mean(epss, na.rm=TRUE))
				dumpdata <- addMeanMaxToDumpdata("epss", dumpdata)
			}
			else if(data$sigma0mode %in% 3:4){
				warning("sigma0mode = 3 or 4 is no longer supported")
			}
			
			# Calculate the backscatter for the voxels of the current radial distance:
			nsig[thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
			out[thesebeams] <- colSums(withoutetaj * etaj)
			outnext[thesebeams] <- colSums(withoutetaj * (1 - etaj))
			
			
			# Dump. Using the mean is ok since there is an equal number of beams for each frequency, as per restriction in the echoIBM package:
			dumpdata <- addMeanMaxToDumpdata("B_T1", dumpdata)
			dumpdata <- addMeanMaxToDumpdata("B_T2", dumpdata)
			dumpdata <- addMeanMaxToDumpdata("B_L", dumpdata)
			dumpdata <- addMeanMaxToDumpdata("etaa", dumpdata)
			dumpdata <- addMeanMaxToDumpdata("withoutetaj", dumpdata)
			
			#dumpdata$B_T1 <- c(dumpdata$B_T1, mean(B_T1, na.rm=TRUE))
			#dumpdata$B_T2 <- c(dumpdata$B_T2, mean(B_T2, na.rm=TRUE))
			#dumpdata$B_L <- c(dumpdata$B_L, mean(B_L, na.rm=TRUE))
			#dumpdata$etaa <- c(dumpdata$etaa, mean(etaa, na.rm=TRUE))
			#dumpdata$withoutetaj <- c(dumpdata$withoutetaj, mean(withoutetaj, na.rm=TRUE))
		} # End of for(i in seq_along(data$uniquek)).
		
		# Dump:
		dumpdata <- lapply(dumpdata, function(x) if(length(x)==0) cbind(Mean=0, Max=0) else x)
		dumpdata <- as.data.frame(lapply(dumpdata, colMeans))
	} # End of if(Nf>0).
	
	# Add the sigma0mode:
	dumpdata <- cbind(dumpdata, sigma0mode=c(Mean=data$sigma0mode, Max=data$sigma0mode))
	
	# Return:
	list(sv=out, svnext=outnext, nsig=nsig, dumpdata=dumpdata)
	##################################################
	##################################################
}
