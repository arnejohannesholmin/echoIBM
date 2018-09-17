#*********************************************
#*********************************************
#' Simulates one radial sphere of one echo sounder observation based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school. One for loop used (for each frequency).
#'
#' @param j  is the sampling interval index.
#' @param data  is the list of data for the simulation, containing all of the following elements:
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
#' @rdname echoIBM.oneping.j1
#'
echoIBM.oneping.j1 <- function(j, data){
	
	############### LOG: ###############
	# Start: 2008-03-10 - Clean version.
	# Update: 2011-02-20 - Changed to use beamPattern.TSD() which takes as input a list of inputs named by the convencional TSD-names. Thus the parameter "mod" as a list element of the beam pattern function is no longer neede.
	# Update: 2012-02-07 - Changes to the treatment of sonar beam patterns, and added the side lobe level factors 'sllf' used for adjusting side lobe levels to the desired level.
	# Last: 2012-11-29 - Expanded the dump.
	############ VARIABLES: ############
	# ---j--- is the sampling interval index.
	# ---data--- is the list of data for the simulation, containing all of the following elements:
	#		"numb"
	#		"lthesel"
	#		"transducerposL"
	#		"thesel"
	#		"rres"
	#		"validr"
	#		"psxf"
	#		"psyf"
	#		"pszf"
	#		"psxv"
	#		"psyv"
	#		"pszv"
	#		"rtzv"
	#		"rtxv"
	#		"rtyv"
	#		"psze"
	#		"ssil"
	#		"sigma0mode"
	#		"uniquek"
	#		"wavenumber"
	#		"dira"
	#		"dire"
	#		"rad1"
	#		"bpt1"
	#		"indi"
	#		"equalbp_em_re"
	#		"rad2"
	#		"bpt2"
	#		"bptf"
	#		"lenl"
	#		"ssif"
	#		"grsf"
	#		"absr"
	#		"fish"
	#		"asps"
	#		"epss"
	#		"epsl"
	
	
	##################################################
	##################################################
	thisdata <- list()
	out <- zeros(data$numb)
	outnext <- zeros(data$numb)
	nsig <- zeros(data$numb)
	
	# Declare the dump data to return:
	dumpdata <- NAs(8)
	names(dumpdata) <- c("etaj", "B_L", "etaa", "epss", "epsl", "chi", "sigma0mode", "withoutetaj")
	
	# If more than 0 fish are to be treated:
	if(data$lthesel>0){
		
		# The radial weighting due to the time interval of the incoming sound (to the receiver):
		thisdata$etaj <- 1 - abs(data$transducerposL[data$thesel,1]/data$rres-data$validr[j] + 1)
		# Dump 'etaj':
		dumpdata[1] <- sum(thisdata$etaj, na.rm=TRUE)
		
		# Extracting the positions of the fish in the coordinate system of the vessel (V):
		thisdata$fishposV <- rotate3D(
			cbind(data$psxf[data$thesel], data$psyf[data$thesel], data$pszf[data$thesel]) - matrix(c(data$psxv,data$psyv,data$pszv), ncol=3, nrow=data$lthesel, byrow=TRUE), 
			by="zxy", 
			ang=c(data$rtzv, data$rtxv, data$rtyv),
			drop.out=TRUE)
		if(data$lthesel==1){
			dim(thisdata$fishposV) <- c(1,3)
			}
		thisdata$fishposV[,3] <- thisdata$fishposV[,3] - data$psze
		
		# Assure that data$ssil has the apropriate dimension:
		if(!is.null(data$ssil) && any(data$sigma0mode %in% 3:4)){
			data$ssil <- matrix(rep(data$ssil, length.out=max(data$thesel) * length(data$uniquek)), nrow=data$lthesel, ncol=length(data$uniquek))
			}
		
		# For loop through the frequencies:
		for(i in seq_along(data$uniquek)){
			
			# Which beams are to be treated:
			thisdata$thesebeams <- which(data$wavenumber==data$uniquek[i])
			
			# Fish positions in (T) for the current beam, giving the fish directions 'data$fishdirT' in the current beam:
			thisdata$fishdirT <- rotate3D(
				thisdata$fishposV, 
				"zx", 
				cbind(data$dira[thisdata$thesebeams]-pi/2,-data$dire[thisdata$thesebeams]), 
				drop.out=FALSE, 
				sph.out=TRUE)
			
			# Beam pattern on emission at the trasducer as a matrix of dimension c(data$lthesel,data$kfreq), where kfreq is the number of beams of the same frequency:
			if(data$pbp1=="circularPiston_ellipticRadius_sidelobefit"){
				# Use as input the 'dira' [lthesel*kfreq], 'dire' [lthesel*kfreq], 'wnsz' [lthesel*kfreq,2] and 'sllf' [lthesel*kfreq,2]:
				thisdata$B_T1_i <- data$bpt1( list(  dira=thisdata$fishdirT[,2,],  dire=thisdata$fishdirT[,3,],  wnsz=cbind(rep(data$wavenumber[thisdata$thesebeams]*data$rad1[thisdata$thesebeams,1],each=data$lthesel),rep(data$wavenumber[thisdata$thesebeams]*data$rad1[thisdata$thesebeams,2],each=data$lthesel)),  sllf=cbind(rep(data$sllf[thisdata$thesebeams,1],each=data$lthesel),rep(data$sllf[thisdata$thesebeams,2],each=data$lthesel)) ) , appr=circularPiston_appr )
				}
			else if(data$pbp1=="circularPiston_ellipticRadius"){
				# Use as input the 'dira' [lthesel*kfreq], 'dire' [lthesel*kfreq] and 'wnsz' [lthesel*kfreq,2]:
				thisdata$B_T1_i <- data$bpt1( list(  dira=thisdata$fishdirT[,2,],  dire=thisdata$fishdirT[,3,],  wnsz=cbind(rep(data$wavenumber[thisdata$thesebeams]*data$rad1[thisdata$thesebeams,1],each=data$lthesel),rep(data$wavenumber[thisdata$thesebeams]*data$rad1[thisdata$thesebeams,2],each=data$lthesel)) ) , appr=circularPiston_appr )
				}
			else if(data$pbp1=="circularPiston"){
				# Use as input the 'dire' [lthesel*kfreq] and 'wnsz' [lthesel*kfreq,2]:
				thisdata$B_T1_i <- data$bpt1( list(  dire=thisdata$fishdirT[,3,],  wnsz=cbind(rep(data$wavenumber[thisdata$thesebeams]*data$rad1[thisdata$thesebeams,1],each=data$lthesel),rep(data$wavenumber[thisdata$thesebeams]*data$rad1[thisdata$thesebeams,2],each=data$lthesel)) ) , appr=circularPiston_appr )
				}
			else{
				# Use as input the 'dira' [lthesel*kfreq], 'dire' [lthesel*kfreq] and 'indi' [lthesel*kfreq,2]:
				thisdata$B_T1_i <- data$bpt1( list(  dira=thisdata$fishdirT[,2,],  dire=thisdata$fishdirT[,3,]  ,indi=rep(data$indi[thisdata$thesebeams],each=data$lthesel)  ) )
				}
			
			# Beam pattern on reception at the trasducer as a matrix of dimension c(data$lthesel,data$kfreq), where kfreq is the number of beams of the same frequency:
			if(data$equalbp_em_re){
				thisdata$B_T2_i <- thisdata$B_T1_i
				}
			else{	
				if(data$pbp2=="circularPiston_ellipticRadius_sidelobefit"){
					# Use as input the 'dira' [lthesel*kfreq], 'dire' [lthesel*kfreq], 'wnsz' [lthesel*kfreq,2] and 'sllf' [lthesel*kfreq,2]:
					thisdata$B_T2_i <- data$bpt2( list(  dira=thisdata$fishdirT[,2,],  dire=thisdata$fishdirT[,3,],  wnsz=cbind(rep(data$wavenumber[thisdata$thesebeams]*data$rad2[thisdata$thesebeams,1],each=data$lthesel),rep(data$wavenumber[thisdata$thesebeams]*data$rad2[thisdata$thesebeams,2],each=data$lthesel)),  sllf=cbind(rep(data$sllf[thisdata$thesebeams,1],each=data$lthesel),rep(data$sllf[thisdata$thesebeams,2],each=data$lthesel)) ) , appr=circularPiston_appr )
					}
				else if(data$pbp2=="circularPiston_ellipticRadius"){
					# Use as input the 'dira' [lthesel*kfreq], 'dire' [lthesel*kfreq] and 'wnsz' [lthesel*kfreq,2]:
					thisdata$B_T2_i <- data$bpt2( list(  dira=thisdata$fishdirT[,2,],  dire=thisdata$fishdirT[,3,],  wnsz=cbind(rep(data$wavenumber[thisdata$thesebeams]*data$rad2[thisdata$thesebeams,1],each=data$lthesel),rep(data$wavenumber[thisdata$thesebeams]*data$rad2[thisdata$thesebeams,2],each=data$lthesel)) ) , appr=circularPiston_appr )
					}
				else if(data$pbp2=="circularPiston"){
					# Use as input the 'dire' [lthesel*kfreq] and 'wnsz' [lthesel*kfreq,2]:
					thisdata$B_T2_i <- data$bpt2( list(  dire=thisdata$fishdirT[,3,],  wnsz=cbind(rep(data$wavenumber[thisdata$thesebeams]*data$rad2[thisdata$thesebeams,1],each=data$lthesel),rep(data$wavenumber[thisdata$thesebeams]*data$rad2[thisdata$thesebeams,2],each=data$lthesel)) ) , appr=circularPiston_appr )
					}
				else{
					# Use as input the 'dira' [lthesel*kfreq], 'dire' [lthesel*kfreq] and 'indi' [lthesel*kfreq,2]:
					thisdata$B_T2_i <- data$bpt2( list(  dira=thisdata$fishdirT[,2,],  dire=thisdata$fishdirT[,3,]  ,indi=rep(data$indi[thisdata$thesebeams],each=data$lthesel)  ) )
					}
				}
			
			# Adding up the contribution to the currently insonified fish from the beams given by 'thisdata$thesebeams' (beams having equal frequency):
			# Take the mean of the beam pattern of all beams of the current frequency for each fish:
			dim(thisdata$B_T1_i) <- c(data$lthesel, data$kfreq)
			thisdata$B_T1_i <- rowSums(thisdata$B_T1_i)
			
			# Beam pattern of the fish:
			# If each fish is exposed to different frequencies, there are different beam patterns from the same fish for the different frequencies:
			thisdata$B_L_i <- data$bptf(list(dira=data$transducerposL[data$thesel,2], dire=data$transducerposL[data$thesel,3], wnsz=data$uniquek[i]*data$lenl[data$thesel]))
			# Dump 'B_L':
			dumpdata[2] <- sum(dumpdata[2], thisdata$B_L_i, na.rm=TRUE)
			
			# chi:
			if(any(data$sigma0mode %in% 3:4)){
				if(is.null(data$ssil)){
					if(length(data$grsf)==1){
						thisdata$chi <- suppressWarnings(approx(data$grsf, data$ssif, data$uniquek[i] * data$lenl[data$thesel], method="constant",rule=2)$y)
						}
					else{
						thisdata$chi <- suppressWarnings(approx(data$grsf, data$ssif, data$uniquek[i] * data$lenl[data$thesel], method="linear",rule=2)$y)
						}
					}
				else{
					thisdata$chi <- data$ssil[,i]
					}
				# Dump 'chi':
				dumpdata[6] <- sum(dumpdata[6], thisdata$chi, na.rm=TRUE)
				}
				
			# Absorption, assuming that the transducer and receiver is located on the same place (the same physical object). The expression differs from the one in the documentation due to the definition of the absorption coefficient used by Simrad, which is absr=alpha/10, where 'alpha' is the absorption coefficient used in the documentation:
			thisdata$thisetaa <- 10^(outer(-0.2*data$transducerposL[data$thesel,1] , data$absr[thisdata$thesebeams],"*"))
			# Dump 'etaa':
			dumpdata[3] <- sum(dumpdata[3], thisdata$thisetaa, na.rm=TRUE)
			dumpdata[7] <- data$sigma0mode
			
			### Add to the output:
			# data$sigma0mode==1: If the optimal backscattering cross section 'sgbs' (sigma_bs) is present, no relation to target size is specified:
			# data$sigma0mode==2: If the coefficient 'epss' linking 'sgbs' to fish size is present, it may be a function of frequency, or simply a numeric vector:
			# data$sigma0mode==3: If the optimal acoustic cross sectional area 'acca' (A_0) is present, no relation to target size is specified:
			# data$sigma0mode==4: If the coefficient 'epsl' linking 'acca' to fish size is present, it may be a function of frequency, or simply a numeric vector:
			# See also echoIBM.oneping.oneschool() and echoIBM.default.oneschool() for explaination of data$sigma0mode:
			if(data$sigma0mode==1){
				# Simply multiply the varialbes 'fish', 'etaj', 'B_T1_i', 'B_L_i', 'B_T2_i' and 'thisetaa':
				withoutetaj <- data$fish[data$thesel] * thisdata$B_T1_i * thisdata$B_L_i * thisdata$B_T2_i * thisdata$thisetaa
				nsig[thisdata$thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
				out[thisdata$thesebeams] <- colSums(withoutetaj * thisdata$etaj)
				outnext[thisdata$thesebeams] <- colSums(withoutetaj * (1-thisdata$etaj))
				}
				
			else if(data$sigma0mode==2){
				# Include the echo ability 'epss' as a function of frequency:
				if(is.function(data$epss)){
					thisdata$epss <- data$epss(data$uniquek[i] * data$asps / (2*pi))
					withoutetaj <- data$fish[data$thesel] * thisdata$B_T1_i * thisdata$B_L_i * thisdata$B_T2_i * thisdata$thisetaa * thisdata$epss
					nsig[thisdata$thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
					out[thisdata$thesebeams] <- colSums(withoutetaj * thisdata$etaj)
					outnext[thisdata$thesebeams] <- colSums(withoutetaj * (1-thisdata$etaj))
					}
				# Include the echo ability 'epss' as a vector:
				else{
					thisdata$epss <- data$epss[data$thesel]
					withoutetaj <- data$fish[data$thesel] * thisdata$B_T1_i * thisdata$B_L_i * thisdata$B_T2_i * thisdata$thisetaa * thisdata$epss
					nsig[thisdata$thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
					out[thisdata$thesebeams] <- colSums(withoutetaj * thisdata$etaj)
					outnext[thisdata$thesebeams] <- colSums(withoutetaj * (1-thisdata$etaj))
					}
				# Dump 'esps':
				dumpdata[4] <- sum(dumpdata[4], thisdata$epss, na.rm=TRUE)
				}
				
			else if(data$sigma0mode==3){
				# Multiply the varialbes 'fish', 'etaj', 'B_T1_i', 'B_L_i', 'B_T2_i' and 'thisetaa', but also include the spherical surface integral 'ssil':
				withoutetaj <- data$fish[data$thesel] * thisdata$B_T1_i * thisdata$B_L_i * thisdata$B_T2_i * thisdata$thisetaa / thisdata$chi
				nsig[thisdata$thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
				out[thisdata$thesebeams] <- colSums(withoutetaj * thisdata$etaj)
				outnext[thisdata$thesebeams] <- colSums(withoutetaj * (1-thisdata$etaj))
				}
				
			else if(data$sigma0mode==4){
				# Include the echo ability 'epsl' as a function of frequency:
				if(is.function(data$epsl)){
					thisdata$epsl <- data$epsl(data$uniquek[i]*data$asps/(2*pi))
					withoutetaj <- data$fish[data$thesel] * thisdata$B_T1_i * thisdata$B_L_i * thisdata$B_T2_i * thisdata$thisetaa * thisdata$epsl / thisdata$chi
					nsig[thisdata$thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
					out[thisdata$thesebeams] <- colSums(withoutetaj * thisdata$etaj)
					outnext[thisdata$thesebeams] <- colSums(withoutetaj * (1-thisdata$etaj))
					}
				# Include the echo ability 'epsl' as a vector:
				else{
					thisdata$epsl <- data$epsl[data$thesel]
					withoutetaj <- data$fish[data$thesel] * thisdata$B_T1_i * thisdata$B_L_i * thisdata$B_T2_i * thisdata$thisetaa * thisdata$epsl / thisdata$chi
					nsig[thisdata$thesebeams] <- apply(withoutetaj, 2, nsignif, na.rm=TRUE)
					out[thisdata$thesebeams] <- colSums(withoutetaj * thisdata$etaj)
					outnext[thisdata$thesebeams] <- colSums(withoutetaj * (1-thisdata$etaj))
					}
				# Dump 'espl':
				dumpdata[5] <- sum(dumpdata[5], thisdata$epsl, na.rm=TRUE)
				}
			dumpdata[8] <- sum(dumpdata[8], withoutetaj, na.rm=TRUE)
			} # End of for(i in seq_along(data$uniquek)).
		# Dump:
		# Detailed information about the school ("etaj","B_L","etaa","epss","epsl","chi"):
		dumpdata[1] <- dumpdata[1] / data$lthesel
		dumpdata[2] <- dumpdata[2] / data$lthesel / length(data$wavenumber)
		dumpdata[3] <- dumpdata[3] / data$lthesel / length(data$wavenumber)
		dumpdata[4] <- dumpdata[4] / length(data$wavenumber)
		dumpdata[5] <- dumpdata[5] / length(data$wavenumber)
		dumpdata[6] <- dumpdata[6] / data$lthesel / length(data$wavenumber)
		dumpdata[8] <- dumpdata[8] / length(data$wavenumber)
		} # End of if(data$lthesel>0).
	
	# Return:
	list(sv=out, svnext=outnext, nsig=nsig, dumpdata=dumpdata)
	##################################################
	##################################################
	}
