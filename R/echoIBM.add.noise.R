#*********************************************
#*********************************************
#' Adds noise to an array of sv-values, given the noise data present in 'data'.
#'
#' @param sv  is an array dimension [lenb,numb] (or a vector of the same length) of sv-values.
#' @param indt  is a the time point index, required to pick out the correct seed.
#' @param noise  is a vector of character strings of length 2, specifying which types of noise to apply to the data:
#' @param data  is the list of beams inputs as returned from read.TSD. For the treatment of noise the following variables are required:
#' @param parlist  is a list of input parameters to the noise generation method specified by 'noise'. Use echoIBM_rexp_defaults() to generate defaul values. Must include a seed vector of at least the same length as the maximum indt value to be used.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR echoIBM.pdns
#' @importFrom TSD arr.ind2ind dim_all strff zeros
#' @importFrom stats rexp
#'
#' @export
#' @rdname echoIBM.add.noise
#'
echoIBM.add.noise<-function(sv, indt, noise=c("bg","pn","ms"), data=NULL, parlist=list()){
	
	############### LOG: ###############
	# Start: 2010-10-16 - First version.
	# Update: 2011-11-26 - Added the option of adding correlated exponential noise (noise="cex").
	# Update: 2011-12-06 - Applied the exponential perturbation for each of "background noise", "close-range noise" and "signal", and summing afterwards, instead of adding all data first and applying the exponential perturbation afterwards.
	# Update: 2012-02-14 - Added the method echoIBM_rexp_MultSines() for generating correlated exponential/Barakat values.
	# Update: 2014-04-05 - Added "cpp" and "cap" in 'noise'.
	# Update: 2014-04-06 - Removed "nr_fun" in 'noise'.
	# Last: 2014-10-14 - Changed to draw exponential values only at the end.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(sum(nchar(noise))==0){
		sv
		}
	# Function used for adding exponential independent or correlated noise:
	addex <- function(sv, noise, indt, parlist){
		# Pick out the seed of the current time step:
		parlist$currentseed <- parlist$seed[indt]
		
		### # Treat the seed so that if NULL is given draw random seeds, if 'indt' is given pick out the indt'th seed, and otherwise use the first element: 
		### if(any(c("ms", "ex", "bk") %in% noise)){
		### 	if(is.null(parlist$seed)){
		### 		parlist$seed <- runif(1, 0, 1e9)
		### 		}
		### 	# If indt is given (as a single integer), pick out the indt'th seed value:
		### 	else if(length(indt)){
		### 		parlist$seed <- parlist$seed[indt]
		### 	}
		### }
		
		# Add correlated exponential noise by multiple sines:
		if(any(c("ms") %in% noise)){
			# Define the overlap array if the phase angles are different for each time step:
			if("pn" %in% noise && is.character(parlist$olpn) && strff("p", parlist$olpn)){
				parlist$olpn <- echoIBM.olpt(data, indt=indt)$olpt
				dim(parlist$olpn) <- c(dim(parlist$olpn)[1:2], parlist$luqf, parlist$nuqf)
				}
			# Set the dimension of the data to a matrix of beams along the columns, and sum over rows and select the rows that have sum > 0, and expand to one row on either side of the selected rows. Do this to save CPU time:
			dim(sv) <- c(parlist$J, parlist$luqf * parlist$nuqf)
			valid <- which(rowSums(sv) > 0)
			if(length(valid)>0){
				valid <- seq(max(1,min(valid)-1), min(parlist$J,max(valid,na.rm=TRUE)+1))
				# Redefine the number of samples along beams:
				parlist$J <- length(valid)
				#thisexp <- echoIBM_rexp_MultSines(parlist)$noise^(1/1.6) / gamma(1+1/1.6)
				thisexp <- echoIBM_rexp_MultSines(parlist)$noise
				#Apply the exponent 'Wshp' to the exponentially distributed variables to obtain Weibull distributed variables:
				#thisexp[,451:475] <- thisexp[,451:475]^(1/1.5) / gamma(1+1/1.5)
				#thisexp <- thisexp^parlist$Wshp[seq_len(parlist$J),] / gamma(1+1/parlist$Wshp[seq_len(parlist$J),])
				
				# Insert the exponentially distributed values in the correct rows:
				sv[valid,] <- sv[valid,]*thisexp
				}
			
			#thisexp <- echoIBM_rexp_MultSines(parlist)$noise
			#sv <- sv*thisexp
			}
		# Add correlated exponential noise by rearrangement method:
		#else if(any(c("acex","aex","cex") %in% noise)){
		#	if("ex" %in% noise){
		#		warning("Correlated exponential noise simulated (not independent exponential noise as indicated by noise = \"ex\")")
		#		}
		#	parlist$rate <- 1/c(sv)
		#	sv <- echoIBM_rexp_Rearr(parlist)
		#	}
		# Add independent exponential noise:
		else if("ex" %in% noise){
			# Set the seed:
			set.seed(parlist$currentseed)
			# Draw from the exponential distribution:
			sv <- rexp(max(data$lenb)*data$numb, 1/c(sv))
			dim(sv) <- c(max(data$lenb), data$numb)
			}
		else if("bk" %in% noise){
			# Set the seed:
			set.seed(parlist$currentseed)
			# Draw from the exponential distribution:
			sv <- rexp(max(data$lenb)*data$numb, 1/c(sv))
			# Apply the Barakat distribution:
			if(length(parlist$Brkt)>0){
				sv[parlist$Brkt] <- rBarakatI(length(parlist$Brkt), nsig=parlist$nsig, magn=sv[parlist$Brkt]/parlist$nsig, N=parlist$N, max.memory=parlist$max.memory)
				}
			dim(sv) <- c(max(data$lenb), data$numb)
			}
		sv
		}
	# If parlist is not given, it is extracted from the method of generating defaults:
	if(length(parlist)==0){
		warning("'parlist' not given and was defaulted using echoIBM_rexp_defaults()")
		parlist <- echoIBM_rexp_defaults(noise=noise, data=data, parlist=parlist)
		}
			
	
	########## Execution ##########
	# Define the expected value of the noise in each voxel:
	beta <- zeros(dim_all(sv))
	
	### (1) Add background noise: ###
	if("bg" %in% noise){
		if(is.null(data$bgns)){
			warning("Background noise data 'bgns' missing. To generate noise data use get.bgns()")
			}
		else{
			# Add background noise to the output:
			beta <- beta + matrix(data$bgns, nrow=max(data$lenb), ncol=data$numb, byrow=TRUE)
			}
		}	
	
	
	### (2) Add periodic noise: ###
	if("pn" %in% noise){
		pdns.present <- all(length(data$sint)>0, length(data$lenb)>0, length(data$pns1)>0, length(data$pns2)>0, length(data$pns3)>0, length(data$acfq)>0, length(data$harm)>0) & (length(data$bgns)>0 | length(data$numb)>0)
		if(!pdns.present){
			warning("Periodic noise data missing. To generate noise data use get.bgns.MS70()")
			}
		else{
			# Add periodic noise to the output:
			if(is.character(parlist$olpn) && strff("p", parlist$olpn)){
				beta <- beta + echoIBM.pdns(data, indt=indt, TVG=FALSE)$pdns
				}
			else if(strff("ek60",data$esnm[1])){
				beta <- beta + echoIBM.pdns(data, indt=NULL, TVG=FALSE)$pdns
				}
			else{
				beta <- beta + echoIBM.pdns(data, indt=NULL, TVG=FALSE)$pdns
				}
			}
		}	
	
	
	### (3) Add high intensity noise: ###
	if("hi" %in% noise){
		# Error if one of 'numb', 'freq' or "lenb" is missing:
		if(any(length(data$numb)==0, length(data$lenb)==0)){
			stop("'numb' and 'lenb' must be present for high intensity noise to be added")
			}
		# Add high intensity noise if present:
		if(all(length(data$hins)>0, length(data$hini)>0)){
			warning("High intensity noise data missing. To generate noise data use get.hins.MS70()")
			}
		else{
			# If 'indt' is missing, choose the first time step that has high intensity noise:
			if(length(indt)==0){
				indt <- min(data$hini[,3])
				}
			theseindt <- data$hini[,3]==indt
			theseind <- arr.ind2ind(data$hini[theseindt,1:2], c(max(data$lenb),data$numb))
			# Add high intensity noise to the output:
			beta[theseind] <- beta[theseind]+data$hins[theseindt]
			}
		}	
	
	
	### (4) Add close-range passive noise: ###
	if(any(c("nr", "np", "cp") %in% noise)){
		if(is.null(data$nr0p)){
			warning("Close-range passive noise data 'nr0p' missing. To generate noise data use get.nrns.MS70()")
			}
		else{
			# Close-range noise may be given as a vector of one value for each sampling interval along the beams (with length equal to the length of the longest beam):
			if(length(dim(data$nr0p))!=2){
				data$nr0p <- matrix(data$nr0p, nrow=max(data$lenb), ncol=data$numb)
				}
			}
		}
	# Or as a matrix/array of time steps:
	else if("cpp" %in% noise){
		if(is.null(data$cnpM)){
			warning("Pingwise close-range passive noise data 'cnpM' missing.")
			}
		else{
			# Close-range noise may be given for each time step as a vector of one value for each sampling interval along the beams (with length equal to the length of the longest beam):
			dimcnpM <- dim(data$cnpM)
			# When 'cnpM' is given, it is assumed that the last dimenstion holds time steps:
			if(length(dimcnpM)==2){
				data$nr0p <- matrix(data$cnpM[,(indt-1) %% dimcnpM[2] + 1], nrow=max(data$lenb), ncol=data$numb)
				}
			else if(length(dimcnpM)==3){
				data$nr0p <- data$cnpM[,,(indt-1) %% dimcnpM[2] + 1]
				}
			else{
				warning("Wrong dimension of 'cnpM'. Should have dimension (length of beams, number of beams, number of pings) or (length of beams, number of pings)")
				}
			}
		}
	
	# Adjust the length of the near-range noise matrix and add close-range to the output::
	if(any(c("nr", "np", "cp", "cpp") %in% noise) && length(data$nr0p)>0){
		if(nrow(data$nr0p)<max(data$lenb)){
			data$nr0p <- rbind(data$nr0p, zeros(max(data$lenb)-nrow(data$nr0p), data$numb))
			}
		else{
			data$nr0p <- data$nr0p[1:max(data$lenb),]
			}
		beta <- beta+data$nr0p
		}
	
	
	### (5) Add close-range passive noise: ###
	if(any(c("na", "ca") %in% noise)){
		if(is.null(data$nr0a)){
			warning("Close-range active noise data 'nr0a' missing. To generate noise data use get.nrns.MS70()")
			}
		else{
			# Close-range noise may be given as a vector of one value for each sampling interval along the beams (with length equal to the length of the longest beam):
			if(length(dim(data$nr0a))!=2){
				data$nr0a <- matrix(data$nr0a, nrow=max(data$lenb), ncol=data$numb)
				}
			}
		}
	# Or as a matrix/array of time steps:
	else if("cap" %in% noise){
		if(is.null(data$cnaM)){
			warning("Pingwise close-range active noise data 'cnaM' missing.")
			}
		else{
			# Close-range noise may be given for each time step as a vector of one value for each sampling interval along the beams (with length equal to the length of the longest beam):
			# When 'cnaM' is given, it is assumed that the last dimenstion holds time steps:
			dimcnaM <- dim(data$cnaM)
			if(length(dimcnaM)==2){
				data$nr0a <- matrix(data$cnaM[,(indt-1) %% dimcnaM[2] + 1], nrow=max(data$lenb), ncol=data$numb)
				}
			else if(length(dimcnpM)==3){
				data$nr0a <- data$cnaM[,,(indt-1) %% dimcnaM[2] + 1]
				}
			else{
				warning("Wrong dimension of 'cnaM'. Should have dimension (length of beams, number of beams, number of pings) or (length of beams, number of pings)")
				}
			}
		}
	
	# Adjust the length of the near-range noise matrix and add close-range to the output:
	if(any(c("na", "ca", "cap") %in% noise) && length(data$nr0a)>0){
		if(NROW(data$nr0a)<max(data$lenb)){
			data$nr0a <- rbind(data$nr0a, zeros(max(data$lenb)-nrow(data$nr0a), data$numb))
			}
		else{
			data$nr0a <- data$nr0a[1:max(data$lenb),]
			}
		beta <- beta+data$nr0a
		}
	### (6) Do not run the Barakat PDF on noise: ###
	else if("bk" %in% noise){
		parlist$Brkt <- NULL
		}
	
	sv <- addex(sv=beta+sv, noise=noise, indt=indt, parlist=c(list(olp0=parlist$input_olpn),parlist))
	
	
	########## Output ##########
	sv
	##################################################
	##################################################
	}
