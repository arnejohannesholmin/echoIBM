#*********************************************
#*********************************************
#' Adds noise to an array of sv-values for the ping 'indt', given the noise data present in 'data'. 
#'
#' @param sv  is an array dimension [lenb,numb] (or a vector of the same length) of sv-values.
#' @param indt  is a the time point index, required to pick out the correct seed.
#' @param noisetypes  is a vector of character strings of length 2, specifying which types of noise to apply to the data:
#' @param data  is the list of beams inputs as returned from read.TSD. For the treatment of noise the following variables are required:
#' @param parlist  is a list of input parameters to the noise generation method specified by 'noisetypes'. Use echoIBM_rexp_defaults() to generate defaul values. Must include a seed vector of at least the same length as the maximum indt value to be used.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR echoIBM.pdns
#' @importFrom TSD arr.ind2ind dim_all strff zeros
#' @importFrom stats rexp
#' @importFrom utils head
#'
#' @export
#' @rdname echoIBM.add.noise
#'
echoIBM.add.noise <- function(sv, indt, noisetypes=c("bg","pn","ms"), data=NULL, parlist=list()){
	
	############### LOG: ###############
	# Start: 2010-10-16 - First version.
	# Update: 2011-11-26 - Added the option of adding correlated exponential noise (noisetypes="cex").
	# Update: 2011-12-06 - Applied the exponential perturbation for each of "background noise", "close-range noise" and "signal", and summing afterwards, instead of adding all data first and applying the exponential perturbation afterwards.
	# Update: 2012-02-14 - Added the method echoIBM_rexp_MultSines() for generating correlated exponential/Barakat values.
	# Update: 2014-04-05 - Added "cpp" and "cap" in 'noisetypes'.
	# Update: 2014-04-06 - Removed "nr_fun" in 'noisetypes'.
	# Last: 2014-10-14 - Changed to draw exponential values only at the end.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Modify lenb and numb:
	if(length(data$lenb)){
		data$lenb <- max(data$lenb)
	}
	if(length(data$numb)){
		data$numb <- data$numb[1]
	}
	
	if(sum(nchar(noisetypes))==0){
		sv
	}
	# Function used for adding exponential independent or correlated noise:
	addex <- function(sv, noisetypes, indt, parlist){
		# Pick out the seed of the current time step:
		parlist$currentseed <- parlist$seed[indt]
		
		### # Treat the seed so that if NULL is given draw random seeds, if 'indt' is given pick out the indt'th seed, and otherwise use the first element: 
		### if(any(c("ms", "ex", "bk") %in% noisetypes)){
		### 	if(is.null(parlist$seed)){
		### 		parlist$seed <- runif(1, 0, 1e9)
		### 	}
		### 	# If indt is given (as a single integer), pick out the indt'th seed value:
		### 	else if(length(indt)){
		### 		parlist$seed <- parlist$seed[indt]
		### }
		### }
		
		# Add correlated exponential noise by multiple sines:
		if(any(c("ms") %in% noisetypes)){
		    # Define the overlap array if the phase angles are different for each time step:
			if("pn" %in% noisetypes && is.character(parlist$olpn) && strff("p", parlist$olpn)){
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
				if(length(thisexp)) {
				    sv[valid,] <- sv[valid,] * thisexp
				}
			}
			
			#thisexp <- echoIBM_rexp_MultSines(parlist)$noise
			#sv <- sv*thisexp
		}
		# Add correlated exponential noise by rearrangement method:
		#else if(any(c("acex","aex","cex") %in% noisetypes)){
		#	if("ex" %in% noisetypes){
		#		warning("Correlated exponential noise simulated (not independent exponential noise as indicated by noisetypes = \"ex\")")
		#	}
		#	parlist$rate <- 1/c(sv)
		#	sv <- echoIBM_rexp_Rearr(parlist)
		#}
		# Add independent exponential noise:
		else if("ex" %in% noisetypes){
			# Set the seed:
			set.seed(parlist$currentseed)
			# Draw from the exponential distribution:
			sv <- rexp(data$lenb * data$numb, 1/c(sv))
			dim(sv) <- c(data$lenb, data$numb)
		}
		else if("bk" %in% noisetypes){
			# Set the seed:
			set.seed(parlist$currentseed)
			# Draw from the exponential distribution:
			sv <- rexp(data$lenb * data$numb, 1/c(sv))
			# Apply the Barakat distribution:
			if(length(parlist$Brkt)>0){
				sv[parlist$Brkt] <- rBarakatI(length(parlist$Brkt), nsig=parlist$nsig, magn=sv[parlist$Brkt]/parlist$nsig, N=parlist$N, max.memory=parlist$max.memory)
			}
			dim(sv) <- c(data$lenb, data$numb)
		}
		sv
	}
	
	# Function for adjusting the length of a noise matrix to the lengths of the beams:
	adjustTolenb <- function(data, var){
		if(NROW(data[[var]]) < data$lenb){
			rbind(data[[var]], zeros(data$lenb - nrow(data[[var]]), data$numb))
		}
		else{
			data[[var]][seq_len(data$lenb),, drop=FALSE]
		}
	}
	
	# Function for generating one particular noise type at one time step:
	getNoiseOnePing <- function(data, var){
		dimnoise <- dim(data[[var]])
		if(is.function(data[[var]])){
			data[[var]] <- data[[var]](indt)
		}
		if(is.null(data[[var]])){
			warning(paste0("Noise data '", var, "' missing."))
		}
		else{
			# Add noise to the output. Support for three dimensional arrays with time along the third dimension, two dimensional arrays specifying noise for all voxels, or one dimensional specifying one value per beam:
			if(length(dimnoise)==3){
				# Use the modulus of time step index 'indt':
				data[[var]] <- data[[var]][,,(indt-1) %% dimnoise[3] + 1]
			}
			else if(length(dimnoise)<2){
				data[[var]] <- matrix(data[[var]], nrow=data$lenb, ncol=data$numb, byrow=TRUE)
			}
		}
		adjustTolenb(data, var)
	}
	
	
	
	# If parlist is not given, it is extracted from the method of generating defaults:
	if(length(parlist)==0){
		warning("'parlist' not given and was defaulted using echoIBM_rexp_defaults()")
		parlist <- echoIBM_rexp_defaults(noisetypes=noisetypes, data=data, parlist=parlist)
	}
			
	
	########## Execution ##########
	# Define the expected value of the noise in each voxel:
	beta <- zeros(dim_all(sv))
	
	### (1) Add background noise: ###
	if("bg" %in% noisetypes){
		beta <- beta + c(getNoiseOnePing(data, "bgns"))
	}
	
	
	### (2) Add periodic noise: ###
	if("pn" %in% noisetypes){
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
	if("hi" %in% noisetypes){
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
			theseind <- arr.ind2ind(data$hini[theseindt,1:2], c(data$lenb, data$numb))
			# Add high intensity noise to the output:
			beta[theseind] <- beta[theseind]+data$hins[theseindt]
		}
	}	
	
	
	### (4) Add close-range noise: ###
	if(any(c("nr", "np", "cp", "cpp", "na", "ca", "cap") %in% noisetypes)){
		var <- head(intersect(c("nr0a", "nr0p", "cnaM", "cnpM"), names(data)), 1)
		beta <- beta + getNoiseOnePing(data, var)
	}
	
	### (6) Do not run the Barakat PDF on noise: ###
	else if("bk" %in% noisetypes){
		parlist$Brkt <- NULL
	}
	
	sv <- addex(sv=beta+sv, noisetypes=noisetypes, indt=indt, parlist=c(list(olp0=parlist$input_olpn),parlist))
	
	
	########## Output ##########
	sv
	##################################################
	##################################################
}
