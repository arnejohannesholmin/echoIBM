#*********************************************
#*********************************************
#' Calls the function rexp_MultSines(), which generates correlated vectors of autocorrelated exponential/Barakat variables. Only one time step is generated.
#'
#' @param parlist  is a list of input parameters:
#' @param data  is a list of beam configuration data, used to extract defaults if 'parlist' is not given.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all ftim2utim read.TSD read.TSDs utim.TSD zeros
#' @importFrom stats runif median
#'
#' @export
#' @rdname echoIBM_rexp_MultSines
#'
echoIBM_rexp_MultSines<-function(parlist=list(), data=list(), Npre=1000, shape=1){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-10-16 - First version.
	# Update: 2012-02-17 - Added the option 'buffer'.
	# Update: 2012-02-17 - Removed the option 'buffer'.
	# Update: 2012-02-28 - Changed to using 'parlist' holding the parameters.
	# Update: 2013-09-24 - Added the option of drawing from pre-generated sets of correlated exponential noise using parlist$pre==TRUE.
	# Last: 2013-11-04 - Added 'olp0', which is the string naming the type of correlation structure to use (one of "c" - constant for all voxels, "b" - variable between beams, "v" - variable between volxels, and "p" variable between voxels and pings, normally used when the phase of the periodic noise is given).
	########### DESCRIPTION: ###########
	# Calls the function rexp_MultSines(), which generates correlated vectors of autocorrelated exponential/Barakat variables. Only one time step is generated.
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---parlist--- is a list of input parameters:
	#		'J' is the number of sample intervals along the beams.
	#		'nuqf' is the number of fans comprising the system (typically >1 for a rectangular grid of beams).
	#		'luqf' is the number of beams in each fan.
	#		'L' is the number of targets pr voxel. Low values results in increasingly non-rayleigh presure and non-exponential intensity (which is what is returned).
	#		'N' is the number of sample points pr voxel. High values increases stability and accuracy in the simulated noise, but increases time demand of the function.
	#		'P' is the number of periods pr voxel, which should be a high enough number, but not too high, to ensure good performance of the funciton.
	#		'w' is a vector of the lengths of the sine waves in units of the time intervals constituting the voxels, one for each fan (constant in MS70, variable in EK60). The autocorrelation will depend on this value.
	#		'olpn' is eihter a character string specifying the type of default correlation structure, or a vector, matrix or array of correlations, representing constant correlation for all voxels, variable correlation between beams, and variable correlation for all voxels. If given as a string, the following values are accepted: "const", for constant correlation for all voxels, "beams", for variable correlation between beams, and "voxels", for variable correlation for all voxels. Abbreviations are accepted. This function assumes the default phase angles. Variable phase angles are obtained using the special function echoIBM.olpt().	
	#		'seed' is the seed vector, one seed for each time step.
	#		'currentseed' is the seed of the current time step.
	#		'pre' can be set to TRUE to draw from pre-generated sets of Npre pings of correlated exponential noise, in which case 'seed' defines the time step translated into 1:Npre.
	# ---data--- is a list of beam configuration data, used to extract defaults if 'parlist' is not given.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	### # Defaults (see "Test_of_rexp_MultSines.R"). This combination is as far as the analysis established, the least time consuming combination resulting in sufficiently exponential variables:
	### #if(length(parlist)==0){
	### #	parlist=echoIBM_rexp_defaults(noise="ms",data=data,parlist=parlist)
	### #	}
	### if(!any(c("L", "N", "P", "olp0", "olpn", "J", "nuqf", "luqf", "rate", "prob", "l", "rho", "w", "wC", "C", "Cind", "buffer") %in% names(parlist))){
	### 	parlist <- c(echoIBM_rexp_defaults(noise="ms", data=data, parlist=parlist), parlist)
	### 	}
	# Error if not all required variables are present:
	#if(!all(c("J", "C", "wC", "w", "rate", "prob", "l", "seed", "Cind", "buffer") %in% names(parlist))){
	#	stop("The list of parameters 'parlist' must contain the following variables: \"J\", \"C\", \"wC\", \"w\", \"rate\", \"prob\", \"l\", \"seed\", \"Cind\", \"buffer\"")
	#	}

	########## Execution ##########
	if(isTRUE(parlist$pre)){
		preGenerated_rexp_MultSines <- file.path(echoIBM_datasets_directory(), "Resources", "rexp_MultSines")
		# Check for compatibility with existing noise data:
		# Read the overlap values for the voxels, and the number of fans and beams:
		rexp_MultSinesDir <- preGenerated_rexp_MultSines
		if(file.exists(rexp_MultSinesDir)){
			rexp_MultSinesFiles <- list.files(rexp_MultSinesDir, full.names=TRUE)
			parlist$parmatch <- logical(length(rexp_MultSinesFiles))
			
			rexp_MultSinesVar_list <- as.data.frame(zeros(length(rexp_MultSinesFiles), 4))
			names(rexp_MultSinesVar_list) <- c("olp0", "nuqf", "luqf", "npsw")
			
			for(i in seq_along(rexp_MultSinesFiles)){
				rexp_MultSinesVar <- read.TSD(rexp_MultSinesFiles[i], var=c("olp0", "nuqf", "luqf", "npsw"))
				rexp_MultSinesVar_list[i,] <- sapply(rexp_MultSinesVar[c("olp0", "nuqf", "luqf", "npsw")], median)
				# Match with the parameters in 'parlist':
				if(length(rexp_MultSinesVar$olp0)==length(parlist$olp0)){
					matcholp0 <- tolower(rexp_MultSinesVar$olp0) == tolower(parlist$olp0)
					}
				else{
					matcholp0 <- FALSE
					}
				if(length(rexp_MultSinesVar$nuqf)==length(parlist$nuqf)){
					matchnuqf <- !any(abs(rexp_MultSinesVar$nuqf - parlist$nuqf)>1e-4)
					}
				else{
					matchnuqf <- FALSE
					}
				if(length(rexp_MultSinesVar$luqf)==length(parlist$luqf)){
					matchluqf <- !any(abs(rexp_MultSinesVar$luqf - parlist$luqf)>1e-4)
					}
				else{
					matchluqf <- FALSE
					}
				if(length(rexp_MultSinesVar$npsw)==length(parlist$w)){
					matchnpsw <- !any(abs(rexp_MultSinesVar$npsw - parlist$w)>1e-4)
					}
				else{
					matchnpsw <- FALSE
					}
				parlist$parmatch[i] <- all(matcholp0, matchnuqf, matchluqf, matchnpsw)
				}
			
			# If several files are matched, select the one closest in time to the reference time point:
			if(sum(parlist$parmatch)>1){
				parlist$utim <- ftim2utim(parlist$ftim)
				# Read the time points of the files:
				utim <- sapply(utim.TSD(read.TSDs(rexp_MultSinesFiles[parlist$parmatch], var="utim", merge=TRUE)), "[[",1)
				if(length(parlist$ftim)==0){
					parlist$parmatch <- rexp_MultSinesFiles[parlist$parmatch][which.min(utim)]
					warning(paste("Several files matching the parameters specified by 'parlist' (specifically, 'olp0', 'nuqf', 'luqf', 'w') were found. In that case the element 'ftim' should be included in parlist, to select the file with time point closest to this value (given as \"yyyymmddHHMMSS\"). Since 'ftim' is not present, the earliest file is chosen:", parlist$parmatch))
					}
				else{
					parlist$parmatch <- rexp_MultSinesFiles[parlist$parmatch][which.min(abs(utim-parlist$utim))]
					}		
				}
			else{
				parlist$parmatch <- rexp_MultSinesFiles[parlist$parmatch]
				}
			
			if(length(parlist$parmatch)==0){
				warning(paste0(
					"'parlist$pre' set to TRUE, but no files containing pre-generated noise data found in the directory ", 
					rexp_MultSinesDir, 
					" for the specified parameter settings: olp0 = ", 
					parlist$olp0, 
					" (", 
					paste(rexp_MultSinesVar_list$olp0, collapse=", "), 
					"), nuqf = ", 
					parlist$nuqf, 
					" (", 
					paste(rexp_MultSinesVar_list$nuqf, collapse=", "), 
					"), luqf = ", 
					parlist$luqf, 
					" (", 
					paste(rexp_MultSinesVar_list$luqf, collapse=", "), 
					"), w = ", 
					parlist$w, 
					" (", 
					paste(rexp_MultSinesVar_list$npsw, collapse=", "), 
					")"))
				parlist$pre <- FALSE
				}
			else{
				#cat("Pre-generated correlated exponentially distributed values read and applied to the simulations\n")
				}
			}
		else{
			warning(paste("Pre-generated correlated exponentially distributed values were requested by parlist$pre=TRUE, but were not found in the directory", preGenerated_rexp_MultSines))
			}
		}
	# If there was a match between the parameters in 'parlist' and the parameters of the pre-generated noise, select a ping in the pre-generated noise:
	if(isTRUE(parlist$pre)){
		if(length(parlist$parmatch)>1){
			warning(paste0("There are ", length(parlist$parmatch), " files containing pre-generated noise data for the specified parameter settings. The first chosen"))
			}
		# Read the pre-generated noise at the time step given by 'parlist$seed' translated into indices 1, ..., Npre:
		noise <- read.TSD(parlist$parmatch[1], var="vbsc", t=(round(parlist$currentseed-1) %% Npre) + 1)$vbsc
		if(parlist$J>nrow(noise)){
			warning("The pre-generated noise values are too short for the specified range of the sonar. Maximum range is 900 m")
			}
		noise <- noise[seq_len(parlist$J),]
		}
	else{
		noise <- zeros(parlist$J, parlist$luqf, parlist$nuqf)
		# Expand the seed to one value per unique frequency, since the multisines method considers each each frequency separately: 
		set.seed(parlist$currentseed)
		parlist$currentseed <- runif(parlist$nuqf, 0, 1e9)
		
		
		# Different dimensions of the overlap information require different array indexes:
		# Constant overlap for all beams and all voxels. Dimension of 'olpn' expanded in rexp_MultSines() to a matrix [length(olpn), number of beams]:
		if(length(dim_all(parlist$olpn))==1){
			for(i in seq_len(parlist$nuqf)){
				noise[,,i] <- rexp_MultSines(J=parlist$J, I=parlist$luqf, L=parlist$L, N=parlist$N, P=parlist$P, w=parlist$w[i], olpn=parlist$olpn, seed=parlist$currentseed[i], shape=shape)
				}
			}
		# Constant overlap for all beams and all voxels within fans, but variable between fans. Dimension of 'olpn' [number of adjacet beams correlated, number of beams]:
		else if(length(dim_all(parlist$olpn))==2){
			for(i in seq_len(parlist$nuqf)){
				noise[,,i] <- rexp_MultSines(J=parlist$J, I=parlist$luqf, L=parlist$L, N=parlist$N, P=parlist$P, w=parlist$w[i], olpn=parlist$olpn[,i], seed=parlist$currentseed[i], shape=shape)
				}
			}
		# Variable overlap for all beams and constant overlap for all voxels within beams. Dimension of 'olpn' [number of adjacet beams correlated, number of voxels along the beams, number of beams]:
		else if(length(dim_all(parlist$olpn))==3){
			for(i in seq_len(parlist$nuqf)){
				noise[,,i] <- rexp_MultSines(J=parlist$J, I=parlist$luqf, L=parlist$L, N=parlist$N, P=parlist$P, w=parlist$w[i], olpn=parlist$olpn[,,i], seed=parlist$currentseed[i], shape=shape)
				}
			}
		# Variable overlap for all voxels, only used with the periodic noise. Dimension of 'olpn': [number of adjacent correlated beams, number of voxels along the beams, number of beams in each fan of constant frequency, number of fans of constant frequency]:
		else if(length(dim_all(parlist$olpn))==4){
			for(i in seq_len(parlist$nuqf)){
				noise[,,i] <- rexp_MultSines(J=parlist$J, I=parlist$luqf, L=parlist$L, N=parlist$N, P=parlist$P, w=parlist$w[i], olpn=parlist$olpn[,seq_len(parlist$J),,i], seed=parlist$currentseed[i],shape=shape)
				}
			}
		else{
			stop("Only overlap matrices up to 4 dimensions are supported")
			}
		dim(noise)=c(parlist$J, parlist$luqf*parlist$nuqf)
		}
			
	
	########## Output ##########
	list(noise=noise, parlist=parlist)
	##################################################
	##################################################
	}
