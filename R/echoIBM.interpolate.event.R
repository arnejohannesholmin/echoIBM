#*********************************************
#*********************************************
#' Runs through segmentation data and interpolates between pings, calculating the probability that voxels in neighboring pings, corresponding to each voxel included in the segmentation mask at a time step, are occupied by the school as well.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param t  is either the number of the ping to be treated, as listed from 1 to the number of pings in the event, or the time point given as a string "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF". Only one time step alowed!
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event.
#' @param adds  is an optional list of variables overriding the variables in 'data'.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. Currently implemented are "EK60", "ME70", "MS70" and "SH80"/"SX80"/"SH90"/"SX90" (may be given in lover case)
#' @param type  is the type of interpolation:
#' @param bw  is the smooting parameter of the interpolation, given as the standard deviation for Gaussian kernel and the base of the exponential decays (bw^x).
#' @param thr  is the threshold for the interpolated probabilities og school, defining the interpolated segmentation masks.
#' @param kernelthr  is the threshold below which the kernel is set to 0.
#' @param ind  is a list of indexes, as typed into the [] of an array, where 0 and NULL denotes all indexes.
#' @param range  is a list of elements with names matching names of 'data', specifying the range of the corresponding elements. If range is given as a list of 3 elements, xyz
#' @param subset  is a numeric or logical vector/expression indicating elements or rows to keep. Missing values are taken as false, and subset=0 or subset=NULL indicates no subsetting.
#' @param segpar  is a list of elements named "bwGp", "lsth"/"rlst", "usth"/"rust", or "sgth"/"sgt0" specifying the parameters of the segmentation data to read.
#' @param spar  is used in smooth.spline() when the interpolated probabilities of school are smoothed.
#' @param TOV  is the time offset of the vessel information, discovered by Holmin and Korneliussen in 2013. The default is found in "/Applications/echoIBM/Documentation/Error in yaw MS70/Error in yaw MS70.R".
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR cm.school event.path find_voxels.TSD read.event rotate3D subset_TSD
#' @importFrom TSD NAs strff write.TSD zeros
#' @importFrom stats dnorm smooth.spline
#'
#' @export
#' @rdname echoIBM.interpolate.event
#'
echoIBM.interpolate.event<-function(event=1, t=1, cruise=2009116, adds=NULL, esnm="MS70", type=c("le","lg","e","g","emax"), bw=0.9, thr=NULL, kernelthr=0.1, ind=list(-(1:30),NULL), range=list(), subset=NULL, segpar=list(), spar=NULL, TOV=0){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-03-22 - Clean version.
	########### DESCRIPTION: ###########
	# Runs through segmentation data and interpolates between pings, calculating the probability that voxels in neighboring pings, corresponding to each voxel included in the segmentation mask at a time step, are occupied by the school as well.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---event--- is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
	# ---t--- is either the number of the ping to be treated, as listed from 1 to the number of pings in the event, or the time point given as a string "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF". Only one time step alowed!
	# ---cruise--- is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event.
	# ---adds--- is an optional list of variables overriding the variables in 'data'.
	# ---esnm--- is the name of the acoustical instrument, given as a four character string. Currently implemented are "EK60", "ME70", "MS70" and "SH80"/"SX80"/"SH90"/"SX90" (may be given in lover case)
	# ---type--- is the type of interpolation:
	#			'le' - exponential decay kernel forward and backward, with the maximum of these taken for each ping, applied in the logarithmic domain.
	#			'lg' -  Gaussian kernel, applied in the logarithmic domain.
	#			'e' -  exponential decay kernel forward and backward, with the maximum of these taken for each ping.
	#			'g' - Gaussian kernel.
	#			'emax' - exponential decay values forward and backward, with the maximum of these taken for each ping.
	# ---bw--- is the smooting parameter of the interpolation, given as the standard deviation for Gaussian kernel and the base of the exponential decays (bw^x).
	# ---thr--- is the threshold for the interpolated probabilities og school, defining the interpolated segmentation masks.
	# ---kernelthr--- is the threshold below which the kernel is set to 0.
	# ---ind--- is a list of indexes, as typed into the [] of an array, where 0 and NULL denotes all indexes.
	# ---range--- is a list of elements with names matching names of 'data', specifying the range of the corresponding elements. If range is given as a list of 3 elements, xyz
	# ---subset--- is a numeric or logical vector/expression indicating elements or rows to keep. Missing values are taken as false, and subset=0 or subset=NULL indicates no subsetting.
	# ---segpar--- is a list of elements named "bwGp", "lsth"/"rlst", "usth"/"rust", or "sgth"/"sgt0" specifying the parameters of the segmentation data to read.
	# ---spar--- is used in smooth.spline() when the interpolated probabilities of school are smoothed.
	# ---TOV--- is the time offset of the vessel information, discovered by Holmin and Korneliussen in 2013. The default is found in "/Applications/echoIBM/Documentation/Error in yaw MS70/Error in yaw MS70.R".
					
	
	##################################################
	##################################################
	##### Preparation #####
	# Read the segmentation data for the entire sequence of pings:
	data = c(adds,read.event(event=event, t=t, cruise=cruise, adds=adds, esnm=esnm, segpar=segpar, var=c("time","sgsc","pr0S","psze","beams","vessel"), TOV=TOV))
	t = data$indt
	numt = length(t)
	
	### Define the kernel function used in the extrapolation: ###
	if(!any(type %in% c("g", "lg", "e", "le", "emax"))){
		stop("'type' must be one of \"g\", \"lg\", \"e\", \"le\", or \"emax\"")
	}
	# Gaussian kernel:
	if(strff("g", type[1]) || strff("lg", type[1])){
		kernelfun = function(x){
			dnorm(x, 0, bw)
		}
	}
	# Exponential decay kernel:
	else if(strff("e", type[1]) || strff("le", type[1])){
		kernelfun = function(x){
			bw^abs(x)
		}
	}
	# Logarithmic kernel requires a 'pSchool' lower than 1 to prevent value=Inf:
	if(strff("l", type[1])){
		pSchool = 0.99999
		value = log(1-pSchool)
	}
	else{
		value = 1
	}
	
	# Define the threshold used to obtain the final interpolated segmentation masks 'sgsi':
	ndetect = 0
	ndetect_default = 2
	if(length(thr)==0){
		ndetect = ndetect_default
	}
	else if(thr>=1 && thr%%1==0){
		ndetect = thr
	}
	else if(thr<=0 || thr>=1){
		warning("'thr' must be a positive integer or a value between 0 and 1. Chosen to correspond to discarding two consecutive identified voxels")
		ndetect = ndetect_default
	}
	if(ndetect>0){
		# We wish to avoid a voxel identified falsely at 'ndetect' consecutive pings to register in the interpolated segmentation mask for a two second difference in time between pings:
		nsecBetweenPings = 2
		maxutimdiff = log(kernelthr)/log(bw)
		k = bw^seq(0, maxutimdiff, nsecBetweenPings)
		thr = sum(k[seq_len(ndetect)])/sum(k)
	}
	
	
	##### Execution #####
	### Read voxel positions for each ping and extract the voxels included in the segmentation masks: ###
	data$psxx = vector("list", numt)
	data$psyx = vector("list", numt)
	data$pszx = vector("list", numt)
	data$cmxS = NAs(numt)
	data$cmyS = NAs(numt)
	data$cmzS = NAs(numt)
	data$scxS = NAs(numt)
	data$scyS = NAs(numt)
	data$sczS = NAs(numt)
	data$utim = NAs(numt)
	
	# Loop through the pings:
	cat("Extracting voxels included in the segmentation masks...\n")
	for(i in seq_len(numt)){
		cat(t[i], " (", i, " of ", numt, ")\n", sep="")
		# Read the voxel position data:
		vox = read.event(event=event, t=t[i], cruise=cruise, esnm=esnm, var=c("psxx","psyx","pszx","time"), TOV=TOV)
		data$utim[i] = vox$utim
		# Include the segmentation of the current ping in the 'subset' used in subset_TSD():
		if(is.logical(subset)){
			thissubset = intersect(which(subset), data$sgsc[[i]])
		}
		else if(length(subset)>0){
			thissubset = intersect(subset, data$sgsc[[i]])
		}
		else{
			thissubset = data$sgsc[[i]]
		}
		# Apply subset_TSD():
		if(length(thissubset)>0){
			vox = subset_TSD(vox, ind=ind, range=range, subset=thissubset, ind.out=TRUE)
			# Extract the subsets:
			data$psxx[[i]] = vox$psxx
			data$psyx[[i]] = vox$psyx
			data$pszx[[i]] = vox$pszx
			data$sgsc[[i]] = vox$subs[[1]]
		}
		# Get the center of mass of the values in the segmentation mask:
		if(length(data$psxx[[i]])>0){
			cmS = cm.school(vox)
			data$cmxS[i] = cmS[1]
			data$cmyS[i] = cmS[2]
			data$cmzS[i] = cmS[3]
		}
	}
	# Clean up memory:
	rm(vox)
	gc()
	
	# Spline smooth the points in each dimension:
	notna = !is.na(data$cmxS)
	data$scxS[notna] = smooth.spline(data$utim[notna], data$cmxS[notna], spar=spar)$y
	data$scyS[notna] = smooth.spline(data$utim[notna], data$cmyS[notna], spar=spar)$y
	data$sczS[notna] = smooth.spline(data$utim[notna], data$cmzS[notna], spar=spar)$y
	
	# Declare the interpolated segmentation data:
	data$sgsI = vector("list", numt)
	data$sgsi = vector("list", numt)
	data$psis = vector("list", numt)
	data$sgsIForeward = vector("list", numt)
	data$psisForeward = vector("list", numt)
	data$sgsIBackward = vector("list", numt)
	data$psisBackward = vector("list", numt)
	
	# Declare the kernels sums used in the smoothing:
	kernelsum = zeros(numt)
	kernelsumForeward = zeros(numt)
	kernelsumBackward = zeros(numt)
	
	
	### Extrapolation (including the reference ping, indicated by this1=i): ###
	cat("Extrapolation...\n")
	for(i in seq_len(numt)){
		cat(t[i], " (", i, " of ", numt, ")\n", sep="")
		if(length(data$psxx[[i]])>0){
			# Define the kernel values in reference to the current ping (dependent on the time difference between the following pings and the current ping):
			kernel = kernelfun(data$utim-data$utim[i])
			# Apply the threshold on 'kernel', avoiding too many included pings:
			add = which(kernel>kernelthr)
			addForeward = add[add>=i]
			addBackward = add[add<=i]
			
			# Get the kernel sums:
			kernelsum[i] = sum(kernel[add])
			kernelsumForeward[i] = sum(kernel[addForeward])
			kernelsumBackward[i] = sum(kernel[addBackward])
						
			# Loop through the pings to extrapolate to:
			for(j in add){
				# Rotate and translate the voxels of the current ping into the coordinate system of ping 'j' (also subtract the change in position of the school, specified in the input data):
				rthetaphi = rotate3D(
					x=	cbind(
							data$psxx[[i]], data$psyx[[i]], data$pszx[[i]]) - 
						matrix(c(
							data$psxv[j] + data$scxS[j] - data$scxS[i], 
							data$psyv[j] + data$scyS[j] - data$scyS[i], 
							data$pszv[j] + data$sczS[j] - data$sczS[i] + data$psze), 
							ncol=3, nrow=length(data$psxx[[i]]), byrow=TRUE), 
					by="z", 
					ang=data$rtzv[j], 
					drop.out=FALSE, 
					sph.out=TRUE)[,,1]
				
				if(length(dim(rthetaphi))==0){
					dim(rthetaphi) = c(1,3)
				}
				# Get the voxel indices in ping 'j':
				extrapolatedVoxels = find_voxels.TSD(rthetaphi, data[c("lenb","dire","dira","sint","asps")])
				# Remove missing voxels:
				extrapolatedVoxels = extrapolatedVoxels[!is.na(extrapolatedVoxels)]
				# Remove duplicate elements (voxels which by chance are located to the same voxel in the ping to which they are extrapolated):
				extrapolatedVoxels = unique(extrapolatedVoxels)
				
				# Insert the extrapolated voxels for the exponential decay kernels, which include the present ping in the backward and the foreward direction:
				if(strff("e", type[1]) || strff("le", type[1])){
					if(j<=i){
						data$sgsIForeward[[j]] = c(data$sgsIForeward[[j]], extrapolatedVoxels)
						data$psisForeward[[j]] = c(data$psisForeward[[j]], value*rep(kernel[j], length(extrapolatedVoxels)))
					}
					if(j>=i){
						data$sgsIBackward[[j]] = c(data$sgsIBackward[[j]], extrapolatedVoxels)
						data$psisBackward[[j]] = c(data$psisBackward[[j]], value*rep(kernel[j], length(extrapolatedVoxels)))
					}
				}
				# Insert the extrapolated voxels for the Gaussian kernels:
				else{
					data$sgsI[[j]] = c(data$sgsI[[j]], extrapolatedVoxels)
					data$psis[[j]] = c(data$psis[[j]], value*rep(kernel[j], length(extrapolatedVoxels)))
				}	
			}
		}
	}	
	
	
	### Apply the smoothing: ###
	cat("Smoothing...\n")
	for(i in seq_len(numt)){
		cat(t[i], " (", i, " of ", numt, ")\n", sep="")
		# Gaussian kernel, should be summed for each voxel of each ping:
		if(strff("g", type[1]) || strff("lg", type[1])){
			if(length(data$psis[[i]])>0){
				b = by(data$psis[[i]], data$sgsI[[i]], sum, na.rm=TRUE)
				data$psis[[i]] = c(b)/kernelsum[i]
				data$sgsI[[i]] = as.numeric(names(b))
			}
		}
		# Logarithmic exponential decay kernel, should be summed for each voxel of each ping separately for foreward and backward extrapolation, and the maximum of these should be used:
		else if(strff("e", type[1]) || strff("le", type[1])){
			# Treat each direciton separately (dividing by the kernel sum), and take the maximum for each voxel afterwards:
			if(length(data$psisForeward[[i]])>0){
				# Get the smoothed values forward:
				b = by(data$psisForeward[[i]], data$sgsIForeward[[i]], sum, na.rm=TRUE)
				data$psisForeward[[i]] = c(b)/kernelsumForeward[i]
				data$sgsIForeward[[i]] = as.numeric(names(b))
			}
			if(length(data$psisBackward[[i]])>0){
				# Get the smoothed values backward:
				b = by(data$psisBackward[[i]], data$sgsIBackward[[i]], sum, na.rm=TRUE)
				data$psisBackward[[i]] = c(b)/kernelsumBackward[i]
				data$sgsIBackward[[i]] = as.numeric(names(b))
			}
			# Get the maximum for each voxel:
			if(length(data$psisForeward[[i]])>0 || length(data$psisBackward[[i]])>0){
				b = by(c(data$psisForeward[[i]], data$psisBackward[[i]]), c(data$sgsIForeward[[i]], data$sgsIBackward[[i]]), max, na.rm=TRUE)
				data$psis[[i]] = c(b)
				data$sgsI[[i]] = as.numeric(names(b))
			}
		}
		# Exponential decay of the closest "School"-voxel, which is obtained by taking the maximum for each voxel of each ping:
		else if(strff("emax", type[1])){
			if(length(data$psisForeward[[i]])>0 && length(data$psisBackward[[i]])>0){
				b = by(c(data$psisForeward[[i]], data$psisBackward[[i]]),  c(data$sgsIForeward[[i]], data$sgsIBackward[[i]]), max, na.rm=TRUE)
				data$psis[[i]] = c(b)
				data$sgsI[[i]] = as.numeric(names(b))
			}
		}
	}
	# Convert to linear values for logarithmic kernels:
	if(strff("l", type[1])){
		for(i in seq_len(numt)){
			if(length(data$psis[[i]])>0){
				data$psis[[i]] = 1-exp(data$psis[[i]])
			}
		}
	}
	# Threshold the interpolated probabilities:
	for(i in seq_len(numt)){
		thissgsi = data$sgsI[[i]][data$psis[[i]]>thr]
		if(length(thissgsi)>0){
			data$sgsi[[i]] = thissgsi
		}
	}
	
	
	##### Output #####
	# Write the interpolated segmentation masks to file:
	towrite = c(data[c("utim","indt","sgsi","sgsI","psis","cmxS","cmyS","cmzS","scxS","scyS","sczS")], list(krnI=kernel,typI=type,sctI=thr,smtI=bw,krtI=kernelthr,pSch=pSchool))
	event = event.path(event=event, cruise=cruise, esnm=esnm, get.name=TRUE)
	if(length(grep("interpolated.seg", data$segfile[1], fixed=TRUE))>0){
		filename = paste(substr(data$segfile[1], 1, nchar(data$segfile[1])-17), "_interpolated.seg", sep="")
	}
	else{
		filename = paste(substr(data$segfile[1], 1, nchar(data$segfile[1])-4), "_interpolated.seg", sep="")
	}
	#filename = file.path(event[1],"seg",paste(event[2],"interpolated.seg",sep="_"))
	write.TSD(towrite, filename)
	
	# Return:
	c(towrite, data["sgsc"])
	##################################################
	##################################################
}
