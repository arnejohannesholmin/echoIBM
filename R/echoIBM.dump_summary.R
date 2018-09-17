#*********************************************
#*********************************************
#' Writes a summary of a variable used by echoIBM() to the dumpfile.
#'
#' @param data  is a variable of which to write summary to the dumpfile.
#' @param dumpfile  is the path to the file to which to write the information.
#' @param type  is the type of variable to be dumped.
#' @param append  is used in write.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all info.TSD
#' @importFrom stats quantile
#' @importFrom utils write.table
#'
#' @export
#' @rdname echoIBM.dump_summary
#'
echoIBM.dump_summary <- function(data, dumpfile, type=c("beams", "ctd", "vessel", "staticschool", "dynschool", "freqschool", "noise"), append=TRUE){
	
	############### LOG: ###############
	# Start: 2011-08-12 - Clean version.
	# Update: 2012-07-26 - Added periodic noise variables.
	# Update: 2012-11-29 - Expanded the dump.
	# Update: 2011-12-01 - Added the parameters 'scls' and 'rand.sel' used to reduce CPU time.
	# Last: 2013-09-27 - Cleaned up and linked to the TSD-description file.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	dumpMaxLength <- 10
	beamsnames <- c("esnm", "numb", "lenb", "freq", "wavenumber", "absr", "psze", "sint", "dira", "dire", "bwtl", "bwtt", "bwtx", "bwty", "pbp1", "pbp2", "rad1", "rad2", "bpt1", "bpt2", "sllf", "cali", "eqba", "gain", "tpow", "sacr")
	staticschoolnames <- c("pbpf", "graf", "gref", "grsf", "grif", "dbpf", "ebpf", "bptf", "size", "lenl", "tilt", "gamw", "gaml", "obln", "ssil", "grff", "epsl_table", "epss_table", "epsl", "epss", "acca", "ssif", "sgbs", "spow", "kLseq", "chivector", "sigma0mode")
	ctdnames1 <- c("lonc", "latc", "rho0", "gacc", "hpr0", "asps")
	ctdnames2 <- c("pszc", "ihpr", "temp", "slty", "isps")
	vesselnames <- c("latv", "lonv", "psxv", "psyv", "pszv", "rtxv", "rtyv", "rtzv", "ispv", "lat0", "lon0")
	dynamicschoolnames <- c("psxf", "psyf", "pszf", "vlxf", "vlyf", "vlzf", "rtxf", "rtyf", "rtzf", "transducerposL", "size", "lenl", "etaC", "etar4", "etaomega", "fish")
	freqschoolnames <- c("etaj", "B_L", "B_T1", "B_T2", "etaa", "epss", "sigma0mode", "withoutetaj", "Nl")
	noisenames <- c("J", "nuqf", "luqf", "seed", "L", "N", "P", "w", "input_olpn", "input_olps", "olpn", "olps", "pdns_scale", "rate", "prob", "l", "rho", "w", "wC", "C", "Cind", "buffer", "acfq", "badb", "harm", "pns1", "pns2", "pns3", "bgns")
	rand.sel_scalenames <- c("rand.sel", "scls")
		
	nonTSD <- c(
	"wavenumber",		"Wavenumbers of the sonar",
	"Nl",				"Number of targets",
	"epsl_table",		"Grid table of the factor linking frequency <grff> and frequency dependent optimal acoustic cross sectional area <acca>",
	"epss_table",		"Grid table of the factor linking frequency <grff> and frequency dependent optimal backscattering cross section <sgbs>",
	"kLseq",			"Grid of the product of wavenumber and length of the acoustic model of the target corresponding to spherical surface integral of the beam pattern",
	"chivector",		"Spherical surface integral of the beam pattern corresponding to kLseq",
	"transducerposL",	"Spherical position of sonar in target coordinate system",
	"etaC",				"Depth compression of the targets",
	"etar4",			"Geometrical spreading of the acoustic wave (r^-4)",
	"etaomega",			"Orientation factor for the acoustic cross sectional receiving area",
	"fish",				"The value transfered to echoIBM.oneping.j1.R or echoIBM.oneping.j2.R",
	"etaj",				"Mean and max radial weighting",
	"B_L",				"Mean and max beam pattern of the fish in the direciton of the transducer",
	"B_T1",				"Mean and max beam pattern of the transducer in the direction of the targets at emission",
	"B_T2",				"Mean and max beam pattern of the transducer in the direction of the targets at reception",
	"etaa",				"Mean and max acoustical absorption factor",
	"epss", 			"Mean and max factor linking frequency <grff> and frequency dependent optimal backscattering cross section <sgbs>",
	"sigma0mode",		"Mode of backscatter: 1 - The optimal backscattering cross section 'sgbs' (sigma_bs) is present, no relation to target size is specified, 2 - the coefficient 'epss' linking 'sgbs' to fish size is present, it may be a function of frequency, or simply a numeric vector",
	"J",				"The number of sample intervals along the beams",
	"seed",				"The seed of the simulation. If 'seed' is NULL, a random seed is used",
	"L",				"(ms) The number of targets pr voxel. Low values results in increasingly non-rayleigh presure and non-exponential intensity",
	"N",				"(ms) The number of sample points pr voxel. High values increases stability and accuracy in the simulated noise, but increases time demand of the function",
	"P",				"(ms) The number of periods pr voxel",
	"w",				"(ms) The length of the sine waves in units of the time intervals constituting the voxels. The autocorrelation will depend on this value",
	"scale",			"The expectation of the exponential variables, estimated in \"Test_of_rexp_MultSines.R\" in the \"extdata\" directory of the echoIBM package",
	"rate",				"(cex) The parameter in the exponential distribution",
	"prob",				"(cex) The probability values assigned to each of the consecutive values sorted by their proximity to the prevoius value. If w = 3 and prob = [0.7,0.5,0,5] the closest value will be assigned to the current position with probability 0.7, and if this fails, the second closest value will be assigned to the current position with probability 0.5, and so on. Used to tune the impact of the rearrangement and reducing autocorrelation",
	"l",				"(cex) The number of values included when randomly rearranging the values not assigned to the current position. We may include more values than 'w' too reduce the effect of nevative autocorrelation at lags ≈ 5",
	"rho",				"(cex) A vector of correlation values defining the dependence of a vector to its neighbor vectors on both sides. The reference vector is in the middle of 'rho' (usually having value 1) and the length of 'rho' must be odd",
	"w",				"(cex) The number of consecutive values including the present value. If w=2, the method choses between the current and the next value depending on the proximity to the prevoius values",
	"wC",				"(cex) The number of significant vectors correlated to each vector, used to deduce Cind. Must be an odd number",
	"C",				"(cex) The number of consecutive values including the present value. If w=2, the method choses between the current and the next value depending on the proximity to the prevoius values",
	"Cind",				"(cex) Optional matrix of index values for the significant correlation values of C for each row",
	"buffer",			"(cex) A value used in the function rexp_Rearr.R to avoid a strange error in the c++ function of the same name",
	"rand.sel",			"The fraction of the targets used in the simulation",
	"scls",				"Scaling for the backscatter of the targets",
	"input_olpn",		"Character giving the type of overlap between voxels for the noise (c = constant for all voxels, b = variable between beams, v = variable between voxels, p = variable between voxels and pings)",
	"input_olps", 		"Character giving the type of overlap between voxels for the signal (c = constant for all voxels, b = variable between beams, v = variable between voxels, p = variable between voxels and pings)", 
	"withoutetaj",		"Mean and max of the product fish * B_T1 * B_L * B_T2 * etaa * epss from echoIBM.oneping.oneschool.onesample()")
	
	nonTSD <- as.data.frame(matrix(nonTSD, ncol=2))
	names(nonTSD) <- c("name", "description")
	
	
	getDescription <- function(name, large=FALSE, full=FALSE){
		out <- info.TSD(name, info.out=FALSE, clean=TRUE)[1 + full]
		#if(is.na(out[1])){
			m <- match(name, nonTSD$name)
			out[m] <- nonTSD$description[m]
		#}
		if(large){
			out <- paste(out," (range)")
		}
		paste0("\n\"",name,"\" (", out, "):")
	}
	
	
	########## Execution and output ##########
	# Sonar variables:
	if(identical(tolower(type[1]), "beams")){
		write("\n\n\n##### ACOUSTICAL PROPERTIES OF THE SONAR: #####", dumpfile, append=TRUE)
		for(j in seq_along(beamsnames)){
			this <- data[[beamsnames[j]]]
			if(length(this)>0){
				if(is.list(this)){
					warning(paste0("The variable ", beamsnames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(beamsnames[j], large=length(this)>dumpMaxLength), dumpfile, append=TRUE)
				
				# Write the function body for funciton elements:
				if(is.function(this)){
					write(format(this), dumpfile, append=append)
				}
				# If 'data' has length > dumpMaxLength write the range of the value:
				else if(length(this)<=dumpMaxLength){
					write(this, dumpfile, append=TRUE, ncolumns=dumpMaxLength)
				}
				else{
					write(range(this), dumpfile, append=TRUE, ncolumns=5, sep="\t")
				}
			}
		}
	}
	# Static school variables
	else if(identical(tolower(type[1]), "staticschool")){
		write("\n\n\n##### ACOUSTIC MODEL OF THE TARGETS: #####", dumpfile, append=TRUE)
		for(j in seq_along(staticschoolnames)){
			this <- data[[staticschoolnames[j]]]
			if(length(this)>0){
				if(is.list(this)){
					warning(paste0("The variable ", staticschoolnames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(staticschoolnames[j]), dumpfile, append=TRUE)
				# Write the function body for funciton elements:
				if(is.function(this)){
					write(format(this), dumpfile, append=append)
				}
				# Else if 'data' has 2 or 3 columns, write summary and histogram of each column:
				else if(NCOL(this) %in% 2:3){
					for(i in seq_len(NCOL(this))){
						h <- hist_simple(this[,i], breaks=10)
						s <- c(quantile(this[,i]), mean(this[,i]))
						s <- s[c(1:3, 6, 4:5)]
						names(s) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
						write(names(s), dumpfile, append=append, ncolumns=10, sep="\t")
						write(s, dumpfile, append=append, ncolumns=10, sep="\t")
						write(c("breaks", "counts"), dumpfile, append=append, ncolumns=10, sep="\t")
						write(t(cbind(h$breaks, c(h$counts, rep("", length(h$breaks)-length(h$counts))))), dumpfile, append=append, ncolumns=2, sep="\t")
					}
				}
				# Else if the length of the variable is < dumpMaxLength, write the value:
				else if(length(this)<=dumpMaxLength){
					write(this, dumpfile, append=append, ncolumns=10)
				}
				# Else collapse and write summary and histogram:
				else{
					if(length(dim(this))>1){
						h <- hist_simple(c(this), breaks=10)
						s <- c(quantile(c(this)), mean(c(this)))
						s <- s[c(1:3, 6, 4:5)]
						names(s) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
					}
					else{
						h <- hist_simple(this, breaks=10)
						s <- c(quantile(this), mean(this))
						s <- s[c(1:3, 6, 4:5)]
						names(s) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
					}
					write(names(s), dumpfile, append=append, ncolumns=10, sep="\t")
					write(s, dumpfile, append=append, ncolumns=10, sep="\t")
					write(c("breaks", "counts"), dumpfile, append=append, ncolumns=10, sep="\t")
					write(t(cbind(h$breaks, c(h$counts, ""))), dumpfile, append=append, ncolumns=2, sep="\t")
				}
			}
		}
	}
	# Vessel variables:
	else if(identical(tolower(type[1]), "vessel")){
		write("\n\n\n##### DYNAMIC VESSEL INFORMATION: #####", dumpfile, append=TRUE)
		# If the length of the variable is < dumpMaxLength, write the value:
		for(j in seq_along(vesselnames)){
			this <- data[[vesselnames[j]]]
			if(length(this)>0){
				if(is.list(this)){
					warning(paste0("The variable ", vesselnames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(vesselnames[j]), dumpfile, append=TRUE)
				write(this, dumpfile, append=append)
			}
		}
	}
		
	# Dynamic school variables:
	else if(identical(tolower(type[1]), "dynschool")){
		write("\n\n\n##### DYNAMIC TARGET INFORMATION: #####", dumpfile, append=TRUE)
		for(j in seq_along(dynamicschoolnames)){
			this <- data[[dynamicschoolnames[j]]]
			if(length(this)>0){
				if(is.list(this)){
					warning(paste0("The variable ", dynamicschoolnames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(dynamicschoolnames[j]), dumpfile, append=TRUE)
				# If 'data' has 2 or 3 columns, write summary and histogram of each column:
				if(NCOL(this) %in% 2:3){
					for(i in seq_len(NCOL(this))){
						h <- hist_simple(this[,i], breaks=10)
						s <- c(quantile(this[,i]), mean(this[,i]))
						s <- s[c(1:3, 6, 4:5)]
						names(s) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
						write(paste0("Column ", i, ":"), dumpfile, append=append)
						write(names(s), dumpfile, append=append, ncolumns=10, sep="\t")
						write(s, dumpfile, append=append, ncolumns=10, sep="\t")
						write(c("breaks", "counts"), dumpfile, append=append, ncolumns=10, sep="\t")
						write(t(cbind(h$breaks, c(h$counts, ""))), dumpfile, append=append, ncolumns=2, sep="\t")
					}
				}
				# Else if the length of the variable is < dumpMaxLength, write the value:
				else if(length(this)<=dumpMaxLength){
					write(this, dumpfile, append=append, ncolumns=10)
				}
				# Else collapse and write summary and histogram:
				else{
					if(length(dim(this))>1){
						h <- hist_simple(c(this), breaks=10)
						s <- c(quantile(c(this)), mean(c(this)))
						s <- s[c(1:3, 6, 4:5)]
						names(s) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
					}
					else{
						h <- hist_simple(this, breaks=10)
						s <- c(quantile(this), mean(this))
						s <- s[c(1:3, 6, 4:5)]
						names(s) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
					}
					write(names(s), dumpfile, append=append, ncolumns=10, sep="\t")
					write(s, dumpfile, append=append, ncolumns=10, sep="\t")
					write(c("breaks", "counts"), dumpfile, append=append, ncolumns=10, sep="\t")
					write(t(cbind(h$breaks, c(h$counts, rep("", length(h$breaks)-length(h$counts))))), dumpfile, append=append, ncolumns=2, sep="\t")
				}
			}
		}
	}
	
	# Dynamic school variables:
	else if(identical(tolower(type[1]), "freqschool")){
		write(paste0("\n\n##### FREQUENCY DEPENDENT TARGET INFORMATION: #####"), dumpfile, append=TRUE)
		for(j in seq_along(freqschoolnames)){
			# Get the current variable:
			this <- data[[freqschoolnames[j]]]
			if(length(this)>0){
				# Apply default significant digits:
				this <- signif(this, 6)
				if(is.list(this)){
					warning(paste0("The variable ", freqschoolnames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(freqschoolnames[j]), dumpfile, append=TRUE)
				if(length(this)==2){
					write(paste(c("Mean:", "Max:"), this), dumpfile, append=append, ncolumns=10, sep="\t")
				}
				else{
					write(paste("Sum:", this), dumpfile, append=append, ncolumns=10, sep="\t")
				}
				
			}
		}
	}
	
	# CTD variables:
	else if(identical(tolower(type[1]), "ctd")){
		write("\n\n\n##### CTD INFORMATION: #####", dumpfile, append=TRUE)
		# Write the ctd data uncompressed:
		for(j in seq_along(ctdnames1)){
			if(length(data[[ctdnames1[j]]])>0){
				write(getDescription(ctdnames1[j]), dumpfile, append=TRUE)
				write(data[[ctdnames1[j]]], dumpfile, append=append)
			}
		}
		ctdnames2present <- sapply(data[ctdnames2], function(x) length(x)>0)
		if(any(ctdnames2present)){
			for(j in seq_along(ctdnames2[ctdnames2present])){
				write(getDescription(ctdnames2[ctdnames2present][j]), dumpfile, append=TRUE)
			}
			suppressWarnings(write.table(round(as.data.frame(data[ctdnames2[ctdnames2present]]), digits=2), dumpfile, append=TRUE, row.names=FALSE, sep="\t"))
		}
	}
	
	# Noise variables:
	else if(identical(tolower(type[1]), "noise")){
		write("\n\n\n##### PARAMETERS IN THE GENERATION OF NOISE AND RANDOMNESS: #####", dumpfile, append=TRUE)
		for(j in seq_along(noisenames)){
			this <- data[[noisenames[j]]]
			if(length(this)>0){
				if(is.list(this)){
					warning(paste0("The variable ", noisenames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(noisenames[j]), dumpfile, append=TRUE)
				if(noisenames[j]=="C"){
					s <- min(25, NCOL(this))
					write(this[s,s], dumpfile, append=append, ncolumns=s)
				}
				else if(noisenames[j]=="olpn"){
					if(length(dim_all(this)) %in% c(1,2)){
						write(this, dumpfile, append=append, ncolumns=10)
					}
					else if(length(dim_all(this))==3){
						olp_fans <- c(1,12,13,14,15,16,17)
						for(i in olp_fans){
							write(paste("Displaying fan",i), dumpfile, append=TRUE)
							write(this[,,i], dumpfile, append=append, ncolumns=dim(this)[1])
						}
					}
					else if(length(dim_all(this))==4){
						olp_beam <- 13
						olp_fans <- c(1,12,13,14,15,16,17)
						for(i in olp_fans){
							write(paste("Displaying 25 voxels of beam", (i-1)*dim(this)[3]+olp_beam), dumpfile, append=TRUE)
							write(this[,1:25,olp_beam,i], dumpfile, append=append, ncolumns=dim(this)[1])
						}
					}
				}
				else if(noisenames[j]=="olps"){
					if(length(dim_all(this)) %in% c(1,2)){
						write(this, dumpfile, append=append, ncolumns=10)
					}
					else if(length(dim_all(this))==3){
						olp_fans <- c(1,12,13,14,15,16,17)
						for(i in olp_fans){
							write(paste("Displaying fan",i), dumpfile, append=TRUE)
							write(this[,,i], dumpfile, append=append, ncolumns=dim(this)[1])
						}
					}
					else if(length(dim_all(this))==4){
						olp_beam <- 13
						olp_fans <- c(1,12,13,14,15,16,17)
						for(i in olp_fans){
							write(paste("Displaying 25 voxels of beam", (i-1)*dim(this)[3]+olp_beam), dumpfile, append=TRUE)
							write(this[,1:25,olp_beam,i], dumpfile, append=append, ncolumns=dim(this)[1])
						}
					}
				}
				else if(noisenames[j]=="seed"){
					write(prettyIntegers(this), dumpfile, append=append, ncolumns=10)
				}
				else{
					write(this, dumpfile, append=append, ncolumns=10)
				}
			}
		}
	}
		
	# Random selection and scaling of the targets:
	else if(identical(tolower(type[1]), "rand.sel_scale")){
		write("\n\n\n##### RANDOM SELECTION AND SCALING OF THE TARGETS: #####", dumpfile, append=TRUE)
		# If the length of the variable is non-zero, write the value:
		for(j in seq_along(rand.sel_scalenames)){
			this <- data[[rand.sel_scalenames[j]]]
			if(length(this)>0){
				if(is.list(this)){
					warning(paste0("The variable ", rand.sel_scalenames[j], " is specified inconsistantly with the time steps of the simulation event"))
				}
				write(getDescription(rand.sel_scalenames[j]), dumpfile, append=TRUE)
				write(this, dumpfile, append=append)
			}
		}
	}
}
