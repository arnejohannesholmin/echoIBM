#*********************************************
#*********************************************
#' Calculates and (optionally) writes to file the spherical surface integral of a beam pattern from acoustic fish models at radius 1 meter, for the grid of the product 'kL' of wavenumber and size of the acoustic model.
#'
#' @param x  is either a list of the data ('grsf' (grid of the product of wavenumber and size of the acoustic model) and 'ebpf' (empirical beam pattern table) and at least one of 'graf' (grid azimuth angle) and 'gref' (grid elevation angle)) or a path to the TSD file holding the empirical beam pattern.
#' @param outfile  is the path to the file to which the data should be written, FALSE if no data should be written to file, or if length(outfile) == 0 the default file name which is the input file added "_ssif" to the name is used (only applied in the case that 'x' is a TSD-file).
#' @param pres  is the desired presition of the integration (used in integrateonsphere()).
#' @param max.cells  is the maximum number of cells in the grid (used in integrateonsphere()).
#' @param onlyssif.out  is TRUE if the output should be only 'ssif', not the entire list read, added 'ssif'.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSD write.TSD zeros
#'
#' @export
#' @rdname write.ssif
#'
write.ssif<-function(x=NULL,outfile=NULL,pres=1e-6,max.cells=1e6,onlyssif.out=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-08-19 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates and (optionally) writes to file the spherical surface integral of a beam pattern from acoustic fish models at radius 1 meter, for the grid of the product 'kL' of wavenumber and size of the acoustic model.
	########## DEPENDENCIES: ###########
	# read.TSDs(), beamPattern.TSD(), integrateonsphere(), write.TSD()
	############ VARIABLES: ############
	# ---x--- is either a list of the data ('grsf' (grid of the product of wavenumber and size of the acoustic model) and 'ebpf' (empirical beam pattern table) and at least one of 'graf' (grid azimuth angle) and 'gref' (grid elevation angle)) or a path to the TSD file holding the empirical beam pattern.
	# ---outfile--- is the path to the file to which the data should be written, FALSE if no data should be written to file, or if length(outfile) == 0 the default file name which is the input file added "_ssif" to the name is used (only applied in the case that 'x' is a TSD-file).
	# ---pres--- is the desired presition of the integration (used in integrateonsphere()).
	# ---max.cells--- is the maximum number of cells in the grid (used in integrateonsphere()).
	# ---onlyssif.out--- is TRUE if the output should be only 'ssif', not the entire list read, added 'ssif'.
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(is.list(x)){
		if(!any(is.null(x$grsf), is.null(x$ebpf)) && !all(is.null(x$graf), is.null(x$gref))){
			warning("the input must contain the elements \"grsf\" (grid of the product of wavenumber and size of the acoustic model) and \"ebpf\" (empirical beam pattern table) and at least one of \"graf\" (grid azimuth angle) and \"gref\" (grid elevation angle)")
			}
		if(length(outfile)==0){
			# Do not write to file if the input is a list and 'outfile' is not specified:
			outfile <- FALSE
			}
		}
	else if(file.exists(x)){
		# Define the output file if missing:
		if(length(outfile)==0){
			outfile <- paste0(substr(x, 1, nchar(x)-4), "_ssif.tsd")
			}
		# Read the relevant x:
		x <- read.TSD(x, header=FALSE, var="all")
		}
	else{
		stop("Invalid input 'x'. Must be a list of the required elements or the path to a TSD-file holding the required variables (\"grsf\" (grid of the product of wavenumber and size of the acoustic model) and \"ebpf\" (empirical beam pattern table) and at least one of \"graf\" (grid azimuth angle) and \"gref\" (grid elevation angle))")
		}
	# Prepare the beam pattern function of the fish:
	bptf <- beamPattern.TSD(x)
	# Put dimension on the empirical beam pattern:
	if(!any(is.null(x$ebp), is.null(x$dbp))){
		dim(x[[which(substr(names(x),1,3)=="ebp")]]) <- x$dbp
		}
	
		
	
	##### Execution #####
	# If empirical beam pattern is present, it may be a function of one or two angles:
	if(!length(x$ebpf)==0){
		x$ssif <- zeros(length(x$grsf))
		cat("kL-range: (",paste0(round(range(x$grsf), digits=3), collapse=", "), ")\n")
			
		for(i in seq_along(x$ssif)){
			cat(x$grsf[i]," ",sep="")
			if(!length(x$graf)==0 && !length(x$gref)==0){
				x$ssif[i] <- integrateonsphere(function(y) bptf(list(dira=y[,1], dire=y[,2], wnsz=x$grsf[i])), ndim=2, pres=pres, max.cells=max.cells, print=FALSE)$out
				}
			else if(!length(x$gref)==0){
				x$ssif[i] <- integrateonsphere(function(y) bptf(list(dire=y, wnsz=x$grsf[i])), ndim=1, pres=pres, max.cells=max.cells, print=FALSE)$out
				}
			else{
				stop("Elevation grid angles \"gref\" is missing")
				}	
			}
		}
	# Else use the parametric beam pattern function:
	else if(!length(x$pbpf)==0){
		integrateonsphere(function(x) bptf(list(dire=x, wnsz=x$wnsz)), ndim=1, pres=pres, max.cells=max.cells, print=FALSE)$out
		}
	else{
		stop("Beam pattern missing")
		}
	
	
	##### Output #####
	# Write the data:
	if(!identical(outfile,FALSE)){
		write.TSD(x, outfile, dimension=TRUE)
		}
	# Return the data:
	if(onlyssif.out){
		x$ssif
		}
	else{
		x
		}
	##################################################
	##################################################
	}
