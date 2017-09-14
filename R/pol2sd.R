#*********************************************
#*********************************************
#' Gets the standard deviation to apply to the function ARsim_school().
#'
#' @param pol  is the desired polarization value, either given as the Huth measure (mean angle deviation between individual and group heading) or the Couzin measure.
#' @param speed  is the speed of the fish.
#' @param dt  is the difference between time steps.
#' @param type  see pol.school().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom akima interp
#' @importFrom TSD NAs read.TSD
#' @importFrom stats approx
#'
#' @export
#' @rdname pol2sd
#'
pol2sd<-function(pol, speed=1, dt=1, gamma=0.8, type="h"){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-03-26 - Clean version
	########### DESCRIPTION: ###########
	# Gets the standard deviation to apply to the function ARsim_school().
	########## DEPENDENCIES: ###########	
	#
	############ VARIABLES: ############
	# ---pol--- is the desired polarization value, either given as the Huth measure (mean angle deviation between individual and group heading) or the Couzin measure.
	# ---speed--- is the speed of the fish.
	# ---dt--- is the difference between time steps.
	# ---type--- see pol.school().
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Read the polarization tables:
	#echoIBM_datadir_ = file.path(echoIBM_frameworks, "R", "Functions", "echoIBM Main", "Utilities")
	#filebasename = "ARsim_school_table.TSD"
	#polfile = file.path(echoIBM_datadir_, filebasename)
	#poltables = read.TSD(polfile)
	poltables = read.TSD(system.file("extdata", "ARsim_school_table.TSD", package="echoIBM"))
	thisvar = paste("tbl", tolower(substr(type,1,1)), sep="")
	
	# If given, interpolate first to the given 'gamma':
	xy = as.matrix(expand.grid(poltables$grns, poltables$grgm))
	poltables[[thisvar]] = interp(xy[,1], xy[,2], poltables[[thisvar]], xo=poltables$grns, yo=gamma)$z
	
	
	########## Execution and output ##########
	sigma = NAs(length(pol))
	for(i in seq_along(pol)){
		sigma[i] = approx(poltables[[thisvar]], poltables$grns, pol[i])$y * speed * dt
		}
	list(sigma=sigma, gamma=gamma)
	##################################################
	##################################################
	}
