#*********************************************
#*********************************************
#' Converts from radial overlap values 'w' in sonar and echosounder systems to factors 'scale' used for calibrating simulations made by echoIBM() using the method of superimposing multiple sine waves for generation of correlated exponentially distributed values.
#'
#' @param w  is the radial overlap value, specifically the duration of the sine waves used to generate the correlated exponentially distributed values, usually calculated as the ratio of pulselength and sampling interval duration multiplied by the factor 3/4 (experienced from previous simulations to be suitable for the MS70 sonar and the SX90 sonar).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.w2scale
#'
echoIBM.w2scale<-function(w){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-09-18 - Clean version.
	########### DESCRIPTION: ###########
	# Converts from radial overlap values 'w' in sonar and echosounder systems to factors 'scale' used for calibrating simulations made by echoIBM() using the method of superimposing multiple sine waves for generation of correlated exponentially distributed values.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---w--- is the radial overlap value, specifically the duration of the sine waves used to generate the correlated exponentially distributed values, usually calculated as the ratio of pulselength and sampling interval duration multiplied by the factor 3/4 (experienced from previous simulations to be suitable for the MS70 sonar and the SX90 sonar).
	
	
	##################################################
	##################################################
	# Assume intersection at 0 length of the sine waves:
	a=0
	b=2.455 # Found using the script "Test_of_rexp_MultSines.R" in the "extdata" directory of the echoIBM package.
	a + b*w
	##################################################
	##################################################
	}
