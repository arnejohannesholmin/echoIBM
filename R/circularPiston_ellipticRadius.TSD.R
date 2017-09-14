#*********************************************
#*********************************************
#' Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius as an elliptical function of azimuth angle 'data$dira' given the radii along the x-axis, data$wnsz[,1], and along the y-axis, data$wnsz[,2], as a function of incidence angle 'data$dire'.
#'
#' @param data  is a list containing the elements:
#' @param appr  is an approximating function used to replace the expression (2*besselJ(ang,1)/ang)^2.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname circularPiston_ellipticRadius.TSD
#'
circularPiston_ellipticRadius.TSD<-function(data,appr=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-01-25 - Clean version.
	# Last: 2012-01-07 - Simplified the preparations and reduced process time.
	########### DESCRIPTION: ###########
	# Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius as an elliptical function of azimuth angle 'data$dira' given the radii along the x-axis, data$wnsz[,1], and along the y-axis, data$wnsz[,2], as a function of incidence angle 'data$dire'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the elements:
	#		'dira' (azimuth angle in radians [0,pi/2]), given as a vector
	#		'dire' (elevarion angle in radians [0,pi/2]), given as a vector
	#		'wnsz' (product of wave number and size), given as a two column matrix
	# ---appr--- is an approximating function used to replace the expression (2*besselJ(ang,1)/ang)^2.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Define the output for negative angles:
	if(min(data$dire)<0){
		data$dire=abs(data$dire)
		}
	# Clamp 'ang' to the absolute range [0,pi/2]:
	if(max(data$dire)>pi/2){
		data$dire[data$dire>pi/2]=pi/2
		}
	
	# The input 'wnsz' needs to be either a non-dimensional or one-dimensional vector or a two column matrix, in which case the transducer element radii and the side lobe specifications differ in the x- and y-direction (elliptical beam patterns):
	if(length(dim(data$wnsz))==0 || identical(ncol(data$wnsz),1L)){
		data$wnsz=cbind(data$wnsz,data$wnsz)
		}
	else if(ncol(data$wnsz)!=2L){
		data$wnsz=c(data$wnsz)
		data$wnsz=cbind(data$wnsz,data$wnsz)
		}
	
	# Define the radius of the ellipse representing the transducer element. An ad hoc method is used where the beam pattern is found by the circular piston of radius as an elliptical function of azimuth angle in the coordinate system of the transducer:
	# The following calculations are based on Rottmann, 7'th norwegian edition, page 50:
	# 'eccentricity2' is the squared eccentricity, as given in Rottmann: 
	eccentricity2=(data$wnsz[,1]^2-data$wnsz[,2]^2)/data$wnsz[,1]^2
	rho=sqrt(data$wnsz[,2]^2/(1-eccentricity2*cos(data$dira)^2))
	# Set the agrument to the Bessel function:
	ang=sin(data$dire) * rho
	
	
	##### Execution and output #####
	if(is.function(appr)){
		out=appr(ang)
		}
	else{
		out=(2*besselJ(ang,1)/ang)^2
		}
	out[ang==0]=1
	out
	##################################################
	##################################################
	}
