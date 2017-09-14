#*********************************************
#*********************************************
#' Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius as an elliptical function of azimuth angle 'ang1' given the radii along the x-axis, kb[,1], and along the y-axis, kb[,2], as a function of incidence angle 'ang2'.
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
#' @rdname circularPiston_ellipticRadius_sidelobefit.TSD
#'
circularPiston_ellipticRadius_sidelobefit.TSD<-function(data,appr=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-01-25 - Clean version.
	# Last: 2012-01-07 - Simplified the preparations and reduced process time.
	########### DESCRIPTION: ###########
	# Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius as an elliptical function of azimuth angle 'ang1' given the radii along the x-axis, kb[,1], and along the y-axis, kb[,2], as a function of incidence angle 'ang2'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the elements:
	#		'dira' (azimuth angle in radians [0,pi/2]), given as a vector
	#		'dire' (elevarion angle in radians [0,pi/2]), given as a vector
	#		'wnsz' (product of wave number and size), given as a two column matrix
	#		'sllf' is a two column matrix of factors adjusting the beam pattern to match the side lobe level to the desired value, given by the horizontal factor in the first column and the vertical in the second column. These factor must be given for each device. And as an example, the MS70 has side lobes at -25dB horizontally and -35 dB vertically. The simulation model based on the circular piston has -17 dB norizontally and -34 dB vertically. Thus the sllf is c(25/17,35/34).
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
	
	# The two inputs 'wnsz' and sllf' need to be either non-dimensional or one-dimensional vectors or two column matrices, in which case the transducer element radii and the side lobe specifications differ in the x- and y-direction (elliptical beam patterns):
	if(length(dim(data$wnsz))==0 || identical(ncol(data$wnsz),1L)){
		data$wnsz=cbind(data$wnsz,data$wnsz)
		}
	else if(ncol(data$wnsz)!=2L){
		data$wnsz=c(data$wnsz)
		data$wnsz=cbind(data$wnsz,data$wnsz)
		}
	if(length(dim(data$sllf))==0 || identical(ncol(data$sllf),1L)){
		data$sllf=cbind(data$sllf,data$sllf)
		}
	else if(ncol(data$sllf)!=2L){
		data$sllf=c(data$sllf)
		data$sllf=cbind(data$sllf,data$sllf)
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
	
	# Apply the side lobe suppression by raising the beam pattern to a power given by sllf (side lobe level factor):
	out[ang==0]=1
	# Apply the sidelobe specification:
	eccentricity2sl=(data$sllf[,1]^2-data$sllf[,2]^2)/data$sllf[,1]^2
	# 'rho' is the distance to the perifery of the ellipse representing the side lobe level: 
	exponent=sqrt(data$sllf[,2]^2/(1-eccentricity2sl*cos(data$dira)^2))
	
	# Output:
	out^exponent
	##################################################
	##################################################
	}
