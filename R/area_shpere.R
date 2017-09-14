#*********************************************
#*********************************************
#' Calculates the surface area of segments of spheres of radius 'r', defined by the azimuth and elevation angles 'theta' and 'phi'. The calculation os based on the expression S = 2 * pi * R * h, given on http://mathworld.wolfram.com/SphericalCap.html. 'h' can be expressed as R * (1-cos(phi)), and for two elevation angles phi2 >= phi1 we get the following surface area of the resulting "disc" (can be illustrated by cutting an orange by two parallel intersections, and calculate the are of the orange skin): S = 2 * pi * R^2 * ( (1-cos(phi2)) - (1-cos(phi1))) = 2 * pi * R^2 * ( cos(phi1) - cos(phi2)). Splitting the resulting surface into parts by azimuth angle is straight foreward: S = R^2 * ( (1-cos(phi2)) - (1-cos(phi1))) * (theta2 - theta1).
#'
#' @param r  is the radii of the spheres.
#' @param theta  is either a vector of two elements representing the lower and upper azimuth angle of the spherical sector, or a two column matrix of the same.
#' @param phi  is either a vector of two elements representing the lower and upper elevation angle of the spherical sector, or a two column matrix of the same.
#' @param deg  is TRUE if the angles are given in degrees.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname area_shpere
#'
area_shpere<-function(r=1,theta=c(0,2*pi),phi=c(0,pi),deg=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-06-19 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the surface area of segments of spheres of radius 'r', defined by the azimuth and elevation angles 'theta' and 'phi'. The calculation os based on the expression S = 2 * pi * R * h, given on http://mathworld.wolfram.com/SphericalCap.html. 'h' can be expressed as R * (1-cos(phi)), and for two elevation angles phi2 >= phi1 we get the following surface area of the resulting "disc" (can be illustrated by cutting an orange by two parallel intersections, and calculate the are of the orange skin): S = 2 * pi * R^2 * ( (1-cos(phi2)) - (1-cos(phi1))) = 2 * pi * R^2 * ( cos(phi1) - cos(phi2)). Splitting the resulting surface into parts by azimuth angle is straight foreward: S = R^2 * ( (1-cos(phi2)) - (1-cos(phi1))) * (theta2 - theta1).
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---r--- is the radii of the spheres.
	# ---theta--- is either a vector of two elements representing the lower and upper azimuth angle of the spherical sector, or a two column matrix of the same.
	# ---phi--- is either a vector of two elements representing the lower and upper elevation angle of the spherical sector, or a two column matrix of the same.
	# ---deg--- is TRUE if the angles are given in degrees.
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(deg){
		theta=theta*pi/180
		phi=phi*pi/180
		}
	if(length(dim(theta))==2){
		theta=theta[,2]-theta[,1]
		}
	else{
		theta=theta[2]-theta[1]
		}
	if(length(dim(phi))==2){
		phi=cos(phi[,1])-cos(phi[,2])
		}
	else{
		phi=cos(phi[1])-cos(phi[2])
		}
	
	
	##### Execution and output #####
	r^2*theta*phi
	##################################################
	##################################################
	}
