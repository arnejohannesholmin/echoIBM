#*********************************************
#*********************************************
#' Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius 'b', as a function of incidence data$direle 'data$dire'.
#'
#' @param data  is a list containing the following elements:
#' @param appr  is an approximating function used to replace the expression (2*besselJ(ang,1)/ang)^2.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname circularPiston.TSD
#'
circularPiston.TSD<-function(data,appr=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ Ldata$direUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-01-25 - Clean version.
	# Last: 2012-01-07 - Simplified the preparations and reduced process time.
	########### DESCRIPTION: ###########
	# Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius 'b', as a function of incidence data$direle 'data$dire'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the following elements:
	#		'dire' (elevarion angle in radians [0,pi/2])
	#		'wnsz' (product of wave number and size) 
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
	
	## The input 'wnsz' needs to be either a non-dimensional or one-dimensional vector or a two column matrix, in which case the transducer element radii and the side lobe specifications differ in the x- and y-direction (elliptical beam patterns):
	#if(length(dim(data$wnsz))==0 || identical(ncol(data$wnsz),1L)){
	#	data$wnsz=cbind(data$wnsz,data$wnsz)
	#	}
	#else if(ncol(data$wnsz)!=2L){
	#	data$wnsz=c(data$wnsz)
	#	data$wnsz=cbind(data$wnsz,data$wnsz)
	#	}
	
	data$dire=sin(data$dire) * data$wnsz
	
	
	##### Execution and output #####
	if(is.function(appr)){
		out=appr(data$dire)
		}
	else{
		out=(2*besselJ(data$dire,1)/data$dire)^2
		}
	out[data$dire==0]=1
	out
	##################################################
	##################################################
	}
