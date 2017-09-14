#*********************************************
#*********************************************
#' The hydrostatic compression factor of a swim bladder. Two cases are supported, combining the compression of length and width, given by L(z) = L(0)(1-(rho*g*z)/p0)^gammaL and W(z) = W(0)(1-(rho*g*z)/p0)^gammaW: (1) Area: The hydrostatic compression factor of the cross sectional area of the swim bladder is found by multiplication of L and W, resulting in (1-(rho*g*z)/p0)^(gammaL+gammaW). (2) The hydrostatic compression factor of the oblongness is found by the ratio L/W, giving (1-(rho*g*z)/p0)^(gammaL-gammaW)
#'
#' @param data  is a list containing the elements "rho0", specifying the mass density of the sea, "gacc", representing the gravitational accelleration, "pszf" which is the z-position of the fish in (G), "hpr0", giving the hydrostatic pressure at sea surface (equal to air pressure at the sea surface) and "gaml" and "gamw" representing the depth compression parameters of the length and the width of the swim bladder model.
#' @param type  is a character string representing the type of the compression factor, where "area" denotes compression of the cross sectional area of the swin bladder and "obln" corresponds to changes in the oblongness of the swin bladder (see DESCRIPTION).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname etaCompression
#'
etaCompression<-function(data=list(),type=c("area","obln")){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-01-18 - Clean version.
	# Last: 2011-10-18 - Added defaults.
	########### DESCRIPTION: ###########
	# The hydrostatic compression factor of a swim bladder. Two cases are supported, combining the compression of length and width, given by L(z) = L(0)(1-(rho*g*z)/p0)^gammaL and W(z) = W(0)(1-(rho*g*z)/p0)^gammaW: (1) Area: The hydrostatic compression factor of the cross sectional area of the swim bladder is found by multiplication of L and W, resulting in (1-(rho*g*z)/p0)^(gammaL+gammaW). (2) The hydrostatic compression factor of the oblongness is found by the ratio L/W, giving (1-(rho*g*z)/p0)^(gammaL-gammaW)
	########## DEPENDENCIES: ###########
	# 	
	############ VARIABLES: ############
	# ---data--- is a list containing the elements "rho0", specifying the mass density of the sea, "gacc", representing the gravitational accelleration, "pszf" which is the z-position of the fish in (G), "hpr0", giving the hydrostatic pressure at sea surface (equal to air pressure at the sea surface) and "gaml" and "gamw" representing the depth compression parameters of the length and the width of the swim bladder model.
	# ---type--- is a character string representing the type of the compression factor, where "area" denotes compression of the cross sectional area of the swin bladder and "obln" corresponds to changes in the oblongness of the swin bladder (see DESCRIPTION).
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Defaults:
	if(length(data$rho0)==0){
		data$rho0=1026
		}
	if(length(data$gacc)==0){
		data$gacc=9.82
		}
	if(length(data$hpr0)==0){
		data$hpr0=101325
		}
	if(length(data$gaml)==0){
		data$gaml=0
		}
	if(length(data$gamw)==0){
		data$gamw=-0.23
		}
	if(length(data$pszf)==0){
		data$pszf=0
		}
	
	
	########## Execution and output ##########
	if(type[1]=="area"){
		exponent=data$gaml+data$gamw
		(1-data$rho0*data$gacc*data$pszf/data$hpr0)^(data$gaml+data$gamw)
		}
	else if(type[1]=="obln"){
		exponent=data$gaml-data$gamw
		(1-data$rho0*data$gacc*data$pszf/data$hpr0)^(data$gaml-data$gamw)
		}
	else{
		stop("Invalid type (needs to be one of \"area\" or \"obln\")")
		}
	(1-data$rho0*data$gacc*data$pszf/data$hpr0)^exponent
	##################################################
	##################################################
	}
