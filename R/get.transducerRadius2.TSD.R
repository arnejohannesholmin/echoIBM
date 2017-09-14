#*********************************************
#*********************************************
#' Calculates the radii of circular transducer elements based on -3 dB beam angles.
#'
#' @param beamwidth  is a vector of -3 dB beamwidths for the beams given in radians (the angle between the directions in which the acoustic intensity is half of the maximum acoustic intensity).
#' @param beampattern  is the beam pattern function having angle in radians as its first argument and the product of wave number and transducer element radius as its second argument, e. g. circularPiston(ang,kb).
#' @param wavenumber  is a vector of wave numbers for the beams.
#' @param max.radius  is an expected maximum radius which should be low enough to avoid exhaustive or inaccurate calculations.
#' @param twoway  is TRUE if the beampattern is to be squared, corresponding to a two way beam forming. I.e., if the beam widths are given as two-way angles bwtl=4 degrees, the one way beam width is actually 4*sqrt(2)=5.656854, and twoway should be TRUE to account for this.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD zeros
#' @importFrom stats uniroot
#'
#' @export
#' @rdname get.transducerRadius2.TSD
#'
get.transducerRadius2.TSD<-function(data){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: #############
	# Start: 2010-01-20 - Clean version.
	# Update: 2010-01-25 - Compatibility changed to only support beam patterns of exactly 2 arguments: angle of incidence and the product of wave number and size of the source.
	# Last: 2010-12-07 - Added the option 'twoway' for two way beam forming (beam pattern squared).
	########### DESCRIPTION: ###########
	# Calculates the radii of circular transducer elements based on -3 dB beam angles.
	########## DEPENDENCIES: ###########
	# beampattern(), zeros()
	############ VARIABLES: ############
	# ---beamwidth--- is a vector of -3 dB beamwidths for the beams given in radians (the angle between the directions in which the acoustic intensity is half of the maximum acoustic intensity).
	# ---beampattern--- is the beam pattern function having angle in radians as its first argument and the product of wave number and transducer element radius as its second argument, e. g. circularPiston(ang,kb).
	# ---wavenumber--- is a vector of wave numbers for the beams.
	# ---max.radius--- is an expected maximum radius which should be low enough to avoid exhaustive or inaccurate calculations.
	# ---twoway--- is TRUE if the beampattern is to be squared, corresponding to a two way beam forming. I.e., if the beam widths are given as two-way angles bwtl=4 degrees, the one way beam width is actually 4*sqrt(2)=5.656854, and twoway should be TRUE to account for this.
	
	
	##################################################
	##################################################
	##### Preparation and execution #####
	# Defaults:
	if(length(data$max.radius)==0){
		data$max.radius=0.2
		}
	if(length(data$twoway)==0){
		data$twoway=FALSE
		}
	
	# Division by 2 because the beamwidth is defined as twice the angle where the intensity is 1/2 of that on the axis of acoustic maximum:
	data$bwty=data$bwty/2
	data$bwtx=data$bwtx/2
	# If the beamwidths are given, the radii 'radius' of the transducer elements are calculated by solving beampattern=1/2:
	Ni=max(length(data$bwty),length(data$bwtx),length(data$wavenumber))
	data$bwty=rep(data$bwty,length.out=Ni)
	data$bwtx=rep(data$bwtx,length.out=Ni)
	data$wavenumber=rep(data$wavenumber,length.out=Ni)
	data$wnsz=data$wavenumber*data$max.radius
	
	# If max.radius is chosen too low, the minimum of 'wavenumber' might cause the beam pattern to stay above 0.5:
	data$dira=0
	# This is the elevation angle in the coordinate system of the transducer, which has z axis along the beam maximum, so elevation angle is the incidence angle:
	data$dire=data$bwtx
	if(any(data$beampattern(data)>=0.5)){
		stop(paste("max.radius=",data$max.radius,"to low"))
		}
	# Similarly for dira:
	data$dira=pi/2
	# This is the elevation angle in the coordinate system of the transducer, which has z axis along the beam maximum, so elevation angle is the incidence angle:
	data$dire=data$bwty
	if(any(data$beampattern(data)>=0.5)){
		stop(paste("max.radius=",data$max.radius,"to low"))
		}
		
	# 'size' is the relative size of the echo sounder elements, to be derived from the beam widths:
	size=zeros(Ni,2)
	if(data$twoway){
		for(i in 1:Ni){
			fx=function(s){
				thisdata=list()
				thisdata$dira=0
				# This is the elevation angle in the coordinate system of the transducer, which has z axis along the beam maximum, so elevation angle is the incidence angle:
				thisdata$dire=data$bwtx[i]
				thisdata$wnsz=s
				thisdata$sllf=data$sllf[i,1]
				thisdata$indi=data$indi[i]
				1/2-data$beampattern(thisdata)^2
				}
			fy=function(s){
				thisdata=list()
				thisdata$dira=pi/2
				# This is the elevation angle in the coordinate system of the transducer, which has z axis along the beam maximum, so elevation angle is the incidence angle:
				thisdata$dire=data$bwty[i]
				thisdata$wnsz=s
				thisdata$sllf=data$sllf[i,2]
				thisdata$indi=data$indi[i]
				1/2-data$beampattern(thisdata)^2
				}
			# Storing the relative size of the echo sounder element:
			size[i,1]=uniroot(fx,c(.Machine$double.eps,data$wavenumber[i]*data$max.radius))$root
			size[i,2]=uniroot(fy,c(.Machine$double.eps,data$wavenumber[i]*data$max.radius))$root
			}
		}
	else{
		for(i in 1:Ni){
			fx=function(s){
				thisdata=list()
				thisdata$dira=0
				# This is the elevation angle in the coordinate system of the transducer, which has z axis along the beam maximum, so elevation angle is the incidence angle:
				thisdata$dire=data$bwtx[i]
				thisdata$wnsz=s
				thisdata$sllf=data$sllf[i,1]
				thisdata$indi=data$indi[i]
				1/2-data$beampattern(thisdata)
				}
			fy=function(s){
				thisdata=list()
				thisdata$dira=pi/2
				# This is the elevation angle in the coordinate system of the transducer, which has z axis along the beam maximum, so elevation angle is the incidence angle:
				thisdata$dire=data$bwty[i]
				thisdata$wnsz=s
				thisdata$sllf=data$sllf[i,2]
				thisdata$indi=data$indi[i]
				1/2-data$beampattern(thisdata)
				}
			# Storing the relative size of the echo sounder element:
			size[i,1]=uniroot(fx,c(.Machine$double.eps,data$wavenumber[i]*data$max.radius))$root
			size[i,2]=uniroot(fy,c(.Machine$double.eps,data$wavenumber[i]*data$max.radius))$root
			}
		}
	
	
	##### Output #####
	size/data$wavenumber
	##################################################
	##################################################
	}
