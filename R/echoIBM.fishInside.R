#*********************************************
#*********************************************
#' Discards fish outside of the sampling region of the sonar/echosounder.
#'
#' @param dynschool  is a list of the dynamic fish information to be subsetted.
#' @param data  is a list holding the variables 'esnm', and vessel spesifications.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom sonR is.sonar rotate3D
#' @importFrom TSD car2sph
#'
#' @export
#' @rdname echoIBM.fishInside
#'
echoIBM.fishInside<-function(dynschool, data, nw=c(r=Inf,az=Inf,el=Inf), dumpfile){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-02-27 - Clean version.
	########### DESCRIPTION: ###########
	# Discards fish outside of the sampling region of the sonar/echosounder.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---dynschool--- is a list of the dynamic fish information to be subsetted.
	# ---data--- is a list holding the variables 'esnm', and vessel spesifications.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(all(is.infinite(nw[2:3]))){
		return(dynschool)
		}
	# Rotate the fish positions into a special coordinate system, which puts echosounders as well as sonars straight out from the port side of the vessel:
	if(is.sonar(bydirs=TRUE,data=data,margin=10*pi/180)){
		beamsRotated=rotate3D(cbind(1,data$dira,data$dire),by="z",ang=pi/2,sph.in=TRUE,sph.out=TRUE)
		data$dira=beamsRotated[,2]
		data$dire=beamsRotated[,3]
		}
	fishposInSonar=car2sph(dynschool$psxf,dynschool$psyf,dynschool$pszf)
	
	maxr=(max(data$lenb)+nw[1]) * soundbeam_range(data, pos="rres")
	maxba=max(data$bwtl)
	maxbe=max(data$bwtt)
	mina=min(data$dira-maxba*nw[2])
	maxa=max(data$dira+maxba*nw[2])
	mine=min(data$dire-maxbe*nw[3])
	maxe=max(data$dire+maxbe*nw[3])
	
	
	########## Execution and output ##########
	# Discard fish outside of the volume at output:
	inside = fishposInSonar[,1]<maxr  &  mina<fishposInSonar[,2] & fishposInSonar[,2]<maxa  &  mine<fishposInSonar[,3] & fishposInSonar[,3]<maxe
	if(length(dumpfile)>0 && nchar(dumpfile)>0){
		write(paste("\n\n# Proportion of fish discarded outside of the observation volume: ",mean(!inside),sep=""),dumpfile,append=TRUE)
		}
	lapply(dynschool,function(x) if(length(x)==length(inside)) x[inside] else x)
	##################################################
	##################################################
	}
