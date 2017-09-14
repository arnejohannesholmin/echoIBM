#*********************************************
#*********************************************
#' Returns a function for the beam pattern of an acoustic source. The function beamPattern() has 3 methods depending on the type of the input object 'data'. Inputs are prioritized as (1) functions, (2) NULL, (3) empirical tables, (4) character strings naming functions.
#'
#' @param data  is a list elements defining the beam pattern (parametric or empirical). Names for the elements of the list adopted from read.TSD(). Asterix "*" is either "f", "1" or "2", representing the school (fish), the echo sounder at emission, or the echo sounder at reception:
#' @param method  is "closest" if the beam pattern value of the closest grid point is to be selected, and "linear" if linear interpolation should be used to extract the beam pattern value (time demanding and only available for 2D grids).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ones zeros
#' @importFrom fBasics linearInterpp
#'
#' @export
#' @rdname beamPattern.TSD
#'
beamPattern.TSD<-function(data,method=c("closest","linear")){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-02-11 - Clean version.
	# Update: 2011-02-20 - Added transformation from degrees to radians, and added normalizing empirical beam patter for each value of the product of wave number and size.
	# Last: 2012-02-07 - Added circularpiston_ellipticradius_sidelobefit.TSD.
	########### DESCRIPTION: ###########
	# Returns a function for the beam pattern of an acoustic source. The function beamPattern() has 3 methods depending on the type of the input object 'data'. Inputs are prioritized as (1) functions, (2) NULL, (3) empirical tables, (4) character strings naming functions.
	########## DEPENDENCIES: ###########
	# circularPiston.TSD(), lineSource.TSD(), pointSource.TSD()
	############ VARIABLES: ############
	# ---data--- is a list elements defining the beam pattern (parametric or empirical). Names for the elements of the list adopted from read.TSD(). Asterix "*" is either "f", "1" or "2", representing the school (fish), the echo sounder at emission, or the echo sounder at reception:
	#
	#		[['pbp*']] - The parametric beam pattern of the source as a function of one or two angle variables and a variable representing the relative size of the source. Given either as a function or as the name of a predefined function (one of "pointSource", "lineSource" or "circularPiston").
	#		[['gra*']] - A vector of arbitrary length representing the grid vector of azimuth angle to the source. Only required if empirical beam pattern 'ebp*' as a function of two angles of direction is given.
	#		[['gre*']] - A vector of arbitrary length representing the grid vector of elevation (incidence) angle to the source. Must be given along with the empirical beam pattern 'ebp*'.
	#		[['grs*']] - A vector of arbitrary length representing the grid vector of size of the source relative to the wavelength, represented by the product of size of the source and wave number.
	#		[['gri*']] - A vector of arbitrary length representing the grid numbering index of the sources.
	#		[['gro*']] - A vector of arbitrary length representing the grid oblongness of the sources.
	#		[['dbp*']] - A vector representing the dimension of the empirical beam pattern.
	#		[['ebp*']] - An array of dimension no less than c(length(data$graf), length(data$gref), length(data$rszf)) representing the empirical beam pattern values of the fish.
	# ---method--- is "closest" if the beam pattern value of the closest grid point is to be selected, and "linear" if linear interpolation should be used to extract the beam pattern value (time demanding and only available for 2D grids).
		

	##################################################
	##################################################
	##### Preparation, execution #####
	# The possible grid vectors "gra" (grid azimuth angle), "gre" (grid elevation angle), "grs" (grid relative size), "gri" (grid index), "gro" (grid oblongness):
	legalgrid=c("gra","gre","grs","gri","gro")
	legalinput=c("dira","dire","wnsz","indi","obln")
					
	nameroot=substr(names(data),1,3)
	object=substr(names(data),4,4)
	if(length(unique(object[nameroot %in% c(legalgrid,"ebp")]))>1){
		stop("Data not from the same category (\"school\", \"1\" (emission), \"2\" (reception))")
		}
	names(data)<-nameroot
	
	# WARNING: Empirical beam patterns overrides parametric beam patterns:
	if(length(data$pbp)>0 && length(data$ebp)==0){
		if(is.function(data$pbp)){
			fun=data$pbp
			}
		else if(tolower(data$pbp) %in% c("p", "ps", "pointsource", "pointsource.tsd")){
			fun=pointSource.TSD
			}
		else if(tolower(data$pbp) %in% c("l", "ls", "linesource", "linesource.tsd")){
			fun=lineSource.TSD
			}
		else if(tolower(data$pbp) %in% c("c", "cp", "circularpiston", "circularpiston.tsd")){
			fun=circularPiston.TSD
			}
		else if(tolower(data$pbp) %in% c("ce", "cpe", "circularpiston_ellipticradius", "circularpiston_ellipticradius.tsd")){
			fun=circularPiston_ellipticRadius.TSD
			}
		else if(tolower(data$pbp) %in% c("ces", "cpes", "circularpiston_ellipticradius_sidelobefit", "circularpiston_ellipticradius_sidelobefit.tsd")){
			fun=circularPiston_ellipticRadius_sidelobefit.TSD
			}
		else{
			warning("Parametric beam pattern not recognized (must be one of \"pointSource\", \"lineSource\", \"circularPiston\", \"circularPiston_ellipticRadius\", \"circularPiston_ellipticRadius_sidelobefit.TSD\")")
			}
		}
	else if(length(data$ebp)>0){
		# If the given parametric beam patter in "pointsource", any empirical beam pattern is ignored:
		if(length(data$pbp)>0){
			if(tolower(data$pbp)=="pointsource"){
				warning("Both parametric and empirical beam pattern present in 'data'. Since pbp=\"pointsource\", the empirical beam pattern is ignored")
				return(pointSource.TSD)
				}
			else{
				warning("Both parametric and empirical beam pattern present in 'data'. Empirical chosen")
				}
			}
		# The grid vectors that are present in the data:
		gridpresentind=which(legalgrid %in% names(data))
		legalgrid=legalgrid[gridpresentind]
		legalinput=legalinput[gridpresentind]
		
		# Grid angle vectors are required to be in radians, and if values>2*pi are found, degrees are assumed in the input and the grid angle vectors are transformed to radians:
		if(length(data$gra)>0 && max(data$gra)>2*pi){
			warning("Grid azimuth angles ($gra) for the empirical beam pattern transformed from degrees to radians")
			data$gra=data$gra*pi/180 
			}
		if(length(data$gre)>0 && max(data$gre)>2*pi){
			warning("Grid elevation angles ($gre) for the empirical beam pattern transformed from degrees to radians")
			data$gre=data$gre*pi/180 
			}
			
		# The dimension of the empirical beam pattern:
		lengths=sapply(data[legalgrid],length)
		if(length(data$dbp)==0){
			dim(data$ebp)=lengths
			lendiff=zeros(length(lengths))
			}
		else{
			if(!identical(dim(data$ebp),data$dbp)){
				dim(data$ebp)=data$dbp
				}
			lendiff=lengths-data$dbp
			}
		
		# Expand the grid elevation vector if given only in the range [0,90] degrees:
		if(all(range(data$gre)==c(0,pi/2)) && "gre" %in% legalgrid && length(dim(data$ebp))==2){
			data$gre=c(data$gre,(pi-rev(data$gre))[-1])
			data$ebp=rbind(data$ebp,data$ebp[rev(seq_len(dim(data$ebp)[1]))[-1],])
			}
		
		# Normalize the empirical beam patterns to max=1:
		# which(gridpresentind==3) gets the position of the relative size grid vector in the empirical beam pattern:
		atgrs=which(gridpresentind==3)
		ndim=length(dim(data$ebp))
		#####################################################################################################################
		##### Get the maximum value of the beam pattern for each relative size value (product of wave number and size): #####
		#####################################################################################################################
		normalization=rep(apply(data$ebp,atgrs,max),prod(dim(data$ebp)[-atgrs]))
		# Organize 'normalization' to match the dimension and order of 'data$ebp':
		dim(normalization)=c(length(data$grs),dim(data$ebp)[-atgrs])
		newpermutation=ones(ndim)
		newpermutation[-atgrs]=2:ndim
		normalization=aperm(normalization,newpermutation)
		# Apply the normalization:
		data$ebp=data$ebp/normalization
		#####################################################################################################################
				
		ndbp=length(data$dbp)
		# The lengths of the explanatory variables should be equal to or one larger than the dimensions of the response variable:
		if(!all(lendiff %in% 0:1)){
			warning("Lengths of explanatory variables should be equal to or one larger than the dimensions of the response variable")
			}
		
		# Execution:
		if(any(sapply(data[legalgrid],is.unsorted))){
			stop("The grid vectors need to be sorted increasingly")
			}
		
		if(method[1]=="closest" || ndbp>3){
			fun=function(l){
				inputpresent=match(names(l),legalinput)
				inputpresent=inputpresent[!is.na(inputpresent)]
				legalgrid=legalgrid[inputpresent]
				inputpresent=legalinput[inputpresent]
				l=l[inputpresent]
				ll=zeros(length(l[[1]]),length(l))
				# Finding the indexes of the input arguments in the grid given by 'data':
				for(i in seq_along(l)){
					# If the length of the grid variable equals the corrsponding dimension of the empirical beampattern, the closest grid point is selected, rounding 0.5 up to 1:
					if(lendiff[i]==0){
						thisgrid=data[[legalgrid[i]]]
						ll[,i]=findInterval(l[[i]],c(-Inf,thisgrid[-length(thisgrid)]+diff(thisgrid)/2,Inf),rightmost.closed=TRUE,all.inside=TRUE)
						}
					else{
						ll[,i]=findInterval(l[[i]],thisgrid,all.inside=TRUE,rightmost.closed=TRUE)
						}
					}
				# Extracting the beam pattern values:
				data$ebp[ll]
				}
			}
		else if(method[1]=="linear"){
			fun=function(l){
				inputpresent=match(names(l),legalinput)
				inputpresent=inputpresent[!is.na(inputpresent)]
				#legalgrid=which(legalgrid[inputpresent] %in% names(data))
				legalgrid=legalgrid[inputpresent]
				inputpresent=legalinput[inputpresent]
				l=l[inputpresent]
				thisgrid=expand.grid(data[legalgrid])
				thislgrid=expand.grid(l[1:2])
				# Finding the indexes of the input arguments in the grid given by 'data':
				fBasics::linearInterpp(thisgrid[,1], thisgrid[,2], data$ebp, thislgrid[,1], thislgrid[,2])
				}
			}
		else{
			stop("Wrong 'method'")
			}
		}
	else{
		if(object[1]=="f"){
			warning("Parametric and empirical beam pattern function missing and defaulted to \"lineSource\"")
			fun=lineSource.TSD
			}
		else if(object[1] %in% c("1","2")){
			warning("Parametric and empirical beam pattern function missing and defaulted to \"circularPiston\"")
			fun=circularPiston.TSD
			}
		}

		
	##### Output #####
	fun
	##################################################
	##################################################
	}
