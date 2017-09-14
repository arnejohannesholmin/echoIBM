#*********************************************
#*********************************************
#' Wrapper for the method of drawing vessel positions based on school centres of mass, given as an option to the fundtion echoIBM (specifying path=TRUE). set.path.vessel.cm() uses the function set.path.vessel() to select the vessel positions. The input variables in this function are especially selected for the single purpose of being used in echoIBM(), and are thus not user friendly for any othr use.
#'
#' @param schooldirs  is a vector of the directories of the schools to be simulated.
#' @param dynschoolfiles  is a list of the dynamic school files for the schools to be simulated.
#' @param recycle  see echoIBM(). Here, if 'recycle' is different from FALSE or does not contains zeros, the first time step is selected.
#' @param t  is the vector of the time steps for which to draw the vesse movement.
#' @param tvessel  is the UNIX time values to be written to the ".vessel" file. If tvessel==NULL the school time points are used.
#' @param numt  is the number of time steps for which to select vessel positions.
#' @param smooth  is the function used for smoothing the vessel path.
#' @param margin  is a vector of the margins on either side of the span of 'x' (recycled if not of length 4).
#' @param scan.volume  is TRUE if the acoustic observation volume schould be scanned for the positions of the fish.
#' @param beams  is a list of the beam configuration.
#' @param ctd  is a list of the CTD-specification of the sea.
#' @param rph  is a matrix of two columns of length 3 representing the mean (column 1) and the standard deviation (column 2) of the roll values (rtxv), pitch values (rtyv) and heave values (przv) of the vessel, used in gaussian simulation of the roll, pitch and heave values.
#' @param event  is a vector of strings representing the paths to the directories and files containing the necessary informtation for the simulation.
#' @param eventname  is the name of the main directory of the simulation, in which to store the simulate acoustical files.
#' @param pathnr  is the number of the vessel path drawn, to be used in the name of the ".vessel" file. If a new ".vessel" file is to be written, 'pathnr' must be different than the pathnr of existing files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom sonR cm.school rotate3D
#' @importFrom TSD car2global ones read.TSD utim.TSD write.TSD zeros
#' @importFrom stats rnorm
#'
#' @export
#' @rdname set.path.vessel.cm
#'
set.path.vessel.cm<-function(schooldirs,dynschoolfiles,recycle,t,tvessel,numt,smooth,margin,scan.volume,beams,ctd,rph,event,eventname,pathnr){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-01-06 - Clean version.
	# Last: 2013-11-07 - Added support for list of function input of 'recycle'.
	########### DESCRIPTION: ###########
	# Wrapper for the method of drawing vessel positions based on school centres of mass, given as an option to the fundtion echoIBM (specifying path=TRUE). set.path.vessel.cm() uses the function set.path.vessel() to select the vessel positions. The input variables in this function are especially selected for the single purpose of being used in echoIBM(), and are thus not user friendly for any othr use.
	########## DEPENDENCIES: ###########
	# read.TSD(), utim.TSD(), zeros(), ones(), cm.school(), set.path.vessel(), car2global(), write.TSD()
	############ VARIABLES: ############
	# ---schooldirs--- is a vector of the directories of the schools to be simulated.
	# ---dynschoolfiles--- is a list of the dynamic school files for the schools to be simulated.
	# ---recycle--- see echoIBM(). Here, if 'recycle' is different from FALSE or does not contains zeros, the first time step is selected.
	# ---t--- is the vector of the time steps for which to draw the vesse movement.
	# ---tvessel--- is the UNIX time values to be written to the ".vessel" file. If tvessel==NULL the school time points are used.
	# ---numt--- is the number of time steps for which to select vessel positions.
	# ---smooth--- is the function used for smoothing the vessel path.
	# ---margin--- is a vector of the margins on either side of the span of 'x' (recycled if not of length 4).
	# ---scan.volume--- is TRUE if the acoustic observation volume schould be scanned for the positions of the fish.
	# ---beams--- is a list of the beam configuration.
	# ---ctd--- is a list of the CTD-specification of the sea.
	# ---rph--- is a matrix of two columns of length 3 representing the mean (column 1) and the standard deviation (column 2) of the roll values (rtxv), pitch values (rtyv) and heave values (przv) of the vessel, used in gaussian simulation of the roll, pitch and heave values.
	# ---event--- is a vector of strings representing the paths to the directories and files containing the necessary informtation for the simulation.
	# ---eventname--- is the name of the main directory of the simulation, in which to store the simulate acoustical files.
	# ---pathnr--- is the number of the vessel path drawn, to be used in the name of the ".vessel" file. If a new ".vessel" file is to be written, 'pathnr' must be different than the pathnr of existing files.
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(length(dynschoolfiles)>1){
		cat("Select the school to use as the reference school when selecting the vessel positions from the list below (type in the number in the list:)\n\n")
		print(schooldirs)
		answer=as.numeric(readline())
		if(!is.na(answer)){
			dynschoolfiles=dynschoolfiles[answer]
			recycle=recycle[[answer]]
			}
		else{
			stop("Wrong input (must be numberic), function terminated")
			}
		}
	else{
		recycle=recycle[[1]]
		}
	# 'dynschoolfiles' is a list and needs to be unlisted:
	dynschoolfiles=dynschoolfiles[[1]]
		
	utim=list()
	# For loop through the school files:
	for(i in seq_along(dynschoolfiles)){
		this=read.TSD(dynschoolfiles[i],var="time",t="all",header=TRUE)
		utim=c(utim,list(utim.TSD(this)))
		}
	
	# 'filet' is a matrix of two columns, where the first column will be the ping number of the file and the second will be the file number
	filet=zeros(numt,2)
	# Getting the file numbers:
	fileind=unlist(lapply(utim,function(x) x[1]))
	
	fileind=order(fileind)
	# Moving through the dynaminc school files and inserting filenumbers and ping numbers:
	
	
	# If one or more of the dynamic data of the school are to be recycled (e. g., for applying the same bottom data to all time steps), 'recycle' specifies the time steps to recycle, and using the modulus function "%%" the appropriate time steps are selected below:
	if(is.function(recycle) || all(as.numeric(recycle))){
		filet=ones(numt,2)
		}
	for(i in seq_along(dynschoolfiles)){
		ind=which(unlist(utim) %in% utim[[i]])
		filet[ind,1]=1:length(utim[[i]])
		filet[ind,2]=fileind[i]
		}
	
	
	##### Execution #####
	# Calculating center of mass of the school
	cm=NULL
	utim=NULL
	cat("Calculating center of mass of the school for time steps\n")
	for(i in t){
		cat(i,", ",sep="")
		# Read the school positions
		schoolt=read.TSD(dynschoolfiles[filet[i,2]],t=filet[i,1],var=c("psxf","psyf","pszf","time"))
		# Store the centers of mass
		cm=rbind(cm,cm.school(schoolt))
		# Stor the unix time:
		utim=c(utim,utim.TSD(schoolt))
		}
	cat("\n")
	if(is.null(tvessel)){
		tvessel=utim
		}
	vessel=set.path.vessel(cm,n=length(t),t=tvessel,smooth=smooth,margin=margin)
	
	# Report the ranges of the school and the ranges of the volume, issuing a warning if any fish are outside of the volume, defined as one beams$bwtl away from the range of theta and phi values for the centers of the beams, in the spherical coordinate system of the echo sounder:
	if(scan.volume){
		cat("Scanning sonar volume...\n")
		rangesfish=zeros(length(t),6)
		j=0
		for(i in t){
			cat(i,", ",sep="")
			j=j+1
			# Read the school positions
			schoolt=read.TSD(dynschoolfiles[filet[i,2]],t=filet[i,1],var=c("psxf","psyf","pszf","time"))
			# Rotate to get the ranges of the school in the spherical coordinate system of the echo sounder:
			fishposV=rotate3D(cbind(schoolt$psxf,schoolt$psyf,schoolt$pszf)-matrix(c(vessel$psxv[i],vessel$psyv[i],vessel$pszv[i]),ncol=3,nrow=length(schoolt$psxf),byrow=TRUE),by="z",ang=vessel$rtzv[i],sph.out=TRUE)
			# Store the ranges:
			if(length(dim(fishposV))==1){
				rangesfish[j,]=c(fishposV[1],fishposV[2] %% (2*pi),fishposV[3])
				}
			else{
				rangesfish[j,]=c(range(fishposV[,1]),range(fishposV[,2] %% (2*pi)),range(fishposV[,3]))
				}
			}
		cat("\n")
		# The max beam width
		maxbw=max(beams$bwtl)
		# The ranges of the volume:
		rangethetaB=range(beams$dira)
		rangephiB=range(beams$dire)
		rangesB=c(0, soundbeam_range(beams, pos="max"))
		# Print a warning if the ranges exceed the allowed ranges:
		cat("Max beam width is: ",maxbw,"\n")
		cat("Ranges of the school compared to the ranges of the echo sounder volume:\n")
		cat("Volume:\n")
		cat("\tr\t\t",rangesB[1],"\t",rangesB[2],"\n",sep="")
		cat("\ttheta\t",rangethetaB[1],"\t",rangethetaB[2],"\t","\n",sep="")
		cat("\tphi\t\t",rangephiB[1],"\t",rangephiB[2],"\n",sep="")
		cat("School:\n")
		cat("\tr\t\t",min(rangesfish[,1]),"\t",max(rangesfish[,2]),"\n",sep="")
		cat("\ttheta\t",min(rangesfish[,3])-maxbw,"\t",max(rangesfish[,4])+maxbw,"\n",sep="")
		cat("\tphi\t\t",min(rangesfish[,5])-maxbw,"\t",max(rangesfish[,6])+maxbw,"\n",sep="")
		}		
	
	# Chose origin of the cartesian coordinate system (G) in the geographical coordinate system of Earth:
	orgn=readline("Set the origin of the cartesian coordinate system (G) in decimal longitude and latitude separated by comma\n")
	orgn=as.numeric(unlist(strsplit(orgn,",")))
	if(identical(orgn,NA) || length(orgn)<2){
		vessel$lon0=0
		vessel$lat0=0
		}
	else{
		vessel$lon0=orgn[1]
		vessel$lat0=orgn[2]
		}
	# Transform to geographic coordinates:
	lonvlatv=car2global(vessel)
	vessel[c("lonv","latv")]=list(lonvlatv[,1],lonvlatv[,2])
	# Remove the origin of the global coordinate system, as all information is available through 'latv' and 'lonv':
	vessel=vessel[setdiff(names(vessel),c("lon0","lat0"))]
	
	rm(schoolt)
	gc()
	
	# Generate 'roll' ('rtxv'), pitch (rtyv) and heave (pszv) from the Gaussian distribution if rph is a matrix of two columns of length 3 representing the mean and the standard deviation:
	if(is.list(rph)){
		vessel$rtxv=rnorm(length(vessel$psxv),rph[1,1],rph[1,2])
		vessel$rtyv=rnorm(length(vessel$psxv),rph[2,1],rph[2,2])
		vessel$pszv=rnorm(length(vessel$psxv),rph[3,1],rph[3,2])
		}
	# Add ping index:
	vessel$indt=as.double(1:length(vessel$psxv))
	
	
	##### output #####
	# Remove elements of 'vessel' having names that do not have 4 characters:
	vessel=vessel[nchar(names(vessel))==4]
	# Set the path to the .vessel file:
	vesselpathname=paste(event[1],"/",eventname,"_",paste(zeros(3-nchar(pathnr)),collapse=""),pathnr,".vessel",sep="")
	# Write the vessel data to file:
	write.TSD(x=vessel,con=vesselpathname,numt=length(vessel$psxv),ts=1)
	# Return the vessel dynamics:
	if(scan.volume){
		return(c(vessel,list(ranges=rangesfish,volume=c(rangesB,rangethetaB,rangephiB))))
		}
	else{
		return(vessel)
		}
	##################################################
	##################################################
	}
