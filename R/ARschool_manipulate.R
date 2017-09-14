#*********************************************
#*********************************************
#' Generates and writes to file, fish positions which are AR(1) but not interacting, by the model in Holmin et al. 2012, "Simulations of multi-beam sonar echos from schooling individual fish in a quiet environment".
#'
#' @param con  is a the path to the directory in which the school files to be modified are located.
#' @param t_offset  is the time offset of the manipulation.
#' @param par  is a list or data frame with vectors c(t0, x0,y0,z0, xI,yI,zI, Rx,Ry,Rz, vx,vy,vz, edger,edget), where
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR rotate3D
#' @importFrom TSD car2sph ftim2utim read.TSD write.TSD
#'
#' @export
#' @rdname ARschool_manipulate
#'
ARschool_manipulate<-function(con, t="all", t_offset=0, par=list(t0=ftim2utim(20140101), x0=0,y0=0,z0=0, xI=0,yI=0,zI=0, Rx=12,Ry=3,Rz=2, vx=1,vy=1,vz=1, edger=1,edget=0.2), dir=c(0,pi/2)){
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-03-11 - Finished.
	# Last: 2014-12-12 - Added 't_offset'.
	########### DESCRIPTION: ###########
	# Generates and writes to file, fish positions which are AR(1) but not interacting, by the model in Holmin et al. 2012, "Simulations of multi-beam sonar echos from schooling individual fish in a quiet environment".
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---con--- is a the path to the directory in which the school files to be modified are located.
	# ---t_offset--- is the time offset of the manipulation.
	# ---par--- is a list or data frame with vectors c(t0, x0,y0,z0, xI,yI,zI, Rx,Ry,Rz, vx,vy,vz, edger,edget), where
	#		[[1]] t0 = the start time of the impulse in UNIX time.
	#		[[2]] x0 = the x-position of the point of origin of the orientation wave.
	#		[[3]] y0 = the y-position of the point of origin of the orientation wave.
	#		[[4]] z0 = the z-position of the point of origin of the orientation wave.
	#		[[5]] xI = the x-position of the point of information of the orientation wave, set at infinite distance to simulate equal turning for all fish affected by the wave.
	#		[[6]] yI = the y-position of the point of information of the orientation wave, set at infinite distance to simulate equal turning for all fish affected by the wave.
	#		[[7]] zI = the z-position of the point of information of the orientation wave, set at infinite distance to simulate equal turning for all fish affected by the wave.
	#		[[8]] Rx = the range of the wave in the x-direction before discipating.
	#		[[9]] Ry = the range of the wave in the y-direction before discipating.
	#		[[10]] Rz = the range of the wave in the z-direction before discipating.
	#		[[11]] vx = the speed of the wave in the x-direction.
	#		[[12]] vy = the speed of the wave in the y-direction.
	#		[[13]] vz = the speed of the wave in the z-direction.
	#		[[14]] edger = the decay of the strength factor of the wave as a function of the range r. edger=10 corresponds 10 meters transition between the fish completely affected by the wave and those unaffected.
	#		[[15]] edget = the decay of the information factor of the wave as a function of the range r.
	#		[[16]] rotx = the rotation of the ellipsoid around the x-axis. NOT IMPLEMENTED!
	#		[[17]] roty = the rotation of the ellipsoid around the y-axis. NOT IMPLEMENTED!
	#		[[18]] rotz = the rotation of the ellipsoid around the z-axis. NOT IMPLEMENTED!
	##################################################
	##################################################
	##### Preparation #####
	nimpulses=length(par$t0)
	# Function for calculating the radius measure of ellipsoids, meaning the fractional distance out to the periphery in any direction:
	r_ellipse=function(x,y,z, rx,ry,rz, maxr){
		xsph=car2sph(x,y,z)
		a1 = (xsph[,1]^2*cos(xsph[,2])^2*sin(xsph[,3])^2/rx^2)
		a2 = (xsph[,1]^2*sin(xsph[,2])^2*sin(xsph[,3])^2/ry^2)
		a3 = (xsph[,1]^2*cos(xsph[,3])^2/rz^2)
		(a1+a2+a3) * maxr
		}
	# The mixing factor between the old and the new rotation angle
	mix=function(data,par,t_offset){
		if(data$utim>=par$t0 && data$utim<=par$tmax){
			# Time diference:
			t=data$utim-par$t0+t_offset
			# Get the longest axis of the ellipsoid:
			maxr=max(par$Rx,par$Ry,par$Rz)
			maxrfromt=max(par$vx*t,par$vy*t,par$vz*t)
			# Get the intensity of the wave and the information at the point of origin of the wave, where intensity of the wave is static in time and only depending on the distance to the point of origin, and the information travels with constant speed outwards:
			p0r = 1/2 + maxr/par$edger
			p0t = 1/2 + maxrfromt/par$edget
			
			# Get the current ellipsoid-radius of the fish positions in the range of the wave and in the current propagation of the wave:
			r=r_ellipse(data$psxf-par$x0,data$psyf-par$y0,data$pszf-par$z0, par$Rx,par$Ry,par$Rz, maxr)
			if(maxrfromt>0){
				rfromt=r_ellipse(data$psxf-par$x0,data$psyf-par$y0,data$pszf-par$z0, par$vx*t,par$vy*t,par$vz*t, maxrfromt)
				}
			else{
				rfromt=rep(Inf,length(r))
				}
			# Calculate the mixing of the previous and new direction:
			pr = p0r - r/par$edger
			pr[pr>1]=1
			pr[pr<0]=0
			pt = p0t - rfromt/par$edget
			pt[pt>1]=1
			pt[pt<0]=0
			p = pr*pt
			list(p=p,r=r,rfromt=rfromt,pr=pr,pt=pt)
			}
		else{
			return(NULL)
			}
		
		}
	maxlengthpar=max(sapply(par,length))
	# Repeat to the max length:
	par=lapply(par,rep,length.out=maxlengthpar)
	
	
	##### Execution and output #####
	utim=read.TSD(con,t=t,var="utim")$utim
	# Move through the time steps and manipulate the school:
	newcon=file.path(dirname(con),sub(".","_mod.",basename(con),fixed=TRUE))
	
	# Define parameters used when plotting the time bar for the generation of bottom points:
	infostring="Processing, inducing orientation responses on the school:"
	cat(infostring,"\n",sep="")
	totalsteps=length(utim)
	stepfact=nchar(infostring)/totalsteps
	oldvalue=0
	for(i in seq_along(utim)){
		# Print a dot if the floor of the new value exceeds the old value in:
		thisvalue=floor(i*stepfact)
		if(thisvalue > oldvalue){
			cat(rep(".",thisvalue-oldvalue),if(i==totalsteps) "\n", sep="")
			oldvalue=thisvalue
			}
		# Read the data:
		towrite=read.TSD(con,t=i)
		velxyz=cbind(towrite$vlxf,towrite$vlyf,towrite$vlzf)
		
		# Get the mixing values:
		thisjs=NULL
		for(j in seq_len(nimpulses)){
			thispar=lapply(par,"[[",j)
			m=mix(towrite,thispar)
			if(length(m)>0){
				thisjs=c(thisjs,j)
				# Get the angle between the inpulse point and the fish positions:
				angletoimpulse=car2sph(towrite$psxf-thispar$xI,towrite$psyf-thispar$yI,towrite$pszf-thispar$zI)
				velxyz=rotate3D(velxyz,by="zx",ang=-cbind(angletoimpulse[,2]-dir[1],angletoimpulse[,3]-dir[2])*m$p,paired=TRUE)
				}
			}
		cat(i,thisjs,"\n",sep=", ")
		# Get the velozity and angle data:
		velxyzsph=car2sph(velxyz)
		towrite$vlxf=velxyz[,1]
		towrite$vlyf=velxyz[,2]
		towrite$vlzf=velxyz[,3]
		towrite$rtzf=velxyzsph[,2]-pi/2
		towrite$rtxf=pi/2-velxyzsph[,3]
		
		write.TSD(towrite,newcon,append=i>1)
		}
	newcon
	##################################################
	##################################################
	}
	
