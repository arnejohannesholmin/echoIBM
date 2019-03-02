#*********************************************
#*********************************************
#' Calculates the (volume) size of the spherical volume segment defined by 'r', 'theta' and 'phi'.
#'
#' @param r  is range of radial coordinates, given as a vector of interval values (in which case length(r)-1 spherical segemnts are treated) or as a two column matrix where the first column represents the lower limit of the sperical segments, and the second column represents the upper limit.
#' @param theta  is the range of azimuth angles in radians, given in the same way as 'r'.
#' @param phi  is the range of elevation angles in radians, given in the same way as 'r.
#' @param var  is a vector of the variables to return. Currently implemented are "volx" for volumes of the voxels and "harx" for horizontal area of the voxels.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD repm strff
#'
#' @export
#' @rdname vol.sph
#'
vol.sph <- function(r=c(0,1), theta=c(0,2*pi), phi=c(0,pi), var=c("volx", "harx")){
		
	############### LOG: ###############
	# Start: 2011-06-26 - Clean version.
	# Last: 2011-09-23 - Added support for multiple volumes, given as two columen matrices or as vectors where the intervals are between consecutive points.
	
	##### Preparation #####
	# Warnings for values not in r>=0, theta = (0,2*pi) or phi = (0,pi):
	if(length(r) && min(r,na.rm=TRUE)<0){
		warning("Negative values detected in 'r', leading to subtraction of the volume for which r<0")
		}
	if(length(theta) && min(theta,na.rm=TRUE)<0){
		warning("Negative values detected in 'theta', leading to subtraction of the volume for which theta<0")
		}
	if(length(theta) && max(theta,na.rm=TRUE)>2*pi){
		warning("Values detected for which theta>2*pi, leading to addition of the volume for which theta>2*pi")
		}
	if(length(phi) && min(phi,na.rm=TRUE)<0){
		warning("Negative values detected in 'phi', leading to subtraction of the volume for which phi<0")
		}
	if(length(phi) && max(phi,na.rm=TRUE)>pi){
		warning("Values detected for which phi>pi, leading to addition of the volume for which phi>pi")
		}
	
	# Accept if only one value is given:
	if(length(r)==1){
		r=c(r,r)
		}
	if(length(theta)==1){
		theta=c(theta,theta)
		}
	if(length(phi)==1){
		phi=c(phi,phi)
		}
	# If given as vectors, transform to matrices of intervals:
	if(length(dim(r))==0){
		r=cbind(r[-length(r)],r[-1])
		}
	if(length(dim(theta))==0){
		theta=cbind(theta[-length(theta)],theta[-1])
		}
	if(length(dim(phi))==0){
		phi=cbind(phi[-length(phi)],phi[-1])
		}
		
	# Ensure that 'r', 'theta' and 'phi' have equal dimension:
	if(diff(range(c(nrow(r),nrow(theta),nrow(phi)))) <= 0){
		m=max(c(nrow(r),nrow(theta),nrow(phi)))
		r=repm(r,length.out=m,byrow=TRUE)
		theta=repm(theta,length.out=m,byrow=TRUE)
		phi=repm(phi,length.out=m,byrow=TRUE)
		}
	
	
	##### Execution #####
	out=list()
	if(strff("volx", var[1])){
		# The volume of the spherical segments, calculated by the following argument:
		# Volume of an open spherical sector (see http://mathworld.wolfram.com/SphericalSector.html), is 2 pi r^3(cos(phi1) - cos(phi2))/3, where 'r' is the radial distance to the edge of the sphere, and 'theta' and 'phi' are the azimuth and elevation angles, respectively. Subtracting for two radial distances and subsetting the resulting strip by theta relative to 2*pi gives: 
		# (theta2-theta1) (cos(phi1) -cos(phi2))(r2^3 - r1^3) / 3:
		out$volx = (theta[,2]-theta[,1])/3 * (cos(phi[,1])-cos(phi[,2])) * (r[,2]^3-r[,1]^3)
		}
	if(strff("harx", var[1])){
		# The horizontal area of the voxels is calculated so that projected onto the surface, the voxels in one layer of horizontally aligned voxels are disjoint and separated at the horizontal lines passing through the acoustic axis of the beams coresponding to the voxels at each sphere separating the voxels radially:
		out$harx = (r[,2]^2-r[,1]^2) * (cos((phi[,1]+phi[,2])/2))^2 * (theta[,2]-theta[,1])/2
		}
	
	
	##### Output #####
	out
}
