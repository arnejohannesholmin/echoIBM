#*********************************************
#*********************************************
#' Numerical calculation of the surface integral of the function 'fun' on (parts of) the unit sphere.
#'
#' @param fun  is the input function to integrate. Must take a column matrix of theta- and phi-values, or a vector of phi-values as input, and return a vector for the paris of 'theta' and 'phi'.
#' @param ndim  is the number of dimensions of the input to 'fun'.
#' @param theta  is a vector of two elements representing the lower and upper bounds of theta (azimuth angle). Only needed if ndim==2.
#' @param phi  is a vector of two elements representing the lower and upper bounds of the phi (elevation angle).
#' @param R  is the radius of the sphere.
#' @param pres  is the desired presition of the integration.
#' @param l  is the initial lengths of the theta and phi grid.
#' @param incriment  is the factor to multiply the lengths of the grids by.
#' @param equal.size  is TRUE if all cells are to have equal size (speeds up function) and FALSE if the anglular incriment is to be equal.
#' @param max.cells  is the maximum number of cells in the grid.
#' @param print  is TRUE if the improvement of the iteration is to be printed.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname integrateonsphere
#'
integrateonsphere<-function(fun,ndim=2,theta=c(0,2*pi),phi=c(0,pi),R=1,pres=1e-6,l=c(100,100),incriment=2,equal.size=TRUE,max.cells=1e6,print=TRUE){
	
	############### LOG: ###############
	# Start: 2010-01-26 - Clean version.
	
	##### Preparation #####
	# 'iter' is the number of iterations:
	iter=0
	# 'outlast' is the last calculated integral value, and 'out' is the present value:
	outlast=0
	out=0
	# 'diffr' is the difference between the present and the last integral value
	diffr=pres+1
	# The total area of the surface of the spherical segment:
	totalarea=area_shpere(theta=theta,phi=phi)
	# Print the improvements:
	if(print){
		cat("Improvement:\n")
		}
		
		
	##### Execution #####
	# One dimensional function of phi (elevation angle):
	if(ndim==1){
		# Setting the bounds of the integration from 'phi':
		from=phi[1]
		to=phi[2]
		# Warning message if integration bounds exceeds the unit sphere:
		if(any(from<0,to>pi)){
			warning("Integration bounds exceeds phi=c(0,pi)")
			}
		# Loop through iterations, terminated if the desired precision is met or if the maximum number of cells is reached:
		while(diffr>pres && l[1]<=max.cells){
			# Update 'iter':
			iter=iter+1
			# Fastest way:
			if(equal.size){
				# The areas of the grid cells are all equal in this case, and the positions at which the function 'fun' is calculated is found by the get.grid_phi() function (see info(get.grid_phi)).
				# 'phi' are the positions at which the function 'fun' is measured:
				phi=get.grid_phi(from,to,area=c(0.5,double(l[1]-1)+1,0.5))[1:l[1]+1]
				deltaA=totalarea/l[1]
				}
			# Less sensitive to complex structure of 'fun' at the poles (phi approx 0,pi):
			else{
				# In this case 'deltaphi' is constant, while the areas 'deltaA' vary:
				deltaphi=(to-from)/l[1]
				# 'phim' is the gridlines between which the area of the grid cells are calculated, and 'phi' are the positions at which the function 'fun' is measured:
				phim=seq(from,to,deltaphi)
				phi=phim[-1] - deltaphi/2
				deltaA=2*pi*R*(cos(phim[-(l[1]+1)])-cos(phim[-1]))
				}
			# Integration:
			out=sum(fun(phi)*deltaA)
			# Updating 'diffr', 'outlast' and 'l':
			diffr=abs(out-outlast)
			if(print){
				cat("--- ",diffr,", number of cells: ",l[1],"\n")
				}
			outlast=out
			l[1]=abs(l[1]*incriment)
			}
		}
	
	# Two dimensional function of theta and phi (azimuth and elevation angle):
	else{
		# Setting the bounds of the integration from 'theta' and 'phi':
		from=c(theta[1],phi[1])		
		to=c(theta[2],phi[2])
		# Warning message if integration bounds exceeds the unit sphere:
		if(any(from<0,to>c(2*pi,pi))){
			warning("Integration bounds exceeds theta=c(0,2*pi) and phi=c(0,pi)")
			}
		# Loop through iterations, terminated if the desired precision is met or the maximum number of cells is reached:
		while(diffr>pres && prod(l)<=max.cells){
			# The grid of 'theta':
			deltatheta=(to[1]-from[1])/l[1]
			theta=seq(from[1],to[1],deltatheta)[-1] - deltatheta/2
			# Update 'iter':
			iter=iter+1
			# Fastest way:
			if(equal.size){
				# The areas of the grid cells are all equal in this case, and the positions at which the function 'fun' is calculated is found by the get.grid_phi() function (see info(get.grid_phi)).
				# 'phi' are the positions at which the function 'fun' is measured:
				phi=get.grid_phi(from[2],to[2],area=c(0.5,double(l[2]-1)+1,0.5))[1:l[2]+1]
				deltaA=totalarea/prod(l)
				}
			# Less sensitive to complex structure of 'fun' at the poles (phi approx 0,pi):
			else{
				# In this case 'deltaphi' is constant, while the areas 'deltaA' vary:
				deltaphi=(to[2]-from[2])/l[2]
				# 'phim' is the gridlines between which the area of the grid cells are calculated, and 'phi' are the positions at which the function 'fun' is measured:
				phim=seq(from[2],to[2],deltaphi)
				phi=phim[-1] - deltaphi/2
				deltaA=outer( rep(deltatheta,l=l[1]) , R*(cos(phim[-(l[1]+1)])-cos(phim[-1])) )
				}
			# Integration:
			out=sum(fun(expand.grid(theta,phi))*deltaA)
			# Updating 'difr', 'outlast', 'xgrid' and 'xgrid':
			diffr=abs(out-outlast)
			if(print){
				cat("--- ",diffr,", number of cells: ",prod(l),"\n")
				}
			outlast=out
			l=abs(l*incriment)
			}
		}
	if(iter==1){
		warning("Only one iteration done. The result may be wrong. Try increasing the initial number of bins of 'theta' and 'phi' by increasing the value of 'l'")
		}
	
	
	##### Output #####
	list(out=out,iter=iter)
}
