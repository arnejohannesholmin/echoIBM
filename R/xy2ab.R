#*********************************************
#*********************************************
#' Parametrizes the piecewise linear trace given by the input. For vertical lines, $Inf$ is returned.
#'
#' @param x  and 'y' are the x-coordinates and y-coordinates of the input trace, given as separate vectors 'x' and 'y', or as a list or matrix 'x' holding both the x-coordinates and the y-coordinates, or if y==NULL 'x' is interpreted as 'y' and x=seq_along(x).
#' @param ang.out  is TRUE if angles are to be returned instead of slope.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname xy2ab
#'
xy2ab<-function(x,y=NULL,ang.out=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-10-29 - Finished.
	# Update:  2009-05-05 - Cleaned up and added support for differing lengths of 'x' and 'y'. Name changed from paramline() to xy2ab().
	# Update:  2009-06-10 - Added support for list input and matrix input, and for angle output instead of slope output.
	# Update:  2009-07-29 - Changed to correspond to functions like car2pol(), where list input gives data.frame output and matrix and vector input give matrix output.
	# Last:  2010-08-26 - Removed the data.frame output.
	########### DESCRIPTION: ###########
	# Parametrizes the piecewise linear trace given by the input. For vertical lines, $Inf$ is returned.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- and 'y' are the x-coordinates and y-coordinates of the input trace, given as separate vectors 'x' and 'y', or as a list or matrix 'x' holding both the x-coordinates and the y-coordinates, or if y==NULL 'x' is interpreted as 'y' and x=seq_along(x).
	# ---ang.out--- is TRUE if angles are to be returned instead of slope.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Return list if input is a list:
	listinputx=FALSE
	# Support for vector, matrix and list input for 'x':
	if(is.list(x)){
		names(x)=tolower(names(x))
		if(!is.null(x$x) && !is.null(x$y)){
			y=x$y
			x=x$x
			}
		else{
			y=x[[2]]
			x=x[[1]]
			}
		listinputx=TRUE
		}
	else if(is.null(y)){
		dimx=dim(x)
		if(length(dimx)==2){
			if(dimx[2]==1){
				y=drop(x)
				x=seq_along(x)
				}
			else{
				y=x[,2]
				x=x[,1]
				}
			}
		# Add zeros for the 'y' values:
		else if(is.null(dimx)){
			y=x
			x=seq_along(x)
			}
		else{
			stop("Invalid input")
			}
		}
	# 'x' and 'y' need to have equal length:
	lx=length(x)
	ly=length(y)
	if(lx!=ly){
		stop("'x' and 'y' lengths differ")
		}
	# Giving 50% longer time usage for length(x)=length(y)=1e6:
	#xy=xy.coords(x,y)
	#x=xy$x
	#y=xy$y
	#lx=length(x)
		
	
	##### Execution #####
	diffx=diff(x)
	diffy=diff(y)
	b=diffy/diffx
	a=y[-lx]-x[-lx]*b
	
	
	##### Output #####
	if(ang.out && listinputx){
		list(a=a,ang=atan2(diffy,diffx))
		}
	else if(ang.out && !listinputx){
		cbind(a=a,ang=atan2(diffy,diffx))
		}
	else if(listinputx){
		list(a=a,b=b)
		}
	else{
		cbind(a=a,b=b)
		}
	##################################################
	##################################################
	}
