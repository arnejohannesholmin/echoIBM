#*********************************************
#*********************************************
#' Returns a list of x-values and y-values of points defined by the fractions 'w' along the partial linear path defined by the points (x,y). Cumulative lengths of the path are also returned.
#'
#' @param x  is the input x-variable.
#' @param y  is the input y-variable.
#' @param w  are the positions along the path (x,y). If 'w' is a single integer, it is set to seq(0,1,length.out=w). If 'speed' is given 'w' is interpreted as time values.
#' @param dw  are the piecewise lengths along the path (x,y), used only if w!=NULL. If 'dw' has length>1 it overrides 'w' by w=cumsum(dw).
#' @param speed  are the speed values along the path (x,y), or along the x-axis if alongx==TRUE. If given it overrides 'w'.
#' @param start  is a vector of two values representing the start point if 'x' and 'y' represent angles and lengths. If 'start' is given, 'x' and 'y' are interpreted as angles and lengths.
#' @param alongx  is TRUE if tracking is to be done along the x-axis, and not along the path (x,y).
#' @param normalized  is TRUE if the positions given by 'w' (and 'speed') are or should be normalized to [0,1], where 1 refers to the end of the path.
#' @param list.out  is TRUE if the output should be a data frame with the names "x", "y" and "l". Else a matrix with the same names is returned.
#' @param plot  is TRUE if a simple plot is to be drawn.
#' @param ...  are parameters to be passed on to plot() if plot==TRUE.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD setrange
#' @importFrom utils tail
#' @importFrom graphics hist points
#'
#' @export
#' @rdname track.line
#'
track.line<-function(x,y=NULL,w=NULL,dw=NULL,speed=NULL,start=NULL,alongx=FALSE,normalized=FALSE,list.out=TRUE,plot=FALSE,...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-10-29 - First clean version.
	# Update: 2009-05-07 - Function cleaned up and effectivized. Changed to accepting/transforming w to the interval [0,1], or if a single value, an equispaced grid of length 'w' in the interval.
	# Update: 2009-06-10 - Added support for start point 'start' and angles and lengths of the line segments represented by 'x' and 'y'. Function restructured and simplified in addition to adding support for list input and matrix input.
	# Update: 2009-08-02 - Fixed bug when alongx=TRUE, in which case the calculation of the output must be a little different (see output section).
	# Update: 2010-02-19 - Added suport for time differences 'dw' if speed is given.
	# Last:  2010-08-26 - Removed the data.frame output.
	########### DESCRIPTION: ###########
	# Returns a list of x-values and y-values of points defined by the fractions 'w' along the partial linear path defined by the points (x,y). Cumulative lengths of the path are also returned.
	########## DEPENDENCIES: ###########
	# xy2ab(), setrange()
	############ VARIABLES: ############
	# ---x--- is the input x-variable.
	# ---y--- is the input y-variable.
	# ---w--- are the positions along the path (x,y). If 'w' is a single integer, it is set to seq(0,1,length.out=w). If 'speed' is given 'w' is interpreted as time values.
	# ---dw--- are the piecewise lengths along the path (x,y), used only if w!=NULL. If 'dw' has length>1 it overrides 'w' by w=cumsum(dw).
	# ---speed--- are the speed values along the path (x,y), or along the x-axis if alongx==TRUE. If given it overrides 'w'.
	# ---start--- is a vector of two values representing the start point if 'x' and 'y' represent angles and lengths. If 'start' is given, 'x' and 'y' are interpreted as angles and lengths.
	# ---alongx--- is TRUE if tracking is to be done along the x-axis, and not along the path (x,y).
	# ---normalized--- is TRUE if the positions given by 'w' (and 'speed') are or should be normalized to [0,1], where 1 refers to the end of the path.
	# ---list.out--- is TRUE if the output should be a data frame with the names "x", "y" and "l". Else a matrix with the same names is returned.
	# ---plot--- is TRUE if a simple plot is to be drawn.
	# ---...--- are parameters to be passed on to plot() if plot==TRUE.
	
	
	##################################################
	##################################################
	# Support for vector, matrix and list input for 'x':
	if(is.list(x) && length(x)>1){
		if(!is.null(x$x) && !is.null(x$y)){
			y=x$y
			x=x$x
			}
		else if(!is.null(x$ang) && !is.null(x$l)){
			y=x$l
			x=x$ang
			}
		else{
			y=x[[2]]
			x=x[[1]]
			}
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
	# If only one point is given as input, there is no line to track along, a list/matrix of max(length(w),1) points identical to (x,y) is returned:
	if(lx==1){
		xywout=cbind(rep(x,length.out=max(length(w),1)),y,0,1)
		colnames(xywout)=c("x","y","l","pos")
		if(list.out){
			return(list(x=xywout[,1],y=xywout[,2],l=xywout[,3],pos=xywout[,4]))
			}
		else{
			return(xywout)
			}
		}
	
		
	##### Execution #####
	# If 'start' is given, 'x' and 'y' are interpreted as angles and lengths.
	if(!is.null(start)){
		ang=x
		diffx=y*cos(ang)
		diffy=y*sin(ang)
		if(alongx){
			lengths=diffx
			}
		else{
			lengths=y
			}
		x=start[1]+c(0,cumsum(diffx))
		y=start[2]+c(0,cumsum(diffy))
		lx=lx+1
		}
	else{
		# Sucsessive differences between x- and y-values:
		diffx=diff(x)
		diffy=diff(y)
		# Parametrization of the line segments, needed to calculate the positions on the line segments:
		ang=xy2ab(list(x,y),ang.out=TRUE)$ang
		# Calculating lengths and cummulative lengths along the path or along the x-axis:
		if(alongx){
			lengths=abs(diffx)
			}
		else{
			lengths=sqrt(diffx^2+diffy^2)
			}
		}
	
	# Cummulative lengths along the track:
	cumlengths=c(0,cumsum(lengths))
	totlengths=cumlengths[lx]
	# Scaling 'w' to match 'cumlengths':
	
	# If the length of 'dw' exceeds 1, it overrides 'w':
	if(length(dw)>1){
		w=cumsum(dw)
		}
	# If w==NULL, regularly spaced points along the path are returned:
	else if(is.null(w)){
		w=seq(0,1,length.out=lx)
		normalized=TRUE
		}
	# If 'w' has length 1 regularly spaced points are returned:
	else if(length(w)==1){
		if(!is.null(dw)){
			w=seq(0,by=dw,length.out=w)
			}
		else{
			w=seq(0,1,length.out=w)
			normalized=TRUE
			}
		}
	# 'w' needs to be sorted:
	else{
		w=sort(w)
		}
	
	# If 'speed' is given, the calculation is done on times not lengths:
	if(!is.null(speed)){
		t=w
		# The length variable to be used in the function:
		times=lengths/speed
		cumtimes=c(0,cumsum(times))
		tottimes=tail(cumtimes,1)
		# Scaling 'w' to match 'cumlengths':
		if(normalized){
			t=setrange(t)*tottimes
			}
		# 'w' must be non-negative:
		if(min(t)<0 || max(t)>tottimes){
			warning("Elements of 'w' (time) exceed the path defined by the given circles, and are ignored")
			t=t[t>=0 & t<=tottimes]
			}
		
		# 'pos' gives in which line segments the points defined by 'w' are positioned. If 'totheend' is TRUE, special attention is needed in the calculation to avoid an error:
		pos=hist(t,breaks=cumtimes,plot=FALSE)$counts
		pos=rep(1:(lx-1),pos)
		# Length of the line segments to the left of the points:
		lefttime=t-cumtimes[pos]
		
		w=cumlengths[pos]+abs(lefttime*speed[pos])
		# The output:
		if(alongx){
			posneg=sign(diffx)
			xywout=cbind(x[pos]+lefttime*speed[pos]*posneg[pos],y[pos]+lefttime*speed[pos]*tan(ang[pos])*posneg[pos],w,t,pos)
			}
		else{
			xywout=cbind(x[pos]+lefttime*speed[pos]*cos(ang[pos]),y[pos]+lefttime*speed[pos]*sin(ang[pos]),w,t,pos)
			}
		colnames(xywout)=c("x","y","l","t","pos")
		}
	else{
		# Scaling 'w' to match 'cumlengths':
		if(normalized){
			w=setrange(w)*totlengths
			}
		# 'w' must be non-negative:
		if(min(w)<0 || max(w)>totlengths){
			warning("Elements of 'w' exceed the path defined by the given circles, and are ignored")
			w=w[w>=0 & w<=totlengths]
			}
		
		# 'pos' gives in which line segments the points defined by 'w' are positioned. If 'totheend' is TRUE, special attention is needed in the calculation to avoid an error:
		pos=hist(w,breaks=cumlengths,plot=FALSE)$counts
		pos=rep(1:(lx-1),pos)
		# Length of the line segments to the left of the points:
		leftdist=w-cumlengths[pos]
		
		# The output:
		if(alongx){
			posneg=sign(diffx)
			xywout=cbind(x[pos]+leftdist*posneg[pos],y[pos]+leftdist*tan(ang[pos])*posneg[pos],w,pos)
			}
		else{
			xywout=cbind(x[pos]+leftdist*cos(ang[pos]),y[pos]+leftdist*sin(ang[pos]),w,pos)
			}
		colnames(xywout)=c("x","y","l","pos")
		}
	
	
	##### Output #####
	# Ploting the result:
	if(plot){
		plot(x,y,type="l",...)
		points(xywout,pch="*",col=2)
		}
	# The output:
	if(list.out){
		if(ncol(xywout)==4){
			list(x=xywout[,1],y=xywout[,2],l=xywout[,3],pos=xywout[,4])
			}
		else if(ncol(xywout)==5){
			list(x=xywout[,1],y=xywout[,2],l=xywout[,3],t=xywout[,4],pos=xywout[,5])
			}
		}
	else{
		xywout
		}
	##################################################
	##################################################
	}
