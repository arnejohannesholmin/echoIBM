#*********************************************
#*********************************************
#' Interactive function for selecting the path and speed of an object based on reference points 'x', intended for use in echoIBM().
#'
#' @param x  are the reference points.
#' @param n  and 't' are the numbers of points and the times for the points along the path respectively. The combination of 'n' and 't', and whether the user selects speed values interactively, defines the mode of the output. There are 7 modes of the output, as liste below (lt=length(t)):
#' @param ref  is a list of reference information from a previous set.path.vessel session (given as the output from set.path.vessel()). The points of ref will be plotted in grey, to form a visual reference for the user.
#' @param margin  is a vector of the margins on either side of the span of 'x' (recycled if not of length 4).
#' @param maxn  maximum number of points entered (to prevent "getting stuck" in the function).
#' @param path.col  color of the entered path.
#' @param grid  is TRUE if a 100 meter grid is to be plotted.
#' @param smooth  has one possible value "spline", smooths the x-values and the y-values of the path separately using the spline function.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD setrange
#'
#' @importFrom utils tail
#' @importFrom graphics plot abline points lines text locator par arrows
#' @importFrom grDevices rainbow col2rgb rgb
#' @importFrom stats spline smooth.spline
#'
#' @export
#' @rdname set.path.vessel
#'
set.path.vessel<-function(x,n=nrow(x),t=NULL,ref=NULL,maxspeed=15,margin=500,maxn=100,path.col="blue",grid=TRUE,smooth="spline"){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-05-08 - Finished.
	########### DESCRIPTION: ###########
	# Interactive function for selecting the path and speed of an object based on reference points 'x', intended for use in echoIBM().
	########## DEPENDENCIES: ###########
	# rounded.rect(), track.line(), xy2ab()
	############ VARIABLES: ############
	# ---x--- are the reference points.
	# ---n--- and 't' are the numbers of points and the times for the points along the path respectively. The combination of 'n' and 't', and whether the user selects speed values interactively, defines the mode of the output. There are 7 modes of the output, as liste below (lt=length(t)):
	#
	#			'n'		't'		'speed'		Output
	#	--------------------------------------------------------------------------------------------------
	#		1	NULL	NULL	not given	The entered points 'pos' that define the path
	#		2	'n'		NULL	not given	'n' equally spaced points along the path defined by 'pos'
	#		3	NULL	NULL	given		nrow(x) points along the path, spaced according to time points 
	#											separated by 1 sec, for equally separated speed values along
	#											the given speed graph
	#		4	'n'		NULL	given		'n' points along the path, spaced according to time points 
	#											separated by 1 sec, for equally separated speed values along
	#											the given speed graph
	#		5	NULL	lt==1	given		nrow(x) points along the path, spaced according to time points 
	#											separated by 't' sec, for equally separated speed values along
	#											the given speed graph
	#		6	'n'		lt==1	given		'n' points along the path, spaced according to time points 
	#											separated by 't' sec, for equally separated speed values along
	#											the given speed graph
	#		7	any		lt>1	given		length(t) points along the path, spaced according to time points 
	#											't', for speed values according to the time points 't' along
	#											the given speed graph
	#	--------------------------------------------------------------------------------------------------
	#
	# ---ref--- is a list of reference information from a previous set.path.vessel session (given as the output from set.path.vessel()). The points of ref will be plotted in grey, to form a visual reference for the user.
	# ---margin--- is a vector of the margins on either side of the span of 'x' (recycled if not of length 4).
	# ---maxn--- maximum number of points entered (to prevent "getting stuck" in the function).
	# ---path.col--- color of the entered path.
	# ---grid--- is TRUE if a 100 meter grid is to be plotted.
	# ---smooth--- has one possible value "spline", smooths the x-values and the y-values of the path separately using the spline function.
		
	
	##################################################
	##################################################
	##### Preparation #####
	# According to the convension on sea, knots are used as the speed measure:
	knot=1.852/3.6
	### Plot the points of 'x' and the enter button: ###
	# How much space that is reserved for the speed plot and the buttons, as a fraction of 'margin':
	speedpanelspace=1/2
	xlim=range(x[,1])+c(-margin,margin)
	ylim=range(x[,2])+c(-(1+speedpanelspace)*margin,margin)
	# Plot the points 'x':
	plot(x[,1:2],type="l",xlim=xlim,ylim=ylim,xlab="x",ylab="y")
	if(grid){
		abline(v=seq.int(ceiling(xlim[1]/100),floor(xlim[2]/100))*100,col="grey")
		abline(h=seq.int(ceiling(ylim[1]/100),floor(ylim[2]/100))*100,col="grey")
		}
	nrowx=nrow(x)
	points(x[,1:2],pch=".",cex=4,xlim=xlim,ylim=ylim,xlab="x",ylab="y",col=rainbow(nrowx,start=0,end=0.8))
	text(x[1,1],x[1,2],labels=1,cex=0.8)
	text(x[nrowx,1],x[nrowx,2],labels=nrowx,cex=0.8)
	# Position and size of the return button:
	o_return=c(xlim[1]+diff(xlim)*0.04,ylim[1]+diff(ylim)*0.04)
	w=diff(xlim)*0.1
	h=diff(ylim)*0.08
	# Plot the return button:
	rounded.rect(o_return,w,h,r=0.3,col="grey44",border="grey55",adjust=TRUE)
	text(xlim[1]+diff(xlim)*0.04,ylim[1]+diff(ylim)*0.04,"return",vfont=c("sans serif","bold"))
	# Position and size of the speed button:
	o_speed=c(xlim[1]+diff(xlim)*0.17,ylim[1]+diff(ylim)*0.04)
	# Plot the speed button:
	rounded.rect(o_speed,w,h,r=0.3,col="grey44",border="grey55",adjust=TRUE)
	text(xlim[1]+diff(xlim)*0.17,ylim[1]+diff(ylim)*0.04,"speed",vfont=c("sans serif","bold"))
	# If reference path is given as an element named 'refpos' in 'ref', this path is plotted in grey:
	if(!is.null(ref$refpos)){
		points(ref$refpos,pch=".",cex=5,col="grey")
		lines(ref$refpos,col="grey")
		}
	
	
	##### Execution #####
	### Positions defining the path: ###
	pos=NULL
	# First point:
	p1=1
	pos=rbind(pos,unlist(locator(1)))
	isenter=o_return[1]-w/2<pos[p1,1] & o_return[1]+w/2>pos[p1,1] & o_return[2]-h/2<pos[p1,2] & o_return[2]+h/2>pos[p1,2]
	if(isenter){
		return(NULL)
		}
	isspeed=FALSE
		
	# Loop through the points entered. If the user presses the speed button or the return button the loop is terminated:
	while(p1<maxn && !isenter && !isspeed){
		points(pos[p1,1],pos[p1,2],col=path.col,pch=".",cex=5)
		if(p1>1){
			lines(pos[p1+c(-1,0),],col=path.col)
			}
		p1=p1+1
		pos=rbind(pos,unlist(locator(1)))
		isenter=o_return[1]-w/2<pos[p1,1] & o_return[1]+w/2>pos[p1,1] & o_return[2]-h/2<pos[p1,2] & o_return[2]+h/2>pos[p1,2]
		isspeed=o_speed[1]-w/2<pos[p1,1] & o_speed[1]+w/2>pos[p1,1] & o_speed[2]-h/2<pos[p1,2] & o_speed[2]+h/2>pos[p1,2]
		}
	# Plot the result in stonger color:
	endcol=col2rgb(path.col)/3
	lines(pos[-p1,,drop=FALSE],col=rgb(endcol[1],endcol[2],endcol[3],maxColorValue=255),lwd=1.5)
		
	### Choose speed: ###
	# Size of the plotting window:
	din=mean(par()$din)
	# First point:
	isspeed=o_speed[1]-w/2<pos[p1,1] & o_speed[1]+w/2>pos[p1,1] & o_speed[2]-h/2<pos[p1,2] & o_speed[2]+h/2>pos[p1,2]
	p2=0
	# The points defining the speed:
	xwin=NULL
	ywin=NULL
	# Because the speed points are transformed later, the original speed points are stored and returned in the output list:
	originalxwin=NULL
	originalywin=NULL
	isenter=FALSE
	if(isspeed){
		# If reference speed is given as an element named 'refspeed' in 'ref', these speed points are is plotted in grey:
		if(!is.null(ref$refspeed)){
			points(ref$refspeed,pch=".",cex=5,col="grey")
			lines(ref$refspeed,col="grey")
			}
		# The speed plotting is terminated if the user presses the return button:
		while(!isenter){
			# Plotting the speed panel:
			if(p2==0){
				x0=xlim[1]+diff(xlim)*0.30
				x1=xlim[1]+diff(xlim)*0.95
				y0=ylim[1]+diff(ylim)*0.015
				y1=ylim[1]+2*margin*speedpanelspace
				arrows(c(x0,x0),c(y0,y0),c(x1,x0),c(y0,y1),length=din/50)
				# Text on the axes:
				text(x0-(x1-x0)/15,y1,"knots",cex=0.6)
				text(x1,y0-(y1-y0)/15,"t",cex=0.6)
				# Ticks on the axes:
				tick=10^(nchar(maxspeed))*c(0.1,0.2,0.4,0.5)
				tickind=which.min(abs(1-maxspeed/tick))
				step=tick[tickind]/10
				nsteps=floor(maxspeed/tick[tickind]*10)
				newmaxspeed=(nsteps+1)*step
				ticks=y0+1:nsteps/(nsteps+1)*(y1-y0)
				for(i in seq_along(ticks)){
					lines(x0+c(-diff(xlim)*0.01,diff(xlim)*0.01),rep(ticks[i],2))
					}
				for(i in 1:(ceiling(nsteps/2)-1)*2){
					text(x0-diff(xlim)*0.03,ticks[i],step*i,cex=0.6)
					}
				}
			point=locator(1)
			isenter=o_return[1]-w/2<point$x & o_return[1]+w/2>point$x & o_return[2]-h/2<point$y & o_return[2]+h/2>point$y
			# If the user clicks on the speed button, nothing happens:
			isspeed0=o_speed[1]-w/2<point$x & o_speed[1]+w/2>point$x & o_speed[2]-h/2<point$y & o_speed[2]+h/2>point$y
			# If the user mouse clicks outside of the speed panel, these points are interpreted as edge points:
			if(!isenter && !isspeed0){
				if(point$x < x0){
					xwin=append(xwin,x0)
					}
				else if(point$x > x1){
					xwin=append(xwin,x1)
					}
				else{
					xwin=append(xwin,point$x)
					}
					
				if(point$y < y0){
					ywin=append(ywin,y0)
					}
				else if(point$y > y1){
					ywin=append(ywin,y1)
					}
				else{
					ywin=append(ywin,point$y)
					}
					
				# Plot the selected point:
				points(xwin[p2+1],ywin[p2+1],col="skyblue2",pch=".",cex=5)
				if(p2>0){
					lines(xwin[p2+0:1],ywin[p2+0:1],col="skyblue2")
					}
				p2=p2+1
				}
			}
		# Sorting the speed values with regard to the x-values, and plotting the end result:
		if(is.null(xwin)){
			isspeed=FALSE
			}
		else{
			orden=order(xwin)
			xwin=xwin[orden]
			ywin=ywin[orden]
			originalxwin=xwin
			originalywin=ywin
			lines(xwin,ywin,col="blue")
			# Transfoming the speed points to knots:
			xwin=setrange(xwin)
			ywin=(ywin-y0)/(y1-y0)*newmaxspeed
			}
		}
	# Strip 'pos' of the last row:
	pos=pos[-nrow(pos),,drop=FALSE]
	# Smooth the path:
	if(smooth[1]=="spline"){
		originalpos=pos
		lspline=nrow(pos)*10
		pos=cbind(spline(pos[,1],n=lspline)$y,spline(pos[,2],n=lspline)$y)
		}
	else if(smooth[1]=="smooth.spline"){
		originalpos=pos
		pos=cbind(smooth.spline(pos[,1])$y,smooth.spline(pos[,2])$y)
		}
	# 't' is missing, it is set to 1 second as a dummy variable, and points spaced accordint to the equally spaced speed values are returned:
	normalized=FALSE
	if(is.null(t)){
		t=1
		normalized=TRUE
		}
		
		
	##### output #####
	# If speed values are not chosen:
	if(!isspeed){
		# Extract the first time difference 'dt' of 't' if 't' has length > 1, and set 'dt' to 1 otherwise:
		if(length(t)>1){
			dt=diff(t)[1]
			}
		else{
			dt=1
			}
		# Mode 1:
		if(is.null(n)){
			out=originalpos
			# If only one point is selected, heading is defaulted to 0:
			if(nrow(originalpos)<=1){
				rtzv=0
				}
			# Subtract pi/2 to account for the difference in definition between polar angle and z-rotation (theta=0 ~ rtzv?-pi/2):
			else{
				rtzv=xy2ab(originalpos,ang.out=TRUE)[,2]-pi/2
				# Add the last element of the rotation vector:
				rtzv=c(rtzv,tail(rtzv,1))
				}
			# Set the 'speed', assuming dt=1 seconds:
			speed=sqrt(diff(originalpos[,1])^2+diff(originalpos[,2])^2)/dt/knot
			speed=c(speed,tail(speed,1))
			}
		# Mode 2:
		else{
			n=n[1]
			out=track.line(pos,w=n)
			# If only one point is selected, heading is defaulted to 0:
			if(nrow(pos)<=1){
				rtzv=0
				}
			# Subtract pi/2 to account for the difference in definition between polar angle and z-rotation (theta=0 ~ rtzv?-pi/2):
			else{
				rtzv=xy2ab(pos,ang.out=TRUE)[out$pos,2]-pi/2
				}
			# Set the 'speed', assuming dt=1 seconds:
			speed=sqrt(diff(out$x)^2+diff(out$y)^2)/dt/knot
			speed=c(speed,tail(speed,1))
			}
		}
	# Else if speed values are chosen:
	else{
		# Mode 3, 5 or 7:
		if(is.null(n)){
			n=nrow(x)
			}
		# Mode 4, 6 or 7:
		else{
			n=n[1]
			}
		# Mode 7:
		if(length(t)>1){
			dt=c(0,diff(t))
			speed=track.line(xwin,ywin,w=t,alongx=TRUE,normalized=TRUE)
			}
		# Mode 3, 4, 5 or 6:
		else{
			speed=track.line(xwin,ywin,w=n,alongx=TRUE,normalized=TRUE)
			dt=rep(t,length.out=nrow(speed))
			}
		# Track the output points along the pathe defined by 'pos':
		out=track.line(pos,w=cumsum(speed$y*knot*dt),normalized=normalized)
		# If only one point is selected, heading is defaulted to 0:
		if(nrow(pos)<=1){
			rtzv=0
			}
		# Subtract pi/2 to account for the difference in definition between polar angle and z-rotation (theta=0 ~ rtzv?-pi/2):
		else{
			rtzv=xy2ab(pos,ang.out=TRUE)[out$pos,2]-pi/2
			}
		}
	# Plot the selected points in red:
	lines(pos,col="green",lwd=1)
	points(out,pch="*",col="red")
	
	### Return: ###
	# If speed is chosen by the user, the speed values are returned as 'ispv' and 'ctim'. 'speed' is simply the reference speeds to be used as reference to the function later:
	null=double(length(out[[1]]))
	if(!isspeed){
		list(psxv=out$x,psyv=out$y,pszv=null,rtxv=null,rtyv=null,rtzv=rtzv,ispv=speed,ctim=rep(t,length.out=length(out[[1]])),refpos=originalpos)
		}
	# Else if speed values are chosen:
	else{
		list(psxv=out$x,psyv=out$y,pszv=null,rtxv=null,rtyv=null,rtzv=rtzv,ispv=speed$y,ctim=cumsum(out$dt),refpos=originalpos,refspeed=cbind(originalxwin,originalywin))
		}
	##################################################
	##################################################
	}
