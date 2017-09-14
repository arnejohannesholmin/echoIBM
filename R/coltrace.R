#*********************************************
#*********************************************
#' Draws a color spectrum on which the user can assign points with the mouse pointer defining a color trace of desired length.
#'
#' @param n  is the length of the output color vector.
#' @param ref  is input from a previous session of coltrace(), if the user wishes to adjust an old session.
#' @param coldens  is the number of color lines to be plotted in the color panel. Large value of 'coldens' is more time demanding but looks nicer.
#' @param ...  methods passed on to sub-functions (not used, but allowing for unused arguments).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR contains
#' @importFrom stats quantile
#' @importFrom grDevices rgb colors
#' @importFrom graphics par plot text points lines polygon locator
#'
#' @export
#' @rdname coltrace
#'
coltrace<-function(n=20, ref=NULL, coldens=20, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-10-30 - Beta version.
	# Udate: 2008-11-03 - Some changes.
	# Update:  2009-02-12 - Visual clean up. Some errors found and corrected. Cooperation with other functions like 'contains' must be checked. Function not working properly!!!!
	# Last: 2009-05-12 - Function body changed to be based on one mouse click at the time, to ensure good performance on every platform. Buttons changed to "undo", "redo", "next"/"prev" and "return".
	########### DESCRIPTION: ###########
	# Draws a color spectrum on which the user can assign points with the mouse pointer defining a color trace of desired length.
	########## DEPENDENCIES: ###########
	# rounded.rect(), track.line(), clamp(), contains()
	############ VARIABLES: ############
	# ---n--- is the length of the output color vector.
	# ---ref--- is input from a previous session of coltrace(), if the user wishes to adjust an old session.
	# ---coldens--- is the number of color lines to be plotted in the color panel. Large value of 'coldens' is more time demanding but looks nicer.
	# ---...--- methods passed on to sub-functions (not used, but allowing for unused arguments).
	
	
	##################################################
	##################################################
	##### Preparation #####
	clamp<-function(x,a=NULL,b=NULL,type=c("l","q","p"),...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-05-13 - Finished.
	# Last: 2009-09-01 - Support for quantile and percentage values for 'a' and 'b' added.
	########### DESCRIPTION: ###########
	# Clamps the input 'x' to the range defined by 'a' (and 'b'), i.e. all values outside of the range are given the value of the closest range value.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# - 'x' is the input object to be clamped.
	# - 'a' is either the lower value or the range to which 'x' is to be clamped.
	# - 'b' is the upper value to which 'x' is to be clamped.
	# - 'type' (only the first element) defines the clamping method, where "l" (linear) denotes that the desired range is given by 'a' and/or 'b', "q" denotes that quantile values are given by 'a' and/or 'b', and "p" is used if percentage og the original range is given by 'a' and/or 'b'.
	# - '...' values passed on to quantile().
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(length(a)>1){
		b=a[2]
		a=a[1]
		}
	type=type[1]
	if(tolower(type)=="p"){
		r=range(x)
		diffr=diff(r)
		if(!is.null(a)){
			a=r[1]+a*diffr
			}
		if(!is.null(b)){
			b=r[1]+b*diffr
			}
		}
	else if(tolower(type)=="q"){
		if(!is.null(a)){
			a=quantile(x,a,...)
			}
		if(!is.null(b)){
			b=quantile(x,b,...)
			}
		}
	
	
	
	##### Execution and output #####
	if(is.null(a) && is.null(b)){
		x
		}
	else if(is.null(b)){
		low=x<a
		x[low]=a
		x
		}
	else if(is.null(a)){
		high=x>b
		x[high]=b
		x
		}
	else{
		low=x<a
		high=x>b
		x[low]=a
		x[high]=b
		x
		}
	}
	# Internal function calculating the colors using a function that transforms values in the color panel to colors:
	plane2col<-function(x,y){
		# Scale input from 'x' in (0,10) and 'y' in (0,10) to 'x' in (0,255*6) and 'y' in (0,255*2), to fit the colorwheel:
		x=ceiling(x*255*6/10)
		y=y*255*2/10
		# Defining the colorwheel:
		# Vectors for convenience:
		v255=rep(255,255)	
		v0=rep(0,255)
		upp=seq(0,254,length.out=255) 
		down=seq(255,1,length.out=255)
		# Colorsequences to darker:
		rseq1=c(v255,down,v0,v0,upp,v255)
		gseq1=c(upp,v255,v255,down,v0,v0)
		bseq1=c(v0,v0,upp,v255,v255,down)
		# Colorsequences to lighter:
		ll=length(rseq1)
		rseq2=c(rseq1[(ll/2+1):ll],rseq1[1:(ll/2)])
		gseq2=c(gseq1[(ll/2+1):ll],gseq1[1:(ll/2)])
		bseq2=c(bseq1[(ll/2+1):ll],bseq1[1:(ll/2)])
		# rgb-functions for 'Colorsequences to darker' and 'Colorsequences to lighter':
		colv1=function(x,y){
			rgb(rseq1[x]*y/255,gseq1[x]*y/255,bseq1[x]*y/255,maxColorValue=255)
			}
		colv2=function(x,y){
			rgb(255-rseq2[x]*y/255,255-gseq2[x]*y/255,255-bseq2[x]*y/255,maxColorValue=255)
			}
		# Creating the output colorsequence:
		colout=double()	
		for(i in 1:length(x)){
			if(y[i]<=255){
				colout[i]=colv1(x[i],y[i])
				}
			else{
				y[i]=510-y[i]
				colout[i]=colv2(x[i],y[i])
				}	
			}
		return(colout)	
		}
	
	# Drawing the plot frame:
	par(mar=double(4),oma=double(4),ann=FALSE,bty="n")
	din=par()$din
	pin=par()$pin
	plot(1,xlim=c(0,10),ylim=c(-5.5,10),axes=FALSE,col="white")
	
	# Plotting of action buttons:
	oundo=c(2,-5.5)
	oredo=c(4,-5.5)
	onext=c(6,-5.5)
	oreturn=c(8,-5.5)
	wbuttons=1.5
	hbuttons=1
	# "undo" button:
	undo.button<-function(bg=colors()[235]){
		rounded.rect(oundo,w=wbuttons,h=hbuttons,col=bg,adjust=TRUE)
		text(2,-5.5,"undo",col="darkorange",cex=sum(din)/10)
		}
	# "redo" button:
	redo.button<-function(bg=colors()[235]){
		rounded.rect(oredo,w=wbuttons,h=hbuttons,col=bg,adjust=TRUE)
		text(4,-5.5,"redo",col="darkorange",cex=sum(din)/10)
		}
	# "next" button:
	next.button<-function(bg=colors()[235],legend="next"){
		rounded.rect(onext,w=wbuttons,h=hbuttons,col=bg,adjust=TRUE)
		text(6,-5.5,legend,col="darkorange",cex=sum(din)/10)
		}
	# "return" button:
	return.button<-function(bg=colors()[235]){
		rounded.rect(oreturn,w=wbuttons,h=hbuttons,col=bg,adjust=TRUE)
		text(8,-5.5,"return",col="darkorange",cex=sum(din)/10)
		}
	# Message to print on the "next"/"prev" button
	message=c("next","prev")
	# Plotting the buttons:
	undo.button()
	redo.button()
	next.button(,message[1])
	return.button()
	
	# Vectors used when plotting the colors:
	v255=rep(255,coldens)	
	v0=rep(0,coldens)
	upp=seq(0,255,length.out=coldens) 
	down=seq(255,0,length.out=coldens)
	v255_0=rep(255,255)	
	v0_0=rep(0,255)
	upp_0=seq(0,254,length.out=255) 
	down_0=seq(255,1,length.out=255)
	
	# The color sequences:
	rseq1=c(v255,down,v0,v0,upp,v255)
	gseq1=c(upp,v255,v255,down,v0,v0)
	bseq1=c(v0,v0,upp,v255,v255,down)
	ll=length(rseq1)
	rseq2=c(rseq1[(ll/2+1):ll],rseq1[1:(ll/2)])
	gseq2=c(gseq1[(ll/2+1):ll],gseq1[1:(ll/2)])
	bseq2=c(bseq1[(ll/2+1):ll],bseq1[1:(ll/2)])
	
	# Plotting the color panel:
	for(i in upp){
		colv1=rgb(rseq1*i/255,gseq1*i/255,bseq1*i/255,maxColorValue=255)
		colv2=rgb(255-rseq2*i/255,255-gseq2*i/255,255-bseq2*i/255,maxColorValue=255)
		points(seq(0,10,length.out=ll),rep(i/255*5,ll),col=colv1,pch=".",cex=350/coldens)
		points(seq(0,10,length.out=ll),10-rep(i/255*5,ll),col=colv2,pch=".",cex=350/coldens)
		}
		
	if(identical(ref,NULL)==FALSE){
		lines(ref$colx,ref$coly,col="grey")
		}
	polygon(c(-0.3,0,0,-0.3),c(-0.3,-0.3,10.3,10.3),col="white",border=FALSE)		
	polygon(c(-0.3,10.3,10.3,-0.3),c(-0.3,-0.3,0,0),col="white",border=FALSE)		
	polygon(c(-0.3,0,0,-0.3)+10.3,c(-0.3,-0.3,10.3,10.3),col="white",border=FALSE)		
	polygon(c(-0.3,0,0,-0.3),c(-0.3,-0.3,10.3,10.3)+10.3,col="white",border=FALSE)
	lines(c(0,10),c(0,0))	
	lines(c(10,10),c(0,10))	
	lines(c(0,10),c(10,10))	
	lines(c(0,0),c(10,0))	

	## Plotting the speed panel:
	speedpanel<-function(lwd=1.5,coll="black",colt="darkblue"){
		lines(c(0,10),c(-3.7,-3.7),lwd=lwd,col=coll)
		lines(c(0,10),c(-0.3,-0.3),lwd=lwd,col=coll)
		lines(c(0,0),c(-3.7,-0.3),lwd=lwd,col=coll)
		lines(c(10,10),c(-3.7,-0.3),lwd=lwd,col=coll)
		h=-3.7+0:9*3.4/9
		for(i in 1:10){
			lines(c(-0.05,0.05),rep(h[i],2),lwd=lwd,col=coll)
			lines(c(-0.05,0.05)+10,rep(h[i],2),lwd=lwd,col=coll)
			text(-0.2,h[i],i,col=colt,cex=pin[2]/6)
			text(0.2+10,h[i],i,col=colt,cex=pin[2]/6)
			}	
		}	
	speedpanel()
	if(identical(ref,NULL)==FALSE){
		lines(ref$weightx,ref$weighty,col="grey")
		}
	
	
	##### Execution and output #####
	# Input variables and temporary input variables:
	xin=yin=xwin=ywin=thisxin=thisyin=thisxwin=thisywin=c()
	
	# Position along the input defined by the undo and redo history:
	inat=winat=0
	# 'mode' is 0 for the color panel and 1 for the speed panel:
	mode=0
	# Logical variable ending the while loop is set to TRUE:
	gotoreturn=FALSE
	
	# The user selects points in the color panel or the speed panle while gotoreturn==FALSE:
	while(!gotoreturn){
		# Updating the temporary inputvariables. These are convenient when appending new values, as they are adjusted for the undo/redo status:
		if(inat==0){
			thisxin=c()
			thisyin=c()
			}
		else{
			thisxin=xin[1:inat]
			thisyin=yin[1:inat]
			}
		if(inat==0){
			thisxwin=c()
			thisywin=c()
			}
		else{
			thisxwin=xwin[1:winat]
			thisywin=ywin[1:winat]
			}
		# Length of the input variables:
		lxin=length(xin)
		lxwin=length(xwin)
		# If the user has pressed the undo button until the start of the input variables, 'canundo' is set to FALSE. similarly, the user cannot redo something that has not been done yet (canredo=FALSE). The same holds for 'wcanundo' and 'wcanredo'
		canundo=inat>0
		canredo=inat<lxin
		wcanundo=winat>0
		wcanredo=winat<lxwin
	
		# The input mouse click:
		point=locator(1)
	
		# If the user clickes on the "return" button, the loop is terminated:
		if(contains(point$x,oreturn[1]+c(-1,1)*wbuttons/2,"b") && contains(point$y,oreturn[2]+c(-1,1)*hbuttons/2,"b")){
			return.button("grey56")
			Sys.sleep(0.15)
			return.button()
			gotoreturn=TRUE
			}
		
		# If the user clickes on the "next"/"prev" button, the mode is changed:
		else if(contains(point$x,onext[1]+c(-1,1)*wbuttons/2,"b") && contains(point$y,onext[2]+c(-1,1)*hbuttons/2,"b")){
			next.button("grey56")
			Sys.sleep(0.15)
			next.button()
			mode=1-mode
			next.button(,message[mode+1])
			}
		
		# If the user clickes on the "undo" button, 'inat' or 'winat' is reduced by 1, depending on the value of 'mode':
		else if(contains(point$x,oundo[1]+c(-1,1)*wbuttons/2,"b") && contains(point$y,oundo[2]+c(-1,1)*hbuttons/2,"b")){
			if(mode==0  && canundo){
				inat=inat-1
				undo.button("grey56")
				Sys.sleep(0.15)
				undo.button()
				# Remove the previous line:
				if(inat==0){
					xcoverline=thisxin[inat+1]
					ycoverline=thisyin[inat+1]
					linecol=plane2col(xcoverline,ycoverline)
					points(xcoverline,ycoverline,col=linecol,pch=".",cex=2)
					}
				else{
					lcoverline=sqrt((thisxin[inat+1]-thisxin[inat])^2+(thisyin[inat+1]-thisyin[inat])^2)*28
					xcoverline=seq(thisxin[inat],thisxin[inat+1],length.out=lcoverline)
					ycoverline=seq(thisyin[inat],thisyin[inat+1],length.out=lcoverline)
					linecol=plane2col(xcoverline,ycoverline)
					points(xcoverline,ycoverline,col=linecol,pch=".",cex=2)
					}
				}
			else if(mode==1 && wcanundo){
				winat=winat-1
				undo.button("grey56")
				Sys.sleep(0.15)
				undo.button()
				lines(thisxwin[winat+0:1],thisywin[winat+0:1],col="white",lwd=2.2)
				}
			}
		
		# If the user clickes on the "redo" button, 'inat' or 'winat' is increased by 1, depending on the value of 'mode':
		else if(contains(point$x,oredo[1]+c(-1,1)*wbuttons/2,"b") && contains(point$y,oredo[2]+c(-1,1)*hbuttons/2,"b")){
			# Add the next line:
			if(mode==0 && canredo){
				inat=inat+1
				redo.button("grey56")
				Sys.sleep(0.15)
				redo.button()
				lines(xin[inat-1:0],yin[inat-1:0])
				}
			else if(mode==1 && wcanredo){
				winat=winat+1
				redo.button("grey56")
				Sys.sleep(0.15)
				redo.button()
				lines(xwin[winat-1:0],ywin[winat-1:0],col="skyblue2",lwd=din[2]/4)
				}
			}
		
		# A point is appended to the input variables:
		else{
			if(mode==0){
				xin=append(thisxin,clamp(point$x,c(0,10)))
				yin=append(thisyin,clamp(point$y,c(0,10)))
				lines(xin[inat+0:1],yin[inat+0:1])
				inat=inat+1
				}
			else if(mode==1){
				# If the user mouse clicks outside of the speed panel, these points are interpreted as edge points:
				if(point$x<0){
					xwin=append(thisxwin,0)
					}
				else if(point$x>10){
					xwin=append(thisxwin,10)
					}
				else{
					xwin=append(thisxwin,point$x)
					}
					
				if(point$y<(-3.7)){
					ywin=append(thisywin,-3.7)
					}
				else if(point$y>(-0.3)){
					ywin=append(thisywin,-0.3)
					}
				else{
					ywin=append(thisywin,point$y)
					}
				lines(xwin[winat+0:1],ywin[winat+0:1],col="skyblue2",lwd=din[2]/4)
				winat=winat+1
				}
			}
		}
	
	# If the user has not selected any points, NULL is returned:
	if(is.null(xin) || inat==0){
		return(NULL)
		}
	# Else if only speed values are missing, these are set to flat:
	else if(is.null(xwin) || winat==0){
		xwin=c(0,10)
		ywin=rep(-3.7,2)
		winat=2
		}
	if(min(xwin)>0){
		xwin=c(0,xwin)
		ywin=c(-3.7,ywin)
		}
	if(max(xwin)<10){
		xwin=c(xwin,10)
		ywin=c(ywin,-3.7)
		}
	
	# Sorting the speed values with regard to the x-values, and plotting the end result:
	orden=order(xwin)
	xwin=xwin[orden]
	ywin=ywin[orden]	
	lines(xwin,ywin,col="blue",lwd=din[2]/4)
	# Calculating colors:
	xin=xin[1:inat]
	yin=yin[1:inat]
	xwin=xwin[1:winat]
	ywin=ywin[1:winat]
	w=track.line(xwin,ywin,seq(0,1,length.out=n),alongx=TRUE,normalized=TRUE)
	w=1+(w$y+3.7)*10/3.4
	cumsumw=cumsum(w)
	ww=cumsumw/max(cumsumw)
	numeric_colorspectrum=track.line(xin,yin,ww,normalized=TRUE)
	colorspectrum=plane2col(numeric_colorspectrum$x,numeric_colorspectrum$y)
	
	# Plotting the colors:
	xpdiff=10/n
	xp=t(array(c(0,xpdiff,xpdiff,0),dim=c(4,n)))+array(0:(n-1)*xpdiff,dim=c(n,4))
	yp=c(-4.8,-4.8,-3.8,-3.8)
	for(i in 1:n){
		polygon(x=xp[i,],y=yp,density=NULL,angle=45,border=FALSE,col=colorspectrum[i],lty=NULL,xpd=NULL)
		}
	
	# Return:
	list(col=colorspectrum,colx=xin,coly=yin,weightx=xwin,weighty=ywin)
	##################################################
	##################################################
	}
