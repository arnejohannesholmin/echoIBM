#*********************************************
#*********************************************
#' Draws a rounded rectangle.
#'
#' @param origin  is the origin of the rounded rectangle.
#' @param w  is the width of the rounded rectangle.
#' @param h  is the height of the rounded rectangle.
#' @param r  is the radius of the rounded corners as a fraction on 'w'.
#' @param n  is the number of points along the rounded corners of the rounded rectangle.
#' @param plot  is whether the rounded rectangle is to be plotted or not.
#' @param adjust  is TRUE if the rounded corners should appear as circles regardless of the frame dimensions.
#' @param by  is along which of width, height or both, the adjustment for reshaped window dimensions of the rounded rectangles are to be done.
#' @param ...  are parameters like 'col' and 'border' used in polygon().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom sonR circle
#' @importFrom graphics par polygon
#'
#' @export
#' @rdname rounded.rect
#'
rounded.rect<-function(origin=0,w=1,h=1,r=0.1,n=100,plot=TRUE,adjust=FALSE,by="hw",...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-05-08 - Finished.
	########### DESCRIPTION: ###########
	# Draws a rounded rectangle.
	########## DEPENDENCIES: ###########
	# circle()
	############ VARIABLES: ############
	# ---origin--- is the origin of the rounded rectangle.
	# ---w--- is the width of the rounded rectangle.
	# ---h--- is the height of the rounded rectangle.
	# ---r--- is the radius of the rounded corners as a fraction on 'w'.
	# ---n--- is the number of points along the rounded corners of the rounded rectangle.
	# ---plot--- is whether the rounded rectangle is to be plotted or not.
	# ---adjust--- is TRUE if the rounded corners should appear as circles regardless of the frame dimensions.
	# ---by--- is along which of width, height or both, the adjustment for reshaped window dimensions of the rounded rectangles are to be done.
	# ---...--- are parameters like 'col' and 'border' used in polygon().
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Default 'r':
	r_default=0.1
	# Origin need to have length 2:
	origin=rep(origin,length.out=2)
	# 'size' is used to calibrating the radius 'r' of the rounded corners to the widht, the height or both width and heigth of the rectangle. 'size' is a vector holding 'w' and/or 'h':
	size=NULL
	if(!identical(grep("h",by),integer(0))){
		size=c(size,h)
		}
	if(!identical(grep("w",by),integer(0))){
		size=c(size,w)
		}
	
	# If adjust==TRUE, the rounded conrners is stretched according to the frame dimensions and the inner coordinate dimensions:
	adjustx=adjusty=1
	if(adjust){
		displayh=(par()$din[1]-sum(par()$mai[c(2,4)]))*diff(par()$usr[3:4])
		displayw=(par()$din[2]-sum(par()$mai[c(1,3)]))*diff(par()$usr[1:2])
		stretch=displayh/displayw
		if(stretch>1){
			# adjustx=1/stretch
			adjustx=1/sqrt(stretch)
			adjusty=sqrt(stretch)
			}
		else{
			# adjusty=1/stretch
			adjusty=1/sqrt(stretch)
			adjustx=sqrt(stretch)
			}
		# Calibrate according to the input 'by'. If by has length
		}
	
	if(r*mean(size)/min(h,w)*min(adjustx,adjusty)>0.5){
		warning("Corneres to widely rounded or illegal value of 'r'")
		}
	r=r*mean(size)*min(adjustx,adjusty)
	# Adjusting inputs:
	w=w/2
	h=h/2
	n=ceiling(n/4)
	
	
	##### Execution #####
	# Origins for the four quarter circles (upperleft, lowerleft, lowerright, upperright):
	o1=origin+c(-w+r,h-r)
	o2=origin+c(-w+r,-h+r)
	o3=origin+c(w-r,-h+r)
	o4=origin+c(w-r,h-r)
	
	# The four circle segments:
	c1=circle(o1,r,seq(0,pi/2,length.out=n)+pi/2)
	c2=circle(o2,r,seq(0,pi/2,length.out=n)+pi)
	c3=circle(o3,r,seq(0,pi/2,length.out=n)+3*pi/2)
	c4=circle(o4,r,seq(0,pi/2,length.out=n)+2*pi)
	
	# 'upperx', 'lowerx', 'uppery' and 'lowery' are the upper and the lower x- and y-coordinates of the rectangle:
	upperx=origin[1]+w
	lowerx=origin[1]-w
	uppery=origin[2]+h
	lowery=origin[2]-h
	c1=cbind(lowerx+(c1[,1]-lowerx)*adjustx,uppery+(c1[,2]-uppery)*adjusty)
	c2=cbind(lowerx+(c2[,1]-lowerx)*adjustx,lowery+(c2[,2]-lowery)*adjusty)
	c3=cbind(upperx+(c3[,1]-upperx)*adjustx,lowery+(c3[,2]-lowery)*adjusty)
	c4=cbind(upperx+(c4[,1]-upperx)*adjustx,uppery+(c4[,2]-uppery)*adjusty)
	
	# Constructing the output matrix:
	p=rbind(c1,c2,c3,c4)
	colnames(p)=c("x","y")
	
	
	##### Return #####
	if(plot){
		polygon(p,...)
		}
	p
	##################################################
	##################################################
	}
