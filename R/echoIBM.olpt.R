#*********************************************
#*********************************************
#' Calculates the one step (to the adjacet voxel of the same sampling interval) overlap matrix given by the background noise and the periodic noise for the given time step index. The output is a matrix with the overlap values alont the columns, and one column for each voxel of each beam.
#'
#' @param indt  is the time index used when constructing the periodic noise.
#' @param adds  is a list of variables required in the function but not read from the data stored in noise.path(). This is used primarily when creating the file of overlap values, in which case 'crb1' and 'crp1' are needed.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.olpt
#'
echoIBM.olpt<-function(data,indt=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-09-25 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the one step (to the adjacet voxel of the same sampling interval) overlap matrix given by the background noise and the periodic noise for the given time step index. The output is a matrix with the overlap values alont the columns, and one column for each voxel of each beam.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---indt--- is the time index used when constructing the periodic noise.
	# ---adds--- is a list of variables required in the function but not read from the data stored in noise.path(). This is used primarily when creating the file of overlap values, in which case 'crb1' and 'crp1' are needed.
	
	
	##################################################
	##################################################
	##### Preparation and execution #####
	dimb=c(data$numb/length(unique(data$freq)),length(unique(data$freq)))

	# Generate the periodic noise and the background noise:
	dist=echoIBM.distToPeak(data,indt=indt)
	
	# Define the correlation matrices for the background noise and the periodic noise, by expanding from one value for each beam to one value for each voxel:
	if(length(data$crb1)==0){
		stop("Correlation data 'crb1' missing in the simulation files/input")
		}
	### REMOVED 2013-09-26 ### crb1M=matrix(data$crb1,nrow=max(data$lenb),ncol=data$numb,byrow=TRUE)
	crtV=matrix(data$crb1,nrow=max(data$lenb),ncol=data$numb,byrow=TRUE)
	if(length(data$crp1)==0){
		stop("Correlation data 'crp1' missing in the simulation files/input")
		}
	### REMOVED 2013-09-26 ### crp1M=matrix(data$crp1,nrow=max(data$lenb),ncol=data$numb,byrow=TRUE)
	
	# Calculate the correlation as the correlation of the background noise, added the difference between the correlation of the periodic noise and the background noise times the size of the periodic noise compared to the backgroun noise:
	### REMOVED 2013-09-26 ### crtV = crb1M
	for(i in which(data$badb==1)){
		crtV[,i] = data$cp1j[i,dist[,i]+1]
		}
	
	# Apply correlation also for the last beam of each fan, and reshape to a four dimensional array at the end, so that one fan is treated at the time, and no correlation is exists between fans:
	for(i2 in seq_len(dimb[2])){
		crtV[,dimb[1]+(i2-1)*dimb[1]] = crtV[,dimb[1]+(i2-1)*dimb[1]-1] 
		} 
	
	# Calculate the overlap values for the correlation matrix of the total correlation, variable for each voxel:
	olpt=cor2overlap(crtV)
	
	
	##### Output #####
	# Return a list of the correlation matrix for all voxels 'crtV' and the overlap array 'olpt'
	list(olpt=olpt,crtV=crtV)
	##################################################
	##################################################
	}
