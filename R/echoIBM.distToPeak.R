#*********************************************
#*********************************************
#' Calculates the distance to the nearest peak of the periodic noise for time step 'indt'.
#'
#' @param data  is a list containing beam configuration, and periodic noise parameters. Speficically the following variables must be included: Beams variables: 'sint', 'lenb'; periodic noise variables: 'pns1', 'pns2', 'pns3', 'acfq', 'harm'; and either 'bgns' or 'numb'
#' @param indt  is the time step index.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD zeros
#'
#' @export
#' @rdname echoIBM.distToPeak
#'
echoIBM.distToPeak<-function(data,indt=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-06-26 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the distance to the nearest peak of the periodic noise for time step 'indt'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing beam configuration, and periodic noise parameters. Speficically the following variables must be included: Beams variables: 'sint', 'lenb'; periodic noise variables: 'pns1', 'pns2', 'pns3', 'acfq', 'harm'; and either 'bgns' or 'numb'
	# ---indt--- is the time step index.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(length(dim(data$bgns))==2){
		data$dist=zeros(dim(data$bgns))
		}
	else if(length(data$lenb)>0 && length(data$numb)>0){
		data$dist=zeros(max(data$lenb),data$numb)
		}
	
	
	########## Execution ##########
	if(!all(c("sint","lenb","numb","acfq","badb","harm","pns3") %in% names(data))){
		warning("Periodic noise not subtracted. Specificaiton of periodic noise not present in the data")
		return(data$dist)
		}
	periodic=which(data$badb==1)
	# Define the frequency of the periodic noise, and the sequence along beams:
	a=2*pi * data$acfq*data$sint * data$harm[periodic]
	r=seq_len(max(data$lenb))
	# Calculate the periodic noise:
	if(length(indt)>0){
		if(indt[1]<1 && indt[1]>NCOL(data$pn3M)){
			warning(paste("'indt' (",indt[1],") chosed outside of the valid range [1,",NCOL(data$pn3M),"]",sep=""))
			}
		data$pns3=data$pn3M[,indt[1]]
		}
	
	# Define the period of the sine wave:
	period=1/(data$acfq*data$sint * data$harm[periodic])
	# Calculate the distance from a voxel to the nearest peak to the right:
	for(i in seq_along(periodic)){
		data$dist[,periodic[i]] = round(((pi/2-data$pns3[periodic[i]])/a[i] - r) %% period[i])
		}
			
	
	########## Output ##########
	data$dist
	##################################################
	##################################################
	}
