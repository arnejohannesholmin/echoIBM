#*********************************************
#*********************************************
#' Merges empirical beam patterns located in separate TSD-files. If duplicated 'grsf' values are present the first is allways chosen.
#'
#' @param files  is a vector of strings representing the paths to the files to merge, or a single string representing the directory in which the files are located.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSDs write.TSD
#'
#' @export
#'
merge_ebpf<-function(files,con=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-06-26 - Clean version.
	########### DESCRIPTION: ###########
	# Merges empirical beam patterns located in separate TSD-files. If duplicated 'grsf' values are present the first is allways chosen.
	########## DEPENDENCIES: ###########
	# read.TSDs()
	############ VARIABLES: ############
	# ---files--- is a vector of strings representing the paths to the files to merge, or a single string representing the directory in which the files are located.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Get the relevant data:
	data=read.TSDs(files,var=c("graf","gref","grsf","ebpf","obln","pbpf"),clean=FALSE)
	# Only merge if more than one beam patterns is present
	if(sum(names(data)=="ebpf")>1 && sum(names(data)=="grsf")>1){
		out=list()
		# Check if all grid azimuth angles 'graf' are of equal length, and discard this angle if not:
		if(any(names(data)=="graf")){
			graf=data[names(data)=="graf"]
			if(all(sapply(graf,length)==length(graf[[1]]))){
				out$graf=graf[[1]]
				}
			}
		# Check if all grid elevation angles 'gref' are of equal length, and discard this angle if not:
		if(any(names(data)=="gref")){
			gref=data[names(data)=="gref"]
			if(all(sapply(gref,length)==length(gref[[1]]))){
				out$gref=gref[[1]]
				}
			}
		
		
		##### Execution and output #####
		# If both grid relative size 'grsf' and empirical beam pattern 'ebpf' is present, merge these by discarding the duplicated values of 'grsf' and the corresponding entries of 'ebpf':
		if(any(names(data)=="grsf") && any(names(data)=="ebpf")){
			grsf=unlist(data[names(data)=="grsf"])
			ebpf=unlist(data[names(data)=="ebpf"])
			# Set the dimension of 'ebpf':
			newdim=c(length(out$graf),length(out$gref),length(grsf))
			newdim=newdim[newdim!=0]
			dim(ebpf)=newdim
			# Remove the duplicated entries:
			dup_grsf=duplicated(grsf)
			out$grsf=grsf[!dup_grsf]
			if(length(newdim)==2){
				out$grsf=grsf[!dup_grsf]
				out$ebpf=ebpf[,!dup_grsf]
				# Order the data:
				ordergrsf=order(out$grsf)
				out$grsf=out$grsf[ordergrsf]
				out$ebpf=out$ebpf[,ordergrsf]
				}
			else if(length(newdim)==3){
				out$grsf=grsf[!dup_grsf]
				out$ebpf=ebpf[,,!dup_grsf]
				# Order the data:
				ordergrsf=order(out$grsf)
				out$grsf=out$grsf[ordergrsf]
				out$ebpf=out$ebpf[,ordergrsf]
				}
			else if(length(newdim)==1){
				warning("Grid elevation angle corresponding to the empirical beam pattern of the targets is missing")
				}
			}
		if(!any(length(con)==0,identical(con,FALSE))){
			out$dbpf=dim(out$ebpf)
			if(!is.null(data$obln)){
				out$obln=data$obln
				}
			if(!is.null(data$pbpf)){
				out$pbpf=data$pbpf
				}
			
			write.TSD(out,con,ts=0)
			}
		out
		}
	else{
		data
		}
	##################################################
	##################################################
	}
