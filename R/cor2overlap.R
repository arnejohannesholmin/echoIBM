#*********************************************
#*********************************************
#' Transforms from correlation values to overlap values used in rexp_MultSines().
#'
#' @param cor  is a vector or matrix of correlation values.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all
#'
#' @export
#' @rdname cor2overlap
#'
cor2overlap<-function(cor){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-09-25 - Clean version.
	########### DESCRIPTION: ###########
	# Transforms from correlation values to overlap values used in rexp_MultSines().
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---cor--- is a vector or matrix of correlation values.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	generateOverlapMatrix=function(overlap){
		N=ceiling(max(overlap))
		M=outer(seq(-N+1,0),overlap,"+")
		M[M<0]=0
		M[M>1]=1
		if(ncol(M)==1){
			rbind(M,1,matrix(M[rev(seq_len(N)),],length(M),1))
			}
		else{
			rbind(M,1,M[rev(seq_len(N)),])
			}
		}
	# Save the old dimension and transform to no dimension:
	olddim=dim_all(cor)
	dim(cor)=NULL
	
	
	########## Execution ##########
	# Calculate the overlap values:
	overlap=findInterval(cor,overlaptable[,1],all.inside=TRUE)
	overlap=overlaptable[overlap,2]
	overlap=generateOverlapMatrix(overlap)
	
	
	########## Output ##########
	# Reset dimensions:
	dim(overlap)=c(dim(overlap)[1],olddim)
	overlap
	##################################################
	##################################################
	}
