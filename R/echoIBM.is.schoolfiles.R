#*********************************************
#*********************************************
#' Simple function for recognizing and returning school files from a vector of paths to files.
#'
#' @param files  is a vector of paths to files from which school files are to be extracted.
#' @param names  is a vector of four character strings representing static and dynamic school valiable names.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSD
#' @importFrom utils tail
#'
#' @export
#' @rdname echoIBM.is.schoolfiles
#'
echoIBM.is.schoolfiles<-function(files, names){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-02-20 - Clean version
	########### DESCRIPTION: ###########
	# Simple function for recognizing and returning school files from a vector of paths to files.
	########## DEPENDENCIES: ###########
	# read.TSD()
	############ VARIABLES: ############
	# ---files--- is a vector of paths to files from which school files are to be extracted.
	# ---names--- is a vector of four character strings representing static and dynamic school valiable names.
	
	
	##################################################
	##################################################
	########## Preparation and execution ##########
	# For loop through the school files:
	schoolfiles=logical(length(files))
	for(i in seq_along(files)){
		ext=tail(unlist(strsplit(files[i],".",fixed=TRUE)),1)
		suppressWarnings(varnames<-read.TSD(files[i],var="",header=TRUE)$labl)
		if(identical(ext,"school") || any(names %in% varnames)){
			schoolfiles[i]=TRUE
			}
		}
	
	
	########## Output ##########
	schoolfiles
	##################################################
	##################################################
	}
