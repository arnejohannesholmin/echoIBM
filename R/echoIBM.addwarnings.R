#*********************************************
#*********************************************
#' Writes warning messages to the end of a echoIBM-dumpfile.
#'
#' @param event  is either the path to a directory containing the dumpfile or the path to the dumpfile itself.
#' @param warn  is a vector of string with warning messages as returned from catchWarnings() or from warnings2character() or generated in some other way.
#' @param pre  is a string to be written on the line preceding the warnings.
#' @param post  is a string to be written on the line following the warnings.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.addwarnings
#'
echoIBM.addwarnings<-function(event,warn=NULL,pre="\n\n",post=""){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-08-12 - Clean version.
	########### DESCRIPTION: ###########
	# Writes warning messages to the end of a echoIBM-dumpfile.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---event--- is either the path to a directory containing the dumpfile or the path to the dumpfile itself.
	# ---warn--- is a vector of string with warning messages as returned from catchWarnings() or from warnings2character() or generated in some other way.
	# ---pre--- is a string to be written on the line preceding the warnings.
	# ---post--- is a string to be written on the line following the warnings.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Allow for a vector of strings:
	event=event[1]
	# Extract the dumpfile if 'event' is a directory:
	if(file.info(event)$isdir){
		files=list.files(event,full.names=TRUE)
		isdumpfile=logical(length(files))
		for(i in seq_along(isdumpfile)){
			suppressWarnings(line1 <- readLines(files[i],n=1))
			ncharline1=nchar(line1,allowNA=TRUE)
			if(length(ncharline1)>0 && !is.na(ncharline1)){
				isdumpfile[i]=length(grep("dump from the simulation experiment built from the files located in the directory", line1, ignore.case=TRUE))>0
				}
			}
		if(sum(isdumpfile)==0){
			stop(paste("No dumpfile starting with the line \"Dump from the simulation experiment built from the files located in the directory (with subdirectories)\" found in the directory ",event,sep=""))
			}
		else if(sum(isdumpfile)>1){
			# Get the newest dumpfile:
			mtime=unclass(file.info(files[isdumpfile])$mtime)
			event=files[isdumpfile][which.max(mtime)]
			warning(paste("More than one dumpfile found. The file ",event," chosen",sep=""))
			}
		else{
			event=files[isdumpfile][1]
			}
		}
	# Else verify the dumpfile:
	else{
		line1=readLines(event,n=1)
		isdumpfile=line1 %in% c("########## DUMP FROM THE SIMULATION EXPERIMENT BUILT FROM THE FILES LOCATED IN THE DIRECTORY (WITH SUBDIRECTORIES): ##########","########## DUMP FROM MERGING AND ADDING NOISE TO SIMULATION FILES LOCATED IN THE DIRECTORY (WITH SUBDIRECTORIES): ##########")
		if(!isdumpfile){
			stop(paste("The file ",event," does not correspond to the fist line assumed for dumpfiles (\"Dump from the simulation experiment built from the files located in the directory (with subdirectories)\")",sep=""))
			}
		}
	
	
	########## Execution and output ##########
	if(length(warn)==0){
		warn=warnings2character(header=FALSE,numbering=FALSE)
		}
	if(is.list(warn)){
		thesecalls=unlist(lapply(warn$warnings,function(x) x$call))
		thesewarn=unlist(lapply(warn$warnings,function(x) x$message))
		warn=paste(paste("In ",thesecalls,sep=""),paste(":\n  ",thesewarn,sep=""),sep="")
		}
	# Remove duplicated warnings:
	warn=unique(warn)
	if(length(warn)>0){
		write(c(pre,"##### OTHER WARNINGS: #####",paste(as.character(seq_along(warn)),": ",warn,sep=""),post),file=event,append=TRUE)
		}
	else{
		write("##### OTHER WARNINGS: #####\nnone",file=event,append=TRUE)
		}	
	##################################################
	##################################################
	}
