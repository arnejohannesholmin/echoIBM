#*********************************************
#*********************************************
#' Set the directory of the acoustic data used by the sonR package
#'
#' @param path	The path to the directory in which to put simulated data, structured as CruiseName/"Events"/EventName/EchosounderName/"raw" for the raw files and substituting "raw" by "tsd" for TSD files.
#'
#' @export
#'
echoIBM_datasets_directory <- function(path=NULL){
	file = "echoIBM_datasets_directory.txt"
	extdata = file.path(find.package("echoIBM"), "extdata")
	file = file.path(extdata, file)
	if(length(path)>0 && nchar(path)>0){
		return(writeLines(path, file))
	}
	extdatafiles = list.files(extdata, full.names=TRUE)

	if(file %in% extdatafiles){
		readLines(file, warn=FALSE)
		}
	else{
		stop(paste0("The file \n\n", file, "\n\nmissing. Use \n\nechoIBM_datasets_directory(PATH_TO_DIRECTORY_HOLDING_ECHOIBM_PROJECTS_AND_RESOURCE_FILES) \n\nto set the path to the directory holding the echoIBM projects and resource files (such as \"~/Data/echoIBM\")"))
		}
	}
