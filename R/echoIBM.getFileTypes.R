#*********************************************
#*********************************************
#' Reads the header of all files in an event and classifies into the following file types: beams (beam configuration), ctd, vessel, pings, school, tsd, other.
#'
#' @param files		list of files for which file type should be detemined.
#' @param ext		Vector of file extensions to test the files against.
#' @param key		Key varialbe names corresponding to the file types given in \code{ext}.
#' @param recursive	Used in list.files() when files is a single directory.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom tools file_ext
#'
#' @export
#'
echoIBM.getFileTypes <- function(files, ext=c("beams", "ctd", "pings", "school", "tsd", "vessel"), key=list(), recursive=TRUE){
	# Get the files, if 'files' is a directory:
	if(length(files)==1 && isTRUE(file.info(files)$isdir)){
		files <- list.files(files, recursive=recursive, full.names=TRUE)
	}
	# Read variable names and file extensions:
	suppressWarnings(labl <- lapply(files, function(x) read.TSD(x, header=TRUE, var=NULL)$labl))
	filesext <- tools::file_ext(files)
	
	# Get the indices at each file extension:
	ind <- lapply(ext, function(x) which(filesext==x | unlist(lapply(labl, function(y) any(y %in% key[[x]])))))
	names(ind) <- ext
	
	# And the files with extension that are not included in 'ext':
	ind <- c(ind, list(other=which(!ext %in% ext)))
	# Return the files:
	lapply(ind, function(x) files[x])
}
