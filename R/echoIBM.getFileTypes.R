#*********************************************
#*********************************************
#' Reads the header of all files in an event and classifies into the following file types: beams (beam configuration), ctd, vessel, pings, school, tsd, other.
#'
#' @param files		list of files for which file type should be detemined.
#' @param ext		Vector of file extensions to test the files against.
#' @param labl		Key varialbe names corresponding to the file types given in \code{ext}.
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
echoIBM.getFileTypes <- function(
	files, 
	ext=list(beams="beams", ctd="ctd", pings="pings", school="school", vessel="vessel", noise=NULL), 
	labl=list(beams=c("cali"), noise=c("bgns", "nrns", "nr0a", "nr0p")), 
	recursive=TRUE)
	{
	
	matchExtLabl <- function(x, filesext, filelabl, ext, labl){
		which( filesext %in% ext[[x]] | (unlist(lapply(filelabl, function(y) any(y %in% labl[[x]])))) )
	}

	# Get the files, if 'files' is a directory:
	if(length(files)==1 && isTRUE(file.info(files)$isdir)){
		files <- list.files(files, recursive=recursive, full.names=TRUE)
	}
	# Read variable names and file extensions:
	suppressWarnings(filelabl <- lapply(files, function(x) read.TSD(x, header=TRUE, var=NULL)$labl))
	filesext <- tools::file_ext(files)
	
	# Get the indices at each file extension:
	if(!is.list(ext)){
		namesext <- ext
		ext <- as.list(ext)
		names(ext) <- namesext
	}
	
	categories <- unique(c(names(ext), names(labl)))
	ind <- lapply(categories, matchExtLabl, filesext=filesext, filelabl=filelabl, ext=ext, labl=labl)
	names(ind) <- categories
	
	# And the files with extension that are not included in 'ext':
	ind <- c(ind, list(other=which(!filesext %in% unlist(ext))))
	# Return the files:
	lapply(ind, function(x) files[x])
}
