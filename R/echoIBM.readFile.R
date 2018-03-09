#*********************************************
#*********************************************
#' Returns a list of the following strings: (1) the path to the event, (2) the event name, (3) the event number, (4) the path to the cruise, and (5) the cruise name.
#'
#' @param path		A list of the following elements: (1) 'path', giving the paths to the sub-events, (2) 'esnm', giving the names of the acoustic instruments in the events (same length as 'path'), and (3) 'name', giving the name of the event.
#' @param ext		A string naming the file extension of the file to read.
#' @param t			The time steps to read.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom tools file_ext
#'
#' @export
#' @rdname echoIBM.setup
#'
echoIBM.readFile <- function(path, ext="vessel", t=1){
	file <- list.files(path, full.names=TRUE)
	file <- file[tools::file_ext(file)==ext]
	#if(length(pattern)){
	#	file <- file[grep(pattern, file)]
	#}
	#file <- head(file, 1)
	read.TSDs(file, t=t)
}