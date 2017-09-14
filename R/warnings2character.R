#*********************************************
#*********************************************
#' Returns the warning messages stored in 'x' which is an object of type "warnings" as retunred from warnings().
#'
#' @param x  is the warnings-object to translate to character string vector.
#' @param numbering  is TRUE if numbering should be adde to the warnings, as is done when warnings() is used.
#' @param header  is TRUE the first string of the output is "Warning message" or "Warning messages".
#' @param width.cutoff  is an integer in [20, 500] determining the cutoff at which line-breaking is tried.
#' @param nlines  is an integer giving the maximum number of lines to produce. Negative values indicate no limit.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname warnings2character
#'
warnings2character<-function(x=warnings(),numbering=TRUE,header=TRUE,width.cutoff=50L,nlines=2L){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-08-04 - First version.
	########### DESCRIPTION: ###########
	# Returns the warning messages stored in 'x' which is an object of type "warnings" as retunred from warnings().
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---x--- is the warnings-object to translate to character string vector.
	# ---numbering--- is TRUE if numbering should be adde to the warnings, as is done when warnings() is used.
	# ---header--- is TRUE the first string of the output is "Warning message" or "Warning messages".
	# ---width.cutoff--- is an integer in [20, 500] determining the cutoff at which line-breaking is tried.
	# ---nlines--- is an integer giving the maximum number of lines to produce. Negative values indicate no limit.
		
	
	##################################################
	##################################################
	if (n <- length(x)) {
		# Define the character variable to return and initiate the heading string:
		warn=character(n)
		
		msgs=names(x)
		# Run through the warnings and transform to strings as displayed using print.warnings():
		for(i in seq_len(n)) {
			# Add numbering if required:
			if(n>1L && numbering){
				ind=paste(i, ": ", sep = "")
				}
			else{
				ind=""
				}
			# Transform to string and insert to the output:
			if (length(x[[i]])){
				temp <- deparse(x[[i]], width.cutoff = width.cutoff, nlines = nlines)
				sm <- strsplit(msgs[i], "\n")[[1L]]
				nl <- if (nchar(ind, "w") + nchar(temp[1L], "w") + nchar(sm[1L], "w") <= 75L) 
				  " "
				else "\n  "
				warn[i]=paste(ind, "In ", temp[1L], if (length(temp) > 1L) {" ..."}, " :", nl, msgs[i], sep = "")
				}
			else{
				warn[i]=paste(ind, msgs[i], sep = "")
				}
			}
		if(header){
			c(ngettext(n, "Warning message:\n", "Warning messages:\n"),warn)
			}
		else{
			warn
			}
		}
	##################################################
	##################################################
	}
