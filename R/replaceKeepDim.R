#*********************************************
#*********************************************
#' Replaces values in a list but keeping dimension
#'
#' @param x				The list to replace values in.
#' @param replacement	A list of the values to replace.
#' @param esnm			The name of the sonar, used when the replacement is a list with \code{esnm} as names.
#' @param ...			Arguments passed on to replacement functions.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname echoIBM.setup
#'
replaceKeepDim <- function(x, replacement, esnm, ...){
	# Get the variables to reaplace:
	inreplacement <- intersect(names(replacement), names(x))
	# Replace only if 'inreplacement' has length>0:
	if(length(inreplacement)){
		# Move through the replacements:
		for(i in seq_along(inreplacement)){
			# If the replacement is a list, assume there is one element per sonar:
			this <- replacement[[inreplacement[i]]]
			if(is.list(this)){
				this <- this[[esnm]]
			}
			if(length(this)){
				# Get the dimension of the original variable:
				thisdim <- dim_all(x[[inreplacement[i]]])
				# If a funciton of the dimension, run this, with additional parameters given in the '...':
				if(is.function(this)){
					x[[inreplacement[i]]]  <- this(thisdim, ...)
				}
				else{
					# Otherwese repeat to a vector or an array:
					if(length(thisdim)<2){
						x[[inreplacement[i]]] <- rep(this, length.out=thisdim)
					}
					else{
						x[[inreplacement[i]]] <- array(this, dim=thisdim)
					}
				}
			}
		}
	}
	return(x)
}
