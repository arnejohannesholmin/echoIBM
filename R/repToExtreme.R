#*********************************************
#*********************************************
#' Repeats all elements of the list \code{x} to the 
#'
#' @param x			A list of arrays.
#' @param ndim		The number of dimensions of the arrays to repeat to the maximum dimension.
#' @param skip.len	The length of arrays to leave be.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname repToExtreme
#'
repToExtreme <- function(x, ndim=2, skip.len=1){
	# Get the extreme dimension:
	d <- dim_all(x)
	# Treat only arrays of 'ndim' dimensions:
	valid <- which(sapply(d, length) == ndim)
	if(length(valid)){
		# Get a matrix of the dimensions of the arrays:
		d_valid <- do.call("rbind", d[valid])
		# And a vector of the lengths:
		l_valid <- apply(d_valid, 1, prod)
		maxd <- apply(d_valid, 2, max)
		# Skip lengths unequal to 'skip.len':
		valid <- valid[l_valid != skip.len]
		# Repeat valid arrays to the maximum dimension:
		for(i in valid){
			x[[i]] <- array(x[[i]], dim=maxd)
		}
	}
	x	
}

