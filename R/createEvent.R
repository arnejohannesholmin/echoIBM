#*********************************************
#*********************************************
#' Creates event skeleton.
#'
#' @param event	The name of the event.
#' @param dir	The directory in which to put the event.
#' @param esnm	The names of the acoustic instruments, such as EK60, SU90, ME70 and MS70.
#'
#' @return Paths and name of the events.
#'
#' @export
#' @rdname createEvent
#' 
createEvent <- function(event, dir, esnm, ow=NULL){
	if(!all(basename(event)=="tsd")){
		event <- file.path(dir, event, esnm, "tsd")
	}
	name <- rev(pathparts(event[1]))[3]
	
	# Create the event skeletons:
	event.create <- function(event, ow){
		if(file.exists(event) && isTRUE(file.info(event)$isdir)){
			if(length(ow)==0){
				ans <- readline(paste0("Event ", event, " already existis. Overwrite? (y/n)"))
				if(ans=="y"){
					unlink(dirname(event), recursive=TRUE, force=TRUE)
					cat("Overwriting\n")
					dir.create(event, recursive=TRUE)
				}
				else{
					cat("Not overwriting\n")
				}
			}
			else if(isTRUE(ow)){
				unlink(dirname(event), recursive=TRUE, force=TRUE)
				cat("Overwriting\n")
				dir.create(event, recursive=TRUE)
			}
			else{
				cat("Not overwriting\n")
			}
		}
		else{
			dir.create(event, recursive=TRUE)
		}
	}
	
	lapply(event, event.create, ow=ow)
	# The returned list element 'event' may have length>1, whereas length(eventName) is always 1:
	return(list(path=event, esnm=esnm, name=name))
}