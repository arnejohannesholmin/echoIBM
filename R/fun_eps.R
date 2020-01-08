#*********************************************
#*********************************************
#' Returns a function for the echo ability 'eps*' equal to the optimal acoustic cross sectional area of a target divided by the length of the target to the desired power (usually squared). The function fun_eps() has 2 methods depending on the type of the input object 'data'.
#'
#' @param data  may be of 2 different types: (1) a function of 1 argument representing the acoustic cross sectional area of the target as a function of frequency, and (2) a list representing the empirical acoustic cross sectional area of a source in intervals or at points of the arguments. Names for the elements of the list adopted from read.TSD(): [['grff']] - A vector of arbitrary length representing the grid vector of acoustic frequency, and [['eps*']] - A vector no shorter than length(grff) holding the echo ability 'eps*' corresponding to the frequency grid vector 'grff'.
#' @param Nl  is the number of targets from which the echo is simulated.
#' @param method  defines the interpolation method.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom stats approxfun
#'
#' @export
#' @rdname fun_eps
#'
fun_eps<-function(data=NULL,Nl=1,method=c("closest","linear","constant")){
	
	############### LOG: ###############
	# Start: 2009-10-19 - Clean version, adapted from beamPattern().
	# Last: 2011-08-06 - Changed to return NULL if data==NULL and a repreate 'data' to the length 'Nl' if 'data' is a numeric vector. Also changed the name to fun_eps, and the function no accepts all variables starting with 'eps' in the list 'data' as the echo ability.
		

	##################################################
	##################################################
	##### Preparation, execution #####
	# (1) If 'data' is given as a function, this funciton is returned with the requirement that the number of arguments to the function is 3:
	if(is.function(data)){
		fun=data
		}	
		
	# (2) If 'data' is is given as a list containing a grid frequency vector 'grff' and the corresponding 'eps' vector:
	else if(is.list(data)){
		
		if(length(data$eps)==0){
			fun=NULL
			}
		# If both 'grff' and 'eps' are present in 'data', an empirical table of values is assumed:
		else if(!is.null(data$grff) && !is.null(data$eps)){
			# If the vetors have equal lengths:
			if(length(data$grff)==length(data$eps)){
				# If the length of the grid variable equals and method[1]="closest", the corrsponding dimension of the empirical beampattern, the closest grid point is selected:
				if(method[1]=="closest"){
					fun=function(x){
						x=findInterval(x,c(-Inf, data$grff[-length(data$grff)] + diff(data$grff)/2, Inf),rightmost.closed=TRUE,all.inside=TRUE)
						# Extracting the values:
						data$eps[x]
						}
					}
				else{
					fun=approxfun(data$grff, data$eps, method=method[1], rule = 2)
					}
				}
			# Else if the response variable is related to intervals of the explanatory variable, linear interpolation is done by the function approxfun() on the points in the middle of the intervals:
			else if(length(data$grff)==length(data$eps)+1){
				if(method[1]=="linear"){
					between=data$grff[-length(data$grff)] + diff(data$grff)/2
					fun=approxfun(between, data$eps, method=method[1], rule = 2)
					}
				else{
					fun=function(x){
						x=findInterval(x,data$grff,all.inside=TRUE,rightmost.closed=TRUE)
						# Extracting the values:
						data$eps[x]
						}
					}
				}
			# Else return the numeric vector data$eps:
			else{
				fun=rep(data$eps,length.out=Nl)
				}
			}
		# Return the function if data$eps is a function:
		else if(is.function(data$eps)){
			fun=data$eps
			}
		# Return the numeric vector if data$eps is a numeric vector:
		else if(is.numeric(data$eps)){
			fun=rep(data$eps,length.out=Nl)
			}
		else{
			stop("Invalid value of data$eps (must be a numeric table along with data$grff, a function of one parameter, or a numeric vector)")
			}
		}
	# 'data' has length 0, return NULL:
	else if(length(data)==0){
		fun=NULL
		}	
	# If 'data' is numeric, return 'data':
	else if(is.numeric(data)){
		fun=rep(data,length.out=Nl)
		}
	else{
		stop("Invalid value of 'data' (must be a function of one parameter or a numeric vector)")
		}
		
		
	##### Output #####
	fun
	##################################################
	##################################################
	}
