#*********************************************
#*********************************************
#' Calculates the cummulative probability function of the total intensity (square of pressure amplitude) from a finite number of sine waves of arbitrary individual amplitudes, given (for the upper probability) in Equation 65 of Barakat 1974, using the terminology of that paper. NOTE: The expression in Equation 65 in Barakat 1974 an error. The argument h_0 should be sqrt(h_0).
#'
#' @param n  is the number of observations to draw.
#' @param beta_n2  is the vector of individual intensities, or a single number giving the number of significant scatterers.
#' @param N  is the number of positive roots of the Bessel function of the first kind.
#' @param max.memory  is the maximum memory occupied by the function before splitting into a loop for each value of 'r'.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD NAs
#' @importFrom stats optimize rnorm
#'
#' @export
#' @rdname rBarakatI
#'
rBarakatI=function(n,beta_n2=NULL,magn=1,nsig=NULL,N=100,betadistr=c("seq","flat"),max.memory=1e9){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-27 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the cummulative probability function of the total intensity (square of pressure amplitude) from a finite number of sine waves of arbitrary individual amplitudes, given (for the upper probability) in Equation 65 of Barakat 1974, using the terminology of that paper. NOTE: The expression in Equation 65 in Barakat 1974 an error. The argument h_0 should be sqrt(h_0).
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---n--- is the number of observations to draw.
	# ---beta_n2--- is the vector of individual intensities, or a single number giving the number of significant scatterers.
	# ---N--- is the number of positive roots of the Bessel function of the first kind.
	# ---max.memory--- is the maximum memory occupied by the function before splitting into a loop for each value of 'r'.
		

	##################################################
	##################################################
	##### Preparation #####
	u=runif(n)
	
	
	##### Execution #####
	out=NAs(n)
	# If 'nsig' is given, this indicates the number of significant scatterers for each simulated value
	if(length(nsig)>0){
		# Repeat 'nsig' if it does not match the 'n':
		if(length(nsig)<n){
			nsig=rep(nsig,length.out=n)
			}
		# Repeat 'magn' if it does not match the 'n':
		if(length(magn)<n){
			magn=rep(magn,length.out=n)
			}
		
		# Simulate one observation at the time:
		for(i in seq_len(n)){
			# Generate the 'beta_n' vector:
			beta_n2=getSeqBarakat(nsig[i],betadistr=betadistr)
			# Return the maximum if only one scatterer is present (degenerage distribution):
			if(length(beta_n2)==1){
				out[i]=magn[i]*beta_n2
				}
			else{
				# Max of h:
				H=sum(sqrt(beta_n2))^2
				# Define the funciton of the Barakat CDF to minimize:			
				fun=function(x){
					abs(pBarakatI(x,beta_n2=nsig[i],N=N,max.memory=max.memory)-u[i])
					}
				out[i]=magn[i] * optimize(fun,lower=0,upper=H)$minimum
				}
			}
		}
	# Else run with the same scatterers for each simulated value:	
	else{
		# Return the maximum if only one scatterer is present (degenerage distribution):
		if(length(beta_n2)==1){
			return(rep(beta_n2,n))
			}
		# Max of h:
		H=sum(sqrt(beta_n2))^2
		# Simulate one observation at the time:
		for(i in seq_len(n)){
			fun=function(x){
				abs(pBarakatI(x,beta_n2=beta_n2,N=N,max.memory=max.memory)-u[i])
				}
			out[i]=magn[i] * optimize(fun,lower=0,upper=H)$minimum
			}
		}
		
		
	##### Output #####
	out
	##################################################
	##################################################
	}
