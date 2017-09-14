#*********************************************
#*********************************************
#' Calculates the cummulative probability function of the pressure amplitude from a finite number of sine waves of arbitrary individual amplitudes, given (for the upper probability) in Equation 59 of Barakat 1974, using the terminology of that paper.
#'
#' @param r_0  is the argument (superimposed pressure amplitude).
#' @param beta_n  is the vector of individual pressure amplitudes, or a single number giving the number of significant scatterers.
#' @param N  is the number of positive roots of the Bessel function of the first kind.
#' @param max.memory  is the maximum memory occupied by the function before splitting into a loop for each value of 'r'.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom gsl bessel_zero_J0
#'
#' @export
#' @rdname pBarakatA
#'
pBarakatA=function(r_0,beta_n,magn=1,N=100,betadistr=c("seq","flat"),max.memory=1e9){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-27 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the cummulative probability function of the pressure amplitude from a finite number of sine waves of arbitrary individual amplitudes, given (for the upper probability) in Equation 59 of Barakat 1974, using the terminology of that paper.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---r_0--- is the argument (superimposed pressure amplitude).
	# ---beta_n--- is the vector of individual pressure amplitudes, or a single number giving the number of significant scatterers.
	# ---N--- is the number of positive roots of the Bessel function of the first kind.
	# ---max.memory--- is the maximum memory occupied by the function before splitting into a loop for each value of 'r'.
		

	##################################################
	##################################################
	##### Preparation #####
	# Function calculating the characteristic function of the sum of scatterers (Barakat 1974, equation 31):
	phi_y=function(omega,beta_n){
		apply(besselJ(outer(beta_n,omega,"*"),0),2,prod)
		}
	# Generate the 'beta_n' vector:
	beta_n=magn*getSeqBarakat(beta_n,betadistr=betadistr)
	# Max of r:
	R=sum(beta_n)
	# Set the lower limit for r as the maximum of 0 and max(beta_n)-sum(beta_n[-which.max(beta_n)])
	R0=max(0,max(beta_n)-sum(beta_n[-which.max(beta_n)]))
	# Orders of the Bessel function of the first kind:
	N=seq_len(N)
	
	
	##### Execution #####
	# Roots of the Bessel function of the first kind and zeroth order:
	gamma_n=bessel_zero_J0(N)
	# NUMERATOR of the expression in Barakat, equation 59 (dependent on 31):
	NUMERATOR=phi_y(gamma_n/R,beta_n)
	# Denominator of the expression in Barakat, equation 59 (dependent on 31):
	DENOMINATOR=gamma_n*(besselJ(gamma_n,1))^2
	# Factor of the expression in Barakat, equation 59 (dependent on 31):
	if(.Platform$r_arch=="x86_64" && length(gamma_n)*length(r_0)*16<max.memory || .Platform$r_arch=="i386" && length(gamma_n)*length(r_0)*8<max.memory){
		FACTOR=besselJ(outer(gamma_n,r_0,"*")/R,1)
		# The expression in Barakat, equation 59 (dependent on 31):
		out=2*r_0/R * colSums(NUMERATOR / DENOMINATOR * FACTOR)
		out[r_0<R0]=0
		out[r_0>R]=1
		}
	else{
		out=r_0
		for(i in seq_along(r_0)){
			# Factor of the expression in Barakat, equation 59 (dependent on 31):
			FACTOR=besselJ(gamma_n*r_0[i]/R,1)
			# The expression in Barakat, equation 59 (dependent on 31):
			out[i]=2*r_0[i]/R * sum(NUMERATOR / DENOMINATOR * FACTOR)
			}
		out[r_0<R0]=0
		out[r_0>R]=1
		}
	
		
	##### Output #####
	out
	##################################################
	##################################################
	}
