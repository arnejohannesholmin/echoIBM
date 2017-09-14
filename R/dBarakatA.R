#*********************************************
#*********************************************
#' Calculates the probability density function of the pressure amplitude from a finite number of sine waves of arbitrary individual amplitudes, given in Equation 56 of Barakat 1974, using the terminology of that paper.
#'
#' @param r  is the argument (superimposed pressure amplitude).
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
#' @rdname dBarakatA
#'
dBarakatA=function(r,beta_n,magn=1,N=100,betadistr=c("seq","flat"),max.memory=1e9){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-27 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the probability density function of the pressure amplitude from a finite number of sine waves of arbitrary individual amplitudes, given in Equation 56 of Barakat 1974, using the terminology of that paper.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---r--- is the argument (superimposed pressure amplitude).
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
	# NUMERATOR of the expression in Barakat, equation 56 (dependent on 55 and 31):
	NUMERATOR=phi_y(gamma_n/R,beta_n)
	# Denominator of the expression in Barakat, equation 56 (dependent on 55 and 31):
	DENOMINATOR=(R*besselJ(gamma_n,1))^2
	# Factor of the expression in Barakat, equation 56 (dependent on 55 and 31):
	if(.Platform$r_arch=="x86_64" && length(gamma_n)*length(r)*16<max.memory || .Platform$r_arch=="i386" && length(gamma_n)*length(r)*8<max.memory){
		FACTOR=besselJ(outer(gamma_n,r,"*")/R,0)
		# The expression in Barakat, equation 56 (dependent on 55 and 31):
		out=2*r*colSums(NUMERATOR / DENOMINATOR * FACTOR)
		out[r<R0 | r>R]=0
		}
	else{
		out=r
		for(i in seq_along(r)){
			# Factor of the expression in Barakat, equation 56 (dependent on 55 and 31):
			FACTOR=besselJ(gamma_n*r[i]/R,0)
			# The expression in Barakat, equation 56 (dependent on 55 and 31):
			out[i]=2*r[i]*sum(NUMERATOR / DENOMINATOR * FACTOR)
			}
		out[r<R0 | r>R]=0
		}
	
		
	##### Output #####
	out
	##################################################
	##################################################
	}
