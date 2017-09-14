#*********************************************
#*********************************************
#' Calculates the probability density function of the total intensity (square of pressure amplitude) from a finite number of sine waves of arbitrary individual amplitudes, given in Equation 64 of Barakat 1974, using the terminology of that paper.
#'
#' @param h  is the argument (superimposed intensity).
#' @param beta_n2  is the vector of individual intensities, or a single number giving the number of significant scatterers.
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
#' @rdname dBarakatI
#'
dBarakatI=function(h,beta_n2,magn=1,N=100,betadistr=c("seq","flat"),max.memory=1e9){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-27 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the probability density function of the total intensity (square of pressure amplitude) from a finite number of sine waves of arbitrary individual amplitudes, given in Equation 64 of Barakat 1974, using the terminology of that paper.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---h--- is the argument (superimposed intensity).
	# ---beta_n2--- is the vector of individual intensities, or a single number giving the number of significant scatterers.
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
	beta_n2=getSeqBarakat(beta_n2,betadistr=betadistr)
	# Apply the magnitude to 'beta_n2':
	beta_n2=magn*beta_n2
	# Transform to the pressure amplitudes used in Barakat 1974:
	beta_n=sqrt(beta_n2)
	rm(beta_n2)
	# Max of r:
	R=sum(beta_n)
	# Set the lower limit for r as the maximum of 0 and max(beta_n)-sum(beta_n[-which.max(beta_n)])
	R0=max(0,max(beta_n)-sum(beta_n[-which.max(beta_n)]))
	# Orders of the Bessel function of the first kind:
	N=seq_len(N)
	
	
	##### Execution #####
	# Roots of the Bessel function of the first kind and zeroth order:
	gamma_n=bessel_zero_J0(N)
	# NUMERATOR of the expression in Barakat, equation 64 (dependent on 31):
	NUMERATOR=phi_y(gamma_n/R,beta_n)
	# Denominator of the expression in Barakat, equation 64 (dependent on 31):
	DENOMINATOR=(R*besselJ(gamma_n,1))^2
	# Factor of the expression in Barakat, equation 64 (dependent on 31):
	FACTOR=besselJ(outer(gamma_n,sqrt(h),"*")/R,0)
	# The expression in Barakat, equation 64 (dependent on 31):
	out=colSums(NUMERATOR / DENOMINATOR * FACTOR)
	out[h<R0^2]=0
	out[h>R^2]=0
	
		
	##### Output #####
	out
	##################################################
	##################################################
	}
