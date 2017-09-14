#*********************************************
#*********************************************
#' Generates 'beta_n' used in dBarakatA(), pBarakatA(), dBarakatI(), pBarakatI(), rBarakatA(), rBarakatI().
#'
#' @param beta_n  is the vector of individual pressure amplitudes, or a single number giving the number of significant scatterers.
#' @param betadistr  is an identifier specifying how to generate 'beta_n' in the case that 'beta_n' is a single numeric.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ones
#'
#' @export
#' @rdname getSeqBarakat
#'
getSeqBarakat<-function(beta_n,betadistr=c("seq","flat")){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-29 - Clean version.
	########### DESCRIPTION: ###########
	# Generates 'beta_n' used in dBarakatA(), pBarakatA(), dBarakatI(), pBarakatI(), rBarakatA(), rBarakatI().
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---beta_n--- is the vector of individual pressure amplitudes, or a single number giving the number of significant scatterers.
	# ---betadistr--- is an identifier specifying how to generate 'beta_n' in the case that 'beta_n' is a single numeric.
		

	##################################################
	##################################################
	##### Preparation and execution #####
	# 'beta_n' may be given as the number of significant scatterers. In this case 'beta_n' is generated  depending on the value of 'betadistr': If betadistr=="seq" a sequence of values from a_1 to 1 is generated which sums to 'beta_n', in such a way that a_1 is as low as possible. If betadistr=="flat" 'beta_n' is generated as a vector of ones followed by the remainder:
	if(length(beta_n)==1){
		if(tolower(substr(betadistr[1],1,1))=="s"){
			n=seq(ceiling(beta_n),floor(2*beta_n))
			a=2*beta_n/n-1
			am=which.min(a)
			beta_n=seq(a[am],1,length.out=n[am])
			beta_n=beta_n[beta_n>0]
			}
		else{
			beta_n=c(ones(beta_n),beta_n%%1)
			}
		}
	
		
	##### Output #####
	beta_n
	##################################################
	##################################################
	}
