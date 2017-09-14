#*********************************************
#*********************************************
#' Expands the background noise into a matrix corresponding to the dimension of the data.
#'
#' @param data  is a list containing beam configuration, and periodic noise parameters. Speficically the following variables must be included: Beams variables: 'sint', 'lenb'; periodic noise variables: 'pns1', 'pns2', 'pns3', 'acfq', 'harm'; and either 'bgns' or 'numb'
#' @param TVG  is TRUE if TVG should be added to the periodic noise.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw apply.TVG
#'
#' @export
#' @rdname echoIBM.bgns
#'
echoIBM.bgns<-function(data,TVG=FALSE,TVG.exp=2){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-09-25 - Clean version.
	########### DESCRIPTION: ###########
	# Expands the background noise into a matrix corresponding to the dimension of the data.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing beam configuration, and periodic noise parameters. Speficically the following variables must be included: Beams variables: 'sint', 'lenb'; periodic noise variables: 'pns1', 'pns2', 'pns3', 'acfq', 'harm'; and either 'bgns' or 'numb'
	# ---TVG--- is TRUE if TVG should be added to the periodic noise.
	# ---TVG.exp--- is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
	
	
	##################################################
	##################################################
	if(TVG){
		data$bgns=apply.TVG(matrix(data$bgns,nrow=max(data$lenb),ncol=data$numb,byrow=TRUE),data,TVG.exp=TVG.exp)
		}
	else{
		data$bgns=matrix(data$bgns,nrow=max(data$lenb),ncol=data$numb,byrow=TRUE)
		}
	##################################################
	##################################################
	}
