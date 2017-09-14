#*********************************************
#*********************************************
#' Calculates the grid of elevation angle phi, used in integrateonsphere() to make all cells have the same area.
#'
#' @param from  is the lower limit (of phi) of the integral.
#' @param to  is the upper limit (of phi) of the integral.
#' @param area  is the length of the grid vector if length(area)==1, and the relative sizes of the cells, if a vector of length > 1 is given (less useful).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname get.grid_phi
#'
get.grid_phi<-function(from=0,to=pi,area=100){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-01-27 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the grid of elevation angle phi, used in integrateonsphere() to make all cells have the same area.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---from--- is the lower limit (of phi) of the integral.
	# ---to--- is the upper limit (of phi) of the integral.
	# ---area--- is the length of the grid vector if length(area)==1, and the relative sizes of the cells, if a vector of length > 1 is given (less useful).
	
	
	##################################################
	##################################################
	##### Preparation #####
	# The size of the areas of the bands on the sphere, between consecutive values of phi. The whole point of this function is that these areas are set to be constant:
	if(length(area)==1){
		constants=(2:area-1)/area * (cos(from)-cos(to))
		}
	else{
		constants=cumsum(area[-length(area)])/sum(area) * (cos(from)-cos(to))
		}
	
	
	##### Execution and output #####
	# Return the grid, with fixed end points and incriments between that are derived from 
	# 2 * pi * ( cos(phi_i+1) - cos(phi_i) ) = 2 * pi * ( cos(from) - cos(to) ) / l
	# Found from integrating bands of width dphi from phi_(i+1) to phi_i, on the surface of the unit sphere => area = 2 * pi * ( cos(phi_(i+1)) - cos(phi_i) ), and equaling this to the desired area 2 * pi * ( cos(from) - cos(to) ) / l:
	c(from,acos(cos(from)-constants),to)
	##################################################
	##################################################
	}
