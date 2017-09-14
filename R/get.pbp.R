get.pbp <- function(x){
	x = tolower(x)
	if(x == "ps"){
		fun=pointSource.TSD
		}
	else if(x == "ls"){
		fun=lineSource.TSD
		}
	else if(x == "PrSph"){
		#fun=lineSource.TSD
		}
	else if(x == "cp"){
		fun=circularPiston.TSD
		}
	else if(x == "cpe"){
		fun=circularPiston_ellipticRadius.TSD
		}
	else if(x == "cpes"){
		fun=circularPiston_ellipticRadius_sidelobefit.TSD
		}
	else{
		warning("Parametric beam pattern not recognized (must be one of \"pointSource\", \"lineSource\", \"circularPiston\", \"circularPiston_ellipticRadius\", \"circularPiston_ellipticRadius_sidelobefit.TSD\")")
		}
	}
