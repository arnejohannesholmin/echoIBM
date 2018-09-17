#*********************************************
#*********************************************
#' Function used for approximating the beam pattern of the beam pattern from a circular piston. The values past the 8 side lobe is approximated by a line (see examples).
#'
#' @param x  is a vector.
#'
#' @return
#'
#' @examples
#' \dontrun{
#' x <- seq(0, 100, length.out=1e5)
#' system.time(c1 <- circularPiston_appr(x))
#' system.time(c2 <- (2 * besselJ(x, 1) / x)^2)
#' plot(10*log10(c1), type="l")
#' lines(10*log10(c2), type="l", col=2)
#'}
#'
#' @export
#' @rdname circularPiston_appr
#'
besselJ1_appr_angmax <- 30
besselJ1_appr_l <- 1e6
besselJ1_appr_d <- besselJ1_appr_angmax / besselJ1_appr_l
besselJ1_appr_ang <- seq(0.5 * besselJ1_appr_d, (besselJ1_appr_l - 0.5) * besselJ1_appr_d, l=besselJ1_appr_l)
besselJ1_appr_b <- (2 * besselJ(besselJ1_appr_ang,1) / besselJ1_appr_ang)^2

circularPiston_appr <- function(x){
	
	############### LOG: ###############
	# Start: 2012-10-24 - Clean version.
	
	# Values used avobe the maximum angle, to approximate the beam pattern:
	intersect_above_angmax <- -32.5
	slope_above_angmax <- -0.275
	# Get indices for the argument in the predefined angle and beam pattern vectors:
	at <- round(x / besselJ1_appr_angmax * besselJ1_appr_l + 0.5)
	# Set index 0 to 1:
	if(any(at==0)){
		at[at==0] <- 1
	}
	out <- double(length(at))
	above_angmax <- x > besselJ1_appr_angmax
	if(any(above_angmax)){
		out[!above_angmax] <- besselJ1_appr_b[at[!above_angmax]]
		out[above_angmax] <- 10^((intersect_above_angmax + slope_above_angmax * x[above_angmax])/10)
	}
	else{
		out <- besselJ1_appr_b[at]
	}
	out
}
