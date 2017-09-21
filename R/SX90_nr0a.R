#*********************************************
#*********************************************
#' Generates close range noise in active mode to be used for fishery sonar. The noise is generated in non-TVG-amplified mode such that when TVG-amplified, the noise follows a pure log r-relationship, i.e., m * log(r). The term "close range noise" is here used somewhat misleading, since it is merely used to mimic the background level of a fishery sonar, which includes surface reverberation and possibly other backscattering which does not fall into the definition of noise in this system (all energy not backscattered from objects).
#'
#' @param beams  is a list of beam configuration data, specifically 
#' @param SNref  is the reference noise level in dB at range 'rref'.
#' @param rref  is the range of the reference noise level.
#' @param SN0  is the noise level at zero range.
#' @param rjoin  is the range to the point where the transient like close range noise (passive mode) joins the active mode close range noise.
#' @param m  specifies the TVG shape of the noise as modeled by 10 m log(r). If m=2, the noise (except absorption) follows the clean background noise, whereas if m=0, the noise (except absorption) is flat (independen of the range form the sonar). The absorption is always present, causing decreasing signal to noise ratio (SNR).
#' @param TVG.exp  a vector of two elements specifying the near field of the noise, where the first element is the noise
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#'
#' @export
#' @rdname SX90_nr0a
#'
SX90_nr0a <- function(beams, SNref=-65, rref=600, SN0=-25, rjoin=20, m=2, TVG.exp=2){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2015-09-02 - Clean version.
	########### DESCRIPTION: ###########
	# Generates close range noise in active mode to be used for fishery sonar. The noise is generated in non-TVG-amplified mode such that whn TVG-amplified, the noise follows a pure log r-relationship, i.e., m * log(r). The term "close range noise" is here used somewhat misleading, since it is merely used to mimic the background level of a fishery sonar, which includes surface reverberation and possibly other backscattering which does not fall into the definition of noise in this system (all energy not backscattered from objects).
	########## DEPENDENCIES: ###########
	#
	############ DETAILS: ############
	# 
	############ VALUE: ############
	# 
	############ REFERENCES: ############
	#
	############ SEAALSO: ############
	#
	############ EXAMPLES: ############
	# event = "/Volumes/Acoustics/S2014119_PG.O.Sars[4174]/Events/S2014119_D200141030_E0022/SX90/tsd"
	# t = 132
	# beams = read.event(event=event, t=t,var="beams")
	# beams$bmmd
	# vbsc = read.event(event=event, t=t,var="vbsc")
	# dimvbsc = dim(vbsc$vbsc)
	# beam1 = 10*log10(vbsc$vbsc[,1])
	# beam21 = 10*log10(vbsc$vbsc[,21])
	# beam51 = 10*log10(vbsc$vbsc[,51])
	# nn = SX90_nr0a(beams, m=0, rjoin=25, SN0=-25, SNref=-62)
	# aa = apply.TVG(nn,beams)
	# beamNoise = 10*log10(aa[seq_len(dimvbsc[1])])
	# ylim = c(beam1, beam21, beam51, beamNoise)
	# ylim = range(ylim[is.finite(ylim)])
	# plotl(beam1, ylim=ylim)
	# lines(beam21, col=2)
	# lines(beam51, col=3)
	# lines(beamNoise, col=4)
	############ VARIABLES: ############
	# ---beams--- is a list of beam configuration data, specifically 
	#	(1) the average speed of sound "asps", 
	#	(2) sampling interval "sint", 
	#	(3) beam length "lenb", 
	#	(4) number of beams "numb", 
	#	(5) and absorption coefficient in units of dB re 1 m^-1 "absr".
	# ---SNref--- is the reference noise level in dB at range 'rref'.
	# ---rref--- is the range of the reference noise level.
	# ---SN0--- is the noise level at zero range.
	# ---rjoin--- is the range to the point where the transient like close range noise (passive mode) joins the active mode close range noise.
	# ---m--- specifies the TVG shape of the noise as modeled by 10 m log(r). If m=2, the noise (except absorption) follows the clean background noise, whereas if m=0, the noise (except absorption) is flat (independen of the range form the sonar). The absorption is always present, causing decreasing signal to noise ratio (SNR).
	# ---TVG.exp--- a vector of two elements specifying the near field of the noise, where the first element is the noise
	
	
	##################################################
	##################################################
	########## Preparation ##########
	#dr = beams$asps*beams$sint/2
	#r = seq(0, max(beams$lenb)-1)*dr
	#r[1] = r[2]/2
	r = soundbeam_range(beams, "mid")
	# Linearize the reference SN:
	SNref = 10^(SNref/10)
	SN0 = 10^(SN0/10)
	# Define 
	# the non-TVG-amplified noise N, and 
	# the TVG-amplified noise 
	# 	SN = N * g, (1)
	# where g = g(r, a, m) = r^m * 10^(r * a/5).
	# Considering the special case that a=0 (no absroption), we model the noise as
	#	SN = N * r^m0
	#	   = C * r^m,
	# yeilding
	#	N = C * r^(m-m0)
	#	N = C * g(r, 0, m-m0)
	# 
	# Given the TVG-amplified noise SNref at range rref, we determine the value of C:
	#	SN(rref) = SNref
	#	       = C * g(rref, 0, m-m0) * g(rref, a, m0)
	#	       = C * rref^(m-m0)  *  rref^m0 * 10^(rref * a/5)
	#	       = C * rref^m * 10^(rref * a/5)
	#	       = C * g(rref, a, m),
	# yealding
	#	C = SNref / g(rref, a, m)
	# 
	# The final model of N is then:
	#	N = SNref * g(r, 0, m-m0)/ g(rref, a, m)
	#
	# Define the range apmlification function (TVG):
	g = function(r, a, m){
		r^m * 10^(r * a/5)
		}
	
	##### Execution and output #####
	noise = list()
	# Get the close range active noise:
	noise$nr0a = SNref * g(r, 0, m-TVG.exp) / g(rref, beams$absr[1], m)
	
	# Add the close range passive noise, given by 'SN0' and 'rjoin':
	atNearFieldEnd = min(which(r>rjoin))-1
	nf = seq_len(atNearFieldEnd)
	noise$nr0a[nf] = seq(SN0, noise$nr0a[atNearFieldEnd]*g(r[atNearFieldEnd], beams$absr[1], TVG.exp), l=atNearFieldEnd) / g(r[nf], beams$absr[1], TVG.exp)
	array(noise$nr0a,dim = c(max(beams$lenb),beams$numb))
	}
