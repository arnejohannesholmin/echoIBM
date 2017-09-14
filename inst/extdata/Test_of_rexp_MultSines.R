############################################################
### The below code tests the performance of the function ###
### rexp_MultSines() with respect to the histogram of the ##
### squared values (which should be exponential), the acf ##
### along beams, and the correlation across beams, which ###
### should reflect the estimated from real data. The test ##
### considers different values of the number 'L' of sine ###
### waves located in each sampling interval (which should ##
### affect the shape of the histogram and result in non- ###
### exponential histogram for the squared values for low ###
### values of 'L', the number of sample points 'N' in the ##
####### sample intervals, the number of periods 'P'  #######
############################################################




# Define the grid values:
L=c(1,5,10,20,50)
N=c(10,100,1000)
P=c(5,10,20)
cputime=NAs(length(L),length(N),length(P))

ks=NAs(length(L),length(N),length(P))
for(i1 in seq_along(L)){
	for(i2 in seq_along(N)){
		for(i3 in seq_along(P)){
			cputime[i1,i2,i3]<-system.time(s<-rexp_MultSines(J=10000,I=1,L=L[i1],N=N[i2],P=P[i3],w=3,olpn=1,seed=runif(1)*1e9)^2)[3]
			if(!any(is.infinite(s),is.na(s))){
				cat(c(i1,i2,i3),"\n")
				ks[i1,i2,i3] = ks.test(s,pexp,1/mean(s))$p
				}
			else{
				cat("\t",c(i1,i2,i3),"\n")
				}
			}
		}
	}

# These are the printouts from three trials, where the last also repported time used by the function. We see that for all three trials the combination L=5, N=100, P=20 gives good results (0.7352132, 0.7898040, 0.93083043). For this combination the time used is only 0.687 seconds, and we thus take this combination as a starting point for fine tuning of the parameters.

# > ks
, , 1

     [,1]      [,2]       [,3]
[1,]    0 0.0000000 0.00000000
[2,]    0 0.0136028 0.61938359
[3,]    0 0.6092815 0.41871558
[4,]    0 0.4093084 0.09469639
[5,]    0 0.3026382 0.15621922

, , 2

     [,1]      [,2]      [,3]
[1,]    0 0.0000000 0.0000000
[2,]    0 0.1365101 0.4400649
[3,]    0 0.5868748 0.1114223
[4,]    0 0.3061418 0.8695362
[5,]    0 0.4816164 0.8172545

, , 3

     [,1]      [,2]       [,3]
[1,]    0 0.0000000 0.00000000
[2,]    0 0.7352132 0.07234405
[3,]    0 0.5378454 0.78467838
[4,]    0 0.7438617 0.61049766
[5,]    0 0.5202529 0.94399112



# > ks
, , 1

     [,1]      [,2]       [,3]
[1,]    0 0.0000000 0.00000000
[2,]    0 0.2638863 0.02702015
[3,]    0 0.4424839 0.37233987
[4,]    0 0.1651657 0.48758985
[5,]    0 0.4848821 0.50813158

, , 2

     [,1]      [,2]       [,3]
[1,]    0 0.0000000 0.00000000
[2,]    0 0.4723880 0.05543447
[3,]    0 0.4413080 0.85943110
[4,]    0 0.5980443 0.94470336
[5,]    0 0.2839800 0.83313533

, , 3

     [,1]         [,2]       [,3]
[1,]    0 1.110223e-16 0.00000000
[2,]    0 7.898040e-01 0.09572924
[3,]    0 8.504530e-01 0.77141444
[4,]    0 9.692905e-01 0.53094481
[5,]    0 9.921026e-01 0.87782566











L=c(1,5,10,20,50)
N=c(10,100,1000)
P=c(5,10,20)
# > ks
, , 1

     [,1]       [,2]      [,3]
[1,]    0 0.00000000 0.0000000
[2,]    0 0.04396717 0.1768777
[3,]    0 0.35463934 0.3293339
[4,]    0 0.63657947 0.3433057
[5,]    0 0.11268029 0.5653756

, , 2

     [,1]      [,2]      [,3]
[1,]    0 0.0000000 0.0000000
[2,]    0 0.2392148 0.1479066
[3,]    0 0.4027709 0.3976411
[4,]    0 0.1708620 0.7835222
[5,]    0 0.9954830 0.9040878

, , 3

     [,1]       [,2]      [,3]
[1,]    0 0.00000000 0.0000000
[2,]    0 0.93083043 0.5602019
[3,]    0 0.09106283 0.7550364
[4,]    0 0.20450995 0.7749663
[5,]    0 0.87120481 0.4183957

# > cputime
, , 1

      [,1]  [,2]   [,3]
[1,] 0.016 0.152  1.594
[2,] 0.072 0.688  6.574
[3,] 0.140 1.361 13.228
[4,] 0.279 2.678 26.670
[5,] 0.683 6.660 65.843

, , 2

      [,1]  [,2]   [,3]
[1,] 0.016 0.152  1.586
[2,] 0.071 0.689  6.626
[3,] 0.138 1.375 13.327
[4,] 0.277 2.678 26.681
[5,] 0.679 6.682 65.911

, , 3

      [,1]  [,2]   [,3]
[1,] 0.016 0.159  1.593
[2,] 0.071 0.687  6.664
[3,] 0.139 1.384 13.392
[4,] 0.276 2.679 26.775
[5,] 0.679 6.728 66.032






# Fine tune the parameters based on L=5, N=100, P=20:
L=c(2,3,4,5,6,7,8)
N=c(20,40,60,80,100,120,140,170,200)
P=20
# Number of runs:
M=10
cputimeFine=NAs(length(L),length(N))

ksFine=NAs(length(L),length(N),M)
for(i1 in seq_along(L)){
	for(i2 in seq_along(N)){
		for(i3 in seq_len(M)){
			cputimeFine[i1,i2]<-system.time(s<-rexp_MultSines(J=10000,I=1,L=L[i1],N=N[i2],P=P,w=3,olpn=1,seed=runif(1)*1e9)^2)[3]
			if(!any(is.infinite(s),is.na(s))){
				cat(c(i1,i2),"\n")
				ksFine[i1,i2,i3] = ks.test(s,pexp,1/mean(s))$p
				}
			else{
				cat("\t",c(i1,i2),"\n")
				}
			}
		}
	}
	
	


L=c(2,3,4,5,6,7,8)
N=c(20,40,60,80,100,120,140,170,200)
# > apply(ksFine,1:2,mean)
     [,1] [,2]         [,3]      [,4]       [,5]        [,6]        [,7]       [,8]       [,9]
[1,]    0    0 1.895573e-03 0.2509040 0.03618491 0.007611678 0.003292011 0.01030188 0.01106903
[2,]    0    0 3.318096e-04 0.8302697 0.31433887 0.252308016 0.141119472 0.06910602 0.18039502
[3,]    0    0 3.375390e-05 0.6562473 0.49043803 0.305746559 0.276714419 0.11860527 0.18285117
[4,]    0    0 1.257896e-05 0.4740522 0.42477367 0.694641839 0.371759373 0.33956520 0.23683710
[5,]    0    0 1.126440e-07 0.4200197 0.65231051 0.638391299 0.496063336 0.41467648 0.56547026
[6,]    0    0 1.175986e-06 0.3762552 0.56494512 0.514027529 0.489464154 0.44690116 0.43477725
[7,]    0    0 3.216059e-06 0.4015749 0.57588331 0.606062064 0.567433905 0.60542098 0.51990038

# > apply(ksFine,1:2,sd)
     [,1] [,2]         [,3]      [,4]       [,5]       [,6]      [,7]       [,8]       [,9]
[1,]    0    0 2.992040e-03 0.2141203 0.05092072 0.01043162 0.0039076 0.01119361 0.01310941
[2,]    0    0 8.605872e-04 0.1213590 0.24576666 0.21433852 0.2024499 0.05998922 0.19358641
[3,]    0    0 6.380919e-05 0.2151500 0.28251376 0.32307021 0.2751465 0.13440806 0.25055917
[4,]    0    0 2.575826e-05 0.3019557 0.26217061 0.18575907 0.2652508 0.28602619 0.22092564
[5,]    0    0 1.518099e-07 0.2732722 0.25260964 0.26893930 0.3181547 0.33927307 0.27921919
[6,]    0    0 3.310748e-06 0.2872471 0.28141639 0.28583488 0.3969625 0.33205969 0.34892372
[7,]    0    0 1.013804e-05 0.2618063 0.27362149 0.18228192 0.3138867 0.31041089 0.27261371

# > cputimeFine
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]
[1,] 0.057 0.112 0.169 0.225 0.280 0.336 0.403 0.480 0.559
[2,] 0.086 0.168 0.253 0.392 0.506 0.503 0.607 0.803 1.206
[3,] 0.203 0.245 0.348 0.533 0.673 1.403 1.242 1.549 1.767
[4,] 0.281 0.513 0.650 0.750 0.785 0.868 0.967 1.180 1.381
[5,] 0.175 0.346 0.572 0.660 0.838 1.024 1.149 1.425 1.654
[6,] 0.192 0.383 0.593 0.764 0.955 1.163 1.334 1.619 1.978
[7,] 0.219 0.441 0.657 0.871 1.094 1.499 1.895 1.855 2.174	










# Quick test of the dependence of N:
L=3
N=c(20,40,60,80,100,120,140,170,200)
P=10
# Number of runs:
M=10
cputimeFine10=NAs(length(L),length(N))

ksFine10=NAs(length(L),length(N),M)
for(i1 in seq_along(L)){
	for(i2 in seq_along(N)){
		for(i3 in seq_len(M)){
			cputimeFine10[i1,i2]<-system.time(s<-rexp_MultSines(J=10000,I=1,L=L[i1],N=N[i2],P=P,w=3,olpn=1,seed=runif(1)*1e9)^2)[3]
			if(!any(is.infinite(s),is.na(s))){
				cat(c(i1,i2),"\n")
				ksFine10[i1,i2,i3] = ks.test(s,pexp,1/mean(s))$p
				}
			else{
				cat("\t",c(i1,i2),"\n")
				}
			}
		}
	}
	
# > apply(ksFine10,1:2,mean)
     [,1]      [,2]      [,3]      [,4]      [,5]       [,6]      [,7]       [,8]      [,9]
[1,]    0 0.4758227 0.1336829 0.1210044 0.1228057 0.06423378 0.1016342 0.08677021 0.1445408

# > apply(ksFine10,1:2,sd)
     [,1]      [,2]      [,3]      [,4]      [,5]      [,6]       [,7]       [,8]      [,9]
[1,]    0 0.2817212 0.2072054 0.1431279 0.1356394 0.1080381 0.09966389 0.07395362 0.2128694

# > cputimeFine10
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]
[1,] 0.086 0.171 0.255 0.356 0.422 0.506 0.584 0.707 0.834	
	



### Aha, it is of some importance where the sampling points are, and if N=4*P, we expect sampling points to end up close to peaks along the sine waves. We move on with this restriction:
L=3
P=c(5,10,20,50)
N=4*P
# Number of runs:
M=100
cputimeFine4P=NAs(length(L),length(N))

ksFine4P=NAs(length(L),length(N),M)
for(i1 in seq_along(L)){
	for(i2 in seq_along(N)){
		for(i3 in seq_len(M)){
			cputimeFine4P[i1,i2]<-system.time(s<-rexp_MultSines(J=10000,I=1,L=L[i1],N=N[i2],P=P[i2],w=3,olpn=1,seed=runif(1)*1e9)^2)[3]
			if(!any(is.infinite(s),is.na(s))){
				cat(c(i1,i2),"\n")
				ksFine4P[i1,i2,i3] = ks.test(s,pexp,1/mean(s))$p
				}
			else{
				cat("\t",c(i1,i2),"\n")
				}
			}
		}
	}

# It seems P=5 and N=20 works:
apply(ksFine4P,1:2,mean)
          [,1]      [,2]     [,3]     [,4]
[1,] 0.5136558 0.5930814 0.637432 0.665327
apply(ksFine4P,1:2,sd)/sqrt(M)
           [,1]       [,2]       [,3]       [,4]
[1,] 0.02408024 0.02740568 0.02661452 0.02615342
cputimeFine4P
      [,1]  [,2]  [,3]  [,4]
[1,] 0.087 0.176 0.351 0.846

# Beautifull histogram of the last value of P:
histexp(s,breaks=400)

# Although the p-value in the Kolmogorov-Smirnov-test seems to increase with increasing P and thus increasing N, the gain seems to be largest from P=5 to P=10. We thus choose P=10, N=40, and L=3, which results in the following time usage for one fan of the MS70:




system.time(s<-rexp_MultSines(J=1324,I=25,L=3,N=40,P=10,w=3,olpn=c(0.5,1,0.5),seed=runif(1)*1e9))
#   user  system elapsed 
#  0.800   0.001   0.798 
  
dim_all(s)

library(fields)
image.plot(olpn(s))
pp(5,5,oma=ones(4),mar=ones(4))
for(i in 1:25){
	acf(s[,i]^2,lag.max=6)
	abline(h=0.4)
	}

pp(5,5,oma=ones(4),mar=ones(4))
for(i in 1:25){
	histexp(s[,i]^2,breaks=100)
	}

# Brilliant, and the function only used 0.8 seconds, which is 16 seconds for one ping and approximately one hour for all the 240 pings of the second paper. Yey!








## Finally tune the correlation to 0.4, as determined to be an apropriate value for all beams in the MS70 in the script "S2009116_PG.O.Sars[4174] - acf and correlation":

lag1=arr.ind2ind(cbind(2:25,1:24),c(25,25))
cor1=c(0.4,0.5,0.6)
corout=zeros(length(cor1))
for(i in seq_along(cor1)){
	system.time(s<-rexp_MultSines(J=10000,I=25,L=3,N=40,P=10,w=3,olpn=c(cor1[i],1,cor1[i]),seed=runif(1)*1e9))
	corout[i]=mean(olpn(s)[lag1])
	}

# Finer grid around 0.5:
lag1=arr.ind2ind(cbind(2:25,1:24),c(25,25))
cor1=seq(0.45,0.5,0.01)
corout=zeros(length(cor1))
for(i in seq_along(cor1)){
	system.time(s<-rexp_MultSines(J=10000,I=25,L=3,N=40,P=10,w=3,olpn=c(cor1[i],1,cor1[i]),seed=runif(1)*1e9))
	corout[i]=mean(olpn(s)[lag1])
	}
> corout
[1] 0.3945725 0.4012673 0.4083585 0.4103687 0.4171318
[6] 0.4238842

# We thus choose c(0.46,1,0.46) to obtain the correlation 0.4:







# Finally we need to scale the results to obtain expectation equal to 1 for the combination of parameters. For full resolution the expectation will be very close to the number of sine waves affecting each sampling interval (product of 'L' and 'w') but for the parameters found above, this is no longer valid:

system.time(s<-rexp_MultSines(J=1e6,I=1,L=3,N=400,P=10,w=3,olpn=1,seed=runif(1)*1e9))
#    user  system elapsed 
# 169.554   0.403 172.988 
# > mean(s)
# [1] 8.963494

system.time(s<-rexp_MultSines(J=1e6,I=1,L=3,N=40,P=10,w=3,olpn=1,seed=runif(1)*1e9))
#    user  system elapsed 
#  17.455   0.042  17.686 
# > mean(s)
# [1] 7.352798

# We spend some computational time on getting an accurate estimate of the expectation for the given combination:
# First get an impression of the variation:
means=zeros(100)
seeds=1:100
for(i in seq_along(means)){
	cat(i,"\ ")
	s=rexp_MultSines(J=1e4,I=1,L=3,N=40,P=10,w=3,olpn=1,seed=seeds[i])
	means[i]=mean(c(s)^2)
	}
# > sd(means)
# [1] 0.09636948
# > cvar(means)
# [1] 0.01310664

# If we increase the number of voxels to 10 runs of 1e7 the coeficient of variation is reduced to 0.01310664/sqrt(10*100) = 0.0001310664, which we feel is acceptable:
lmeans=10
means=zeros(lmeans)
seeds=seq_len(lmeans)
for(i in seq_along(means)){
	print(system.time(s<-rexp_MultSines(J=1e7,I=1,L=3,N=40,P=10,w=3,olpn=1,seed=seeds[i])))
	means[i]=mean(s)
	}

mean(means)
# > mean(means)
# [1] 7.350175

# This value is added to the function for the MS70





########## However, for different values of 'w' we need a table linking 'w' and 'scale': ##########
w=c(1, 5, 10, 20, 40, 75, 130, 200)
lw=length(w)
nruns=1
means=zeros(lw,nruns)
cputime=zeros(lw)

seeds=array(seq_len(lw*nruns),dim=c(lw,nruns))
for(i in seq_len(lw)){
	cat("w = ",w[i]," (",i,"/",lw,")\n",sep="")
	ppp=proc.time()[3]
	for(j in seq_len(nruns)){
		s<-rexp_MultSines(J=1e7,I=1,L=3,N=40,P=10,w=w[i],olpn=1,seed=seeds[i,j])
		means[i,j]=mean(s)
		}
	cputime[i]=proc.time()[3]-ppp
	print(cputime[i])
	}

# Observe that the scale increases linearly by increases in 'w':
ploto(w,means)
l=lm(means ~ 0 + w)
abline(l,col=3)
# Coefficients:
#     w  
# 2.455 

# Write the values to file:
w_scale=cbind(w,means)
colnames(w_scale)=c("w","scale")
dump("w_scale",file.path("/Applications/echoIBM/Frameworks/R/Functions/echoIBM Main/echoIBM","Table of w and scale.R"))

