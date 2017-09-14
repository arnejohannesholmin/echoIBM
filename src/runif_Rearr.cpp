#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdlib>

extern "C" {
	
	//////////////////////// AUTHOR(S): ////////////////////////
	// Arne Johannes Holmin
	//////////////////////// LANGUAGE: /////////////////////////
	// English
	/////////////////////////// LOG: ///////////////////////////
	// Start: 2011-11-15 - Clean version.
	////////////////////// DESCRIPTION: ////////////////////////
	// Generates correlated uniform vectors that are internally autocorrelated, based on the Rearrangement method given by Holmin.
	//////////////////////// VARIABLES: ////////////////////////
	// int *n		[1]		The number of valiables to generate in each vector.
	// int *p		[1]		The number of vectors.
	// double *C	[p,p]	A "correlation" matrix determining the degree of correlation between vectors. These values do not fully correspond to correlations, and trial and error is required to obtain the desired correlation. The method is subject to limitations concerning the acheivable autocorrelation and correlation.
	// int *wC		[1]		The number of significant vectors correlated to each vector. This value is used for reducing the operations of correlating insignificantly correlated vectors.
	// int *Cind	[n,wC]	A matrix of indeces soecifying which vectors are significantly correlated. For simplicity the number of significant vectors is constant for all vectors, so that indexes for zeros-correlation vectors may occur.
	// int *w		[1]		The number of consecutive values including the present value. If w=2, the method choses between the current and the next value depending on the proximity to the prevoius values.
	// double *prob	[w]		Probability values assigned to each of the consecutive values sorted by their proximity to the prevoius value. If w = 3 and prob = [0.7,0.5,0,5] the closest value will be assigned to the current position with probability 0.7, and if this fails, the second closest value will be assigned to the current position with probability 0.5, and so on.
	// int *l		[1]		The number of values included when randomly rearranging the values not assigned to the current position. We may include more values than 'w' too reduce the effect of nevative autocorrelation at lags ≈ 5.
	// int *sample	[n,l]	A matrix of prevoiusly generated sampled (without replacement) index vectors used when rearranging the values not assigned to the current position. There is one index vector of 'l' indices for each step in the rearrangement iteration.
	// int *seed	[1]		The random seed of the utilization of 'prob'.
	// double *Z	[n,p]	The matrix of the independent uniform vectors.
	// double *U	[n,p]	The output matrix of autocorrelated and correlated uniform vectors.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array indexes and the number (n) of rows of the array:
	int runif_Rearr_ind(int i, int j, int n)
	{
		return i + j * n;
	}
	
	// The main funciton for generating correlated vectors of autocorrelated uniform variables. Parameters to the function are:
	void runif_Rearr(int *n, int *p, double *C, int *wC, int *Cind, int *w, double *prob, int *l, int *sample, int *seed, double *Z, double *U)
	{
		// Define variables used in the function ('f' is the position of the closest value to the prevoius and 'thisrand' is the random number in [0,1]):
		int f = 0;
		double thisrand;
		// Define the sum of the difference between the proceding elements along the beam and the current elements across the beams, weighted by C. 'd' is the sum of the differences and 's' will be the sorted version of 'd':
		double d[*w], s[*w];
		
		// Set the random seed:
		srand(*seed);
		
		// Initiate the output and assign the independent uniform values to 'U':	
		for(int i = 0; i < *n * *p; i++)
		{
			U[i] = 0.0;
		}
		for(int i = 0; i < *n * *p; i++)
		{
			U[i] = Z[i];
		}
		
		// Run through the voxels and rearrange independent uniformly distributied variables:
		for (int i1 = 1; i1 < *n - *w; i1++) // along beams
		{
			for (int i2 = 0; i2 < *p; i2++) // across beams
			{
				// Clean up 'd' and 's':
				for (int i3 = 0; i3 < *w; i3++) // along the window of voxels along the current beam
				{
					s[i3] = 0.0;
					d[i3] = 0.0;
				} // End of for i3
				
				// Determime the summed distance to the previous values:
				for (int i3 = 0; i3 < *w; i3++) // along the window of voxels along the current beam
				{
					for (int i4 = 0; i4 < *wC; i4++) // across beams having positive correlation to the current beam 
					{
						int thisind = Cind[ runif_Rearr_ind(i2,i4,*p) ];
						if(thisind >= 0 && thisind < *p)
						{
							d[i3] += C[runif_Rearr_ind(i2, Cind[ runif_Rearr_ind(i2,i4,*p) ], *p)] * abs( U[runif_Rearr_ind(i1+i3,i2,*n)] - U[runif_Rearr_ind(i1-1, Cind[ runif_Rearr_ind(i2,i4,*p) ], *n)] );
						}
					} // End of for i4
				} // End of for i3
				
				// Copy 'd' to 's' and sort 's':
				for (int i3 = 0; i3 < *w; i3++) // along the window of voxels along the current beam
				{
					s[i3] = d[i3];
				} // End of for i3
				sort(s, s + *w);
				
				// Move through 'd' and switch if 'prob' is satified:
				int pos = 0;
				while(pos < *w)
				{
					// Get the current value in the sorted version of 'd' (named 's'):
					int posf = 0;
					while(1)
					{
						if(d[posf] == s[pos])
						{
							f=posf;
							break;
						}
						posf += 1;
					}
					
					// Generate random number and compare to 'prob':
					thisrand = rand() / double(RAND_MAX);
					
					// If f>0 the current value should be replaced by value number 'f':
					if(f==0 && thisrand < prob[pos])
					{
						break;
					}
					else if(thisrand < prob[pos])
					{
						// Create a temporary array of length 'l' to hold the uniform values and insert from this array below:
						double temp[*l];
						for (int i5 = 0; i5 < *l; i5++){
							temp[i5] = U[runif_Rearr_ind(i1+i5,i2,*n)];
						}
						
						// Move through the 'l' values and rearrange if the current sample value is not 'f':
						// 'sample' is a matrix of rows representing sample realizations of the vector 0:(l-1)
						int posl = 0;
						for (int i5 = 0; i5 < *l; i5++)
						{
							// Define the current sample value for simplicity:
							int thissample = sample[runif_Rearr_ind(i1,i5,*n)];
							if(thissample != f  &&  i1+2+thissample <= *n)
							{
								U[runif_Rearr_ind(i1+posl+1,i2,*n)] = temp[thissample];
								posl += 1;
							} // End of if sample[i1* *l + k] != f
						} // End of for k
						// Insert the closest accepted value to the current position:
						U[runif_Rearr_ind(i1,i2,*n)] = temp[f];
						// Break the while-loop:
						pos = *w;
					} // End of if f>0  &&  thisrand < prob[pos]
					pos += 1;	
				} // End of while(pos<w)
			} // End of for i2
		} // End of for i1
	} // End of void
} // End of extern "C"
