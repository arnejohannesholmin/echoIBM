#include <iostream>
using namespace std;
#include <valarray>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string>

extern "C" {
	
	//////////////////////// AUTHOR(S): ////////////////////////
	// Arne Johannes Holmin
	//////////////////////// LANGUAGE: /////////////////////////
	// English
	/////////////////////////// LOG: ///////////////////////////
	// Start: 2012-01-22 - Clean version.
	// Update: 2012-09-21 - Added the option of variable correlation between the beams, consistent with the need for simulating the MS70 data with high correlation between beams in the center of fans with periodic noise.
	// Last: 2012-09-21 - Added the option of variable correlation between and along the beams.
	////////////////////// DESCRIPTION: ////////////////////////
	// Generates correlated rayleigh vectors that are internally autocorrelated, based on superimposed sine waves.
	//////////////////////// VARIABLES: ////////////////////////
	// int *n		[1]		The number of valiables to generate in each vector.
	// int *p		[1]		The number of vectors.
	// double *U	[n,p]	The output matrix of autocorrelated and correlated uniform vectors.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int rexp_MultSines_VaribleCorrelation_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	int rexp_MultSines_VaribleCorrelation_ind3d(int i, int j, int k, int Li, int Lj)
	{
		return i + j * Li + k * Li * Lj;
	}
	
	// A modified sine function returning 0 outside the limits of the set on which the function is defined:
	double rexp_MultSines_VaribleCorrelation_sin0(double x, double omega, double maxx)
	{
		double out;
		if(x<0 || x>maxx)
		{
			out=0.0;
		}
		else
		{
			out = sin(omega*x);
		}
		return out;
	}
	
		
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	
	// 'J' is the number of voxels along beams.
	// 'I' is the number of beams.
	// 'L' is the number of targets in each voxels.
	// 'N' is the number of samples in each voxels.
	// 'P' is the frequency of the sound, where the length of sampling intervals is 1.
	// 'w' is the length of the scattered pulse in lengths of voxels.
	// 'ci' is the number of beam affecting the current beam on each side of the beam.
	// 'olp' is a vector of length I * (2*ci+1) giving the degree of overlap between beams.
	// 'X' is the vector to be outputted.
	// 'wM' is the minimum distance between maxima in the superimposed wave, used for extracting only the proper peaks of the wave, and not all the minor peaks (if any).
	// 'seed' is the seed used in the simulation.
	
	
	void rexp_MultSines_VaribleCorrelation(int *J, int *I, int *L, int *N, double *P, double *w, int *ci, double *olp, double *X, int *wM, int *seed)
	{
		
		// Set the random seed:
		srand(*seed);
		
		// Declare the constants and variables used in the calculations:
		// 'cw' is the number of voxels from which sine waves affect the current voxel, found by rounding up 'w':
		int cw = ceil(*w);
		// 'W' is the number of voxels of which sine waves are considered, one more than 'cw':
		int W = cw + 1;
		// 'fullI' is the number of beams used in the calculations ('ci' on each side of every beam):
		int fullI = *I + 2 * *ci;
		// 'ncor' is the number of rows of 'olp', corresponding to the number of beams summed over in the function:
		int ncor = 2 * *ci + 1;
		// 'lX' is the length of the output vector (including the 'cw' first voxels of each beam):
		int lX = *J * *I;
		// 'lS' is the length of the superimposed sine vector, new for each voxel:
		int lS = *N * *I;
		// The number of targets cosidered for each voxel:
		int WL = W * *L;
		// 'lphis' is the total number of targets considered for each sample point (the number of targets considered for each sample point 'WL' * the number of beams considered 'fullI'):
		int lphis = WL * fullI;
		// 'add' is the value to add when generating uniformly distributed phases:
		double add = 0.0;
		// 'npos' is the current position in the curren voxel (see below):
		double npos = 0.0;
		// Reserve in memory the last maximum position of the superimposed sine wave, which is used as when applying the criterion of not to closely positioned maxima:
		int lastMaxPos = 0;
		// Reserve in memory the last maximum value:
		double lastMax = 0.0;
		// Reserve in memory the position of the currently located closest maximum position:
		int MaxPosClosest = 0;
		// Reserve in memory the position of the previously located closest maximum position:
		int lastMaxPosClosest = 0;
		// Store the value og 'P' to denot the midpoint of voxels:
		int midP = ceil(*N/2);
		
		// Modify the frequency P:
		const double PI = 3.141592;
		*P = *P * 2 * PI;
		
		// Initialize the output vector [J,I]:
		for(int i = 0; i < lX; i++)
		{
			X[i] = 0.0;
		}
		// Initialize the vector of the summed echo for each voxel [N,I]:
		double S[lS];
		for(int i = 0; i < lS; i++)
		{
			S[i] = 0.0;
		}
		// Initialize the vector of phases [W*L,I]:
		double phis[lphis];
		for(int i = 0; i < lphis; i++)
		{
			phis[i] = -1000.0;
		}
		// Initialize the vector of sine values [N,J]:
		double s[lphis];
		for(int i = 0; i < lphis; i++)
		{
			s[i] = 0.0;
		}
		
		
		// First sampling interval (W voxels):	
		// Generate the phases:
		for(int l = 0; l < W * *L; l++)
		{
			for(int i = 0; i < fullI; i++)
			{
				add = floor(l / *L) - W + 1;
				phis[rexp_MultSines_VaribleCorrelation_ind2d(l, i, W * *L)] = rand() / double(RAND_MAX) + add;
			}
		}
		
		
		// Loop through the sampling points of the first voxel:
		for(int n = 0; n < *N; n++)
		{
			npos = (double) (n+1) / *N;
			// Calculate the sine values:
			for(int p = 0; p < lphis; p++)
			{
				s[p] = rexp_MultSines_VaribleCorrelation_sin0(npos - phis[p], *P, *w);
			}
			
			// Loop through the beams:
			for(int i = 0; i < *I; i++)
			{
				// Loop through the neighborhood of beams:
				for(int ii = 0; ii < ncor; ii++)
				{
					// Only add the contribution from the targets in the beam if the value of 'olp' is positive:
					if(olp[rexp_MultSines_VaribleCorrelation_ind3d(ii, 0, i, ncor, *J)]>0)
					{
						// Loop through the 'WL' relevant targets:
						for(int l = 0; l < WL; l++)
						{
							S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)] += s[rexp_MultSines_VaribleCorrelation_ind2d(l, i+ii, WL)] * olp[rexp_MultSines_VaribleCorrelation_ind3d(ii, 0, i, ncor, *J)];
						}
					}
				}
			}
		}
		
		// Extract the local maxima in each voxel of each beam of sample interval j=0;
		int j = 0;
		for(int i = 0; i < *I; i++)
		{
			// Reinitialize:
			lastMaxPos = 0;
			lastMax = 0.0;
			MaxPosClosest = 0;
			lastMaxPosClosest = 0;
			
			// Move through the voxel:
			for(int n = 1; n < *N - 1; n++)
			{
				// Is the current position a local maximum (larger than the neighbors)?:
				if( S[rexp_MultSines_VaribleCorrelation_ind2d(n-1, i, *N)]<S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)] && S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)]>S[rexp_MultSines_VaribleCorrelation_ind2d(n+1, i, *N)] )
				{
					// If the current position is a local maximum, and the position differs sufficiently from the last position of maximum, and in addition the new position of maximum is closer to the mid point than the previous, accept this position as a peak but store the old position in case the new position was a false peak, and the proper peak is found to be farther away from the mid point than the pervious:
					if(n-lastMaxPos > *wM)
					{
						// If the position of the currently located maximum is closer to the mid point of the voxel than the previously located maximum, assign the current maximum to 'x':
						if(abs(n-midP) < abs(MaxPosClosest-midP))
						{
							// Assign the currently located maximum to 'x':
							X[rexp_MultSines_VaribleCorrelation_ind2d(j, i, *J)] = pow(S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)],2);
							// Update the last closest maximum position 'lastMaxPosClosest' before updating the closest maximum position 'MaxPosClosest':
							lastMaxPosClosest = MaxPosClosest;
							MaxPosClosest = n;
							lastMaxPos = n;
						}
						// Else there cannot be any closer local maxima than the current:
						else
						{
							break;
						}
					}
					// If the local maximum is too close to the previous, but larger, it is a candidate for the closest local maximum:
					else if(S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)] > lastMax)
					{
						// If the currently located maximum is closer to the mid point of the voxel than the previously located closest maximum at position 'lastMaxPosClosest', the current value is assigned to 'x', and the positions of the closest maximum position 'MaxPosClosest' and the maximum position 'lastMaxPos' are updated:
						if(abs(n-midP) < abs(lastMaxPosClosest-midP))
						{
							X[rexp_MultSines_VaribleCorrelation_ind2d(j, i, *J)] = pow(S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)],2);
							MaxPosClosest = n;
							lastMaxPos = n;
						}
						// Else if the currently located maximum is farther away from the mid point of the voxel than the previously located closest maximum at position 'lastMaxPosClosest', the old closest position is reassignied to 'x', and no maxima can be closer so the loop is breaked:
						else
						{
							X[rexp_MultSines_VaribleCorrelation_ind2d(j, i, *J)] = pow(S[rexp_MultSines_VaribleCorrelation_ind2d(lastMaxPosClosest, i, *N)],2);
							break;
						}
					}
				} // End of if
			} // End of for n
		}
		
		
		// Move on with the rest of the voxels:
		for(int j = 1; j < *J; j++)
		{
			// Clean up the array 'S':
			double S[lS];
			for(int i = 0; i < lS; i++)
			{
				S[i] = 0.0;
			}
			
			// Shift the positions of the phases and insert new phases on the last L positions:
			for(int l = 0; l < (W-1) * *L; l++)
			{
				for(int i = 0; i < fullI; i++)
				{
					phis[rexp_MultSines_VaribleCorrelation_ind2d(l, i, W * *L)] = phis[rexp_MultSines_VaribleCorrelation_ind2d(l + *L, i, W * *L)];
				}
			}
			for(int l = (W-1) * *L; l < W * *L; l++)
			{
				for(int i = 0; i < fullI; i++)
				{
					phis[rexp_MultSines_VaribleCorrelation_ind2d(l, i, W * *L)] = rand() / double(RAND_MAX) + j;
				}
			}
			
			
			// Loop through the sampling points of the j'th voxel:
			for(int n = 0; n < *N; n++)
			{
				npos = (double) (n+1) / *N + j;
				// Calculate the sine values:
				for(int p = 0; p < lphis; p++)
				{
					s[p] = rexp_MultSines_VaribleCorrelation_sin0(npos - phis[p], *P, *w);
				}
				// Loop through the beams:
				for(int i = 0; i < *I; i++)
				{
					// Loop through the neighborhood of beams:
					for(int ii = 0; ii < ncor; ii++)
					{
						// Only add the contribution from the targets in the beam if the value of 'olp' is positive:
						if(olp[rexp_MultSines_VaribleCorrelation_ind3d(ii, j, i, ncor, *J)]>0)
						{
							// Loop through the 'WL' relevant targets:
							for(int l = 0; l < WL; l++)
							{
								S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)] += s[rexp_MultSines_VaribleCorrelation_ind2d(l, i+ii, WL)] * olp[rexp_MultSines_VaribleCorrelation_ind3d(ii, j, i, ncor, *J)];
							}
						}
					}
				}
			}
			
			// Extract the local maxima in each voxel of each beam
			for(int i = 0; i < *I; i++)
			{
				// Reinitialize:
				lastMaxPos = 0;
				lastMax = 0.0;
				MaxPosClosest = 0;
				lastMaxPosClosest = 0;
				
				// Move through the voxel:
				for(int n = 1; n < *N - 1; n++)
				{
					// Is the current position a local maximum (larger than the neighbors)?:
					if( S[rexp_MultSines_VaribleCorrelation_ind2d(n-1, i, *N)]<S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)] && S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)]>S[rexp_MultSines_VaribleCorrelation_ind2d(n+1, i, *N)] )
					{
						// If the current position is a local maximum, and the position differs sufficiently from the last position of maximum, and in addition the new position of maximum is closer to the mid point than the previous, accept this position as a peak but store the old position in case the new position was a false peak, and the proper peak is found to be farther away from the mid point than the pervious:
						if(n-lastMaxPos > *wM)
						{
							// If the position of the currently located maximum is closer to the mid point of the voxel than the previously located maximum, assign the current maximum to 'x':
							if(abs(n-midP) < abs(MaxPosClosest-midP))
							{
								// Assign the currently located maximum to 'x':
								X[rexp_MultSines_VaribleCorrelation_ind2d(j, i, *J)] = pow(S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)],2);
								// Update the last closest maximum position 'lastMaxPosClosest' before updating the closest maximum position 'MaxPosClosest':
								lastMaxPosClosest = MaxPosClosest;
								MaxPosClosest = n;
								lastMaxPos = n;
							}
							// Else there cannot be any closer local maxima than the current:
							else
							{
								break;
							}
						}
						// If the local maximum is too close to the previous, but larger, it is a candidate for the closest local maximum:
						else if(S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)] > lastMax)
						{
							// If the currently located maximum is closer to the mid point of the voxel than the previously located closest maximum at position 'lastMaxPosClosest', the current value is assigned to 'x', and the positions of the closest maximum position 'MaxPosClosest' and the maximum position 'lastMaxPos' are updated:
							if(abs(n-midP) < abs(lastMaxPosClosest-midP))
							{
								X[rexp_MultSines_VaribleCorrelation_ind2d(j, i, *J)] = pow(S[rexp_MultSines_VaribleCorrelation_ind2d(n, i, *N)],2);
								MaxPosClosest = n;
								lastMaxPos = n;
							}
							// Else if the currently located maximum is farther away from the mid point of the voxel than the previously located closest maximum at position 'lastMaxPosClosest', the old closest position is reassignied to 'x', and no maxima can be closer so the loop is breaked:
							else
							{
								X[rexp_MultSines_VaribleCorrelation_ind2d(j, i, *J)] = pow(S[rexp_MultSines_VaribleCorrelation_ind2d(lastMaxPosClosest, i, *N)],2);
								break;
							}
						}
					} // End of if
				} // End of for n
			} // End of for i
		} // End of for j
	} // End of void
} // End of extern "C"
