#include "pch-seek.h"

#include "toolkit.h"
#include <math.h>

#include <stdlib.h>

float** pch_seek_search_spectrum(float* amplitudes, int ndat, float xscale, float threshold, int* ncand, float rms,float xoff){
	
	float** results;
	int arrsize;
	int nres;

	results = (float**) malloc(sizeof(float*) * 2);

	arrsize = 128;
	results[0] = (float*) malloc(sizeof(float)*arrsize); // freq
	results[1] = (float*) malloc(sizeof(float)*arrsize); // SNR

	
	nres=0;
	for (int i = 0; i < ndat; i++){
		if(xscale*(float)i < 0.1) continue;
		if (amplitudes[i]/rms > threshold){
			results[0][nres] = xscale*(float)i+xoff;
			results[1][nres] = amplitudes[i]/rms;
			nres++;
			if(nres >= arrsize){
				arrsize*=2;
				results[0] = (float*) realloc(results[0],sizeof(float)*arrsize);
				results[1] = (float*) realloc(results[1],sizeof(float)*arrsize);
			}
		}
	}

	(*ncand) = nres;
	return results;
}
