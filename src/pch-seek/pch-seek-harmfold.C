#include "pch-seek.h"

#include <stdlib.h>
#include <math.h>


/**
 * This is a 'stupid' harmonic folder that does each harmonic fold independantly.
 */
float** pch_seek_harmfold_simple(float* amplitudes, int ndat, int* folds, int nfolds){
	int ifold, foldval,f;
	int samp;
	float** hfolds;
	float fac;

	hfolds = (float**) malloc(sizeof(float*) * nfolds);
	

	for (ifold=0; ifold < nfolds; ifold++){
		foldval = folds[ifold];
		hfolds[ifold] = (float*) calloc(ndat,sizeof(float));
		for (samp = 0; samp < ndat; samp++){
			for (f=1; f <= foldval; f++){
				hfolds[ifold][samp] += amplitudes[(int)((samp*f)/foldval)];
			}
		}
		
	}

	return hfolds;
}
