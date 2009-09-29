#include "pch-seek.h"

#include <stdlib.h>
#include <math.h>

void pch_seek_harmfold_stupid(float* amplitudes, int ndat, int* folds, int ifolds, float** hfolds);

/**
 * This is a 'stupid' harmonic folder that does each harmonic fold independantly.
 */
float** pch_seek_harmfold_simple(float* amplitudes, int ndat, int* folds, int nfolds){
	int ifold;
	float** hfolds;

	hfolds = (float**) malloc(sizeof(float*) * nfolds);

	// Loop over each harmonic fold
	for (ifold=0; ifold < nfolds; ifold++){
		// Make the array that stores the output
		hfolds[ifold] = (float*) calloc(ndat,sizeof(float));
		if(hfolds[ifold] == NULL){
			printf("Could not allocate enough memory for harmfold %d\n",folds[ifold]);
			exit(1);
		}
		pch_seek_harmfold_stupid(amplitudes,ndat,folds,ifold,hfolds);
	}
	return hfolds;
}

float** pch_seek_harmfold_smart(float* amplitudes, int ndat, int* folds, int nfolds){
	int ifold,max_clever_folds;
	float** hfolds;

	hfolds = (float**) malloc(sizeof(float*) * nfolds);
	// Make the array that stores the output
	for (ifold=0; ifold < nfolds; ifold++){
		hfolds[ifold] = (float*) calloc(ndat,sizeof(float));
		if(hfolds[ifold] == NULL){
			fprintf(stderr,"Could not allocate enough memory for harmfold %d\n",folds[ifold]);
			exit(1);
		}
	}
		

	// See if we have 2,4,8,16,... etc
	for (ifold=0; ifold < nfolds; ifold++){
		if (folds[ifold] != (int)pow(2,ifold+1) || ifold > 6)break;
	}
	max_clever_folds=ifold;
	if(max_clever_folds > 0){
		// if we did have powers of 2, we can do the 'clever' folding
		// of AGL et al.
		//printf("Smart_folding < %d\n",max_clever_folds);
		pch_seek_harmfold_pow_two(amplitudes,ndat,hfolds,max_clever_folds);
	}

	// Any left over folds, we have to do the 'stupid' way
	for (ifold=max_clever_folds; ifold < nfolds; ifold++){
		//printf("Stupid_folding %d\n",ifold);
		pch_seek_harmfold_stupid(amplitudes,ndat,folds,ifold,hfolds);
	}
	return hfolds;
}




void pch_seek_harmfold_stupid(float* amplitudes, int ndat, int* folds, int ifold, float** hfolds){
	int foldval,f;
	int samp,idx,last_idx;
	float fac;

	foldval = folds[ifold];
// Loop over each sample in output
	for (samp = 0; samp < ndat; samp++){
		last_idx=0;
		for (f=1; f <= foldval; f++){
			// add samp/f 2samp/f 3samp/f 4samp/f ... fsamp/f
			idx=(int)((samp*f)/(float)foldval + 0.5);
			if (idx > 1 && idx != last_idx){
				hfolds[ifold][samp] +=
					amplitudes[idx];
			}
			last_idx=idx;
		}
	}
}
