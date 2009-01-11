#include "pch-seek.h"

#include "TKfit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.14159265

float pch_seek_recon_add(float* amplitudes, float* phases, int ndat, int foldval, float freq, float xscale){

	int samp,bin;
	float phase,amp,re,im,s_re,s_im;
	double *ph_arr, *am_arr, *idx;
	double p[2];
	int m;
	float bac=0;
	double prevphase=-10000;


	ph_arr= (double*)malloc(sizeof(double)*foldval);
	am_arr= (double*)malloc(sizeof(double)*foldval);
	idx= (double*)malloc(sizeof(double)*foldval);



	samp = (int)(freq/xscale+0.5); // the bin in the hfolded fold

	s_re=0;
	s_im=0;
	for (int f=1; f <= foldval; f++){
		bin=(int)((samp*f)/foldval);
		phase = phases[bin];
		amp=amplitudes[bin];
		// if the amplitude is less than 0
		// then it can cause a +ve increase 
		// if the phase is correct.
		// I have chosen to subtract all -ve
		// amplitudes from the final amplitude
		if (amp < 0){
			bac+=amp;
			amp=0;
		}
		while (phase < prevphase) phase+=2*PI;
		ph_arr[f-1] = phase;
		am_arr[f-1] = amp;
		idx[f-1]=f;
		prevphase=phase;
	}

	// fit a straight line to the data
	m=1;
	TKleastSquares_svd_noErr(idx,ph_arr,foldval, p, m, TKfitPolyOrigin);

	// subtract the straight line
	
	for (int f=1; f <= foldval; f++){
		ph_arr[f-1]-=(f*p[0]);
	}

	for (int f=1; f <= foldval; f++){
		re=am_arr[f-1]*cos(ph_arr[f-1]);
		im=am_arr[f-1]*sin(ph_arr[f-1]);
		s_re+=re;
		s_im+=im;
	}

	free(ph_arr);
	free(am_arr);
	free(idx);

	return (sqrt(s_re*s_re + s_im*s_im)+bac) / sqrt(foldval);



}
