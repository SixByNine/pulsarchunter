#include "pch-seek.h"

#include "TKfit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#define PI 3.14159265

float pch_seek_recon_ralph(float* amplitudes, float* phases, int ndat, int foldval, float freq, float xscale){

	float bin0;
	float peak;
	int harmmax,i,ifun;
	fftwf_complex* fft;
	fftwf_plan plan;

	peak = -10000.0f;
	bin0 = ((freq/xscale)/foldval);

	harmmax = ndat/(int)(bin0);
	
	if (foldval > harmmax) foldval = harmmax;

	fft = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*(foldval*2+1));

	fft[0][0] = 0;
	fft[0][1] = 0;
	for (i=0; i < foldval; i++){
		ifun = (int)((bin0*(i+1)+0.5));

		fft[i+1][0] = amplitudes[ifun]*cos(phases[ifun]);
		fft[i+1][1] = amplitudes[ifun]*sin(phases[ifun]);
		// make the mirror the congrugate to make the transform real.
		fft[(foldval*2-i)][0] = fft[i+1][0];
		fft[(foldval*2-i)][1] = -fft[i+1][1];

	}
	plan = fftwf_plan_dft_1d(foldval*2+1, fft, fft, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftwf_execute(plan);


	for (i=0; i < foldval*2; i++){
	//	printf("%f %f\n",fft[i][0],fft[i][1]);
		if(fft[i][0] > peak)peak=fft[i][0];
	}
	//printf("\n");




	fftwf_destroy_plan(plan);


	fftwf_free(fft);

	return peak/sqrt(foldval);

}

float pch_seek_recon_add(float* amplitudes, float* phases, int ndat, int foldval, float freq, float xscale, float spectral_snr){

	int samp,bin;
	float phase,amp,re,im,s_re,s_im,sum;
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
	sum=0;
	for (int f=1; f <= foldval; f++){
		bin=(int)((samp*f)/(float)foldval + 0.5);
		phase = phases[bin];
		amp=amplitudes[bin];
		sum+=amp;
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

	float ret = ((sqrt(s_re*s_re + s_im*s_im)+bac) / sqrt(foldval));
	if(ret < 0) ret=0;
	return ret;
//	return spectral_snr*((sqrt(s_re*s_re + s_im*s_im)) / foldval);



}
