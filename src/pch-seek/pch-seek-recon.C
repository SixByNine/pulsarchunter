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
	float phase,amp,re,peak,pulsep;
	double *ph_arr, *am_arr, *idx;
	double p;
	int b;
	int m;
	float rms,mean;
	double prevphase=-10000;
	ph_arr= (double*)malloc(sizeof(double)*foldval);
	am_arr= (double*)malloc(sizeof(double)*foldval);

/*
	switch(foldval){
	case 1:
			mean=0;
			rms=1;
			break;
	case 2:
			mean=0.956;
			rms=0.2779;
			break;
	case 4:
			mean=1.08;
			rms=0.2215;
			break;
	case 8:
			mean=1.1837;
			rms=0.1803;
			break;
	case 16:
			mean=1.2715;
			rms=0.1493;
			break;
	}
*/
	samp = (int)(freq/xscale+0.5); // the bin in the hfolded fold

	peak=0;
	for (int f=1; f <= foldval; f++){
		bin=(int)((samp*f)/(float)foldval + 0.5);
		phase = phases[bin];
		amp=amplitudes[bin];
		ph_arr[f-1] = phase;
		am_arr[f-1] = amp;
	//	printf("%d\t%f\t%f\t%f\t%f\n",f,1000.0/freq,phase,amp,spectral_snr);
	}
	int nphases=128;
	for (int b = 0 ; b < nphases; b++){
		re=0;
		pulsep= 2.0*PI*(float)b / (float)nphases;
		for (int f=1; f <= foldval; f++){
			p=ph_arr[f-1]+f*pulsep;
			re+=am_arr[f-1]*cos(p);
		}
//		printf("%f\t%f\t%f\t%f\t%f\t%d\n",pulsep,re/sqrt(foldval),am_arr[0],1000.0/freq,spectral_snr,foldval);
		if(re > peak)peak=re;
	}
	return peak/sqrt(foldval);
//	else return (sqrt(peak/sqrt(foldval))-mean)/rms;
}
