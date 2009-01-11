#include <stdlib.h>
#include <complex>

#include <fftw3.h>

#include "pch-seek.h"

/**
 *
 * Do an in place fourier transform on the input
 * real time series.
 *
 * produces a complex spectrum.
 *
 * Input and output pointers are the same.
 *
 * performs a simple fftw 32_bit r2c transform.
 */
void pch_seek_fourier_r2c(float* timeseries, int npoints){
	fftwf_complex* output;
	fftwf_plan transform_plan;

	output = (fftwf_complex*)timeseries;

	transform_plan = fftwf_plan_dft_r2c_1d(npoints, timeseries, output, FFTW_ESTIMATE);
	// hopefully we have wisdom already, so we will be better than ESTIMATE!

	fftwf_execute(transform_plan);
	
	fftwf_destroy_plan(transform_plan);

	return;
}
	

void pch_seek_form_phase_amp(fftwf_complex* complex_spec, float* amp_spec, float* phase_spec, int ncomplex, char twiddle){
	int i;
	std::complex<float>* cmp;
	cmp = reinterpret_cast<std::complex<float>* >(complex_spec);
	if(twiddle){
		// this is the strange way of computing the amplitudes
		// taken from seek (formspec.f)
		float ar,ai,arl,ail,a1,a2;
		for (i=0; i < ncomplex; i++){
			ar=real(cmp[i]);
			ai=imag(cmp[i]);
			a1=ar*ar+ai*ai;
			arl=ar-arl;
			ail=ai-ail;
			a2=(arl*arl+ail*ail)/2.0;
			if(a1>a2)amp_spec[i]=sqrt(a1);
			else amp_spec[i]=sqrt(a2);
			arl=ar;
			ail=ai;
			phase_spec[i] = arg(cmp[i]);

		}
	}else{
		for (i=0; i < ncomplex; i++){
			amp_spec[i] = abs(cmp[i]);
			phase_spec[i] = arg(cmp[i]);
		}
	}
}


