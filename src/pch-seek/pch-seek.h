#ifndef PCHFOLD_H_
#define PCHFOLD_H_
#include <psrxml.h>
#include <fftw3.h>


typedef struct pch_seek_operations{
	char dump_tim;
	char fft_input;
	char form_amplitudes;
	char twiddle_amplitudes;
	char dump_amplitudes;
	char dump_phases;
	char dump_normalised;
	char dump_harmfolds;
	char phase_fit;
	char fscrunch;
	char normalise_median;
	char normalise_agl;
	char normalise_powerlaw;
	char harmfold_simple;
	char search_amplitudes;
	char recon_add;
	char write_prd;

	// some options
	float* dmtrials;
	int ndm;
	int* harmfolds;
	int nharms;
	char prdfile[256];
	float amp_thresh;
} pch_seek_operations_t;


/* pch-seek-do-search.C */
void pch_seek_do_search(pch_seek_operations_t* operations, psrxml* header, float** data_arr);

/* pch-seek-read-file.C */
float** pch_seek_read_file(psrxml* header);
void pch_seek_init_operations(pch_seek_operations_t* ops);

/* pch-seek-output.C */
void pch_seek_dump(float* data, int len, float xscale, char* filename);
void pch_seek_write_prd(char* filename, float** freq, float** spec, float** recon, int* ncand, int* harms, int nharm, psrxml* header);

/* pch-seek-fourier.C */
void pch_seek_fourier_r2c(float* timeseries, int npoints);
void pch_seek_form_phase_amp(fftwf_complex* complex_spec, float* amp_spec, float* phase_spec, int nsamp,char twiddle);

/* pch-seek-phase-fit.C */
void pch_seek_phase_fit_simple(float** phases, int nsamp, int nchan, float* dmtrials, float ndm, float fch1, float foff, float tobs);

/* pch-seek-normalise.C */
void pch_seek_normalise_median(float* amplitudes, int ndat, int nrun);
void pch_seek_normalise_agl_mean(float* amplitudes, int ndat, int nrun);
void pch_seek_normalise_powerlaw(float* amplitudes, int ndat);

/* pch-seek-harmfold.C */
float** pch_seek_harmfold_simple(float* amplitudes, int ndat, int* folds, int nfolds);


/* pch-seek-search-spec.C */
float** pch_seek_search_spectrum(float* amplitudes, int ndat, float xscale, float threshold, int* ncand, float rms);

/* pch-seek-recon.C */
float pch_seek_recon_add(float* amplitudes, float* phases, int ndat, int foldval, float freq, float xscale);



#endif /*PCHFOLD_H_*/

