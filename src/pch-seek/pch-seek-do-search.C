#include "pch-seek.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

bool pch_seek_sanity_check(pch_seek_operations_t* operations, psrxml* header);
void pch_seek_search_flat_spec(pch_seek_operations_t* operations, psrxml* header, float* amplitude_fscrunch, float* phase_spectum, float xoff);

/**
 *
 * Calls the appropriate operations in the correct sequence
 *
 * Checks that the required operations are valid.
 *
 * On success, the operation is left set to 1.
 * On failure, the operation is set to -1.
 *
 * This function should do as little work as possible, with the hard work
 * done inside other subroutines.
 *
 * To save memory, it may overwrite the same memory area many times.
 *
 * Therefore there may be many pointers that point to the same memory
 * area, although they should be set to NULL once the data is no longer
 * valid for that pointer.
 * 
 * More or less the order is:
 *
 * time_arr -> complex_spectra -> amplitude_spectra/phase_spectra
 *
 * For the amplitude search: 
 * amplitdude_spectra -> amplitude_fscrunch (+ amplitude_harmfolds)
 * suspects are stored in amplitude_harmfold_freqs, amplitude_harmfold_spectral, amplitude_harmfold_recon
 *
 */
void pch_seek_do_search(pch_seek_operations_t* operations, psrxml* header, float** time_arr){
	int samp,chan,nchan,nsamp,ncomplex,ifold,icand;
	char filename[1024];
	char remeber_filename[1024];
	fftwf_complex** complex_spectra; // the real-imag pairs
	float** amplitude_spectra; // the amplitudes
	float** phase_spectra; // the phases
	float*  amplitude_fscrunch; // the freq flattened phases
	float** amplitude_harmfolds; // the harmonic folds
	float** amplitude_harmfold_freqs; // the best frequencies
	float** amplitude_harmfold_spectral; // the best snrs
	float** amplitude_harmfold_recon; // the best recon snrs
	float xoff; // we might offset the x (freq) scale if we twiddle amps, for instance.
	float tobs;

	tobs=header->actualObsTime;
	nchan=header->numberOfChannels;
	nsamp=header->numberOfSamples;
	ncomplex=nsamp/2;
	amplitude_fscrunch=NULL;
	amplitude_harmfolds=NULL;
	amplitude_spectra = (float**)malloc(sizeof(float*)*nchan);
	phase_spectra     = (float**)malloc(sizeof(float*)*nchan);
	complex_spectra   = (fftwf_complex**)malloc(sizeof(fftwf_complex*)*nchan);
	xoff=0;
	amplitude_harmfold_freqs=amplitude_harmfold_spectral=amplitude_harmfold_recon=NULL;
	for (chan=0; chan < nchan; chan++){
		amplitude_spectra[chan] = NULL;
		phase_spectra[chan] = NULL;
		complex_spectra[chan] = NULL;
	}




	if (pch_seek_sanity_check(operations,header)){
		fprintf(stderr,"Sanity check failed :(\nCannot continue.\n");
	} else {
		/*
		 * Now do the work.
		 *
		 * Perhaps it would be better code if this wasnt all inside
		 * a giant 'else' statement.
		 *
		 * Perhaps it can be made nicer later.
		 *
		 * This part of the method is quiet long, but it doesn't do any
		 * actual work. All the operations are in external functions
		 * and the operations to do are determined by 'operations'
		 *
		 */
		printf("Operations:\n");
		if(operations->dump_tim)printf(" - Dump timeseries\n");
		if(operations->hist_tim)printf(" - Histogram timeseries\n");
		if(operations->fft_input)printf(" - Fourier transform\n");
		if(operations->write_presto_fft)printf(" - Write raw spectrum to presto file ('%s.fft')\n",operations->presto_fft_file);
		if(operations->form_amplitudes){
			printf(" - Form phases/amplitudes");
			if(operations->twiddle_amplitudes)printf(" (twiddle mode)\n");
			else printf(" (natural mode)\n");
		}
		if(operations->dump_amplitudes)printf(" - Dump amplitudes\n");
		if(operations->hist_amplitudes)printf(" - Histogram amplitudes\n");
		if(operations->giant_search)printf(" - Single Pulse Search (GSearch)\n");
		if(operations->dump_phases)printf(" - Dump phases\n");
		if(operations->phase_fit)printf(" - Phase fit search\n");
		if(operations->search_chans)printf(" - Search all %d channels in turn\n",nchan);
		if(operations->fscrunch && nchan > 1)printf(" - Flatten in frequency\n");
		if(operations->normalise_median)printf(" - Normalising amplitudes (median method)\n");
		if(operations->normalise_agl)printf(" - Normalising amplitudes (Lyne et al. mean/rms method)\n");
		if(operations->normalise_powerlaw)printf(" - Normalising amplitudes (powerlaw fitting method)\n");
		if(operations->dump_normalised)printf(" - Dump normalised amplitudes\n");
		if(operations->hist_normalised)printf(" - Histogram normalised amplitudes\n");
		if(operations->harmfold_smart && operations->nharms > 0)printf(" - Harmonicly folding %d times (SMART method)\n",operations->nharms);
		if(operations->harmfold_simple && operations->nharms > 0)printf(" - Harmonicly folding %d times (simple method)\n",operations->nharms);
		if(operations->dump_harmfolds)printf(" - Dump harmonicaly folded amplitudes\n");
		if(operations->hist_harmfolds)printf(" - Histogram harmonicaly folded amplitudes\n");
		if(operations->search_amplitudes)printf(" - Search harmonicaly folded amplitudes\n");
		if(operations->recon_add)printf(" - Compute recon SNR (addition method)\n");

		if(operations->write_prd)printf(" - Write prd file '%s'\n",operations->prdfile);




		if (operations->dump_tim){
			/* The 'Dump Timeseries' operation acts before any other
			 * operations as it just ascii dumps the data as loaded.
			 */
			printf("OP[START]: Dump Timeseries\n");
			for (chan=0; chan < nchan; chan++){
				sprintf(filename,"timeseries_ch%03d.ascii",chan);
				pch_seek_dump(time_arr[chan],header->numberOfSamples,header->currentSampleInterval, filename);
			}
			printf("OP[END]: Dump Timeseries\n");
		}
		if (operations->hist_tim){
			printf("OP[START]: Histogram Timeseries\n");
			for (chan=0; chan < nchan; chan++){
				sprintf(filename,"timeseries_ch%03d.hist",chan);
				pch_seek_histogram(time_arr[chan],header->numberOfSamples,256, filename);
			}
			printf("OP[END]: Histogram Timeseries\n");
		}


		if (operations->giant_search){
			/* 
			*/
			printf("OP[START]: Single Pulse Search (GSearch) file: '%s'\n",operations->giantfile);
			for (chan=0; chan < nchan; chan++){
				strcpy(filename,operations->giantfile);
				for (chan=0; chan < nchan; chan++){
					if(nchan > 1) sprintf(filename,"%s_ch%03d",operations->giantfile,chan);
					pch_seek_singlepulse(time_arr[chan],header->numberOfSamples,5.0,header->referenceDm,8,header->receiverBeamNumber, filename);
				}

			}
			printf("OP[END]: Single Pulse Search (GSearch)\n");
		}



		if (operations->fft_input){
			/* We just fft the input data for later work.
			*/
			printf("OP[START]: FFT Timeseries\n");
			for (chan=0; chan < header->numberOfChannels; chan++){
				pch_seek_fourier_r2c(time_arr[chan],header->numberOfSamples);
				complex_spectra[chan] = (fftwf_complex*)(time_arr[chan]);
				time_arr[chan] = NULL;
			}
			printf("OP[END]: FFT Timeseries\n");

		}
		if (operations->write_presto_fft){
			/* This operation writes out the fft'd data in a form that
			 * can be read in by S.Ransom's presto software.
			 */
			printf("OP[START]: Write raw spectrum to presto file\n");
			if(nchan > 1){
				// If we have multiple channels, we have to write out many files as
				// I don't think presto supports multi-channel spectra
				for (chan=0; chan < nchan; chan++){
					sprintf(filename,"%s%03d",operations->presto_fft_file,chan);
					pch_seek_write_presto_fft(filename,complex_spectra[chan],ncomplex,header);
				}
			} else {
				// otherwise just call the write function on the data.
				pch_seek_write_presto_fft(operations->presto_fft_file,complex_spectra[0],ncomplex,header);
			}
			printf("OP[END]: Write raw spectrum to presto file\n");
		}

		if (operations->form_amplitudes){
			printf("OP[START]: Form phase-amplitude spectra");
			if (operations->twiddle_amplitudes){
				printf(" ('Twiddle' amplitude mode)\n");
				xoff=-0.25/tobs;
			}
			else printf(" (Natural mode)\n");

			// form the phase/amplitude arrays from the real-imag pairs.
			// Uses (mostly) the same arrays
			// we in fact shift over one block of memory so that we dont
			// overwrite the data we wanted to use.
			// Probably very slow compared to in-place, but we will have
			// to see.
			amplitude_spectra[0] = (float*) fftwf_malloc(sizeof(float)*nsamp);
			phase_spectra[0] =  amplitude_spectra[0]+ncomplex;
			pch_seek_form_phase_amp(complex_spectra[0],amplitude_spectra[0],phase_spectra[0],ncomplex,operations->twiddle_amplitudes);
			// remove the bins 'longer' than 20 seconds, as this is likely to be 'bunk'
			int zz=0;
			for(samp=0 ; (samp/tobs) < 0.05; samp++){
				amplitude_spectra[0][samp]=0;
				zz++;
			}
			printf("Removed leading %d bins\n",zz);
			for (chan=1; chan < nchan; chan++){
				amplitude_spectra[chan] = ((float*)complex_spectra[chan-1]);
				phase_spectra[chan] = ((float*)complex_spectra[chan-1]) + ncomplex;
				pch_seek_form_phase_amp(complex_spectra[chan],amplitude_spectra[chan],phase_spectra[chan],ncomplex,operations->twiddle_amplitudes);
				for(samp=0 ; (samp/tobs) < 0.05; samp++){
					amplitude_spectra[0][samp]=0;
				}

			}
			fftwf_free(complex_spectra[nchan-1]);
			// we have overwritten the complex spectra, so set the pointers to null.
			for (chan=0; chan < nchan; chan++){
				complex_spectra[chan]=NULL;
			}
			printf("OP[END]: Form phase-amplitude spectra\n");

		}

		if (operations->hist_amplitudes){
			printf("OP[START]: Histogram amplitudes\n");
			for (chan=0; chan < nchan; chan++){
				sprintf(filename,"amplitudes_ch%03d.hist",chan);
				pch_seek_histogram(amplitude_spectra[chan],ncomplex, 256, filename);
			}
			printf("OP[END]: Histogram amplitudes\n");
		}


		if (operations->dump_amplitudes){
			printf("OP[START]: Dump amplitudes\n");
			for (chan=0; chan < nchan; chan++){
				sprintf(filename,"amplitudes_ch%03d.ascii",chan);
				pch_seek_dump(amplitude_spectra[chan],ncomplex, 1.0/tobs, filename);
			}
			printf("OP[END]: Dump amplitudes\n");
		}
		if (operations->dump_phases){
			printf("OP[START]: Dump phases\n");
			for (chan=0; chan < nchan; chan++){
				sprintf(filename,"phases_ch%03d.ascii",chan);
				pch_seek_dump(phase_spectra[chan],ncomplex, 1.0/tobs, filename);
			}
			printf("OP[END]: Dump phases\n");
		}
		if (operations->phase_fit){

			printf("OP[START]: Phase fit search\n");
			pch_seek_phase_fit_simple(phase_spectra,ncomplex,nchan,operations->dmtrials,operations->ndm,header->centreFreqCh1,header->freqOffset,tobs); 
			printf("OP[END]: phase fit search\n");
		}

		if (operations->search_chans){
			// search each channel for periodic signals
			// currently destroys the contents of the arrays!
                        strcpy(remeber_filename,operations->prdfile);
			for (chan=0; chan < nchan; chan++){
				sprintf(operations->prdfile,"%s.ch%04d",remeber_filename,chan);
				pch_seek_search_flat_spec(operations,header,amplitude_spectra[chan],phase_spectra[chan],xoff);
			}
			strcpy(operations->prdfile,remeber_filename);
		}

		if (operations->fscrunch){
			// fscrunch the amplitudes.
			if(nchan>0)printf("OP[START]: Flattening in frequency\n");

			amplitude_fscrunch=amplitude_spectra[0];
			for (chan=1; chan < nchan; chan++){
				for (samp=0; samp<ncomplex; samp++){
					amplitude_fscrunch[samp]+=amplitude_spectra[chan][samp];
				}
			}
			if(nchan>0)printf("OP[END]: Flattening in frequency\n");
			pch_seek_search_flat_spec(operations,header,amplitude_fscrunch,phase_spectra[0],xoff);
		}
		/*
		 * Clean up the arrays that are left.
		 *
		 */
		for (chan=0; chan < header->numberOfChannels; chan++){
			if (time_arr[chan] != NULL)fftwf_free(time_arr[chan]);
			if (complex_spectra[chan] != NULL)fftwf_free(complex_spectra[chan]);
			if (amplitude_spectra[chan] != NULL)fftwf_free(amplitude_spectra[chan]);
		}

		free(time_arr);
		free(complex_spectra);
		free(amplitude_spectra);
		free(phase_spectra);

	}
}
void pch_seek_search_flat_spec(pch_seek_operations_t* operations, psrxml* header, float* amplitude_fscrunch, float* phase_spectrum,float xoff)
{

	int samp,chan,nchan,nsamp,ncomplex,ifold,icand;
	char filename[1024];
	fftwf_complex** complex_spectra; // the real-imag pairs
	float** amplitude_harmfolds; // the harmonic folds
	float** amplitude_harmfold_freqs; // the best frequencies
	float** amplitude_harmfold_spectral; // the best snrs
	float** amplitude_harmfold_recon; // the best recon snrs
	int* ncand_amp;
	float tobs;

	tobs=header->actualObsTime;
	nchan=header->numberOfChannels;
	nsamp=header->numberOfSamples;
	ncomplex=nsamp/2;
	amplitude_harmfolds=NULL;
	ncand_amp = (int*)malloc(sizeof(int)*(operations->nharms+1));
	amplitude_harmfold_freqs=amplitude_harmfold_spectral=amplitude_harmfold_recon=NULL;

	if (operations->normalise_median){
		// normalise the spectra using the reatough/agl median method
		printf("OP[START]: Normalising using median method\n");
		pch_seek_normalise_median(amplitude_fscrunch,ncomplex,128);
		printf("OP[END]: Normalising using median method\n");
	}

	if (operations->normalise_agl){
		// normalise the spectra using the reatough/agl median method
		printf("OP[START]: Normalising using Lyne et al. mean/rms method\n");
		pch_seek_normalise_agl_mean(amplitude_fscrunch,ncomplex,128);
		printf("OP[END]: Normalising using Lyne et al. mean/rms method\n");
	}

	if (operations->normalise_powerlaw){
		// normalise the spectra using the reatough/agl median method
		printf("OP[START]: Normalising using powerlaw fitting\n");
		pch_seek_normalise_powerlaw(amplitude_fscrunch,ncomplex);
		amplitude_fscrunch[0]=0;
		printf("OP[END]: Normalising using powerlaw fitting\n");
	}



	if (operations->dump_normalised){
		printf("OP[START]: Dump normalised amplitudes\n");
		sprintf(filename,"normalised_amp.ascii",chan);
		pch_seek_dump(amplitude_fscrunch,ncomplex, 1.0/tobs, filename);
		printf("OP[END]: Dump normalised amplitudes\n");
	}
	if (operations->hist_normalised){
		printf("OP[START]: Histogram normalised amplitudes\n");
		sprintf(filename,"normalised_amp.hist",chan);
		pch_seek_histogram(amplitude_fscrunch,ncomplex,256, filename);
		printf("OP[END]: Histogram normalised amplitudes\n");
	}

	if (operations->harmfold_smart && operations->nharms > 0){
		printf("OP[START]: Harmonic folding (smart method) (H =");
		for(int i=0; i <  operations->nharms; i++){
			printf(" %d",operations->harmfolds[i]);
		}
		printf(")\n");
		amplitude_harmfolds = pch_seek_harmfold_smart(amplitude_fscrunch,ncomplex,operations->harmfolds,operations->nharms);
		printf("OP[END]: Harmonic folding (smart method)");
		printf("\n");
	}

	if (operations->harmfold_simple && operations->nharms > 0){
		printf("OP[START]: Harmonic folding (simple method) (H =");
		for(int i=0; i <  operations->nharms; i++){
			printf(" %d",operations->harmfolds[i]);
		}
		printf(")\n");
		amplitude_harmfolds = pch_seek_harmfold_simple(amplitude_fscrunch,ncomplex,operations->harmfolds,operations->nharms);
		printf("OP[END]: Harmonic folding (simple method)");
		printf("\n");
	}

	if (operations->dump_harmfolds){
		printf("OP[START]: Dump harmonicaly folded amplitudes\n");
		for (ifold=0; ifold < operations->nharms; ifold++){
			sprintf(filename,"harmfold_%03d.ascii",ifold);
			pch_seek_dump(amplitude_harmfolds[ifold],ncomplex, 1.0/(operations->harmfolds[ifold]*tobs), filename);
		}
		printf("OP[END]: Dump harmonicaly folded amplitudes\n");
	}
	if (operations->hist_harmfolds){
		printf("OP[START]: Histogram harmonicaly folded amplitudes\n");
		for (ifold=0; ifold < operations->nharms; ifold++){
			sprintf(filename,"harmfold_%03d.hist",ifold);
			pch_seek_histogram(amplitude_harmfolds[ifold],ncomplex, 256, filename);
		}
		printf("OP[END]: Histogram harmonicaly folded amplitudes\n");
	}



	if (operations->search_amplitudes){
		printf("OP[START]: Search amplitudes for SNR > %f\n",operations->amp_thresh);

		amplitude_harmfold_freqs = (float**)malloc(sizeof(float*)*(operations->nharms+1));
		amplitude_harmfold_spectral = (float**)malloc(sizeof(float*)*(operations->nharms+1));

		float rms=0;

		// work out an estimate of the 'actual' rms of the normalised data
		// this method does the same thing as 'seek'
		// perhaps it isn't the best thing to do...
		int n=0;
		for(samp=0; samp < ncomplex; samp+=1024){
			if(fabs(amplitude_fscrunch[samp]) < 3.0){
				rms+= amplitude_fscrunch[samp]*amplitude_fscrunch[samp];
				n++;
			}
		}
		rms = sqrt(rms/(float)n);
		printf("Spectral RMS: %f\n",rms);


		float** res;
		// search the non-folded data
		//float** pch_seek_search_spectrum(float* amplitudes, int ndat, float tobs, float threshold, int* ncand);
		res=pch_seek_search_spectrum(amplitude_fscrunch,ncomplex,1.0/tobs,operations->amp_thresh,&(ncand_amp[0]),rms,xoff);

		amplitude_harmfold_freqs[0]=res[0];
		amplitude_harmfold_spectral[0]=res[1];
		// search the folded data
		for (ifold=0; ifold < operations->nharms; ifold++){
			res=pch_seek_search_spectrum(amplitude_harmfolds[ifold],ncomplex,1.0/(operations->harmfolds[ifold]*tobs),operations->amp_thresh,&(ncand_amp[ifold+1]),rms*sqrt(operations->harmfolds[ifold]),xoff/operations->harmfolds[ifold]);
			// Important note! If we 'twiddled' the amplitudes, our bins are now off by 1/4 of a bin
			amplitude_harmfold_freqs[ifold+1]=res[0];
			amplitude_harmfold_spectral[ifold+1]=res[1];
		}
		printf("OP[END]: Search amplitudes for SNR > %f\n",operations->amp_thresh);
	}

	if(operations->recon_ralph){
		printf("OP[START]: Computing Recon SNR (R.Eatough method)\n");
		amplitude_harmfold_recon = (float**)malloc(sizeof(float*)*(operations->nharms+1));
		amplitude_harmfold_recon[0] = (float*)malloc(sizeof(float)*ncand_amp[0]);
		for (icand = 0; icand < ncand_amp[0]; icand++){
			amplitude_harmfold_recon[0][icand] = pch_seek_recon_ralph(amplitude_fscrunch, phase_spectrum, ncomplex, 1, amplitude_harmfold_freqs[0][icand],1.0/tobs);
		}
		for (ifold=0; ifold < operations->nharms; ifold++){
			amplitude_harmfold_recon[ifold+1] = (float*)malloc(sizeof(float)*ncand_amp[ifold+1]);
			for (icand = 0; icand < ncand_amp[ifold+1]; icand++){
				amplitude_harmfold_recon[ifold+1][icand] = pch_seek_recon_ralph(amplitude_fscrunch, phase_spectrum, ncomplex, operations->harmfolds[ifold], amplitude_harmfold_freqs[ifold+1][icand],1.0/(operations->harmfolds[ifold]*tobs));
			}
		}
		printf("OP[END]: Computing Recon SNR (R.Eatough method)\n");
	}

	if(operations->recon_add){
		printf("OP[START]: Computing Recon SNR (addition method)\n");
		amplitude_harmfold_recon = (float**)malloc(sizeof(float*)*(operations->nharms+1));
		amplitude_harmfold_recon[0] = (float*)malloc(sizeof(float)*ncand_amp[0]);
		for (icand = 0; icand < ncand_amp[0]; icand++){
			amplitude_harmfold_recon[0][icand] = pch_seek_recon_add(amplitude_fscrunch, phase_spectrum, ncomplex, 1, amplitude_harmfold_freqs[0][icand],1.0/tobs,amplitude_harmfold_spectral[0][icand]);
		}
		for (ifold=0; ifold < operations->nharms; ifold++){
			amplitude_harmfold_recon[ifold+1] = (float*)malloc(sizeof(float)*ncand_amp[ifold+1]);
			for (icand = 0; icand < ncand_amp[ifold+1]; icand++){
				amplitude_harmfold_recon[ifold+1][icand] = pch_seek_recon_add(amplitude_fscrunch, phase_spectrum, ncomplex, operations->harmfolds[ifold], amplitude_harmfold_freqs[ifold+1][icand],1.0/(operations->harmfolds[ifold]*tobs),amplitude_harmfold_spectral[ifold+1][icand]);
			}
		}
		printf("OP[END]: Computing Recon SNR (addition method)\n");
	}

	if (operations->nharms > 0 && operations->hfold_bonus_factor != 0.0){
		printf("OP[START]: Adding 'bonus' sqrt(hfold)*%f to SNR for each fold\n",operations->hfold_bonus_factor);
		for (icand = 0; icand < ncand_amp[0]; icand++){
			amplitude_harmfold_spectral[0][icand] += operations->hfold_bonus_factor;
			if (amplitude_harmfold_recon!=NULL)amplitude_harmfold_recon[0][icand] += operations->hfold_bonus_factor;

		}
		for (ifold=0; ifold < operations->nharms; ifold++){
			for (icand = 0; icand < ncand_amp[ifold+1]; icand++){
				amplitude_harmfold_spectral[ifold+1][icand] += sqrt(operations->harmfolds[ifold])*operations->hfold_bonus_factor;
				if (amplitude_harmfold_recon!=NULL)amplitude_harmfold_recon[ifold+1][icand] += operations->hfold_bonus_factor;
			}
		}
		printf("OP[END]: Adding 'bonus' sqrt(hfold)*%f to SNR for each fold\n",operations->hfold_bonus_factor);
	}

	if (operations->write_prd){
		printf("OP[START]: Writing 'prd' file: '%s'",operations->prdfile);
		if(operations->append_output)printf(" (Appending)");
		printf("\n");
		int* mod_harmfolds = (int*)malloc(sizeof(int)*(operations->nharms+1));
		// we have to add the '1' harmonic fold.
		mod_harmfolds[0]=1;
		memcpy(mod_harmfolds+1,operations->harmfolds,sizeof(int)*operations->nharms);

		pch_seek_write_prd(operations->prdfile, amplitude_harmfold_freqs, amplitude_harmfold_spectral,
				amplitude_harmfold_recon, ncand_amp, mod_harmfolds, operations->nharms+1, header,operations->append_output);
		printf("OP[END]: Writing 'prd' file: '%s'\n",operations->prdfile);
	}

	/*
	 *Clean up remaining arrays.
	 */
	if(amplitude_harmfolds!=NULL){
		for(ifold = 0; ifold < operations->nharms; ifold++){
			free(amplitude_harmfolds[ifold]);
		}
		free(amplitude_harmfolds);
	}
	if(amplitude_harmfold_recon!=NULL){
		for(ifold = 0; ifold < operations->nharms+1; ifold++){
			free(amplitude_harmfold_recon[ifold]);
		}
		free(amplitude_harmfold_recon);
	}
	if(amplitude_harmfold_spectral!=NULL){
		for(ifold = 0; ifold < operations->nharms+1; ifold++){
			free(amplitude_harmfold_spectral[ifold]);
			free(amplitude_harmfold_freqs[ifold]);
		}
		free(amplitude_harmfold_spectral);
		free(amplitude_harmfold_freqs);
	}

}


/**
 *
 * Check that the specified options are compatable.
 *
 * This should be improved in future when more options are known.
 *
 *
 */
bool pch_seek_sanity_check(pch_seek_operations_t* operations, psrxml* header){
	/*
	 * Check we can work with this file!
	 */
	if (header->numberOfSamples < 1024){
		// too small!
		fprintf(stderr, "Cannot work with this small number of samples\n");
		return 1;
	}

	
	/*
	 * Check for sane combinations of options
	 */

	if (operations->write_presto_fft){
		operations->fft_input=1;
	}
	if (operations->dump_amplitudes || operations->dump_phases || operations->hist_amplitudes){
		operations->fft_input = 1;
		operations->form_amplitudes=1;
	}

	if (operations->phase_fit){
		if (header->numberOfChannels < 4){
			fprintf(stderr, "Need at least four channels to do a phase fit search\n");
			return 2;
		}
		if (operations->ndm < 1 || operations->dmtrials == NULL){
			fprintf(stderr, "Need at least one DM trial to do a phase fit search\n");
			return 2;
		}
		operations->fft_input = 1;
		operations->form_amplitudes=1;

	}


	if (operations->search_chans){
		operations->fft_input = 1;
		operations->form_amplitudes=1;
		if (operations->fscrunch){
			fprintf(stderr, "Sorry, cannot fscrunch and search all channels :(\n");
			return 3;
		}
	} else {
		operations->fscrunch=1;
	}


	if (operations->dump_normalised || operations->hist_normalised){
		operations->fft_input = 1;
		operations->form_amplitudes=1;
		if(!(operations->normalise_agl || operations->normalise_powerlaw || operations->normalise_median))operations->normalise_powerlaw=1;
	}

	if (operations->dump_harmfolds || operations->hist_harmfolds){
		if(!operations->harmfold_simple)
			operations->harmfold_smart=1;
	}

	if(operations->recon_add || operations->recon_ralph){
		if (header->numberOfChannels > 1){
			fprintf(stderr, "Cannot compute recon SNR with multi-channel data\n");
			return 4;

		}
		if (operations->harmfold_simple)
			operations->harmfold_simple = 1;
		operations->search_amplitudes=1;
	}


	if (operations->write_prd){
		if (operations->harmfold_simple)
			operations->harmfold_simple = 1;
		operations->search_amplitudes=1;
	}

	if (operations->harmfold_simple || operations->harmfold_smart){
		operations->fft_input = 1;
		operations->form_amplitudes=1;
		if(!(operations->normalise_agl || operations->normalise_powerlaw || operations->normalise_median))operations->normalise_powerlaw=1;
	}




	return 0;

}
