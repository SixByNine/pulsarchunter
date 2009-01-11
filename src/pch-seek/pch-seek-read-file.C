#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>

#include "pch-seek.h"

/**
 * float** pch_seek_read_file(psrxml*)
 *
 * returns a 2-dimentional float array with nchans x nsamps data points
 * takes a psrxml header that describes the data to be read
 *
 * Uses the psrxml IO routines to read the data file.
 * Output data are created using fftwf_malloc so are alligned for MMX/SSE SIMD operations
 * i.e. this makes FFTW faster! Note that you should then free with fftw_free.
 *
 * 2008-12-22	M.Keith		Initial version
 *
 */
float** pch_seek_read_file(psrxml* header){

	unsigned char* byte_array;
	float** output;
	float* float_array;
	float** block_chan_array;
	int data_size;
	int chan;
	int start,end,nsamps;
	char swap_chans;
	int nbytes_to_read;
	int read,blocks_read;
	int update_on_block;
	dataFile* file;

	file=header->files[0];
	swap_chans = header->freqOffset > 0;

	fprintf(stdout,"Reading data file '%s'.\n",file->filename);

	//	the byte_array stores one block of byte data
	byte_array = (unsigned char*)malloc(file->blockLength);
	//	float_array stores byte_array in floats, channel by channel (i.e. PFT order)
	float_array = (float*) malloc(file->blockLength*sizeof(float));
	// block_chan_array are pointers to the start of each channel in float_array
	block_chan_array = (float**)malloc(sizeof(float*) * header->numberOfChannels);
	//      the output store the channel time-series
	output = (float**) malloc(sizeof(float*)*header->numberOfChannels);

	// now try and malloc enough memory for the entire file!
	// we malloc for an extra two samples to cover for the later transform to a complex spectrum
	fprintf(stdout,"Trying to malloc %d bytes for the data.\n", header->numberOfChannels*sizeof(float)*(header->numberOfSamples+2));
	for (chan=0; chan < header->numberOfChannels; chan++){
		output[chan] = (float*)malloc(sizeof(float)*(header->numberOfSamples+2));
		if(output[chan]==NULL){
			// we probably ran out of memory...
			fprintf(stderr,"Could not allocate enough memory\n");
			return NULL;
		}
	}


	// prepare for reading!
	readPsrXMLPrepDataFile(file, file->filename);

	start = 0;
	nbytes_to_read = header->numberOfSamples * (file->bitsPerSample / 8) * header->numberOfChannels;
	blocks_read=0;
	
	update_on_block = (int)(1024*1024*10/file->blockLength);
	if(update_on_block < 1)update_on_block=1;
	// we now try and read the entire file, 'block' at a time into the float arrays.
	while (nbytes_to_read > 0) {
		read = readPsrXmlNextDataBlockIntoExistingArray(file, byte_array);
		if (read <= 0) {
			// we have reached the end of the file...
			fprintf(stdout, "\n\nReached the end of file...\n");

			nbytes_to_read = 0;
			break;
		}
		nbytes_to_read -= read;

		nsamps = (read * 8/file->bitsPerSample) / header->numberOfChannels;
		end = start+nsamps;
		unpackDataChunk(byte_array, float_array,header,0,nsamps,0,nsamps,swap_chans);
		unpackToChannels(float_array,block_chan_array,header->numberOfChannels,nsamps);	

		// Now I do a memcopy. this is probably not the most efficient way to do it but I don't want to have to transform the entire file.
		for (chan=0; chan < header->numberOfChannels; chan++){
			memcpy(output[chan]+start,block_chan_array[chan], nsamps*sizeof(float));
		}
		start+=nsamps;
		blocks_read++;
		if(!(blocks_read % update_on_block)){	
			fprintf(stdout,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			fprintf(stdout,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			fprintf(stdout,"Read %04d data blocks, %09d samples",blocks_read, start);
			fflush(stdout);
		}
	}
	fprintf(stdout,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	fprintf(stdout,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

	fprintf(stdout,"Read %04d data blocks, %09d samples\n",blocks_read, start);


	// these frees should match the ones above except for the 'output' array.
	free(byte_array);
	free(block_chan_array);
	free(float_array);



	return output;

}

void pch_seek_init_operations(pch_seek_operations_t* operations){
	operations->dump_tim=0;
	operations->fft_input=0;
	operations->form_amplitudes=0;
	operations->twiddle_amplitudes=0;
	operations->dump_amplitudes=0;
	operations->dump_phases=0;
	operations->dump_normalised=0;
	operations->dump_harmfolds=0;
	operations->phase_fit=0;
	operations->fscrunch=0;
	operations->normalise_median=0;
	operations->normalise_powerlaw=0;
	operations->normalise_agl=0;

	operations->harmfold_simple=0;
	operations->search_amplitudes=0;
	operations->write_prd=0;
	operations->recon_add=0;

	operations->ndm=0;
	operations->dmtrials=NULL;
	operations->harmfolds=NULL;
	operations->nharms=0;
	strcpy(operations->prdfile,"out.prd");
	operations->amp_thresh=5;

}
