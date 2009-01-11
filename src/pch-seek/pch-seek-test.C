#include "pch-seek.h"
#include "toolkit.h"

#include <stdlib.h>
#include <stdio.h>
#include <psrxml.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#define PI 3.14159265


int test_single_chan();

/**
 *
 * Run some test cases
 *
 *
 */
int main(int argc, char** argv){
	return test_single_chan();
}

int test_single_chan(){

	float* data;
	float* data_run;
	float** chans;
	int npts;
	float tsamp;
	long seed;
	float ffac;
	float sinfac;
	float cperiod;
	float mean,rms;
	int hfolds[4];
	psrxml* header;
	pch_seek_operations_t ops;

	hfolds[0]=2;
	hfolds[1]=4;
	hfolds[2]=8;
	hfolds[3]=16;


	seed = 1024;
	tsamp = 250/1e6; //250 us
	cperiod = 0.1; //s

	npts = (int) (pow(2,23));
	data = (float*) fftwf_malloc(sizeof(float)*npts);
	ffac = 10/(sqrt(npts));
	sinfac = tsamp*2*PI/cperiod;

	header = (psrxml*) malloc(sizeof(psrxml));

	header->numberOfSamples=npts;
	header->currentSampleInterval=tsamp;
	header->actualObsTime=tsamp*(float)npts;
	header->numberOfChannels = 1;


	pch_seek_init_operations(&ops);

	ops.write_prd=1;
	ops.normalise_agl=1;
	ops.harmfolds=hfolds;
	ops.nharms=4;
	ops.twiddle_amplitudes=1;
	ops.dump_amplitudes=1;
	ops.dump_normalised=1;
		




	// initialise the array with gausian distributed data
	for(int i = 0; i < npts; i++){
		data[i] = TKgaussDev(&seed);
	}
	meanrms(data,npts,&mean,&rms);
	printf("mean: %f rms: %f\n",mean,rms);

	pch_seek_fourier_r2c(data,npts);

	        data_run=(float*) fftwf_malloc(sizeof(float)*npts);

	for(int i = 0; i < npts; i++){
		data_run[i] = TKgaussDev(&seed);
	}

	pch_seek_fourier_r2c(data_run,npts);

	for(int i = 0; i < npts; i++){
//		data[i] = data[i] + data_run[i];
	}


	for(int i = 0; i < npts/2-1; i++){
		data[2*i] = sqrt(data[2*i]*data[2*i] + data[2*i+1]*data[2*i+1]);
	}
	meanrms(data,npts/2,&mean,&rms);
	printf("mean: %f rms: %f\n",mean,rms);



	return 0;



	// initialise the array with gausian distributed data
	for(int i = 0; i < npts; i++){
		data[i] = TKgaussDev(&seed);
	}



	strcpy(ops.prdfile,"test_1ch_noise.prd");
	data_run=(float*) fftwf_malloc(sizeof(float)*npts);
	memcpy(data_run,data,sizeof(float)*npts);
	chans = (float**) malloc(sizeof(float*));
	chans[0]=data_run;
	pch_seek_do_search(&ops,header,chans);
	

	return 0;



	// now add a sine wave
	for(int i = 0; i < npts; i++){
		data[i] += ffac*sin((float)i*sinfac);
	}

	strcpy(ops.prdfile,"test_1ch_1h.prd");
	data_run=(float*) fftwf_malloc(sizeof(float)*npts);
	memcpy(data_run,data,sizeof(float)*npts);
	chans = (float**) malloc(sizeof(float*));
	chans[0]=data_run;
	pch_seek_do_search(&ops,header,chans);


        // now add a 2nd harmonic sine wave
        for(int i = 0; i < npts; i++){
                data[i] += (ffac/sqrt(1.0))*sin((float)i*sinfac*2);
        }

        strcpy(ops.prdfile,"test_1ch_2h.prd");
        data_run=(float*) fftwf_malloc(sizeof(float)*npts);
        memcpy(data_run,data,sizeof(float)*npts);
        chans = (float**) malloc(sizeof(float*));
        chans[0]=data_run;
        pch_seek_do_search(&ops,header,chans);


        // now add a 3rd and 4th harmonic sine wave
        for(int i = 0; i < npts; i++){
                data[i] += (ffac/sqrt(1.0))*sin((float)i*sinfac*3);
                data[i] += (ffac/sqrt(1.0))*sin((float)i*sinfac*4);
        }

        strcpy(ops.prdfile,"test_1ch_4h.prd");
        data_run=(float*) fftwf_malloc(sizeof(float)*npts);
        memcpy(data_run,data,sizeof(float)*npts);
        chans = (float**) malloc(sizeof(float*));
        chans[0]=data_run;
        pch_seek_do_search(&ops,header,chans);


	return 0;
}
