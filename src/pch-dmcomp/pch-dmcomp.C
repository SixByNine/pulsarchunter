#include "pch-dmcomp.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <psrxml.h>
#include <getopt.h>
#include <math.h>

#define KDM 8.287616

int set_options(struct option* long_opt, int* opt_flag);
void pch_dmcomp_print_usage();

/**
 * PulsarCHunter SEEK
 *
 * This is a program for finding pulsars.
 *
 * This 'main' method doesn't really do anything,
 * just read the file and calls the 'do_search' function 
 * which calls the various search commands.
 *
 * 2008-12-22	M.Keith		Initial version
 *
 */
int main(int argc, char** argv){
	psrxml* header;
	const char* args = "2d:FhSN:l";
	char c;
	int long_opt_idx = 0;
	int opt_flag = 0;
	float dm =0;
	struct option long_opt[256];
	bool pow2only=false;
	bool compute_fourier_size=false, compute_scrunch_factor=false;
	bool compute_dm_band_dly=false;
	bool compute_ndm_steps=false;
	int ndm=1;

	header = (psrxml*) malloc(sizeof (psrxml));

	// arrange for the options struct to be filled out
	// this uses the GNU standard 'getopt.h' and so might not work well
	// on a Solaris/UNIX machine.
	long_opt_idx = set_options(long_opt,&opt_flag);

	// now we read the options
	// This is a rather large switch block... with another switch inside!
	// try and keep up!
	while ((c = getopt_long(argc, argv, args, long_opt, &long_opt_idx)) != -1) {
		// handle args
		switch (opt_flag) {

			default:
				opt_flag = 0;
				switch (c) {
					case 'h':
						pch_dmcomp_print_usage();
						exit(0);
						break;
					case 'd':
						dm = atof(optarg);
						break;
					case 'F':
						compute_fourier_size=true;
						break;
					case '2':
						pow2only=true;
						break;
					case 'S':
						compute_scrunch_factor=true;
						break;
					case 'l':
						compute_dm_band_dly=true;
						break;

					case 'N':
						compute_ndm_steps=true;
						ndm=atoi(optarg);
						break;
				}
		}
		opt_flag=0;
	}

	// Options are now read!
	// Now we read the psrxml file

	if (optind >= argc) {
		printf("Error: Need a psrxml file to work on... (use --help for options)\n");
		exit(1);
	}

	readPsrXml(header, argv[optind]);
	float cfreq = header->centreFreqCh1 + header->freqOffset*header->numberOfChannels/2;
	cfreq /= 1000.0; // GHz
	float dm_t = KDM *1e-6 * fabs(header->freqOffset) * pow(cfreq,-3) * dm; // seconds
	int nsamp=header->numberOfSamples;
	

	if(compute_dm_band_dly){
		float f1 = header->centreFreqCh1;
		float f2 = header->centreFreqCh1 + header->freqOffset*header->numberOfChannels;
		f1 /= 1000.0; // GHz
		f2 /= 1000.0; // GHz
		float dm_band_dly = KDM/2.0 * 1e-3 * fabs(pow(f1,-2)-pow(f2,-2)) * dm; // seconds
		printf("%f %d\n",dm_band_dly,(int)(dm_band_dly/header->currentSampleInterval));
	}
	if (compute_ndm_steps){
		float maxdm = pch_getDMtable_ndm(dm, ndm, header->currentSampleInterval*1e6, 40.0, header->freqOffset, cfreq, header->numberOfChannels, 1.25);
		printf("%f\n",maxdm);
	}

	if (compute_scrunch_factor){
		float tsamp=header->currentSampleInterval;
		int scrfactor=1;
		// check for dm smeering
		if (dm > 0){
			while ( dm_t > 2*tsamp ){
				tsamp*=2;
				scrfactor*=2;
			}
			nsamp /= scrfactor;
		}

		printf("%d\n",scrfactor);


	}


	if (compute_fourier_size){

		int size=pch_seek_fourier_size(nsamp,pow2only);
		printf("%d\n",size);
	}

	return 0;
}


/*
 *
 * Configures the options array for getopt
 * 
 *
 * Flag values are:
 * 00000001	dump timeseries
 * 00000002	???
 *
 */
int set_options(struct option* long_opt, int* opt_flag){
	int long_opt_idx = 0;
	long_opt[long_opt_idx].name = "help";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'h';

	long_opt[long_opt_idx].name = "dm";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'd';

	long_opt[long_opt_idx].name = "fourier-size";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'F';

	long_opt[long_opt_idx].name = "power-of-two";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = '2';

	long_opt[long_opt_idx].name = "scrunch-factor";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'S';

	long_opt[long_opt_idx].name = "ndm-steps";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'N';




	// The TERMINATOR
	long_opt[long_opt_idx].name = 0;
	long_opt[long_opt_idx].has_arg = 0;
	long_opt[long_opt_idx].flag = 0;
	long_opt[long_opt_idx++].val = 0;

	return long_opt_idx;
}

void pch_dmcomp_print_usage(){

	printf("-F\tSensible Fourier size\n");
	printf("-2\tOnly power of 2\n");
	printf("-N\tWhat DM after N DM Trials\n");
	printf("-D\tSet DM used\n");
}

