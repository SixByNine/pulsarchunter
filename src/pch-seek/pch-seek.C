#include "pch-seek.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <psrxml.h>
#include <fftw3.h>
#include <getopt.h>
#include <math.h>

int set_options(struct option* long_opt, int* opt_flag);
void pch_seek_print_usage();

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

	float** data_array;
	psrxml* header;
	const char* args = "Adhprt:T:G:H:";
	char c;
	FILE *dmfile;
	int dmtrials_size;
	int threads=0;
	int long_opt_idx = 0;
	int opt_flag = 0;
	int scrunch_factor=1;
	int fft_size = -1;
	struct option long_opt[256];
	pch_seek_operations_t operations;

	header = (psrxml*) malloc(sizeof (psrxml));

	dmfile=NULL;

	pch_seek_init_operations(&operations);

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
			case 0x00000001: // dump-tim
				operations.dump_tim=1;
				break;
			case 0x00000002: // dump-amplitudes
				operations.dump_amplitudes=1;
				break;
			case 0x00000003: // dump-phases
				operations.dump_phases=1;
				break;
			case 0x00000004: // dump-normalised
				operations.dump_normalised=1;
				break;
			case 0x00000005: // dump-harmfolds
				operations.dump_harmfolds=1;
				break;

			case 0x00000008: // hist-tim
				operations.hist_tim=1;
				break;
			case 0x00000009: // hist-amplitudes
				operations.hist_amplitudes=1;
				break;
			case 0x0000000A: // hist-tim
                                operations.hist_normalised=1;
                                break;
                        case 0x0000000B: // hist-amplitudes
                                operations.hist_harmfolds=1;
                                break;





			case 0x00000011: // normaise-median
				operations.normalise_median=1;
				break;
			case 0x00000012: // normaise-agl
				operations.normalise_agl=1;
				break;
			case 0x00000013: // normalise-powerlaw
				operations.normalise_powerlaw=1;
				break;


			case 0x00000021: //harmfold-simple
				operations.harmfold_simple=1;
				break;

			case 0x00000031:
				operations.recon_ralph=1;
				break;
			case 0x00000101: //twiddle-amplitudes
				operations.twiddle_amplitudes=1;
				break;
			case 0x00000111: // write-prd
				operations.write_prd=1;
				strcpy(operations.prdfile,optarg);
				break;
			case 0x00000112: // write-presto-fft
				operations.write_presto_fft=1;
				strcpy(operations.presto_fft_file,optarg);
				break;

			case 0x00000121: // amp-thresh
				operations.amp_thresh=atof(optarg);
				break;
			case 0x00000122: // hfold-bonus-factor
				operations.hfold_bonus_factor=atof(optarg);
				break;



			case 0x00000201: //seach-chans
				operations.search_chans=1;
				break;
			case 0x00000202: //tscrunch
				scrunch_factor = atol(optarg);
				break;



			default:
				opt_flag = 0;
				switch (c) {
					case 'h':
						pch_seek_print_usage();
						exit(0);
						break;
					case 'p':
						operations.phase_fit=1;
						break;
					case 'T':
						threads=atoi(optarg);
#ifndef FFTW_THREADS
						fprintf(stderr,"Warning: FFTW was not compiled in multi-threaded mode.\n");
#endif
						break;
					case 't': //fft-size
						fft_size = atol(optarg);
						break;


					case 'd': // dm file
						dmfile = fopen(optarg,"r");
						break;

					case 'r': //recon-add
						operations.recon_add=1;
						break;
					case 'G':
						operations.giant_search=1;
						strcpy(operations.giantfile,optarg);
						break;

					case 'H': //select harmonic folds
						{
							// this is a little difficult to work out!
							// (because C is for monkeys)
							char* end;
							char* str;
							int val;
							int i,arrlen;
							arrlen=8;
							i=0;
							operations.harmfolds=(int*)malloc(sizeof(int)*arrlen);
							str=optarg;
							end=str+strlen(str);
							while (str<end){
								if (i > arrlen){
									arrlen*=2;
									operations.harmfolds=(int*)realloc(operations.harmfolds,(sizeof(int)*arrlen));
								}
								sscanf(str,"%d",&val);
								operations.harmfolds[i]=val;
								i++;
								while(str[0] != ' ' && str[0] !='\0')str++;
								str++;
							}
							operations.nharms=i;
						}

						break;
					case 'A':
						operations.append_output=1;
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

	if(fft_size > 1){
		// we pretend the file is shorter than it is to force the size to be sensible
		int read_length = (int)pow(2,fft_size)*scrunch_factor;
		if (read_length < header->numberOfSamples){
			header->numberOfSamples = read_length;
			header->actualObsTime=read_length*header->currentSampleInterval;
		}

	}

	/*
	 * If we are using a 'dm file' to specify some dispersion measures to use
	 * read them in now.
	 *
	 * This is used by the phase-search algorithm only at the moment.
	 *
	 */
	if(dmfile!=NULL){
		dmtrials_size = 10;
		operations.dmtrials = (float*)malloc(sizeof(float)*dmtrials_size);
		operations.ndm=0;
		while(!feof(dmfile)){
			if (operations.ndm > dmtrials_size){	
				dmtrials_size*=2;
				operations.dmtrials=(float*)realloc(operations.dmtrials,sizeof(float)*dmtrials_size);
			}
			fscanf(dmfile,"%f\n",operations.dmtrials+operations.ndm);
			operations.ndm++;
		}
		printf("Using %d DM trials\n",operations.ndm);
	}

	// Now we read the data file into memory, if it fits!
	data_array = pch_seek_read_file(header,scrunch_factor);
	if(data_array == NULL){
		// we cannot contiune... data couldn't be read or memory couldn't be allocated.
		fprintf(stderr,"Error, could not read data\n");
		return 1;
	}

#ifdef FFTW_THREADS
	// If we are using FFTW in multi-threaded mode, we initialise it here.
	if (threads){
		printf("Using %d threads\n",threads);
		fftwf_init_threads();
		fftwf_plan_with_nthreads(threads);
	}
#endif

	// call this function to do the work. 
	// The arrays passed are destroyed and free'd in the operation
	// of the search, so we don't need to free the data array.
	pch_seek_do_search(&operations,header,data_array);

#ifdef FFTW_THREADS
	if(threads){
		fftwf_cleanup_threads();
	}
#endif

	free(operations.dmtrials);
	if(operations.harmfolds!=NULL)free(operations.harmfolds);

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

	long_opt[long_opt_idx].name = "dump-tim";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000001;

	long_opt[long_opt_idx].name = "dump-amplitudes";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000002;

	long_opt[long_opt_idx].name = "dump-phases";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000003;

	long_opt[long_opt_idx].name = "dump-normalised";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000004;

	long_opt[long_opt_idx].name = "dump-harmfolds";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000005;

	long_opt[long_opt_idx].name = "hist-tim";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000008;

	long_opt[long_opt_idx].name = "hist-amplitudes";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000009;

	long_opt[long_opt_idx].name = "hist-normalised";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x0000000A;

	long_opt[long_opt_idx].name = "hist-harmfolds";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x0000000B;






	long_opt[long_opt_idx].name = "phase-fit";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'p';

	long_opt[long_opt_idx].name = "threads";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'T';

	long_opt[long_opt_idx].name = "dm-trials";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'd';

	long_opt[long_opt_idx].name = "normalise-median";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000011;

	long_opt[long_opt_idx].name = "normalise-agl";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000012;

	long_opt[long_opt_idx].name = "normalise-powerlaw";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000013;



	long_opt[long_opt_idx].name = "harmfold-simple";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000021;

	long_opt[long_opt_idx].name = "harm-folds";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'H';

	long_opt[long_opt_idx].name = "recon-add";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'r';
	long_opt[long_opt_idx].name = "recon-ralph";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000031;



	long_opt[long_opt_idx].name = "twiddle-amps";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000101;

	long_opt[long_opt_idx].name = "enable-interbin";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000101;



	long_opt[long_opt_idx].name = "write-prd";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000111;

	long_opt[long_opt_idx].name = "write-presto-fft";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000112;
	long_opt[long_opt_idx].name = "append-output";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'A';


	long_opt[long_opt_idx].name = "amp-thresh";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000121;

	long_opt[long_opt_idx].name = "hfold-bonus-factor";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000122;



	long_opt[long_opt_idx].name = "giant-search";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'G';


	long_opt[long_opt_idx].name = "search-chans";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000201;


	long_opt[long_opt_idx].name = "tscrunch";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = opt_flag;
	long_opt[long_opt_idx++].val = 0x00000202;


	long_opt[long_opt_idx].name = "fft-size";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 't';




	// The TERMINATOR
	long_opt[long_opt_idx].name = 0;
	long_opt[long_opt_idx].has_arg = 0;
	long_opt[long_opt_idx].flag = 0;
	long_opt[long_opt_idx++].val = 0;

	return long_opt_idx;
}

void pch_seek_print_usage(){

}

