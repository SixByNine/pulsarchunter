#include "pch-seek.h"
#include <config.h>

#include "toolkit.h"

#include <psrxml.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/**
 * A 'dummy' operation that dumps spectra etc.
 */
void pch_seek_dump(float* data, int len, float xscale, char* filename){

	FILE* file;
	int i;

	file = fopen(filename,"w");

	for(i=0; i < len; i++){
		fprintf(file,"%f %f\n",xscale*i,data[i]);
	}

	fclose(file);

	return;
}

void pch_seek_write_prd(char* filename, float** freq, float** spec, float** recon, int* ncand, int* harms, int nharm, psrxml* header){

	int maxncand=0;
	int ifold;
	int icand;
	int** indexes;
	float** sorter;

	FILE *file;

	file = fopen(filename, "w");

	indexes = (int**)malloc(sizeof(int*)*nharm);

	if(recon!=NULL)sorter=recon;
	else sorter=spec;

	for(ifold = 0; ifold < nharm; ifold++){
		// sort the candidates
		indexes[ifold] = (int*) malloc(sizeof(int)*ncand[ifold]);

		for (icand=0; icand < ncand[ifold]; icand++){
			indexes[ifold][icand]=icand;
		}
		quicksort_index(sorter[ifold],indexes[ifold],ncand[ifold]);
		// reverse the order of the array
		{
			int swap;
			int i,j;
			j=ncand[ifold];
			for(i=0;i<--j;i++){
				swap=indexes[ifold][i];
				indexes[ifold][i]=indexes[ifold][j];
				indexes[ifold][j]=swap;
			}
		}

		// remove freqs that are very close
		float* existing_freqs=(float*) malloc(sizeof(float)*ncand[ifold]);
		int nex=0;
		char found=0;
		int remd=0;
		for (icand=0; icand < ncand[ifold] ; icand++){
			found = 0;

			for(int x=0; x < nex; x++){
				if(fabs(existing_freqs[x]/freq[ifold][indexes[ifold][icand]]-1) < 0.001){
					// swap index to the end
					for(int y=icand+1; y < ncand[ifold]; y++){
						indexes[ifold][y-1]=indexes[ifold][y];
					}
					// we have one less cand
					ncand[ifold]--;
					// because we have changed which cand is at 'icand' we want to do it again
					icand--;
					// end the loop
					found=1;
					break;
				}
			}
			if(!found){
				existing_freqs[nex] = freq[ifold][indexes[ifold][icand]];
				nex++;
			}
		}
		if(maxncand < ncand[ifold])maxncand=ncand[ifold];

	}

	/*
	 * Add a header
	 *##BEGIN HEADER##
	 SOURCEID = GEGJ1027K
	 FREF =   3796.5 MHz
	 TSTART =   54562.2848
	 TELESCOPE = Parkes
	 RAJ = 10:28:20.2
	 DECJ = -58:20:11.1
	 TSAMP =    250.0 us
	 PROGRAM = SEEK
	 VERSION = 4.3                             
	 HARM_FOLDS =    1   2   4   8  16
	 COLS = SNR_SPEC SNR_RECON PERIOD
##END HEADER##

*/

	fprintf(file,"##BEGIN HEADER##\n");
	fprintf(file,"SOURCEID = %s\n",header->sourceName);
	fprintf(file,"FREF = %f MHz\n",header->centreFreqCh1);
	fprintf(file,"TSTART = %f\n",header->mjdObs);
	fprintf(file,"TELESCOPE = %s\n",header->telescope.name);

	fprintf(file,"RAJ = %f\n",header->startCoordinate.ra);
	fprintf(file,"DECJ = %f\n",header->startCoordinate.dec);

	fprintf(file,"TSAMP = %f us\n",header->currentSampleInterval);
	fprintf(file,"PROGRAM = %s\n",PACKAGE_NAME);
	fprintf(file,"VERSION = %f\n",PACKAGE_VERSION);

	fprintf(file,"HARM_FOLDS =");
	for(ifold = 0; ifold < nharm; ifold++){
		fprintf(file," %d",harms[ifold]);
	}
	fprintf(file,"\n");

	if(recon==NULL){
		fprintf(file,"COLS = SNR_SPEC PERIOD\n");

	} else {
		fprintf(file,"COLS = SNR_SPEC SNR_RECON PERIOD\n");

	}



	fprintf(file,"##END HEADER##\n");


	fprintf(file," DM:% 11.6f      AC:% 11.6f      AD:% 11.6f\n",header->referenceDm,0,0);

	for (icand=0; icand < maxncand; icand++){
		for(ifold = 0; ifold < nharm; ifold++){
			if(icand >= ncand[ifold]){
				// no more cands in this fold
				fprintf(file,"     0.0");
				if(recon!=NULL)fprintf(file,"     0.0");
				fprintf(file,"    1.00000000");
			} else {
				fprintf(file,"% 8.1f",spec[ifold][indexes[ifold][icand]]);
				if(recon!=NULL)fprintf(file,"% 8.1f",recon[ifold][indexes[ifold][icand]]);
				fprintf(file,"% 13.8f",1000.0/freq[ifold][indexes[ifold][icand]]);
			}
		}
		fprintf(file,"\n");
	}

	fclose(file);

	for(ifold = 0; ifold < nharm; ifold++){
		free(indexes[ifold]);
	}
	free(indexes);

}
