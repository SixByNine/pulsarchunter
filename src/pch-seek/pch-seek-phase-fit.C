#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pch-seek.h"

#define KDM 4.148808e-3
#define TWOPI 6.28318531
#define PI 3.14159265


void pch_seek_phase_fit_simple(float** phases, int nsamp, int nchan, float* dmtrials, float ndm, float fch1, float foff, float tobs){

	int samp,chan,idm,i;
	float dm;
	float* sample_phases, *best_sample_phases;
	float** delays;
	float phase_off;
	float bin_freq,freq_mult;
	float fch,fhi;
	float sum,sumsq;
	float rms,mean,cmean;
	float st_sum,st_sumsq;
	float best_sum,best_sumsq,best_rms;
	int ifold,nfold,fold[8],meanrms[8];
	int nsubfold,foldbin;
	int navg;

	FILE* file;

	nfold=8;
	fold[0]=1;
	fold[1]=2;
	fold[2]=3;
	fold[3]=4;
	fold[4]=6;
	fold[5]=8;
	fold[6]=16;
	fold[7]=32;
	for (i=0; i<nfold; i++)meanrms[i]=0;

	file = fopen("phase_test","w");

	freq_mult=1.0/tobs;
	sample_phases=(float*)malloc(sizeof(float)*nchan);
	best_sample_phases=(float*)malloc(sizeof(float)*nchan);
	delays = (float**)malloc(sizeof(float*)*nchan);
	fhi = fch1/1000.0;
	for (chan = 0; chan < nchan; chan++){
		delays[chan] = (float*)malloc(sizeof(float)*ndm);
		fch=(fch1 + chan*foff)/1000.0;
		
		for (idm = 0 ; idm < ndm ; idm++){
			delays[chan][idm] = dmtrials[idm] * (KDM * (1.0 / (fhi * fhi) - 1.0
						/ (fch * fch)));
		}
	}

	/**
	 *
	 * Here we begin the work.
	 *
	 * There are a lot of nested loops here:
	 *  over 'samples' (i.e. phases from the FFT)
	 *  over trial DM (should be relplaced with a DM fitting algorithm!)
	 *  over 'harmonic folds'
	 *  over 'sub-folds' because each bin maps to more than one harmonicaly related bin
	 *  over spectral channels
	 */
	for (samp = 1 ; samp < nsamp ; samp++){
		bin_freq = samp*freq_mult;
		// loop over each dm trial
		for (idm = 0 ; idm < ndm ; idm++){
			dm = dmtrials[idm];
			sumsq=0;
			sum=0;
			st_sum=sum;
			st_sumsq=sumsq;
			cmean=0;
			for (ifold=0; ifold<nfold; ifold++){
				nsubfold=fold[ifold]/2;
				best_rms=100000;
				if(ifold==0)nsubfold=0;
				for(i= -nsubfold; i <= nsubfold ; i++){
					foldbin = fold[ifold]*samp + i;
					if (foldbin >= nsamp)continue;

					//if(samp==10801 && fold[ifold]==6)foldbin=fold[ifold]*samp-2;
					sum=st_sum;
					sumsq=st_sumsq;
					for (chan = 0; chan < nchan; chan++){
						phase_off = delays[chan][idm]*(bin_freq*fold[ifold]) * TWOPI;
						phase_off -= (fold[ifold]-1)*cmean; // rotate the harmonics if cmean is set
//						if (foldbin==64804)printf("phaseoff %d %d %f %f %f\n",samp,chan,phase_off,cmean,phases[chan][foldbin]);
						sample_phases[chan] = phases[chan][foldbin] - phase_off;
						//@TODO: investigage the problems where points fall close to +/- PI.
						while(sample_phases[chan] > PI) sample_phases[chan] -= TWOPI;
						while(sample_phases[chan] < -PI) sample_phases[chan] += TWOPI;
						sum  += sample_phases[chan];
						sumsq+= sample_phases[chan]*sample_phases[chan];
					}
					navg = (ifold+1)*nchan;
					mean = sum/navg;
					if(ifold==0)cmean=mean;
					rms  = sqrt((sumsq/navg - mean*mean) / (ifold+1));
					if(rms < best_rms){
						best_rms =rms;
						best_sum=sum;
						best_sumsq=sumsq;
						memcpy(best_sample_phases,sample_phases,sizeof(float)*nchan);
					}
				}
				st_sum=best_sum;
				st_sumsq=best_sumsq;
				meanrms[ifold] += best_rms;
				if(best_rms < 1.0){
					fprintf(file,"%d %f %f %f %d %f *",fold[ifold],(float)foldbin/tobs/fold[ifold],dm,1.0/best_rms,foldbin,mean);
					for (chan = 0; chan < nchan; chan++){
						fprintf(file," %f",best_sample_phases[chan]);
					}
					fprintf(file,"\n");
				}

			}
		}
	}

	for (i=0; i<nfold; i++)printf("%d %f\n",i,meanrms[i]/nsamp);
	

	free(sample_phases);
	free(delays);

	fclose(file);

}
