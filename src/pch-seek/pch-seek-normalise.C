#include "pch-seek.h"

#include "TKfit.h"
#include "toolkit.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

#define median_NMAX 8192



/**
 * Try and normalise the data.
 *
 * This aproach fits a power law to the data and then divides by this
 *
 * Then we do the AGL mean/rms fitting with a large block size, to remove any residual effects.
 */
void pch_seek_normalise_powerlaw(float* amplitudes, int ndat){
	double* log_x;
	double* log_y;
	double p[2];
	int m,x;
	float mult,power;
	float mean,rms;

	log_x = (double*)malloc(sizeof(double)*(int)(ndat));
	log_y = (double*)malloc(sizeof(double)*(int)(ndat));



	x=0;
	for (int i = 1; i < ndat; i+=10){
		if(amplitudes[i] < 0.1) continue;
		log_x[x] = log((double)i);
		log_y[x] = log((double)amplitudes[i]);
//		printf("%lf %lf\n",log_x[x],log_y[x]);
		x++;
	}


	m=2;
	TKleastSquares_svd_noErr(log_x,log_y,(int)(ndat/1024), p, m, TKfitPoly);

	mult=exp(p[0]);
	power=p[1];

//	printf("--\n%f %f\n",mult,power);

	free(log_x);
	free(log_y);

	for (int i = 0; i < ndat; i++){
		if(amplitudes[i] > 0.1){
			amplitudes[i] /= mult*pow((float)i,power);
		}
	}

/*	meanrms(amplitudes, ndat,&mean,&rms);

	for (int i = 0; i < ndat; i++){
		amplitudes[i] = (amplitudes[i]-mean)/rms;
	}
*/

	pch_seek_normalise_agl_mean(amplitudes,ndat,16384);


}

/**
 *
 * Use the sigproc 'AGL' whitening method.
 *
 * This does mean and rms of data blocks, zapping things that are more than 3 sigma away.
 *
 *
 */
void pch_seek_normalise_agl_mean(float* amplitudes, int ndat, int nrun){
	int j,l,h,k;
	float mean,rms;
	float* small_block_array;
	float* rmean,*rrms;
	int zapped;

	rmean = (float*) malloc(sizeof(float)*(int)(ndat/nrun+1.5));
	rrms = (float*) malloc(sizeof(float)*(int)(ndat/nrun+1.5));

	small_block_array = (float*) malloc(sizeof(float)*nrun);

	l=0;
	j=0;
	k=0;
	// This code isn't very optimised as I copied it directly from ralph/andrew's code.
	for ( int i = 1; i < ndat; i++){
		if(fabs(amplitudes[i]) > 0.000001){
			l++;
			small_block_array[l]=amplitudes[i];
		}
		j++;
		if(j==nrun || i>=(ndat-1)){
			mean=0;
			rms=0;

			int nrej=0;
			int nrejold=-1;
			float sum=0;
			float sumsq=0;
			for(int x=0; x < l; x++){
				sum+=small_block_array[x];
				sumsq+=small_block_array[x]*small_block_array[x];
			}

			int counter=0;
			float s=sum;
			float ss=sumsq;
			while(counter < 6 && nrej != nrejold){
				nrejold=nrej;
				counter++;
				float vv=(sumsq-sum*sum/(float)(l-nrej))/(float)(l-nrej-1.0);
				mean=sum/(float)(l-nrej);
				if(vv < 0)vv=0;
				rms = sqrt(vv);
				sum=s;
				sumsq=ss;
				nrej=0;
				for(int x=0; x < l; x++){
					if(fabs(small_block_array[x]-mean) > 3*rms){
						sum-=small_block_array[x];
						sumsq-=small_block_array[x]*small_block_array[x];
						nrej++;
					}
				}
			}

			/*
			zapped=0;
			for(int x = 0; x < 3; x++){
				if(rms> 0){
					// zap things more than 3 sigma away.
					for(int y=1; y < l; y++){
						if((small_block_array[y]-mean) > 3*rms){
							zapped++;
							small_block_array[y]=0;
						}
					}
				}
				//				if(zapped > 0.1*nrun) break;
				mean=0;
				rms=0;
				if(l > 128){
					meanrms(small_block_array,l,&mean,&rms);
				} else {
					if(k > 1){
						mean=rmean[k-1];
						rms=rrms[k-1];
					}
				}
			}
			*/
			rmean[k]=mean;
			rrms[k]=rms;
			k++;
			l=0;
			j=0;
		}
	}


	// Now do the 'agl' whitening
	h=-1;
	for(int i = 0; i < ndat; i++){
		if (!(i%nrun))h++; // every nrun, increment h.
		if(h>=k) h=k-1;;
		if (rrms[h] < 0.00001 || rmean[h] == 0)  amplitudes[i] = 0.0;
		else amplitudes[i] = (amplitudes[i] - rmean[h]) / rrms[h];
	}

	free(small_block_array);
	free(rmean);
	free(rrms);
}


/**
 *
 * Normalise the data by computing the median every nrun bins
 *
 */
void pch_seek_normalise_median(float* amplitudes, int ndat, int nrun){
	int k,h;
	float* small_block_array;
	float* rmed;

	rmed = (float*) malloc(sizeof(float)*(int)(ndat/nrun+0.5));
	small_block_array = (float*) malloc(sizeof(float)*nrun);

	k=0;

	for ( int i = 0; i < ndat; i+=nrun){
		memcpy(small_block_array,amplitudes+i,nrun*sizeof(float));
		quicksort_inplace(small_block_array,nrun);
		rmed[k]=small_block_array[(int)(nrun/2)];
		k++;
	}


	// Now do the 'eatough-lyne' median whitening
	h=-1;
	for(int i = 0; i < ndat; i++){
		if (!(i%nrun))h++; // every nrun, increment h.
		if(h>=k) h=k-1;;
		if (fabs(rmed[h]) < 0.001)  amplitudes[i] = 0.0;
		else amplitudes[i] = (amplitudes[i]/rmed[h])-1.0;
	}

	free(small_block_array);
	free(rmed);
}

