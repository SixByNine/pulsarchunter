#include <config.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "TKfit.h"
#include "toolkit.h"

#define median_NMAX 8192


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
	for ( int i = 0; i < ndat; i++){
		if(fabs(amplitudes[i]) > 0.000001){
			small_block_array[l]=amplitudes[i];
			l++;
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



/*
 * This is the M.Keith attempt at spectral whitening
 *
 * It subtracts an interpolated median of increasing block sizes (as in Presto),
 * Then divides by the interpolated RMS of the resultant -ve half of the data.
 */
void pch_seek_normalise_median_smoothed(float* amplitudes, int ndat, int nrun){

	int binnum,npts;
	int nrunmax,nrunstart,nruno,nrunoo;
	float *buf_old,*buf_new,*realbuffer, *buf_old_old;
	float median,medlow,new_median,new_medlow,slope_median,slope_medlow,new_rms,slope_rms;
	float ssq,rms;
	float u_quartile,l_quartile;


	nrunstart=nrun;
	nrunmax=1024;
	realbuffer=(float*)malloc(sizeof(float)*nrunmax);

	buf_old=amplitudes;
	binnum=0;
	memcpy(realbuffer,buf_old,nrun*sizeof(float));
	quicksort_inplace(realbuffer,nrun);
	npts=nrun;
	while(realbuffer[npts-1] < 1e-12){
		npts--;
	}

	
	median=realbuffer[(int)(npts/2)];
//	medlow=median - realbuffer[(int)(5*npts/6)]; // need some conversion factor?
	u_quartile=realbuffer[(int)(npts/4)];
	l_quartile=realbuffer[(int)(3*npts/4)];

	nruno=nrun;
	nrun*=2;
	binnum+=nruno;
	buf_new=NULL;
	while (binnum+nrun < ndat){
		buf_new = buf_old+nruno;

		memcpy(realbuffer,buf_new,nrun*sizeof(float));
		quicksort_inplace(realbuffer,nrun);
		npts=nrun;
		while(realbuffer[npts-1] < 1e-12){
			npts--;
		}

		new_median=realbuffer[(int)(npts/2)];
//		new_medlow=(realbuffer[(int)(npts/6)] - realbuffer[(int)(5*npts/6)])/2.0; // need some conversion factor?

		slope_median=(new_median-median)/nrun;
//		slope_medlow=(new_medlow-medlow)/nrun;
		slope_rms=0;//(new_rms-rms)/nrun;

		for (int i = 0; i < nruno; i++){
			if(buf_old[i] < 1e-12){
				buf_old[i]=0;
			} else {
				buf_old[i] -= median + slope_median*i;
			}
		}
		ssq=0;
		int count=0;
		for (int i = 0; i < nruno; i++){
//			if(buf_old[i] > l_quartile && buf_old[i] < u_quartile && fabs(buf_old[i]) > 1e-12){
			if(buf_old[i] < 0){
				ssq+=buf_old[i]*buf_old[i];
				count++;
			}
		}

		new_rms=sqrt(ssq/count);

		slope_rms=(new_rms-rms)/nruno;

		if(buf_old_old!=NULL){
			for (int i = 0; i < nrunoo; i++){
				buf_old_old[i] /= rms + slope_rms*i;
			}
		}



		u_quartile=realbuffer[(int)(npts/4)];
		l_quartile=realbuffer[(int)(3*npts/4)];

		median=new_median;
		rms=new_rms;
//		medlow=new_medlow;


		buf_old_old=buf_old;
		buf_old=buf_new;
		binnum+=nrun;
		nrunoo=nruno;
		nruno=nrun;
		if(nrun < nrunmax) nrun = nrunstart + (int)(binnum/1000);
	}

	// finish off the last little bit
	for (int i = binnum-nrunoo-nruno; i < binnum-nruno; i++){
		amplitudes[i]/=rms;
	}

	for (int i = binnum-nruno; i < ndat; i++){
		amplitudes[i]-= median;
		amplitudes[i]/=rms;
	}

}
