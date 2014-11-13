#include <stdio.h>
void pch_seek_fourier_size_test_slack(int npts, int val, int *bestsize, int *bestslack){
	int slack = npts - val;
	if (slack < 0)slack=npts;
	if (slack < *bestslack){
		*bestsize=val;
		*bestslack=slack;
	}
}

int pch_seek_fourier_size(int npts, bool pow2only){

	int bestslack=npts;
	int bestsize;
	int val,x,y;
	if(!pow2only){
		// try some things
		// powers of 3
		x=1;
		while (x < npts){
			y=1;
			//powers of 5
			while (y<npts){
				val = y*x*pch_seek_fourier_size(npts/(x*y),true);
				pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
				y *= 5;
			}
			x*=3;
		}
	}
	val = 128;
	while ( val <= npts ){
		val *=2 ;
	}
	val /=2;
	pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
	return bestsize;
}


