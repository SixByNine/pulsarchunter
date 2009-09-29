#include <stdio.h>
void pch_seek_fourier_size_test_slack(int npts, int val, int *bestsize, int *bestslack){
	int slack = npts - val;
	if (slack < *bestslack){
		*bestsize=val;
		*bestslack=slack;
	}
}

int pch_seek_fourier_size(int npts, bool pow2only){

	int bestslack=npts;
	int bestsize;
	int val;
	if(!pow2only){
		// try some things
		// 3
		val = 3*pch_seek_fourier_size(npts/3,true);
		pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
		// 3^2
		val = 9*pch_seek_fourier_size(npts/9,true);
		pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
		// 3^3
		val = 27*pch_seek_fourier_size(npts/27,true);
		pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
		// 3^4
		val = 81*pch_seek_fourier_size(npts/81,true);
		pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
		// 3^5
		val = 243*pch_seek_fourier_size(npts/243,true);
		pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
		// 3^6
		val = 729*pch_seek_fourier_size(npts/729,true);
		pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
	}
	val = 128;
	while ( val <= npts ){
		val *=2 ;
	}
	val /=2;
	pch_seek_fourier_size_test_slack(npts,val,&bestsize,&bestslack);
	return bestsize;
}


