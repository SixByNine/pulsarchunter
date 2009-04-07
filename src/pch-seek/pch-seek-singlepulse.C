#include "gtools.h"
#include <vector>
#include <stdio.h>

void pch_seek_singlepulse(float* timeseries, int npts, float thresh, float dm, int max_scrunch, int beam, char* filename){
	FILE* file;
	vector<Gpulse> pulses = findgiants(npts, timeseries, thresh, 30.0, max_scrunch, dm, 1);
	file = fopen(filename,"a");
	
	std::vector<Gpulse>::iterator itr;
	for ( itr = pulses.begin(); itr != pulses.end(); ++itr ){
		//
		//filename   amplitude    SNR   startbin   peakbin   width    tscrunch   dm    beam
		fprintf(file,"%s\t%f\t%f\t%d\t%d\t%d\t%d\t%f\t%d\n","xx",itr->amp,itr->SNR, itr->start, itr->width, itr->tscrfac,itr->dm,beam);
	}

	fclose(file);
}
