#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include "pch-seek.h"



void pch_seek_zap_sigproc(char* zapfilename, float* spectrum,int npts,float xscale){
	FILE* zapfile;
	char line[1024];
	float blo,bhi,bf,df;
	int nh_start,nh,nzap,count;
	char stat; // static widths


	zapfile = fopen(zapfilename,"r");
	if (!zapfile){
		fprintf(stderr,"Error: Cannot open zapfile '%s'\n",zapfilename);
		return;
	}

	count=0;
	while(!feof(zapfile)){
		fgets(line,1024,zapfile);
		count++;
		if (line[0]=='#' || line[0]=='\n')continue;
		int nentries = sscanf(line,"%f %f %d %d",&blo,&bhi,&nh_start,&nh);
		float blo_o=blo;
		float bhi_o=bhi;
		if(nentries == 3){
			// old 3 number format
			nh=nh_start;
			nh_start=1;
			nentries=4;
		}
		if(nentries!=4){
			printf("Bad line (#%d) in zap file!",count);
			continue;
		}

		bf=0.5*(bhi+blo);
		df=0.5*(bhi-blo);
		if (nh < 0 ){
			stat=1;
			nh=-nh;
		} else stat=0;
		nzap=0;
		for (int i=nh_start; i <= nh; i++){
			if(stat){
				blo=bf*(float)i-df;
				bhi=bf*(float)i+df;
			} else {
				blo=(bf-df)*(float)i;
				bhi=(bf+df)*(float)i;
			}
			//convert freq to bin number.
			int nlo=(int)(blo/xscale+0.5);
			int nhi=(int)(bhi/xscale+0.5);
			// we could be zapping aliased harmonics
			if (nlo > npts){
				nlo = npts-nlo;
			}
			if (nhi > npts){
				nhi = npts-nhi;
			}
			if (nlo > nhi){
				int dd=nlo;
				nlo=nhi;
				nhi=dd;
			}

			// now do the zapping
			for(int j=nlo; j < nhi; j++){
				spectrum[j]=0;
				nzap++;
			}
		}
		// Now we are ready to start zapping...
		printf("Filter: %f -> %f Hz %d->%d harmonics (zeroed %d bins)\n"
				,blo_o,bhi_o,nh,nh_start,nzap);



	}
	fclose(zapfile);
}
