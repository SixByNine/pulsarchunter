#include "cbarycentre.h"
/*
UNKNOWN(0),
    JODRELL(1),
    LOVELL(2),
    MK2(3),
    FOURTYTWOFT(4),
    WARDLE(5),
    KNOCKIN(6),
    DEFFORD(7),
    TABLEY(8),
    DARNHALL(9),
    CAMBRIDGE(10),
    PARKES(11),
    ARECIBO(12),
    EFFELSBERG(13),
    GBT(14),
    NANCAY(15),
    GMRT(16);
*/ 
int get_telid(char* telname){
	if(strcmp(telname,"JODRELL")==0)return 1;
	if(strcmp(telname,"LOVELL")==0)return 2;
	if(strcmp(telname,"MK2")==0)return 3;
	if(strcmp(telname,"FOURTYTWOFT")==0)return 4;
        if(strcmp(telname,"WARDLE")==0)return 5;
        if(strcmp(telname,"KNOCKIN")==0)return 6;
        if(strcmp(telname,"DEFFORD")==0)return 7;
        if(strcmp(telname,"TABLEY")==0)return 8;
        if(strcmp(telname,"DARNHALL")==0)return 9;
        if(strcmp(telname,"CAMBRIDGE")==0)return 10;
        if(strcmp(telname,"PARKES")==0)return 11;
        if(strcmp(telname,"ARECIBO")==0)return 12;
        if(strcmp(telname,"EFFELSBERG")==0)return 13;
        if(strcmp(telname,"GBT")==0)return 14;
        if(strcmp(telname,"NANCAY")==0)return 15;
        if(strcmp(telname,"GMRT")==0)return 16;

	return 0;
}


double barycentre_doppler_factor(double epoch, int telid,double ra,double dec){
	double tobs,pobs,xma,btdb;
	barycentre(epoch,telid,ra,dec,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,&tobs,&pobs,&xma,&btdb);
	return pobs;
}

void barycentre(double EPOCH, int ITELNO, double RA20,
double DEC20, double PBEPOCH, double PB,
double PBDOT, double EPBIN, double PBIN, double ASINI, double WBIN,
double ECC, double* TOBS, double* POBS, double* XMA, double* BTDB){

	psrephb_(&EPOCH, &ITELNO, &RA20, &DEC20, &PBEPOCH, &PB, &PBDOT,
&EPBIN, &PBIN, &ASINI, &WBIN, &ECC, TOBS, POBS, XMA, BTDB);
}


