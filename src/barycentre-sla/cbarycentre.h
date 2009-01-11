#ifndef CBARYCENTRE_H_
#define CBARYCENTRE_H_

#ifdef __cplusplus
extern "C" {
#endif

double barycentre_doppler_factor(double epoch, int telid,double ra,double dec);

void barycentre(double EPOCH, int ITELNO, double RA20,
		double DEC20, double PBEPOCH, double PB,
		double PBDOT, double EPBIN, double PBIN, double ASINI, double WBIN,
		double ECC, double* TOBS, double* POBS, double* XMA, double* BTDB);
int get_telid(char* telname);
#ifdef __cplusplus
}
#endif

#endif
