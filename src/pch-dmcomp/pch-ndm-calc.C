// Creates a list of accelerations for each DM 
// Lina Levin
// 2008-06-25
// Modified by Mkeith for a nefarious purpose! Don't use it if you want the origianal functionality!

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

double pch_getDMtable_ndm(float DMstart, int max_ndms, double tsamp, double ti, double bw, double cfreq, int Nchan, double tol){

   int Ndms=0;
   float DM = DMstart;
   int j;

   cout << "DMstart is " << DMstart << " Ndms is " << Ndms << endl;

   cout << "tsamp is " << tsamp << endl;
   cout << "ti is " << ti << endl;
   cout << "bw is " << bw << endl;
   cout << "cfreq is " << cfreq << endl;
   cout << "Nchan is " << Nchan << endl;
   cout << "tol is " << tol << endl;

   while(Ndms < max_ndms){
	  double oldDM = DM;
	  double t00 = sqrt(tsamp*tsamp + ti*ti + pow(8.3*bw*DM/pow(cfreq,3.0),2));
	  double a = 8.3*bw/pow(cfreq,3);
	  double b = 8.3*Nchan*bw/4/pow(cfreq,3);
	  double c = tol*tol*t00*t00 - tsamp*tsamp - ti*ti;
	  double newDM = (b*b*DM + sqrt(-a*a*b*b*DM*DM + a*a*c + b*b*c))/(a*a+b*b);
	  //DM = oldDM + 2.0*(newDM-oldDM);
	  DM = newDM;
	  Ndms++;
   }
   return DM;


}

