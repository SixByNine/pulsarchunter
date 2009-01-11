/* $Source: /cvsroot/pulsarhunter/pulsarhunter/src/native/sla/f772c.h,v $
   $Revision: 1.1.1.1 $
   $Date: 2007/04/10 10:06:40 $
   $Author: sixbynine $ */

#ifndef F772C_H
#define F772C_H

/* redwards --- just a simple macro to transform fortran symbols into
   C symbols. NOTE --- you should also always specify the fortran
   symbol in lower case. 
   The ## is required to make the preprocessor recognize x as a token
   on its own and to the substitution */

/* For linux, if name contains '_' (e.g. sla_galeq), have to add '__'
   RNM March 2000 
   Russell found a GNU compiler flag that makes the __ unnecessary
   WVS February 2001
*/

#define F772C(x) x##_

#endif

