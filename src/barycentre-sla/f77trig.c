
/* redwards 4 Feb 01 -- adding these routines since g77 doesn't
   support them as intrinsics */

#include <math.h>
#include "f772c.h"

double F772C(sind)(double x) {return sin(x);}
double F772C(cosd)(double x) {return cos(x);}
double F772C(atand)(double x) {return atan(x);}
double F772C(tand)(double x) {return tan(x);}
