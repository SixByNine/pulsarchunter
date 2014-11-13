dnl @synopsis SWIN_LIB_FFTWF
dnl 
AC_DEFUN([SWIN_LIB_FFTWf],
[
  AC_PROVIDE([SWIN_LIB_FFTWf])

  SWIN_PACKAGE_OPTIONS([fftw3f])

  AC_MSG_CHECKING([for single-precision FFTW-3 library])

  if test "$have_fftw3f" != "user disabled"; then

    SWIN_PACKAGE_FIND([fftw3f],[fftw3.h])
    SWIN_PACKAGE_TRY_COMPILE([fftw3f],[#include <fftw3.h>])

    if test $have_fftw3f = yes; then
      SWIN_PACKAGE_FIND([fftw3f],[libfftw3f.*])
      SWIN_PACKAGE_TRY_LINK([fftw3f],[#include <fftw3.h>],
                            [fftwf_plan_dft_1d(0,0,0,FFTW_FORWARD,FFTW_ESTIMATE);],
                            [-lfftw3f -lm])
    fi

  else
    have_fftw3f=no
  fi

  AC_MSG_RESULT([$have_fftw3f])

  if test $have_fftw3f = yes; then
    AC_DEFINE(HAVE_FFTW3_f,1,[Define if the FFTW3 single-precision library is installed])
    FFTWf_LIBS="$fftw3f_LIBS $FFTW3F_LIBS"
    FFTWf_CFLAGS="$fftw3f_CFLAGS $FFTW3F_CFLAGS"
  fi

  AC_SUBST(FFTWf_LIBS)
  AC_SUBST(FFTWf_CFLAGS)

  AM_CONDITIONAL(HAVE_FFTW3f,[test "$have_fftw3f" = yes])

])

