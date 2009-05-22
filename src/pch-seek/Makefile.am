bin_PROGRAMS=pch-seek pch-seek-test
INCLUDES=@PSRXML_CFLAGS@ @FFTW_CFLAGS@

core_src=pch-seek.h toolkit.h TKfit.h gtools.h pch-seek-init.C pch-seek-read-file.C pch-seek-output.C pch-seek-do-search.C pch-seek-fourier.C pch-seek-phase-fit.C pch-seek-normalise.C pch-seek-harmfold.C pch-seek-search-spec.C pch-seek-recon.C  pch-seek-singlepulse.C TKfit.C gtools.C toolkit.C

pch_seek_SOURCES=pch-seek.C $(core_src)
pch_seek_LDADD=@PSRXML_LIBS@ @FFTW_LIBS@ @FFTW_CFLAGS@

pch_seek_test_SOURCES=pch-seek-test.C $(core_src)
pch_seek_test_LDADD=@PSRXML_LIBS@ @FFTW_LIBS@ @FFTW_CFLAGS@

