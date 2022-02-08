ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = chromopainter

#OPTIMIZATION = -O1 -g -Wall
OPTIMIZATION = -O3 -Wall -mfpmath=sse -msse -msse2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations  -lgsl -lgslcblas
chromopainter_CXXFLAGS = $(GSL_CFLAGS) $(OPTIMIZATION) -Wall 
chromopainter_LDFLAGS = $(OPTIMIZATION) 
chromopainter_LDADD = $(GSL_LIBS)
chromopainter_SOURCES = ChromoPainterConstants.h ChromoPainterData.h ChromoPainterDonors.h ChromoPainterError.h ChromoPainterInfiles.h ChromoPainterMutEM.h ChromoPainterOutfiles.h ChromoPainterPar.h ChromoPainterReading.h ChromoPainterRecmap.h ChromoPainterSampler.h ChromoPainterData.c ChromoPainterDonors.c ChromoPainterError.c ChromoPainterInfiles.c ChromoPainterMutEM.c ChromoPainterOutfiles.c ChromoPainterPar.c ChromoPainterReading.c ChromoPainterRecmap.c ChromoPainterSampler.c ChromoPainterMain.c

check_SCRIPTS = profilerun.sh
TESTS = $(check_SCRIPTS)


profilerun.sh:
	echo './chromopainter -s 532772982 -x 100 -y 50 testout' > profilerun.sh
	chmod +x profilerun.sh
