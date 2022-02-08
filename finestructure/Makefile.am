ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = finestructure

#OPTIMIZATION = -O1 -g -Wall
OPTIMIZATION = -O3 -Wall -mfpmath=sse -msse -msse2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations -I ./eigen/ -lgsl -lgslcblas
finestructure_CXXFLAGS = $(GSL_CFLAGS) $(OPTIMIZATION) -Wall 
finestructure_LDFLAGS = $(OPTIMIZATION) 
finestructure_LDADD = $(GSL_LIBS)
finestructure_SOURCES = safegetline.h data.h pcadata.h fsxml.h state.h node.h inf1.h infadmixture.h infextract.h infextract2.h infextract3.h infextract4.h infextract5.h infmcmc.h fines.h finesfunctions.h rng.h safegetline.cpp data.cpp pcadata.cpp fsxml.cpp state.cpp node.cpp inf1.cpp infadmixture.cpp infextract.cpp infextract2.cpp infextract3.cpp infextract4.cpp infextract5.cpp infmcmc.cpp finesfunctions.cpp fines.cpp rng.cpp

check_SCRIPTS = profilerun.sh
TESTS = $(check_SCRIPTS)


profilerun.sh:
	echo './finestructure -s 532772982 -x 100 -y 50 testout' > profilerun.sh
	chmod +x profilerun.sh
