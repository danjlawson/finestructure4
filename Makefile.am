## CC = /opt/homebrew/opt/llvm/bin/clang
## CXX = $(CC)++

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = fs mixPainter

##OPTIMIZATION = -O3 -lm -Wall -Wno-write-strings -mfpmath=sse -msse -msse2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations -fopenmp -Wno-unused-variable $(OPENMP_CXXFLAGS)
OPTIMIZATION = -O3 -Wall -Wno-write-strings -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations -Wno-unused-variable -fopenmp $(OPENMP_CXXFLAGS)
## For mac: may need to add -Wl,-ld_classic
## For others: may need to remove them?
fs_CXXFLAGS = $(OPTIMIZATION) -Wall -Wl,-ld_classic ##  -fsanitize=leak
fs_LDFLAGS = $(OPTIMIZATION) -lm -static-libstdc++

mixPainter_CXXFLAGS = $(OPTIMIZATION) -Wall
mixPainter_LDFLAGS = $(OPTIMIZATION) -static-libstdc++

all_SOURCES = fsproject.h fsproject.cpp fssettings.h fssettings.cpp fsutils.h fsutils.cpp fsparam.h fsparam.cpp fscmds.h fscmds.cpp fsdonor.cpp fsdonor.h finestructure/fines.h alphanum.hpp fsconstants.h cp/ChromoPainterReading.h cp/ChromoPainterConstants.h cp/ChromoPainterError.h cp/ChromoPainterData.h  cp/ChromoPainterMutEM.h cp/ChromoPainterPar.h cp/ChromoPainterRecmap.h cp/ChromoPainterDonors.h cp/ChromoPainterOutfiles.h cp/ChromoPainterInfiles.h cp/ChromoPainterSampler.h cp/ChromoPainterReading.c cp/ChromoPainterError.c cp/ChromoPainterData.c  cp/ChromoPainterMutEM.c cp/ChromoPainterPar.c cp/ChromoPainterRecmap.c cp/ChromoPainterDonors.c cp/ChromoPainterOutfiles.c cp/ChromoPainterInfiles.c cp/ChromoPainterSampler.c  cp/adler32.c cp/compress.c cp/crc32.c cp/crc32.h cp/deflate.c cp/deflate.h cp/gzclose.c cp/gzguts.h cp/gzlib.c cp/gzread.c cp/gzwrite.c cp/infback.c cp/inffast.c cp/inffast.h cp/inffixed.h cp/inflate.c cp/inflate.h cp/inftrees.c cp/inftrees.h cp/trees.c cp/trees.h cp/uncompr.c cp/zconf.h cp/zlib.h cp/zutil.c cp/zutil.h finestructure/safegetline.h finestructure/data.h finestructure/fsxml.h finestructure/state.h finestructure/prior.h finestructure/node.h finestructure/inf1.h finestructure/infadmixture.h finestructure/infconcordance.h  finestructure/infextract.h finestructure/infextract2.h finestructure/infextract3.h finestructure/infextract4.h finestructure/infextract5.h finestructure/infextractdonors.h finestructure/infmcmc.h finestructure/fines.h finestructure/finesfunctions.h finestructure/rng.h finestructure/safegetline.cpp finestructure/data.cpp finestructure/fsxml.cpp finestructure/state.cpp finestructure/prior.cpp finestructure/node.cpp finestructure/inf1.cpp finestructure/infadmixture.cpp finestructure/infconcordance.cpp finestructure/infextract.cpp finestructure/infextract2.cpp finestructure/infextract3.cpp finestructure/infextract4.cpp finestructure/infextract5.cpp finestructure/infextractdonors.cpp finestructure/infmcmc.cpp finestructure/finesfunctions.cpp finestructure/fines.cpp finestructure/rng.cpp chromocombine/ChromoCombine.cpp chromocombine/ChromoCombineFuns.cpp chromocombine/ChromoCombineFuns.h chromocombine/ChromoCombine.h
fs_SOURCES = fs.cpp $(all_SOURCES)
mixPainter_SOURCES = badpainter.cpp $(all_SOURCES)

dist_noinst_SCRIPTS = scripts/finestructuregreedy.sh scripts/qsub_array.sh scripts/sbatch_array.sh scripts/makeuniformrecfile.pl scripts/convertrecfile.pl scripts/idsFromChunkcounts.sh scripts/phasescreen.pl scripts/phasesubsample.pl scripts/impute2chromopainter.pl scripts/beagle2chromopainter.pl scripts/msms2cp.pl Manual/manual.pdf examples/example1 examples/example2 examples/example3 examples/example4 LICENSE-FS LICENSE-CP
