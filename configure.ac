AC_INIT([fs], [4.1.0], [dan.lawson@bristol.ac.uk],
             [fs], [http://www.paintmychromosomes.com])
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([foreign 1.10 -Wall no-define])
AC_CONFIG_HEADERS([config.h])

AC_MSG_NOTICE([Looking for a GNU C++ compiler])

case "$host_os" in
darwin*)

    AC_MSG_NOTICE([Searching for an OpenMP-capable GNU compiler])

    if test -x "$srcdir/whereisglob"; then

        for c in `$srcdir/whereisglob 'g++*' | sort -V -r`; do
            AC_MSG_NOTICE([Trying g++ compiler $c])

            if echo 'int main(){return 0;}' |
                $c -fopenmp -x c++ - -o conftest >/dev/null 2>&1
            then
                CXX="$c"
                AC_MSG_NOTICE([Using C++ compiler: $c])
                break
            fi
        done

        for c in `$srcdir/whereisglob 'gcc*' | sort -V -r`; do
            AC_MSG_NOTICE([Trying gcc compiler $c])

            if echo 'int main(){return 0;}' |
                $c -fopenmp -x c - -o conftest >/dev/null 2>&1
            then
                CC="$c"
                AC_MSG_NOTICE([Using C compiler: $c])
                break
            fi
        done
    fi
    ;;
esac

AC_PROG_CXX
AC_PROG_CC # Perhaps causes errors in some places, but should be fine
CFLAGS="$CFLAGS -std=gnu11"

LT_INIT
AC_CONFIG_MACRO_DIR([m4])
#AM_PROG_CC_C_O
#AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# Checks for libraries.
# find gsl, wherever it's hiding on the system
#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
#AX_PATH_GSL(,,AC_MSG_ERROR(could not find required version of GSL))
#AC_CHECK_LIB([z],[gzread],,AC_MSG_ERROR(Could not find zlib))

AC_OPENMP
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")

# The Apple linker (Xcode 15+) needs -ld_classic; GNU ld rejects it (it reads
# -ld_classic as a request to link libd_classic and fails to find it). Set the
# flag only on macOS so non-Apple hosts link out of the box.

AC_OPENMP(
  [AC_DEFINE(HAVE_OPENMP,1,[OpenMP available])],
  [AC_MSG_WARN([OpenMP not found; building single-threaded version])]
)

case "$host_os" in
  darwin*) LD_CLASSIC="-Wl,-ld_classic" ;;
  *)       LD_CLASSIC="" ;;
esac
AC_SUBST([LD_CLASSIC])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
