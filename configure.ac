AC_INIT([fs], [3.0.0], [dan.lawson@bristol.ac.uk],
             [fs], [http://www.paintmychromosomes.com])
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([foreign 1.10 -Wall no-define])
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
AC_PROG_CXX
AC_PROG_CC # Perhaps causes errors in some places, but should be fine

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

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
