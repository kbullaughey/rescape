#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([regevscape], [0.9], [kbullaughey@gmail.com])
AC_CONFIG_SRCDIR([array.cpp])
#AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [log])
AC_CHECK_LIB([gslcblas], [cblas_sgemm])
AC_CHECK_LIB([gsl], [gsl_ran_multinomial])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Set up OpenMP
AC_OPENMP

# Checks for library functions.
AC_FUNC_LSTAT
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([pow strtol strtoul])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
