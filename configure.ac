AC_INIT([boot-ana],[1.0],[francesco.sanfilippo@infn.it])

AC_CONFIG_SRCDIR([lib])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(include/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_RANLIB

AC_PROG_CC
AC_PROG_CXX
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX(17,noext,mandatory)

AX_CXXFLAGS_WARN_ALL

#silent automake
AM_SILENT_RULES([yes])

AC_CONFIG_FILES(Makefile bin/Makefile)

AC_OUTPUT
