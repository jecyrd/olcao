AC_INIT([OLCAO],[0.6],[rulisp@umkc.edu])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])


dnl Link to the appropriate library.  The dtsevx and zherk are just test cases
dnl   to make sure the libraries exist.  If they are found, then they are
dnl   included in the link stage.
dnl AC_CHECK_LIB(lapack, dstevx)
dnl AC_CHECK_LIB([blas], [zherk])
dnl AC_CHECK_LIB([veclib], [zherk])
dnl AC_CHECK_LIB([xt-libsci], [zherk])


dnl ######################################################################
dnl
dnl Require the user to specify the location of the lapack and blas
dnl
dnl ######################################################################


dnl ######################################################################
dnl
dnl Confirm that we have access to all the necessary files.
dnl
dnl ######################################################################
dnl AC_CHECK_FILE(${HDF5_LIBPATH}/libhdf5_fortran.so,LIBHDF5_A=${HDF5_LIBPATH}/hdf5_fortran.h,LIBHDF5_A="")
dnl if test -n "LIBHDF5_A"; then
dnl    ac_cv_have_hdf5=yes
dnl else
dnl    ac_cv_have_hdf5=no
dnl fi

dnl ######################################################################
dnl
dnl Allow perl modules to be installed as desired.
dnl
dnl ######################################################################

AC_ARG_ENABLE([perllibdir],
   [AS_HELP_STRING([--enable-perllibdir=DIR],
   [perl module install directory ($(bindir)/perl5)])],
   [perllibdir=$enableval],
   [perllibdir='$(bindir)/perl5'])
AC_SUBST([perllibdir])

dnl ######################################################################
dnl
dnl Set up the fortran compile flags.
dnl
dnl ######################################################################

AC_ARG_VAR(FC,Fortran compiler)
AC_ARG_VAR(FCFLAGS,Fortran compiler flags)


dnl ######################################################################
dnl
dnl Define for whether MPI found
dnl
dnl ######################################################################

AC_MSG_NOTICE
AC_MSG_NOTICE([Checks for MPI])
AC_MSG_NOTICE
AX_MPI()
if test "$HAVE_MPI" = ""; then
  AC_MSG_ERROR([Unable to find MPI installation])
fi

dnl ######################################################################
dnl
dnl Define for whether hdf5 fortran found.
dnl Look at macro for explanation of decisions made below.
dnl
dnl ######################################################################

AC_MSG_NOTICE
AC_MSG_NOTICE([Checks for HDF5])
AC_MSG_NOTICE
AX_LIB_HDF5([parallel])
if test "$with_hdf5_fortran" = "yes"; then
  FC="$HDF5_FC"
else
  AC_MSG_ERROR([Unable to find parallel HDF5, required for olcao compilation])
fi

dnl ######################################################################
dnl
dnl Set up correct use of the pre-processor (C or Fortran).
dnl
dnl ######################################################################

dnl if test -z "$CPP"; then
dnl    AC_PATH_PROG(CPP,fpp)
dnl else
   AC_PROG_CPP
dnl fi

dnl if test -z "$CPPFLAGS"; then
dnl   CPPFLAGS="-Wp,-C,-P"
dnl fi

AC_SUBST(CPP)
dnl AC_SUBST(CPPFLAGS)



dnl ######################################################################
dnl
dnl All Done!   Output the results.
dnl
dnl ######################################################################
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
