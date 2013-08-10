dnl -*- m4 -*-
dnl ============================================================================================
dnl
dnl	This file is part of LiDIA --- a library for computational number theory
dnl
dnl	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
dnl
dnl	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
dnl
dnl --------------------------------------------------------------------------------------------
dnl
dnl	$Id$
dnl
dnl	Author	: Safuat Hamdy
dnl	Changes	: see CVS log
dnl
dnl ============================================================================================


dnl
dnl Process this file with autoconf-2.52 (or later) to produce a
dnl configure script.
dnl


AC_INIT([LiDIA],[2.3.0],[lidia@cdc.informatik.tu-darmstadt.de])
AC_REVISION($Id$)
AC_PREREQ(2.52)
AC_CONFIG_SRCDIR(.LIDIA_SOURCE_DIR)
AM_INIT_AUTOMAKE([tar-ustar])

# Update version info in some files for later distribution
AC_CONFIG_FILES([LiDIA.spec])

# This disables the use of Autoheader for updating *.h.in
AUTOHEADER=:

# <config.h> holds information about available headers, types, functions, etc.
# <LiDIA/LiDIA.h> holds information relevant for binary compatibility issues.
# <LiDIA/path.h> holds database installation paths.
dnl Due to a bug in Automake 1.5, all headers must be declared in one
dnl single AM_CONFIG_HEADER, within one source line.
AM_CONFIG_HEADER(include/config.h:config.h.in include/LiDIA/LiDIA.h:LiDIA.h.in include/LiDIA/path.h:path.h.in)

AC_LANG(C++)

# Some definitions for <LiDIA/LiDIA.h>
AC_DEFINE(LIDIA_MAJOR_VERSION, 2)
AC_DEFINE(LIDIA_MINOR_VERSION, 3)
AC_DEFINE(LIDIA_PATCH_VERSION, 0)

# Automake defines `pkgdatadir' in the generated `Makefile.in's (though it
# does not set a corresponding shell variable here).
# Expand the (to-be-assigned) value of `pkgdatadir' recursively
# and use the result for defining the installation paths.
pkgdatadir_x='${datadir}/${PACKAGE}'
LIDIA_EXPAND_DIRVAR(pkgdatadir_x)

# Definitions for <LiDIA/path.h>
AC_DEFINE_UNQUOTED([LIDIA_GF2N_DB],"${pkgdatadir_x}/GF2n.database",
	[The filename of LiDIA's GF(2^n) database])
AC_DEFINE_UNQUOTED([LIDIA_PRIMELIST_NAME],"${pkgdatadir_x}/LIDIA_PRIMES",
	[The filename of LiDIA's primes list for remainder arithmetic])
AC_DEFINE_UNQUOTED([LIDIA_ECO_DB],"${pkgdatadir_x}/MOD_EQ",
	[The directory holding LiDIA ECO's GF(p) database files])
AC_DEFINE_UNQUOTED([LIDIA_ECO2_DB],"${pkgdatadir_x}/MOD_EQ_2",
	[The directory holding LiDIA ECO's GF(2^n) database files])
AC_DEFINE_UNQUOTED([LIDIA_DISCRIMINANTS_DB],"${pkgdatadir_x}/Discriminants",
	[The directory holding LiDIA GEC's imaginary quadratic discriminant database files])
AC_DEFINE_UNQUOTED([LIDIA_PRIME_DIFFS],"${pkgdatadir_x}/LIDIA_PRIME_DIFFS",
	[The file holding the differences of subsequential primes used for LiDIA GEC's prime_proof])



#
#	CHECKING FOR PROGRAMS
#

dnl AC_PROG_MAKE_SET	-- done by AM_INIT_AUTOMAKE
dnl AC_PROG_INSTALL	-- done by AM_INIT_AUTOMAKE

dnl  after AC_PROG_CC and AC_PROG_CXX we can't say anymore whether
dnl  the user set CGLAGS and CXXFLAGS, respectively, outside of
dnl  configure. We therefore store the original values in
dnl  ac_USER_DEFINED_CFLAGS and ac_USER_DEFINED_CXXFLAGS

ac_USER_DEFINED_CFLAGS="${CFLAGS}"
test -z "${CFLAGS}" && CFLAGS="-O2"
AC_PROG_CC
LIDIA_CHECK_CC_FLAGS

ac_USER_DEFINED_CXXFLAGS="${CXXFLAGS}"
test -z "${CXXFLAGS}" && CXXFLAGS="-O2"
AC_PROG_CXX
LIDIA_CHECK_CXX_FLAGS
AC_PROG_CXXCPP


#
#	check for some ISO C++ features --- they're absolutely mandatory!
#

LIDIA_CHECK_CXX_ISO_BOOL
LIDIA_CHECK_CXX_ISO_INLINE
LIDIA_CHECK_CXX_ISO_CAST
LIDIA_CHECK_CXX_ISO_EXPLICIT
LIDIA_CHECK_CXX_ISO_MUTABLE
LIDIA_CHECK_CXX_ISO_TEMPLATE_SPEC
LIDIA_CHECK_CXX_ISO_TEMPLATE_INST
dnl LIDIA_CHECK_CXX_ISO_TEMPLATE_EXPORT

AC_ARG_ENABLE(namespaces,
[  --enable-namespaces     use namespaces (default = yes)],
[
	if test "$enable_namespaces" = "yes"; then
		LIDIA_CHECK_CXX_ISO_NAMESPACES
	fi
],
[
	LIDIA_CHECK_CXX_ISO_NAMESPACES
])


dnl AC_PROG_RANLIB	-- done by LIDIA_PROG_LIBTOOL
AC_DISABLE_SHARED
LIDIA_PROG_LIBTOOL


LIDIA_CHECK_PROG_LATEX
LIDIA_CHECK_PROG_GUNZIP
LIDIA_CHECK_PROG_BUNZIP2



AC_ARG_WITH(extra-includes,
[  --with-extra-includes=DIR
                          adds non standard include path],
[
	EXTRA_INCLUDES="${EXTRA_INCLUDES} -I${with_extra_includes}"
],
)


AC_ARG_WITH(extra-libs,
[  --with-extra-libs=DIR   adds non standard library path],
[
	LIBS="${LIBS} -L${with_extra_libs}"
],
)


#
#	CHECKING FOR LIBRARIES
#

AC_ARG_WITH(arithmetic,
[  --with-arithmetic=ARG   select multi-precision arithmetic:
                          gmp: GNU MP (default)
                          cln
                          piologie
                          libI],
[
	case $with_arithmetic in
	gmp)
		LIDIA_CHECK_LIB_GMP
		if test "$lidia_cv_gmp" = "yes"; then
			lidia_arithmetic=gmp
			LIBS="${LIBS} -lgmp"
			AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_GMP)
		else
			AC_MSG_ERROR([
you don't appear to have GNU MP.  Either add --with-extra-includes and --with-extra-libs
or get a copy of the last release of GNU MP from http://www.swox.com/gmp/])
		fi
		;;
	cln)
		LIDIA_CHECK_LIB_CLN
		if test "$lidia_cv_cln" = "yes"; then
			lidia_arithmetic=cln
			LIBS="${LIBS} -lcln"
			AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_CLN)
		else
			AC_MSG_ERROR([
you don't appear to have CLN.  Either add --with-extra-includes and --with-extra-libs
or get a copy of the last release of CLN from http://clisp.cons.org/~haible/])
		fi
		;;
	piologie)
		LIDIA_CHECK_LIB_PIOLOGIE
		if test "$lidia_cv_piologie" = "yes"; then
			lidia_arithmetic=piologie
			LIBS="${LIBS} -lpiologie"
			AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_PIOLOGIE)
		else
			AC_MSG_ERROR([
you don't appear to have Piologie.  Either add --with-extra-includes and --with-extra-libs
or get a copy of the last release of Piologie from http://www.hipilib.de/piologie.htm])
		fi
		;;
	libI)
		LIDIA_CHECK_LIB_LIBI
		if test "$lidia_cv_libI" = "yes"; then
			lidia_arithmetic=libI
			LIBS="${LIBS} -lI"
			AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_LIBI)
		fi
		;;
	*)
		AC_MSG_ERROR([Invalid value for --with-arithmetic])
		;;
	esac
],
[
	AC_MSG_NOTICE(checking which multi-precision arithmetic library to use)
	LIDIA_CHECK_LIB_GMP
	if test "$lidia_cv_gmp" = "yes"; then
		lidia_arithmetic="gmp"
		LIBS="${LIBS} -lgmp"
		AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_GMP)
	else
	LIDIA_CHECK_LIB_CLN
	if test "$lidia_cv_cln" = "yes"; then
		lidia_arithmetic="cln"
		LIBS="${LIBS} -lcln"
		AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_CLN)
	else
	LIDIA_CHECK_LIB_PIOLOGIE
	if test "$lidia_cv_piologie" = "yes"; then
		lidia_arithmetic="piologie"
		LIBS="${LIBS} -lpiologie"
		AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_PIOLOGIE)
	else
	LIDIA_CHECK_LIB_LIBI
	if test "$lidia_cv_libI" = "yes"; then
		lidia_arithmetic="libI"
		LIBS="${LIBS} -lI"
		AC_DEFINE(LIDIA_ARITHMETIC, LIDIA_ARITH_LIBI)
	else
	AC_MSG_ERROR(could not find a multi-precision integer arithmetic.)
	fi fi fi fi
])
dnl AC_SUBST(lidia_arithmetic)	-- no longer needed outside configure



#
#	CHECKING FOR HEADER FILES
#

AC_HEADER_STDC

# check for some POSIX headers
AC_CHECK_HEADERS(sys/times.h sys/stat.h unistd.h signal.h,,
	AC_MSG_WARN([
It seems that you don't have some required POSIX C headers.
])
)

# check for some C++ style headers
AC_CHECK_HEADERS(iostream cstddef cstdio cstdlib cmath cctype cstring ctime cassert climits,,
	AC_MSG_ERROR([
It seems that you don't have some required ISO C++ style headers.
])
)

# check for new style string streams
AC_CHECK_HEADERS(sstream,,
	[AC_CHECK_HEADERS(strstream,,
		AC_MSG_ERROR([
Your system lacks both <sstream> and <strstream> headers.
At least one of these standard C++ headers is required.
])
	)
	AC_MSG_WARN([
Your system lacks the standard C++ header <sstream>.
LiDIA will be built using the deprecated char*-Streams in <strstream>.
])
	AC_DEFINE(OSSTREAM_USE_STRSTREAM,1)]
)


#
#	CHECKING FOR TYPEDEFS
#



#
#	CHECKING FOR STRUCTURES
#



#
#	CHECKING FOR COMPILER CHARACTERISTICS
#

LIDIA_CHECK_CXX_TYPES
LIDIA_CHECK_GCC_INSTANTIATION



#
#	CHECKING FOR LIBRARY FUNCTIONS
#

LIDIA_CHECK_STLPORT
LIDIA_CHECK_FUNC_RANDOM
LIDIA_CHECK_FUNC_MKSTEMP



#
#	CHECKING FOR SYSTEM SERVICES
#

LIDIA_CHECK_SYS_POSIX_SIGNALS
LIDIA_CHECK_SYS_POSIX_TIMES
LIDIA_CHECK_SYS_POSIX_TIME



#
#	CHECKING FOR FEATURES
#

dnl do not cache the enable-variables!


AC_ARG_ENABLE(inline,
[  --enable-inline         inline multiprecision arithmetic functions (default = yes)],
[
	if test "$enable_inline" = "yes"; then
		AC_DEFINE(INLINE_INTERFACE, 1)
	else
		AC_DEFINE(INLINE_INTERFACE, 0)
	fi
],
[
	AC_DEFINE(INLINE_INTERFACE, 1)
])


dnl Check whether the user wants to disable exceptions.
dnl If we know how to tell the compiler to do so, build the library
dnl without code for exception handling.
AC_ARG_ENABLE(exceptions,
[  --enable-exceptions     compile with error handling via exceptions activated (default = yes)],
[
   if test "$enable_exceptions" = "yes"; then
     AC_DEFINE(LIDIA_EXCEPTIONS, 1)
   else
     LIDIA_DISABLE_EXCEPTIONS
   fi
],
[
     AC_DEFINE(LIDIA_EXCEPTIONS, 1)
])


AC_ARG_ENABLE(assert,
[  --enable-assert         compile with assert macros activated (default = no)],
[
	if test "$enable_assert" != "yes" ; then
		AC_DEFINE(NDEBUG, 1)
	fi
],
[
	AC_DEFINE(NDEBUG, 1)
])


lidia_base="yes"


AC_ARG_ENABLE(ff,
[  --enable-ff             compile the finite fields package (default = yes)],
[
	if test "$enable_ff" = "yes" ; then
		lidia_ff="yes"
	else
		lidia_ff="no"
	fi
],
[
	lidia_ff="yes"
])


AC_ARG_ENABLE(la,
[  --enable-la             compile the linear algebra package (default = yes)],
[
	if test "$enable_la" = "yes" ; then
		lidia_la="yes"
	else
		lidia_la="no"
	fi
],
[
	lidia_la="yes"
])


AC_ARG_ENABLE(lt,
[  --enable-lt             compile the lattices package (default = yes)],
[
	if test "$enable_lt" = "yes" ; then
		lidia_lt="yes"
	else
		lidia_lt="no"
	fi
],
[
	lidia_lt="yes"
])


AC_ARG_ENABLE(nf,
[  --enable-nf             compile the number fields package (default = yes)],
[
	if test "$enable_nf" = "yes" ; then
		lidia_nf="yes"
	else
		lidia_nf="no"
	fi
],
[
	lidia_nf="yes"
])


AC_ARG_ENABLE(ec,
[  --enable-ec             compile the elliptic curves package (default = yes)],
[
	if test "$enable_ec" = "yes" ; then
		lidia_ec="yes"
	else
		lidia_ec="no"
	fi
],
[
	lidia_ec="yes"
])


AC_ARG_ENABLE(eco,
[  --enable-eco            compile the elliptic curve order package (default = yes)],
[
	if test "$enable_eco" = "yes" ; then
		lidia_eco="yes"
	else
		lidia_eco="no"
	fi
],
[
	lidia_eco="yes"
])


AC_ARG_ENABLE(gec,
[  --enable-gec            compile the elliptic curve generation package (default = yes)],
[
	if test "$enable_gec" = "yes" ; then
		lidia_gec="yes"
	else
		lidia_gec="no"
	fi
],
[
	lidia_gec="yes"
])


LIDIA_LTCLIBS="\${LIDIA_BASE_LTCLIB}"
PKG_DIRS="base"
AC_CONFIG_FILES(library/base/Makefile examples/base/Makefile)


AC_MSG_CHECKING(whether to compile the finite fields package)
if test "$lidia_ff" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/finite_fields ; then
	AC_MSG_RESULT(no (not present))
	lidia_ff="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_FF_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} finite_fields"
	AC_CONFIG_FILES(library/finite_fields/Makefile examples/finite_fields/Makefile)
fi

AC_MSG_CHECKING(whether to compile the linear algebra package)
if test "$lidia_la" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/linear_algebra ; then
	AC_MSG_RESULT(no (not present))
	lidia_la="no"
elif test "$lidia_ff" != "yes"; then
	AC_MSG_RESULT(no (requires finite fields package))
	lidia_la="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_LA_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} linear_algebra"
	AC_CONFIG_FILES(library/linear_algebra/Makefile examples/linear_algebra/Makefile)
fi


AC_MSG_CHECKING(whether to compile the lattice package)
if test "$lidia_lt" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/lattices ; then
	AC_MSG_RESULT(no (not present))
	lidia_lt="no"
elif test "$lidia_la" != "yes"; then
	AC_MSG_RESULT(no (requires linear algebra package))
	lidia_lt="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_LT_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} lattices"
	AC_CONFIG_FILES(library/lattices/Makefile examples/lattices/Makefile)
fi


AC_MSG_CHECKING(whether to compile the number fields package)
if test "$lidia_nf" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/number_fields ; then
	AC_MSG_RESULT(no (not present))
	lidia_nf="no"
elif test "$lidia_lt" != "yes"; then
	AC_MSG_RESULT(no (requires lattice package))
	lidia_nf="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_NF_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} number_fields"
	AC_CONFIG_FILES(library/number_fields/Makefile examples/number_fields/Makefile)
fi


AC_MSG_CHECKING(whether to compile the elliptic curve package)
if test "$lidia_ec" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/elliptic_curves ; then
	AC_MSG_RESULT(no (not present))
	lidia_ec="no"
elif test "$lidia_la" != "yes"; then
	AC_MSG_RESULT(no (requires linear algebra package))
	lidia_ec="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_EC_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} elliptic_curves"
	AC_CONFIG_FILES(library/elliptic_curves/Makefile examples/elliptic_curves/Makefile)
fi


AC_MSG_CHECKING(whether to compile the elliptic curve order package)
if test "$lidia_eco" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/ec_order ; then
	AC_MSG_RESULT(no (not present))
	lidia_eco="no"
elif test "$lidia_ec" != "yes"; then
	AC_MSG_RESULT(no (requires elliptic curves package))
	lidia_eco="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_ECO_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} ec_order"
	AC_CONFIG_FILES(library/ec_order/Makefile examples/ec_order/Makefile)
fi


AC_MSG_CHECKING(whether to compile the elliptic curve generation package)
if test "$lidia_gec" != "yes"; then
	AC_MSG_RESULT(no)
elif test ! -d $srcdir/src/gec ; then
	AC_MSG_RESULT(no (not present))
	lidia_gec="no"
elif test "$lidia_eco" != "yes"; then
	AC_MSG_RESULT(no (requires elliptic curve order package))
	lidia_gec="no"
elif test "$lidia_nf" != "yes"; then
	AC_MSG_RESULT(no (requires number field package))
	lidia_gec="no"
else
	AC_MSG_RESULT(yes)
	LIDIA_LTCLIBS="\${LIDIA_GEC_LTCLIB} ${LIDIA_LTCLIBS}"
	PKG_DIRS="${PKG_DIRS} gec"
	AC_CONFIG_FILES(library/gec/Makefile examples/gec/Makefile)
fi


AC_SUBST(LIDIA_LTCLIBS)
AC_SUBST(PKG_DIRS)
test -z "${EXTRA_INCLUDES}" || CPPFLAGS="${CPPFLAGS} ${EXTRA_INCLUDES}"



#
#	Prepare to create links according to selected multi precision library
#

AC_CONFIG_LINKS(src/base/interfaces/bigint/bigint.cc:src/interfaces/bigint/${lidia_arithmetic}/bigint.cc)
AC_CONFIG_LINKS(src/base/interfaces/udigit/udigit.cc:src/interfaces/udigit/${lidia_arithmetic}/udigit.cc)
AC_CONFIG_LINKS(src/base/interfaces/udigit/udigit_appl.cc:src/interfaces/udigit/${lidia_arithmetic}/udigit_appl.cc)
AC_CONFIG_LINKS(include/LiDIA/kernel/bigint_def.h:src/interfaces/bigint/${lidia_arithmetic}/bigint_def.h)
AC_CONFIG_LINKS(include/LiDIA/kernel/bigint_interface.h:src/interfaces/bigint/${lidia_arithmetic}/bigint_interface.h)
AC_CONFIG_LINKS(include/LiDIA/kernel/udigit_def.h:src/interfaces/udigit/${lidia_arithmetic}/udigit_def.h)
AC_CONFIG_LINKS(include/LiDIA/kernel/udigit_interface.h:src/interfaces/udigit/${lidia_arithmetic}/udigit_interface.h)

if test "$lidia_cv_func_random" != "yes" ; then
	AC_CONFIG_LINKS(include/LiDIA/random.h:src/portability/random.h)
fi



#
#	create output files
#

AC_CONFIG_FILES(Makefile library/Makefile examples/Makefile doc/Makefile)
AC_OUTPUT
