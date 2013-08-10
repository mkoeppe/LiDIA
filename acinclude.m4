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
dnl	Macro to expand directory variables recursively.
dnl	Example usage:
dnl		var='${datadir}/${PACKAGE}'
dnl		LIDIA_EXPAND_DIRVAR(var)
dnl
AC_DEFUN([LIDIA_EXPAND_DIRVAR],
[
# Caution: prefix=NONE exec_prefix=NONE are not resolved before AC_OUTPUT.
restore_NONE=""
case $prefix in NONE)
	prefix=$ac_default_prefix
	restore_NONE="$restore_NONE prefix=NONE"
	;;
esac
case $exec_prefix in NONE)
	exec_prefix='${prefix}'
	restore_NONE="$restore_NONE exec_prefix=NONE"
	;;
esac
while :; do
	case $$1 in
	*\${*}*)
		eval "lidia_expand_dirvar_tmp=\"$$1\"" 2>/dev/null || break
		case $lidia_expand_dirvar_tmp in
		"" | "$$1") break ;;
		*) $1=$lidia_expand_dirvar_tmp ;;
		esac ;;
	*) break ;;
	esac
done
eval "$restore_NONE"
])


dnl
dnl
dnl
AC_DEFUN([LIDIA_CHECK_CXX_TYPES],
[
	AC_CHECK_SIZEOF(short, 2)
	AC_CHECK_SIZEOF(int, 4)
	AC_CHECK_SIZEOF(long, 4)
	AC_CHECK_SIZEOF(double, 8)

	AC_C_BIGENDIAN
])



dnl
dnl	find good optimization flags for target
dnl
dnl	for the time being, set OPTFLAGS simply to -O2
dnl	though some compilers accept -fast (SUNpro CC, HP aCC) or -Ofast (MIPSpro CC)
dnl	for maximum performance
dnl
AC_DEFUN([LIDIA_CHECK_CC_FLAGS],
[
  if test "${GCC}" != "yes" && test -z "${ac_USER_DEFINED_CFLAGS}"; then
    AC_MSG_WARN([
Currently, gcc is the only C compiler for which we know
how to set the compiler switches, but we found ${CC}.
We will assume CFLAGS = ${CFLAGS} and go on, but we strongly
recommend you set CFLAGS on the command line.
    ])
  fi
])


AC_DEFUN([LIDIA_CHECK_CXX_FLAGS],
[
  if test "$GXX" = "yes"; then
    case ${CXXFLAGS} in
      *)
        ;;
    esac
  elif test -z "${ac_USER_DEFINED_CXXFLAGS}"; then
    AC_MSG_WARN([
Currently, g++ is the only C++ compiler for which we know
how to set the compiler switches, but we found ${CXX}.
We will assume CXXFLAGS = \"${CXXFLAGS}\" and go on, but
we strongly recommend you set CXXFLAGS on the command
line.
    ])
  fi
])


dnl     disable generation of exception handling code
AC_DEFUN([LIDIA_DISABLE_EXCEPTIONS],
[
  if test "$GXX" = "yes"; then
    case " ${CXXFLAGS} " in
      *" -fno-exceptions "*)
        ;;
      *)
        CXXFLAGS="${CXXFLAGS} -fno-exceptions"
        ;;
    esac
  fi
])

dnl
dnl	check for LaTeX 2e
dnl
AC_DEFUN([LIDIA_CHECK_PROG_LATEX],
[
	AC_REQUIRE([AM_MISSING_HAS_RUN])
	AC_CHECK_PROG(LATEX, latex, latex,,)
	case $LATEX in
	"")	LATEX="$MISSING latex"
	;;
	*)
		AC_MSG_CHECKING(for LaTeX 2e)
		AC_CACHE_VAL(lidia_cv_latex2e,
		[
			cat > conftest.tex <<EOF
\documentclass{book}
\usepackage{amsmath,amssymb}
\usepackage{graphics}
\usepackage{xspace,url}
\usepackage{geometry}
\usepackage{fontenc}
\usepackage{fancyhdr}
\newcommand{\LiDIA}{\textsf{LiDIA}\xspace}
\begin{document}
  \LiDIA: \url{http://www.informatik.tu-darmstadt.de/TI/LiDIA/}
\end{document}
EOF
			if AC_TRY_COMMAND(${LATEX} conftest.tex) < /dev/null > /dev/null 2>&1; then
				lidia_cv_latex2e="yes"
			else
				lidia_cv_latex2e="no"
			fi
		])
		AC_MSG_RESULT($lidia_cv_latex2e)
	;;
	esac

	if test "$lidia_cv_latex2e" = "yes" ; then
		AC_CHECK_PROG(MAKEINDEX, makeindex, makeindex, :,)
		AC_CHECK_PROG(BIBTEX, bibtex, bibtex, :,)
		AC_CHECK_PROG(DVIPS, dvips, dvips,
			$MISSING dvips,)
		AC_CHECK_PROG(LATEX2HTML, latex2html, latex2html,
			$MISSING latex2html,)
		AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex,
			$MISSING pdflatex,)
		AC_CHECK_PROGS(TEXI2DVI, texi2dvi texi2dvi4a2ps,
			$MISSING texi2dvi,)
	else
		AC_MSG_WARN([
You don't appear to have LaTeX 2e or some required LaTeX 2e packages, therefore you
will not be able to build the manual.  You can obtain the manual in PostScript format
from http://www.informatik.tu-darmstadt.de/TI/LiDIA/])
	fi
])


dnl
dnl	check for GNU MP
dnl
AC_DEFUN([LIDIA_CHECK_LIB_GMP],
[
	AC_MSG_CHECKING(for GNU MP (version >= 3.1))
	AC_CACHE_VAL(lidia_cv_gmp,
	[
		cat > conftest.C <<EOF
#include	<gmp.h>

#if ((__GNU_MP_VERSION < 3) || \
    ((__GNU_MP_VERSION == 3 ) && (__GNU_MP_VERSION_MINOR < 1)))
#error Need GMP >= 3.1!
#endif

int main ()
{
	mpz_t	x, y, z;

	mpz_init(x);
	mpz_init(y);
	mpz_init(z);

	mpz_gcd(z, x, y);

	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(z);
	return 0;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest ${CPPFLAGS} ${EXTRA_INCLUDES} conftest.C ${LIBS} -lgmp) > /dev/null 2>&1; then
			lidia_cv_gmp="yes"
		else
			lidia_cv_gmp="no"
		fi
	])
	AC_MSG_RESULT($lidia_cv_gmp)
])




dnl
dnl	check for CLN
dnl
AC_DEFUN([LIDIA_CHECK_LIB_CLN],
[
	AC_MSG_CHECKING(for CLN)
	AC_CACHE_VAL(lidia_cv_cln,
	[
		cat > conftest.C <<EOF
#include	<cln/integer.h>

int main ()
{
	cln::cl_I	x, y, z;

	z = gcd(x, y);

	return 0;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest ${CPPFLAGS} ${EXTRA_INCLUDES} conftest.C ${LIBS} ${LDFLAGS} -lcln) > /dev/null 2>&1; then
			lidia_cv_cln="yes"
		else
			lidia_cv_cln="no"
		fi
	])
	AC_MSG_RESULT($lidia_cv_cln)
])



dnl
dnl	check for libI
dnl
AC_DEFUN([LIDIA_CHECK_LIB_LIBI],
[
	AC_MSG_CHECKING(for libI)
	AC_CACHE_VAL(lidia_cv_libI,
	[
		cat > conftest.C <<EOF
#include	<iint.h>

int main ()
{
	Integer	x, y, z;

	cI(&x);
	cI(&y);
	cI(&z);

	Idgcd(&z, &x, &y);

	dI(&x);
	dI(&y);
	dI(&z);
	return 0;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest ${CPPFLAGS} ${EXTRA_INCLUDES} conftest.C ${LIBS} -lI) > /dev/null 2>&1; then
			lidia_cv_libI="yes"
		else
			lidia_cv_libI="no"
		fi
	])
	AC_MSG_RESULT($lidia_cv_libI)
])


dnl
dnl	check for Piologie
dnl
AC_DEFUN([LIDIA_CHECK_LIB_PIOLOGIE],
[
	AC_MSG_CHECKING(for Piologie)
	AC_CACHE_VAL(lidia_cv_piologie,
	[
		cat > conftest.C << EOF
#include	<integer.h>

int main ()
{

	Integer x, y, z;

	z = gcd(x, y);

	return 0;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest ${CPPFLAGS} ${EXTRA_INCLUDES} conftest.C ${LIBS} -lpiologie) > /dev/null 2>&1; then
			lidia_cv_piologie="yes"
		else
			lidia_cv_piologie="no"
		fi
	])
	AC_MSG_RESULT($lidia_cv_piologie)
])



dnl
dnl	check for type bool
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_BOOL],
[
	AC_CACHE_CHECK(for type bool,
		lidia_cv_iso_bool,
	[
		cat > conftest.C << EOF 
void frob(bool flag)
{
	if (flag == true) {
		flag = false;
	}
}
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_iso_bool="yes"
		else
			lidia_cv_iso_cast="no"
			AC_MSG_ERROR(Your compiler does not support the ISO C++ bool type.)
		fi
	])
])



dnl
dnl	check whether the keyword inline is supported
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_INLINE],
[
	AC_CACHE_CHECK(for ISO C++ inline,
		lidia_cv_iso_inline,
	[
		cat > conftest.C << EOF
inline long s(long i)
{
	return i * i;
}

static long g(long i, long j)
{
	return s(i) + s(j);
}
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_iso_inline="yes"
		else
			lidia_cv_iso_inline="no"
			AC_MSG_ERROR(Your C++ compiler does not support the ISO C++ inline keyword.)
		fi
	])
])



dnl
dnl	check whether the ISO C++ keyword mutable is supported
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_MUTABLE],
[
	AC_CACHE_CHECK(for mutable class members by ISO C++,
		lidia_cv_iso_mutable,
	[
		cat > conftest.C << EOF
class A {
private:
	mutable long	_m;

public:
	long m() const;
};

long A::m() const
{
	long m = _m;
	_m = 0;
	return m;
}

long f(const A& a)
{
	return a.m();
}
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_iso_mutable="yes"
		else
			lidia_cv_iso_mutable="no"
			AC_MSG_ERROR(Your C++ compiler does not support the ISO C++ mutable keyword.)
		fi
	])
])



dnl
dnl	check whether ISO C++ casting is supported
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_CAST],
[
	AC_CACHE_CHECK(for ISO C++ casting,
		lidia_cv_iso_cast,
	[
		cat > conftest.C << EOF
struct B { virtual void f(); };
struct A : public B { void f(); };

void g (A*);

void f(const A* a, B* b)
{
	const B* c = static_cast<const B*>(a);
	g(const_cast<A*>(a));
}
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_iso_cast="yes"
		else
			lidia_cv_iso_cast="no"
			AC_MSG_ERROR(Your C++ compiler doesn't support ISO C++ casts.)
		fi
	])
])



dnl
dnl	check whether ISO C++ explicit constructors are supported
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_EXPLICIT],
[
	AC_CACHE_CHECK(for explicit constructors by ISO C++,
		lidia_cv_iso_explicit,
	[
		cat > conftest.C << EOF
class A {
public:
	A();
	explicit A(int);
};

void f()
{
	A a(0);
}
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_iso_explicit="yes"
		else
			lidia_cv_iso_explicit="no"
			AC_MSG_ERROR(Your C++ compiler doesn't support explicit constructors by ISO C++.)
		fi
	])
])



dnl
dnl	check whether ISO C++ template specialization works.
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_TEMPLATE_SPEC],
[
	AC_CACHE_CHECK(for working template<>,
		lidia_cv_explicit_templates,
	[
		cat > conftest.C << EOF
template <class T> class c { T t; };
template <> class c<int> { int x; };
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_explicit_templates="yes"
		else
			lidia_cv_explicit_templates="no"
			AC_MSG_ERROR(Your C++ compiler doesn't support ISO C++ template specialization.)
		fi
	])
])




dnl
dnl	check whether ISO C++ explicit template instantiation works.
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_TEMPLATE_INST],
[
	AC_CACHE_CHECK(for explicit template instantiation by ISO C++,
		lidia_cv_iso_templates,
	[
		cat > conftest.C << EOF
template <class T> class X { T t; };

template class X<int>;
EOF
		if AC_TRY_COMMAND(${CXX} -c conftest.C) > /dev/null 2>&1; then
			lidia_cv_iso_templates="yes"
		else
			lidia_cv_iso_templates="no"
			AC_MSG_ERROR(Your C++ compiler doesn't support explicit template instatiation by ISO C++.)
		fi
	])
])



dnl
dnl	check whether ISO C++ namespaces are supported
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_NAMESPACES],
[
	AC_CACHE_CHECK(for ISO C++ namespaces,
		lidia_cv_iso_namespaces,
	[
		cat > conftest.C << EOF
int a();
int a(int);

namespace A {
	int a(int);
}

namespace A {
	int a();
}

int a()
{
	return 0;
}

int a(int i)
{
	return -i;
}

namespace A {
int a(int i)
{
	return i;
}
}

namespace B {
int a(int i)
{
	return 0;
}
}

int A::a()
{
	return 1;
}

int main(int, char**)
{
	if ((a() == 0) && (A::a() == 1) && (a(1) == -1) && (A::a(1) == 1) && (B::a(1) == 0)) {
		return 0;
	}
	return 4;
}
EOF
		lidia_cv_iso_namespaces="no"
		if AC_TRY_COMMAND(${CXX} -o conftest conftest.C) > /dev/null 2>&1 ; then
			if ./conftest ; then
				lidia_cv_iso_namespaces="yes"
			fi
		fi
	])
	if test "$lidia_cv_iso_namespaces" = "yes" ; then
		AC_DEFINE(LIDIA_NAMESPACE, 1)
	fi
])



dnl
dnl	check whether ISO C++ export of templates is supported --- not used, yet.
dnl
AC_DEFUN([LIDIA_CHECK_CXX_ISO_TEMPLATE_EXPORT],
[
	AC_CACHE_CHECK(whether ISO C++ export is supported,
		lidia_cv_iso_export,
	[
		cat > conftest1.C << EOF
export template <class T>T twice(T t) { return t + t; }
EOF
		cat > conftest2.C << EOF
template <class T> T twice(T t);

int f (int i) { return twice(i); }

int main () { return 0; }
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest conftest1.C conftest2.C) > /dev/null 2>&1; then
			lidia_cv_iso_export="yes"
		else
			lidia_cv_iso_export="no"
		fi
	])
])



dnl
dnl	Prepare some object lists
dnl
# LIDIA_SETUP_OBJECTLIST(OBJVAR, LTOBJVAR)
# -- clear the object list OBJVAR and prepare its translation to LTOBJVAR
AC_DEFUN([LIDIA_SETUP_OBJECTLIST],
	[$1=
	AC_SUBST($1)
	AC_CONFIG_COMMANDS_PRE(
		[$2=`echo "$$1 " | sed 's:\.[[^.]]* :.lo :g; s: $::'`
		AC_SUBST($2)])])

AC_DEFUN([LIDIA_SETUP_LIDIA_BASE_XLIBOBJS],
	[LIDIA_SETUP_OBJECTLIST(LIDIA_BASE_XLIBOBJS, LIDIA_BASE_XLTLIBOBJS)])

AC_DEFUN([LIDIA_SETUP_LIDIA_LT_XLIBOBJS],
	[LIDIA_SETUP_OBJECTLIST(LIDIA_LT_XLIBOBJS, LIDIA_LT_XLTLIBOBJS)])



dnl
dnl	Check for applicability of STLport template instantiation requests
dnl
AC_DEFUN([LIDIA_CHECK_STLPORT],
[
	AC_REQUIRE([LIDIA_SETUP_LIDIA_BASE_XLIBOBJS])
	AC_CACHE_CHECK(for use of STLport,
		lidia_cv_use_stlport,
		[AC_EGREP_CPP([_STLP_IOSTREAM],
			[#include <iostream>
			_STLP_IOSTREAM],
			[lidia_cv_use_stlport=no],
			[lidia_cv_use_stlport=yes])])
	case $lidia_cv_use_stlport in yes)
		AC_CACHE_CHECK(for some STLport templates,
			lidia_cv_instantiate_stlport,
			[AC_TRY_COMPILE(
				[#include "${srcdir}/src/portability/stlport.cc"],,
				[lidia_cv_instantiate_stlport=yes],
				[lidia_cv_instantiate_stlport=no])])
		case $lidia_cv_instantiate_stlport in yes)
			LIDIA_BASE_XLIBOBJS="$LIDIA_BASE_XLIBOBJS stlport.$ac_objext"
		;;
		esac
	;;
	esac
])



dnl
dnl	check whether random and srandom are defined in cstdlib.
dnl
AC_DEFUN([LIDIA_CHECK_FUNC_RANDOM],
[
	AC_REQUIRE([LIDIA_SETUP_LIDIA_BASE_XLIBOBJS])
	AC_CACHE_CHECK(for random and srandom,
		lidia_cv_func_random,
	[
		cat > conftest.C << EOF
#include	<cstdlib>

int main (int argc, char** argv) { srandom(random()); return 0; }
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest conftest.C) > /dev/null 2>&1; then
			lidia_cv_func_random="yes"
		else
			lidia_cv_func_random="no"
		fi
	])
	if test "$lidia_cv_func_random" != "yes" ; then
		LIDIA_BASE_XLIBOBJS="$LIDIA_BASE_XLIBOBJS random.$ac_objext"
	else
		AC_DEFINE(HAVE_RANDOM,1)
	fi
])



dnl
dnl	Check for instantiation objects that are used with GCC only
dnl
AC_DEFUN([LIDIA_CHECK_GCC_INSTANTIATION],
[	AC_REQUIRE([LIDIA_SETUP_LIDIA_LT_XLIBOBJS])
	case $GXX in yes)
		LIDIA_LT_XLIBOBJS="$LIDIA_LT_XLIBOBJS lattice_modules_instant.$ac_objext" ;;
	esac
])



dnl
dnl    check whether mkstemp is defined in stdlib.h
dnl
AC_DEFUN([LIDIA_CHECK_FUNC_MKSTEMP],
[
	AC_CACHE_CHECK(for mkstemp,
		lidia_cv_func_mkstemp,
	[
		cat > conftest.C << EOF
#include <stdlib.h>
[
int main (void) {
	char f[ ] = "/tmp/xxXXXXXX";
	return mkstemp(f);
}]
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest conftest.C) > /dev/null 2>&1; then
			lidia_cv_func_mkstemp="yes"
		else
			lidia_cv_func_mkstemp="no"
		fi
	])
	if test "$lidia_cv_func_mkstemp" = "yes" ; then
		AC_DEFINE(HAVE_MKSTEMP,1)
	fi
])



dnl
dnl	check whether POSIX signals work (i.e. sigaction)
dnl
AC_DEFUN([LIDIA_CHECK_SYS_POSIX_SIGNALS],
[
	AC_CACHE_CHECK(whether POSIX signals work,
		lidia_cv_posix_signals,
	[
		lidia_cv_posix_signals="no"
		cat > conftest.C << EOF
#include	<signal.h>
#include	<stdio.h>

static int ret;

void sig_handler (int signo)
{
	if (signo == SIGUSR1) ret = 0;
}

int main (int, char**)
{
	struct sigaction sact;

	ret = 2;
	sact.sa_handler = &sig_handler;
	sact.sa_flags = 0;
	sigemptyset(&sact.sa_mask);
	if (sigaction(SIGUSR1, &sact, NULL) != 0) {
		perror(NULL);
		return 1;
	}
	if (raise(SIGUSR1) != 0) {
		perror(NULL);
		return 1;
	}

	return ret;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest conftest.C) > /dev/null 2>&1; then
			if ./conftest ; then
				lidia_cv_posix_signals="yes"
			fi
		fi
	])
	if test "$lidia_cv_posix_signals" = "yes" ; then
		AC_DEFINE(HAVE_POSIX_SIGNALS, 1)
	fi
])



AC_DEFUN([LIDIA_CHECK_SYS_POSIX_TIMES],
[
	AC_CACHE_CHECK(whether POSIX times work,
		lidia_cv_posix_times,
	[
		lidia_cv_posix_times="no"
		cat > conftest.C << EOF
#include	<sys/times.h>
#include	<stdio.h>

int main (int, char**)
{
	struct tms t_info;
	if (times(&t_info) == static_cast<clock_t>(-1)) {
		perror(NULL);
		return 1;
	}
	return 0;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest conftest.C) > /dev/null 2>&1; then
			if ./conftest ; then
				lidia_cv_posix_times="yes"
			fi
		fi
	])
	if test "$lidia_cv_posix_times" = "yes" ; then
		AC_DEFINE(HAVE_POSIX_TIMES, 1)
	fi
])



AC_DEFUN([LIDIA_CHECK_SYS_POSIX_TIME],
[
	AC_CACHE_CHECK(whether POSIX time work,
		lidia_cv_posix_time,
	[
		lidia_cv_posix_time="no"
		cat > conftest.C << EOF
#include	<time.h>
#include	<stdio.h>

int main (int, char**)
{
	if (time(NULL) == static_cast<time_t>(-1)) {
		perror(NULL);
		return 1;
	}
	return 0;
}
EOF
		if AC_TRY_COMMAND(${CXX} -o conftest conftest.C) > /dev/null 2>&1; then
			if ./conftest ; then
				lidia_cv_posix_time="yes"
			fi
		fi
	])
	if test "$lidia_cv_posix_time" = "yes" ; then
		AC_DEFINE(HAVE_POSIX_TIME, 1)
	fi
])



dnl
dnl    check for gunzip
dnl
AC_DEFUN([LIDIA_CHECK_PROG_GUNZIP],
[
	AC_CHECK_PROG(GUNZIP, gunzip, gunzip)
	if test -n "$GUNZIP" ; then
		AC_DEFINE(HAVE_GUNZIP,1)
	fi
])


dnl
dnl    check for bunzip2
dnl
AC_DEFUN([LIDIA_CHECK_PROG_BUNZIP2],
[
	AC_CHECK_PROG(BUNZIP2, bunzip2, bunzip2)
	if test -n "$BUNZIP2" ; then
		AC_DEFINE(HAVE_BUNZIP2,1)
	fi
])



dnl
dnl LIDIA_PROG_LIBTOOL
dnl	-- a wrapper for AC_PROG_LIBTOOL that works around some bugs
dnl
dnl AC_PROG_LIBTOOL does its actual work in an AC_REQUIRE, therefore the
dnl wrapper must be put into AC_REQUIREs too, in order to be placed right.
dnl
dnl Excuse the non-indentation, it makes the `configure' code look better.
dnl

AC_DEFUN([LIDIA_PROG_LIBTOOL_PRE],
[# Libtool 1.4 uses $CC for linking and manipulates CFLAGS for its tests :-(
lidia_CC=$CC		CC=$CXX
lidia_CFLAGS=$CFLAGS	CFLAGS=$CXXFLAGS
AC_LANG_PUSH(C)])

AC_DEFUN([LIDIA_PROG_LIBTOOL_POST],
[# Restore the variables that were manipulated for the Libtool tests
AC_LANG_POP(C)
CFLAGS=$lidia_CFLAGS
CC=$lidia_CC])

AC_DEFUN([LIDIA_PROG_LIBTOOL],
[dnl
AC_REQUIRE([LIDIA_PROG_LIBTOOL_PRE])dnl
AC_REQUIRE([AC_PROG_LIBTOOL])dnl
AC_REQUIRE([LIDIA_PROG_LIBTOOL_POST])dnl
])


