##	=============================================================
##
##	Automake definitions for the LiDIA EC package's library stuff
##
##	=============================================================

AUTOMAKE_OPTIONS = 1.5 nostdinc

ETAGS_ARGS = -C -l c++

include $(top_srcdir)/Paths.def


##	-------------
##	INCLUDE PATHS
##	-------------

AM_CPPFLAGS = \
	-I$(LIDIA_CONFIG_INCLUDEDIR) \
	-I$(LIDIA_BASE_INCLUDEDIR) -I$(srcdir)/$(LIDIA_BASE_INCLUDEDIR) \
	-I$(LIDIA_FF_INCLUDEDIR) -I$(srcdir)/$(LIDIA_FF_INCLUDEDIR) \
	-I$(LIDIA_LA_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LA_INCLUDEDIR) \
	-I$(LIDIA_EC_INCLUDEDIR) -I$(srcdir)/$(LIDIA_EC_INCLUDEDIR)


##	---------
##	LIBRARIES
##	---------

noinst_LTLIBRARIES = libLiDIA-EC.la

# We don't want to have $(LIBS) in the convenience library.
# FIXME: This is not safe against command-line overrides.
LIBS =


##	---------------
##	LIBRARY SOURCES
##	---------------

libLiDIA_EC_la_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_finite_fields/elliptic_curve_gf_element.cc \
	$(LIDIA_EC_SRCDIR)/ec_finite_fields/elliptic_curve_rep_gf_element.cc \
	$(LIDIA_EC_SRCDIR)/ec_finite_fields/point_gf_element.cc \
	$(LIDIA_EC_SRCDIR)/ec_finite_fields/point_operations_gf_element.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/complex_periods.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/ec_arith.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/elliptic_curve_bigint.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/elliptic_curve_rep_bigint.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/instant_curve_isomorphism.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/instant_ec_bigint.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/kodaira_code.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/minimal_model.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/point_bigint.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/point_operations_bigint.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/quartic.cc \
	$(LIDIA_EC_SRCDIR)/ec_rationals/reduction_type.cc \
	$(LIDIA_EC_SRCDIR)/instantiate/bv_point_bigint.cc \
	$(LIDIA_EC_SRCDIR)/instantiate/ec_bigrational.cc \
	$(LIDIA_EC_SRCDIR)/instantiate/ec_gf_element.cc


##	-------
##	HEADERS
##	-------

instantincludedir = $(pkgincludedir)/instantiate
ecincludedir = $(pkgincludedir)/elliptic_curves

instantinclude_HEADERS = \
	$(LIDIA_EC_PKGINCLUDEDIR)/instantiate/ec_finite_fields.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/instantiate/ec_rationals.cc

ecinclude_HEADERS = \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/base_elliptic_curve.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/base_elliptic_curve.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/base_elliptic_curve_rep.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/base_elliptic_curve_rep.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/base_point.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/base_point.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/ec_arith.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/elliptic_curve_rep.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/elliptic_curve_rep.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/elliptic_curve_rep_bigint.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/point_operations.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/point_operations.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curves/point_operations_bigint.h

pkginclude_HEADERS = \
	$(LIDIA_EC_PKGINCLUDEDIR)/complex_periods.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/curve_isomorphism.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/curve_isomorphism.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curve.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curve.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curve_bigint.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/elliptic_curve_flags.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/kodaira_code.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/minimal_model.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/point.cc \
	$(LIDIA_EC_PKGINCLUDEDIR)/point.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/point_bigint.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/quartic.h \
	$(LIDIA_EC_PKGINCLUDEDIR)/reduction_type.h


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST =


##	-------------
##	SPECIAL RULES
##	-------------

# Remove installation subdirectories after uninstalling, unless left non-empty.
uninstall-hook:
	rmdir	$(DESTDIR)$(instantincludedir) \
		$(DESTDIR)$(ecincludedir) \
		$(DESTDIR)$(pkgincludedir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

