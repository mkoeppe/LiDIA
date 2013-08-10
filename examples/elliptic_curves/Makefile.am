##	========================================================
##
##	Automake definitions for the LiDIA EC package's examples
##
##	========================================================

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


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	complex_periods_appl \
	curve_isomorphism_appl \
	ec_bigint_appl \
	ec_bigint_test_appl \
	ec_bigrational_appl \
	elliptic_curve_bigint_appl \
	minimal_model_appl \
	point_bigint_appl \
	quartic_appl \
	ec_gf_element_appl \
	ec_point_gf_element_appl

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

complex_periods_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/complex_periods_appl.cc

curve_isomorphism_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/curve_isomorphism_appl.cc

ec_bigint_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/ec_bigint_appl.cc

ec_bigint_test_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/ec_bigint_test_appl.cc

ec_bigrational_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/ec_bigrational_appl.cc

elliptic_curve_bigint_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/elliptic_curve_bigint_appl.cc

minimal_model_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/minimal_model_appl.cc

point_bigint_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/point_bigint_appl.cc

quartic_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_rationals/quartic_appl.cc

ec_gf_element_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_finite_fields/ec_gf_element_appl.cc

ec_point_gf_element_appl_SOURCES = \
	$(LIDIA_EC_SRCDIR)/ec_finite_fields/ec_point_gf_element_appl.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/README \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out1 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out2 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out3 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out4 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out5 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out6 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/complex_periods_appl.out7 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/curve_isomorphism_appl.out \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data1 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data2 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data3 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data4 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data5 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data6 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/data7 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out1 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out2 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out3 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out4 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out5 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out6 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_appl.out7 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/elliptic_curve_bigint_test.in \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out1 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out2 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out3 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out4 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out5 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out6 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/minimal_model_appl.out7 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out1 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out2 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out3 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out4 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out5 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out6 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/point_bigint_appl.out7 \
	${LIDIA_EC_SRCDIR}/ec_rationals/appl_data/quartic_appl.out \
	${LIDIA_EC_SRCDIR}/ec_rationals/testit.sh


##	-------------
##	SPECIAL RULES
##	-------------

# Remove installation subdirectories after uninstalling, unless left non-empty.
uninstall-hook:
	rmdir	$(DESTDIR)$(applexecdir) \
		$(DESTDIR)$(pkglibdir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

