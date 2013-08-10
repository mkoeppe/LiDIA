##	=========================================================
##
##	Automake definitions for the LiDIA ECO package's examples
##
##	=========================================================

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
	-I$(LIDIA_EC_INCLUDEDIR) -I$(srcdir)/$(LIDIA_EC_INCLUDEDIR) \
	-I$(LIDIA_ECO_INCLUDEDIR) -I$(srcdir)/$(LIDIA_ECO_INCLUDEDIR)


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	eco_2_appl \
	eco_appl \
	find_curve_appl \
	meq_convert_appl \
	verify_order_appl \
	ff1_appl \
	ff2_appl \
	wep_rat_function_appl \
	eco2_appl \
	find_curve_gf2n_appl \
	verify_strong_curve_appl \
	mv_poly_appl \
	mv_term_appl

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

eco_2_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/eco_gfp/eco_2_appl.cc

eco_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/eco_gfp/eco_appl.cc

find_curve_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/eco_gfp/find_curve_appl.cc

meq_convert_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/eco_gfp/meq_convert_appl.cc

verify_order_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/eco_gfp/verify_order_appl.cc

ff1_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/ffx/ff1_appl.cc

ff2_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/ffx/ff2_appl.cc

wep_rat_function_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gfp/wep_rat_function/wep_rat_function_appl.cc

eco2_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gf2n/eco_gf2n/eco2_appl.cc

find_curve_gf2n_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gf2n/eco_gf2n/find_curve_gf2n_appl.cc

verify_strong_curve_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gf2n/eco_gf2n/verify_strong_curve_appl.cc

mv_poly_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gf2n/mv_pol/mv_poly_appl.cc

mv_term_appl_SOURCES = \
	$(LIDIA_ECO_SRCDIR)/eco_gf2n/mv_pol/mv_term_appl.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	${LIDIA_ECO_SRCDIR}/eco_gf2n/mv_pol/mv_poly_timing.cc \
	${LIDIA_ECO_SRCDIR}/eco_gfp/wep_rat_function/wep_rat_function_timing.cc


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

