##	========================================================
##
##	Automake definitions for the LiDIA NF package's examples
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
	-I$(LIDIA_LT_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LT_INCLUDEDIR) \
	-I$(LIDIA_NF_INCLUDEDIR) -I$(srcdir)/$(LIDIA_NF_INCLUDEDIR)


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	alg_number_appl \
	ideal_appl \
	ideal_factorization_appl \
	module_appl \
	findh_appl \
	qi_class_appl \
	qi_class_real_appl \
	quadratic_form_appl \
	quadratic_ideal_appl \
	quadratic_ideal_appl2 \
	quadratic_ideal_appl3 \
	quadratic_number_power_product_appl \
	quadratic_number_standard_appl \
	quadratic_order_appl \
	quadratic_order_appl2

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

alg_number_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/algebraic_num/alg_number_appl.cc

ideal_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/algebraic_num/ideal_appl.cc

ideal_factorization_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/algebraic_num/ideal_factorization_appl.cc

module_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/algebraic_num/module_appl.cc

findh_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/findh_appl.cc

qi_class_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qi_class_appl.cc

qi_class_real_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qi_class_real_appl.cc

quadratic_form_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_form_appl.cc

quadratic_ideal_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_ideal_appl.cc

quadratic_ideal_appl2_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_ideal_appl2.cc

quadratic_ideal_appl3_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_ideal_appl3.cc

quadratic_number_power_product_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_number_power_product_appl.cc

quadratic_number_standard_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_number_standard_appl.cc

quadratic_order_appl_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_order_appl.cc

quadratic_order_appl2_SOURCES = \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_order_appl2.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST =


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

