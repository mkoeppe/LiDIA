##	=============================================================
##
##	Automake definitions for the LiDIA NF package's library stuff
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
	-I$(LIDIA_LT_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LT_INCLUDEDIR) \
	-I$(LIDIA_NF_INCLUDEDIR) -I$(srcdir)/$(LIDIA_NF_INCLUDEDIR)


##	---------
##	LIBRARIES
##	---------

noinst_LTLIBRARIES = libLiDIA-NF.la

# We don't want to have $(LIBS) in the convenience library.
# FIXME: This is not safe against command-line overrides.
LIBS =


##	---------------
##	LIBRARY SOURCES
##	---------------

libLiDIA_NF_la_SOURCES = \
	$(LIDIA_NF_SRCDIR)/algebraic_num/alg_factor.cc \
	$(LIDIA_NF_SRCDIR)/algebraic_num/alg_number.cc \
	$(LIDIA_NF_SRCDIR)/algebraic_num/factorization/sf_alg_ideal.cc \
	$(LIDIA_NF_SRCDIR)/algebraic_num/ideal.cc \
	$(LIDIA_NF_SRCDIR)/algebraic_num/module.cc \
	$(LIDIA_NF_SRCDIR)/algebraic_num/order.cc \
	$(LIDIA_NF_SRCDIR)/algebraic_num/prime_ideal.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bpf_quadratic_number_standard.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bpp_quadratic_ideal.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_bpp_qi_bigint.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_p_bigint.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_partial_relation.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_pp_sf_alg_ideal.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_qi_class.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_qi_class_real.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_quadratic_form.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_quadratic_ideal.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_quadratic_number_power_product.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_quadratic_number_standard.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/bv_xbigfloat.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/fact_alg_ideal.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/ht_ideal_node.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/ht_qi_class_real.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/ht_quadratic_form.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/iht_ideal_node.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/iht_qi_class_real.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/lpf_quadratic_number_standard.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/mv_quadratic_number_standard.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/pair_bigint.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/rcl_quadratic_number_power_product_basis.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/sv_p_bigint.cc \
	$(LIDIA_NF_SRCDIR)/instantiate/sv_pp_sf_alg_ideal.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/compute_class_group.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/matrix/bigint_matrix_qo.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qi_class.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qi_class_real.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qo_list.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qo_sieve.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qo_subexp.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/qo_util.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_form.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_ideal.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_ideal_power_product.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_number_logarithm.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_number_power_product.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_number_power_product_basis.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_number_standard.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_order1.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_order2.cc \
	$(LIDIA_NF_SRCDIR)/quadratic_order/quadratic_order3.cc


##	-------
##	HEADERS
##	-------

instantincludedir = $(pkgincludedir)/instantiate
specialincludedir = $(pkgincludedir)/specialization
nfincludedir = $(pkgincludedir)/number_fields

instantinclude_HEADERS = \
	$(LIDIA_NF_PKGINCLUDEDIR)/instantiate/power_product.cc

specialinclude_HEADERS = \
	$(LIDIA_NF_PKGINCLUDEDIR)/specialization/alg_ideal.special

nfinclude_HEADERS = \
	$(LIDIA_NF_PKGINCLUDEDIR)/number_fields/partial_relation.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/number_fields/qo_list.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/number_fields/qo_sieve.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/number_fields/qo_util.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/number_fields/sf_alg_ideal.h

pkginclude_HEADERS = \
	$(LIDIA_NF_PKGINCLUDEDIR)/alg_number.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/prime_ideal.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/qi_class.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/qi_class_real.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_form.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_ideal.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_ideal_power_product.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_number_logarithm.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_number_power_product.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_number_power_product_basis.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_number_standard.h \
	$(LIDIA_NF_PKGINCLUDEDIR)/quadratic_order.h


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
		$(DESTDIR)$(specialincludedir) \
		$(DESTDIR)$(nfincludedir) \
		$(DESTDIR)$(pkgincludedir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

