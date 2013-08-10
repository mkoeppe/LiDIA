##	=============================================================
##
##	Automake definitions for the LiDIA LA package's library stuff
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
	-I$(LIDIA_LA_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LA_INCLUDEDIR)


##	---------
##	LIBRARIES
##	---------

noinst_LTLIBRARIES = libLiDIA-LA.la

# We don't want to have $(LIBS) in the convenience library.
# FIXME: This is not safe against command-line overrides.
LIBS =


##	---------------
##	LIBRARY SOURCES
##	---------------

libLiDIA_LA_la_SOURCES = \
	$(LIDIA_LA_SRCDIR)/bigfloat_matrix.cc \
	$(LIDIA_LA_SRCDIR)/bigint_matrix.cc \
	$(LIDIA_LA_SRCDIR)/bigmod_matrix.cc \
	$(LIDIA_LA_SRCDIR)/crt_and_prime_handling.cc \
	$(LIDIA_LA_SRCDIR)/file_adjoint.cc \
	$(LIDIA_LA_SRCDIR)/hnf_instant.cc \
	$(LIDIA_LA_SRCDIR)/instant_bigint_matrix.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigfloat_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigfloat_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigfloat_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigfloat_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigrational_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigrational_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigrational_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bm_bigrational_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bv_bigfloat_int.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bv_bigint_matrix.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bv_bpp_matrix_bigint.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/bv_bpp_matrix_long.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigfloat_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigfloat_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigfloat_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigfloat_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigrational_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigrational_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigrational_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_bigrational_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_double_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_double_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_double_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/fm_double_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigfloat_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigfloat_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigfloat_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigfloat_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigint_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigint_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigint_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigint_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigrational_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigrational_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigrational_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_bigrational_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_double_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_double_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_double_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_double_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_int_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_int_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_int_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_int_sub.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_long_dense.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_long_normal.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_long_sparse.cc \
	$(LIDIA_LA_SRCDIR)/instantiate/rm_long_sub.cc \
	$(LIDIA_LA_SRCDIR)/m_multiple_gcd.cc \
	$(LIDIA_LA_SRCDIR)/modular_instant.cc \
	$(LIDIA_LA_SRCDIR)/trans_matrix.cc


##	-------
##	HEADERS
##	-------

matrixincludedir = $(pkgincludedir)/matrix

matrixinclude_HEADERS = \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/bigint_matrix_algorithms.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/bigint_matrix_algorithms.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/crt_and_prime_handling.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/dense_bigint_matrix_kernel.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/dense_bigint_matrix_kernel.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/dense_bigint_matrix_modules.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/dense_bigint_matrix_modules.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/dense_fp_matrix_kernel.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/dense_fp_matrix_kernel.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/fp_matrix_algorithms.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/fp_matrix_algorithms.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/hnf_conf.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/hnf_conf.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/hnf_kernel.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/hnf_kernel.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/modular_arithmetic.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/modular_arithmetic.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/modular_functions.inl \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/normalize_kernel.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/normalize_kernel.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_bigint_matrix_kernel.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_bigint_matrix_kernel.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_bigint_matrix_modules.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_bigint_matrix_modules.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_fp_matrix_algorithms.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_fp_matrix_algorithms.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_fp_matrix_kernel.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix/sparse_fp_matrix_kernel.h

pkginclude_HEADERS = \
	$(LIDIA_LA_PKGINCLUDEDIR)/base_power_product.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/base_power_product.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/base_ppair.cc \
	$(LIDIA_LA_PKGINCLUDEDIR)/base_ppair.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/bigfloat_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/bigint_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/bigmod_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/dense_bigint_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/dense_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/long_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/math_matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/matrix_intern.h \
	$(LIDIA_LA_PKGINCLUDEDIR)/sparse_matrix.h


##	----
##	DATA
##	----

dist_pkgdata_DATA = \
	$(LIDIA_PKGDATADIR)/LIDIA_PRIMES


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST =


##	-------------
##	SPECIAL RULES
##	-------------

# Remove installation subdirectories after uninstalling, unless left non-empty.
uninstall-hook:
	rmdir	$(DESTDIR)$(matrixincludedir) \
		$(DESTDIR)$(pkgincludedir) \
		$(DESTDIR)$(pkgdatadir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

