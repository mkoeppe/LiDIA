##	=============================================================
##
##	Automake definitions for the LiDIA LT package's library stuff
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
	-I$(LIDIA_LT_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LT_INCLUDEDIR)


##	---------
##	LIBRARIES
##	---------

noinst_LTLIBRARIES = libLiDIA-LT.la

# Add the configured supplementary object files
libLiDIA_LT_la_LIBADD = @LIDIA_LT_XLTLIBOBJS@
libLiDIA_LT_la_DEPENDENCIES = @LIDIA_LT_XLTLIBOBJS@

# We don't want to have $(LIBS) in the convenience library.
# FIXME: This is not safe against command-line overrides.
LIBS =


##	---------------
##	LIBRARY SOURCES
##	---------------

libLiDIA_LT_la_SOURCES = \
	$(LIDIA_LT_SRCDIR)/bf_basis.cc \
	$(LIDIA_LT_SRCDIR)/bf_gsys_alg1.cc \
	$(LIDIA_LT_SRCDIR)/bf_gsys_alg2.cc \
	$(LIDIA_LT_SRCDIR)/bf_gsys_base.cc \
	$(LIDIA_LT_SRCDIR)/bfl_alg1.cc \
	$(LIDIA_LT_SRCDIR)/bfl_base.cc \
	$(LIDIA_LT_SRCDIR)/bi_basis.cc \
	$(LIDIA_LT_SRCDIR)/bi_gsys_alg1.cc \
	$(LIDIA_LT_SRCDIR)/bi_gsys_alg2.cc \
	$(LIDIA_LT_SRCDIR)/bi_gsys_alg3.cc \
	$(LIDIA_LT_SRCDIR)/bi_gsys_base.cc \
	$(LIDIA_LT_SRCDIR)/bil_alg1.cc \
	$(LIDIA_LT_SRCDIR)/bil_alg2.cc \
	$(LIDIA_LT_SRCDIR)/bil_alg3.cc \
	$(LIDIA_LT_SRCDIR)/bil_base.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_bfl1.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_bfl2.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_bfl3.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_bfl4.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_dbl1.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_dbl2.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_dbl3.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_dbl4.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_xdbl1.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_xdbl2.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_xdbl3.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bfl_xdbl4.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl1.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl2.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl3.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl4.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl5.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl6.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl7.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_bfl8.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl1.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl2.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl3.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl4.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl5.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl6.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl7.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_dbl8.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl1.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl2.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl3.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl4.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl5.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl6.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl7.cc \
	$(LIDIA_LT_SRCDIR)/instantiate/bil_xdbl8.cc \
	$(LIDIA_LT_SRCDIR)/lattice_basis.cc \
	$(LIDIA_LT_SRCDIR)/lattice_gensys.cc

EXTRA_libLiDIA_LT_la_SOURCES = \
	$(LIDIA_LT_SRCDIR)/instantiate/lattice_modules_instant.cc


##	-------
##	HEADERS
##	-------

instantincludedir = $(pkgincludedir)/instantiate
ltincludedir = $(pkgincludedir)/lattices

instantinclude_HEADERS = \
	$(LIDIA_LT_PKGINCLUDEDIR)/instantiate/lattices.cc

ltinclude_HEADERS = \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/bf_lattice_basis.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/bf_lattice_gensys.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/bi_lattice_basis.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/bi_lattice_gensys.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/lattice_defs.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/lattice_kernel.cc \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/lattice_kernel.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/lattice_modules.cc \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/lattice_modules.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/p_vector.cc \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattices/p_vector.h

pkginclude_HEADERS = \
	$(LIDIA_LT_PKGINCLUDEDIR)/bigfloat_lattice.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/bigint_lattice.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattice_basis.h \
	$(LIDIA_LT_PKGINCLUDEDIR)/lattice_gensys.h


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
		$(DESTDIR)$(ltincludedir) \
		$(DESTDIR)$(pkgincludedir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

