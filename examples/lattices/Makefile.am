##	========================================================
##
##	Automake definitions for the LiDIA LT package's examples
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
	-I$(LIDIA_LT_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LT_INCLUDEDIR)


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	lattice_appl

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

lattice_appl_SOURCES = \
	$(LIDIA_LT_SRCDIR)/lattice_appl.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	${LIDIA_LT_SRCDIR}/Lattices/d10 \
	${LIDIA_LT_SRCDIR}/Lattices/d10.gp \
	${LIDIA_LT_SRCDIR}/Lattices/d15 \
	${LIDIA_LT_SRCDIR}/Lattices/d15.gp \
	${LIDIA_LT_SRCDIR}/Lattices/d20 \
	${LIDIA_LT_SRCDIR}/Lattices/d20.gp \
	${LIDIA_LT_SRCDIR}/Lattices/d25 \
	${LIDIA_LT_SRCDIR}/Lattices/d25.gp \
	${LIDIA_LT_SRCDIR}/Lattices/d30 \
	${LIDIA_LT_SRCDIR}/Lattices/d30.gp \
	${LIDIA_LT_SRCDIR}/Lattices/d5 \
	${LIDIA_LT_SRCDIR}/Lattices/d5.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t10 \
	${LIDIA_LT_SRCDIR}/Lattices/t10.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t15 \
	${LIDIA_LT_SRCDIR}/Lattices/t15.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t20 \
	${LIDIA_LT_SRCDIR}/Lattices/t20.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t25 \
	${LIDIA_LT_SRCDIR}/Lattices/t25.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t30 \
	${LIDIA_LT_SRCDIR}/Lattices/t30.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t40 \
	${LIDIA_LT_SRCDIR}/Lattices/t40.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t5 \
	${LIDIA_LT_SRCDIR}/Lattices/t5.gp \
	${LIDIA_LT_SRCDIR}/Lattices/t50 \
	${LIDIA_LT_SRCDIR}/Lattices/t50.gp


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

