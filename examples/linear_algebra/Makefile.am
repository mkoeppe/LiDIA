##	========================================================
##
##	Automake definitions for the LiDIA LA package's examples
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
	-I$(LIDIA_LA_INCLUDEDIR) -I$(srcdir)/$(LIDIA_LA_INCLUDEDIR)


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	bigint_matrix_appl \
	field_matrix_appl \
	ring_matrix_appl 

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

bigint_matrix_appl_SOURCES = \
	$(LIDIA_LA_SRCDIR)/bigint_matrix_appl.cc

field_matrix_appl_SOURCES = \
	$(LIDIA_LA_SRCDIR)/field_matrix_appl.cc

ring_matrix_appl_SOURCES = \
	$(LIDIA_LA_SRCDIR)/ring_matrix_appl.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	${LIDIA_LA_SRCDIR}/bigint_matrix_appl.dat \
	${LIDIA_LA_SRCDIR}/ring_matrix_appl_bigint.dat \
	${LIDIA_LA_SRCDIR}/ring_matrix_appl_int.dat \
	${LIDIA_LA_SRCDIR}/ring_matrix_appl_long.dat 


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

