##	=========================================================
##
##	Automake definitions for the LiDIA gec package's examples
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
	-I$(LIDIA_NF_INCLUDEDIR) -I$(srcdir)/$(LIDIA_NF_INCLUDEDIR) \
	-I$(LIDIA_EC_INCLUDEDIR) -I$(srcdir)/$(LIDIA_EC_INCLUDEDIR) \
	-I$(LIDIA_ECO_INCLUDEDIR) -I$(srcdir)/$(LIDIA_ECO_INCLUDEDIR) \
	-I$(LIDIA_GEC_INCLUDEDIR) -I$(srcdir)/$(LIDIA_GEC_INCLUDEDIR)

##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	gec_appl \
	prime_proof_appl

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

gec_appl_SOURCES = \
	$(LIDIA_GEC_SRCDIR)/gec_appl.cc 


prime_proof_appl_SOURCES = \
	$(LIDIA_GEC_SRCDIR)/prime_proof/prime_proof_appl.cc


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

