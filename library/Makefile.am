##	================================================
##
##	Automake definitions for the LiDIA library stuff
##
##	================================================

AUTOMAKE_OPTIONS = 1.5 nostdinc

ETAGS_ARGS = -C -l c++

include $(top_srcdir)/Paths.def


##	--------------
##	SUBDIRECTORIES
##	--------------

SUBDIRS = $(PKG_DIRS)
DIST_SUBDIRS = $(PKG_DIRS)


##	-------------
##	WHAT TO BUILD
##	-------------

## Specifications for `libLiDIA.la'.
## All source files are compiled in the $(SUBDIRS), so we just link
## the resulting libtool convenience archives $(LIDIA_LTCLIBS) here.
## Do not add $(LIBS) to libLiDIA_la_LIBADD, they are linked in anyway.

lib_LTLIBRARIES = libLiDIA.la
libLiDIA_la_LIBADD = $(LIDIA_LTCLIBS)
libLiDIA_la_DEPENDENCIES = $(LIDIA_LTCLIBS)
libLiDIA_la_LDFLAGS = -version-info $(LIDIA_LTVERSION)
libLiDIA_la_SOURCES =

## Tell Automake to use the C++ linker.
## This is necessary because no source files are specified here.

CXXLD = $(CXX)
libLiDIA_la_LINK = \
	$(LIBTOOL) --mode=link $(CXXLD) $(AM_LDFLAGS) $(LDFLAGS) -o $@


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST =


##	-------------
##	SPECIAL RULES
##	-------------

# # Some supplemental rules to ensure that "make libLiDIA.la" updates
# # the sub-libraries like "make [all]" does.  Unfortunately, they cause
# # "make [all]" to recurse into $(SUBDIRS) twice, with the second pass
# # re-checking all dependencies just for finding "nothing to be done" :-(
#
# $(LIDIA_LTCLIBS): .FORCE
#	cd $(@D) && $(MAKE) $(AM_MAKEFLAGS) $(@F)
#
# .PHONY: .FORCE
# .FORCE:

# When removing $(DISTFILES) in $(distdir), we need to make sure
# that no files outside $(distdir) are deleted.
RELDISTFILES = `echo ' $(DISTFILES)' | sed \
	-e 's, $(srcdir)[\\/], ,g' \
	-e 's, $(top_srcdir)\([\\/]\), $(top_builddir)\1,g'`

# This directory is not specific to a certain sub-package, therefore all its
# stuff goes into the base (or full) distribution and is removed from other
# distributions.
dist-hook:
	case ' $(PKG_DIRS) ' in \
	  *' base '*) ;; \
	  *) (cd $(distdir) && rm -rf $(RELDISTFILES)) ;; \
	esac

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

