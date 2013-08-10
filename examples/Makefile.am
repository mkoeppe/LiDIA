##	=======================================================
##
##	Automake definitions for the LiDIA application examples
##
##	=======================================================

AUTOMAKE_OPTIONS = 1.5 nostdinc

ETAGS_ARGS = -C -l c++

## include $(top_srcdir)/Paths.def


##	--------------
##	SUBDIRECTORIES
##	--------------

SUBDIRS = $(PKG_DIRS)
DIST_SUBDIRS = $(PKG_DIRS)


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST =


##	-------------
##	SPECIAL RULES
##	-------------

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

