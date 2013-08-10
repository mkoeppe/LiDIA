##	=======================================================
##
##	Top-level Makefile.am for LiDIA
##
##	=======================================================

## Option `foreign' disables warnings about missing files (AUTHORS, ChangeLog)
## and prevents "automake --force-missing --add-missing" from overwriting
## the `COPYING' file with the GPL.  However, it provides no INSTALL file.

AUTOMAKE_OPTIONS = 1.5 foreign check-news dist-bzip2 dist-zip nostdinc


##	--------------
##	SUBDIRECTORIES
##	--------------

# The complete subdirectory list, used when making `distclean',
# `maintainer-clean', or a standard `dist*' target
DIST_SUBDIRS = library examples doc

# To get finer control over Automake's recursion, we set SUBDIRS=. and do the
# actual recursion by adding dependencies to Automake's `*-recursive' targets.
# Note that `*-local' targets should not be used for extra recursions because
# these targets are subject to subtle internal dependencies which would cause
# multiple recursions or even more trouble.
SUBDIRS = .

.PHONY all-recursive: all-subdirs
all-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) all
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) all
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) all

.PHONY check-recursive: check-subdirs
check-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) check	# updates the lib
	cd examples && $(MAKE) $(AM_MAKEFLAGS) check	# builds examples
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) check

.PHONY dvi-recursive: dvi-subdirs
dvi-subdirs:
##	cd library  && $(MAKE) $(AM_MAKEFLAGS) dvi
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) dvi
	cd doc      && $(MAKE) $(AM_MAKEFLAGS) dvi

.PHONY install-recursive: install-subdirs
install-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) install
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) install
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) install

.PHONY install-data-recursive: install-data-subdirs
install-data-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) install-data
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) install-data
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) install-data

.PHONY install-exec-recursive: install-exec-subdirs
install-exec-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) install-exec
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) install-exec
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) install-exec

# Note that there is no recursive target for `install-strip'
# which uses `install' internally.

.PHONY installdirs-recursive: installdirs-subdirs
installdirs-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) installdirs
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) installdirs
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) installdirs

.PHONY installcheck-recursive: installcheck-subdirs
installcheck-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) installcheck
##	cd examples && $(MAKE) $(AM_MAKEFLAGS) installcheck
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) installcheck

.PHONY uninstall-recursive: uninstall-subdirs
uninstall-subdirs:
	cd doc      && $(MAKE) $(AM_MAKEFLAGS) uninstall	# just in case
	cd examples && $(MAKE) $(AM_MAKEFLAGS) uninstall
	cd library  && $(MAKE) $(AM_MAKEFLAGS) uninstall

.PHONY clean-recursive: clean-subdirs
clean-subdirs:
	cd doc      && $(MAKE) $(AM_MAKEFLAGS) clean
	cd examples && $(MAKE) $(AM_MAKEFLAGS) clean
	cd library  && $(MAKE) $(AM_MAKEFLAGS) clean

.PHONY mostlyclean-recursive: mostlyclean-subdirs
mostlyclean-subdirs:
	cd doc      && $(MAKE) $(AM_MAKEFLAGS) mostlyclean
	cd examples && $(MAKE) $(AM_MAKEFLAGS) mostlyclean
	cd library  && $(MAKE) $(AM_MAKEFLAGS) mostlyclean

.PHONY tags-recursive: tags-subdirs
tags-subdirs:
	cd library  && $(MAKE) $(AM_MAKEFLAGS) tags
	cd examples && $(MAKE) $(AM_MAKEFLAGS) tags
##	cd doc      && $(MAKE) $(AM_MAKEFLAGS) tags

# "make tags" searches $(SUBDIRS) for `TAGS' files.  That's not enough here...
ETAGS_ARGS = -C -l c++ \
	`test -f  library/TAGS && echo \ -i  library/TAGS; \
	 test -f examples/TAGS && echo \ -i examples/TAGS; \
	 test -f      doc/TAGS && echo \ -i      doc/TAGS`

##	----------------
##	ADDITIONAL STUFF
##	----------------

## Compensate for the low strictness level by including some Readmes explicitly
## Added `depcomp' here because Automake 1.5 initially "forgets" it
EXTRA_DIST = \
	LiDIA.spec \
	COPYING NEWS INSTALL README \
	depcomp \
	.LIDIA_SOURCE_DIR \
	Paths.def


##	-------------
##	SPECIAL RULES
##	-------------

# Automake 1.5 does not remove links or stamp files when making distclean
# Autoconf >= 2.54 creates a cache dir which Automake < 1.7 is not aware of
distclean-local:
	rm -rf include src/base/interfaces autom4te.cache

# Convenience targets (some of which are actually no-ops)
.PHONY: \
library			examples		doc			\
all-library		all-examples		all-doc			\
check-library		check-examples		check-doc		\
install-library		install-examples	install-doc		\
install-data-library	install-data-examples	install-data-doc	\
install-exec-library	install-exec-examples	install-exec-doc	\
install-strip-library	install-strip-examples	install-strip-doc	\
installdirs-library	installdirs-examples	installdirs-doc		\
installcheck-library	installcheck-examples	installcheck-doc	\
uninstall-library	uninstall-examples	uninstall-doc		\
clean-library		clean-examples		clean-doc		\
mostlyclean-library	mostlyclean-examples	mostlyclean-doc		\
tags-library		tags-examples		tags-doc		\
pdf dsc ps psgz psbz2

# These rules are for makes that do not respect .PHONY: declarations
library  : all-library  ;
examples : all-examples ;
doc      : all-doc      ;

pdf dsc ps psgz psbz2:
	cd doc && $(MAKE) $(AM_MAKEFLAGS) $@

all-library		all-examples		all-doc			\
check-library		check-examples		check-doc		\
install-library		install-examples	install-doc		\
install-data-library	install-data-examples	install-data-doc	\
install-exec-library	install-exec-examples	install-exec-doc	\
install-strip-library	install-strip-examples	install-strip-doc	\
installdirs-library	installdirs-examples	installdirs-doc		\
installcheck-library	installcheck-examples	installcheck-doc	\
uninstall-library	uninstall-examples	uninstall-doc		\
clean-library		clean-examples		clean-doc		\
mostlyclean-library	mostlyclean-examples	mostlyclean-doc		\
tags-library		tags-examples		tags-doc		\
:
	@target=`expr $@ : '\(.*\)-[^-]*$$'`; \
	subdir=`expr $@ : '.*-\([^-]*\)$$'`; \
	echo "Making $$target in $$subdir"; \
	cd "$$subdir" && $(MAKE) $(AM_MAKEFLAGS) "$$target"

all-examples check-examples install-examples \
install-exec-examples install-strip-examples: library

# Make sure that the configuration headers are up to date
all-subdirs		all-library		all-examples		\
check-subdirs		check-library		check-examples		\
install-subdirs		install-library		install-examples	\
install-data-subdirs	install-data-library	install-data-examples	\
install-exec-subdirs	install-exec-library	install-exec-examples	\
			install-strip-library	install-strip-examples	\
installcheck-subdirs	installcheck-library	installcheck-examples	\
: include/config.h include/LiDIA/LiDIA.h include/LiDIA/path.h


## ------------------------------------------------------------------------
## RULES FOR MAKING SUB-PACKAGE DISTRIBUTIONS
##
## Essentially, these work by calling "make dist" with a modified PKG_DIRS
## variable, thereafter renaming the generated archive.  Doing so for each
## single package and then for all packages together gives the desired
## partial and full distribution archives.  To avoid partial overlaps
## between individual package archives, some `Makefile.am's need a dist-hook
## that removes non-package-specific stuff from any distribution that does
## not include the base package.
##
## Additionally, some checks have been integrated to ensure that the file
## lists of the individual package archives are disjoint and that their
## union equals the file list of the full distribution.  Discrepancies are
## diagnosed and listed in error log files.  The generated list files are
## left around as a possibly interesting by-product.
## ------------------------------------------------------------------------

# When removing $(DISTFILES) in $(distdir), we need to make sure
# that no files outside $(distdir) are deleted.
RELDISTFILES = `echo ' $(DISTFILES)' | sed \
	-e 's, $(srcdir)[\\/], ,g' \
	-e 's, $(top_srcdir)\([\\/]\), $(top_builddir)\1,g'`

# This directory is not specific to an individual sub-package, therefore
# all its stuff goes into the base (or full) distribution and is removed
# from other distributions by this dist-hook.
dist-hook:
	case ' $(PKG_DIRS) ' in \
	  *' base '*) ;; \
	  *) (cd $(distdir) && rm -rf $(RELDISTFILES)) ;; \
	esac
# Additionally, a distribution file list can be produced for consistency
# checks.  However, without further precautions, "make distcheck" would
# complain about the extra list file generated by its internal "make dist".
# To prevent this, one may cause "make distclean" to delete the list file.
# Another approach is to create the list file only on explicit request.
# But how do we safely keep it up-to-date then?  The solution is to update
# the list file iff it already exists.
#
# The "sleep 3" shall ensure that the distribution list gets an earlier
# timestamp than the generated archives, to avoid tie cases (see below).
# The timestamp of files has on some machines a granularity of 2 seconds.
# Therefore, it is not enough to sleep one second. In order to be on the
# safe side, we sleep for three seconds.
	if test -f $(distdir).lst; then \
	  find $(distdir) ! -type d | sort -o $(distdir).lst; \
	  error=$$?; sleep 3; exit $$error; \
	else :; fi			# need "else :" for broken shells

# Filenames for storing lists of duplicated or missing files
dupfiles_lst = dupfiles.lst
fullonly_lst = fullonly.lst
partonly_lst = partonly.lst

# Archives that can be created by "make dist-all"
dist_archives = \
	$(distdir).tar.gz \
	$(distdir).tar.bz2 \
	$(distdir).zip

# What to put into the env.vars GZIP/BZIP2/ZIP when using gzip/bzip2/zip
GZIP_ENV = -9
BZIP_ENV =
BZIP2_ENV = -9
ZIP_ENV = -9

# Note the use of the environment variable AM_MAKEFLAGS to pass changes
# of the PKG_DIRS variable to nested sub-makes.  GNU make would do that
# by itself, but Automake tries to be portable, and so do we.
#
# Automake versions later than 1.5 changed the meaning of `dist' from
# a `dist-gzip' equivalent to `dist-all'.  Therefore we determine the
# set of generated archives dynamically, using the timestamp of
# `$(distdir).lst'.
#
# Solaris sh's test builtin does not suppport -ot or -nt.
# We work around this by using "ls -t" and sed.  The results are
# indeterminate in case of equal timestamps.  That's why the "sleep 3"
# has been added to the dist-hook.

.PHONY: dist-pkgs dist-gzip-pkgs dist-bzip2-pkgs dist-zip-pkgs dist-all-pkgs

dist-pkgs dist-gzip-pkgs dist-bzip2-pkgs dist-zip-pkgs dist-all-pkgs:
	@GZIP='$(GZIP_ENV)'; \
	BZIP='$(BZIP_ENV)'; \
	BZIP2='$(BZIP2_ENV)'; \
	ZIP='$(ZIP_ENV)'; \
	export GZIP BZIP BZIP2 ZIP; \
	dist=`expr $@ : '\(.*\)-pkgs$$'`; \
	trap \
	  'if test -f $(distdir).lst; then \
	    rm -f `ls -t $(distdir).lst $(dist_archives) 2>/dev/null | \
	      sed "/\\.lst\$$/q"`; \
	  fi' 0; \
	partlists=; \
	for p in $(PKG_DIRS); do \
	  case $$p in \
	    finite_fields) part=FF ;; \
	    linear_algebra) part=LA ;; \
	    lattices) part=LT ;; \
	    number_fields) part=NF ;; \
	    elliptic_curves) part=EC ;; \
	    ec_order) part=ECO ;; \
	    gec) part=GEC ;; \
	    *) part=$$p ;; \
	  esac; \
	  archive=$(PACKAGE)-$$part-$(VERSION); \
	  echo "Making $$dist of $$archive"; \
	  echo "" >$(distdir).lst; \
	  AM_MAKEFLAGS="$$AM_MAKEFLAGS PKG_DIRS='$$p'" \
	  $(MAKE) $(AM_MAKEFLAGS) PKG_DIRS="$$p" $$dist || exit $$?; \
	  ls -t $(distdir).lst $(dist_archives) 2>/dev/null | sed -n \
	  -e "s:^$(distdir)\\.\\(.*\\)\$$:echo mv -f & $$archive.\\1:p" \
	  -e 's:^echo ::p' \
	  -e '/\.lst$$/q' | $(SHELL) -e || exit $$?; \
	  partlists="$$partlists $$archive.lst"; \
	done; \
	archive=$(PACKAGE)-$(VERSION); \
	echo "Making full $$dist of $$archive"; \
	echo "" >$(distdir).lst; \
	$(MAKE) $(AM_MAKEFLAGS) $$dist || exit $$?; \
	trap 0; \
	case $$partlists in "") exit 0 ;; esac; \
	error=0; \
	echo "Checking disjointness of the partial $$dist archives:"; \
	if sort -m $$partlists | uniq -d | \
	   tee $(dupfiles_lst) | grep .; then \
	  echo "Error: The above files occur in several $$dist parts."; \
	  echo "Duplicate files list stored in \`$(dupfiles_lst)'."; \
	  error=1; \
	else \
	  echo "OK: No duplicate files among different $$dist parts."; \
	  rm -f $(dupfiles_lst); \
	fi; \
	echo "Checking completeness of the partial $$dist archives:"; \
	if sort -m -u $$partlists | comm -23 $(distdir).lst - | \
	   tee $(fullonly_lst) | grep .; then \
	  echo "Error: The above files are only in the full $$dist archive(s)."; \
	  echo "Extra files list stored in \`$(fullonly_lst)'."; \
	  error=1; \
	else \
	  echo "OK: In total, the parts contain the full $$dist stuff."; \
	  rm -f $(fullonly_lst); \
	fi; \
	echo "Checking completeness of the full $$dist archive(s):"; \
	if sort -m -u $$partlists | comm -23 - $(distdir).lst | \
	   tee $(partonly_lst) | grep .; then \
	  echo "Error: The above files are only in some partial $$dist archive."; \
	  echo "Missing files list stored in \`$(partonly_lst)'."; \
	  error=1; \
	else \
	  echo "OK: The full $$dist archive(s) contain all parts."; \
	  rm -f $(partonly_lst); \
	fi; \
	echo "Check summary:"; \
	if test $$error -ne 0; then \
	  echo "There are missing or duplicated files in the generated $$dist archives."; \
	  echo "Check the produced file lists and review the \`Makefile.am's."; \
	else \
	  echo "OK: The full $$dist archive(s) are the disjoint union of their parts."; \
	fi; \
	exit $$error

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

