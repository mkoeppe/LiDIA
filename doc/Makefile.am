##	================================================
##
##	Automake definitions for the LiDIA documentation
##
##	================================================

AUTOMAKE_OPTIONS = 1.5 nostdinc

## Because Automake cannot automatically handle LaTeX sources yet,
## we have to provide the compilation rules explicitly.
## Nevertheless, Automake generates the usual maintainer rules, like
## those used by "make dist".

MOSTLYCLEANFILES = \
	*.aux *.toc *.lof *.lot *.bbl \
	*.log *.[bdi]lg *.idx *.ind *.gl[os]

CLEANFILES = \
	LiDIA.dvi	LiDIA.ps	LiDIA.pdf	LiDIA.dsc \
	LiDIA.dvi.gz	LiDIA.ps.gz	LiDIA.pdf.gz	LiDIA.dsc.gz \
	LiDIA.dvi.bz2	LiDIA.ps.bz2	LiDIA.pdf.bz2	LiDIA.dsc.bz2


##	-------------
##	WHAT TO BUILD
##	-------------

noinst_DATA = LiDIA.dvi

LiDIA_dvi_TEXSRCS = \
	LiDIA.tex \
		graphics.cfg \
		color.cfg \
		LiDIA.def \
		copyright.tex \
		overview.tex \
			lidia_model.tex \
		use.tex \
		LiDIA_base.tex \
			bigint.tex \
			udigit_mod.tex \
			bigmod.tex \
			multi_bigmod.tex \
			bigrational.tex \
			xdouble.tex \
			bigfloat.tex \
			xbigfloat.tex \
			bigcomplex.tex \
			base_vector.tex \
			math_vector.tex \
			sort_vector.tex \
			file_vector.tex \
			lidia_vector.tex \
			hash_table.tex \
			indexed_hash_table.tex \
			base_matrix.tex \
			math_matrix.tex \
			matrix_GL2Z.tex \
			polynomial.tex \
			big_polynomial.tex \
			fact.tex \
			single_factor.tex \
			factorization.tex \
			sf_bigint.tex \
			modular_functions.tex \
			nmbrthry_functions.tex \
			crt.tex \
			prime_list.tex \
			prime_proof.tex \
			power_functions.tex \
			timer.tex \
			handlers.tex \
			basic_error.tex \
			generic_error.tex \
			cast_error.tex \
			precondition_error.tex \
			index_out_of_bounds_error.tex \
			gmm.tex \
			random_generator.tex \
			lidia_signal.tex \
		LiDIA_FF.tex \
			gf2n.tex \
			galois_field.tex \
			galois_field_iterator.tex \
			gf_element.tex \
			gf_polynomial.tex \
			Fp_polynomial.tex \
			Fp_poly_modulus.tex \
			gf_poly_modulus.tex \
			Fp_rational_function.tex \
			power_series.tex \
			sf_Fp_polynomial.tex \
			sf_gf_polynomial.tex \
			fft_prime.tex \
		LiDIA_LA.tex \
			bigint_matrix.tex \
			bigmod_matrix.tex \
		LiDIA_LT.tex \
			bigint_lattice.tex \
			bigfloat_lattice.tex \
		LiDIA_NF.tex \
			quadratic_form.tex \
			quadratic_order.tex \
			qi_class.tex \
			qi_class_real.tex \
			quadratic_ideal.tex \
			quadratic_number_standard.tex \
			quadratic_number_power_product_basis.tex \
			quadratic_number_power_product.tex \
			field.tex \
			order.tex \
			alg_number.tex \
			module.tex \
			prime_ideal.tex \
			sf_alg_ideal.tex \
		LiDIA_EC.tex \
			elliptic_curve_flags.tex \
			elliptic_curve.tex \
			point.tex \
			elliptic_curve_bigint.tex \
			point_bigint.tex \
			curve_isomorphism.tex \
			quartic.tex \
			kodaira_code.tex \
		LiDIA_ECO.tex \
			eco_prime.tex \
			meq_prime.tex \
			trace_mod.tex \
			trace_list.tex \
		LiDIA_GEC.tex \
			gec.tex \
			gec_complex_multiplication.tex



##	-----------
##	EXTRA STUFF
##	-----------

EXTRA_DIST = \
	${LiDIA_dvi_TEXSRCS} \
	LiDIA.bib \
	LiDIA.ist

# Unused CVS stuff
extra_stuff = \
	FIG/lidia_model.fig \
	.latex2html-init \
	doc2html.awk \
	filt.c \
	make_lidia_index \
	make_lidia_index.cmd \
	man_pages.tex \
	unused/bi_bf_lattice_basis.tex \
	unused/bi_bf_lattice_gensys.tex \
	unused/bigfloat_int.tex \
	unused/double_int.tex \
	unused/lattice_basis.tex \
	unused/lattice_gensys.tex \
	unused/use_vecmat.tex


##	-------------
##	SPECIAL RULES
##	-------------

# Convenience targets

.PHONY: dvi pdf dsc ps psgz psbz2

dvi: LiDIA.dvi
pdf: LiDIA.pdf
dsc: LiDIA.dsc
ps: LiDIA.ps
psgz: LiDIA.ps.gz
psbz2: LiDIA.ps.bz2


# Environment settings, Programs, and Flags
TEX_ENV = TEXINPUTS="${srcdir}${PATH_SEPARATOR}$$TEXINPUTS"
BIBTEX_ENV = BIBINPUTS="${srcdir}"
MAKEINDEX_ENV = INDEXSTYLE="${srcdir}"
DVIPS_ENV = TEXINPUTS="${srcdir}"

MAKEINDEX_OPTIONS = -s LiDIA.ist
DVIPS_OPTIONS = -ta4

# GNU texi2dvi supports LaTeX/bibtex/makeindex loops and adds the source
# dir to the search paths (except for BibTeX), thus is very practical here.
#
# don't capitalize LaTeX to work around Debian bug 472758
LATEX2DVI = $(TEXI2DVI) -l latex
LATEX2PDF = ${LATEX2DVI} -p

LATEX2DVI_ENV = \
	TEX="$(LATEX)" PDFTEX="$(PDFLATEX)" $(TEX_ENV) \
	BIBTEX="$(BIBTEX)" $(BIBTEX_ENV) \
	MAKEINDEX="$(MAKEINDEX) $(MAKEINDEX_OPTIONS)" $(MAKEINDEX_ENV)
LATEX2PDF_ENV = ${LATEX2DVI_ENV}

# What to put into the environment's GZIP/BZIP2 vars when using gzip/bzip2
GZIP_ENV = -9
BZIP_ENV =
BZIP2_ENV = -9

## Build rules

# For compatibility, $< is used in implicit rules only,
# and GNU Make's implicit pattern rules (%) are used only for
# targets that usually do not need to be (re)made.
# Whenever possible, specifications of implicit rules use the
# less flexible but more portable suffix form.

SUFFIXES = .tex .dvi .pdf .ps .dsc

.tex.dvi:
	${LATEX2DVI_ENV} ${LATEX2DVI} $<

.tex.pdf:
	${LATEX2PDF_ENV} ${LATEX2PDF} $<

.dvi.ps:
	${DVIPS_ENV} ${DVIPS} ${DVIPS_OPTIONS} -o $@ $<

.pdf.dsc:
	pdf2dsc $< $@

LiDIA.dvi LiDIA.pdf: ${LiDIA_dvi_TEXSRCS}

# Pattern rules for producing compressed targets.
# The given rules preserve the uncompressed file, thus unnecessary
# "make" invocations do not rebuild anything.
%.gz: %
	GZIP='$(GZIP_ENV)' gzip -c $< >$@

%.bz2: %
	BZIP='$(BZIP_ENV)' BZIP2='$(BZIP2_ENV)' bzip2 -c $< >$@

# The explicit compression rules for `LiDIA.ps', for non-GNU Makes
LiDIA.ps.gz: LiDIA.ps
	GZIP='$(GZIP_ENV)' gzip -c LiDIA.ps >$@

LiDIA.ps.bz2: LiDIA.ps
	BZIP='$(BZIP_ENV)' BZIP2='$(BZIP2_ENV)' bzip2 -c LiDIA.ps >$@

# This directory is not specific to a certain sub-package, therefore all its
# stuff goes into the base (or full) distribution and is removed from other
# distributions.
dist-hook:
	case ' $(PKG_DIRS) ' in \
	  *' base '*) ;; \
	  *) rm -rf $(distdir) ;; \
	esac

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

