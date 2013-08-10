##	=============================================================
##
##	Automake definitions for the LiDIA FF package's library stuff
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
	-I$(LIDIA_FF_INCLUDEDIR) -I$(srcdir)/$(LIDIA_FF_INCLUDEDIR)


##	---------
##	LIBRARIES
##	---------

noinst_LTLIBRARIES = libLiDIA-FF.la

# We don't want to have $(LIBS) in the convenience library.
# FIXME: This is not safe against command-line overrides.
LIBS =


##	---------------
##	LIBRARY SOURCES
##	---------------

libLiDIA_FF_la_SOURCES = \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/Fp_pol_arithmetic.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/Fp_pol_tools.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/Fp_poly_modulus.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/Fp_poly_multiplier.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/Fp_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/compose.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/crossover_class.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/factoring.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/fft_arith.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/fft_reps.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/fractions.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/gcd.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/plain_arith.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/poly_argument.cc \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/poly_matrix.cc \
	$(LIDIA_FF_SRCDIR)/Fp_rational_function/Fp_rational_function.cc \
	$(LIDIA_FF_SRCDIR)/discrete_log/pohlig_hellman/dlp.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/Fp_factor.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/Fp_factor_generic.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/berlekamp.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/binomial.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/canzass.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/ddf.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/edf.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/find_factors.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/old_ddf.cc \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/squarefree.cc \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/checked_min_poly.cc \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/gf_berlekamp.cc \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/gf_can_zass.cc \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/gf_fact_util.cc \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/gf_factor.cc \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/gf_squarefree.cc \
	$(LIDIA_FF_SRCDIR)/factorization/sf_Fp_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/factorization/sf_gf_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/fft_prime/bit_reverse_table.cc \
	$(LIDIA_FF_SRCDIR)/fft_prime/fft_prime.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2nIO.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n_init.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n_invert.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n_lib.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n_mult.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n_reduce.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2nin.cc \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2nout.cc \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_find_root.cc \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_poly_karatzuba.cc \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_poly_modulus.cc \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_poly_remainder.cc \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/gf2n_rational_function/gf2n_rational_function.cc \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_pol_arith.cc \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_pol_instant.cc \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_pol_util.cc \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_poly_argument.cc \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_poly_modulus.cc \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/galois_field.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/galois_field_iterator.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/galois_field_rep.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/gf_element.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/gf_rep_GF2n.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/gf_rep_GFp.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/gf_rep_GFpn.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/info_gf2n.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/info_gf2n_inv.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/info_gf2n_mul.cc \
	$(LIDIA_FF_SRCDIR)/gfpn/info_gf2n_red.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bpf_gf2n.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bpf_gf_element.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_Fp_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_gf2n.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_gf_element.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_pp_sf_Fp_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_pp_sf_gf_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_spc_bigint.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_spc_bigmod.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_spc_bigrational.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/bv_spc_udigit_mod.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/dpsr_bigint.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/dpsr_bigrational.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/fact_Fp_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/fact_gf_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/lpf_gf2n.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/spsr_bigint.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/spsr_bigrational.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/spsr_udigit_mod.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/sv_pp_sf_Fp_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/sv_pp_sf_gf_polynomial.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/sv_spc_bigint.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/sv_spc_bigmod.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/sv_spc_bigrational.cc \
	$(LIDIA_FF_SRCDIR)/instantiate/sv_spc_udigit_mod.cc \
	$(LIDIA_FF_SRCDIR)/powser/dpsr_bigmod.cc \
	$(LIDIA_FF_SRCDIR)/powser/dpsr_udigit_mod.cc \
	$(LIDIA_FF_SRCDIR)/powser/spsr_bigmod.cc


##	-------
##	HEADERS
##	-------

instantincludedir = $(pkgincludedir)/instantiate
specialincludedir = $(pkgincludedir)/specialization
ffincludedir = $(pkgincludedir)/finite_fields

instantinclude_HEADERS = \
	$(LIDIA_FF_PKGINCLUDEDIR)/instantiate/power_series.cc

specialinclude_HEADERS = \
	$(LIDIA_FF_PKGINCLUDEDIR)/specialization/Fp_polynomial.special \
	$(LIDIA_FF_PKGINCLUDEDIR)/specialization/dense_power_series.special \
	$(LIDIA_FF_PKGINCLUDEDIR)/specialization/gf_element.special \
	$(LIDIA_FF_PKGINCLUDEDIR)/specialization/gf_polynomial.special

ffinclude_HEADERS = \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/Fp_polynomial_fft.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/Fp_polynomial_util.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/base_dense_power_series.cc \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/base_dense_power_series.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/base_sparse_power_series.cc \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/base_sparse_power_series.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/bit_reverse_table.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/coeff_sparse_power_series.cc \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/coeff_sparse_power_series.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/dense_power_series_bigmod.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/dense_power_series_udigit_mod.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/fft_mul_mod.inl \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/galois_field_rep.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/gf2nIO.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/info_gf2n.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/sf_Fp_polynomial.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/finite_fields/sf_gf_polynomial.h

pkginclude_HEADERS = \
	$(LIDIA_FF_PKGINCLUDEDIR)/Fp_pol_crossover.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/Fp_poly_modulus.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/Fp_poly_multiplier.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/Fp_polynomial.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/Fp_rational_function.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/dense_power_series.cc \
	$(LIDIA_FF_PKGINCLUDEDIR)/dense_power_series.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/dlp.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/fft_prime.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/galois_field.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/galois_field_iterator.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/gf2n.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/gf2n_poly_modulus.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/gf2n_polynomial.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/gf2n_rational_function.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/gf_element.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/gf_polynomial.h \
	$(LIDIA_FF_PKGINCLUDEDIR)/sparse_power_series.cc \
	$(LIDIA_FF_PKGINCLUDEDIR)/sparse_power_series.h


##	----
##	DATA
##	----

dist_pkgdata_DATA = \
	$(LIDIA_PKGDATADIR)/GF2n.database


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	${LIDIA_FF_SRCDIR}/Fp_polynomial/crossover.tbl \
	${LIDIA_FF_SRCDIR}/Fp_polynomial/crossover.tbl.default \
	${LIDIA_FF_SRCDIR}/Fp_polynomial/make_crossover.cc \
	${LIDIA_FF_SRCDIR}/gf2n_polynomial/crossover.tbl \
	${LIDIA_FF_SRCDIR}/gf2n_polynomial/crossover.tbl.default \
	${LIDIA_FF_SRCDIR}/gf2n_polynomial/make_crossover.cc


##	-------------
##	SPECIAL RULES
##	-------------

# Remove installation subdirectories after uninstalling, unless left non-empty.
uninstall-hook:
	rmdir	$(DESTDIR)$(instantincludedir) \
		$(DESTDIR)$(specialincludedir) \
		$(DESTDIR)$(ffincludedir) \
		$(DESTDIR)$(pkgincludedir) \
		$(DESTDIR)$(pkgdatadir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

