##	========================================================
##
##	Automake definitions for the LiDIA FF package's examples
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
	-I$(LIDIA_FF_INCLUDEDIR) -I$(srcdir)/$(LIDIA_FF_INCLUDEDIR)


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	factorization_appl \
	Fp_polynomial_test_appl \
	Fp_rational_function_appl \
	dlp_appl \
	fpc_appl \
	Fp_pol_factor_appl \
	berlekamp_appl \
	can_zass_appl \
	gf_pol_factor_appl \
	bit_reverse_table_appl \
	fft_prime_appl \
	gf2n_appl \
	gf2n_find_root_appl \
	gf2n_poly_modulus_appl \
	gf2n_polynomial_appl \
	gf2n_rational_function_appl \
	gf_pol_appl \
	gfpn_appl \
	galois_field_iterator_appl \
	dense_power_series_appl \
	dpsr_udigit_mod_appl \
	sparse_power_series_appl

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

factorization_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/factorization/factorization_appl.cc

Fp_polynomial_test_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/Fp_polynomial/Fp_polynomial_test_appl.cc

Fp_rational_function_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/Fp_rational_function/Fp_rational_function_appl.cc

dlp_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/discrete_log/pohlig_hellman/dlp_appl.cc

fpc_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/discrete_log/pohlig_hellman/fpc_appl.cc

Fp_pol_factor_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/Fp_pol_factor_appl.cc

berlekamp_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/berlekamp_appl.cc

can_zass_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/factorization/Fp_polynomial/can_zass_appl.cc

gf_pol_factor_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/factorization/gf_polynomial/gf_pol_factor_appl.cc

bit_reverse_table_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/fft_prime/bit_reverse_table_appl.cc

fft_prime_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/fft_prime/fft_prime_appl.cc

gf2n_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gf2n/gf2n_appl.cc

gf2n_find_root_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_find_root_appl.cc

gf2n_poly_modulus_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_poly_modulus_appl.cc

gf2n_polynomial_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gf2n_polynomial/gf2n_polynomial_appl.cc

gf2n_rational_function_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gf2n_rational_function/gf2n_rational_function_appl.cc

gf_pol_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gf_polynomial/gf_pol_appl.cc

gfpn_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gfpn/gfpn_appl.cc

galois_field_iterator_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/gfpn/galois_field_iterator_appl.cc

dense_power_series_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/powser/dense_power_series_appl.cc

dpsr_udigit_mod_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/powser/dpsr_udigit_mod_appl.cc

sparse_power_series_appl_SOURCES = \
	$(LIDIA_FF_SRCDIR)/powser/sparse_power_series_appl.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	${LIDIA_FF_SRCDIR}/Fp_polynomial/Fp_poly_modulus_timing.cc \
	${LIDIA_FF_SRCDIR}/fft_prime/test_mul_mod.cc \
	${LIDIA_FF_SRCDIR}/gf2n/gf2n_timing.cc \
	${LIDIA_FF_SRCDIR}/gf2n_polynomial/gf2n_poly_modulus_appl.cc.asli \
	${LIDIA_FF_SRCDIR}/gf2n_polynomial/gf2n_poly_modulus_timing.cc \
	${LIDIA_FF_SRCDIR}/gf2n_polynomial/gf2n_polynomial_timing.cc \
	${LIDIA_FF_SRCDIR}/gfpn/gfpn_timings_appl.cc


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

