##	==========================================================
##
##	Automake definitions for the LiDIA base package's examples
##
##	==========================================================

AUTOMAKE_OPTIONS = 1.5 nostdinc

ETAGS_ARGS = -C -l c++

include $(top_srcdir)/Paths.def


##	-------------
##	INCLUDE PATHS
##	-------------

AM_CPPFLAGS = \
	-I$(LIDIA_CONFIG_INCLUDEDIR) \
	-I$(LIDIA_BASE_INCLUDEDIR) -I$(srcdir)/$(LIDIA_BASE_INCLUDEDIR)


##	--------
##	PROGRAMS
##	--------

applexecdir = $(pkglibdir)/examples

applexec_PROGRAMS = \
	udigit_appl \
	random_generator_appl \
	signal_lidia_appl \
	timer_appl \
	bigcomplex_appl \
	bigfloat_appl \
	bigfloat_int_appl \
	bigint_appl \
	cornacchia_appl \
	is_power_appl \
	is_prime_appl \
	is_square_appl \
	ressol_appl \
	bigmod_appl \
	multi_bigmod_appl \
	bigrational_appl \
	crt_appl \
	rational_factorization_appl \
	divisors_appl \
	prime_list_appl \
	prime_list_bench_appl \
	udigit_mod_appl \
	xbigfloat_appl \
	xdouble_appl \
	bigint_factor_appl \
	jenkins_traub_appl \
	polynomial_cast_appl \
	resultant_check_appl \
	hash_table_appl \
	indexed_hash_table_appl \
	vector_appl \
	base_matrix_appl \
	power_functions_appl

$(applexec_PROGRAMS): $(LIDIA_LTLIB)
LDADD = $(LIDIA_LTLIB)


##	---------------
##	PROGRAM SOURCES
##	---------------

nodist_udigit_appl_SOURCES = \
	$(LIDIA_CONFIG_LINKDIR)/udigit/udigit_appl.cc

random_generator_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/system/random_generator_appl.cc

signal_lidia_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/system/signal_lidia_appl.cc

timer_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/system/timer_appl.cc

bigcomplex_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigcomplex/bigcomplex_appl.cc

bigfloat_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_appl.cc

bigfloat_int_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat_int/bigfloat_int_appl.cc

bigint_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/bigint_appl.cc

cornacchia_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/cornacchia_appl.cc

is_power_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/is_power_appl.cc

is_prime_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/is_prime_appl.cc

is_square_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/is_square_appl.cc

ressol_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/ressol_appl.cc

bigmod_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigmod/bigmod_appl.cc

multi_bigmod_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigmod/multi_bigmod_appl.cc

bigrational_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigrational/bigrational_appl.cc

crt_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/chinese_rem/crt_appl.cc

rational_factorization_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/rational_factorization_appl.cc

divisors_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/nmbrthry_functions/divisors_appl.cc

prime_list_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/prime_list/prime_list_appl.cc

prime_list_bench_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/prime_list/prime_list_bench_appl.cc

udigit_mod_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/udigit_mod/udigit_mod_appl.cc

xbigfloat_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/xbigfloat/xbigfloat_appl.cc

xdouble_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/simple_classes/xdouble/xdouble_appl.cc

bigint_factor_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/bigint_factor_appl.cc

jenkins_traub_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/polynomial/jenkins_traub_appl.cc

polynomial_cast_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/polynomial/polynomial_cast_appl.cc

resultant_check_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/polynomial/resultant_check_appl.cc

hash_table_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/hash_table/hash_table_appl.cc

indexed_hash_table_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/hash_table/indexed_hash_table_appl.cc

vector_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/vector/vector_appl.cc

base_matrix_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/matrix/base_matrix_appl.cc

power_functions_appl_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/misc/power_functions_appl.cc


##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	$(LIDIA_CONFIG_SRCDIR)/udigit/cln/udigit_appl.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/gmp/udigit_appl.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/libI/udigit_appl.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/piologie/udigit_appl.cc \
	\
	${LIDIA_BASE_SRCDIR}/hash_table/hash_table_appl_int.dat \
	${LIDIA_BASE_SRCDIR}/hash_table/indexed_hash_table_appl_int.dat \
	${LIDIA_BASE_SRCDIR}/matrix/dense_input_bigint \
	${LIDIA_BASE_SRCDIR}/matrix/dense_input_int \
	${LIDIA_BASE_SRCDIR}/matrix/dense_input_long \
	${LIDIA_BASE_SRCDIR}/matrix/dense_ring_matrix_appl_bigint.dat \
	${LIDIA_BASE_SRCDIR}/matrix/dense_ring_matrix_appl_int.dat \
	${LIDIA_BASE_SRCDIR}/matrix/dense_ring_matrix_appl_long.dat \
	${LIDIA_BASE_SRCDIR}/matrix/input_bigfloat \
	${LIDIA_BASE_SRCDIR}/matrix/input_bigint \
	${LIDIA_BASE_SRCDIR}/matrix/input_double \
	${LIDIA_BASE_SRCDIR}/matrix/input_int \
	${LIDIA_BASE_SRCDIR}/matrix/input_long \
	${LIDIA_BASE_SRCDIR}/matrix/matrix_inst \
	${LIDIA_BASE_SRCDIR}/matrix/sparse_input_bigint \
	${LIDIA_BASE_SRCDIR}/matrix/sparse_input_int \
	${LIDIA_BASE_SRCDIR}/matrix/sparse_input_long \
	${LIDIA_BASE_SRCDIR}/matrix/sparse_ring_matrix_appl_bigint.dat \
	${LIDIA_BASE_SRCDIR}/matrix/sparse_ring_matrix_appl_int.dat \
	${LIDIA_BASE_SRCDIR}/matrix/sparse_ring_matrix_appl_long.dat \
	${LIDIA_BASE_SRCDIR}/polynomial/test_roots.1 \
	${LIDIA_BASE_SRCDIR}/polynomial/test_roots.1.cmp \
	${LIDIA_BASE_SRCDIR}/polynomial/test_roots_appl.cc \
	${LIDIA_BASE_SRCDIR}/simple_classes/bigint/test_bigint_io.cc


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

