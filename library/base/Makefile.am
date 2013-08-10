##	===============================================================
##
##	Automake definitions for the LiDIA base package's library stuff
##
##	===============================================================

AUTOMAKE_OPTIONS = 1.5 nostdinc

ETAGS_ARGS = -C -l c++

include $(top_srcdir)/Paths.def


##	-------------
##	INCLUDE PATHS
##	-------------

AM_CPPFLAGS = \
	-I$(LIDIA_CONFIG_INCLUDEDIR) \
	-I$(LIDIA_BASE_INCLUDEDIR) -I$(srcdir)/$(LIDIA_BASE_INCLUDEDIR)


##	---------
##	LIBRARIES
##	---------

noinst_LTLIBRARIES = libLiDIA-base.la

# Add the configured portability object files
libLiDIA_base_la_LIBADD = @LIDIA_BASE_XLTLIBOBJS@
libLiDIA_base_la_DEPENDENCIES = @LIDIA_BASE_XLTLIBOBJS@

# We don't want to have $(LIBS) in the convenience library.
# FIXME: This is not safe against command-line overrides.
LIBS =


##	---------------
##	LIBRARY SOURCES
##	---------------

nodist_libLiDIA_base_la_SOURCES = \
	$(LIDIA_CONFIG_LINKDIR)/bigint/bigint.cc \
	$(LIDIA_CONFIG_LINKDIR)/udigit/udigit.cc

libLiDIA_base_la_SOURCES = \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_check_gcd.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_ecm.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_factor.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_fermat_method.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_mpqs.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_pollard_p-1.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_pollard_rho.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_trialdiv.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/bigint/sf_williams_p+1.cc \
	$(LIDIA_BASE_SRCDIR)/factorization/sf_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_bigint_dense.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_bigint_normal.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_bigint_sparse.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_bigint_sub.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_double_dense.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_double_normal.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_double_sparse.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_double_sub.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_int_dense.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_int_normal.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_int_sparse.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_int_sub.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_long_dense.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_long_normal.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_long_sparse.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bm_long_sub.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bpf_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bpf_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_bigcomplex.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_bigfloat.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_bigrational.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_char.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_double.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_int.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_long.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_param_desc_ptr.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_pp_sf_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_rf_single_factor.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_udigit_mod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/bv_xdouble.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/fact_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/fv_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/fv_bigrational.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/fv_double.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/fv_int.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/fv_long.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/ht_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/ht_int.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/iht_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/iht_int.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/lpf_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/lpf_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_bigfloat.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_bigrational.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_double.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_int.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_long.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_udigit_mod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/mv_xdouble.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/rcl_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/spf_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/spf_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_bigrational.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_double.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_int.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_long.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_pp_sf_bigint.cc \
	$(LIDIA_BASE_SRCDIR)/instantiate/sv_rf_single_factor.cc \
	$(LIDIA_BASE_SRCDIR)/matrix/base_matrix_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/matrix/matrix_error_msg.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/bigcomplex_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/bigfloat_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/bigint_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/bigrational_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/cast_bigcomplex_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/cast_bigfloat_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/cast_bigrational_polynomial.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/jenkins_traub.cc \
	$(LIDIA_BASE_SRCDIR)/polynomial/poly_error_msg.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigcomplex/bigcomplex.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigcomplex/bigcomplex_io.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigcomplex/bigcomplex_trans.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat2.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_const.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_exp.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_hyp.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_io.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_miscfunc.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat/bigfloat_trig.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat_int/bigfloat_int.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigfloat_int/bigfloat_int_trans.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/bigint_share.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/chinese_remainder.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/cornacchia.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/fermat.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/is_power.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/is_prime.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/is_square.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/jacobi.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/multiple_gcd.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/nearest.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/newton_root.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/next_prime.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/power_mod.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/previous_prime.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigint/ressol.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigmod/base_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigmod/bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigmod/bigmod_lib.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigmod/multi_bigmod.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigrational/bigrational.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/bigrational/bigrational_io.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/chinese_rem/crt.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/chinese_rem/crt_table.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/ec_point.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/ecm.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/ecm_primes.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/file_f.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/mpqs.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/rational_factorization.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/factorization/trialdiv.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/lanczos/lanczos.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/lanczos/lanczos_small_matrix.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/lanczos/lanczos_spmatrix.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/lanczos/lanczos_spvector.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/lanczos/lanczos_vblock.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/lanczos/process.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/matrix_GL2Z.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/nmbrthry_functions/divisors.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/nmbrthry_functions/modular_functions.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/prime_list/prime_list.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/single_precision/single_precision.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/udigit_mod/udigit_mod.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/xbigfloat/xbigfloat.cc \
	$(LIDIA_BASE_SRCDIR)/simple_classes/xdouble/xdouble.cc \
	$(LIDIA_BASE_SRCDIR)/system/debug.cc \
	$(LIDIA_BASE_SRCDIR)/system/error.cc \
	$(LIDIA_BASE_SRCDIR)/system/file_lidia.cc \
	$(LIDIA_BASE_SRCDIR)/system/info.cc \
	$(LIDIA_BASE_SRCDIR)/system/memory.cc \
	$(LIDIA_BASE_SRCDIR)/system/osstream.cc \
	$(LIDIA_BASE_SRCDIR)/system/param_desc.cc \
	$(LIDIA_BASE_SRCDIR)/system/param_desc_instant.cc \
	$(LIDIA_BASE_SRCDIR)/system/precondition_error.cc \
	$(LIDIA_BASE_SRCDIR)/system/precondition_error_instant.cc \
	$(LIDIA_BASE_SRCDIR)/system/random_generator.cc \
	$(LIDIA_BASE_SRCDIR)/system/signal_lidia.cc \
	$(LIDIA_BASE_SRCDIR)/system/timer.cc \
	$(LIDIA_BASE_SRCDIR)/system/warning.cc \
	$(LIDIA_BASE_SRCDIR)/vector/vector_error_msg.cc

EXTRA_libLiDIA_base_la_SOURCES = \
	$(LIDIA_PORTAB_SRCDIR)/stlport.cc \
	$(LIDIA_PORTAB_SRCDIR)/random.c \
	$(LIDIA_PORTAB_SRCDIR)/random.h


##	-------
##	HEADERS
##	-------

kernelincludedir = $(pkgincludedir)/kernel
instantincludedir = $(pkgincludedir)/instantiate
specialincludedir = $(pkgincludedir)/specialization
matrixincludedir = $(pkgincludedir)/matrix
baseincludedir = $(pkgincludedir)/base

instantinclude_HEADERS = \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/factorization.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/hash_table.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/pair.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/power_functions.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/residue_class_list.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/instantiate/vector.cc

specialinclude_HEADERS = \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/bigcomplex.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/bigfloat.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/bigint.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/bigmod.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/bigrational.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/polynomial.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/single_factor.special \
	$(LIDIA_BASE_PKGINCLUDEDIR)/specialization/udigit_mod.special

matrixinclude_HEADERS = \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/base_matrix_algorithms.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/base_matrix_algorithms.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/dense_base_matrix_kernel.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/dense_base_matrix_kernel.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/dense_field_matrix_kernel.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/dense_field_matrix_kernel.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/dense_ring_matrix_kernel.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/dense_ring_matrix_kernel.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/field_matrix_algorithms.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/field_matrix_algorithms.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/matrix_representation.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/ring_matrix_algorithms.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/ring_matrix_algorithms.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/sparse_base_matrix_kernel.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/sparse_base_matrix_kernel.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/sparse_field_matrix_kernel.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/sparse_field_matrix_kernel.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/sparse_ring_matrix_kernel.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix/sparse_ring_matrix_kernel.h

baseinclude_HEADERS = \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/array_functions.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/b_value.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/base_bigmod.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/base_factor.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/base_factor.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/ecm_primes.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/file_io_class.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/interface_lib.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/poly_intern.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/poly_intern.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/ppair.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/ppair.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/residue_class_list.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/residue_class_list.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base/sf_bigint.h

pkginclude_HEADERS = \
	$(LIDIA_BASE_PKGINCLUDEDIR)/arith.inl \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base_vector.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/base_vector.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigcomplex.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigcomplex_polynomial.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigfloat.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigfloat_config.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigfloat_int.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigfloat_polynomial.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigint.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigint_polynomial.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigmod.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigrational.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/bigrational_polynomial.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/comparator.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/crt.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/crt_table.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/debug.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/dense_base_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/dense_base_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/dense_field_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/dense_field_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/dense_ring_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/dense_ring_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/error.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/factorization.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/factorization.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/field_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/field_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/field_polynomial.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/field_polynomial.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/file_vector.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/file_vector.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/hash_table.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/hash_table.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/indexed_hash_table.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/indexed_hash_table.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/info.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/isstream.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/lanczos.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/lidia_defines.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/lidia_file.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/lidia_reference_counter.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/lidia_signal.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/lidia_vector.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/math_vector.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/math_vector.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix_GL2Z.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/matrix_flags.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/memory.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/modular_functions.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/modular_operations.inl \
	$(LIDIA_BASE_PKGINCLUDEDIR)/mpqs_timing.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/multi_bigmod.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/nmbrthry_functions.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/osstream.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/pair.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/pair.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/polynomial.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/power_functions.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/power_functions.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/param_desc.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/param_desc.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/precondition_error.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/precondition_error.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/prime_list.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/random_generator.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/rational_factorization.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/ring_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/ring_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/single_factor.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/single_factor.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sort_vector.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sort_vector.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sparse_base_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sparse_base_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sparse_field_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sparse_field_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sparse_ring_matrix.cc \
	$(LIDIA_BASE_PKGINCLUDEDIR)/sparse_ring_matrix.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/timer.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/udigit.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/udigit_mod.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/vector_flags.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/vector_representation.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/warning.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/xbigfloat.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/xdouble.h \
	$(LIDIA_BASE_PKGINCLUDEDIR)/xdouble.inl \
	$(LIDIA_BASE_PKGINCLUDEDIR)/xerror.h


# The following file lists are configuration-dependent.
# Therefore, wildcards are used.  These don't work with VPATH,
# but this stuff really is in the build tree.  That's why it works.

# Common members:
#	$(LIDIA_CONFIG_PKGINCLUDEDIR)/LiDIA.h \
#	$(LIDIA_CONFIG_PKGINCLUDEDIR)/path.h
nodist_pkginclude_HEADERS = \
	$(LIDIA_CONFIG_PKGINCLUDEDIR)/*.h

# Common members:
#	$(LIDIA_CONFIG_PKGINCLUDEDIR)/kernel/bigint_def.h \
#	$(LIDIA_CONFIG_PKGINCLUDEDIR)/kernel/bigint_interface.h \
#	$(LIDIA_CONFIG_PKGINCLUDEDIR)/kernel/udigit_def.h \
#	$(LIDIA_CONFIG_PKGINCLUDEDIR)/kernel/udigit_interface.h
nodist_kernelinclude_HEADERS = \
	$(LIDIA_CONFIG_PKGINCLUDEDIR)/kernel/*.h

##      ------------------
##      BUILT DEPENDENCIES
##      ------------------

BUILT_SOURCES = $(LIDIA_CONFIG_INCLUDEDIR)/bytes_to_int_flag.h
CLEANFILES = $(LIDIA_CONFIG_INCLUDEDIR)/bytes_to_int_flag.h

noinst_HEADERS = $(LIDIA_CONFIG_INCLUDEDIR)/bytes_to_int_flag.h

$(LIDIA_CONFIG_INCLUDEDIR)/bytes_to_int_flag.h: bytes_to_int_flag_generator$(EXEEXT)
	./bytes_to_int_flag_generator > $@

noinst_PROGRAMS = bytes_to_int_flag_generator$(EXEEXT)
bytes_to_int_flag_generator_SOURCES = \
	$(LIDIA_PORTAB_SRCDIR)/bytes_to_int_flag_generator.cc

##	----------------
##	ADDITIONAL STUFF
##	----------------

EXTRA_DIST = \
	$(LIDIA_CONFIG_SRCDIR)/bigint/cln/bigint.cc \
	$(LIDIA_CONFIG_SRCDIR)/bigint/cln/bigint_def.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/cln/bigint_interface.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/gmp/bigint.cc \
	$(LIDIA_CONFIG_SRCDIR)/bigint/gmp/bigint_def.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/gmp/bigint_interface.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/libI/bigint.cc \
	$(LIDIA_CONFIG_SRCDIR)/bigint/libI/bigint_def.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/libI/bigint_interface.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/piologie/bigint.cc \
	$(LIDIA_CONFIG_SRCDIR)/bigint/piologie/bigint_def.h \
	$(LIDIA_CONFIG_SRCDIR)/bigint/piologie/bigint_interface.h \
	\
	$(LIDIA_CONFIG_SRCDIR)/udigit/cln/udigit.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/cln/udigit_def.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/cln/udigit_interface.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/gmp/udigit.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/gmp/udigit128.gmp.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/gmp/udigit_def.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/gmp/udigit_interface.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/libI/udigit.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/libI/udigit_def.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/libI/udigit_interface.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/piologie/udigit.cc \
	$(LIDIA_CONFIG_SRCDIR)/udigit/piologie/udigit_def.h \
	$(LIDIA_CONFIG_SRCDIR)/udigit/piologie/udigit_interface.h \
	\
	$(LIDIA_CONFIG_SRCDIR)/memory/boehm/gmm.h \
	$(LIDIA_CONFIG_SRCDIR)/memory/boehm/makefile \
	$(LIDIA_CONFIG_SRCDIR)/memory/boehm/new.c \
	$(LIDIA_CONFIG_SRCDIR)/memory/default/makefile \
	$(LIDIA_CONFIG_SRCDIR)/memory/dentzer/gmm.h \
	$(LIDIA_CONFIG_SRCDIR)/memory/dentzer/makefile \
	$(LIDIA_CONFIG_SRCDIR)/memory/dentzer/new.c \
	$(LIDIA_CONFIG_SRCDIR)/memory/malloc/gmm.h \
	$(LIDIA_CONFIG_SRCDIR)/memory/malloc/makefile \
	$(LIDIA_CONFIG_SRCDIR)/memory/malloc/new.c \
	\
	${LIDIA_BASE_SRCDIR}/simple_classes/factorization/test_machine.cc

# Unused CVS stuff
extra_stuff = \
	${LIDIA_BASE_SRCDIR}/matrix/ring_matrix_bigmod.cc.unfinished


##	-------------
##	SPECIAL RULES
##	-------------

# Remove installation subdirectories after uninstalling, unless left non-empty.
uninstall-hook:
	rmdir	$(DESTDIR)$(kernelincludedir) \
		$(DESTDIR)$(instantincludedir) \
		$(DESTDIR)$(specialincludedir) \
		$(DESTDIR)$(matrixincludedir) \
		$(DESTDIR)$(baseincludedir) \
		$(DESTDIR)$(pkgincludedir) \
		2>/dev/null; :

# Delete half-made targets when their commands fail
.DELETE_ON_ERROR:

