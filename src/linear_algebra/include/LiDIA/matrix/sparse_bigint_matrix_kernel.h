// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_H_GUARD_
#define LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_H_GUARD_


#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
# include	"LiDIA/random_generator.h"
#endif
#ifndef LIDIA_SPARSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_ring_matrix_kernel.h"
#endif
#ifndef LIDIA_CRT_H_GUARD_
# include	"LiDIA/crt.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class MATRIX_TYPE >
class sparse_bigint_matrix_kernel : public SRMK< bigint >
{
	//
	// modul definitions
	//

	const BMA< bigint, SBMK < bigint >, SBMK< bigint > > SS_base_modul;

	//
	// constructor
	//

public:

	sparse_bigint_matrix_kernel() {}

	//
	// destructor
	//

public:

	~sparse_bigint_matrix_kernel() {}

	//
	// divide
	//

	void divide(MATRIX_TYPE &RES, const MATRIX_TYPE &A, const bigint &k) const;
	void compwise_divide(MATRIX_TYPE &RES, const MATRIX_TYPE &A, const MATRIX_TYPE &B) const;

	//
	// remainder
	//

	void remainder(MATRIX_TYPE &, const MATRIX_TYPE &, const bigint &) const;

	//
	// norms and bounds
	//

	void max(const MATRIX_TYPE &, bigint &) const;
	void max_abs(const MATRIX_TYPE &, bigint &) const;
	void max_pos(const MATRIX_TYPE &, bigint &, lidia_size_t &, lidia_size_t &) const;
	void max_abs_pos(const MATRIX_TYPE &, bigint &, lidia_size_t &, lidia_size_t &) const;

	void min(const MATRIX_TYPE &, bigint &) const;
	void min_abs(const MATRIX_TYPE &, bigint &) const;
	void min_pos(const MATRIX_TYPE &, bigint &, lidia_size_t &, lidia_size_t &) const;
	void min_abs_pos(const MATRIX_TYPE &, bigint &, lidia_size_t &, lidia_size_t &) const;

	void hadamard(const MATRIX_TYPE &, bigint &) const;
	void binary_hadamard(const MATRIX_TYPE &, lidia_size_t &) const;

	void row_norm(const MATRIX_TYPE &, bigint &, lidia_size_t, long) const;

	void column_norm(const MATRIX_TYPE &, bigint &, lidia_size_t, long) const;

	//
	// randomize
	//

	void randomize(MATRIX_TYPE &, const bigint &) const;
	void randomize_with_det(MATRIX_TYPE &, const bigint &, const bigint &) const;
	void randomize(MATRIX_TYPE &, const bigint &, const long) const;

	//
	// regular expansion
	//

	void regexpansion(MATRIX_TYPE &, const lidia_size_t *) const;

	///////////////////////////
	// BEGIN: Linear algebra //
	// PART 2                //
	///////////////////////////

	//
	// Hermite normal form
	//

	void hnfmod_dkt(MATRIX_TYPE &, const bigint &) const;

	void hnfmod_cohen(MATRIX_TYPE &, const bigint &) const;

	void hnfmod_mueller(MATRIX_TYPE &, MATRIX_TYPE &, bigint &) const;

	void hnf_storjohann(MATRIX_TYPE &RES) const;
	void hnf_storjohann(MATRIX_TYPE &RES, MATRIX_TYPE &TR, MATRIX_TYPE &C, MATRIX_TYPE &Q) const;

	//
	// Kernel
	//

	void kernel1(MATRIX_TYPE &, const MATRIX_TYPE &) const;
	void kernel2(MATRIX_TYPE &, const MATRIX_TYPE &) const;

	//
	// regular InvImage
	//

	void reginvimage1(MATRIX_TYPE &, const MATRIX_TYPE &, const MATRIX_TYPE &) const;
	void reginvimage2(MATRIX_TYPE &, const MATRIX_TYPE &, const MATRIX_TYPE &) const;

	//
	// Image
	//

	void image1(MATRIX_TYPE &, const MATRIX_TYPE &) const;
	void image2(MATRIX_TYPE &, const MATRIX_TYPE &) const;

	//
	// InvImage
	//

	void invimage(MATRIX_TYPE &, const MATRIX_TYPE &, const bigint *) const;
	void invimage(MATRIX_TYPE &, const MATRIX_TYPE &, const math_vector< bigint > &) const;

	//
	// Smith normal form
	//

	void snf_hartley(MATRIX_TYPE &) const;
	void snf_hartley(MATRIX_TYPE &, MATRIX_TYPE &, MATRIX_TYPE &) const;

	void snf_simple(MATRIX_TYPE &) const;
	void snf_simple(MATRIX_TYPE &, MATRIX_TYPE &, MATRIX_TYPE &) const;

	void snf_havas(MATRIX_TYPE &) const;
	void snf_havas(MATRIX_TYPE &, MATRIX_TYPE & T1, MATRIX_TYPE & T2) const;

	void snf_mult(MATRIX_TYPE &, long) const;
	void snf_mult(MATRIX_TYPE &, MATRIX_TYPE &, MATRIX_TYPE &, long) const;

	void snf_add(MATRIX_TYPE &, long) const;
	void snf_add(MATRIX_TYPE &, MATRIX_TYPE &, MATRIX_TYPE &, long) const;

	void snf_new(MATRIX_TYPE &, long) const;
	void snf_new(MATRIX_TYPE &, MATRIX_TYPE &, MATRIX_TYPE &, long) const;

	void snfmod_dkt(MATRIX_TYPE &, const bigint &) const;

	void snfmod_cohen(MATRIX_TYPE &, const bigint &) const;

	/////////////////////////
	// END: Linear algebra //
	// PART 2              //
	/////////////////////////

	void gauss(MATRIX_TYPE &) const;

	bigint *mgcd2(MATRIX_TYPE &, const bigint *, lidia_size_t) const;

	void basis_completion(MATRIX_TYPE &, bigint *, lidia_size_t) const;
	void simple_basis_completion(MATRIX_TYPE &, bigint *, lidia_size_t) const;
	lidia_size_t cond_matrix(MATRIX_TYPE &, bigint *, lidia_size_t) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/sparse_bigint_matrix_kernel.cc"
#endif



#endif	// LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_H_GUARD_
