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


#ifndef LIDIA_BIGINT_MATRIX_ALGORITHMS_H_GUARD_
#define LIDIA_BIGINT_MATRIX_ALGORITHMS_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class REP1, class REP2, class REP3 >
class bigint_matrix_algorithms
{

public:

	const REP1 rep_modul1;
	const REP2 rep_modul2;
	const REP3 rep_modul3;

	//
	// constructor
	//

	bigint_matrix_algorithms() {}

	//
	// destructor
	//

	~bigint_matrix_algorithms() {}

	//
	// divide
	//

	void divide(matrix< bigint > &RES, const matrix< bigint > &A,
		    const bigint &k) const;
	void compwise_divide(matrix< bigint > &RES, const matrix< bigint > &A,
			     const matrix< bigint > &B) const;

	//
	// remainder
	//

	void remainder(matrix< bigint > &, const matrix< bigint > &,
		       const bigint &) const;
	void remainder(matrix< long > &, const matrix< bigint > &, long) const;

	void trans_remainder(matrix< bigint > &, const matrix< bigint > &,
			     const bigint &) const;
	void trans_remainder(matrix< long > &, const matrix< bigint > &,
			     long) const;
};



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
class modular_bigint_matrix_algorithms
	: public bigint_matrix_algorithms< REP, REP, REP >
{

public:

	const REP rep_modul;
	const SINGLE_MODUL long_modul;
	const MULTI_MODUL bigint_modul;

	const BMA< bigint, DBMK < bigint >, DBMK< bigint > > dmodul;

	//
	// constructor
	//

	modular_bigint_matrix_algorithms() {}

	//
	// destructor
	//

	~modular_bigint_matrix_algorithms() {}

	//
	// Kernel
	//

	void kernel1(matrix< bigint > &, const matrix< bigint > &, const bigint &) const;
	void kernel2(matrix< bigint > &, const matrix< bigint > &) const;

	//
	// regular InvImage
	//

	void reginvimage1(matrix< bigint > &, const matrix< bigint > &,
			  const matrix< bigint > &, const bigint &) const;
	void reginvimage2(matrix< bigint > &, const matrix< bigint > &,
			  const matrix< bigint > &, const bigint &) const;

	//
	// Image
	//

	void image1(matrix< bigint > &, const matrix< bigint > &, const bigint &) const;
	void image2(matrix< bigint > &, const matrix< bigint > &) const;

	//
	// InvImage
	//

	void invimage(matrix< bigint > &, const matrix< bigint > &, const bigint *) const;
	void invimage(matrix< bigint > &, const matrix< bigint > &, const math_vector< bigint > &) const;

	//
	// Smith normal form
	//

	void snf_hartley(matrix< bigint > &) const;
	void snf_hartley(matrix< bigint > &, matrix< bigint > &, matrix< bigint > &) const;

	void snf_simple(matrix< bigint > &) const;
	void snf_simple(matrix< bigint > &, matrix< bigint > &, matrix< bigint > &) const;

	void snf_havas(matrix< bigint > &) const;
	void snf_havas(matrix< bigint > &, matrix< bigint > & T1, matrix< bigint > & T2) const;

	void snf_mult(matrix< bigint > &, long) const;
	void snf_mult(matrix< bigint > &, matrix< bigint > &, matrix< bigint > &, long) const;

	void snf_add(matrix< bigint > &, long) const;
	void snf_add(matrix< bigint > &, matrix< bigint > &, matrix< bigint > &, long) const;

	void snf_new(matrix< bigint > &, long) const;
	void snf_new(matrix< bigint > &, matrix< bigint > &, matrix< bigint > &, long) const;

	void snfmod_dkt(matrix< bigint > &, const bigint &) const;

	void snfmod_cohen(matrix< bigint > &, const bigint &) const;

	/////////////////////////
	// END: Linear algebra //
	// PART 2              //
	/////////////////////////

	void gauss(matrix< bigint > &) const;

	bigint *mgcd2(matrix< bigint > &, const bigint *, lidia_size_t) const;

	void basis_completion(matrix< bigint > &, bigint *, lidia_size_t) const;
	void cond_matrix(matrix< bigint > &, bigint *, lidia_size_t) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/bigint_matrix_algorithms.cc"
#endif



#endif	// LIDIA_BIGINT_MATRIX_ALGORITHMS_H_GUARD_
