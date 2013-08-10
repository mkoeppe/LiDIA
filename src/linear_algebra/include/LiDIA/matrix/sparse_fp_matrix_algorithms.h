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


#ifndef LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_H_GUARD_
#define LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_H_GUARD_


#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
# include	"LiDIA/random_generator.h"
#endif
#ifndef LIDIA_POLYNOMIAL_H_GUARD_
# include	"LiDIA/polynomial.h"
#endif
#ifndef LIDIA_ARITH_INL_GUARD_
# include	"LiDIA/arith.inl"
#endif
#ifndef LIDIA_MODULAR_OPERATIONS_INL_GUARD_
# include	"LiDIA/modular_operations.inl"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_BASE_POWER_PRODUCT_H_GUARD_
# include	"LiDIA/base_power_product.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T, class MATRIX_TYPE >
class sparse_fp_matrix_algorithms
{
protected:

	//
	// internal functions
	//

	T * berlekamp_massay(T *s, lidia_size_t n, const T &mod) const;
	const T multiply_special(const math_vector< T > &w,
				 const math_vector< T > &v,
				 const T &mod) const;
	void multiply_special(math_vector< T > &w,
			      const base_power_product< ring_matrix < T >, lidia_size_t > &B,
			      const T &mod) const;
	void multiply_special(math_vector< T > &w,
			      const MATRIX_TYPE &B,
			      const T &mod) const;
	void multiply_special(math_vector< T > &v,
			      const math_vector< T > &w,
			      const MATRIX_TYPE &B,
			      const T &mod) const;

public:

	//
	// constructor
	//

	sparse_fp_matrix_algorithms() {}

	//
	// destructor
	//

	~sparse_fp_matrix_algorithms() {}

	//
	// column step form
	//

	int STF(MATRIX_TYPE &, const T &) const;

	//
	// extended column step form
	//

	const T * STF_extended(MATRIX_TYPE &, const T &) const
	{
		return NULL;
	}

	//
	// rank
	//

	lidia_size_t rank(MATRIX_TYPE &, const T &) const;

	//
	// rank and linearly independent rows or columns
	//

	lidia_size_t *lininr(MATRIX_TYPE &, const T &) const;
	lidia_size_t *lininc(MATRIX_TYPE &, const T &) const;

	//
	// adjoint matrix
	//

	void adj(MATRIX_TYPE &, const T &) const;

	//
	// determinant
	//

	const T det(MATRIX_TYPE &A, const T &mod) const;
	const T det(const base_power_product< ring_matrix < T >, lidia_size_t > &A, const T &mod) const;

	//
	// Hessenberg form
	//

	void HBF(MATRIX_TYPE &, const T &) const;

	//
	// characteristic polynomial
	//

	T *charpoly(MATRIX_TYPE &, const T &) const;
	T *charpoly(const base_power_product< ring_matrix < T >, lidia_size_t > &A, const T &mod) const;

	//
	// lanczos algorithm
	//

	bool lanczos(const MATRIX_TYPE &A, math_vector< T > &x, const math_vector< T > &b, T &mod) const;
	T lanczos_ZmZ(const MATRIX_TYPE &A, math_vector< T > &x, const math_vector< T > &b, T &mod) const;

	//
	// wiedemann algorithm
	//

	bool wiedemann(const ring_matrix< T > &A,
		       math_vector< T > &x,
		       const math_vector< T > &b,
		       T &mod) const;
	bool wiedemann(const base_power_product< ring_matrix < T >, lidia_size_t > &A,
		       math_vector< T > &x,
		       const math_vector< T > &b, T &mod) const;

	//
	// conjugate_gradient
	//

	bool conjugate_gradient(const ring_matrix< T > &A,
				math_vector< T > &x,
				const math_vector< T > &b,
				T &mod) const;

};



template< class T >
inline T * Fp_polynomial_convert(const Fp_polynomial &res, const T&)
{
	T *bigres = new T[res.degree() + 1];
	for (lidia_size_t i = 0; i <= res.degree(); i++)
		bigres[i] = res[i];
	return bigres;
}



inline long *Fp_polynomial_convert(const Fp_polynomial &res, const long &)
{
	long *bigres = new long[res.degree() + 1];
	for (lidia_size_t i = 0; i <= res.degree(); i++)
		res[i].longify(bigres[i]);
	return bigres;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_H_GUARD_
