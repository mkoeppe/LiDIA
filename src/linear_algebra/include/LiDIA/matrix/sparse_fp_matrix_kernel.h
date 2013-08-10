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


#ifndef LIDIA_SPARSE_FP_MATRIX_KERNEL_H_GUARD_
#define LIDIA_SPARSE_FP_MATRIX_KERNEL_H_GUARD_


#ifndef LIDIA_SPARSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T, class MATRIX_TYPE >
class sparse_fp_matrix_kernel
{

	const SBMK< T > S_base_modul;

public:

	//
	// constructor
	//

	sparse_fp_matrix_kernel() {}

	//
	// destructor
	//

	~sparse_fp_matrix_kernel() {}

	//
	// column step form
	//

	int STF(MATRIX_TYPE &, const T &) const;

	const T * STF_extended(MATRIX_TYPE &, const T &) const;

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

	//
	// Hessenberg form
	//

	void HBF(MATRIX_TYPE &, const T &) const;

	//
	// characteristic polynomial
	//

	T *charpoly(MATRIX_TYPE &, const T &) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SPARSE_FP_MATRIX_KERNEL_H_GUARD_
