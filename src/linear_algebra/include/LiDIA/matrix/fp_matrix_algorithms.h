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
//	$Id: fp_matrix_algorithms.h,v 2.3 2002/06/24 09:42:04 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FP_MATRIX_ALGORITHMS_H_GUARD_
#define LIDIA_FP_MATRIX_ALGORITHMS_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T, class MATRIX_TYPE >
class fp_matrix_algorithms
{

public:

	//
	// constructor
	//

	fp_matrix_algorithms() {}

	//
	// destructor
	//

	~fp_matrix_algorithms() {}

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



#endif	// LIDIA_FP_MATRIX_ALGORITHMS_H_GUARD_
