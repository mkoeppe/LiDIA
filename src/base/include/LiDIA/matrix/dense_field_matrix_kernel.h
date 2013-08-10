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


#ifndef LIDIA_DENSE_FIELD_MATRIX_KERNEL_H_GUARD_
#define LIDIA_DENSE_FIELD_MATRIX_KERNEL_H_GUARD_



#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define dense_field_matrix_kernel DFMK



template< class T >
class dense_field_matrix_kernel : public DRMK< T >
{

public:

	//
	// constructor
	//

	dense_field_matrix_kernel () {}

	//
	// destructor
	//

	~dense_field_matrix_kernel () {}

	//
	// divide
	//

	void divide(MR< T > &RES, const MR< T > &A, const T &k) const;

	void compwise_divide(MR< T > &RES, const MR< T > &A, const MR< T > &B) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/dense_field_matrix_kernel.cc"
#endif



#undef dense_field_matrix_kernel



#endif	// LIDIA_DENSE_FIELD_MATRIX_KERNEL_H_GUARD_
