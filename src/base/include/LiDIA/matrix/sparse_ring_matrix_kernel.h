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


#ifndef LIDIA_SPARSE_RING_MATRIX_KERNEL_H_GUARD_
#define LIDIA_SPARSE_RING_MATRIX_KERNEL_H_GUARD_



#ifndef LIDIA_SPARSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define sparse_ring_matrix_kernel SRMK



template< class T >
class sparse_ring_matrix_kernel : public SBMK< T >
{

public:

	//
	// constructor
	//

	sparse_ring_matrix_kernel () {}

	//
	// destructor
	//

	~sparse_ring_matrix_kernel () {}

	using SBMK< T >::assign;

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// addition
	//

	void add(MR< T > &RES, const MR< T > &M, const MR< T > &N) const;
	void add(MR< T > &RES, const MR< T > &M, const T &a) const;
	void add(MR< T > &RES, const T &a, const MR< T > &M) const;

	//
	// subtraction
	//

	void subtract(MR< T > &RES, const MR< T > &M, const MR< T > &N) const;
	void subtract(MR< T > &RES, const MR< T > &M, const T &a) const;
	void subtract(MR< T > &RES, const T &a, const MR< T > &M) const;

	//////////////////////////////////
	// END: arithmetical procedures //
	//////////////////////////////////
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/sparse_ring_matrix_kernel.cc"
#endif


#undef sparse_ring_matrix_kernel



#endif	// LIDIA_SPARSE_RING_MATRIX_KERNEL_H_GUARD_
