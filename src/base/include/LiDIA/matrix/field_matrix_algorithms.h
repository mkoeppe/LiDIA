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


#ifndef LIDIA_FIELD_MATRIX_ALGORITHMS_H_GUARD_
#define LIDIA_FIELD_MATRIX_ALGORITHMS_H_GUARD_



#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
# include	"LiDIA/matrix/matrix_representation.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define field_matrix_algorithms FMA



template< class T , class ARG1, class ARG2, class ARG3 >
class field_matrix_algorithms
{

private:

	const ARG1 modul1;
	const ARG2 modul2;
	const ARG3 modul3;

public:

	//
	// constructor
	//

	field_matrix_algorithms() {}

	//
	// destructor
	//

	~field_matrix_algorithms() {}

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



#undef field_matrix_algorithms



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/field_matrix_algorithms.cc"
#endif



#endif	// LIDIA_FIELD_MATRIX_ALGORITHMS_H_GUARD_
