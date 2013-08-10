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


#ifndef LIDIA_FIELD_MATRIX_ALGORITHMS_CC_GUARD_
#define LIDIA_FIELD_MATRIX_ALGORITHMS_CC_GUARD_



#ifndef LIDIA_ERROR_H_GUARD_
# include	"LiDIA/error.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define field_matrix_algorithms FMA

//
// debug defines / error defines
//



#define DVALUE LDBL_MATRIX + 10   // Debug value
#define DMESSAGE "field_matrix_algorithms"      // Debug message / Error message

//
// divide
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void
field_matrix_algorithms< T, ARG1, ARG2, ARG3 >::divide (MR< T > &RES,
							const MR< T > &A,
							const T &k) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"divide(const MATRIX_TYPE &, const T &)", DVALUE + 2);

	register lidia_size_t j, i;
	T TMP;

	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < A.columns; i++) {
			LiDIA::divide(TMP, modul2.member(A, j, i), k);
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void
field_matrix_algorithms< T, ARG1, ARG2, ARG3 >::compwise_divide(MR< T > &RES,
								const MR< T > &A,
								const MR< T > &B) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"compwise_divide(const MATRIX_TYPE &, const MATRIX_TYPE &)", DVALUE + 3);

	register lidia_size_t j, i;
	T TMP;

	for (j = 0; j < RES.rows; j++) {
		for (i = 0; i < RES.columns; i++) {
			LiDIA::divide(TMP, modul2.member(A, j, i), modul3.member(B, j, i));
			modul1.sto(RES, j, i, TMP);
		}
	}
}



#undef DVALUE
#undef DMESSAGE


#undef field_matrix_algorithms



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_FIELD_MATRIX_ALGORITHMS_CC_GUARD_
