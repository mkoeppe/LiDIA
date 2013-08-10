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


#ifndef LIDIA_SPARSE_FIELD_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_SPARSE_FIELD_MATRIX_KERNEL_CC_GUARD_



#ifndef LIDIA_ERROR_H_GUARD_
# include	"LiDIA/error.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_SPARSE_FIELD_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_field_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define sparse_field_matrix_kernel SFMK

//
// debug defines / error defines
//


extern const char *PRT;
extern const char *matrix_error_msg[];

#define DV_MM LDBL_MATRIX + 10   // Debug value
#define DM_MM "MATRIX_TYPE"      // Debug message / Error message
#define LMM_ERROR matrix_error_msg



//
// divide
//

template< class T >
inline void
sparse_field_matrix_kernel< T >::divide (MR< T > &RES,
					 const MR< T > &A,
					 const T &k) const
{
	debug_handler_l(DM_MM, "in member - function "
			"divide(const MR< T > &, const T &)", DV_MM + 2);

	register lidia_size_t j, i;
	T *REStmp, *Atmp;

	for (j = 0; j < A.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		for (i = 0; i < A.columns; i++)
			LiDIA::divide(REStmp[i], Atmp[i], k);
	}
}



template< class T >
inline void
sparse_field_matrix_kernel< T >::compwise_divide (MR< T > &RES,
						  const MR< T > &A,
						  const MR< T > &B) const
{
	debug_handler_l(DM_MM, "in member - function "
			"compwise_divide(const MR< T > &, const MR< T > &)", DV_MM + 3);

	register lidia_size_t j, i;
	T *REStmp, *Atmp, *Btmp;

	for (j = 0; j < RES.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		Btmp = B.value[j];
		for (i = 0; i < RES.columns; i++)
			LiDIA::divide(REStmp[i], Atmp[i], Btmp[i]);
	}
}




#undef DV_MM
#undef DM_MM
#undef LMM_ERROR



#undef sparse_field_matrix_kernel



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_FIELD_MATRIX_KERNEL_CC_GUARD_
