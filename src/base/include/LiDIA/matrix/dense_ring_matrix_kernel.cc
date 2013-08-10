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


#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_DENSE_RING_MATRIX_KERNEL_CC_GUARD_



#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define dense_ring_matrix_kernel DRMK

//
// debug defines / error defines
//


#define DV_MM LDBL_MATRIX + 10   // Debug value
#define DM_MM "MATRIX_TYPE"      // Debug message / Error message

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T >
inline void dense_ring_matrix_kernel< T >::
add(MR< T > &RES, const MR< T > &M, const MR< T > &N) const
{
	register lidia_size_t i, j;

	T *Mtmp, *Ntmp, *REStmp;
	for (i = 0; i < RES.rows; i++) {
		Ntmp = N.value[i];
		Mtmp = M.value[i];
		REStmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			LiDIA::add(REStmp[j], Ntmp[j], Mtmp[j]);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
add(MR< T > &RES, const MR< T > &M, const T &a) const
{
	register lidia_size_t i, j;
	T *REStmp, *Mtmp;

	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < RES.columns; j++)
			LiDIA::add(REStmp[j], Mtmp[j], a);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
add(MR< T > &RES, const T &a, const MR< T > &M) const
{
	register lidia_size_t i, j;
	T *REStmp, *Mtmp;

	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < RES.columns; j++)
			LiDIA::add(REStmp[j], Mtmp[j], a);
	}
}



//
// subtraction
//

template< class T >
inline void dense_ring_matrix_kernel< T >::
subtract(MR< T > &RES, const MR< T > &M, const MR< T > &N) const
{
	register lidia_size_t i, j;

	T *REStmp, *Mtmp, *Ntmp;

	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		Ntmp = N.value[i];
		for (j = 0; j < RES.columns; j++)
			LiDIA::subtract(REStmp[j], Mtmp[j], Ntmp[j]);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
subtract(MR< T > &RES, const MR< T > &M, const T &a) const
{
	register lidia_size_t i, j;
	T *REStmp, *Mtmp;

	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < RES.columns; j++)
			LiDIA::subtract(REStmp[j], Mtmp[j], a);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
subtract(MR< T > &RES, const T &a, const MR< T > &M) const
{
	register lidia_size_t i, j;
	T *REStmp, *Mtmp;

	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < RES.columns; j++)
			LiDIA::subtract(REStmp[j], a, Mtmp[j]);
	}
}



//
// multiplication
//

template< class T >
inline void dense_ring_matrix_kernel< T >::
multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B) const
{
	register lidia_size_t j, i, z;
	T TMP, TMP1, *Atmp, *REStmp;

	for (j = 0; j < A.rows; j++) {
		Atmp = A.value[j];
		REStmp = RES.value[j];
		for (i = 0; i < B.columns; i++) {
			TMP = 0; // TMP.assign_zero();
			for (z = 0; z < B.rows; z++) {
				LiDIA::multiply(TMP1, Atmp[z], B.value[z][i]);
				LiDIA::add(TMP, TMP, TMP1);
			}
			REStmp[i] = TMP;
		}
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
multiply(MR< T > &RES, const MR< T > &A, const T &k) const
{
	register lidia_size_t j, i;
	T *REStmp, *Atmp;

	for (j = 0; j < A.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		for (i = 0; i < A.columns; i++)
			LiDIA::multiply(REStmp[i], Atmp[i], k);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
multiply(MR< T > &RES, const T &k, const MR< T > &A) const
{
	register lidia_size_t j, i;
	T *REStmp, *Atmp;
	for (j = 0; j < A.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		for (i = 0; i < A.columns; i++)
			LiDIA::multiply(REStmp[i], Atmp[i], k);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
compwise_multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B) const
{
	register lidia_size_t j, i;
	T *REStmp, *Atmp, *Btmp;

	for (j = 0; j < RES.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		Btmp = B.value[j];
		for (i = 0; i < RES.columns; i++)
			LiDIA::multiply(REStmp[i], Atmp[i], Btmp[i]);
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
multiply_right(const MR< T > &RES, T *&c, const T *v) const
{
	register lidia_size_t i, j;
	T TMP, TMP1, *tmp;

	if (c == v || c == NULL) {
		c = new T[RES.rows];
		memory_handler(c, DM_MM, "multiply_right :: "
			       "Error in memory allocation (c)");
	}

	for (i = 0; i < RES.rows; i++) {
		TMP = 0;
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++) {
			LiDIA::multiply(TMP1, tmp[j], v[j]);
			LiDIA::add(TMP, TMP, TMP1);
		}
		c[i] = TMP;
	}
}



template< class T >
inline void dense_ring_matrix_kernel< T >::
multiply_left(const MR< T > &RES, T *&c, const T *v) const
{
	register lidia_size_t i, j;
	T TMP, TMP1;

	if (c == v || c == NULL) {
		c = new T[RES.columns];
		memory_handler(c, DM_MM, "multiply_left :: "
			       "Error in memory allocation (c)");
	}

	for (i = 0; i < RES.columns; i++) {
		TMP = 0;
		for (j = 0; j < RES.rows; j++) {
			LiDIA::multiply(TMP1, RES.value[j][i], v[j]);
			LiDIA::add(TMP, TMP, TMP1);
		}
		c[i] = TMP;
	}
}



//
// negation
//

template< class T >
inline void dense_ring_matrix_kernel< T >::
negate(MR< T > &RES, const MR< T > &B) const
{
	register lidia_size_t i, j;

	T *Btmp, *REStmp;

	for (i = 0; i < B.rows; i++) {
		Btmp = B.value[i];
		REStmp = RES.value[i];
		for (j = 0; j < B.columns; j++)
			LiDIA::negate(REStmp[j], Btmp[j]);
	}
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////

//
// comparisons
//

template< class T >
inline bool dense_ring_matrix_kernel< T >::
equal(const MR< T > &RES, const MR< T > &N) const
{
	register lidia_size_t i, j;
	if (RES.rows != N.rows || RES.columns != N.columns)
		return false;

	T *tmp, *Ntmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		Ntmp = N.value[i];
		for (j = 0; j < RES.columns; j++)
			if (tmp[j] != Ntmp[j])
				return false;
	}
	return true;
}



//
// trace
//

template< class T >
inline void dense_ring_matrix_kernel< T >::
trace(const MR< T > &RES, T &tr) const
{
	register lidia_size_t i;
	tr = RES.value[0][0];
	for (i = 1; i < RES.rows; i++)
		LiDIA::add(tr, tr, RES.value[i][i]);
}



#undef DV_MM
#undef DM_MM

#undef dense_ring_matrix_kernel



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_DENSE_RING_MATRIX_KERNEL_CC_GUARD_
