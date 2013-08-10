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
//	$Id: sparse_ring_matrix_kernel.cc,v 2.8 2004/06/15 10:19:20 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SPARSE_RING_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_SPARSE_RING_MATRIX_KERNEL_CC_GUARD_



#define sparse_ring_matrix_kernel SRMK


#ifndef LIDIA_ERROR_H_GUARD_
#include "LiDIA/error.h"
#endif

#ifndef LIDIA_SPARSE_RING_MATRIX_KERNEL_H_GUARD_
#include "LiDIA/matrix/sparse_ring_matrix_kernel.h"
#endif


#include	"LiDIA/matrix/sparse_base_matrix_kernel.cc"

#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif


//
// debug defines / error defines
//





extern const char *PRT;
extern const char *matrix_error_msg[];

#define DV_MM LDBL_MATRIX + 10        // Debug value
#define DM_MM "MATRIX_TYPE"           // Debug message / Error message
#define LMM_ERROR matrix_error_msg

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T >
void
sparse_ring_matrix_kernel< T >::add (MR< T > &RES, const MR< T > &M, const MR< T > &N) const
{
	register lidia_size_t i, j, l1, l2, l3, l;

	T TMP;
	T *Mtmp, *Ntmp, *REStmp;
	lidia_size_t *Mindex, *Nindex, *RESindex;
	lidia_size_t counter;

	for (i = 0; i < RES.rows; i++) {
		// Memory allocation
		l = M.allocated[i] + N.allocated[i];
		if (l > RES.columns)
			l = RES.columns;

		if (l == 0) {
			REStmp = NULL;
			RESindex = NULL;
		}
		else {
			REStmp = new T[l];
			RESindex = new lidia_size_t[l];
		}

		l1 = 0;
		l2 = 0;
		l3 = 0;
		Ntmp = N.value[i];
		Nindex = N.index[i];

		Mtmp = M.value[i];
		Mindex = M.index[i];

		RES.allocated[i] = l;
		counter = 0;

		for (j = 0; j < RES.columns; j++) {
			if (l1 < M.value_counter[i] && Mindex[l1] == j) {
				if (l2 < N.value_counter[i] && Nindex[l2] == j) {
					LiDIA::add(TMP, Mtmp[l1], Ntmp[l2]);
					if (TMP != RES.Zero) {
						REStmp[l3] = TMP;
						RESindex[l3] = j;
						l3++;
						counter++;
					}
					l2++;
					l1++;
				}
				else {
					REStmp[l3] = Mtmp[l1];
					RESindex[l3] = j;
					l3++;
					l1++;
					counter++;
				}
			}
			else {
				if (l2 < N.value_counter[i] && Nindex[l2] == j) {
					REStmp[l3] = Ntmp[l2];
					RESindex[l3] = j;
					l3++;
					l2++;
					counter++;
				}
			}
		}

		delete[] RES.value[i];
		delete[] RES.index[i];

		RES.value[i] = REStmp;
		RES.index[i] = RESindex;
		RES.value_counter[i] = counter;
	}
}



template< class T >
void
sparse_ring_matrix_kernel< T >::add (MR< T > &RES, const MR< T > &M, const T &a) const
{
	register lidia_size_t i, j, l1, l3, l;

	T TMP;
	T *Mtmp, *REStmp;
	lidia_size_t *Mindex, *RESindex;
	lidia_size_t counter;

	if (a == RES.Zero)
//Changed by G.A.	
		//assign(RES, M);
		this->assign(RES, M);
		//sparse_base_matrix_kernel< T >::assign(RES, M);

	else
		for (i = 0; i < RES.rows; i++) {
			// Memory allocation
			l = M.columns;
			if (l > RES.columns)
				l = RES.columns;

			if (l == 0) {
				REStmp = NULL;
				RESindex = NULL;
			}
			else {
				REStmp = new T[l];
				RESindex = new lidia_size_t[l];
			}

			l1 = l3 = 0;
			Mtmp = M.value[i];
			Mindex = M.index[i];

			RES.allocated[i] = l;
			counter = 0;

			for (j = 0; j < RES.columns; j++) {
				if (l1 < M.value_counter[i] && Mindex[l1] == j) {
					LiDIA::add(TMP, Mtmp[l1], a);
					if (TMP != RES.Zero) {
						REStmp[l3] = TMP;
						RESindex[l3] = j;
						l3++;
						counter++;
					}
					l1++;
				}
				else {
					REStmp[l3] = a;
					RESindex[l3] = j;
					l3++;
					counter++;
				}
			}

			delete[] RES.value[i];
			delete[] RES.index[i];

			RES.value[i] = REStmp;
			RES.index[i] = RESindex;
			RES.value_counter[i] = counter;
		}
}



template< class T >
void
sparse_ring_matrix_kernel< T >::add (MR< T > &RES, const T &a, const MR< T > &M) const
{

	register lidia_size_t i, j, l1, l3, l;

	T TMP;
	T *Mtmp, *REStmp;
	lidia_size_t *Mindex, *RESindex;
	lidia_size_t counter;

	if (a == RES.Zero)
//Changed by G.A.	
		//assign(RES, M);
		this->assign(RES, M);
		//sparse_base_matrix_kernel< T >::assign(RES, M);
	
	else
		for (i = 0; i < RES.rows; i++) {
			// Memory allocation
			l = M.columns;
			if (l > RES.columns)
				l = RES.columns;

			if (l == 0) {
				REStmp = NULL;
				RESindex = NULL;
			}
			else {
				REStmp = new T[l];
				RESindex = new lidia_size_t[l];
			}

			l1 = l3 = 0;
			Mtmp = M.value[i];
			Mindex = M.index[i];

			RES.allocated[i] = l;
			counter = 0;

			for (j = 0; j < RES.columns; j++) {
				if (l1 < M.value_counter[i] && Mindex[l1] == j) {
					LiDIA::add(TMP, a, Mtmp[l1]);
					if (TMP != RES.Zero) {
						REStmp[l3] = TMP;
						RESindex[l3] = j;
						l3++;
						counter++;
					}
					l1++;
				}
				else {
					REStmp[l3] = a;
					RESindex[l3] = j;
					l3++;
					counter++;
				}
			}

			delete[] RES.value[i];
			delete[] RES.index[i];

			RES.value[i] = REStmp;
			RES.index[i] = RESindex;
			RES.value_counter[i] = counter;
		}
}



//
// subtraction
//

template< class T >
void
sparse_ring_matrix_kernel< T >::subtract (MR< T > &RES, const MR< T > &M, const MR< T > &N) const
{
	register lidia_size_t i, j, l1, l2, l3, l;

	T TMP;
	T *Mtmp, *Ntmp, *REStmp;
	lidia_size_t *Mindex, *Nindex, *RESindex;
	lidia_size_t counter;

	for (i = 0; i < RES.rows; i++) {
		// Memory allocation
		l = M.allocated[i] + N.allocated[i];
		if (l > RES.columns)
			l = RES.columns;

		if (l == 0) {
			REStmp = NULL;
			RESindex = NULL;
		}
		else {
			REStmp = new T[l];
			RESindex = new lidia_size_t[l];
		}

		l1 = 0;
		l2 = 0;
		l3 = 0;
		Ntmp = N.value[i];
		Nindex = N.index[i];

		Mtmp = M.value[i];
		Mindex = M.index[i];

		RES.allocated[i] = l;
		counter = 0;

		for (j = 0; j < RES.columns; j++) {
			if (l1 < M.value_counter[i] && Mindex[l1] == j) {
				if (l2 < N.value_counter[i] && Nindex[l2] == j) {
					LiDIA::subtract(TMP, Mtmp[l1], Ntmp[l2]);
					if (TMP != RES.Zero) {
						REStmp[l3] = TMP;
						RESindex[l3] = j;
						l3++;
						counter++;
					}
					l2++;
					l1++;
				}
				else {
					REStmp[l3] = Mtmp[l1];
					RESindex[l3] = j;
					l3++;
					l1++;
					counter++;
				}
			}
			else {
				if (l2 < N.value_counter[i] && Nindex[l2] == j) {
					REStmp[l3] = -Ntmp[l2];
					RESindex[l3] = j;
					l3++;
					l2++;
					counter++;
				}
			}
		}

		delete[] RES.value[i];
		delete[] RES.index[i];

		RES.value[i] = REStmp;
		RES.index[i] = RESindex;
		RES.value_counter[i] = counter;
	}
}



template< class T >
void
sparse_ring_matrix_kernel< T >::subtract (MR< T > &RES, const MR< T > &M, const T &a) const
{
	register lidia_size_t i, j, l1, l3, l;

	T TMP;
	T *Mtmp, *REStmp;
	lidia_size_t *Mindex, *RESindex;
	lidia_size_t counter;

	if (a == RES.Zero)
//Changed by G.A.	
		//assign(RES, M);
		this->assign(RES, M);
		//sparse_base_matrix_kernel< T >::assign(RES, M);
	
	
	else
		for (i = 0; i < RES.rows; i++) {
			// Memory allocation
			l = M.columns;
			if (l > RES.columns)
				l = RES.columns;

			if (l == 0) {
				REStmp = NULL;
				RESindex = NULL;
			}
			else {
				REStmp = new T[l];
				RESindex = new lidia_size_t[l];
			}

			l1 = l3 = 0;
			Mtmp = M.value[i];
			Mindex = M.index[i];

			RES.allocated[i] = l;
			counter = 0;

			for (j = 0; j < RES.columns; j++) {
				if (l1 < M.value_counter[i] && Mindex[l1] == j) {
					LiDIA::subtract(TMP, Mtmp[l1], a);
					if (TMP != RES.Zero) {
						REStmp[l3] = TMP;
						RESindex[l3] = j;
						l3++;
						counter++;
					}
					l1++;
				}
				else {
					REStmp[l3] = -a;
					RESindex[l3] = j;
					l3++;
					counter++;
				}
			}

			delete[] RES.value[i];
			delete[] RES.index[i];

			RES.value[i] = REStmp;
			RES.index[i] = RESindex;
			RES.value_counter[i] = counter;
		}
}



template< class T >
void
sparse_ring_matrix_kernel< T >::subtract (MR< T > &RES, const T &a, const MR< T > &M) const
{
	register lidia_size_t i, j, l1, l3, l;

	T TMP;
	T *Mtmp, *REStmp;
	lidia_size_t *Mindex, *RESindex;
	lidia_size_t counter;

	if (a == RES.Zero)
//Changed by G.A.	
		//assign(RES, M);
		this->assign(RES, M);
		//sparse_base_matrix_kernel< T >::assign(RES, M);
	else
		for (i = 0; i < RES.rows; i++) {
			// Memory allocation
			l = M.columns;
			if (l > RES.columns)
				l = RES.columns;

			if (l == 0) {
				REStmp = NULL;
				RESindex = NULL;
			}
			else {
				REStmp = new T[l];
				RESindex = new lidia_size_t[l];
			}

			l1 = l3 = 0;
			Mtmp = M.value[i];
			Mindex = M.index[i];

			RES.allocated[i] = l;
			counter = 0;

			for (j = 0; j < RES.columns; j++) {
				if (l1 < M.value_counter[i] && Mindex[l1] == j) {
					LiDIA::subtract(TMP, a, Mtmp[l1]);
					if (TMP != RES.Zero) {
						REStmp[l3] = TMP;
						RESindex[l3] = j;
						l3++;
						counter++;
					}
					l1++;
				}
				else {
					REStmp[l3] = a;
					RESindex[l3] = j;
					l3++;
					counter++;
				}
			}

			delete[] RES.value[i];
			delete[] RES.index[i];

			RES.value[i] = REStmp;
			RES.index[i] = RESindex;
			RES.value_counter[i] = counter;
		}
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////



#undef DV_MM
#undef DM_MM
#undef LMM_ERROR


#undef sparse_ring_matrix_kernel



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_RING_MATRIX_KERNEL_CC_GUARD_
