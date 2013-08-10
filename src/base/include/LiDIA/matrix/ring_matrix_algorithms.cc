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


#ifndef LIDIA_RING_MATRIX_ALGORITHMS_CC_GUARD_
#define LIDIA_RING_MATRIX_ALGORITHMS_CC_GUARD_



#ifndef LIDIA_ERROR_H_GUARD_
# include	"LiDIA/error.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#include	"LiDIA/modular_operations.inl"



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define ring_matrix_algorithms RMA

//
// debug defines / error defines
//

#define DM_MM "ring_matrix_algorithms"  // Debug message / Error message

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
add(MR< T > &RES, const MR< T > &M, const MR< T > &N) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < N.rows; i++) {
		for (j = 0; j < N.columns; j++) {
			LiDIA::add(TMP, modul2.member(M, i , j), modul3.member(N, i, j));
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
add(MR< T > &RES, const MR< T > &M, const T &a) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::add(TMP, modul2.member(M, i, j), a);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
add(MR< T > &RES, const T &a, const MR< T > &M) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::add(TMP, modul2.member(M, i, j), a);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



//
// subtraction
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
subtract(MR< T > &RES, const MR< T > &M, const MR< T > &N) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::subtract(TMP, modul2.member(M, i, j), modul3.member(N, i, j));
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
subtract(MR< T > &RES, const MR< T > &M, const T &a) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::subtract(TMP, modul2.member(M, i, j), a);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
subtract(MR< T > &RES, const T &a, const MR< T > &M) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::subtract(TMP, a, modul2.member(M, i, j));
			modul1.sto(RES, i, j, TMP);
		}
	}
}



//
// multiplication
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B) const
{
	register lidia_size_t j, i, z;
	T TMP, TMP1;

	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < B.columns; i++) {
			TMP = 0; // TMP.assign_zero();
			for (z = 0; z < B.rows; z++) {
				LiDIA::multiply(TMP1, modul2.member(A, j, z), modul3.member(B, z, i));
				LiDIA::add(TMP, TMP, TMP1);
			}
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply(MR< T > &RES, const MR< T > &A, const T &k) const
{
	register lidia_size_t j, i;
	T TMP;

	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < A.columns; i++) {
			LiDIA::multiply(TMP, modul2.member(A, j, i), k);
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply(MR< T > &RES, const T &k, const MR< T > &A) const
{
	register lidia_size_t j, i;
	T TMP;
	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < A.columns; i++) {
			LiDIA::multiply(TMP, modul2.member(A, j, i), k);
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
compwise_multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B) const
{
	register lidia_size_t j, i;
	T TMP;

	for (j = 0; j < RES.rows; j++) {
		for (i = 0; i < RES.columns; i++) {
			LiDIA::multiply(TMP, modul2.member(A, j, i), modul3.member(B, j, i));
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_right(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v) const
{
	register lidia_size_t i, j, l = res.size();
	T TMP, TMP1;

	T *tmp1 = res.get_data_address(), *tmp2 = v.get_data_address();

	if (tmp1 == tmp2) {
		math_vector< T > res2(RES.rows, RES.rows);
		tmp1 = res2.get_data_address();

		for (i = 0; i < RES.rows; i++) {
			TMP = 0;
			for (j = 0; j < RES.columns; j++) {
				LiDIA::multiply(TMP1, modul1.member(RES, i, j), tmp2[j]);
				LiDIA::add(TMP, TMP, TMP1);
			}
			tmp1[i] = TMP;
		}
		res = res2;
	}
	else {
		if (l != RES.rows) {
			if (res.capacity() < RES.rows)
				res.set_capacity(RES.rows);
			if (l != RES.rows)
				res.set_size(RES.rows);
		}
		tmp1 = res.get_data_address();

		for (i = 0; i < RES.rows; i++) {
			TMP = 0;
			for (j = 0; j < RES.columns; j++) {
				LiDIA::multiply(TMP1, modul1.member(RES, i, j), tmp2[j]);
				LiDIA::add(TMP, TMP, TMP1);
			}
			tmp1[i] = TMP;
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_right(const MR< T > &RES, T *&c, const T *v) const
{
	register lidia_size_t i, j;
	T TMP, TMP1;

	if (c == v || c == NULL) {
		c = new T[RES.rows];
		memory_handler(c, DM_MM, "multiply_right :: "
			       "Error in memory allocation (c)");
	}

	for (i = 0; i < RES.rows; i++) {
		TMP = 0;
		for (j = 0; j < RES.columns; j++) {
			LiDIA::multiply(TMP1, modul1.member(RES, i, j), v[j]);
			LiDIA::add(TMP, TMP, TMP1);
		}
		c[i] = TMP;
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_left(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v) const
{
	register lidia_size_t i, j, l = res.size();
	T TMP, TMP1;

	T *tmp1 = res.get_data_address(), *tmp2 = v.get_data_address();

	if (tmp1 == tmp2) {
		math_vector< T > res2(RES.columns, RES.columns);
		tmp1 = res2.get_data_address();
		for (i = 0; i < RES.columns; i++) {
			TMP = 0;
			for (j = 0; j < RES.rows; j++) {
				LiDIA::multiply(TMP1, modul1.member(RES, j, i), tmp2[j]);
				LiDIA::add(TMP, TMP, TMP1);
			}
			tmp1[i] = TMP;
		}
		res.assign(res2);
	}
	else {
		if (l != RES.columns) {
			if (res.capacity() < RES.columns)
				res.set_capacity(RES.columns);
			if (l != RES.columns)
				res.set_size(RES.columns);
		}
		tmp1 = res.get_data_address();

		for (i = 0; i < RES.columns; i++) {
			TMP = 0;
			for (j = 0; j < RES.rows; j++) {
				LiDIA::multiply(TMP1, modul1.member(RES, j, i), tmp2[j]);
				LiDIA::add(TMP, TMP, TMP1);
			}
			tmp1[i] = TMP;
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
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
			LiDIA::multiply(TMP1, modul1.member(RES, j, i), v[j]);
			LiDIA::add(TMP, TMP, TMP1);
		}
		c[i] = TMP;
	}
}



//
// negation
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
negate(MR< T > &RES, const MR< T > &B) const
{
	register lidia_size_t i, j;

	T TMP;

	for (i = 0; i < B.rows; i++) {
		for (j = 0; j < B.columns; j++) {
			LiDIA::negate(TMP, modul2.member(B, i, j));
			modul1.sto(RES, i, j, TMP);
		}
	}
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////

//
// comparisons
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline bool ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
equal(const MR< T > &RES, const MR< T > &N) const
{
	register lidia_size_t i, j;
	if (RES.rows != N.rows || RES.columns != N.columns)
		return false;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++)
			if (modul1.member(RES, i, j) != modul2.member(N, i, j))
				return false;
	}
	return true;
}



//
// trace
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
trace(const MR< T > &RES, T &tr) const
{
	register lidia_size_t i;
	tr = modul1.member(RES, 0, 0);
	for (i = 1; i < RES.rows; i++)
		LiDIA::add(tr, tr, modul1.member(RES, i, i));
}



//
// size reduction
//

template< class T, class ARG1, class ARG2, class ARG3 >
void ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
size_reduction(MR< T > &RES) const
{
	lidia_size_t *no_of_elements = new lidia_size_t[RES.rows];
	memory_handler(no_of_elements, DM_MM, "size_reduction :: "
		       "Error in memory allocation (no_of_elements)");

	lidia_size_t i, j;
	for (i = 0; i < RES.rows; i++) {
		no_of_elements[i] = 0;
		for (j = 0; j < RES.columns; j++)
			if (modul1.member(RES, i, j) != RES.Zero)
				no_of_elements[i]++;
	}

	// main loop
	lidia_size_t rpos = RES.rows - 1;
	lidia_size_t cpos = RES.columns - 1;

	for (i = RES.rows - 1; i >= 0; i--) {
		//      std::cout << " POSION: " << i << " =  >" << no_of_elements[i] << std::endl;
		if (no_of_elements[i] == 1) {
			modul1.swap_rows(RES, i, rpos);
			LiDIA::swap(no_of_elements[i], no_of_elements[rpos]);
			for (j = 0; j <= cpos; j++)
				if (modul1.member(RES, rpos, j) != RES.Zero)
					modul1.swap_columns(RES, j, cpos);

			// update of array no_of_elements
			for (j = 0; j < rpos; j++)
				if (modul1.member(RES, j, cpos) != RES.Zero)
					no_of_elements[j]--;
			i = rpos;
			rpos--;
			cpos--;
			//std::cout << "rpos = " << rpos << " cpos = " << cpos << std::endl;
		}
	}
}



#define modular_ring_matrix_algorithms MRMA

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
add(MR< T > &RES, const MR< T > &M, const MR< T > &N, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < N.rows; i++) {
		for (j = 0; j < N.columns; j++) {
			LiDIA::add_mod(TMP, modul2.member(M, i , j), modul3.member(N, i, j), mod);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
add(MR< T > &RES, const MR< T > &M, const T &a, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::add_mod(TMP, modul2.member(M, i, j), a, mod);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
add(MR< T > &RES, const T &a, const MR< T > &M, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::add_mod(TMP, modul2.member(M, i, j), a, mod);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



//
// subtraction
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
subtract(MR< T > &RES, const MR< T > &M, const MR< T > &N, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::sub_mod(TMP, modul2.member(M, i, j), modul3.member(N, i, j), mod);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
subtract(MR< T > &RES, const MR< T > &M, const T &a, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::sub_mod(TMP, modul2.member(M, i, j), a, mod);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
subtract(MR< T > &RES, const T &a, const MR< T > &M, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			LiDIA::sub_mod(TMP, a, modul2.member(M, i, j), mod);
			modul1.sto(RES, i, j, TMP);
		}
	}
}



//
// multiplication
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B, const T &mod) const
{
	register lidia_size_t j, i, z;
	T TMP, TMP1;

	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < B.columns; i++) {
			TMP = 0; // TMP.assign_zero();
			for (z = 0; z < B.rows; z++) {
				LiDIA::mult_mod(TMP1, modul2.member(A, j, z), modul3.member(B, z, i), mod);
				LiDIA::add_mod(TMP, TMP, TMP1, mod);
			}
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply(MR< T > &RES, const MR< T > &A, const T &k, const T &mod) const
{
	register lidia_size_t j, i;
	T TMP;

	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < A.columns; i++) {
			LiDIA::mult_mod(TMP, modul2.member(A, j, i), k, mod);
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply(MR< T > &RES, const T &k, const MR< T > &A, const T &mod) const
{
	register lidia_size_t j, i;
	T TMP;
	for (j = 0; j < A.rows; j++) {
		for (i = 0; i < A.columns; i++) {
			LiDIA::mult_mod(TMP, modul2.member(A, j, i), k, mod);
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
compwise_multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B, const T &mod) const
{
	register lidia_size_t j, i;
	T TMP;

	for (j = 0; j < RES.rows; j++) {
		for (i = 0; i < RES.columns; i++) {
			LiDIA::mult_mod(TMP, modul2.member(A, j, i), modul3.member(B, j, i), mod);
			modul1.sto(RES, j, i, TMP);
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_right(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v, const T &mod) const
{
	register lidia_size_t i, j, l = res.size();
	T TMP, TMP1;

	T *tmp1 = res.get_data_address(), *tmp2 = v.get_data_address();

	if (tmp1 == tmp2) {
		math_vector< T > res2(RES.rows, RES.rows);
		tmp1 = res2.get_data_address();

		for (i = 0; i < RES.rows; i++) {
			TMP = 0;
			for (j = 0; j < RES.columns; j++) {
				LiDIA::mult_mod(TMP1, modul1.member(RES, i, j), tmp2[j], mod);
				LiDIA::add_mod(TMP, TMP, TMP1, mod);
			}
			tmp1[i] = TMP;
		}
		res = res2;
	}
	else {
		if (l != RES.rows) {
			if (res.capacity() < RES.rows)
				res.set_capacity(RES.rows);
			if (l != RES.rows)
				res.set_size(RES.rows);
		}
		tmp1 = res.get_data_address();

		for (i = 0; i < RES.rows; i++) {
			TMP = 0;
			for (j = 0; j < RES.columns; j++) {
				LiDIA::mult_mod(TMP1, modul1.member(RES, i, j), tmp2[j], mod);
				LiDIA::add_mod(TMP, TMP, TMP1, mod);
			}
			tmp1[i] = TMP;
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_right(const MR< T > &RES, T *&c, const T *v, const T &mod) const
{
	register lidia_size_t i, j;
	T TMP, TMP1;

	if (c == v || c == NULL) {
		c = new T[RES.rows];
		memory_handler(c, DM_MM, "multiply_right :: "
			       "Error in memory allocation (c)");
	}

	for (i = 0; i < RES.rows; i++) {
		TMP = 0;
		for (j = 0; j < RES.columns; j++) {
			LiDIA::mult_mod(TMP1, modul1.member(RES, i, j), v[j], mod);
			LiDIA::add_mod(TMP, TMP, TMP1, mod);
		}
		c[i] = TMP;
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_left(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v, const T &mod) const
{
	register lidia_size_t i, j, l = res.size();
	T TMP, TMP1;

	T *tmp1 = res.get_data_address(), *tmp2 = v.get_data_address();

	if (tmp1 == tmp2) {
		math_vector< T > res2(RES.columns, RES.columns);
		tmp1 = res2.get_data_address();
		for (i = 0; i < RES.columns; i++) {
			TMP = 0;
			for (j = 0; j < RES.rows; j++) {
				LiDIA::mult_mod(TMP1, modul1.member(RES, j, i), tmp2[j], mod);
				LiDIA::add_mod(TMP, TMP, TMP1, mod);
			}
			tmp1[i] = TMP;
		}
		res.assign(res2);
	}
	else {
		if (l != RES.columns) {
			if (res.capacity() < RES.columns)
				res.set_capacity(RES.columns);
			if (l != RES.columns)
				res.set_size(RES.columns);
		}
		tmp1 = res.get_data_address();

		for (i = 0; i < RES.columns; i++) {
			TMP = 0;
			for (j = 0; j < RES.rows; j++) {
				LiDIA::mult_mod(TMP1, modul1.member(RES, j, i), tmp2[j], mod);
				LiDIA::add_mod(TMP, TMP, TMP1, mod);
			}
			tmp1[i] = TMP;
		}
	}
}



template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
multiply_left(const MR< T > &RES, T *&c, const T *v, const T &mod) const
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
			LiDIA::mult_mod(TMP1, modul1.member(RES, j, i), v[j], mod);
			LiDIA::add_mod(TMP, TMP, TMP1, mod);
		}
		c[i] = TMP;
	}
}



//
// negation
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
negate(MR< T > &RES, const MR< T > &B, const T &mod) const
{
	register lidia_size_t i, j;

	T TMP;

	for (i = 0; i < B.rows; i++) {
		for (j = 0; j < B.columns; j++) {
			LiDIA::negate(TMP, modul2.member(B, i, j));
			modul1.sto(RES, i, j, TMP);
		}
	}
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////

//
// trace
//

template< class T, class ARG1, class ARG2, class ARG3 >
inline void modular_ring_matrix_algorithms< T, ARG1, ARG2, ARG3 >::
trace(const MR< T > &RES, T &tr, const T &mod) const
{
	register lidia_size_t i;
	tr = modul1.member(RES, 0, 0);
	for (i = 1; i < RES.rows; i++)
		LiDIA::add_mod(tr, tr, modul1.member(RES, i, i));
}



#undef DM_MM

#undef modular_ring_matrix_algorithms
#undef ring_matrix_algorithms



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_RING_MATRIX_ALGORITHMS_CC_GUARD_
