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
//	$Id: base_matrix_algorithms.cc,v 2.4 2002/06/24 09:41:37 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_MATRIX_ALGORITHMS_CC_GUARD_
#define LIDIA_BASE_MATRIX_ALGORITHMS_CC_GUARD_



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define base_matrix_algorithms BMA

//
// debug defines / error defines
//

extern const char *PRT;
extern const char *matrix_error_msg[];

#define DVALUE LDBL_MATRIX                 // Debug value
#define DMESSAGE "base_matrix_algorithms"  // Debug messages
#define EMESSAGE matrix_error_msg          // Error messages



//
// constructors
//

template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::constructor (MR< T > &A) const
{
	A.rows = A.columns = 0;
	A.sparse_rows = A.sparse_columns = 0;
	A.value = NULL;
	A.index = NULL;
	A.allocated = NULL;
	A.value_counter = NULL;
}



template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::constructor (MR< T > &A,
						     const base_vector< T > &v) const
{
	modul.constructor(A, v.get_size(), 1);

	for (register lidia_size_t i = 0; i < A.rows; i++)
		modul.sto(A, i, 0, v[i]);
}



template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::constructor (MR< T > &A,
						     const MR< T > &M) const
{
	modul.constructor(A, M.rows, M.columns);

	for (register lidia_size_t i = 0; i < A.rows; i++)
		for (register lidia_size_t j = 0; j < A.columns; j++)
			modul.sto(A, i, j, modul1.member(M, i, j));
}



//
// exchange functions / swap functions
//

template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::swap_columns (MR< T > &A,
						      lidia_size_t i,
						      lidia_size_t j) const
{
	T TMP1, TMP2;


	for (register lidia_size_t k = A.rows - 1; k >= 0; k--) {
		TMP1 = modul.member(A, k, j);
		TMP2 = modul.member(A, k, i);
		modul.sto(A, k, j, TMP2);
		modul.sto(A, k, i, TMP1);
	}
}



template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::swap_rows (MR< T > &A,
						   lidia_size_t i,
						   lidia_size_t j) const
{
	T TMP;

	for (register lidia_size_t k = A.columns - 1; k >= 0; k--) {
		TMP = modul.member(A, j, k);
		modul.sto(A, j, k, modul.member(A, i, k));
		modul.sto(A, i, k, TMP);
	}
}



//
// insert_at / copy_from
//

template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::insert_at (MR< T > &M,
						   lidia_size_t r,
						   lidia_size_t c,
						   const MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t len1,
						   lidia_size_t len2) const
{
	lidia_size_t i, j;

	for (i = r; i < r + len1 && i < M.rows; i++) {
		for (j = c; j < c + len2 && j < M.columns; j++)
			modul.sto(M, i, j, modul1.member(A, startr+i-r, startc+j-c));
	}
}



//
// structur functions
//

template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::resize (MR< T > &A,
						lidia_size_t r,
						lidia_size_t c) const
{
	if (A.columns != c)
		modul.set_no_of_columns(A, c);
	if (A.rows != r)
		modul.set_no_of_rows(A, r);
}



template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::kill (MR< T > &A) const
{
	modul.set_no_of_columns(A, 0);
	modul.set_no_of_rows(A, 0);
}



//
// assignments
//

template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::assign (MR< T > &A,
						const MR< T > &M) const
{
	register lidia_size_t i, j;
	for (i = 0; i < A.rows; i++)
		for (j = 0; j < A.columns; j++)
			modul.sto(A, i, j, modul1.member(M, i, j));

	A.bitfield = M.bitfield;
}



//
// diagonal function
//

template< class T, class ARG, class ARG1 >
inline void
base_matrix_algorithms< T, ARG, ARG1 >::diag (MR< T > &A,
					      const T &a,
					      const T &b) const
{
	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++)
			modul.sto(A, i, j, ((i == j) ? a : b));
	}
}



//
// transpose function
//

template< class T, class ARG, class ARG1 >
void
base_matrix_algorithms< T, ARG, ARG1 >::trans (MR< T > &R,
					       const MR< T > &A) const
{
	register lidia_size_t i, j;

	if (A.value != R.value) {
		if (R.columns != A.rows)
			modul.set_no_of_columns(R, A.rows);
		if (R.rows != A.columns)
			modul.set_no_of_rows(R, A.columns);

		for (i = 0; i < A.rows; i++) {
			for (j = 0; j < A.columns; j++)
				modul.sto(R, j, i, modul1.member(A, i, j));
		}
	}
	else {
		T TMP;
		register lidia_size_t oldrows = R.rows, oldcolumns = R.columns;
		if (R.rows != R.columns)
			if (R.columns > R.rows) {
				modul.set_no_of_rows(R, oldcolumns);
				for (i = 0; i < R.rows; i++)
					for (j = 0; j < i; j++) {
						TMP = modul.member(R, i, j);
						modul.sto(R, i, j, modul.member(R, j, i));
						modul.sto(R, j, i, TMP);
					}
				modul.set_no_of_columns(R, oldrows);
			}
			else {
				modul.set_no_of_columns(R, oldrows);
				for (i = 0; i < R.rows; i++)
					for (j = 0; j < i; j++) {
						TMP = modul.member(R, i, j);
						modul.sto(R, i, j, modul.member(R, j, i));
						modul.sto(R, j, i, TMP);
					}
				modul.set_no_of_rows(R, oldcolumns);
			}
		else {
			for (i = 0; i < R.rows; i++)
				for (j = 0; j < i; j++) {
					TMP = modul.member(R, i, j);
					modul.sto(R, i, j, modul.member(R, j, i));
					modul.sto(R, j, i, TMP);
				}
		}
	}
}



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#undef base_matrix_algorithms



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_MATRIX_ALGORITHMS_CC_GUARD_
