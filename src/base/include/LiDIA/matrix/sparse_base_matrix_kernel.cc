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


#ifndef LIDIA_SPARSE_BASE_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_SPARSE_BASE_MATRIX_KERNEL_CC_GUARD_



#ifndef LIDIA_SPARSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define sparse_base_matrix_kernel SBMK

//
// debug defines / error defines
//


extern const char *PRT;
extern const char *matrix_error_msg[];

#define DMESSAGE "sparse_base_matrix_kernel"  // Debug messages
#define EMESSAGE matrix_error_msg             // Error messages

//
// Expansion factor
//

#define DELTA 2

//
// constructor kernel
//

template< class T >
void
sparse_base_matrix_kernel< T >::constructor (MR< T > &A,
					     lidia_size_t r,
					     lidia_size_t c) const
{
	A.rows = r;
	A.columns = c;
	A.sparse_rows = A.sparse_columns = 0;

	A.value = new T *[r];
	memory_handler(A.value, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (A.value)");

	A.index = new lidia_size_t *[r];
	memory_handler(A.index, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (A.index)");

	A.value_counter = new lidia_size_t[r];
	memory_handler(A.value_counter, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (A.value_counter)");

	A.allocated = new lidia_size_t[r];
	memory_handler(A.allocated, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (allocated)");

	for (register lidia_size_t i = 0; i < r; i++) {
		A.value[i] = NULL;
		A.value_counter[i] = 0;
		A.allocated[i] = 0;
		A.index[i] = NULL;
	}
}



template< class T >
void
sparse_base_matrix_kernel< T >::constructor (MR< T > &A,
					     lidia_size_t r,
					     lidia_size_t c,
					     const T **B) const
{
	A.rows = r;
	A.columns = c;

	A.sparse_rows = A.sparse_columns = 0;

	A.value = new T*[r];
	memory_handler(A.value, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) ::"
		       "Error in memory allocation (A.value)");

	A.index = new lidia_size_t*[r];
	memory_handler(A.index, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) ::"
		       "Error in memory allocation (A.index)");

	A.value_counter = new lidia_size_t[r];
	memory_handler(A.value_counter, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) ::"
		       "Error in memory allocation (A.value_counter)");

	A.allocated = new lidia_size_t[r];
	memory_handler(A.allocated, DMESSAGE,
		       "constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) ::"
		       "Error in memory allocation (A.allocated)");

	for (register lidia_size_t i = 0; i < r; i++) {
		register lidia_size_t size = c;
		for (register lidia_size_t j = 0; j < c; j++)
			if (B[i][j] == A.Zero)
				size--;
		register lidia_size_t p = 0;
		if (size) {
			A.index[i] = new lidia_size_t[size];
			memory_handler(A.index[i], DMESSAGE,
				       "constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) ::"
				       "Error in memory allocation (A.index[i])");
			A.value[i] = new T[size];
			memory_handler(A.value[i], DMESSAGE,
				       "constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) ::"
				       "Error in memory allocation (A.value[i])");
			A.value_counter[i] = A.allocated[i] = size;
			for (register lidia_size_t k = 0; k < c; k++)
				if (B[i][k] != A.Zero) {
					A.value[i][p] = B[i][k];
					A.index[i][p] = k;
					p++;
				}
		}
	}
}



template< class T >
void
sparse_base_matrix_kernel< T >::constructor (MR< T > &A,
					     const MR< T > &B) const
{
	A.rows = B.rows;
	A.columns = B.columns;

	A.sparse_rows = B.sparse_rows;
	A.sparse_columns = B.sparse_columns;

	A.value = new T*[A.rows];
	memory_handler(A.value, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.value)");

	A.index = new lidia_size_t*[A.rows];
	memory_handler(A.index, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.index)");

	A.value_counter = new lidia_size_t[A.rows];
	memory_handler(A.value_counter, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.value_counter)");

	A.allocated = new lidia_size_t[A.rows];
	memory_handler(A.allocated, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.allocated)");

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		register lidia_size_t size = B.allocated[i];
		if (size) {
			A.index[i] = new lidia_size_t[size];
			memory_handler(A.index[i], DMESSAGE,
				       "constructor((MR< T > &, const MR< T > &) ::"
				       "Error in memory allocation (index[i])");
			A.value[i] = new T[size];
			memory_handler(A.value[i], DMESSAGE,
				       "constructor((MR< T > &, const MR< T > &) ::"
				       "Error in memory allocation (A.value[i])");
			A.value_counter[i] = B.value_counter[i];
			A.allocated[i] = size;
			for (register lidia_size_t p = 0; p < A.value_counter[i]; p++) {
				A.value[i][p] = B.value[i][p];
				A.index[i][p] = B.index[i][p];
			}
		}
		else {
		  A.value[i] = 0;
		  A.index[i] = 0;
		  A.allocated[i] = A.value_counter[i] = 0;
               }	
	}
}



//
// destructor kernel
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::destructor (MR< T > &A) const
{
	lidia_size_t l = ((A.bitfield.get_orientation() == matrix_flags::column_oriented) ? A.columns : A.rows);

	for (--l; l >= 0; l--) {
		if (A.allocated[l]) {
			delete[] A.value[l];
			delete[] A.index[l];
		}
	}

	if (A.index)
		delete[] A.index;
	if (A.value)
		delete[] A.value;
	if (A.allocated)
		delete[] A.allocated;
	if (A.value_counter)
		delete[] A.value_counter;
}



//
// column access kernel
//

template< class T >
void
sparse_base_matrix_kernel< T >::sto_column (MR< T > &A,
					    const T *v,
					    lidia_size_t l,
					    lidia_size_t j,
					    lidia_size_t from) const
{
	register lidia_size_t p, f;
	for (register lidia_size_t r = from; r < from + l; r++) { // for every row
		for (p = A.value_counter[r] - 1; p >= 0; p--)
			if (A.index[r][p] == j)
				break;

		if (v[r - from] != A.Zero) {
			if (p >= 0) {        // if not zero : replace
				A.value[r][p] = v[r - from];
				continue;
			}
			if (A.allocated[r] == A.value_counter[r]) { // if every place is used:
				// resize
				lidia_size_t len = (A.allocated[r] == 0) ?
					1 : A.allocated[r] * DELTA;
				T *tmp1 = new T[len];
				memory_handler(tmp1, DMESSAGE,
					       "sto_column(MR< T > &, const T *, "
					       "lidia_size_t, lidia_size_t, lidia_size_t) ::"
					       "Error in memory allocation (tmp1)");
				lidia_size_t *tmp2 = new lidia_size_t[len];
				memory_handler(tmp2, DMESSAGE,
					       "sto_column(MR< T > &, const T *, "
					       "lidia_size_t, lidia_size_t, lidia_size_t) ::"
					       "Error in memory allocation (tmp2)");

				for (register lidia_size_t c = 0; c < A.allocated[r]; c++) {
					tmp1[c] = A.value[r][c];
					tmp2[c] = A.index[r][c];
				}

				delete[] A.value[r];
				delete[] A.index[r];
				A.value[r] = tmp1;
				A.index[r] = tmp2;
				A.allocated[r] = len;
			}

			for (p = 0; p < A.value_counter[r]; p++)
				if (A.index[r][p] > j)
					break;

			for (f = A.value_counter[r] - 1; f >= p; f--) {
				A.index[r][f+1] = A.index[r][f];
				A.value[r][f+1] = A.value[r][f];
			}
			A.value[r][p] = v[r - from]; // and insert
			A.index[r][p] = j;
			A.value_counter[r]++;
		}
		else {
			if (p == -1)         // if already 0 do nothing
				continue;
			else {
				for (register lidia_size_t i = p + 1; i < A.value_counter[r]; i++) {
					// close the gap
					A.value[r][i-1] = A.value[r][i];
					A.index[r][i-1] = A.index[r][i];
				}
				A.value_counter[r]--;
			}
		}
	}
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::get_column (const MR< T > &A,
					    T *res,
					    lidia_size_t i) const
{
	for (register lidia_size_t j = 0; j < A.rows; j++) {
		lidia_size_t k;
		for (k = 0; k< A.value_counter[j] && i > A.index[j][k]; k++);
		if (k < A.value_counter[j] && i == A.index[j][k])
			res[j] = A.value[j][k];
		else {
			res[j] = A.Zero;
		}
	}
}



//
// row access kernel
//

template< class T >
void
sparse_base_matrix_kernel< T >::sto_row (MR< T > &A,
					 const T *v,
					 lidia_size_t l,
					 lidia_size_t i,
					 lidia_size_t from) const
{
	lidia_size_t k, j, c = 0;
	for (j = 0; j < l; j++)
		if (v[j] != A.Zero)
			c++;

	T *tmp = new T[c + A.value_counter[i]];
	memory_handler(tmp, DMESSAGE,
		       "sto_row(MR< T > &A, const T *v, lidia_size_t l, "
		       "lidia_size_t i, lidia_size_t from) :: "
		       "Error in memory allocation (tmp)");

	lidia_size_t *tmp1 = new lidia_size_t[c + A.value_counter[i]];
	memory_handler(tmp1, DMESSAGE,
		       "sto_row(MR< T > &A, const T *v, lidia_size_t l, "
		       "lidia_size_t i, lidia_size_t from) :: "
		       "Error in memory allocation (tmp1)");

	j = 0;
	while (A.index[i][j] < from) {
		tmp[j] = A.value[i][j];
		tmp1[j] = A.index[i][j];
		j++;
	}

	c = j;
	for (k = 0; k < l; k++) {
		if (v[k] != A.Zero) {
			tmp[j] = v[k];
			tmp1[j] = from + k;
			j++;
		}
	}

	while (c < A.value_counter[i]) {
		if (A.index[i][c] >= from + l) {
			tmp[j] = A.value[i][c];
			tmp1[j] = A.index[i][c];
			j++;
		}
		c++;
	}

	delete[] A.value[i];
	A.value[i] = tmp;
	delete[] A.index[i];
	A.index[i] = tmp1;

	A.allocated[i] = j;
	A.value_counter[i] = j;
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::get_row (const MR< T > &A,
					 T *res,
					 lidia_size_t i) const
{
	lidia_size_t l = 0, j;
	for (j = 0; j < A.columns; j++) {
		if (l < A.value_counter[i] && j == A.index[i][l]) {
			res[j] = A.value[i][l];
			l++;
		}
		else {
			res[j] = A.Zero;
		}
	}
}



//
// insert and remove kernel
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::insert_columns (MR< T > &A,
						lidia_size_t *ind,
						const T **news) const
{
	ind = NULL;
	news = NULL;
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::insert_column_at (MR< T > &M,
						  lidia_size_t j,
						  const T *v,
						  lidia_size_t len) const
{
	lidia_size_t i;

	for (i = 0; i < len; i++)
		sto(M, i, j, v[i]);
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::remove_columns (MR< T > &A,
						lidia_size_t *rem) const
{
	lidia_size_t i, j, l = 0, l1 = 1;
	for (i = 0; i < A.rows; i++) {
		l = 0;
		l1 = 1;

		for (j = 0; j < A.value_counter[i]; j++) {
			if (l1 > rem[0] || rem[l1] > A.index[i][j]) {
				LiDIA::swap(A.value[i][j], A.value[i][l]);
				A.index[i][l] = A.index[i][j] - l1 + 1;
				l++;
			}
			else {
				if (rem[l1] != A.index[i][j])
					j--;
				l1++;
			}
		}
		A.value_counter[i] = l;
	}
	A.columns -= rem[0];
}



template< class T >
void
sparse_base_matrix_kernel< T >::insert_rows (MR< T > &A,
					     lidia_size_t *ind,
					     const T **news) const
{
	ind = NULL;
	news = NULL;

	ind = NULL;
	news = NULL;

}



template< class T >
inline void
sparse_base_matrix_kernel< T >::insert_row_at (MR< T > &M,
					       lidia_size_t j,
					       const T *v,
					       lidia_size_t len) const
{
	lidia_size_t i;

	for (i = 0; i < len; i++)
		sto(M, j, i, v[i]);
}



template< class T >
void
sparse_base_matrix_kernel< T >::remove_rows (MR< T > &A,
					     lidia_size_t *rem) const
{
	lidia_size_t len = A.rows;
	A.rows -= rem[0];

	T** new_value = new T*[A.rows];
	lidia_size_t **index = new lidia_size_t *[A.rows];
	lidia_size_t *allocated = new lidia_size_t[A.rows];
	lidia_size_t *value_counter = new lidia_size_t[A.rows];

	register lidia_size_t i, l = 0, l1 = 1;
	for (i = 0; i < len; i++) {
		if (l1 > rem[0] || rem[l1] != i) {
			new_value[l] = A.value[i];
			index[l] = A.index[i];
			allocated[l] = A.allocated[i];
			value_counter[l] = A.value_counter[i];
			l++;
		}
		else {
			if (A.allocated[i] != 0) {
				delete[] A.value[i];
				delete[] A.index[i];
			}

			l1++;
		}
	}

	delete[] A.value;
	delete[] A.index;
	delete[] A.value_counter;
	delete[] A.allocated;

	A.value = new_value;
	A.index = index;
	A.value_counter = value_counter;
	A.allocated = allocated;

}



//
// exchange functions / swap functions
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::swap_columns (MR< T > &RES,
					      lidia_size_t i,
					      lidia_size_t j) const
{
	T TMP;
	for (lidia_size_t l = 0; l < RES.rows; l++) {
		TMP = member(RES, l, i);
		sto(RES, l, i, member(RES, l, j));
		sto(RES, l, j, TMP);
	}
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::swap_rows (MR< T > &RES,
					   lidia_size_t i,
					   lidia_size_t j) const
{
	// value array
	T *tmp = RES.value[i];
	RES.value[i] = RES.value[j];
	RES.value[j] = tmp;

	// index array
	lidia_size_t *tmp1 = RES.index[i];
	RES.index[i] = RES.index[j];
	RES.index[j] = tmp1;

	// and the rest
	LiDIA::swap(RES.value_counter[i], RES.value_counter[j]);
	LiDIA::swap(RES.allocated[i], RES.allocated[j]);
}



//
// assignment kernel
//

template< class T >
void
sparse_base_matrix_kernel< T >::assign (MR< T > &A,
					const MR< T > &B) const
{
	lidia_size_t l = A.rows;

	if (l > 0) {
		for (--l; l >= 0; l--) {
			if (A.allocated[l]) {
				if (A.value[l])
					delete[] A.value[l];
				if (A.index[l])
					delete[] A.index[l];
			}
		}
		if (A.index)
			delete[] A.index;
		if (A.value)
			delete[] A.value;
		if (A.allocated)
			delete[] A.allocated;
		if (A.value_counter)
			delete[] A.value_counter;
	}

	A.rows = B.rows;
	A.columns = B.columns;

	A.sparse_rows = B.sparse_rows;
	A.sparse_columns = B.sparse_columns;

	A.bitfield = B.bitfield;

	A.value = new T*[A.rows];
	memory_handler(A.value, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.value)");

	A.index = new lidia_size_t*[A.rows];
	memory_handler(A.index, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.index)");

	A.value_counter = new lidia_size_t[A.rows];
	memory_handler(A.value_counter, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.value_counter)");

	A.allocated = new lidia_size_t[A.rows];
	memory_handler(A.allocated, DMESSAGE,
		       "constructor((MR< T > &, const MR< T > &) ::"
		       "Error in memory allocation (A.allocated)");

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		register lidia_size_t size = B.allocated[i];
		if (size) {
			A.index[i] = new lidia_size_t[size];
			memory_handler(A.index[i], DMESSAGE,
				       "constructor((MR< T > &, const MR< T > &) ::"
				       "Error in memory allocation (index[i])");
			A.value[i] = new T[size];
			memory_handler(A.value[i], DMESSAGE,
				       "constructor((MR< T > &, const MR< T > &) ::"
				       "Error in memory allocation (A.value[i])");
			A.value_counter[i] = B.value_counter[i];
			A.allocated[i] = size;
			for (register lidia_size_t p = 0; p < A.value_counter[i]; p++) {
				A.value[i][p] = B.value[i][p];
				A.index[i][p] = B.index[i][p];
			}
		}
	}
}



//
// stream handling - LIDIA
//

template< class T >
void
sparse_base_matrix_kernel< T >::write_to_beauty (const MR< T > &A,
						 std::ostream &out) const
{
	if (A.bitfield.get_orientation() == matrix_flags::row_oriented) {
		register lidia_size_t i, j, l;
		for (i = 0; i < A.rows; i++) {
			l = 0;
			out << std::endl << "(";
			for (j = 0; j < A.columns; j++)
				if (l < A.value_counter[i] && A.index[i][l] == j) {
					out << A.value[i][l] << " ";
					l++;
				}
				else
					out << A.Zero << " ";
			out << ")";
		}
		out << std::endl << std::flush;
	}
	else {
		register lidia_size_t i, j;
		lidia_size_t *l = new lidia_size_t[A.columns];
		for (i = 0; i < A.columns; i++)
			l[i] = 0;

		for (i = 0; i < A.rows; i++) {
			out << std::endl << "(";
			for (j = 0; j < A.columns; j++)
				if (l[j] < A.value_counter[j] && A.index[j][l[j]] == i) {
					out << A.value[j][l[j]] << " ";
					l[j]++;
				}
				else
					out << A.Zero << " ";
			out << ")";
		}
		out << std::endl << std::flush;
		delete[] l;
	}
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::write_to_stream (const MR< T > &A,
						 std::ostream &out) const
{
	out << A.rows << std::endl;
	out << A.columns << std::endl;

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		out << A.value_counter[i] << " ";
		for (register lidia_size_t j = 0; j < A.value_counter[i]; j++)
			out << A.index[i][j] << " " << A.value[i][j] << " ";
		out << std::endl;
	}
}



template< class T >
void
sparse_base_matrix_kernel< T >::read_from_stream (MR< T > &A,
						  std::istream &in) const
{
    if(in.good()) {
	lidia_size_t new_rows;
	in >> std::ws >> new_rows;
	if(new_rows < 0) {
	    in.setstate(std::ios::failbit);
	}
	if(in.fail()) {
	    return;
	}

	lidia_size_t new_columns;
	in >> std::ws >> new_columns;
	if(new_columns < 0) {
	    in.setstate(std::ios::failbit);
	}
	if(in.fail()) {
	    return;
	}

	lidia_size_t l = A.rows;
	if (l > 0) {
		for (--l; l >= 0; l--) {
			if (A.allocated[l]) {
				if (A.value[l])
					delete[] A.value[l];
				if (A.index[l])
					delete[] A.index[l];
			}
		}
		if (A.index)
			delete[] A.index;
		if (A.value)
			delete[] A.value;
		if (A.allocated)
			delete[] A.allocated;
		if (A.value_counter)
			delete[] A.value_counter;
	}

	A.rows = new_rows;
	A.columns = new_columns;
	A.value = new T*[A.rows];
	memory_handler(A.value, DMESSAGE,
		       "read_from_stream ::"
		       "Error in memory allocation (A.value)");

	A.index = new lidia_size_t*[A.rows];
	memory_handler(A.index, DMESSAGE,
		       "constructor(T **&, lidia_size_t &, lidia_size_t &, lidia_size_t **&, "
		       "lidia_size_t *&, lidia_size_t *&, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (A.index)");

	A.value_counter = new lidia_size_t[A.rows];
	memory_handler(A.value_counter, DMESSAGE,
		       "constructor(T **&, lidia_size_t &, lidia_size_t &, lidia_size_t **&, "
		       "lidia_size_t *&, lidia_size_t *&, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (A.value_counter)");

	A.allocated = new lidia_size_t[A.rows];
	memory_handler(A.allocated, DMESSAGE,
		       "constructor(T **&, lidia_size_t &, lidia_size_t &, lidia_size_t **&, "
		       "lidia_size_t *&, lidia_size_t *&, lidia_size_t, lidia_size_t) ::"
		       "Error in memory allocation (allocated)");

	lidia_size_t size;
	for (register lidia_size_t i = 0; i < A.rows; i++) {
		in >> std::ws >> size;
		if(size < 0) {
		    in.setstate(std::ios::failbit);
		    size = 0;
		}
		A.allocated[i] = A.value_counter[i] = size;

		A.index[i] = new lidia_size_t[size];
		memory_handler(A.index[i], DMESSAGE,
			       "constructor(T **&, lidia_size_t &, lidia_size_t &, lidia_size_t **&, "
			       "lidia_size_t *&, lidia_size_t *&, lidia_size_t, lidia_size_t, const T **) ::"
			       "Error in memory allocation (A.index[i])");

		A.value[i] = new T[size];
		memory_handler(A.value[i], DMESSAGE,
			       "constructor(T **&, lidia_size_t &, lidia_size_t &, lidia_size_t **&, "
			       "lidia_size_t *&, lidia_size_t *&, lidia_size_t, lidia_size_t, const T **) ::"
			       "Error in memory allocation (A.value[i])");

		for (register lidia_size_t j = 0; j < size; j++) {
			in >> std::ws >> A.index[i][j]
			   >> std::ws >> A.value[i][j];
		}
	}
    }
}



//
// stream handling - MATHEMATICA
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::write_to_mathematica (const MR< T > &A,
						      std::ostream &out) const
{

	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1, l1;

	out << "{";
	for (i = 0; i < A.rows; i++) {
		l1 = 0;
		out << "{";
		for (j = 0; j < A.columns; j++) {
			if (l1 < A.value_counter[i] && A.index[i][l1] == j) {
				out << A.value[i][l1];
				l1++;
			}
			else
				out << A.Zero;
			if (j != l)
				out << ", ";
			else
				if (i != k)
					out << "}, " << std::endl;
		}
	}
	out << "}}" << std::endl << std::flush;
}



template< class T >
void
sparse_base_matrix_kernel< T >::read_from_mathematica (MR< T > &A,
						       std::istream &dz) const
{
	char c = 0;

	dz >> std::ws >> c;
	if (c != '{')
		lidia_error_handler("void MR< T >::"
				    "read_from_mathematica(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	base_vector< T > buffer(0, vector_flags(vector_flags::expand));
	lidia_size_t len = 0, local_rows = 0, local_columns = 0;
	bool END = false;

	dz >> std::ws >> c;
	while (!dz.eof() && END == false) {
		if (c == '{') {
			local_rows++;
			while (c != '}') {
				dz >> buffer[len];
				len++;
				dz >> std::ws >> c;
			}
		}
		dz >> std::ws >> c;
		if (c == '}')
			END = true;
	}

	local_columns = len / local_rows;

	len = 0;
	if (A.rows != local_rows)
		set_no_of_rows(A, local_rows);
	if (A.columns != local_columns)
		set_no_of_columns(A, local_columns);

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++) {
			sto(A, i, j, buffer[len]);
			len++;
		}
	}
}



//
// stream handling - MAPLE
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::write_to_maple (const MR< T > &A,
						std::ostream &out) const
{
	register lidia_size_t i, j, l;

	out << "array(1 .. " << A.rows << ", 1 .. " << A.columns << ", [";

	for (i = 0; i < A.rows; i++) {
		l = 0;
		for (j = 0; j < A.columns; j++)
			if (l < A.value_counter[i]) {
				if (A.index[i][l] == j) {
					out << "(" << i+1 << ", " << j+1 << ") = " << A.value[i][l];
					l++;
				}
				else
					out << "(" << i+1 << ", " << j+1 << ") = " << A.Zero;

				if (!(j == A.columns - 1 && i == A.rows - 1))
					out << ", ";
			}
	}
	out << "]); " << std::endl << std::flush;
}



template< class T >
void
sparse_base_matrix_kernel< T >::read_from_maple (MR< T > &A,
						 std::istream &dz) const
{
	char c = 0;
	while (c != '(' && dz.good())
		dz >> std::ws >> c;

	if (dz.fail())
		lidia_error_handler("void MR< T >::"
				    "read_from_maple(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	register lidia_size_t i, j;
	lidia_size_t startr, startc, endr, endc, x, y;
	T TMP;

	// read row dimension
	dz >> startr >> std::ws >> c >> std::ws >> c >> endr >> std::ws >> c;

	// read column dimension
	dz >> startc >> std::ws >> c >> std::ws >> c >> endc >> std::ws >> c;

	endr = endr - startr + 1;
	if (endr < 0)
		endr = -endr;
	endc = endc - startc + 1;
	if (endc < 0)
		endc = -endc;

	if (A.columns != endc)
		set_no_of_columns(A, endc);
	if (A.rows != endr)
		set_no_of_rows(A, endr);

	dz >> std::ws >> c;

	for (i = 0; i < endc; i++)
		for (j = 0; j < endr; j++) {
			dz >> std::ws >> c >> x >> std::ws >> c >> y;
			dz >> std::ws >> c >> std::ws >> c >> std::ws >> TMP >> std::ws >> c;
			sto(A, x - startr, y - startc, TMP);
		}
	dz >> std::ws >> c >> std::ws >> c;
}



//
// stream handling - PARI
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::write_to_gp (const MR< T > &A,
					     std::ostream &out) const
{
	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1, l1;

	out << "[";
	for (i = 0; i < A.rows; i++) {
		l1 = 0;
		for (j = 0; j < A.columns; j++) {
			if (l1 < A.value_counter[i] && A.index[i][l1] == j) {
				out << A.value[i][l1];
				l1++;
			}
			else
				out << A.Zero;
			if (j != l)
				out << ", ";
		}
		if (i != k)
			out << "; ";
	}
	out << "]" << std::flush;
}



template< class T >
void
sparse_base_matrix_kernel< T >::read_from_gp (MR< T > &A,
					      std::istream &dz) const
{
	char c = 0;

	dz >> std::ws >> c;
	if (c != '[')
		lidia_error_handler("void MR< T >::"
				    "read_from_gp(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	base_vector< T > buffer(0, vector_flags(vector_flags::expand));
	lidia_size_t len = 0, local_rows = 0, local_columns = 0;
	bool END = false;

	while (!dz.eof() && END == false) {
		local_rows++;
		do {
			dz >> buffer[len];
			len++;
			dz >> std::ws >> c;
		}
		while (c != ';' && c != ']');

		if (c == ']')
			END = true;
	}

	local_columns = len / local_rows;

	len = 0;
	if (A.rows != local_rows)
		set_no_of_rows(A, local_rows);
	if (A.columns != local_columns)
		set_no_of_columns(A, local_columns);

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++) {
			sto(A, i, j, buffer[len]);
			len++;
		}
	}
}



//
// stream handling - KANT
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::write_to_kash (const MR< T > &A,
					       std::ostream &out) const
{
	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1, l1;

	out << "LIDIA: = Mat(Z, [";
	for (i = 0; i < A.rows; i++) {
		l1 = 0;
		out << "[";
		for (j = 0; j < A.columns; j++) {
			if (l1 < A.value_counter[i] && A.index[i][l1] == j) {
				out << A.value[i][l1];
				l1++;
			}
			else
				out << A.Zero;

			if (j != l)
				out << ", ";
			else
				if (i != k)
					out << "], " << std::endl;
		}
	}
	out << "]]); " << std::endl << std::flush;
}



template< class T >
void
sparse_base_matrix_kernel< T >::read_from_kash (MR< T > &A,
						std::istream &dz) const
{
	char c = 0;
	do {
		dz >> std::ws >> c;
	}
	while (c != 'M' || dz.eof());
	if (c != 'M')
		lidia_error_handler("void MR< T >::"
				    "read_from_kash(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	dz >> std::ws >> c >> std::ws >> c >> std::ws >> c >> std::ws >> c >> std::ws >> c >> std::ws >> c;

	base_vector< T > buffer(0, vector_flags(vector_flags::expand));
	lidia_size_t len = 0, local_rows = 0, local_columns = 0;
	bool END = false;

	dz >> std::ws >> c;
	while (!dz.eof() && END == false) {
		if (c == '[') {
			local_rows++;
			while (c != ']') {
				dz >> buffer[len];
				len++;
				dz >> std::ws >> c;
			}
		}
		dz >> std::ws >> c;
		if (c == ']')
			END = true;
	}

	local_columns = len / local_rows;

	len = 0;
	if (A.rows != local_rows)
		set_no_of_rows(A, local_rows);
	if (A.columns != local_columns)
		set_no_of_columns(A, local_columns);

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++) {
			sto(A, i, j, buffer[len]);
			len++;
		}
	}
	dz >> std::ws >> c >> std::ws >> c;
}



//
// stream handling - LATEX
//

template< class T >
inline void
sparse_base_matrix_kernel< T >::write_to_latex (const MR< T > &A,
						std::ostream &out) const
{
	out << A.value[0][0];
}



//
// stream handling - MAGMA
//

template< class T >
void
sparse_base_matrix_kernel< T >::write_to_magma (const MR< T > &A,
						std::ostream &out) const
{
	lidia_size_t i, j;
	out << "L : = RMatrixSpace(Integers(), " << A.rows << ", " << A.columns << ") ! [" << std::endl;
	for (i = 0; i < A.rows; i++) {
		for (j = 0; j < A.columns; j++) {
			out << member(A, i, j) << std::flush;
			if (j < A.columns - 1)
				out << ", " << std::flush;
			else if (j == A.columns - 1 && i < A.rows - 1)
				out << ", " << std::endl;
			else if (j == A.columns - 1 && i == A.rows - 1)
				out << "]; " << std::endl;
		}
	}
}



//
// structur functions
//

template< class T >
void
sparse_base_matrix_kernel< T >::set_no_of_rows (MR< T > &A,
						lidia_size_t r) const
{
	if (r == A.rows)
		return;

	T **tmp = A.value;
	lidia_size_t **tmp1 = A.index;
	lidia_size_t i;

	if (r < A.rows) {
		for (i = r; i < A.rows; i++) {
			delete[] A.value[i];
			delete[] A.index[i];
		}

		A.value = new T *[r];
		memory_handler(A.value, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.value)");
		A.index = new lidia_size_t *[r];
		memory_handler(A.index, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.index)");

		for (i = 0; i < r; i++) {
			A.value[i] = tmp[i];
			A.index[i] = tmp1[i];
		}
	}
	else {
		A.value = new T *[r];
		memory_handler(A.value, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.value)");

		A.index = new lidia_size_t *[r];
		memory_handler(A.index, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.index)");

		for (i = 0; i < A.rows; i++) {
			A.value[i] = tmp[i];
			A.index[i] = tmp1[i];
		}

		lidia_size_t *tmp2 = A.allocated;
		lidia_size_t *tmp3 = A.value_counter;

		A.allocated = new lidia_size_t[r];
		memory_handler(A.allocated, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.allocated)");

		A.value_counter = new lidia_size_t[r];
		memory_handler(A.value_counter, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.value_counter)");

		for (i = 0; i < A.rows; i++) {
			A.allocated[i] = tmp2[i];
			A.value_counter[i] = tmp3[i];
		}

		for (i = A.rows; i < r; i++) {
			A.index[i] = NULL;
			A.value[i] = NULL;
			A.allocated[i] = 0;
			A.value_counter[i] = 0;
		}

		delete[] tmp2;
		delete[] tmp3;
	}

	A.rows = r;
	delete[] tmp;
	delete[] tmp1;
}



template< class T >
inline void
sparse_base_matrix_kernel< T >::set_no_of_columns (MR< T > &A,
						   lidia_size_t c) const
{
	if (c < A.columns)
		for (register lidia_size_t i = 0; i < A.rows; i++)
			while (A.value_counter[i] > 0 && A.index[i][A.value_counter[i] - 1] >= c)
				A.value_counter[i]--;
	A.columns = c;
}



//
// element access kernel
//

template< class T >
void
sparse_base_matrix_kernel< T >::sto (MR< T > &A,
				     lidia_size_t x,
				     lidia_size_t y,
				     const T &e) const
{
	T wert = e;
	register lidia_size_t p, i;
	if (e != A.Zero) {
		for (p = A.value_counter[x] - 1; p >= 0 && A.index[x][p] >= y; p--) {
			if (A.index[x][p] == y) {
				A.value[x][p] = wert;
				return;
			}
		}

		if (A.allocated[x] == A.value_counter[x]) {	// if every place is used : resize

			bool SW = true;
			if (A.allocated[x] == 0) {
				SW = false;
				A.allocated[x] = 1;
			}

			lidia_size_t len = ((A.allocated[x]*DELTA < A.columns) ? A.allocated[x]*DELTA : A.columns);
			register T *tmp1 = new T[len];
			memory_handler(tmp1, DMESSAGE,
				       "sto(MR< T > &, lidia_size_t, "
				       "lidia_size_t, const T &) :: "
				       "Error in memory allocation (tmp1)");

			register lidia_size_t *tmp2 = new lidia_size_t[len];
			memory_handler(tmp2, DMESSAGE,
				       "sto(MR< T > &, lidia_size_t, "
				       "lidia_size_t, const T &) :: "
				       "Error in memory allocation (tmp2)");

			for (i = A.value_counter[x] - 1; i >= 0; i--) {
				LiDIA::swap(tmp1[i], A.value[x][i]);
				LiDIA::swap(tmp2[i], A.index[x][i]);
			}

			A.allocated[x] = len;

			if (SW) {
				delete[] A.value[x];
				delete[] A.index[x];
			}

			A.value[x] = tmp1;
			A.index[x] = tmp2;

		}

		register T *tmp1 = A.value[x];
		register lidia_size_t *tmp2 = A.index[x];

		for (i = A.value_counter[x]; i - 1 > p; i--) {
			LiDIA::swap(tmp1[i], tmp1[i - 1]);
			LiDIA::swap(tmp2[i], tmp2[i - 1]);
		}

		tmp1[p + 1] = wert;
		tmp2[p + 1] = y;

		A.value_counter[x]++;
	}
	else {
		for (p = A.value_counter[x] - 1; p >= 0 && A.index[x][p] >= y; p--) {
			if (A.index[x][p] == y) {
				for (i = p; i + 1 < A.value_counter[x]; i++) {
					LiDIA::swap(A.value[x][i], A.value[x][i + 1]);
					LiDIA::swap(A.index[x][i], A.index[x][i + 1]);
				}
				A.value_counter[x]--;
			}
		}
	}
}



template< class T >
inline T & sparse_base_matrix_kernel< T >::member (const MR< T > &A,
						   lidia_size_t x,
						   lidia_size_t y) const
{
	register lidia_size_t p;
	for (p = A.value_counter[x]-1; p >= 0; p--)
		if (A.index[x][p] == y)
			return A.value[x][p];
	return (T &)A.Zero;
}



//
// boolean functions
//

template< class T >
inline bool sparse_base_matrix_kernel< T >::is_column_zero (const MR< T > &RES,
							    lidia_size_t c) const
{
	register lidia_size_t i, j;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.value_counter[i] && RES.index[i][j] < c; j++);
		if (j != RES.value_counter[i] && c == RES.index[i][j])
			return false;
	}
	return true;
}



template< class T >
inline bool sparse_base_matrix_kernel< T >::is_row_zero (const MR< T > &RES,
							 lidia_size_t r) const
{
	register lidia_size_t i;
	register T *tmp = RES.value[r];

	for (i = 0; i < RES.value_counter[r]; i++) {
		if (tmp[i] != RES.Zero)
			return false;
	}
	return true;
}



template< class T >
inline bool sparse_base_matrix_kernel< T >::is_matrix_zero (const MR< T > &RES) const
{
	register lidia_size_t i, j;
	register T *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.value_counter[i]; j++)
			if (tmp[j] != RES.Zero)
				return false;
	}
	return true;
}



//
// change orientation
//

template< class T >
void
sparse_base_matrix_kernel< T >::change_orientation (MR< T > &A,
						    unsigned long mode) const
{
	lidia_size_t index1 = 0, index2 = 0;
	if (mode != A.bitfield.get_orientation()) {
		if (mode == matrix_flags::column_oriented) {
			index1 = A.columns;
			index2 = A.rows;
			A.bitfield.set_orientation(matrix_flags::column_oriented);
		}
		if (mode == matrix_flags::row_oriented) {
			index1 = A.rows;
			index2 = A.columns;
			A.bitfield.set_orientation(matrix_flags::row_oriented);
		}

		T **oldvalue = A.value;
		A.value = new T*[index1];
		memory_handler(A.value, DMESSAGE,
			       "change_orientation(MR< T > &, unsigned long) ::"
			       "Error in memory allocation (A.value)");

		lidia_size_t **oldindex = A.index;
		A.index = new lidia_size_t*[index1];
		memory_handler(A.index, DMESSAGE,
			       "change_orientation(MR< T > &, unsigned long) ::"
			       "Error in memory allocation (A.index)");

		lidia_size_t *oldvalue_counter = A.value_counter;
		A.value_counter = new lidia_size_t[index1];
		memory_handler(A.value_counter, DMESSAGE,
			       "change_orientation(MR< T > &, unsigned long) ::"
			       "Error in memory allocation (A.value_counter)");

		lidia_size_t *oldallocated = A.allocated;
		A.allocated = new lidia_size_t[index1];
		memory_handler(A.allocated, DMESSAGE,
			       "change_orientation(MR< T > &, unsigned long) ::"
			       "Error in memory allocation (allocated)");

		register lidia_size_t j; // MM for HP UX 10.20 CC
		for (j = 0; j < index1; j++) {
			A.allocated[j] = 0;
			A.value_counter[j] = 0;
		}

		for (j = index2 - 1; j >= 0; j--) {
			for (register lidia_size_t i = 0; i < oldvalue_counter[j]; i++)
				A.value_counter[oldindex[j][i]]++;
		}

		for (j = 0; j < index1; j++) {
			if (A.value_counter[j] > 0) {
				A.value[j] = new T[A.value_counter[j]];
				memory_handler(A.value[j], DMESSAGE,
					       "change_orientation(MR< T > &, unsigned long) ::"
					       "Error in memory allocation (A.value[i])");

				A.index[j] = new lidia_size_t[A.value_counter[j]];
				memory_handler(A.index[j], DMESSAGE,
					       "change_orientation(MR< T > &, unsigned long) ::"
					       "Error in memory allocation (A.index[i])");
			}
			else {
				A.value[j] = NULL;
				A.index[j] = NULL;
			}
		}
		register lidia_size_t l;
		for (j = 0; j < index2; j++) {
			for (register lidia_size_t i = 0; i < oldvalue_counter[j]; i++) {
				l = oldindex[j][i];
				A.value[l][A.allocated[l]] = oldvalue[j][i];
				A.index[l][A.allocated[l]] = j;
				A.allocated[l]++;
			}
			delete[] oldvalue[j];
			delete[] oldindex[j];
		}
		delete[] oldvalue;
		delete[] oldindex;
		delete[] oldvalue_counter;
		delete[] oldallocated;
	}
}



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE

#undef sparse_base_matrix_kernel



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_BASE_MATRIX_KERNEL_CC_GUARD_
