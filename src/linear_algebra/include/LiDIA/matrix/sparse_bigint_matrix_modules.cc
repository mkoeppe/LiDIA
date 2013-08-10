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


#ifndef LIDIA_SPARSE_BIGINT_MATRIX_MODULES_CC_GUARD_
#define LIDIA_SPARSE_BIGINT_MATRIX_MODULES_CC_GUARD_


#include	"LiDIA/arith.inl"
#include	"LiDIA/modular_operations.inl"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#ifdef LIDIA_NAMESPACE
using std::abs;
#endif



#define column_oriented_sparse_matrix_modules COSMM

#define DELTA 2

//
// member
//

template< class T >
const T &
column_oriented_sparse_matrix_modules< T >::member (const MR< T > &A,
						    lidia_size_t x,
						    lidia_size_t y) const
{
	for (register lidia_size_t i = 0; i < A.value_counter[y]; i++) {
		if (A.index[y][i] == x)
			return A.value[y][i];
	}
	return A.Zero;
}



//
// sto
//

template< class T >
void
column_oriented_sparse_matrix_modules< T >::sto (MR< T > &A,
						 lidia_size_t y,
						 lidia_size_t x,
						 const T &e) const
{
	register lidia_size_t p, i;
	if (e != A.Zero) {
		for (p = A.value_counter[x] - 1; p >= 0 && A.index[x][p] >= y; p--) {
			if (A.index[x][p] == y) {
				A.value[x][p] = e;
				return;
			}
		}

		if (A.allocated[x] == A.value_counter[x]) {
			// if every place is used : resize

			bool SW = true;
			if (A.allocated[x] == 0) {
				SW = false;
				A.allocated[x] = 1;
			}

			lidia_size_t len = ((A.allocated[x]*DELTA < A.columns) ? A.allocated[x]*DELTA : A.rows);
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

		tmp1[p + 1] = e;
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



//
// swap rows
//

template< class T >
void
column_oriented_sparse_matrix_modules< T >::swap_columns (MR< T > &A,
							  lidia_size_t pos1,
							  lidia_size_t pos2) const
{
	T *tmp = A.value[pos1];
	A.value[pos1] = A.value[pos2];
	A.value[pos2] = tmp;

	lidia_size_t *tmp1 = A.index[pos1];
	A.index[pos1] = A.index[pos2];
	A.index[pos2] = tmp1;

	LiDIA::swap(A.allocated[pos1], A.allocated[pos2]);
	LiDIA::swap(A.value_counter[pos1], A.value_counter[pos2]);
}



//
// max
//

template< class T >
inline void
column_oriented_sparse_matrix_modules< T >::max_abs (const MR< T > &RES,
						     T &MAX) const
{
	register lidia_size_t i, j;
	bool SW = false;

	for (i = 0; i < RES.columns; i++) {
		for (j = 0; j < RES.value_counter[i]; j++)
			if (SW == false)
				MAX = abs(RES.value[i][j]);
			else
				if (MAX < abs(RES.value[i][j]))
					MAX = abs(RES.value[i][j]);
	}
}



//
// subtract
//

template< class T >
void
column_oriented_sparse_matrix_modules< T >::subtract_multiple_of_column (MR< T > &A,
									 lidia_size_t index1,
									 const T &q,
									 lidia_size_t index2,
									 lidia_size_t length) const
{
	lidia_size_t len = (A.value_counter[index1] + A.value_counter[index2] < length + 1 ?
			    A.value_counter[index1] + A.value_counter[index2] : length + 1);


	T *tmp = new T[len];
	lidia_size_t *tmp2 = new lidia_size_t[len];
	A.allocated[index1] = len;

	T TMP;

	register lidia_size_t l1, l2, l3;
	for (l1 = l2 = l3 = 0; l2 < A.value_counter[index1] && l1 < A.value_counter[index2];) {
		if (A.index[index2][l1] == A.index[index1][l2]) {
			LiDIA::multiply(TMP, q, A.value[index2][l1]);
			LiDIA::subtract(tmp[l3], A.value[index1][l2], TMP);

			tmp2[l3] = A.index[index2][l1];

			l1++;
			l2++;
			if (tmp[l3] != 0)
				l3++;
		}
		else if (A.index[index2][l1] < A.index[index1][l2]) {
			LiDIA::multiply(tmp[l3], -q, A.value[index2][l1]);
			tmp2[l3] = A.index[index2][l1];

			l1++;
			if (tmp[l3] != 0)
				l3++;
		}
		else if (A.index[index2][l1] > A.index[index1][l2]) {
			tmp[l3] = A.value[index1][l2];

			tmp2[l3] = A.index[index1][l2];

			l2++;
			if (tmp[l3] != 0)
				l3++;
		}
	}

	while (l2 < A.value_counter[index1]) {
		tmp[l3] = A.value[index1][l2];
		tmp2[l3] = A.index[index1][l2];
		l2++;
		l3++;
	}

	while (l1 < A.value_counter[index2]) {
		LiDIA::multiply(tmp[l3], -q, A.value[index2][l1]);
		tmp2[l3] = A.index[index2][l1];
		l1++;
		l3++;
	}


	if (A.value[index1] != NULL)
		delete[] A.value[index1];
	A.value[index1] = tmp;

	if (A.index[index1] != NULL)
		delete[] A.index[index1];
	A.index[index1] = tmp2;

	A.value_counter[index1] = l3;

#if 0
	tmp = new T[l3];
	tmp2 = new lidia_size_t[l3];

	// memory adaption
	for (l1 = 0; l1 < l3; l1++) {
		tmp2[l1] = A.index[index1][l1];
		tmp[l1] = A.value[index1][l1];
	}

	delete[] A.value[index1];
	A.value[index1] = tmp;

	delete[] A.index[index1];
	A.index[index1] = tmp2;
#endif
}



template< class T >
void
column_oriented_sparse_matrix_modules< T >::subtract_multiple_of_column_mod (MR< T > &A,
									     lidia_size_t index1,
									     const T &q,
									     lidia_size_t index2,
									     lidia_size_t length,
									     const T &mod) const
{

	lidia_size_t len = (A.value_counter[index1] + A.value_counter[index2] < length + 1 ?
			    A.value_counter[index1] + A.value_counter[index2] : length + 1);

	T *tmp = new T[len];
	lidia_size_t *tmp2 = new lidia_size_t[len];
	A.allocated[index1] = len;
	T TMP;

	register lidia_size_t l1, l2, l3;
	for (l1 = l2 = l3 = 0; l2 < A.value_counter[index1] && l1 < A.value_counter[index2];) {
		if (A.index[index2][l1] == A.index[index1][l2]) {
			LiDIA::mult_mod(TMP, q, A.value[index2][l1], mod);
			LiDIA::sub_mod(tmp[l3], A.value[index1][l2], TMP, mod);
			tmp2[l3] = A.index[index2][l1];

			l1++;
			l2++;
			if (tmp[l3] != 0)
				l3++;
		}
		else if (A.index[index2][l1] < A.index[index1][l2]) {
			LiDIA::mult_mod(tmp[l3], -q, A.value[index2][l1], mod);
			tmp2[l3] = A.index[index2][l1];

			l1++;
			if (tmp[l3] != 0)
				l3++;
		}
		else if (A.index[index2][l1] > A.index[index1][l2]) {
			tmp[l3] = A.value[index1][l2];

			tmp2[l3] = A.index[index1][l2];

			l2++;
			if (tmp[l3] != 0)
				l3++;
		}
	}

	while (l2 < A.value_counter[index1]) {
		tmp[l3] = A.value[index1][l2];
		tmp2[l3] = A.index[index1][l2];
		l2++;
		l3++;
	}

	while (l1 < A.value_counter[index2]) {

		LiDIA::mult_mod(tmp[l3], q, A.value[index2][l1], mod);
		tmp2[l3] = A.index[index2][l1];
		l1++;
		l3++;
	}

	delete[] A.value[index1];
	A.value[index1] = tmp;

	delete[] A.index[index1];
	A.index[index1] = tmp2;

	A.value_counter[index1] = l3;

	tmp = new T[l3];
	tmp2 = new lidia_size_t[l3];

	// memory adaption
	for (l1 = 0; l1 < l3; l1++) {
		tmp2[l1] = A.index[index1][l1];
		tmp[l1] = A.value[index1][l1];
	}

	delete[] A.value[index1];
	A.value[index1] = tmp;

	delete[] A.index[index1];
	A.index[index1] = tmp2;
}



template< class T >
void
column_oriented_sparse_matrix_modules< T >::normalize_column_mod (MR< T > &A,
								  lidia_size_t index1,
								  const T &q,
								  lidia_size_t index2,
								  lidia_size_t length,
								  const T &mod) const
{
	lidia_size_t len = (A.value_counter[index1] + A.value_counter[index2] < length + 1 ?
			    A.value_counter[index1] + A.value_counter[index2] : length + 1);

	T *tmp = new T[len];
	lidia_size_t *tmp2 = new lidia_size_t[len];
	A.allocated[index1] = len;
	T TMP;

	register lidia_size_t l1, l2, l3;
	for (l1 = l2 = l3 = 0; l2 < A.value_counter[index1] && l1 < A.value_counter[index2];) {
		if (A.index[index2][l1] == A.index[index1][l2]) {
			LiDIA::mult_mod(TMP, q, A.value[index2][l1], mod);
			LiDIA::sub_mod(tmp[l3], A.value[index1][l2], TMP, mod);
			if (tmp[l3] < A.Zero)
				tmp[l3] += mod;
			tmp2[l3] = A.index[index2][l1];

			l1++;
			l2++;
			if (tmp[l3] != 0)
				l3++;
		}
		else if (A.index[index2][l1] < A.index[index1][l2]) {
			LiDIA::mult_mod(tmp[l3], -q, A.value[index2][l1], mod);
			if (tmp[l3] < A.Zero)
				tmp[l3] += mod;
			tmp2[l3] = A.index[index2][l1];

			l1++;
			if (tmp[l3] != 0)
				l3++;
		}
		else if (A.index[index2][l1] > A.index[index1][l2]) {
			tmp[l3] = A.value[index1][l2];
			if (tmp[l3] < A.Zero)
				tmp[l3] += mod;

			tmp2[l3] = A.index[index1][l2];

			l2++;
			if (tmp[l3] != 0)
				l3++;
		}
	}

	while (l2 < A.value_counter[index1]) {
		tmp[l3] = A.value[index1][l2];
		if (tmp[l3] < A.Zero)
			tmp[l3] += mod;

		tmp2[l3] = A.index[index1][l2];
		l2++;
		l3++;
	}

	while (l1 < A.value_counter[index2]) {

		LiDIA::mult_mod(tmp[l3], q, A.value[index2][l1], mod);
		if (tmp[l3] < A.Zero)
			tmp[l3] += mod;

		tmp2[l3] = A.index[index2][l1];
		l1++;
		l3++;
	}

	delete[] A.value[index1];
	A.value[index1] = tmp;

	delete[] A.index[index1];
	A.index[index1] = tmp2;

	A.value_counter[index1] = l3;

	tmp = new T[l3];
	tmp2 = new lidia_size_t[l3];

	// memory adaption
	for (l1 = 0; l1 < l3; l1++) {
		tmp2[l1] = A.index[index1][l1];
		tmp[l1] = A.value[index1][l1];
	}

	delete[] A.value[index1];
	A.value[index1] = tmp;

	delete[] A.index[index1];
	A.index[index1] = tmp2;
}



template< class T >
void
column_oriented_sparse_matrix_modules< T >::negate_column (MR< T > &A,
							   lidia_size_t pos,
							   lidia_size_t len) const
{
	for (register lidia_size_t i = 0; i < A.value_counter[pos] && A.index[pos][i] <= len; i++)
		A.value[pos][i] = -A.value[pos][i];
}



template< class T >
void
column_oriented_sparse_matrix_modules< T >::negate_column_mod (MR< T > &A,
							       lidia_size_t pos,
							       lidia_size_t len,
							       const T &mod) const
{
	for (register lidia_size_t i = 0; i < A.value_counter[pos] && A.index[pos][i] <= len; i++) {
		A.value[pos][i] = -A.value[pos][i];
		LiDIA::best_remainder(A.value[pos][i], A.value[pos][i], mod);
	}
}



//
// hadamard bound
//

template< class T >
void
column_oriented_sparse_matrix_modules< T >::hadamard (const MR< T > &A,
						      bigint &H) const
{
	register lidia_size_t min, i, j;

	if (A.columns < A.rows) {
		min = A.columns;
	}
	else {
		min = A.rows;
	}

	bigint TMP;

	// Step 1 - 11
	register bigint *hcolumns = new bigint[A.columns];
	register bigint *hrows = new bigint[A.rows];

	for (j = 0; j < A.columns; j++) {
		for (i = 0; i < A.value_counter[j]; i++) {
			LiDIA::square(TMP, bigint(A.value[j][i]));
			LiDIA::add(hcolumns[j], hcolumns[j], TMP);
			LiDIA::add(hrows[A.index[j][i]], hrows[A.index[j][i]], TMP);
		}
	}

	// Step 12 - 15
	j = 0;
	if (A.rows < A.columns)
		while (j < A.columns - 1)
			if (hcolumns[j] < hcolumns[j + 1]) {
				LiDIA::swap(hcolumns[j], hcolumns[j + 1]);
				if (j > 0)
					j--;
			}
			else
				j++;


	// Step 16 - 20
	bigint COLUMNS = 1;
	for (i = 0; i < min; i++) {
		if (hcolumns[i] != 0) {
			LiDIA::multiply(COLUMNS, COLUMNS, hcolumns[i]);
		}
	}

	j = 0;
	if (A.rows > A.columns)
		while (j < A.rows - 1)
			if (hrows[j] < hrows[j + 1]) {
				LiDIA::swap(hrows[j], hrows[j + 1]);
				if (j > 0)
					j--;
			}				// h.sort(DOWN);
			else
				j++;


	// Step 35 - 39
	bigint ROWS;
	ROWS = 1;
	for (i = 0; i < min; i++)
		if (hrows[i] != 0) {
			LiDIA::multiply(ROWS, ROWS, hrows[i]);
		}

	delete[] hcolumns;
	delete[] hrows;

	// Step 40 - 45
	bigint M;
	if (COLUMNS < ROWS)
		M.assign(COLUMNS);
	else
		M.assign(ROWS);

	register lidia_size_t B = M.bit_length() - 1;
	bigint E = (bigint(B) / bigint(2)) + bigint(2);
	power(H, bigint(2), E);
	dec(H);
}



//
// Init and update
//

template< class T >
void
column_oriented_sparse_matrix_modules< T >::update_max_array (const MR< T > &A,
							      lidia_size_t i,
							      T *max_array) const
{
	for (register lidia_size_t j = 0; j < A.value_counter[i]; j++)
		if (max_array[i] < abs(A.value[i][j]))
			max_array[i] = abs(A.value[i][j]);
}



template< class T >
inline T*
column_oriented_sparse_matrix_modules< T >::init_max_array (const MR< T > &A) const
{
	T *max_array = new T[A.columns];
	memory_handler(max_array, DMESSAGE, "init :: "
		       "Error in memory allocation (max_array)");
	for (lidia_size_t i = 0; i < A.columns; i++)
		max_array[i] = 0;

	for (register lidia_size_t i = 0; i < A.columns; i++)
		for (register lidia_size_t j = 0; j < A.value_counter[i]; j++)
			if (max_array[i] < abs(A.value[i][j]))
				max_array[i] = abs(A.value[i][j]);

	return max_array;
}



//
// Pivot search
//

template< class T >
inline lidia_size_t
column_oriented_sparse_matrix_modules< T >::normal_pivot (const MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t &index) const
{
	lidia_size_t FOUND = 0;

	for (lidia_size_t i = 0; i <= startc; i++)
		if (member(A, startr, i) != A.Zero) {
			FOUND++;
			index = i;
			for (lidia_size_t j = 0; j <= startc; j++)
				if ((member(A, startr, j) != A.Zero) && (j != index)) {
					FOUND++;
					if (j != index) {
						if (abs(member(A, startr, j)) < abs(member(A, startr, i)))
							index = j;
						break;
					}
				}
			break;
		}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::min_abs_of_row (const MR< T > &A,
							    lidia_size_t startr,
							    lidia_size_t startc,
							    lidia_size_t &index) const
{
	T pivot = 0;
	register lidia_size_t COUNT = 0;
	lidia_size_t j;

	for (register lidia_size_t i = startc; i >= 0; i--) {
		register T *tmp = A.value[i];


		// Position
		j = A.value_counter[i] - 1;

		if (j >= 0 && A.index[i][j] > startr) {
			std::cout << "Indexfehler " << A.index[i][j] << std::endl;
			exit(1);
		}

		if (j >= 0 && A.index[i][j] == startr)
			if (COUNT == 0 || abs(pivot) > abs(tmp[j])) {
				pivot = tmp[j];
				index = i;
				COUNT++;
			}
			else {
				COUNT++;
				if (abs(pivot) == abs(tmp[j])) {
					if (A.value_counter[index] > A.value_counter[i])
						index = i;
				}
			}
	}
	return COUNT;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::pivot_sorting_gcd (const MR< T > &A,
							       lidia_size_t startr,
							       lidia_size_t startc,
							       lidia_size_t &index) const
{
	lidia_size_t FOUND = 0;
	lidia_size_t index_largest = 0;
	T val, maxval = 0;
	index = 0;

	for (lidia_size_t i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (FOUND == 0) {
				maxval = val;
				index_largest = i;
			}
			else if (maxval < val) {
				maxval = val;
				index_largest = i;
			}

			FOUND++;
		}
	}

	bool SW = false;

	if (FOUND > 1) {
		for (lidia_size_t i = 0; i <= startc; i++) {
			val = abs(member(A, startr, i));
			if (val != A.Zero && index_largest != i) {
				if (SW == false) {
					maxval = val;
					index = i;
				}
				else
					if (maxval <= val) {
						maxval = val;
						index = i;
					}
				SW = true;
			}
		}
	}
	else
		index = index_largest;

	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::minimal_norm (const MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t &index,
							  lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	bool Pair = false;
	T max_val = 0, val;

	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				max_val = val;
				SW = true;
				index = i;
			}
			else {
				if (val > max_val) {
					max_val = val;
					Pair = false;
				}
				else if (val == max_val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != max_val) {
					min_norm = column_norm(A, i, norm);
					index = i;
					SW = true;
				}
			}
			else {
				if (Pair) {
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_norm = nnorm;
						index = i;
					}
				}
				else {
					if (val != max_val) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_norm = nnorm;
							index = i;
						}
					}
				}
			}
		}
	}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::min_abs_of_row_plus_minimal_norm (const MR< T > &A,
									      lidia_size_t startr,
									      lidia_size_t startc,
									      lidia_size_t &index,
									      lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T minval = 0, val;
	bool SW = false;

	for (i = startc; i >= 0; i--) {
		val = abs(member(A, startr, i));

		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				index = i;
				minval = val;
				SW = true;
			}
			else {
				if (minval > val) {
					index = i;
					minval = val;
				}
			}
		}
	}

	if (FOUND > 0) {
	  min_norm = column_norm(A, index, norm);
	  
	  // norm computation
	  for (i = index - 1; i >= 0; i--) {
	    if (minval == abs(member(A, startr, i))) {
	      nnorm = column_norm(A, i, norm);
	      if (nnorm < min_norm) {
		index = i;
		min_norm = nnorm;
	      }
	    }
	  }
	}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::minimal_norm_plus_min_abs_of_row (const MR< T > &A,
									      lidia_size_t startr,
									      lidia_size_t startc,
									      lidia_size_t &index,
									      lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	bool Pair = false;
	T max_val = 0, val = 0, min_val = 0;

	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				max_val = val;
				SW = true;
				index = i;
			}
			else {
				if (val > max_val) {
					max_val = val;
					Pair = false;
				}
				else if (val == max_val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != max_val) {
					min_norm = column_norm(A, i, norm);
					index = i;
					min_val = val;
					SW = true;
				}
			}
			else {
				if (Pair) {
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_norm = nnorm;
						min_val = val;
						index = i;
					}
					else if (min_norm == nnorm)
						if (val < min_val) {
							min_norm = nnorm;
							min_val = val;
							index = i;
						}
				}
				else {
					if (val != max_val) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_norm = nnorm;
							min_val = val;
							index = i;
						}
						else if (min_norm == nnorm)
							if (val < min_val) {
								min_norm = nnorm;
								min_val = val;
								index = i;
							}
					}
				}
			}
		}
	}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::minimal_norm_plus_sorting_gcd (const MR< T > &A,
									   lidia_size_t startr,
									   lidia_size_t startc,
									   lidia_size_t &index,
									   lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	bool Pair = false;
	T max_val = 0, val = 0, min_val = 0;

	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				max_val = val;
				SW = true;
				index = i;
			}
			else {
				if (val > max_val) {
					max_val = val;
					Pair = false;
				}
				else if (val == max_val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != max_val) {
					min_norm = column_norm(A, i, norm);
					index = i;
					min_val = val;
					SW = true;
				}
			}
			else {
				if (Pair) {
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_norm = nnorm;
						min_val = val;
						index = i;
					}
					else if (min_norm == nnorm)
						if (val > min_val) {
							min_norm = nnorm;
							min_val = val;
							index = i;
						}
				}
				else {
					if (val != max_val) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_norm = nnorm;
							min_val = val;
							index = i;
						}
						else if (min_norm == nnorm)
							if (val > min_val) {
								min_norm = nnorm;
								min_val = val;
								index = i;
							}
					}
				}
			}
		}
	}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::minimal_norm_plus_min_no_of_elements (const MR< T > &A,
										  lidia_size_t startr,
										  lidia_size_t startc,
										  lidia_size_t &index,
										  lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint no, min_no;
	bigint min_norm, nnorm, TMP;
	bool Pair = false;
	T max_val = 0, val = 0;

	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				max_val = val;
				SW = true;
				index = i;
			}
			else {
				if (val > max_val) {
					max_val = val;
					Pair = false;
				}
				else if (val == max_val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != max_val) {
					min_norm = column_norm(A, i, norm);
					index = i;
					min_no = column_norm(A, i, 0);
					SW = true;
				}
			}
			else {
				if (Pair) {
					no = column_norm(A, i, 0);
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_norm = nnorm;
						min_no = no;
						index = i;
					}
					else if (min_norm == nnorm)
						if (no < min_no) {
							min_norm = nnorm;
							min_no = no;
							index = i;
						}
				}
				else {
					no = column_norm(A, i, 0);
					if (val != max_val) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_norm = nnorm;
							min_no = no;
							index = i;
						}
						else if (min_norm == nnorm)
							if (no < min_no) {
								min_norm = nnorm;
								min_no = no;
								index = i;
							}
					}
				}
			}
		}
	}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::sorting_gcd_plus_minimal_norm (const MR< T > &A,
									   lidia_size_t startr,
									   lidia_size_t startc,
									   lidia_size_t &index,
									   lidia_size_t norm) const
{
	lidia_size_t FOUND = 0;
	lidia_size_t index_largest = 0;
	T val, maxval = 0;
	index = 0;

	for (lidia_size_t i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (FOUND == 0) {
				maxval = val;
				index_largest = i;
			}
			else if (maxval < val) {
				maxval = val;
				index_largest = i;
			}

			FOUND++;
		}
	}

	bool SW = false;

	if (FOUND > 1) {
		for (lidia_size_t i = 0; i <= startc; i++) {
			val = abs(member(A, startr, i));
			if (val != A.Zero && index_largest != i) {
				if (SW == false) {
					maxval = val;
					index = i;
				}
				else
					if (maxval <= val) {
						maxval = val;
						index = i;
					}
				SW = true;
			}
		}
	}
	else {
		index = index_largest;
		maxval = abs(member(A, startr, index));
	}

	bigint nnorm, min_norm = column_norm(A, index, norm);

	for (lidia_size_t i = 0; i <= startc; i++)
		if (maxval == abs(member(A, startr, i))) {
			nnorm = column_norm(A, i, norm);
			if (nnorm < min_norm) {
				min_norm = nnorm;
				index = i;
			}
		}

	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_sparse_matrix_modules< T >::min_no_of_elements_plus_minimal_norm (const MR< T > &A,
										  lidia_size_t startr,
										  lidia_size_t startc,
										  lidia_size_t &index,
										  lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint no, min_no;
	bigint min_norm, nnorm, TMP;
	bool Pair = false;
	T max_val = 0, val = 0;

	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				max_val = val;
				SW = true;
				index = i;
			}
			else {
				if (val > max_val) {
					max_val = val;
					Pair = false;
				}
				else if (val == max_val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(member(A, startr, i));
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != max_val) {
					min_norm = column_norm(A, i, 0);
					index = i;
					min_no = column_norm(A, i, norm);
					SW = true;
				}
			}
			else {
				if (Pair) {
					no = column_norm(A, i, norm);
					nnorm = column_norm(A, i, 0);
					if (min_norm > nnorm) {
						min_norm = nnorm;
						min_no = no;
						index = i;
					}
					else if (min_norm == nnorm)
						if (no < min_no) {
							min_norm = nnorm;
							min_no = no;
							index = i;
						}
				}
				else {
					no = column_norm(A, i, norm);
					if (val != max_val) {
						nnorm = column_norm(A, i, 0);
						if (min_norm > nnorm) {
							min_norm = nnorm;
							min_no = no;
							index = i;
						}
						else if (min_norm == nnorm)
							if (no < min_no) {
								min_norm = nnorm;
								min_no = no;
								index = i;
							}
					}
				}
			}
		}
	}
	return FOUND;
}



//
// Norm
//

template< class T >
inline bigint
column_oriented_sparse_matrix_modules< T >::column_norm (const MR< T > &A,
							 lidia_size_t pos,
							 lidia_size_t norm) const
{
	lidia_size_t j;
	bigint Norm = A.Zero;
	T * tmp = A.value[pos];
	bigint TMP;

	if (norm == 0) {
		for (j = 0; j < A.value_counter[pos]; j++)
			if (tmp[j] != A.Zero)
				Norm++;
	}
	else if (norm == 1) {
		for (j = 0; j < A.value_counter[pos]; j++)
			LiDIA::add(Norm, Norm, tmp[j]);
	}
	else {
		for (j = 0; j < A.value_counter[pos]; j++) {
			power(TMP, tmp[j], norm);
			LiDIA::add(Norm, Norm, abs(TMP));
		}
	}
	return Norm;
}



template< class T >
void
column_oriented_sparse_matrix_modules< T >::kennwerte (MR< T > &RES,
						       T &MAX,
						       lidia_size_t &no_of_elements,
						       T &Durch) const
{
	register lidia_size_t i, j;
	bool SW = false;
	no_of_elements = 0;
	Durch = 0;

	for (i = 0; i < RES.columns; i++) {
		no_of_elements += RES.value_counter[i];
		for (j = 0; j < RES.value_counter[i]; j++) {
			Durch += abs(RES.value[i][j]);
			if (SW == false) {
				SW = true;
				MAX = abs(RES.value[i][j]);
			}
			else
				if (MAX < abs(RES.value[i][j]))
					MAX = abs(RES.value[i][j]);
		}
	}
	if (no_of_elements != 0)
		Durch /= no_of_elements;
	else
		Durch = 0;
}



template< class T >
inline void
column_oriented_sparse_matrix_modules< T >::max (MR< T > &RES, T &MAX) const
{
	register lidia_size_t i, j;
	bool SW = false;

	for (i = 0; i < RES.columns; i++)
		for (j = 0; j < RES.value_counter[i]; j++)
			if (SW == false) {
				SW = true;
				MAX = abs(RES.value[i][j]);
			}
			else
				if (MAX < abs(RES.value[i][j]))
					MAX = abs(RES.value[i][j]);
}



//
// Hermite normal form: elimination
//

template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));
			if (q != 0)
				subtract_multiple_of_column(A, i, q, index, len);
			break;
		}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
									  matrix< bigint > &TR,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));
			if (q != 0) {
				subtract_multiple_of_column(A, i, q, index, len);
				// => subtract_multiple_of_column(TR, i, q, index, TR.get_no_of_rows() - 1);
			}
			break;
		}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
									  math_vector< bigfloat > &vec,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	T q, TMP;
	bigfloat rtemp;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));
			if (q != 0) {
				subtract_multiple_of_column(A, i, q, index, len);
				multiply(rtemp, bigfloat(q), vec.member(index));
				subtract(vec[i], vec[i], rtemp);
			}
			break;
		}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_linear (MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;

	// Init
	lidia_size_t index;
	for (index = startc; index >= 0 && member(A, startr, index) == A.Zero; index--);

	if (index < 0)
		return false;

	if (index != startc)
		swap_columns(A, startc, index);

	for (lidia_size_t i = index - 1; i >= 0; i--) {
		if (member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, startc));
			x = member(A, startr, startc) / RES2;
			y = member(A, startr, i) / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, startc);

				sto(A, l, i, (TMP*x - TMP1*y));
				sto(A, l, startc, (TMP*RES0 + TMP1*RES1));
			}
		}
	}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_linear (MR< T > &A,
							 matrix< bigint > &TR,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	bigint TMP_1, TMP_2;

	// Init
	lidia_size_t index;
	for (index = startc; index >= 0 && member(A, startr, index) == A.Zero; index--);
	if (index < 0)
		return false;
	if (index != startc) {
		swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	for (lidia_size_t i = index - 1; i >= 0; i--)
		if (member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, startc));
			x = member(A, startr, startc) / RES2;
			y = member(A, startr, i) / RES2;

			lidia_size_t l;
			for (l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, startc);

				sto(A, l, i, (TMP*x - TMP1*y));
				sto(A, l, startc, (TMP*RES0 + TMP1*RES1));
			}

			for (l = 0; l < TR.get_no_of_rows(); l++) {
				TMP_1 = TR.member(l, i);
				TMP_2 = TR.member(l, startc);

				TR.sto(l, i, (TMP_1*x - TMP_2*y));
				TR.sto(l, startc, (TMP_1*RES0 + TMP_2*RES1));
			}
		}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_bradley (MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t len) const
{
	lidia_size_t n = startc + 1;
	register lidia_size_t i, j;

	matrix< T > Tr;
	Tr.set_representation(matrix_flags::sparse_representation);
	Tr.resize(A.columns, A.columns);

	Tr.diag(1, 0);

	T *a = new T[n];

	for (i = 0; i < n; i++)
		a[i] = member(A, startr, i);

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		T y, z, TMP, TMP1, t1, t2;
		T g_old = a[0];
		T g = g_old;
		for (i = 1; i < n; i++) {
			if (a[i] != 0) {
				g = xgcd(y, z, g_old, a[i]);
				TMP = -a[i]/g;
				TMP1 = g_old/g;
				for (j = 0; j < n; j++) {
					t1 = Tr(j, i-1);
					t2 = Tr(j, i);
					Tr.sto(j, i-1, TMP*t1+TMP1*t2);
					Tr.sto(j, i, y*t1+z*t2);
				}
				g_old = g;
			}
			else {
				Tr.swap_columns(i-1, i);
			}
		}
		if (g == 0)
			return false;

		// Balaciertes Rueckwaertseinsetzen
		T q, r;
		lidia_size_t l, k;
		for (i = n - 1, j = n - 2; i >= 0 && j >= 0; i--, j--) {
			for (k = j + 1; k < n; k++) {
				if (Tr(i, j) != 0) {
					pos_div_rem(q, r, Tr(i, k), Tr(i, j));

					Tr.sto(i, k, Tr(i, k) - q * Tr(i, j));
					TMP = Tr(i, k) + Tr(i, j);
					TMP1 = Tr(i, k) - Tr(i, j);
					if (abs(TMP) < abs(Tr(i, k))) {
						q--;
						Tr.sto(i, k, TMP);
					}
					else if (abs(TMP1) < abs(Tr(i, k))) {
						q++;
						Tr.sto(i, k, TMP1);
					}

					if (abs(q) > 0)
						for (l = 0; l < i; l++)
							Tr.sto(l, k, Tr(l, k) - q*Tr(l, j));
				}
			}
		}
	}
	//multiply(A, A, Tr);
	T TMP;
	T *tmp2 = new T[A.columns];
	std::cout << "Vor multiply" << std::endl;
	for (i = 0; i < A.rows; i++) {
		std::cout << i << std::endl;
		for (lidia_size_t p = 0; p < A.columns; p++)
			tmp2[p] = member(A, i, p);
		for (j = 0; j < A.columns; j++) {
			TMP = A.Zero;
			for (lidia_size_t l = 0; l < A.columns; l++) {
				if (tmp2[l] != A.Zero && Tr.member(l, j) != A.Zero)
					TMP += T(tmp2[l] * Tr.member(l, j));
			}
			sto(A, i, j, TMP);
		}
	}
	delete[] tmp2;

	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_bradley (MR< T > &A,
							  matrix< bigint > &TR,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t len) const
{
	lidia_size_t n = startc + 1;
	register lidia_size_t i, j;


	matrix< T > Tr;
	Tr.set_representation(matrix_flags::sparse_representation);
	Tr.resize(A.columns, A.columns);

	Tr.diag(1, 0);

	T *a = new T[n];

	for (i = 0; i < n; i++)
		a[i] = member(A, startr, i);

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		T y, z, TMP, TMP1, t1, t2;
		T g_old = a[0];
		T g = g_old;
		for (i = 1; i < n; i++) {
			if (a[i] != 0) {
				g = xgcd(y, z, g_old, a[i]);
				TMP = -a[i]/g;
				TMP1 = g_old/g;
				for (j = 0; j < n; j++) {
					t1 = Tr(j, i-1);
					t2 = Tr(j, i);
					Tr.sto(j, i-1, TMP*t1+TMP1*t2);
					Tr.sto(j, i, y*t1+z*t2);
				}
				g_old = g;
			}
			else {
				Tr.swap_columns(i-1, i);
			}
		}

		if (g == 0)
			return false;

		// Balaciertes Rueckwaertseinsetzen
		T q, r;
		lidia_size_t l, k;
		for (i = n - 1, j = n - 2; i >= 0 && j >= 0; i--, j--)
			for (k = j + 1; k < n; k++) {
				if (Tr(i, j) != 0) {
					pos_div_rem(q, r, Tr(i, k), Tr(i, j));

					Tr.sto(i, k, Tr(i, k) - q*Tr(i, j));
					TMP = Tr(i, k) + Tr(i, j);
					TMP1 = Tr(i, k) - Tr(i, j);
					if (abs(TMP) < abs(Tr(i, k))) {
						q--;
						Tr.sto(i, k, TMP);
					}
					else if (abs(TMP1) < abs(Tr(i, k))) {
						q++;
						Tr.sto(i, k, TMP1);
					}

					if (abs(q) > 0)
						for (l = 0; l < i; l++)
							Tr.sto(l, k, Tr(l, k) - q*Tr(l, j));
				}
			}
	}
	//multiply(A, A, Tr);
	T TMP;
	T *tmp2 = new T[A.columns];
	for (i = 0; i < A.rows; i++) {
		for (lidia_size_t p = 0; p < A.columns; p++)
			tmp2[p] = member(A, i, p);
		for (j = 0; j < A.columns; j++) {
			TMP = A.Zero;
			for (lidia_size_t l = 0; l < A.columns; l++)
				if (tmp2[l] != A.Zero && Tr.member(l, j) != A.Zero)
					TMP += T(tmp2[l] * Tr.member(l, j));
			sto(A, i, j, TMP);
		}
	}
	delete[] tmp2;

	LiDIA::multiply(TR, TR, (matrix< bigint > )Tr);
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_ilio (MR< T > &A,
						       lidia_size_t startr,
						       lidia_size_t startc,
						       lidia_size_t len) const
{
	lidia_size_t n = startc + 1;

	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	T g = A.Zero;
	if (n > 1) {
		T y, z, TMP, TMP1, t1, t2;

		lidia_size_t step = 1, step2 = 2;
		while (step < n) {
			for (i = 0; i < n; i += step2) {
				if (i + step < n) {
					g = xgcd(y, z, member(A, startr, i), member(A, startr, i+step));
					if (g != A.Zero) {
						TMP = -member(A, startr, i+step)/g;
						TMP1 = member(A, startr, i)/g;

						for (j = 0; j < A.rows; j++) {

							t1 = member(A, j, i);
							t2 = member(A, j, i+step);
							sto(A, j, i, y*t1+z*t2);
							sto(A, j, i+step, TMP*t1+TMP1*t2);
						}
					}
				}
			}
			step = step2;
			step2 = 2*step2;
		}
		swap_columns(A, 0, n-1);
		if (member(A, startr, startc) < 0)
		  for (i = 0; i <= startr; i++)
		    sto(A, i, startc, -member(A, i, startc));
	}
	if (g == A.Zero)
		return false;
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_ilio (MR< T > &A,
						       matrix< bigint > &Tr,
						       lidia_size_t startr,
						       lidia_size_t startc,
						       lidia_size_t len) const
{
	lidia_size_t n = startc + 1;

	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	T g = A.Zero;
	if (n > 1) {
		T y, z, TMP, TMP1, t3, t4;
		bigint t1, t2;

		lidia_size_t step = 1, step2 = 2;
		while (step < n) {
			for (i = 0; i < n; i += step2) {
				if (i + step < n) {
					g = xgcd(y, z, member(A, startr, i), member(A, startr, i+step));
					if (g != A.Zero) {
						TMP = -member(A, startr, i+step)/g;
						TMP1 = member(A, startr, i)/g;
						for (j = 0; j < A.columns; j++) {
							t1 = Tr.member(j, i);
							t2 = Tr.member(j, i + step);
							Tr.sto(j, i, y*t1+z*t2);
							Tr.sto(j, i + step, TMP*t1+TMP1*t2);
						}

						for (j = 0; j < A.rows; j++) {
							t3 = member(A, j, i);
							t4 = member(A, j, i+step);
							sto(A, j, i, y*t3+z*t4);
							sto(A, j, i+step, TMP*t3+TMP1*t4);
						}
					}
				}
			}
			step = step2;
			step2 = 2*step2;
		}

		swap_columns(A, 0, n-1);
		Tr.swap_columns(0, n-1);
		if (member(A, startr, startc) < 0)
		  {
		    for (i = 0; i <= startr; i++)
		      sto(A, i, startc, -member(A, i, startc));
		    for (i = 0; i < A.columns; i++)
		      Tr.sto(i, startc, -Tr.member(i, startc));
		  }
	}
	if (g == A.Zero)
		return false;
	return true;
}



template< class T >
bool column_oriented_sparse_matrix_modules< T >::mgcd_opt (MR< T > &A,
								  lidia_size_t startr,
								  lidia_size_t startc,
								  lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;

	// Init
	lidia_size_t index;
	for (index = startc; member(A, startr, index) == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc)
		swap_columns(A, startc, index);

	// sort
	if (startc >= 2) {
		lidia_size_t gcdindex = startc - 1;
		T W = member(A, startr, startc);
		T gcd_min = gcd(W, member(A, startr, startc - 1)), gcd_tmp;
		for (lidia_size_t i = startc - 2; i >= 0; i--) {
			gcd_tmp = gcd(W, member(A, startr, i));
			if (gcd_tmp < gcd_min)
				gcdindex = i;
		}
		swap_columns(A, gcdindex, startc - 1);
	}

	for (lidia_size_t i = startc - 1; i >= 0; i--)
		if (member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, startc));
			x = member(A, startr, startc) / RES2;
			y = member(A, startr, i) / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, startc);

				sto(A, l, i, (TMP*x - TMP1*y));
				sto(A, l, startc, (TMP*RES0 + TMP1*RES1));
			}
		}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_opt (MR< T > &A,
						      matrix< bigint > &TR,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	bigint TMP_1, TMP_2;

	// Init
	lidia_size_t index;
	for (index = startc; index >= 0 && member(A, startr, index) == A.Zero; index--);
	if (index < 0)
		return false;
	if (index != startc) {
		swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	// sort
	if (startc >= 2) {
		lidia_size_t gcdindex = startc - 1;
		T W = member(A, startr, startc);
		T gcd_min = gcd(W, member(A, startr, startc - 1)), gcd_tmp;
		for (lidia_size_t i = startc - 2; i >= 0; i--) {
			gcd_tmp = gcd(W, member(A, startr, i));
			if (gcd_tmp < gcd_min)
				gcdindex = i;
		}
		swap_columns(A, gcdindex, startc - 1);
		TR.swap_columns(gcdindex, startc - 1);
	}

	for (lidia_size_t i = startc - 1; i >= 0; i--)
		if (member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, startc));
			x = member(A, startr, startc) / RES2;
			y = member(A, startr, i) / RES2;

			lidia_size_t l;
			for (l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, startc);

				sto(A, l, i, (TMP*x - TMP1*y));
				sto(A, l, startc, (TMP*RES0 + TMP1*RES1));
			}

			for (l = 0; l < TR.get_no_of_rows(); l++) {
				TMP_1 = TR.member(l, i);
				TMP_2 = TR.member(l, startc);

				TR.sto(l, i, (TMP_1*x - TMP_2*y));
				TR.sto(l, startc, (TMP_1*RES0 + TMP_2*RES1));
			}
		}
	return true;
}



template< class T >
bool column_oriented_sparse_matrix_modules< T >::mgcd_storjohann (MR< T > &A,
									 lidia_size_t startr,
									 lidia_size_t startc,
									 lidia_size_t len) const
{
	T TMP;
	T g = 0, N = 0, a = 0, qa = 0, ra = 0, qb = 0, rb = 0, t = 0;
	T A3 = 0, A2 = 0, A1 = 0, A0 = 0, x = 0, y = 0;
	T TMP1, TMP2, TMP3;

	register lidia_size_t i, j, l;
	lidia_size_t index;

	//
	// init
	//

	for (index = startc; index >= 0 && member(A, startr, index) == A.Zero; index--);

	if (index < 0)
		return false;
	else {
		//
		// conditioning
		//

		N = member(A, startr, index);
		for (j = index - 1; j >= 0 && member(A, startr, j) == A.Zero; j--);
		if (j >= 0) {
			a = member(A, startr, j);
			t = 0;

			for (i = j - 1; i >= 0; i--) {
				if (member(A, startr, i) != 0) {
					g = gcd(member(A, startr, i), a);

					pos_div_rem(qa, ra, a/g, N);
					pos_div_rem(qb, rb, member(A, startr, i)/g, N);

					for (t = 0; gcd((ra + t*rb), N) != 1; t++);

					a = a + t * member(A, startr, i); // NEW

					//M.sto(1, i, t);
					for (l = 0; l <= startr; l++)
						sto(A, l, j, member(A, l, j) + t * member(A, l, i));
				}
			}

			//
			// elimination
			//

			A2 = xgcd(A0, A1, member(A, startr, j), member(A, startr, index));
			pos_div_rem(x, A3, member(A, startr, index), A2);
			pos_div_rem(y, A3, member(A, startr, j), A2);

			for (l = 0; l <= startr; l++) {
				TMP = member(A, l, j);
				TMP1 = member(A, l, index);


				// Atmp1[j] = ((TMP * x) + (TMP1 * y))
				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				sto(A, l, j, TMP2 - TMP3);

				// Atmp1[n-m+i] = ((TMP * A0) + (TMP1 * A1))
				LiDIA::multiply(TMP2, TMP, A0);
				LiDIA::multiply(TMP3, TMP1, A1);
				sto(A, l, index, TMP2 + TMP3);
			}

			for (i = j - 1; i >= 0; i--)
				if (member(A, startr, i) != A.Zero) {
					pos_div_rem(TMP, TMP1, member(A, startr, i), member(A, startr, index));
					for (l = 0; l <= startr; l++)
						sto(A, l, i, member(A, l, i) - TMP * member(A, l, index));
				}

		}

		if (member(A, startr, index) < 0)
			for (i = 0; i <= startr; i++)
				sto(A, i, index, -member(A, i, index));
		if (index != startc)
			swap_columns(A, startc, index);
	}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::mgcd_storjohann (MR< T > &A,
							     matrix< bigint > &TR,
							     lidia_size_t startr,
							     lidia_size_t startc,
							     lidia_size_t len) const
{
	T TMP;
	T g, N, a, qa, ra, qb, rb, t;
	T A3, A2, A1, A0, x, y;
	T TMP1, TMP2, TMP3;

	register lidia_size_t i, j, l;
	lidia_size_t index;

	//
	// init
	//

	for (index = startc; index >= 0 && member(A, startr, index) == A.Zero; index--);

	if (index < 0)
		return false;
	else {
		//
		// conditioning
		//

		N = member(A, startr, index);
		for (j = index - 1; j >= 0 && member(A, startr, j) == A.Zero; j--);
		if (j >= 0) {
			a = member(A, startr, j);
			t = 0;

			for (i = j - 1; i >= 0; i--) {
				if (member(A, startr, i) != 0) {
					g = gcd(member(A, startr, i), a);

					pos_div_rem(qa, ra, a/g, N);
					pos_div_rem(qb, rb, member(A, startr, i)/g, N);

					for (t = 0; gcd((ra + t*rb), N) != 1; t++);

					a = a + t * member(A, startr, i); // NEW

					//M.sto(1, i, t);
					for (l = 0; l <= startr; l++)
						sto(A, l, j, member(A, l, j) + t * member(A, l, i));
					for (l = 0; l < A.columns; l++)
						TR.sto(l, j, TR.member(l, j) + t * TR.member(l, i));
				}
			}

			//
			// elimination
			//

			A2 = xgcd(A0, A1, member(A, startr, j), member(A, startr, index));
			pos_div_rem(x, A3, member(A, startr, index), A2);
			pos_div_rem(y, A3, member(A, startr, j), A2);

			for (l = 0; l <= startr; l++) {
				TMP = member(A, l, j);
				TMP1 = member(A, l, index);


				// Atmp1[j] = ((TMP * x) + (TMP1 * y))
				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				sto(A, l, j, TMP2 - TMP3);

				// Atmp1[n-m+i] = ((TMP * A0) + (TMP1 * A1))
				LiDIA::multiply(TMP2, TMP, A0);
				LiDIA::multiply(TMP3, TMP1, A1);
				sto(A, l, index, TMP2 + TMP3);
			}

			bigint TTMP, TTMP1, TTMP2, TTMP3;
			for (l = 0; l < A.columns; l++) {
				TTMP = TR.member(l, j);
				TTMP1 = TR.member(l, index);


				// Atmp1[j] = ((TMP * x) + (TMP1 * y))
				LiDIA::multiply(TTMP2, TTMP, x);
				LiDIA::multiply(TTMP3, TTMP1, y);
				TR.sto(l, j, TTMP2 - TTMP3);

				// Atmp1[n-m+i] = ((TMP * A0) + (TMP1 * A1))
				LiDIA::multiply(TTMP2, TTMP, A0);
				LiDIA::multiply(TTMP3, TTMP1, A1);
				TR.sto(l, index, TTMP2 + TTMP3);
			}

			for (i = j - 1; i >= 0; i--)
				if (member(A, startr, i) != A.Zero) {
					pos_div_rem(TMP, TMP1, member(A, startr, i), member(A, startr, index));
					for (l = 0; l <= startr; l++)
						sto(A, l, i, member(A, l, i) - TMP * member(A, l, index));
					for (l = 0; l < A.columns; l++)
						TR.sto(l, i, TR.member(l, i) - TMP * TR.member(l, index));
				}

		}

		if (member(A, startr, index) < 0) {
			for (i = 0; i <= startr; i++)
				sto(A, i, index, -member(A, i, index));
			for (i = 0; i < A.columns; i++)
				TR.sto(i, index, -TR.member(i, index));
		}
		if (index != startc) {
			swap_columns(A, startc, index);
			TR.swap_columns(startc, index);
		}
	}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::xgcd_elimination (MR< T > &A,
							      lidia_size_t startr,
							      lidia_size_t startc,
							      lidia_size_t index,
							      lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	T TMP2, TMP3;

	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, index));
			x = member(A, startr, index) / RES2;
			y = member(A, startr, i) / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, index);

				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				LiDIA::subtract(TMP2, TMP2, TMP3);
				sto(A, l, i, TMP2);

				LiDIA::multiply(TMP2, TMP, RES0);
				LiDIA::multiply(TMP3, TMP1, RES1);
				LiDIA::add(TMP2, TMP2, TMP3);
				sto(A, l, index, TMP2);
			}
		}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::xgcd_elimination (MR< T > &A,
							      matrix< bigint > &TR,
							      lidia_size_t startr,
							      lidia_size_t startc,
							      lidia_size_t index,
							      lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	T TMP2, TMP3;

	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, index));
			x = member(A, startr, index) / RES2;
			y = member(A, startr, i) / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, index);

				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				LiDIA::subtract(TMP2, TMP2, TMP3);
				sto(A, l, i, TMP2);

				LiDIA::multiply(TMP2, TMP, RES0);
				LiDIA::multiply(TMP3, TMP1, RES1);
				LiDIA::add(TMP2, TMP2, TMP3);
				sto(A, l, index, TMP2);
			}

			bigint TTMP, TTMP1, TTMP2, TTMP3, TTMP4;
			for (lidia_size_t l = 0; l < A.columns; l++) {
				TTMP = TR.member(l, i);
				TTMP1 = TR.member(l, index);

				LiDIA::multiply(TTMP2, TTMP, x);
				LiDIA::multiply(TTMP3, TTMP1, y);
				LiDIA::subtract(TTMP2, TTMP2, TTMP3);
				TR.sto(l, i, TTMP2);

				LiDIA::multiply(TTMP2, TTMP, RES0);
				LiDIA::multiply(TTMP3, TTMP1, RES1);
				LiDIA::add(TTMP2, TTMP2, TTMP3);
				TR.sto(l, index, TTMP2);
			}
		}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::xgcd_elimination_mod (MR< T > &A,
								  lidia_size_t startr,
								  lidia_size_t startc,
								  lidia_size_t index,
								  lidia_size_t len,
								  const T &mod) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	T TMP2, TMP3;

	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			RES2 = xgcd(RES0, RES1, member(A, startr, i), member(A, startr, index));
			x = member(A, startr, index) / RES2;
			y = member(A, startr, i) / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = member(A, l, i);
				TMP1 = member(A, l, index);

				LiDIA::mult_mod(TMP2, TMP, x, mod);
				LiDIA::mult_mod(TMP3, TMP1, y, mod);
				LiDIA::sub_mod(TMP2, TMP2, TMP3, mod);
				sto(A, l, i, TMP2);

				LiDIA::mult_mod(TMP2, TMP, RES0, mod);
				LiDIA::mult_mod(TMP3, TMP1, RES1, mod);
				LiDIA::add_mod(TMP2, TMP2, TMP3, mod);
				sto(A, l, index, TMP2);
			}
		}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_row (MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len) const
{
	T q, TMP;
	if (member(A, startr, index) != A.Zero)
		for (lidia_size_t i = 0; i <= startc; i++)
			if ((i != index) && member(A, startr, i) != A.Zero) {
				pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));
				if (q != 0)
					subtract_multiple_of_column(A, i, q, index, len);
			}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_row_mod (MR< T > &A,
							       lidia_size_t startr,
							       lidia_size_t startc,
							       lidia_size_t index,
							       lidia_size_t len,
							       const T &mod) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));
			if (q != 0)
				subtract_multiple_of_column_mod(A, i, q, index, len, mod);
		}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_row (MR< T > &A,
							   matrix< bigint > &TR,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len) const
{
	T q, TMP;
	if (member(A, startr, index) != A.Zero)
		for (lidia_size_t i = 0; i <= startc; i++)
			if ((i != index) && member(A, startr, i) != A.Zero) {
				pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));

				// Check
				if (q != 0) {
					subtract_multiple_of_column(A, i, q, index, len);

					for (lidia_size_t j = 0; j < A.columns; j++)
						TR.sto(j, i, TR.member(j, i) - q*TR.member(j, index));
				}
			}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_row (MR< T > &A,
							   math_vector< bigfloat > &vec,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len) const
{
	T q, TMP;
	bigfloat rtemp;
	if (member(A, startr, index) != A.Zero)
		for (lidia_size_t i = 0; i <= startc; i++)
			if ((i != index) && member(A, startr, i) != A.Zero) {
				pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));

				// Check
				if (q != 0) {
					subtract_multiple_of_column(A, i, q, index, len);
					multiply(rtemp, bigfloat(q), vec.member(index));
					subtract(vec[i], vec[i], rtemp);
				}
			}
	return true;
}



template< class T >
inline bool
column_oriented_sparse_matrix_modules< T >::normalize_row (MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len,
							   T *max_array,
							   const T &BOUND) const
{
	T q, TMP;
	if (member(A, startr, index) != A.Zero)
		for (lidia_size_t i = 0; i <= startc; i++)
			if ((i != index) && member(A, startr, i) != A.Zero) {
				pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));

				// Check
				if (q != 0) {
					if (abs(bigint(bigint(max_array[index]) * bigint(q)) + bigint(max_array[i])) > bigint(BOUND))
						return false;

					subtract_multiple_of_column(A, i, q, index, len);
					update_max_array(A, i, max_array);
				}
			}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::normalize_row (MR< T > &A,
							   matrix< bigint > &TR,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len,
							   T *max_array,
							   const T &BOUND) const
{
	T q, TMP;
	if (member(A, startr, index) != A.Zero)
		for (lidia_size_t i = 0; i <= startc; i++)
			if ((i != index) && member(A, startr, i) != A.Zero) {
				pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));

				// Check
				if (q != 0) {
					if (abs(bigint(bigint(max_array[index]) * q) + max_array[i]) > bigint(BOUND))
						return false;

					subtract_multiple_of_column(A, i, q, index, len);

					for (lidia_size_t j = 0; j < startr; j++)
						TR.sto(j, i, TR.member(j, i) - q * TR.member(j, index));

					update_max_array(A, i, max_array);
				}
			}
	return true;
}



template< class T >
bool
column_oriented_sparse_matrix_modules< T >::normalize_row (MR< T > &A,
							   math_vector< bigfloat > &vec,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len,
							   T *max_array,
							   const T &BOUND) const
{
	T q, TMP;
	bigfloat rtemp;
	if (member(A, startr, index) != A.Zero)
		for (lidia_size_t i = 0; i <= startc; i++)
			if ((i != index) && member(A, startr, i) != A.Zero) {
				pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));

				// Check
				if (q != 0) {
					if (abs(bigint(bigint(max_array[index]) * q) + max_array[i]) > bigint(BOUND))
						return false;
					subtract_multiple_of_column(A, i, q, index, len);
					multiply(rtemp, bigfloat(q), vec.member(index));
					subtract(vec[i], vec[i], rtemp);
					update_max_array(A, i, max_array);
				}
			}
	return true;
}



#undef column_oriented_sparse_matrix_modules



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_BIGINT_MATRIX_MODULES_CC_GUARD_
