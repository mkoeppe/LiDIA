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


#ifndef LIDIA_DENSE_BIGINT_MATRIX_MODULES_CC_GUARD_
#define LIDIA_DENSE_BIGINT_MATRIX_MODULES_CC_GUARD_


#include	"LiDIA/modular_operations.inl"



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define row_oriented_dense_matrix_modules RODMM
#define column_oriented_dense_matrix_modules CODMM



//
// max
//

template< class T >
inline void
row_oriented_dense_matrix_modules< T >::max_abs (const MR< T > &RES, T &MAX) const
{
	MAX.assign(abs(RES.value[0][0]));

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MAX < abs(tmp[j]))
				MAX.assign(abs(tmp[j]));
	}
}



//
// special functions
//

template< class T >
inline void
row_oriented_dense_matrix_modules< T >::subtract_multiple_of_column (MR< T > &A,
								     lidia_size_t l1,
								     const T &q,
								     lidia_size_t l2,
								     lidia_size_t l) const
{
	for (register lidia_size_t i = 0; i <= l; i++)
		LiDIA::subtract(A.value[i][l1], A.value[i][l1], q*A.value[i][l2]);
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::subtract_multiple_of_column_mod (MR< T > &A,
									 lidia_size_t l1,
									 const T &q,
									 lidia_size_t l2,
									 lidia_size_t l,
									 const T &mod) const
{
	T TMP;
	for (register lidia_size_t i = 0; i <= l; i++) {
		LiDIA::mult_mod(TMP, q, A.value[i][l2], mod);
		LiDIA::sub_mod(A.value[i][l1], A.value[i][l1], TMP, mod);
	}
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::normalize_column_mod (MR< T > &A,
							      lidia_size_t l1,
							      const T &q,
							      lidia_size_t l2,
							      lidia_size_t l,
							      const T &mod) const
{
	T TMP;
	for (register lidia_size_t i = 0; i <= l; i++) {
		LiDIA::mult_mod(TMP, q, A.value[i][l2], mod);
		LiDIA::sub_mod(A.value[i][l1], A.value[i][l1], TMP, mod);
		if (A.value[i][l1] < A.Zero)
			A.value[i][l1] += mod;
	}
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::negate_column (MR< T > &A,
						       lidia_size_t index,
						       lidia_size_t l) const
{
	for (register lidia_size_t i = 0; i <= l; i++)
		A.value[i][index] = -A.value[i][index];
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::negate_column_mod (MR< T > &A,
							   lidia_size_t index,
							   lidia_size_t l,
							   const T &mod) const
{
	for (register lidia_size_t i = 0; i <= l; i++) {
		A.value[i][index] = -A.value[i][index];
		LiDIA::best_remainder(A.value[i][index], A.value[i][index], mod);
	}
}



template< class T >
inline T *
row_oriented_dense_matrix_modules< T >::init_max_array (const MR< T > &A) const
{
	T *max_array = new T[A.columns];
	memory_handler(max_array, DMESSAGE, "init :: "
		       "Error in memory allocation (max_array)");
	T *tmp;
	for (register lidia_size_t i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		for (register lidia_size_t j = 0; j < A.columns; j++)
			if (max_array[j] < abs(tmp[j]))
				max_array[j] = abs(tmp[j]);
	}
	return max_array;
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::update_max_array (const MR< T > &A,
							  lidia_size_t i,
							  T *max_array) const
{
	for (register lidia_size_t j = 0; j < A.rows; j++)
		if (max_array[i] < abs(A.value[j][i]))
			max_array[i] = abs(A.value[j][i]);
}



//
// Hermite normal form: pivot search
//

template< class T >
inline lidia_size_t
row_oriented_dense_matrix_modules< T >::normal_pivot (const MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t &index) const
{
	lidia_size_t FOUND = 0;

	T *tmp = A.value[startr];

	for (lidia_size_t i = 0; i <= startc; i++)
		if (tmp[i] != A.Zero) {
			FOUND++;
			index = i;
			for (lidia_size_t j = i + 1; j <= startc; j++)
				if (tmp[j] != A.Zero && (j != index)) {
					FOUND++;
					if (j != index) {
						if (abs(tmp[j]) < abs(tmp[i]))
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
row_oriented_dense_matrix_modules< T >::min_abs_of_row (const MR< T > &A,
							lidia_size_t startr,
							lidia_size_t startc,
							lidia_size_t &pos) const
{
	lidia_size_t i, FOUND = 0;
	T *tmp = A.value[startr];
	T val, minval;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;

			if (SW == false) {
				minval = val;
				pos = i;
				SW = true;
			}
			else
				if (minval > val) {
					minval = val;
					pos = i;
				}
		}
	}
	return FOUND;
}



template< class T >
lidia_size_t
row_oriented_dense_matrix_modules< T >::pivot_sorting_gcd (const MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t &index) const
{
	lidia_size_t FOUND = 0;
	lidia_size_t index_largest = 0;
	T *tmp = A.value[startr];
	T maxval, val;

	for (lidia_size_t i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (FOUND == 0) {
				maxval = val;
				index_largest = i;
			}
			else
				if (maxval < val) {
					maxval = val;
					index_largest = i;
				}
			FOUND++;
		}
	}

	bool SW = false;

	if (FOUND > 1) {
		for (lidia_size_t i = 0; i <= startc; i++) {
			val = abs(tmp[i]);
			if (val != A.Zero && index_largest != i) {
				if (SW == false) {
					maxval = val;
					index = i;
					SW = true;
				}
				else
					if (maxval <= val) {
						maxval = val;
						index = i;
					}
			}
		}
	}
	else
		index = index_largest;

	return FOUND;
}



template< class T >
lidia_size_t
row_oriented_dense_matrix_modules< T >::minimal_norm (MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t &index,
						      lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T val, maxval;

	T *tmp = A.value[startr];
	bool Pair = false;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (val == maxval)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != maxval) {
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
					if (val != maxval) {
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
row_oriented_dense_matrix_modules< T >::min_abs_of_row_plus_minimal_norm (MR< T > &A,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t &index,
									  lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T minval = 0, val;
	T *tmp = A.value[startr];
	bool SW = false;

	for (i = startc; i >= 0; i--) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;

			if (SW == false) {
				minval = val;
				index = i;
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
	    if (minval == abs(tmp[i])) {
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
row_oriented_dense_matrix_modules< T >::minimal_norm_plus_min_abs_of_row (MR< T > &A,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t &index,
									  lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T val, maxval, min_val;

	T *tmp = A.value[startr];
	bool Pair = false;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (val == maxval)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != maxval) {
					min_val = val;
					min_norm = column_norm(A, i, norm);
					index = i;
					SW = true;
				}
			}
			else {
				if (Pair) {
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_val = val;
						min_norm = nnorm;
						index = i;
					}
					else if (min_norm == nnorm)
						if (min_val > val) {
							min_val = val;
							min_norm = nnorm;
							index = i;
						}
				}
				else {
					if (val != maxval) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_val = val;
							min_norm = nnorm;
							index = i;
						}
						else if (min_norm == nnorm)
							if (min_val > val) {
								min_val = val;
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
row_oriented_dense_matrix_modules< T >::minimal_norm_plus_sorting_gcd (MR< T > &A,
								       lidia_size_t startr,
								       lidia_size_t startc,
								       lidia_size_t &index,
								       lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T val, maxval, min_val;

	T *tmp = A.value[startr];
	bool Pair = false;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (val == maxval)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != maxval) {
					min_val = val;
					min_norm = column_norm(A, i, norm);
					index = i;
					SW = true;
				}
			}
			else {
				if (Pair) {
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_val = val;
						min_norm = nnorm;
						index = i;
					}
					else if (min_norm == nnorm)
						if (min_val < val) {
							min_val = val;
							min_norm = nnorm;
							index = i;
						}
				}
				else {
					if (val != maxval) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_val = val;
							min_norm = nnorm;
							index = i;
						}
						else if (min_norm == nnorm)
							if (min_val < val) {
								min_val = val;
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
row_oriented_dense_matrix_modules< T >::minimal_norm_plus_min_no_of_elements (MR< T > &A,
									      lidia_size_t startr,
									      lidia_size_t startc,
									      lidia_size_t &index,
									      lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_no, no;
	bigint min_norm, nnorm, TMP;
	T val, maxval;

	T *tmp = A.value[startr];
	bool Pair = false;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (val == maxval)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != maxval) {
					min_no = column_norm(A, i, 0);
					min_norm = column_norm(A, i, norm);
					index = i;
					SW = true;
				}
			}
			else {
				if (Pair) {
					no = column_norm(A, i, 0);
					nnorm = column_norm(A, i, norm);
					if (min_norm > nnorm) {
						min_no = no;
						min_norm = nnorm;
						index = i;
					}
					else if (min_norm == nnorm)
						if (min_no > no) {
							min_no = no;
							min_norm = nnorm;
							index = i;
						}
				}
				else {
					if (val != maxval) {
						no = column_norm(A, i, 0);
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_no = no;
							min_norm = nnorm;
							index = i;
						}
						else if (min_norm == nnorm)
							if (min_no > no) {
								min_no = no;
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
row_oriented_dense_matrix_modules< T >::sorting_gcd_plus_minimal_norm (MR< T > &A,
								       lidia_size_t startr,
								       lidia_size_t startc,
								       lidia_size_t &index,
								       lidia_size_t norm) const
{
	lidia_size_t FOUND = 0;
	lidia_size_t index_largest = 0;
	T *tmp = A.value[startr];
	T maxval, val;

	for (lidia_size_t i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (FOUND == 0) {
				maxval = val;
				index_largest = i;
			}
			else
				if (maxval < val) {
					maxval = val;
					index_largest = i;
				}
			FOUND++;
		}
	}

	bool SW = false;

	if (FOUND > 1) {
		for (lidia_size_t i = 0; i <= startc; i++) {
			val = abs(tmp[i]);
			if (val != A.Zero && index_largest != i) {
				if (SW == false) {
					maxval = val;
					index = i;
					SW = true;
				}
				else
					if (maxval <= val) {
						maxval = val;
						index = i;
					}
			}
		}
	}
	else {
		index = index_largest;
		maxval = abs(tmp[index]);
	}

	bigint nnorm, min_norm = column_norm(A, index, norm);

	for (lidia_size_t i = 0; i <= startc; i++)
		if (abs(tmp[i]) == maxval) {
			nnorm = column_norm(A, index, norm);
			if (nnorm < min_norm) {
				min_norm = nnorm;
				index = i;
			}
		}

	return FOUND;
}



template< class T >
lidia_size_t
row_oriented_dense_matrix_modules< T >::min_no_of_elements_plus_minimal_norm (MR< T > &A,
									      lidia_size_t startr,
									      lidia_size_t startc,
									      lidia_size_t &index,
									      lidia_size_t norm) const
{
	lidia_size_t i, FOUND = 0;
	bigint min_no, no;
	bigint min_norm, nnorm, TMP;
	T val, maxval;

	T *tmp = A.value[startr];
	bool Pair = false;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (val == maxval)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(tmp[i]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != maxval) {
					min_no = column_norm(A, i, norm);
					min_norm = column_norm(A, i, 0);
					index = i;
					SW = true;
				}
			}
			else {
				if (Pair) {
					no = column_norm(A, i, norm);
					nnorm = column_norm(A, i, 0);
					if (min_norm > nnorm) {
						min_no = no;
						min_norm = nnorm;
						index = i;
					}
					else if (min_norm == nnorm)
						if (min_no > no) {
							min_no = no;
							min_norm = nnorm;
							index = i;
						}
				}
				else {
					if (val != maxval) {
						no = column_norm(A, i, norm);
						nnorm = column_norm(A, i, 0);
						if (min_norm > nnorm) {
							min_no = no;
							min_norm = nnorm;
							index = i;
						}
						else if (min_norm == nnorm)
							if (min_no > no) {
								min_no = no;
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



//
// Norm
//

template< class T >
inline bigint
row_oriented_dense_matrix_modules< T >::column_norm (MR< T > &A,
						     lidia_size_t pos,
						     lidia_size_t norm) const
{
	bigint Norm = A.Zero;
	bigint TMP;
	lidia_size_t j;


	if (norm == 0) {
		for (j = 0; j < A.rows; j++)
			if (A.value[j][pos] != A.Zero)
				Norm++;
	}
	else if (norm == 1) {
		for (j = 0; j < A.rows; j++)
			LiDIA::add(Norm, Norm, A.value[j][pos]);
	}
	else {
		for (j = 0; j < A.rows; j++) {
			power(TMP, A.value[j][pos], norm);
			LiDIA::add(Norm, Norm, abs(TMP));
		}
	}
	return Norm;
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::kennwerte (MR< T > &RES,
						   T &MAX,
						   lidia_size_t &no_of_elements,
						   T &Durch) const
{
	register lidia_size_t i, j;
	bool SW = false;
	no_of_elements = 0;
	Durch = 0;

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			if (RES.value[i][j] != RES.Zero) {
				no_of_elements++;
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
	}
	Durch /= no_of_elements;
}



template< class T >
inline void
row_oriented_dense_matrix_modules< T >::max (MR< T > &RES, T &MAX) const
{
	register lidia_size_t i, j;
	bool SW = false;

	for (i = 0; i < RES.rows; i++)
		for (j = 0; j < RES.columns; j++)
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
row_oriented_dense_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
								      lidia_size_t startr,
								      lidia_size_t startc,
								      lidia_size_t index,
								      lidia_size_t len) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));
			if (q != 0)
				subtract_multiple_of_column(A, i, q, index, len);
			break;
		}
	return true;
}



template< class T >
inline bool
row_oriented_dense_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
								      matrix< bigint > &TR,
								      lidia_size_t startr,
								      lidia_size_t startc,
								      lidia_size_t index,
								      lidia_size_t len) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));
			if (q != 0) {
				subtract_multiple_of_column(A, i, q, index, len);
				subtract_multiple_of_column(TR, i, q, index, TR.get_no_of_rows() - 1);
			}
			break;
		}
	return true;
}



template< class T >
inline bool
row_oriented_dense_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
								      math_vector< bigfloat > &vec,
								      lidia_size_t startr,
								      lidia_size_t startc,
								      lidia_size_t index,
								      lidia_size_t len) const
{
	T q, TMP;
	bigfloat rtemp;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));
			if (q != 0) {
				subtract_multiple_of_column(A, i, q, index, len);
				LiDIA::multiply(rtemp, bigfloat(q), vec.member(index));
				LiDIA::subtract(vec[i], vec[i], rtemp);
			}
			break;
		}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_linear (MR< T > &A,
						     lidia_size_t startr,
						     lidia_size_t startc,
						     lidia_size_t len) const
{
	T TMP, RES0, RES1, RES2, TMP1, x, y;

	T *tmp = A.value[startr];

	// Init
	lidia_size_t index;
	for (index = startc; tmp[index] == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc)
		this->swap_columns(A, startc, index);

	for (lidia_size_t i = index - 1; i >= 0; i--)
		if (tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[startc]);
			x = tmp[startc] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][startc];

				A.value[l][i] = (TMP*x - TMP1*y);
				A.value[l][startc] = (TMP*RES0 + TMP1*RES1);
			}
		}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_linear (MR< T > &A,
						     matrix< bigint > &TR,
						     lidia_size_t startr,
						     lidia_size_t startc,
						     lidia_size_t len) const
{
	T TMP, RES0, RES1, RES2, TMP1, x, y;

	bigint TMP_1, TMP_2;
	T *tmp = A.value[startr];

	// Init
	lidia_size_t index;
	for (index = startc; tmp[index] == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc) {
		this->swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	for (lidia_size_t i = index - 1; i >= 0; i--)
		if (tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[startc]);
			x = tmp[startc] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][startc];

				A.value[l][i] = TMP*x - TMP1*y;
				A.value[l][startc] = TMP*RES0 + TMP1*RES1;
			}

			for (lidia_size_t l = 0; l < TR.get_no_of_rows(); l++) {
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
row_oriented_dense_matrix_modules< T >::mgcd_bradley (MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t len) const
{
	lidia_size_t n = startc + 1;
	register lidia_size_t i, j;


	matrix< bigint > Tr(A.columns, A.columns);
	Tr.diag(1, 0);

	T *a = A.value[startr];

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		bigint y, z, TMP, TMP1, t1, t2;
		bigint g_old = a[0];
		bigint g = g_old;
		for (i = 1; i < n; i++) {
			if (a[i] != A.Zero) {
				g = xgcd(y, z, g_old, a[i]);
				TMP = -a[i]/g;
				TMP1 = g_old/g;
				for (j = 0; j < Tr.rows; j++) {
					t1 = Tr.value[j][i-1];
					t2 = Tr.value[j][i];
					Tr.value[j][i-1] = TMP*t1+TMP1*t2;
					Tr.value[j][i] = y*t1+z*t2;
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
		bigint q, r;
		lidia_size_t l, k;
		for (i = n - 1, j = n - 2; i >= 0 && j >= 0; i--, j--)
			for (k = j + 1; k < n; k++) {
				if (Tr.value[i][j] != Tr.Zero) {
					div_rem(q, r, Tr.value[i][k], Tr.value[i][j]);

					Tr.value[i][k] = Tr.value[i][k] - q*Tr.value[i][j];
					TMP = Tr.value[i][k] + Tr.value[i][j];
					TMP1 = Tr.value[i][k] - Tr.value[i][j];
					if (abs(TMP) < abs(Tr.value[i][k])) {
						q--;
						Tr.value[i][k] = TMP;
					}
					else if (abs(TMP1) < abs(Tr.value[i][k])) {
						q++;
						Tr.value[i][k] = TMP1;
					}

					if (abs(q) > 0)
						for (l = 0; l < i; l++)
							Tr.value[l][k] = Tr.value[l][k] - q*Tr.value[l][j];
				}
			}
	}
	//multiply(A, A, Tr);
	T TMP;
	T *tmp2 = new T[A.columns];
	for (i = 0; i < A.rows; i++) {

		for (lidia_size_t p = 0; p < A.columns; p++)
			tmp2[p] = this->member(A, i, p);

		for (j = 0; j < A.columns; j++) {
			TMP = A.Zero;
			for (lidia_size_t l = 0; l < A.columns; l++)
				if (tmp2[l] != A.Zero && Tr.member(l, j) != A.Zero)
					TMP += T(tmp2[l] * Tr.member(l, j));

			this->sto(A, i, j, TMP);
		}
	}
	delete[] tmp2;

	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_bradley (MR< T > &A,
						      matrix< bigint > &TR,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t len) const
{
	lidia_size_t n = startc + 1;
	register lidia_size_t i, j;


	matrix< bigint > Tr(A.columns, A.columns);
	Tr.diag(1, 0);

	T *a = A.value[startr];

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		bigint y, z, TMP, TMP1, t1, t2;
		bigint g_old = a[0];
		bigint g = g_old;
		for (i = 1; i < n; i++) {
			if (a[i] != A.Zero) {
				g = xgcd(y, z, g_old, a[i]);
				TMP = -a[i]/g;
				TMP1 = g_old/g;
				for (j = 0; j < n; j++) {
					t1 = Tr.value[j][i-1];
					t2 = Tr.value[j][i];
					Tr.value[j][i-1] = TMP*t1+TMP1*t2;
					Tr.value[j][i] = y*t1+z*t2;
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
		bigint q, r;
		lidia_size_t l, k;
		for (i = n - 1, j = n - 2; i >= 0 && j >= 0; i--, j--)
			for (k = j + 1; k < n; k++) {
				if (Tr.value[i][j] != Tr.Zero) {
					div_rem(q, r, Tr.value[i][k], Tr.value[i][j]);

					Tr.value[i][k] = Tr.value[i][k] - q*Tr.value[i][j];
					TMP = Tr.value[i][k] + Tr.value[i][j];
					TMP1 = Tr.value[i][k] - Tr.value[i][j];
					if (abs(TMP) < abs(Tr.value[i][k])) {
						q--;
						Tr.value[i][k] = TMP;
					}
					else if (abs(TMP1) < abs(Tr.value[i][k])) {
						q++;
						Tr.value[i][k] = TMP1;
					}

					if (abs(q) > 0)
						for (l = 0; l < i; l++)
							Tr.value[l][k] = Tr.value[l][k] - q*Tr.value[l][j];
				}
			}
	}
	//multiply(A, A, Tr);
	T TMP;
	T *tmp2 = new T[A.columns];
	for (i = 0; i < A.rows; i++) {

		for (lidia_size_t p = 0; p < A.columns; p++)
			tmp2[p] = this->member(A, i, p);

		for (j = 0; j < A.columns; j++) {
			TMP = A.Zero;
			for (lidia_size_t l = 0; l < A.columns; l++)
				if (tmp2[l] != A.Zero && Tr.member(l, j) != A.Zero)
					TMP += T(tmp2[l] * Tr.member(l, j));

			this->sto(A, i, j, TMP);
		}
	}
	delete[] tmp2;

	LiDIA::multiply(TR, TR, Tr);
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_ilio (MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t len) const
{
	lidia_size_t n = startc + 1;

	bigint *a = A.value[startr];

	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	T g = A.Zero;
	if (n > 1) {
		bigint y, z, TMP, TMP1, t1, t2;

		lidia_size_t step = 1, step2 = 2;
		while (step < n) {
			for (i = 0; i < n; i += step2) {
				if (i + step < n) {
					g = xgcd(y, z, a[i], a[i+step]);
					if (g != A.Zero) {
						TMP = -a[i+step]/g;
						TMP1 = a[i]/g;

						for (j = 0; j < A.rows; j++) {
							t1 = A.value[j][i];
							t2 = A.value[j][i+step];
							A.value[j][i] = y*t1+z*t2;
							A.value[j][i+step] = TMP*t1+TMP1*t2;
						}
					}
				}
			}
			step = step2;
			step2 = 2*step2;
		}
		for (j = 0; j < A.rows; j++)
		  LiDIA::swap(A.value[j][0], A.value[j][n-1]);
		if (A.value[startr][startc] < 0)
		  for (j = 0; j < A.rows; j++)
		    A.value[j][startc].negate();
	}
	if (g == A.Zero)
		return false;
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_ilio (MR< T > &A,
						   matrix< bigint > &Tr,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t len) const
{
	lidia_size_t n = startc + 1;

	bigint *a = A.value[startr];

	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	T g = A.Zero;
	if (n > 1) {
		bigint y, z, TMP, TMP1, t1, t2;

		lidia_size_t step = 1, step2 = 2;
		while (step < n) {
			for (i = 0; i < n; i += step2) {
				if (i + step < n) {
					g = xgcd(y, z, a[i], a[i+step]);
					if (g != A.Zero) {
						TMP = -a[i+step]/g;
						TMP1 = a[i]/g;
						for (j = 0; j < A.columns; j++) {
							t1 = Tr.member(j, i);
							t2 = Tr.member(j, i + step);
							Tr.sto(j, i, y*t1+z*t2);
							Tr.sto(j, i + step, TMP*t1+TMP1*t2);
						}

						for (j = 0; j < A.rows; j++) {
							t1 = A.value[j][i];
							t2 = A.value[j][i+step];
							A.value[j][i] = y*t1+z*t2;
							A.value[j][i+step] = TMP*t1+TMP1*t2;
						}
					}
				}
			}
			step = step2;
			step2 = 2*step2;
		}
		for (j = 0; j < A.rows; j++)
			LiDIA::swap(A.value[j][0], A.value[j][n-1]);
		Tr.swap_columns(0, n-1);
		if (A.value[startr][startc] < 0)
		  {
		    for (j = 0; j < A.rows; j++)
		      A.value[j][startc].negate();
		    for (j = 0; j < A.columns; j++)
		      Tr.sto(j, startc, -Tr.member(j, startc));
		  }
	}
	if (g == A.Zero)
		return false;
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_opt (MR< T > &A,
						  lidia_size_t startr,
						  lidia_size_t startc,
						  lidia_size_t len) const
{
	T TMP, RES0, RES1, RES2, TMP1, x, y;

	T *tmp = A.value[startr];

	// Init
	lidia_size_t index;
	for (index = startc; tmp[index] == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc)
		this->swap_columns(A, startc, index);

	// sort
	if (startc >= 2) {
		lidia_size_t gcdindex = startc - 1;
		T gcd_min = gcd(tmp[startc], tmp[startc - 1]), gcd_tmp;
		for (lidia_size_t i = startc - 2; i >= 0; i--) {
			gcd_tmp = gcd(tmp[startc], tmp[i]);
			if (gcd_tmp < gcd_min)
				gcdindex = i;
		}
		this->swap_columns(A, gcdindex, startc - 1);
	}

	for (lidia_size_t i = startc - 1; i >= 0; i--)
		if (tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[startc]);
			x = tmp[startc] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][startc];

				A.value[l][i] = (TMP*x - TMP1*y);
				A.value[l][startc] = (TMP*RES0 + TMP1*RES1);
			}
		}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_opt (MR< T > &A,
						  matrix< bigint > &TR,
						  lidia_size_t startr,
						  lidia_size_t startc,
						  lidia_size_t len) const
{
	T TMP, RES0, RES1, RES2, TMP1, x, y;

	bigint TMP_1, TMP_2;
	T *tmp = A.value[startr];

	// Init
	lidia_size_t index;
	for (index = startc; tmp[index] == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc) {
		this->swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	// sort
	if (startc >= 2) {
		lidia_size_t gcdindex = startc - 1;
		T gcd_min = gcd(tmp[startc], tmp[startc - 1]), gcd_tmp;
		for (lidia_size_t i = startc - 2; i >= 0; i--) {
			gcd_tmp = gcd(tmp[index], tmp[i]);
			if (gcd_tmp < gcd_min)
				gcdindex = i;
		}
		this->swap_columns(A, gcdindex, startc - 1);
		TR.swap_columns(gcdindex, startc - 1);
	}


	for (lidia_size_t i = startc - 1; i >= 0; i--)
		if (tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[startc]);
			x = tmp[startc] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][startc];

				A.value[l][i] = TMP*x - TMP1*y;
				A.value[l][startc] = TMP*RES0 + TMP1*RES1;
			}

			for (lidia_size_t l = 0; l < TR.get_no_of_rows(); l++) {
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
row_oriented_dense_matrix_modules< T >::mgcd_storjohann (MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t len) const
{
	T q, res, TMP;
	T *tmp;
	T g, N, a, qa, ra, qb, rb, t;
	T A3, A2, A1, A0, x, y;
	T TMP1, TMP2, TMP3;

	register lidia_size_t i, j, l;
	lidia_size_t index;

	//
	// init
	//

	tmp = A.value[startr];
	for (index = startc; index >= 0 && tmp[index].is_zero(); index--);

	if (index < 0)
		return false;
	else {
		//
		// conditioning
		//

		N = tmp[index];
		for (j = index - 1; j >= 0 && tmp[j].is_zero(); j--);
		if (j >= 0) {
			a = tmp[j];
			t = 0;

			for (i = j - 1; i >= 0; i--) {
				if (tmp[i] != 0) {
					g = gcd(tmp[i], a);

					div_rem(qa, ra, a/g, N);
					div_rem(qb, rb, tmp[i]/g, N);

					for (t = 0; gcd((ra + t*rb), N) != 1; t++);

					a = a + t * tmp[i]; // NEW

					//M.sto(1, i, t);
					for (l = 0; l <= startr; l++)
						LiDIA::add(A.value[l][j], A.value[l][j], t*A.value[l][i]);
				}
			}

			//
			// elimination
			//

			A2 = xgcd(A0, A1, tmp[j], tmp[index]);
			div_rem(x, A3, tmp[index], A2);
			div_rem(y, A3, tmp[j], A2);

			for (l = 0; l <= startr; l++) {
				TMP = A.value[l][j];
				TMP1 = A.value[l][index];


				// Atmp1[j] = ((TMP * x) + (TMP1 * y))
				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				LiDIA::subtract(A.value[l][j], TMP2, TMP3);

				// Atmp1[n-m+i] = ((TMP * A0) + (TMP1 * A1))
				LiDIA::multiply(TMP2, TMP, A0);
				LiDIA::multiply(TMP3, TMP1, A1);
				LiDIA::add(A.value[l][index], TMP2, TMP3);
			}

			for (i = j - 1; i >= 0; i--)
				if (!tmp[i].is_zero()) {
					div_rem(TMP, TMP1, tmp[i], tmp[index]);
					for (l = 0; l <= startr; l++)
						LiDIA::subtract(A.value[l][i], A.value[l][i], TMP*A.value[l][index]);
				}

		}

		if (tmp[index].is_lt_zero())
			for (i = 0; i <= startr; i++)
				A.value[i][index].negate();
		if (index != startc)
			this->swap_columns(A, startc, index);
	}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::mgcd_storjohann (MR< T > &A,
							 matrix< bigint > &TR,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t len) const
{
	T q, res, TMP;
	T *tmp;
	T g, N, a, qa, ra, qb, rb, t;
	T A3, A2, A1, A0, x, y;
	T TMP1, TMP2, TMP3;

	register lidia_size_t i, j, l;
	lidia_size_t index;

	//
	// init
	//

	tmp = A.value[startr];
	for (index = startc; index >= 0 && tmp[index].is_zero(); index--);

	if (index < 0)
		return false;
	else {
		//
		// conditioning
		//

		N = tmp[index];
		for (j = index - 1; j >= 0 && tmp[j].is_zero(); j--);
		if (j >= 0) {
			a = tmp[j];
			t = 0;

			for (i = j - 1; i >= 0; i--) {
				if (tmp[i] != 0) {
					g = gcd(tmp[i], a);

					div_rem(qa, ra, a/g, N);
					div_rem(qb, rb, tmp[i]/g, N);

					for (t = 0; gcd((ra + t*rb), N) != 1; t++);

					a = a + t * tmp[i]; // NEW

					//M.sto(1, i, t);
					for (l = 0; l <= startr; l++)
						LiDIA::add(A.value[l][j], A.value[l][j], t*A.value[l][i]);
					for (l = 0; l < A.columns; l++)
						TR.sto(l, j, TR.member(l, j) + t * TR.member(l, i));
				}
			}

			//
			// elimination
			//

			A2 = xgcd(A0, A1, tmp[j], tmp[index]);
			div_rem(x, A3, tmp[index], A2);
			div_rem(y, A3, tmp[j], A2);

			for (l = 0; l <= startr; l++) {
				TMP = A.value[l][j];
				TMP1 = A.value[l][index];


				// Atmp1[j] = ((TMP * x) + (TMP1 * y))
				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				LiDIA::subtract(A.value[l][j], TMP2, TMP3);

				// Atmp1[n-m+i] = ((TMP * A0) + (TMP1 * A1))
				LiDIA::multiply(TMP2, TMP, A0);
				LiDIA::multiply(TMP3, TMP1, A1);
				LiDIA::add(A.value[l][index], TMP2, TMP3);
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
				if (!tmp[i].is_zero()) {
					div_rem(TMP, TMP1, tmp[i], tmp[index]);
					for (l = 0; l <= startr; l++)
						LiDIA::subtract(A.value[l][i], A.value[l][i], TMP*A.value[l][index]);
					for (l = 0; l < A.columns; l++)
						TR.sto(l, i, TR.member(l, i) - TMP * TR.member(l, index));
				}

		}

		if (tmp[index].is_lt_zero()) {
			for (i = 0; i <= startr; i++)
				A.value[i][index].negate();
			for (i = 0; i < A.columns; i++)
				TR.sto(i, index, -TR.member(i, index));
		}


		if (index != startc) {
			this->swap_columns(A, startc, index);
			TR.swap_columns(startc, index);
		}
	}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::xgcd_elimination (MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t index,
							  lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	T TMP2, TMP3;

	T *tmp = A.value[startr];

	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[index]);
			x = tmp[index] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][index];

				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				LiDIA::subtract(A.value[l][i], TMP2, TMP3);

				LiDIA::multiply(TMP2, TMP, RES0);
				LiDIA::multiply(TMP3, TMP1, RES1);
				LiDIA::add(A.value[l][index], TMP2, TMP3);
			}
		}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::xgcd_elimination (MR< T > &A,
							  matrix< bigint > &TR,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t index,
							  lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	T TMP2, TMP3;

	T *tmp = A.value[startr];

	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[index]);
			x = tmp[index] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][index];

				LiDIA::multiply(TMP2, TMP, x);
				LiDIA::multiply(TMP3, TMP1, y);
				LiDIA::subtract(A.value[l][i], TMP2, TMP3);

				LiDIA::multiply(TMP2, TMP, RES0);
				LiDIA::multiply(TMP3, TMP1, RES1);
				LiDIA::add(A.value[l][index], TMP2, TMP3);
			}

			bigint TTMP, TTMP1, TTMP2, TTMP3, TTMP4;
			for (lidia_size_t l = 0; l < A.columns; l++) {
				TTMP = TR.member(l, i);
				TTMP1 = TR.member(l, index);

				LiDIA::multiply(TTMP2, TTMP, x);
				LiDIA::multiply(TTMP3, TTMP1, y);
				TR.sto(l, i, TTMP2 - TTMP3);

				LiDIA::multiply(TTMP2, TTMP, RES0);
				LiDIA::multiply(TTMP3, TTMP1, RES1);
				TR.sto(l, index, TTMP2 + TTMP3);
			}
		}
	return true;
}



template< class T >
bool
row_oriented_dense_matrix_modules< T >::xgcd_elimination_mod (MR< T > &A,
							      lidia_size_t startr,
							      lidia_size_t startc,
							      lidia_size_t index,
							      lidia_size_t len,
							      const T &mod) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	T TMP2, TMP3;

	T *tmp = A.value[startr];

	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && tmp[i] != A.Zero) {
			RES2 = xgcd(RES0, RES1, tmp[i], tmp[index]);
			x = tmp[index] / RES2;
			y = tmp[i] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[l][i];
				TMP1 = A.value[l][index];

				LiDIA::mult_mod(TMP2, TMP, x, mod);
				LiDIA::mult_mod(TMP3, TMP1, y, mod);
				LiDIA::sub_mod(A.value[l][i], TMP2, TMP3, mod);

				LiDIA::mult_mod(TMP2, TMP, RES0, mod);
				LiDIA::mult_mod(TMP3, TMP1, RES1, mod);
				LiDIA::add_mod(A.value[l][index], TMP2, TMP3, mod);
			}
		}
	return true;
}



template< class T >
inline bool
row_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
						       lidia_size_t startr,
						       lidia_size_t startc,
						       lidia_size_t index,
						       lidia_size_t len) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));
			if (q != 0)
				subtract_multiple_of_column(A, i, q, index, len);
		}
	return true;
}



template< class T >
inline bool
row_oriented_dense_matrix_modules< T >::normalize_row_mod (MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t index,
							   lidia_size_t len,
							   const T &DET) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));
			if (q != 0)
				subtract_multiple_of_column_mod(A, i, q, index, len, DET);
		}
	return true;
}



template< class T >
inline bool
row_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
						       matrix< bigint > &TR,
						       lidia_size_t startr,
						       lidia_size_t startc,
						       lidia_size_t index,
						       lidia_size_t len) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));

			// Check
			if (q != 0) {
				subtract_multiple_of_column(A, i, q, index, len);

				for (lidia_size_t j = 0; j < A.columns; j++)
					TR.sto(j, i, TR.member(j, i) - q * TR.member(j, index));
			}
		}
	return true;
}



template< class T >
inline bool
row_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
						       lidia_size_t startr,
						       lidia_size_t startc,
						       lidia_size_t index,
						       lidia_size_t len,
						       T *max_array,
						       const T &BOUND) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));

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
row_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
						       matrix< bigint > &TR,
						       lidia_size_t startr,
						       lidia_size_t startc,
						       lidia_size_t index,
						       lidia_size_t len,
						       T *max_array,
						       const T &BOUND) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && this->member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, this->member(A, startr, i), this->member(A, startr, index));

			// Check
			if (q != 0) {
				if (abs(bigint(bigint(max_array[index]) * q) + max_array[i]) > bigint(BOUND))
					return false;
				subtract_multiple_of_column(A, i, q, index, len);

				for (lidia_size_t j = 0; j < A.columns; j++)
					TR.sto(j, i, TR.member(j, i) - q* TR.member(j, index));

				update_max_array(A, i, max_array);
			}
		}
	return true;
}



//
// special functions
//

template< class T >
inline const T &
column_oriented_dense_matrix_modules< T >::member (const MR< T > &A,
						   lidia_size_t x,
						   lidia_size_t y) const
{
	return A.value[y][x];
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::swap_columns (const MR< T > &A,
							 lidia_size_t l1,
							 lidia_size_t l2) const
{
	T *tmp = A.value[l1];
	A.value[l1] = A.value[l2];
	A.value[l2] = tmp;
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::subtract_multiple_of_column (MR< T > &A,
									lidia_size_t l1,
									const T &q,
									lidia_size_t l2,
									lidia_size_t l) const
{
	for (register lidia_size_t i = 0; i <= l; i++)
		LiDIA::subtract(A.value[l1][i], A.value[l1][i], q*A.value[l2][i]);
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::subtract_multiple_of_column_mod (MR< T > &A,
									    lidia_size_t l1,
									    const T &q,
									    lidia_size_t l2,
									    lidia_size_t l,
									    const T &mod) const
{
	T TMP;
	for (register lidia_size_t i = 0; i <= l; i++) {
		LiDIA::mult_mod(TMP, q, A.value[l2][i], mod);
		LiDIA::sub_mod(A.value[l1][i], A.value[l1][i], TMP, mod);
	}
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::normalize_column_mod (MR< T > &A,
								 lidia_size_t l1,
								 const T &q,
								 lidia_size_t l2,
								 lidia_size_t l,
								 const T &mod) const
{
	T TMP;
	for (register lidia_size_t i = 0; i <= l; i++) {
		LiDIA::mult_mod(TMP, q, A.value[l2][i], mod);
		LiDIA::sub_mod(A.value[l1][i], A.value[l1][i], TMP, mod);
		if (A.value[l1][i] < A.Zero)
			A.value[l1][i] += mod;
	}
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::negate_column (MR< T > &A,
							  lidia_size_t index,
							  lidia_size_t l) const
{
	for (register lidia_size_t i = 0; i <= l; i++)
		A.value[index][i] = -A.value[index][i];
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::negate_column_mod (MR< T > &A,
							      lidia_size_t index,
							      lidia_size_t l,
							      const T &mod) const
{
	for (register lidia_size_t i = 0; i <= l; i++) {
		A.value[index][i] = -A.value[index][i];
		LiDIA::best_remainder(A.value[index][i], A.value[index][i], mod);
	}
}



template< class T >
inline T *
column_oriented_dense_matrix_modules< T >::init_max_array (const MR< T > &A) const
{
	T *max_array = new T[A.columns];
	memory_handler(max_array, DMESSAGE, "init :: "
		       "Error in memory allocation (max_array)");
	T *tmp;
	for (register lidia_size_t i = 0; i < A.columns; i++) {
		tmp = A.value[i];
		for (register lidia_size_t j = 0; j < A.rows; j++)
			if (max_array[j] < abs(tmp[j]))
				max_array[j] = abs(tmp[j]);
	}
	return max_array;
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::update_max_array (const MR< T > &A,
							     lidia_size_t i,
							     T *max_array) const
{
	for (register lidia_size_t j = 0; j < A.rows; j++)
		if (max_array[i] < abs(A.value[i][j]))
			max_array[i] = abs(A.value[i][j]);
}



//
// Pivot search
//

template< class T >
inline lidia_size_t
column_oriented_dense_matrix_modules< T >::normal_pivot (const MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t &index) const
{
	lidia_size_t FOUND = 0;

	for (lidia_size_t i = 0; i <= startc; i++)
		if (A.value[i][startr] != A.Zero) {
			FOUND++;
			index = i;
			for (lidia_size_t j = 0; j <= startc; j++)
				if (A.value[j][startr] != A.Zero && (j != index)) {
					FOUND++;
					if (j != index) {
						if (abs(A.value[j][startr]) < abs(A.value[i][startr]))
							index = j;
						break;
					}
				}
			break;
		}
	return FOUND;
}



template< class T >
inline lidia_size_t
column_oriented_dense_matrix_modules< T >::min_abs_of_row (const MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc,
							   lidia_size_t &pos) const
{
	lidia_size_t FOUND = 0;
	for (; startc >= 0 && !FOUND; startc--)
		if (A.value[startc][startr] != A.Zero) {
			pos = startc;
			FOUND++;
		}

	for (; startc >= 0; startc--)
		if (A.value[startc][startr] != A.Zero) {
			FOUND++;
			if (abs(A.value[pos][startr]) >= abs(A.value[startc][startr]))
				pos = startc;
		}

	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_dense_matrix_modules< T >::pivot_sorting_gcd (const MR< T > &A,
							      lidia_size_t startr,
							      lidia_size_t startc,
							      lidia_size_t &index) const
{
	lidia_size_t FOUND = 0;
	lidia_size_t index_largest = 0;
	T val, maxval;

	for (lidia_size_t i = 0; i <= startc; i++) {
		val = abs(A.value[i][startr]);
		if (val != A.Zero) {
			if (FOUND == 0) {
				maxval = val;
				index_largest = i;
			}
			else
				if (maxval < val) {
					maxval = val;
					index_largest = i;
				}

			FOUND++;
		}
	}

	bool SW = false;
	if (FOUND > 1) {
		for (lidia_size_t i = 0; i <= startc; i++) {
			val = abs(A.value[i][startr]);
			if (val != A.Zero && index_largest != i) {
				if (SW == false) {
					maxval = val;
					index = i;
				}
				else
					if (maxval < val) {
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
column_oriented_dense_matrix_modules< T >::minimal_norm (const MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t &index,
							 lidia_size_t norm) const
{
	lidia_size_t i, j, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T val, maxval, minval;

	bool Pair = false;
	bool SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(A.value[i][startr]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (maxval == val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(A.value[i][startr]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || maxval != val) {
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
					if (val != maxval) {
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
column_oriented_dense_matrix_modules< T >::min_abs_of_row_plus_minimal_norm (const MR< T > &A,
									     lidia_size_t startr,
									     lidia_size_t startc,
									     lidia_size_t &index,
									     lidia_size_t norm) const
{
	lidia_size_t i, j, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T min_val = 0, val;
	bool SW = false;
	bool Pair = false;

	for (i = startc; i >= 0; i--) {
		val = abs(A.value[i][startr]);

		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				index = i;
				min_val = val;
				SW = true;
			}
			else {
				if (min_val > val) {
					index = i;
					min_val = val;
				}
			}
		}
	}

	min_norm = column_norm(A, index, norm);

	// norm computation
	for (i = index - 1; i >= 0; i--) {
		if (min_val == abs(A.value[i][startr])) {
			nnorm = column_norm(A, i, norm);
			if (nnorm < min_norm) {
				index = i;
				min_norm = nnorm;
			}
		}
	}
	return FOUND;
}



template< class T >
lidia_size_t
column_oriented_dense_matrix_modules< T >::minimal_norm_plus_min_abs_of_row (const MR< T > &A,
									     lidia_size_t startr,
									     lidia_size_t startc,
									     lidia_size_t &index,
									     lidia_size_t norm) const
{
	lidia_size_t i, j, FOUND = 0;
	bigint min_norm, nnorm, TMP;
	T val, maxval, minval;

	bool SW = false;
	bool Pair = false;

	for (i = 0; i <= startc; i++) {
		val = abs(A.value[i][startr]);
		if (val != A.Zero) {
			FOUND++;
			if (SW == false) {
				maxval = val;
				SW = true;
				index = i;
			}
			else {
				if (val > maxval) {
					maxval = val;
					Pair = false;
				}
				else if (maxval == val)
					Pair = true;
			}
		}
	}

	SW = false;

	for (i = 0; i <= startc; i++) {
		val = abs(A.value[i][startr]);
		if (val != A.Zero) {
			if (SW == false) {
				if (Pair || val != maxval) {
					minval = val;
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
					else if (min_norm == nnorm)
						if (minval > val) {
							minval = val;
							min_norm = nnorm;
							index = i;
						}
				}
				else {
					if (val != maxval) {
						nnorm = column_norm(A, i, norm);
						if (min_norm > nnorm) {
							min_norm = nnorm;
							index = i;
						}
						else if (min_norm == nnorm)
							if (minval > val) {
								minval = val;
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
column_oriented_dense_matrix_modules< T >::sorting_gcd_plus_minimal_norm (const MR< T > &A,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t &index,
									  lidia_size_t norm) const
{
	lidia_size_t FOUND = 0;
	lidia_size_t index_largest = 0;
	T val, maxval;

	for (lidia_size_t i = 0; i <= startc; i++) {
		val = abs(A.value[i][startr]);
		if (val != A.Zero) {
			if (FOUND == 0) {
				maxval = val;
				index_largest = i;
			}
			else
				if (maxval < val) {
					maxval = val;
					index_largest = i;
				}

			FOUND++;
		}
	}

	bool SW = false;
	if (FOUND > 1) {
		for (lidia_size_t i = 0; i <= startc; i++) {
			val = abs(A.value[i][startr]);
			if (val != A.Zero && index_largest != i) {
				if (SW == false) {
					maxval = val;
					index = i;
				}
				else
					if (maxval < val) {
						maxval = val;
						index = i;
					}

				SW = true;
			}
		}
	}
	else {
		index = index_largest;
		maxval = abs(A.value[index][startr]);
	}

	bigint nnorm, min_norm = column_norm(A, index, norm);

	for (lidia_size_t i = 0; i <= startc; i++)
		if (abs(A.value[i][startr]) == maxval) {
			nnorm = column_norm(A, i, norm);
			if (nnorm < min_norm) {
				min_norm = nnorm;
				index = i;
			}
		}
	return FOUND;
}



//
// Norm
//

template< class T >
inline bigint
column_oriented_dense_matrix_modules< T >::column_norm (const MR< T > &A,
							lidia_size_t pos,
							lidia_size_t norm) const
{
	bigint Norm = A.Zero;
	T * tmp = A.value[pos];
	bigint TMP;

	if (norm == 0) {
		for (lidia_size_t j = 0; j < A.rows; j++)
			if (tmp[j] != A.Zero)
				Norm++;
	}
	else if (norm == 1) {
		for (lidia_size_t j = 0; j < A.rows; j++)
			LiDIA::add(Norm, Norm, tmp[j]);
	}
	else {
		for (lidia_size_t j = 0; j < A.rows; j++) {
			power(TMP, tmp[j], norm);
			LiDIA::add(Norm, Norm, abs(TMP));
		}
	}
	return Norm;
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::kennwerte (MR< T > &RES,
						      T &MAX,
						      lidia_size_t &no_of_elements,
						      T &Durch) const
{
	register lidia_size_t i, j;
	bool SW = false;
	no_of_elements = 0;
	Durch = 0;

	for (i = 0; i < RES.columns; i++) {
		for (j = 0; j < RES.rows; j++) {
			if (RES.value[i][j] != RES.Zero) {
				no_of_elements++;
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
	}
	Durch /= no_of_elements;
}



template< class T >
inline void
column_oriented_dense_matrix_modules< T >::max (MR< T > &RES, T &MAX) const
{
	register lidia_size_t i, j;
	bool SW = false;

	for (i = 0; i < RES.columns; i++)
		for (j = 0; j < RES.rows; j++)
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
column_oriented_dense_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
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
column_oriented_dense_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
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
				subtract_multiple_of_column(TR, i, q, index, TR.get_no_of_rows() - 1);
			}
			break;
		}
	return true;
}



template< class T >
inline bool
column_oriented_dense_matrix_modules< T >::normalize_one_element_of_row (MR< T > &A,
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
column_oriented_dense_matrix_modules< T >::mgcd_linear (MR< T > &A,
							lidia_size_t startr,
							lidia_size_t startc,
							lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;

	// Init
	lidia_size_t index;
	for (index = startc; A.value[index][startr] == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc)
		swap_columns(A, startc, index);

	for (lidia_size_t i = index - 1; i >= 0; i--)
		if (A.value[i][startr] != A.Zero) {
			RES2 = xgcd(RES0, RES1, A.value[i][startr], A.value[startc][startr]);
			x = A.value[startc][startr] / RES2;
			y = A.value[i][startr] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[i][l];
				TMP1 = A.value[startc][l];

				A.value[i][l] = TMP*x - TMP1*y;
				A.value[startc][l] = TMP*RES0 + TMP1*RES1;
			}
		}
	return true;
}



template< class T >
bool
column_oriented_dense_matrix_modules< T >::mgcd_linear (MR< T > &A,
							matrix< bigint > &TR,
							lidia_size_t startr,
							lidia_size_t startc,
							lidia_size_t len) const
{
	T TMP;
	T RES0, RES1, RES2, TMP1, TMP2, x, y;

	// Init
	lidia_size_t index;
	for (index = startc; A.value[index][startr] == A.Zero && index >= 0; index--);
	if (index < 0)
		return false;
	if (index != startc) {
		swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	for (lidia_size_t i = index - 1; i >= 0; i--)
		if (A.value[i][startr] != A.Zero) {
			RES2 = xgcd(RES0, RES1, A.value[i][startr], A.value[startc][startr]);
			x = A.value[startc][startr] / RES2;
			y = A.value[i][startr] / RES2;

			for (lidia_size_t l = 0; l <= len; l++) {
				TMP = A.value[i][l];
				TMP1 = A.value[startc][l];

				A.value[i][l] = TMP*x - TMP1*y;
				A.value[startc][l] = TMP*RES0 + TMP1*RES1;
			}

			for (lidia_size_t l = 0; l < TR.get_no_of_rows(); l++) {
				TMP1 = TR(l, i);
				TMP2 = TR(l, startc);

				TR.sto(l, i, (TMP1*x - TMP2*y));
				TR.sto(l, startc, (TMP1*RES0 + TMP2*RES1));
			}
		}
	return true;
}



template< class T >
bool
column_oriented_dense_matrix_modules< T >::xgcd_elimination_mod (MR< T > &A,
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
				LiDIA::sub_mod(A.value[i][l], TMP2, TMP3, mod);

				LiDIA::mult_mod(TMP2, TMP, RES0, mod);
				LiDIA::mult_mod(TMP3, TMP1, RES1, mod);
				LiDIA::add_mod(A.value[index][l], TMP2, TMP3, mod);
			}
		}
	return true;
}



template< class T >
inline bool
column_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
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
		}
	return true;
}



template< class T >
inline bool
column_oriented_dense_matrix_modules< T >::normalize_row_mod (MR< T > &A,
							      lidia_size_t startr,
							      lidia_size_t startc,
							      lidia_size_t index,
							      lidia_size_t len,
							      const T &DET) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));
			if (q != 0)
				subtract_multiple_of_column_mod(A, i, q, index, len, DET);
		}
	return true;
}



template< class T >
inline bool
column_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
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

			// Check
			if (q != 0) {
				subtract_multiple_of_column(A, i, q, index, len);
				subtract_multiple_of_column(TR, i, q, index, TR.get_no_of_rows() - 1);
			}
		}
	return true;
}



template< class T >
inline bool
column_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t index,
							  lidia_size_t len,
							  T *max_array,
							  const T &BOUND) const
{
	T q, TMP;
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
inline bool
column_oriented_dense_matrix_modules< T >::normalize_row (MR< T > &A,
							  matrix< bigint > &TR,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t index,
							  lidia_size_t len,
							  T *max_array,
							  const T &BOUND) const
{
	T q, TMP;
	for (lidia_size_t i = 0; i <= startc; i++)
		if ((i != index) && member(A, startr, i) != A.Zero) {
			pos_div_rem(q, TMP, member(A, startr, i), member(A, startr, index));

			// Check
			if (q != 0) {
				if (abs(bigint(bigint(max_array[index]) * q) + max_array[i]) > bigint(BOUND))
					return false;
				subtract_multiple_of_column(A, i, q, index, len);
				subtract_multiple_of_column(TR, i, q, index, TR.get_no_of_rows() - 1);
				update_max_array(A, i, max_array);
			}
		}
	return true;
}



#undef row_oriented_dense_matrix_modules
#undef column_oriented_dense_matrix_modules



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_DENSE_BIGINT_MATRIX_MODULES_CC_GUARD_
