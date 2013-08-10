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


#ifndef LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_CC_GUARD_


#ifndef LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_bigint_matrix_kernel.h"
#endif
#ifndef LIDIA_MODULAR_OPERATIONS_INL_GUARD_
# include	"LiDIA/modular_operations.inl"
#endif
#ifndef LIDIA_MODULAR_FUNCTIONS_INL_GUARD_
# include	"LiDIA/matrix/modular_functions.inl"
#endif
#ifndef LIDIA_CRT_AND_PRIME_HANDLING_H_GUARD_
# include	"LiDIA/matrix/crt_and_prime_handling.h"
#endif

#include "LiDIA/precondition_error.h"


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


#define DMESSAGE "sparse_bigint_matrix_kernel"  // Debug message
#define EMESSAGE matrix_error_msg               // Error message

//
// divide
//

template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::divide (MATRIX_TYPE &RES,
						    const MATRIX_TYPE &A,
						    const bigint &k) const
{
	register lidia_size_t j, i;
	bigint *REStmp, *Atmp;

	for (j = 0; j < A.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		for (i = 0; i < A.value_counter[j]; i++)
			LiDIA::divide(REStmp[i], Atmp[i], k);
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::compwise_divide (MATRIX_TYPE &RES,
							     const MATRIX_TYPE &A,
							     const MATRIX_TYPE &B) const
{
	register lidia_size_t j, i;
	bigint *REStmp, *Atmp, *Btmp;

	for (j = 0; j < RES.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		Btmp = B.value[j];
		for (i = 0; i < RES.columns; i++)
			LiDIA::divide(REStmp[i], Atmp[i], Btmp[i]);
	}
}



//
// remainder
//

template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::remainder (MATRIX_TYPE &RES,
						       const MATRIX_TYPE & M,
						       const bigint & mod) const
{
	lidia_size_t i, j;
	bigint *REStmp, *Mtmp;

	for (i = 0; i < M.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < M.value_counter[i]; j++)
			LiDIA::best_remainder(REStmp[j], Mtmp[j], mod);
	}
}



//
// norms and bounds
//

template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::max (const MATRIX_TYPE &RES,
						 bigint &MAX) const
{
	MAX.assign(RES.value[0][0]);

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MAX < tmp[j])
				MAX.assign(tmp[j]);
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::max_abs (const MATRIX_TYPE &RES,
						     bigint &MAX) const
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



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::max_pos (const MATRIX_TYPE &RES,
						     bigint & MAX,
						     lidia_size_t & x,
						     lidia_size_t & y) const
{
	MAX.assign(RES.value[0][0]);
	x = 0;
	y = 0;

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MAX < tmp[j]) {
				MAX.assign(tmp[j]);
				x = i;
				y = j;
			}
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::max_abs_pos (const MATRIX_TYPE &RES,
							 bigint & MAX,
							 lidia_size_t & x,
							 lidia_size_t & y) const
{
	MAX.assign(abs(RES.value[0][0]));
	x = 0;
	y = 0;

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MAX < abs(tmp[j])) {
				MAX.assign(abs(tmp[j]));
				x = i;
				y = j;
			}
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::min (const MATRIX_TYPE &RES,
						 bigint &MIN) const
{
	MIN.assign(RES.value[0][0]);

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MIN > tmp[j])
				MIN.assign(tmp[j]);
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::min_abs (const MATRIX_TYPE &RES,
						     bigint &MIN) const
{
	MIN.assign(abs(RES.value[0][0]));

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MIN > abs(tmp[j]))
				MIN.assign(abs(tmp[j]));
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::min_pos (const MATRIX_TYPE &RES,
						     bigint & MIN,
						     lidia_size_t & x,
						     lidia_size_t & y) const
{
	MIN.assign(RES.value[0][0]);
	x = 0;
	y = 0;

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MIN > tmp[j]) {
				MIN.assign(tmp[j]);
				x = i;
				y = j;
			}
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::min_abs_pos (const MATRIX_TYPE &RES,
							 bigint & MIN,
							 lidia_size_t & x,
							 lidia_size_t & y) const
{
	MIN.assign(abs(RES.value[0][0]));
	x = 0;
	y = 0;

	register lidia_size_t i, j;
	register bigint *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (MIN > abs(tmp[j])) {
				MIN.assign(abs(tmp[j]));
				x = i;
				y = j;
			}
	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::hadamard (const MATRIX_TYPE &RES,
						      bigint & H) const
{
	register lidia_size_t min = (RES.columns < RES.rows) ? RES.columns : RES.rows;
	register lidia_size_t i, j;

	// Computation of row and column norms
	register bigint *Wrows = new bigint[RES.rows];
	memory_handler(Wrows, DMESSAGE, "hadamard :: "
		       "Error in memory allocation (Wrows)");

	register bigint *Wcolumns = new bigint[RES.columns];
	memory_handler(Wcolumns, DMESSAGE, "hadamard :: "
		       "Error in memory allocation (Wcolumns)");

	register bigint *tmp;
	bigint TMP;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.value_counter[i]; j++) {
			LiDIA::multiply(TMP, tmp[j], tmp[j]);
			LiDIA::add(Wrows[i], Wrows[i], TMP);
			LiDIA::add(Wcolumns[RES.index[i][j]], Wcolumns[RES.index[i][j]], TMP);
		}
	}

	// Column evaluation
	j = 0;
	if (RES.rows < RES.columns)
		while (j < RES.columns - 1)
			if (Wcolumns[j] < Wcolumns[j + 1]) {
				LiDIA::swap(Wcolumns[j], Wcolumns[j + 1]);
				j = 0;
			}				// h.sort(DOWN);
			else
				j++;

	bigint COLUMNS;
	COLUMNS.assign_one();
	for (j = 0; j < min; j++)
		if (!Wcolumns[j].is_zero())
			LiDIA::multiply(COLUMNS, COLUMNS, Wcolumns[j]);
	delete[] Wcolumns;

	// Row evaluation
	j = 0;
	if (RES.rows > RES.columns)
		while (j < RES.rows - 1)
			if (Wrows[j] < Wrows[j + 1]) {
				LiDIA::swap(Wrows[j], Wrows[j + 1]);
				j = 0;
			}				// h.sort(DOWN);
			else
				j++;

	bigint ROWS;
	ROWS.assign_one();
	for (i = 0; i < min; i++)
		if (!Wrows[i].is_zero())
			LiDIA::multiply(ROWS, ROWS, Wrows[i]);
	delete[] Wrows;

	// Hadamard
	register lidia_size_t B = ((COLUMNS < ROWS) ? COLUMNS.bit_length() : ROWS.bit_length()) - 1;
	bigint E = (bigint(B) / bigint(2)) + bigint(2);
	power(H, bigint(2), E);
	dec(H);
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::binary_hadamard (const MATRIX_TYPE &RES,
							     lidia_size_t &H) const
{
	register lidia_size_t min = (RES.columns < RES.rows) ? RES.columns : RES.rows;
	register lidia_size_t i, j;

	// Computation of row and column norms
	register bigint *Wrows = new bigint[RES.rows];
	memory_handler(Wrows, DMESSAGE, "hadamard :: "
		       "Error in memory allocation (Wrows)");

	register bigint *Wcolumns = new bigint[RES.columns];
	memory_handler(Wcolumns, DMESSAGE, "hadamard :: "
		       "Error in memory allocation (Wcolumns)");

	register bigint *tmp;
	bigint TMP;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++) {
			LiDIA::multiply(TMP, tmp[j], tmp[j]);
			LiDIA::add(Wrows[i], Wrows[i], TMP);
			LiDIA::add(Wcolumns[j], Wcolumns[j], TMP);
		}
	}

	// Column evaluation
	j = 0;
	if (RES.rows < RES.columns)
		while (j < RES.columns - 1)
			if (Wcolumns[j] < Wcolumns[j + 1]) {
				LiDIA::swap(Wcolumns[j], Wcolumns[j + 1]);
				j = 0;
			}				// h.sort(DOWN);
			else
				j++;

	bigint COLUMNS;
	COLUMNS.assign_one();
	for (j = 0; j < min; j++)
		LiDIA::multiply(COLUMNS, COLUMNS, Wcolumns[j]);
	delete[] Wcolumns;

	// Row evaluation
	j = 0;
	if (RES.rows > RES.columns)
		while (j < RES.rows - 1)
			if (Wrows[j] < Wrows[j + 1]) {
				LiDIA::swap(Wrows[j], Wrows[j + 1]);
				j = 0;
			}				// h.sort(DOWN);
			else
				j++;

	bigint ROWS;
	ROWS.assign_one();
	for (i = 0; i < min; i++)
		LiDIA::multiply(ROWS, ROWS, Wrows[i]);
	delete[] Wrows;

	// Hadamard
	register lidia_size_t B = ((COLUMNS < ROWS) ? COLUMNS.bit_length() : ROWS.bit_length()) - 1;
	H = B / 2 + 2;
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::row_norm (const MATRIX_TYPE &M,
						      bigint & RES,
						      lidia_size_t pos,
						      long art) const
{
	RES.assign_zero();

	bigint TMP;
	register bigint *tmp = M.value[pos];

	register lidia_size_t i;

	if (art != 0)
		for (i = 0; i < M.value_counter[pos]; i++) {
			power(TMP, tmp[i], art);
			LiDIA::add(RES, RES, abs(TMP));
		}
	else {
		for (i = 0; i < M.value_counter[pos]; i++) {
			if (!tmp[i].is_zero())
				RES++;
		}
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::column_norm (const MATRIX_TYPE &M,
							 bigint & RES,
							 lidia_size_t pos,
							 long art) const
{
	RES.assign_zero();

	bigint TMP;

	register lidia_size_t i;
	if (art != 0)
		for (i = 0; i < M.rows; i++) {
			power(TMP, M(i, pos), art);
			LiDIA::add(RES, RES, abs(TMP));
		}
	else
		for (i = 0; i < M.rows; i++) {
			if (!M(i, pos).is_zero())
				RES++;
		}
}



//
// randomize
//

template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::randomize (MATRIX_TYPE &RES, const bigint & S) const
{
	register lidia_size_t i, j;
	register bigint *tmp;

	bigint::seed(S);

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			tmp[j].assign(LiDIA::randomize(S));
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::randomize (MATRIX_TYPE &RES,
						       const bigint & S,
						       const long d) const
{
	register lidia_size_t i, j;
	long TMP;
	random_generator gen;

	bigint::seed(S);

	for (i = 0; i < RES.rows; i++) {
		for (j = 0; j < RES.columns; j++) {
			gen >> TMP;
			if (TMP % 100 <= d)
				RES.sto(i, j, LiDIA::randomize(S));
			else
				RES.sto(i, j, RES.Zero);
		}
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::randomize_with_det (MATRIX_TYPE &RES,
								const bigint & S,
								const bigint & DET) const
{
	lidia_size_t i, j;

	MATRIX_TYPE T(RES.rows, RES.columns), T1(RES.rows, RES.columns);
	bigint::seed(S);

	register bigint *tmp = NULL, *tmp1;
	for (i = 0; i < RES.rows; i++) {
		tmp = T.value[i];
		tmp1 = T1.value[i];
		for (j = i+1; j < RES.columns; j++)
			tmp[j].assign(LiDIA::randomize(S));
		for (j = 0; j < i; j++)
			tmp1[j].assign(LiDIA::randomize(S));
		tmp[i].assign_one();
		tmp1[i].assign_one();
	}
	random_generator GEN;
	GEN >> i;
	i %= RES.rows;
	T.value[i][i].assign(DET);
	LiDIA::multiply(RES, T, T1);
}



//
// regular expansion
//

template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::regexpansion (MATRIX_TYPE &RES,
							  const lidia_size_t * v) const
{
	register lidia_size_t k = v[0];
	if (RES.columns > k) {
		register lidia_size_t i = 0, j = 0;
		lidia_size_t diff = RES.columns - RES.rows;
		MATRIX_TYPE A(RES);
		RES.set_no_of_rows(RES.columns);
		RES.set_no_of_columns(RES.columns);
		MATRIX_TYPE ERG(diff, RES.columns);
		while (i< RES.columns && k > 0 && j < diff) {
			if (i != v[k]) {
				ERG.value[j][i].assign_one();
				j++;
			}
			else
				k--;
			i++;
		}
		//dmodul2.compose_v(RES, ERG, A);
		SS_base_modul.insert_at(RES, 0, 0, ERG, 0, 0, ERG.rows, ERG.columns);
		SS_base_modul.insert_at(RES, RES.rows - A.rows, 0, A, 0, 0, A.rows, A.columns);
	}
}



//
// BEGIN: Linear algebra
// PART 2
//

//
// Hermite normal form
//

template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::hnfmod_dkt (MATRIX_TYPE &RES, const bigint &mod) const
{
	//
	// Task: A.hnfmod_dkt(mod);
	// =  > A in Hermite normal form
	// =  > h = lattice determinant of lattice formed
	//              by the columns of matrix A
	// ERROR: rank != rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "in member - function "
			"hnfmod_dkt(const bigint &)", DVALUE + 8);

	// bigint part
	register long i, j, z, diff = RES.columns - RES.rows;

	bigint RES0, RES1, RES2, RES3; // 0 = lggT, 1 = rggt, 2 = ggt
	bigint x, y;
	bigint TMP, TMP1, TMP2, TMP3;
	bigint *Atmp, *Atmp1 = NULL;

	// Step 1, 2
	for (i = 0; i < RES.rows; i++) {
		Atmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			best_remainder(Atmp[j], Atmp[j], mod);
	}

	// Step 3 - 5
	for (i = RES.rows - 1; i >= 0; i--) {
		Atmp = RES.value[i];

		// Step 6 -8
		for (j = diff + i - 1; j >= 0; j--) {
			if (!Atmp[j].is_zero()) {
				// Step 9 - 12
				if (Atmp[diff + i].is_zero())
					Atmp[diff + i].assign(mod);

				// Step 13 - 18
				RES2 = xgcd(RES0, RES1, Atmp[j], Atmp[diff + i]);
				div_rem(x, RES3, Atmp[diff + i], RES2);
				div_rem(y, RES3, Atmp[j], RES2);

				// Step 19 - 28
				for (z = 0; z <= i; z++) {
					Atmp1 = RES.value[z];
					TMP.assign(Atmp1[j]);
					TMP1.assign(Atmp1[diff + i]);

					// Atmp1[j] = ((TMP * x) + (TMP1 * y)) % h;
					mult_mod(TMP2, TMP, x, mod);
					mult_mod(TMP3, TMP1, y, mod);
					sub_mod(Atmp1[j], TMP2, TMP3, mod);

					// Atmp1[n-m+i] = ((TMP * RES0) + (TMP1 * RES1)) % h;
					mult_mod(TMP2, TMP, RES0, mod);
					mult_mod(TMP3, TMP1, RES1, mod);
					add_mod(Atmp1[diff+i], TMP2, TMP3, mod);
				}
			}
		}
	}

	// Step 29 - 32
	bigint D = mod;
	for (i = RES.rows - 1; i >= 0; i--) {
		Atmp = RES.value[i];

		// Step 33 - 36
		if (Atmp[diff + i].is_zero())
			Atmp[diff + i].assign(D);

		// Step 37 - 39
		RES2 = xgcd(RES0, RES1, Atmp[diff + i], D);

		// Step 40 - 46
		for (z = 0; z < RES.rows; z++) {
			Atmp1 = RES.value[z];
			mult_mod(Atmp1[diff + i], Atmp1[diff + i], RES0, D);
			if (Atmp1[diff+i].is_lt_zero())
				LiDIA::add(Atmp1[diff+i], Atmp1[diff+i], D);
		}
		Atmp[RES.columns - RES.rows + i].assign(RES2);

		// Step 47 - 49
		div_rem(D, TMP, D, RES2);
	}

	// Step 50 - 52
	for (i = RES.rows - 1; i >= 0; i--) {
		Atmp = RES.value[i];

		// Step 53 - 56
		if (Atmp[diff + i].is_zero())
			Atmp[diff + i].assign(mod);

		// Step 57 - 59
		for (j = diff + i + 1; j < RES.columns; j++) {
			pos_div_rem(TMP, TMP1, Atmp[j], Atmp[diff + i]);

			// Step 60 - 66
			for (z = 0; z <= i; z++) {
				Atmp1 = RES.value[z];
				LiDIA::multiply(TMP1, Atmp1[diff + i], TMP);
				LiDIA::subtract(Atmp1[j], Atmp1[j], TMP1);
			}
		}
	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::hnfmod_cohen (MATRIX_TYPE &RES,
							  const bigint & D) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"hnfmod_cohen(const bigint &)", DVALUE + 8);

	matrix< bigint > W(RES.rows, RES.columns);
	long i = RES.rows - 1;
	long j = RES.columns - 1;
	long k = RES.columns - 1;
	long z;

	bigint R = D;
	bigint u, v, d, q;
	bigint TMP, TMP1, TMP2;
	bigint *tmp, *tmp1;

	for (i = RES.rows - 1; i >= 0; i--) {
		tmp = RES.value[i];
		if (tmp[k].is_zero())
			tmp[k].assign(R);
		while (j != 0) {
			j--;
			if (!tmp[j].is_zero()) {
				d = xgcd(u, v, tmp[k], tmp[j]);
				div_rem(TMP1, TMP, tmp[j], d);
				div_rem(TMP2, TMP, tmp[k], d);
				for (z = 0; z < RES.rows; z++) {
					tmp1 = RES.value[z];
					TMP.assign(u * tmp1[k] + v * tmp1[j]);
					best_remainder(tmp1[j], TMP2 * tmp1[j] - TMP1 * tmp1[k], R);
					best_remainder(tmp1[k], TMP, R);
				}
			}
		}
		d = xgcd(u, v, tmp[k], R);
		for (z = 0; z < RES.rows; z++) {
			best_remainder(W.value[z][k], u * RES.value[z][k], R);
			if (W.value[z][k] < 0)
				LiDIA::add(W.value[z][k], W.value[z][k], R);
		}
		if (W.value[i][k].is_zero())
			W.value[i][k].assign(R);
		for (j = k + 1; j < RES.columns; j++) {
			pos_div_rem(q, TMP, W.value[i][j], W.value[i][k]);
			for (z = 0; z < RES.rows; z++)
				LiDIA::subtract(W.value[z][j], W.value[z][j], q * W.value[z][k]);
		}
		div_rem(R, TMP, R, d);
		k--;
		j = k;
	}
	RES.assign(W);
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::hnfmod_mueller (MATRIX_TYPE &RES,
							    MATRIX_TYPE & TRANS,
							    bigint &H) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"hnfmod_mueller(MATRIX_TYPE &, const bigint &)", DVALUE + 8);

	register long i, j;
	bigint TMP, *TRANStmp;

	// Step 1
	lidia_size_t *linuz = RES.lininc(H);
	if (linuz[0] != RES.rows)
		precondition_error_handler(linuz[0], "rank", "rank == rows",
				    RES.rows, "rows", "rank == rows",
				    "void MATRIX_TYPE::"
				    "hnfmod_mueller(MATRIX_TYPE & TRANS)",
				    DMESSAGE, EMESSAGE[10]);

	// Step 2, 3
	if (RES.rows == RES.columns) {
		MATRIX_TYPE A1(RES);
		bigint DET;
		RES.det(DET);
		if (DET.is_lt_zero())
			hnfmod_dkt(RES, -DET);
		else
			hnfmod_dkt(RES, DET);
		MATRIX_TYPE A2(A1.rows, A1.columns);
		bigint h;
		hadamard(A1, h);
		A2.adj(A1, h);
		LiDIA::multiply(TRANS, A2, RES);

		for (i = 0; i < RES.rows; i++) {
			TRANStmp = TRANS.value[i];
			for (j = 0; j < RES.rows; j++)
				div_rem(TRANStmp[j], TMP, TRANStmp[j], DET);
		}
	}
	else {
		lidia_size_t diff = RES.columns - RES.rows;

		MATRIX_TYPE EAR(RES);
		regexpansion(EAR, linuz);
		MATRIX_TYPE R(RES.rows, RES.rows), F(RES.rows, diff);
		MATRIX_TYPE E(RES.columns, RES.columns), P(RES.columns, RES.columns);
		E.diag(bigint(1), bigint(0));

		long k = linuz[0], l1 = 0, l2 = 0;
		for (i = 0; i < RES.columns; i++) {
			if (linuz[k] == i) {
				k--;
				// R.sto_column(column(i), rows, l1);
				for (j = 0; j < RES.rows; j++)
					R.value[j][l1].assign(RES.value[j][i]);
				// P.sto_column(E.column(i), columns, l1);
				for (j = 0; j < RES.columns; j++)
					P.value[j][l1].assign(E.value[j][i]);
				l1++;
			}
			else {
				//F.sto_column(column(i), rows, l2);
				for (j = 0; j < RES.rows; j++)
					F.value[j][l2].assign(RES.value[j][i]);
				//P.sto_column(E.column(i), columns, rows + l2);
				for (j = 0; j < RES.columns; j++)
					P.value[j][RES.rows+l2].assign(E.value[j][i]);
				l2++;
			}
		}
		delete[] linuz;

		// Step 4, 5
		hnfmod_dkt(EAR, H);
		MATRIX_TYPE ADJ;
		hadamard(R, H);
		ADJ.adj(R, H);
		bigint DET;
		R.det(DET);
		MATRIX_TYPE U(diff, diff), V(diff, RES.rows);
		MATRIX_TYPE EHNF(RES.rows, RES.rows), Null(RES.rows, diff);

		//dmodul2.split_t(EAR, U, V, Null, EHNF);
		SS_base_modul.insert_at(U, 0, 0, EAR, 0, 0, U.rows, U.columns);
		SS_base_modul.insert_at(V, 0, 0, EAR, 0, EAR.columns - V.columns, V.rows, V.columns);
		SS_base_modul.insert_at(Null, 0, 0, EAR, EAR.rows - Null.rows, 0, Null.rows, Null.columns);
		SS_base_modul.insert_at(EHNF, 0, 0, EAR, EAR.rows - EHNF.rows, EAR.columns - EHNF.columns,
					EHNF.rows, EHNF.columns);

		MATRIX_TYPE M(RES.columns, RES.columns);

		// Step 6
		MATRIX_TYPE T1 = (-ADJ) * (F * U);
		MATRIX_TYPE T2 = (-ADJ) * (F * V) + (ADJ * EHNF);
		MATRIX_TYPE T3 = (U * DET);
		MATRIX_TYPE T4 = (V * DET);
		//dmodul2.compose_t(M, T1, T2, T3, T4);
		SS_base_modul.insert_at(M, 0, 0, T1, 0, 0, T1.rows, T1.columns);
		SS_base_modul.insert_at(M, 0, M.columns - T2.columns, T2, 0, 0, T2.rows, T2.columns);
		SS_base_modul.insert_at(M, M.rows - T3.rows, 0, T3, 0, 0, T3.rows, T3.columns);
		SS_base_modul.insert_at(M, M.rows - T4.rows, M.columns - T4.columns, T4, 0, 0, T4.rows, T4.columns);

		// Step 7
		LiDIA::multiply(TRANS, P, M);
		for (i = 0; i < RES.columns; i++) {
			TRANStmp = TRANS.value[i];
			for (j = 0; j < RES.columns; j++)
				div_rem(TRANStmp[j], TMP, TRANStmp[j], DET);
		}
		LiDIA::multiply(RES, RES, TRANS);
	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::hnf_storjohann (MATRIX_TYPE &RES) const
{
	bigint q, res, TMP;
	bigint *tmp;
	bigint g, N, a, qa, ra, qb, rb, t;
	bigint RES3, RES2, RES1, RES0, x, y;
	bigint TMP1, TMP2, TMP3;

	register lidia_size_t startr, startc, i, j, l;
	lidia_size_t index;

	//
	// Elimination
	//

	for (startc = RES.columns - 1, startr = RES.rows - 1; startr >= 0 && startc >= 0; startc--, startr--) {
		//
		// init
		//

		tmp = RES.value[startr];
		for (index = startc; index >= 0 && member(RES, startr, index) == RES.Zero; index--);

		if (index < 0)
			startc++;
		else {
			//
			// conditioning
			//

			N = member(RES, startr, index);
			for (j = index - 1; j >= 0 && member(RES, startr, j) == 0; j--);
			if (j >= 0) {
				a = member(RES, startr, j);
				t = 0;

				for (i = j - 1; i >= 0; i--) {
					if (member(RES, startr, i) != 0) {
						g = gcd(member(RES, startr, i), a);

						div_rem(qa, ra, a/g, N);
						div_rem(qb, rb, member(RES, startr, i)/g, N);

						for (t = 0; gcd((ra + t*rb), N) != 1; t++);

						a = a + t * member(RES, startr, i); // NEW

						//M.sto(1, i, t);
						for (l = 0; l <= startr; l++) {
							LiDIA::add(TMP, member(RES, l, j), t*member(RES, l, i));
							sto(RES, l, j, TMP);
						}
					}
				}

				//
				// elimination
				//

				RES2 = xgcd(RES0, RES1, member(RES, startr, j), member(RES, startr, index));
				div_rem(x, RES3, member(RES, startr, index), RES2);
				div_rem(y, RES3, member(RES, startr, j), RES2);

				for (l = 0; l <= startr; l++) {
					TMP = member(RES, l, j);
					TMP1 = member(RES, l, index);


					// Atmp1[j] = ((TMP * x) + (TMP1 * y))
					LiDIA::multiply(TMP2, TMP, x);
					LiDIA::multiply(TMP3, TMP1, y);
					LiDIA::subtract(TMP2, TMP2, TMP3);
					sto(RES, l, j, TMP2);

					// Atmp1[n-m+i] = ((TMP * RES0) + (TMP1 * RES1))
					LiDIA::multiply(TMP2, TMP, RES0);
					LiDIA::multiply(TMP3, TMP1, RES1);
					LiDIA::add(TMP2, TMP2, TMP3);
					sto(RES, l , index, TMP2);
				}

				for (i = j - 1; i >= 0; i--)
					if (member(RES, startr, i) != RES.Zero) {
						div_rem(TMP, TMP1, member(RES, startr, i), member(RES, startr, index));
						for (l = 0; l <= startr; l++) {
							LiDIA::subtract(TMP2, member(RES, l, i), TMP*member(RES, l, index));
							sto(RES, l, i, TMP2);
						}
					}

			}

			if (member(RES, startr, index) < RES.Zero)
				for (i = 0; i <= startr; i++) {
					TMP = - member(RES, i , index);
					sto(RES, i, index, TMP);

				}
			if (index != startc)
				RES.swap_columns(startc, index);
		}
	}

	for (startc = RES.columns - 2, startr = RES.rows - 2; startr >= 0 && startc >= 0; startr--, startc--) {
		tmp = RES.value[startr];

		if (member(RES, startr, startc) != RES.Zero) {
			for (i = startc + 1; i < RES.columns; i++) {
				pos_div_rem(q, TMP, member(RES, startr, i), member(RES, startr, startc));
				for (j = 0; j <= startr; j++) {
					LiDIA::multiply(TMP, q, member(RES, j , startc));
					LiDIA::subtract(TMP, member(RES, j, i), TMP);
					sto(RES, j, i, TMP);
				}
			}
		}
		else
			startc++;
	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::hnf_storjohann (MATRIX_TYPE &RES,
							    MATRIX_TYPE &TR,
							    MATRIX_TYPE &C,
							    MATRIX_TYPE &Q) const
{
	MATRIX_TYPE A = RES;

	bigint q, res, TMP;
	bigint *tmp, *tmp1;
	bigint g, N, a, qa, ra, qb, rb, t;
	bigint RES3, RES2, RES1, RES0, x, y;
	bigint TMP1, TMP2, TMP3;

	register lidia_size_t startr, startc, i, j, l;
	lidia_size_t index;

	MATRIX_TYPE Q1(RES.columns, RES.columns);
	MATRIX_TYPE C1(RES.columns, RES.columns);

	Q.diag(1, 0);
	C.diag(1, 0);

	//
	// Elimination
	//

	for (startc = RES.columns - 1, startr = RES.rows - 1; startr >= 0 && startc >= 0; startc--, startr--) {
		//
		// init
		//

		Q1.diag(1, 0);
		C1.diag(1, 0);

		tmp = RES.value[startr];
		for (index = startc; index >= 0 && tmp[index].is_zero(); index--);

		if (index < 0)
			startc++;
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

						for (l = 0; l <= startr; l++)
							LiDIA::add(RES.value[l][j], RES.value[l][j], t*RES.value[l][i]);

						C1.value[i][j] = t;
					}
				}

				//
				// elimination
				//

				RES2 = xgcd(RES0, RES1, tmp[j], tmp[index]);
				div_rem(x, RES3, tmp[index], RES2);
				div_rem(y, RES3, tmp[j], RES2);

				for (l = 0; l <= startr; l++) {
					TMP = RES.value[l][j];
					TMP1 = RES.value[l][index];

					// Atmp1[j] = ((TMP * x) - (TMP1 * y))
					LiDIA::multiply(TMP2, TMP, x);
					LiDIA::multiply(TMP3, TMP1, y);
					LiDIA::subtract(RES.value[l][j], TMP2, TMP3);

					// Atmp1[n-m+i] = ((TMP * RES0) + (TMP1 * RES1))
					LiDIA::multiply(TMP2, TMP, RES0);
					LiDIA::multiply(TMP3, TMP1, RES1);
					LiDIA::add(RES.value[l][index], TMP2, TMP3);

				}

				Q1.value[j][j] = x;
				Q1.value[index][j] = -y;

				Q1.value[index][index] = RES1;
				Q1.value[j][index] = RES0;


				for (i = j - 1; i >= 0; i--) {
					if (!tmp[i].is_zero()) {
						div_rem(TMP, TMP1, tmp[i], tmp[index]);
						for (l = 0; l <= startr; l++)
							LiDIA::subtract(RES.value[l][i], RES.value[l][i], TMP*RES.value[l][index]);

						LiDIA::subtract(Q1.value[index][i], Q1.value[index][i], TMP*Q1.value[index][index]);
						LiDIA::subtract(Q1.value[j][i], Q1.value[j][i], TMP*Q1.value[j][index]);
					}
				}

			}

			if (tmp[index].is_lt_zero()) {
				for (i = 0; i <= startr; i++)
					RES.value[i][index].negate();
				if (j != -1)
					Q1.value[j][index].negate();
				Q1.value[index][index].negate();
			}
			if (index != startc) {
				RES.swap_columns(startc, index);
				Q1.swap_columns(startc, index);
			}
		}

		//
		// Update of the transformation matrix
		C = C*C1;
		Q = Q*LiDIA::adj(C1)*Q1;
		TR = TR * C1 * Q1;
	}

	for (startc = RES.columns - 2, startr = RES.rows - 2; startr >= 0 && startc >= 0; startr--, startc--) {
		tmp = RES.value[startr];
		for (i = startc + 1; i < RES.columns; i++) {
			pos_div_rem(q, TMP, tmp[i], tmp[startc]);
			for (j = 0; j <= startr; j++) {
				tmp1 = RES.value[j];
				LiDIA::multiply(TMP, q, tmp1[startc]);
				LiDIA::subtract(tmp1[i], tmp1[i], TMP);
			}

			for (j = 0; j <= TR.rows; j++) {
				LiDIA::multiply(TMP, q, TR.value[j][startc]);
				LiDIA::subtract(TR.value[j][i], TR.value[j][i], TMP);
			}
		}
	}
}



//
// Kernel
//

template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::kernel1 (MATRIX_TYPE &RES,
						     const MATRIX_TYPE & A) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"kernel1(const MATRIX_TYPE &)", DVALUE + 8);

	bigint *ZBAtmp, *Atmp;
	register long i, j;
	lidia_size_t c = A.columns;

	// Step 1
	bigint H;
	hadamard(A, H);
	lidia_size_t *linuz = A.lininr(H);
	lidia_size_t r = linuz[0];
	if (r == c) {
		MATRIX_TYPE RES1(c, 1);
		RES.assign(RES1);
		delete[] linuz;
		return;
	}

	// Step 2
	MATRIX_TYPE ZBA(r, c);
	for (i = 1; i <= r; i++) {
		ZBAtmp = ZBA.value[i - 1];
		Atmp = A.value[linuz[r - i + 1]];
		for (j = 0; j < c; j++)
			ZBAtmp[j].assign(Atmp[j]);
	}
	delete[] linuz;

	// Step 3
	MATRIX_TYPE TRANS(c, c);
	ZBA.hnf(TRANS);

	// Step 4
	MATRIX_TYPE PART2(c, r);
	if (RES.rows != c)
		RES.set_no_of_rows(c);
	if (RES.columns != c - r)
		RES.set_no_of_columns(c - r);
	//dmodul2.split_h(TRANS, RES, PART2);
	SS_base_modul.insert_at(RES, 0, 0, TRANS, 0, 0, RES.rows, RES.columns);
	SS_base_modul.insert_at(PART2, 0, 0, TRANS, 0, TRANS.columns - PART2.columns, PART2.rows, PART2.columns);
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::kernel2 (MATRIX_TYPE &RES,
						     const MATRIX_TYPE & A) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"kernel2(const MATRIX_TYPE &)", DVALUE + 8);

	register lidia_size_t i;
	MATRIX_TYPE B = A;
	//hnf_havas_cont(B, RES);

	for (i = 0; i < A.columns && B.is_column_zero(i); i++);

	if (i == 0) {
		MATRIX_TYPE C(A.columns, 1);
		RES.assign(C);
	}
	else
		RES.set_no_of_columns(i);
}



//
// regular InvImage
//

template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::reginvimage1 (MATRIX_TYPE &RES,
							  const MATRIX_TYPE & A,
							  const MATRIX_TYPE & B) const
{
	//
	// Task: C.reginvimage1(A, B);
	// =  > A * C.column(j) = g(j)*B.column(j), j = 0, ..., B.columns
	// =  > g(j) minimal
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "in member - function "
			"reginvimage1(const MATRIX_TYPE &, const MATRIX_TYPE &", DVALUE + 8);

	register long i, j;
	bigint TMP, TMP1;

	// Step 1
	bigint H;
	hadamard(A, H);
	bigint DET;
	A.det(DET, H);
	if (DET == 0) {
		precondition_error_handler(DET, "det(A)", "det(A) != 0",
				    "void MATRIX_TYPE::"
				    "reginvimage1(const MATRIX_TYPE & A, const MATRIX_TYPE & B)",
				    DMESSAGE, EMESSAGE[11]);
		return;
	}
	MATRIX_TYPE ADJ;
	ADJ.adj(A, H);

	// Step 2
	MATRIX_TYPE PROD = ADJ * B;
	bigint *u = new bigint[B.rows];
	memory_handler(u, DMESSAGE, "reginvimage1 :: "
		       "Error in memory allocation (u)");
	bigint *g, phi;

	// Step 3
	if (RES.rows != B.rows + 1)
		RES.set_no_of_rows(B.rows + 1);
	if (RES.columns != B.columns)
		RES.set_no_of_columns(B.columns);
	for (i = 0; i < B.columns; i++) {
		for (j = 0; j < PROD.rows; j++)
			u[j].assign(PROD.value[j][i]);
		g = LiDIA::mgcd2(u, PROD.rows);
		div_rem(phi, TMP, DET, gcd(g[0], DET));
		if (phi.is_lt_zero())
			phi.negate();

		// Step 4
		for (j = 0; j < PROD.rows; j++) {
			LiDIA::multiply(TMP, phi, u[j]);
			div_rem(RES.value[j][i], TMP1, TMP, DET);
		}
		RES.value[PROD.rows][i].assign(phi);
		delete[] g;
	}
	delete[] u;
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::reginvimage2 (MATRIX_TYPE &RES,
							  const MATRIX_TYPE & A,
							  const MATRIX_TYPE & B) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"reginvimage2(const MATRIX_TYPE &, const MATRIX_TYPE &", DVALUE + 8);

	register lidia_size_t i, j, len, oldlen;
	bigint TMP, TMP1;

	// Step 1
	bigint H;
	hadamard(A, H);
	bigint DET;
	A.det(DET, H);
	if (DET == 0) {
		precondition_error_handler(DET, "det(A)", "det(A) != 0",
				    "void MATRIX_TYPE::"
				    "reginvimage2(const MATRIX_TYPE & A, const MATRIX_TYPE & B)",
				    DMESSAGE, EMESSAGE[11]);
		return;
	}

	// Step 2
	oldlen = B.rows;
	len = B.rows + 1;
	LiDIA::multiply(RES, LiDIA::adj(A), B);

	bigint *u = new bigint[len];
	memory_handler(u, DMESSAGE, "reginvimage :: "
		       "Error in memory allocation (u)");
	bigint phi;

	// Step 3
	RES.set_no_of_rows(len);
	for (i = 0; i < B.columns; i++) {
		for (j = 0; j < oldlen; j++)
			u[j].assign(RES.value[j][i]);
		u[oldlen].assign(DET);
		LiDIA::mgcd2(TMP1, u, len);
		div_rem(phi, TMP, DET, TMP1);
		if (phi.is_lt_zero())
			phi.negate();

		// Step 4
		for (j = 0; j < oldlen; j++) {
			LiDIA::multiply(TMP, phi, u[j]);
			div_rem(RES.value[j][i], TMP1, TMP, DET);
		}
		RES.value[oldlen][i].assign(phi);
	}
	delete[] u;
}



//
// Image
//

template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::image1 (MATRIX_TYPE &RES,
						    const MATRIX_TYPE & A) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"image1(const MATRIX_TYPE &)", DVALUE + 8);

	bigint *ZBAtmp, *Atmp;
	register long i, j;

	// Step 1
	bigint H;
	hadamard(A, H);
	lidia_size_t *v = A.lininr(H);
	lidia_size_t RANG = v[0];

	// Step 2, 3
	MATRIX_TYPE ZBA(RANG, A.columns);
	for (i = 1; i <= RANG; i++) {
		ZBAtmp = ZBA.value[i - 1];
		Atmp = A.value[v[RANG - i + 1]];
		for (j = 0; j < A.columns; j++)
			ZBAtmp[j].assign(Atmp[j]);
	}
	delete[] v;

	// Step 4
	MATRIX_TYPE TRANS(A.columns, A.columns);
	ZBA.hnf(TRANS);

	// Step 5
	if (RES.rows != A.rows)
		RES.set_no_of_rows(A.rows);
	if (RES.columns != RANG)
		RES.set_no_of_columns(RANG);

	if (A.columns == RANG)
		LiDIA::multiply(RES, A, TRANS);
	else {
		MATRIX_TYPE M(A.rows, A.columns);
		LiDIA::multiply(M, A, TRANS);
		MATRIX_TYPE PART1(RANG, A.columns - RANG);
		//dmodul2.split_h(M, PART1, RES);
		SS_base_modul.insert_at(PART1, 0, 0, M, 0, 0, PART1.rows, PART1.columns);
		SS_base_modul.insert_at(RES, 0, 0, M, 0, M.columns - RES.columns, RES.rows, RES.columns);
	}
}



template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::image2 (MATRIX_TYPE &RES,
						    const MATRIX_TYPE & A) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"image2(const MATRIX_TYPE &)", DVALUE + 8);

	register lidia_size_t i;
	RES.assign(A);
	// hnf_havas_cont(RES);

	for (i = 0; i < RES.columns && RES.is_column_zero(i); i++);

	if (i != 0)
		if (i == RES.columns)
			RES.set_no_of_columns(1);
		else {
			MATRIX_TYPE M(RES);
			MATRIX_TYPE PART1(RES.rows, i);
			RES.set_no_of_columns(RES.columns-i);
			//dmodul2.split_h(M, PART1, RES);
			SS_base_modul.insert_at(PART1, 0, 0, M, 0, 0, PART1.rows, PART1.columns);
			SS_base_modul.insert_at(RES, 0, 0, M, 0, M.columns - RES.columns, RES.rows, RES.columns);
		}
}



//
// InvImage
//

template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::invimage (MATRIX_TYPE &RES,
						      const MATRIX_TYPE & B,
						      const bigint * b) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"invimage(const MATRIX_TYPE &, const bigint *)", DVALUE + 8);

	if (b == NULL)
		precondition_error_handler(PRT, "b", "b != NULL",
				    "void MATRIX_TYPE::"
				    "invimage(const MATRIX_TYPE & B, const bigint * b)",
				    DMESSAGE, EMESSAGE[1]);

	register long i;
	bigint *tmp;

	// Step 1
	MATRIX_TYPE A = B;
	A.set_no_of_columns(B.columns + 1);
	for (i = 0; i < B.rows; i++)
		A.value[i][B.columns].assign(-b[i]);
	kernel1(RES, A);

	// Step 2
	if (RES.is_column_zero(0) || RES.is_row_zero(B.columns)) {
		MATRIX_TYPE C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	tmp = new bigint[RES.columns];
	RES.get_row(tmp, B.columns);
	bigint *g = LiDIA::mgcd2(tmp, RES.columns);
	delete[] tmp;
	if (g[0] > 1) {
		MATRIX_TYPE C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	// Step 3, 4 
	bigint *x = (RES) * &(g[1]);
	delete[] g;

	// Step 5
	kernel1(RES, B);
	RES.set_no_of_columns(RES.columns + 1);
	for (i = 0; i < RES.rows; i++)
		RES.value[i][RES.columns-1].assign(x[i]);
	delete[] x;
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::invimage (MATRIX_TYPE &RES,
						      const MATRIX_TYPE & B,
						      const math_vector< bigint > &b) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"invimage(const MATRIX_TYPE &, const math_vector< bigint > &)", DVALUE + 8);

	register long i;
	// Step 1
	MATRIX_TYPE A = B;
	A.set_no_of_columns(B.columns + 1);

	if (b.size() != B.rows)
		precondition_error_handler(b.size(), "b.size", "b.size == B.rows",
				    B.rows, "B.rows", "b.size == B.rows",
				    "void MATRIX_TYPE::"
				    "invimage(const MATRIX_TYPE & B, const math_vector< bigint > &b)",
				    DMESSAGE, EMESSAGE[1]);

	bigint *tmp = b.get_data_address();
	for (i = 0; i < B.rows; i++)
		A.value[i][B.columns].assign(-tmp[i]);
	kernel1(RES, A);

	// Step 2
	if (RES.is_column_zero(0) || RES.is_row_zero(B.columns)) {
		MATRIX_TYPE C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	tmp = new bigint[RES.columns];
	RES.get_row(tmp, B.columns);
	bigint *g = LiDIA::mgcd2(tmp, RES.columns);

	delete[] tmp;
	if (g[0] > 1) {
		MATRIX_TYPE C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	// Step 3, 4 
	bigint *x = (RES) * &(g[1]);
	delete[] g;

	// Step 5
	kernel1(RES, B);
	RES.set_no_of_columns(RES.columns + 1);
	for (i = 0; i < RES.rows; i++)
		RES.value[i][RES.columns-1].assign(x[i]);
	delete[] x;
}



//
// Smith normal form
//

template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_hartley (MATRIX_TYPE &RES) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_hartley()", DVALUE + 8);

	register lidia_size_t startr, startc, TEILBARKEIT;
	bigint TMP1, TMP2;
	bigint *tmp;
	register lidia_size_t xpivot, ypivot, i, j, z;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivot: first non zero

		xpivot = -1;
		ypivot = -1;
		for (i = startr; i < RES.rows; i++) {
			tmp = RES.value[i];
			for (j = startc; j < RES.columns; j++)
				if (!tmp[j].is_zero()) {
					xpivot = i;
					ypivot = j;
					i = RES.rows;
					j = RES.columns;
				}
		}

		if (xpivot != -1) {
			// swap to diagonalposition
			RES.swap_rows(startr, xpivot);
			RES.swap_columns(startc, ypivot);

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// mgcd computation for row
				for (i = startc + 1; i < RES.columns; i++)
					if (RES.value[startr][i] % RES.value[startr][startc] != 0) {
						div_rem(TMP1, TMP2, RES.value[startr][i], RES.value[startr][startc]);
						for (j = startr; j < RES.rows; j++) {
							tmp = RES.value[j];
							LiDIA::multiply(TMP2, tmp[startc], TMP1);
							LiDIA::subtract(tmp[i], tmp[i], TMP2);
						}
						RES.swap_columns(startc, i);
						i = startc;
					}

				// mgcd computation for column
				for (i = startr + 1; i < RES.rows; i++)
					if (RES.value[i][startc] % RES.value[startr][startc] != 0) {
						div_rem(TMP1, TMP2, RES.value[i][startc], RES.value[startr][startc]);
						for (j = startc; j < RES.columns; j++) {
							LiDIA::multiply(TMP2, RES.value[startr][j], TMP1);
							LiDIA::subtract(RES.value[i][j], RES.value[i][j], TMP2);
						}
						TEILBARKEIT = 0; //perhaps
						RES.swap_rows(i, startr);
						i = startr;
					}
			}

			// row elimination
			for (i = startc + 1; i < RES.columns; i++) {
				div_rem(TMP1, TMP2, RES.value[startr][i], RES.value[startr][startc]);
				for (j = startr; j < RES.rows; j++) {
					tmp = RES.value[j];
					LiDIA::multiply(TMP2, tmp[startc], TMP1);
					LiDIA::subtract(tmp[i], tmp[i], TMP2);
				}
			}

			// column elimination
			for (i = startr + 1; i < RES.rows; i++) {
				div_rem(TMP1, TMP2, RES.value[i][startc], RES.value[startr][startc]);
				for (j = startc; j < RES.columns; j++) {
					LiDIA::multiply(TMP2, RES.value[startr][j], TMP1);
					LiDIA::subtract(RES.value[i][j], RES.value[i][j], TMP2);
				}
			}

			// modulo test
			for (i = startr + 1; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						for (z = 0; z < RES.columns; z++)
							LiDIA::add(RES.value[startr][z], RES.value[startr][z], RES.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc = 0;
						startr = 0;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.columns && i < RES.rows; i++)
		if (RES.value[i][i].is_lt_zero())
			RES.value[i][i].negate();
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_hartley (MATRIX_TYPE &RES,
							 MATRIX_TYPE & T1,
							 MATRIX_TYPE & T2) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_hartley(MATRIX_TYPE &, MATRIX_TYPE &)", DVALUE + 8);

	register lidia_size_t startr, startc, TEILBARKEIT;
	bigint TMP1, TMP2;
	lidia_size_t xpivot, ypivot;
	register lidia_size_t i, j, z;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivot: first non zero
		xpivot = -1;
		ypivot = -1;
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++)
				if (RES.value[i][j] != 0) {
					xpivot = i;
					ypivot = j;
					i = RES.rows;
					j = RES.columns;
				}
		if (xpivot != -1) {
			// swap to diagonalposition
			RES.swap_rows(startr, xpivot);
			T1.swap_rows(startr, xpivot);
			RES.swap_columns(startc, ypivot);
			T2.swap_columns(startc, ypivot);

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// mgcd computation for row
				for (i = startc + 1; i < RES.columns; i++)
					if (RES.value[startr][i] % RES.value[startr][startc] != 0) {
						div_rem(TMP1, TMP2, RES.value[startr][i], RES.value[startr][startc]);
						for (j = startr; j < RES.rows; j++) {
							LiDIA::multiply(TMP2, RES.value[j][startc], TMP1);
							LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
						}
						for (j = 0; j < RES.columns; j++) {
							LiDIA::multiply(TMP2, T2.value[j][startc], TMP1);
							LiDIA::subtract(T2.value[j][i], T2.value[j][i], TMP2);
						}
						RES.swap_columns(startc, i);
						T2.swap_columns(startc, i);
						i = startc;
					}

				// mgcd computation for column
				for (i = startr + 1; i < RES.rows; i++)
					if (RES.value[i][startc] % RES.value[startr][startc] != 0) {
						div_rem(TMP1, TMP2, RES.value[i][startc], RES.value[startr][startc]);
						for (j = startc; j < RES.columns; j++) {
							LiDIA::multiply(TMP2, RES.value[startr][j], TMP1);
							LiDIA::subtract(RES.value[i][j], RES.value[i][j], TMP2);
						}
						for (j = 0; j < RES.rows; j++) {
							LiDIA::multiply(TMP2, T1.value[startr][j], TMP1);
							LiDIA::subtract(T1.value[i][j], T1.value[i][j], TMP2);
						}
						TEILBARKEIT = 0;
						RES.swap_rows(i, startr);
						T1.swap_rows(i, startr);
						i = startr;
					}
			}

			// row elimination
			for (i = startc + 1; i < RES.columns; i++) {
				div_rem(TMP1, TMP2, RES.value[startr][i], RES.value[startr][startc]);
				for (j = startr; j < RES.rows; j++) {
					LiDIA::multiply(TMP2, RES.value[j][startc], TMP1);
					LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
				}
				for (j = 0; j < RES.columns; j++) {
					LiDIA::multiply(TMP2, T2.value[j][startc], TMP1);
					LiDIA::subtract(T2.value[j][i], T2.value[j][i], TMP2);
				}
			}

			// column elimination
			for (i = startr + 1; i < RES.rows; i++) {
				div_rem(TMP1, TMP2, RES.value[i][startc], RES.value[startr][startc]);
				for (j = startc; j < RES.columns; j++) {
					LiDIA::multiply(TMP2, RES.value[startr][j], TMP1);
					LiDIA::subtract(RES.value[i][j], RES.value[i][j], TMP2);
				}
				for (j = 0; j < RES.rows; j++) {
					LiDIA::multiply(TMP2, T1.value[startr][j], TMP1);
					LiDIA::subtract(T1.value[i][j], T1.value[i][j], TMP2);
				}
			}

			// modulo test
			for (i = startr + 1; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						for (z = 0; z < RES.columns; z++)
							LiDIA::add(RES.value[startr][z], RES.value[startr][z], RES.value[i][z]);
						for (z = 0; z < RES.rows; z++)
							LiDIA::add(T1.value[startr][z], T1.value[startr][z], T1.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc = 0;
						startr = 0;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0) {
			RES.value[i][i].negate();
			for (z = 0; z < RES.columns; z++)
				T2.value[z][i].negate();
		}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_simple (MATRIX_TYPE &RES) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_simple()", DVALUE + 8);
	bigint PIVOT, TMP1, TMP2;
	bigint *tmp = NULL, *deltmp;

	MATRIX_TYPE TR1(RES.rows, RES.rows);
	MATRIX_TYPE TR2(RES.columns, RES.columns);

	bigint *REM;
	register lidia_size_t startr, startc, pivot, i, j, z, TEILBARKEIT;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivot: first non zero
		pivot = -1;
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++)
				if (RES.value[i][j] != 0) {
					pivot = i;
					i = RES.rows;
					j = RES.columns;
				}

		if (pivot != -1) {
			// swap pivot in actual row
			RES.swap_rows(startr, pivot);

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// mgcd computation and row elimination
				deltmp = new bigint[RES.columns];
				RES.get_row(deltmp, startr);
				REM = mgcd2(TR2, deltmp, RES.columns);
				delete[] deltmp;

				delete[] REM;
				TR2 = TR2.trans();
				LiDIA::multiply(RES, RES, TR2);

				tmp = RES.value[startr];
				for (i = 0; tmp[i].is_zero() && i < RES.columns; i++);
				RES.swap_columns(startc, i);

				// mgcd computation and column elimination
				deltmp = new bigint[RES.rows];
				RES.get_column(deltmp, startc);
				REM = mgcd2(TR1, deltmp, RES.rows);
				delete[] deltmp;

				delete[] REM;
				LiDIA::multiply(RES, TR1, RES);

				for (i = 0; RES.value[i][startc].is_zero() && i < RES.rows; i++);
				RES.swap_rows(startr, i);

				// control: row == 0
				tmp = RES.value[startr];
				for (i = startc+1; tmp[i].is_zero() && i < RES.columns; i++);
				if (i != RES.columns)
					TEILBARKEIT = 0;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr)
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0)
			RES.value[i][i].negate();
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_simple (MATRIX_TYPE &RES,
							MATRIX_TYPE & T1,
							MATRIX_TYPE & T2) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_simple(MATRIX_TYPE &, MATRIX_TYPE &)", DVALUE + 8);
	bigint PIVOT, TMP1, TMP2;
	bigint *tmp = NULL, *deltmp;

	MATRIX_TYPE TR1 = T1;
	MATRIX_TYPE TR2 = T2;
	bigint *REM;
	register lidia_size_t startr, startc, pivot, i, j, z, TEILBARKEIT;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivot: first non zero
		pivot = -1;
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++)
				if (RES.value[i][j] != 0) {
					pivot = i;
					i = RES.rows;
					j = RES.columns;
				}

		if (pivot != -1) {
			// swap pivot in actual row
			RES.swap_rows(startr, pivot);
			T1.swap_rows(startr, pivot);

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// mgcd computation and row elimination
				deltmp = new bigint[RES.columns];
				RES.get_row(deltmp, startr);
				REM = mgcd2(TR2, deltmp, RES.columns);
				delete[] deltmp;

				delete[] REM;
				TR2 = TR2.trans();
				LiDIA::multiply(RES, RES, TR2);
				LiDIA::multiply(T2, T2, TR2);

				tmp = RES.value[startr];
				for (i = 0; tmp[i].is_zero() && i < RES.columns; i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd computation and column elimination
				deltmp = new bigint[RES.rows];
				RES.get_column(deltmp, startc);
				REM = mgcd2(TR1, deltmp, RES.rows);
				delete[] deltmp;

				delete[] REM;
				LiDIA::multiply(RES, TR1, RES);
				LiDIA::multiply(T1, TR1, T1);

				for (i = 0; RES.value[i][startc].is_zero() && i < RES.rows; i++);
				RES.swap_rows(startr, i);
				T1.swap_rows(startr, i);

				// control: row == 0
				tmp = RES.value[startr];
				for (i = startc+1; tmp[i].is_zero() && i < RES.columns; i++);
				if (i != RES.columns)
					TEILBARKEIT = 0;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr) {
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
							for (z = 0; z < RES.rows; z++)
								LiDIA::add(T1.value[startr][z], T1.value[startr][z], T1.value[i][z]);
						}
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0) {
			RES.value[i][i].negate();
			for (z = 0; z < RES.columns; z++)
				T2.value[z][i].negate();
		}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_havas (MATRIX_TYPE &RES) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_havas()", DVALUE + 8);

//  std::cout << "snf_havas() " << std::endl;
	register lidia_size_t i, j, z, index;
	bigint PIVOT;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, SW, TEILBARKEIT;
	bigint TMP1, TMP2;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {
//      std::cout << " - " << startr << std::flush;
		// pivot: first non zero
		xpivot = -1;
		ypivot = -1;
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++)
				if (!RES.value[i][j].is_zero()) {
					xpivot = i;
					ypivot = j;
					i = RES.rows;
					j = RES.columns;
				}

		if (xpivot != -1) {
			// swap to actual row
			RES.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;

					for (i = 0; i < RES.columns; i++)
						if (i != index && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero() && i < RES.columns; i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr; // no index search
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero() && i < RES.rows; i++);
				RES.swap_rows(startr, i);

				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;
				index = startc;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr)
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0)
			RES.value[i][i].negate();
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_havas (MATRIX_TYPE &RES,
						       MATRIX_TYPE & T1,
						       MATRIX_TYPE & T2) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_havas(MATRIX_TYPE &, MATRIX_TYPE &)", DVALUE + 8);

//  std::cout << "snf_havas2() " << std::endl;
	register lidia_size_t i, j, z, index;
	bigint PIVOT;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, SW, TEILBARKEIT;
	bigint TMP1, TMP2;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {
//      std::cout << " - " << startc << std::endl;
		// pivot: first non zero
		xpivot = -1;
		ypivot = -1;
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++)
				if (!RES.value[i][j].is_zero()) {
					xpivot = i;
					ypivot = j;
					i = RES.rows;
					j = RES.columns;
				}

		if (xpivot != -1) {
			// swap to actual row
			RES.swap_rows(startr, xpivot);
			T1.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;

					for (i = 0; i < RES.columns; i++)
						if (i != index && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, T2.value[j][index], TMP1);
								LiDIA::subtract(T2.value[j][i], T2.value[j][i], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero() && i < RES.columns; i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr; // no index search
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, T1.value[index][j], TMP1);
								LiDIA::subtract(T1.value[i][j], T1.value[i][j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero() && i < RES.rows; i++);
				RES.swap_rows(startr, i);
				T1.swap_rows(startr, i);

				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;
				index = startc;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr) {
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
							for (z = 0; z < RES.rows; z++)
								LiDIA::add(T1.value[startr][z], T1.value[startr][z], T1.value[i][z]);
						}
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0) {
			RES.value[i][i].negate();
			for (z = 0; z < RES.columns; z++)
				T2.value[z][i].negate();
		}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_mult (MATRIX_TYPE &RES, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_mult(MATRIX_TYPE &, MATRIX_TYPE &)", DVALUE + 8);
	register lidia_size_t i, j, z, index, SW;
	bigint TMP1, TMP2;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, TEILBARKEIT;
	bigint ROW, COLUMN, PIVOT, NORM;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivotselection: minimale C * R norm
		xpivot = -1;
		ypivot = -1;
		PIVOT.assign_zero();
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++) {
				if (PIVOT == abs(RES.value[i][j])) {
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::multiply(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if ((PIVOT > abs(RES.value[i][j]) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::multiply(NORM, ROW, COLUMN);
					xpivot = i;
					ypivot = j;
				}
			}

		if (!PIVOT.is_zero()) {

			// swap to actual row
			RES.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((i != index) && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
						}

					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero(); i++);
				RES.swap_rows(startr, i);

				tmp = RES.value[startr];
				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;

				index = startr;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr)
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0)
			RES.value[i][i].negate();
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_mult (MATRIX_TYPE &RES,
						      MATRIX_TYPE & T1,
						      MATRIX_TYPE & T2,
						      long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_mult(MATRIX_TYPE &, MATRIX_TYPE &)", DVALUE + 8);

	register lidia_size_t i, j, z, index, SW;
	bigint TMP1, TMP2;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, TEILBARKEIT;
	bigint ROW, COLUMN, PIVOT, NORM;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivotselection: minimale C * R norm
		xpivot = -1;
		ypivot = -1;
		PIVOT.assign_zero();
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++) {
				if (PIVOT == abs(RES.value[i][j])) {
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::multiply(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if ((PIVOT > abs(RES.value[i][j]) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::multiply(NORM, ROW, COLUMN);
					xpivot = i;
					ypivot = j;
				}
			}

		if (!PIVOT.is_zero()) {

			// swap to actual row
			RES.swap_rows(startr, xpivot);
			T1.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((i != index) && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, T2.value[j][index], TMP1);
								LiDIA::subtract(T2.value[j][i], T2.value[j][i], TMP2);
							}
						}

					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, T1.value[index][j], TMP1);
								LiDIA::subtract(T1.value[i][j], T1.value[i][j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero(); i++);
				RES.swap_rows(startr, i);
				T1.swap_rows(startr, i);

				tmp = RES.value[startr];
				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;

				index = startr;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr) {
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
							for (z = 0; z < RES.rows; z++)
								LiDIA::add(T1.value[startr][z], T1.value[startr][z], T1.value[i][z]);
						}
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0) {
			RES.value[i][i].negate();
			for (z = 0; z < RES.columns; z++)
				T2.value[z][i].negate();
		}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_add (MATRIX_TYPE &RES, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_add(long)", DVALUE + 8);

	register lidia_size_t i, j, z, index, SW;
	bigint TMP1, TMP2;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, TEILBARKEIT;
	bigint ROW, COLUMN, PIVOT, NORM;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivotselection: minimale C * R norm
		xpivot = -1;
		ypivot = -1;
		PIVOT.assign_zero();
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++) {
				if (PIVOT == abs(RES.value[i][j])) {
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::add(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if ((PIVOT > abs(RES.value[i][j]) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::add(NORM, ROW, COLUMN);
					xpivot = i;
					ypivot = j;
				}
			}

		if (!PIVOT.is_zero()) {

			// swap to actual row
			RES.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((i != index) && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
						}

					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero(); i++);
				RES.swap_rows(startr, i);

				tmp = RES.value[startr];
				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;

				index = startr;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr)
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0)
			RES.value[i][i].negate();
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_add (MATRIX_TYPE &RES,
						     MATRIX_TYPE & T1,
						     MATRIX_TYPE & T2,
						     long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_add(MATRIX_TYPE &, MATRIX_TYPE &, long)", DVALUE + 8);

	register lidia_size_t i, j, z, index, SW;
	bigint TMP1, TMP2;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, TEILBARKEIT;
	bigint ROW, COLUMN, PIVOT, NORM;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {

		// pivotselection: minimale C * R norm
		xpivot = -1;
		ypivot = -1;
		PIVOT.assign_zero();
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++) {
				if (PIVOT == abs(RES.value[i][j])) {
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::add(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if ((PIVOT > abs(RES.value[i][j]) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					row_norm(RES, ROW, i, art);
					column_norm(RES, COLUMN, j, art);
					LiDIA::add(NORM, ROW, COLUMN);
					xpivot = i;
					ypivot = j;
				}
			}

		if (!PIVOT.is_zero()) {

			// swap to actual row
			RES.swap_rows(startr, xpivot);
			T1.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((i != index) && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, T2.value[j][index], TMP1);
								LiDIA::subtract(T2.value[j][i], T2.value[j][i], TMP2);
							}
						}

					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, T1.value[index][j], TMP1);
								LiDIA::subtract(T1.value[i][j], T1.value[i][j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero(); i++);
				RES.swap_rows(startr, i);
				T1.swap_rows(startr, i);

				tmp = RES.value[startr];
				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;

				index = startr;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr) {
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
							for (z = 0; z < RES.rows; z++)
								LiDIA::add(T1.value[startr][z], T1.value[startr][z], T1.value[i][z]);
						}
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0) {
			RES.value[i][i].negate();
			for (z = 0; z < RES.columns; z++)
				T2.value[z][i].negate();
		}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_new (MATRIX_TYPE &RES, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_new(long)", DVALUE + 8);

	register lidia_size_t i, j, z, index, SW;
	bigint TMP1, TMP2;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, TEILBARKEIT;
	bigint ROW, COLUMN, PIVOT, NORM;

	bigint *RO = new bigint[RES.rows];
	bigint *CO = new bigint[RES.columns];

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {
		// norm computation
		for (i = 0; i < RES.rows; i++)
			row_norm(RES, RO[i], i, art);
		for (i = 0; i < RES.columns; i++)
			column_norm(RES, CO[i], i, art);

		// pivotselection: new minimale C * R norm
		xpivot = -1;
		ypivot = -1;
		PIVOT.assign_zero();
		NORM.assign_zero();
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++) {
				LiDIA::multiply(TMP1, RO[i], CO[j]);

				if (!RES.value[i][j].is_zero() && (NORM > TMP1 || PIVOT.is_zero())) {
					NORM.assign(TMP1);
					PIVOT.assign(abs(RES.value[i][j]));
					xpivot = i;
					ypivot = j;
				}

				if (NORM == TMP1 && !PIVOT.is_zero() && !RES.value[i][j].is_zero()) {
					if (PIVOT > abs(RES.value[i][j])) {
						PIVOT.assign(abs(RES.value[i][j]));
						NORM.assign(TMP1);
						xpivot = i;
						ypivot = j;
					}
				}
			}

		if (!PIVOT.is_zero()) {

			// swap to actual row
			RES.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((i != index) && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
						}

					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero(); i++);
				RES.swap_rows(startr, i);

				tmp = RES.value[startr];
				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;

				index = startr;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr)
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0)
			RES.value[i][i].negate();
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snf_new (MATRIX_TYPE &RES,
						     MATRIX_TYPE & T1,
						     MATRIX_TYPE & T2,
						     long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_new(MATRIX_TYPE &, MATRIX_TYPE &, long)", DVALUE + 8);
	register lidia_size_t i, j, z, index, SW;
	bigint TMP1, TMP2;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, TEILBARKEIT;
	bigint ROW, COLUMN, PIVOT, NORM;

	bigint *RO = new bigint[RES.rows];
	bigint *CO = new bigint[RES.columns];

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {
		// norm computation
		for (i = 0; i < RES.rows; i++)
			row_norm(RES, RO[i], i, art);
		for (i = 0; i < RES.columns; i++)
			column_norm(RES, CO[i], i, art);

		// pivotselection: new minimale C * R norm
		xpivot = -1;
		ypivot = -1;
		PIVOT.assign_zero();
		NORM.assign_zero();
		for (i = startr; i < RES.rows; i++)
			for (j = startc; j < RES.columns; j++) {
				//row_norm(ROW, i, art);
				//column_norm(RES, COLUMN, j, art);
				LiDIA::multiply(TMP1, RO[i], CO[j]);

				if (!RES.value[i][j].is_zero() && (NORM > TMP1 || PIVOT.is_zero())) {
					NORM.assign(TMP1);
					PIVOT.assign(abs(RES.value[i][j]));
					xpivot = i;
					ypivot = j;
				}

				if (NORM == TMP1 && !PIVOT.is_zero() && !RES.value[i][j].is_zero()) {
					if (PIVOT > abs(RES.value[i][j])) {
						PIVOT.assign(abs(RES.value[i][j]));
						NORM.assign(TMP1);
						xpivot = i;
						ypivot = j;
					}
				}
			}

		if (!PIVOT.is_zero()) {

			// swap to actual row
			RES.swap_rows(startr, xpivot);
			T1.swap_rows(startr, xpivot);

			index = ypivot;

			TEILBARKEIT = 0;

			while (TEILBARKEIT == 0) {
				TEILBARKEIT = 1;

				// gcd2(row(startr), columns, TR2);
				tmp = RES.value[startr];
				do {
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((i != index) && !tmp[i].is_zero()) {
							SW = 1;
							div_rem(TMP1, TMP2, tmp[i], tmp[index]);
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, RES.value[j][index], TMP1);
								LiDIA::subtract(RES.value[j][i], RES.value[j][i], TMP2);
							}
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, T2.value[j][index], TMP1);
								LiDIA::subtract(T2.value[j][i], T2.value[j][i], TMP2);
							}
						}

					for (i = 0; i < RES.columns; i++)
						if ((abs(tmp[index]) > abs(tmp[i])) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do {
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs(RES.value[index][startc]) > abs(RES.value[i][startc])) && !RES.value[i][startc].is_zero())
							index = i;

					for (i = 0; i < RES.rows; i++)
						if ((i != index) && !RES.value[i][startc].is_zero()) {
							SW = 1;
							tmp = RES.value[i];
							div_rem(TMP1, TMP2, tmp[startc], RES.value[index][startc]);
							for (j = 0; j < RES.columns; j++) {
								LiDIA::multiply(TMP2, RES.value[index][j], TMP1);
								LiDIA::subtract(tmp[j], tmp[j], TMP2);
							}
							for (j = 0; j < RES.rows; j++) {
								LiDIA::multiply(TMP2, T1.value[index][j], TMP1);
								LiDIA::subtract(T1.value[i][j], T1.value[i][j], TMP2);
							}
						}
				}
				while (SW == 1);

				for (i = 0; RES.value[i][startc].is_zero(); i++);
				RES.swap_rows(startr, i);
				T1.swap_rows(startr, i);

				tmp = RES.value[startr];
				for (index = startc+1; index < RES.columns && tmp[index].is_zero(); index++);
				if (index != RES.columns)
					TEILBARKEIT = 0;

				index = startr;
			}

			// modulo test
			for (i = startr; i < RES.rows; i++)
				for (j = startc + 1; j < RES.columns; j++)
					if (RES.value[i][j] % RES.value[startr][startc] != 0) {
						if (i != startr) {
							for (z = 0; z < RES.columns; z++)
								LiDIA::add(tmp[z], tmp[z], RES.value[i][z]);
							for (z = 0; z < RES.rows; z++)
								LiDIA::add(T1.value[startr][z], T1.value[startr][z], T1.value[i][z]);
						}
						i = RES.rows;
						j = RES.columns;
						startc--;
						startr--;
					}
		}
	}

	// diagonal >= 0
	for (i = 0; i < RES.rows && i < RES.columns; i++)
		if (RES.value[i][i] < 0) {
			RES.value[i][i].negate();
			for (z = 0; z < RES.columns; z++)
				T2.value[z][i].negate();
		}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snfmod_dkt (MATRIX_TYPE &RES,
							const bigint &mod) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snfmod_dkt(const bigint &)", DVALUE + 8);

	register lidia_size_t diagindex, j, z, l;

	bigint RES0, RES1, RES2, RES3; // 0 = lggT, 1 = rggt, 2 = ggt
	bigint x, y;
	bigint TMP, TMP1, TMP2, TMP3;
	bigint *Atmp, *Atmp1 = NULL;

	// A = A % mod
	for (z = 0; z < RES.rows; z++) {
		Atmp = RES.value[z];
		for (j = 0; j < RES.columns; j++)
			best_remainder(Atmp[j], Atmp[j], mod);
	}

	// Step 3 - 5
	for (diagindex = 0; diagindex < RES.rows; diagindex++) {
		Atmp = RES.value[diagindex];

		// if diagonalelement == 0 then diagonalelement = mod
		if (Atmp[diagindex].is_zero())
			Atmp[diagindex].assign(mod);

		// Step 6 -8
		for (j = diagindex+1; j < RES.columns; j++)
			if (!Atmp[j].is_zero()) {
				// columns reduction
				RES2 = xgcd(RES0, RES1, Atmp[j], Atmp[diagindex]);
				div_rem(x, RES3, Atmp[diagindex], RES2);
				div_rem(y, RES3, Atmp[j], RES2);

				for (z = 0; z < RES.rows; z++) {
					Atmp1 = RES.value[z];
					TMP.assign(Atmp1[j]);
					TMP1.assign(Atmp1[diagindex]);

					mult_mod(TMP2, TMP, x, mod);
					mult_mod(TMP3, TMP1, y, mod);
					sub_mod(Atmp1[j], TMP2, TMP3, mod);

					mult_mod(TMP2, TMP, RES0, mod);
					mult_mod(TMP3, TMP1, RES1, mod);
					add_mod(Atmp1[diagindex], TMP2, TMP3, mod);
				}
			}

		for (j = diagindex+1; j < RES.rows; j++)
			if (!RES.value[j][diagindex].is_zero()) {
				// row reduction
				RES2 = xgcd(RES0, RES1, RES.value[j][diagindex], Atmp[diagindex]);
				div_rem(x, RES3, Atmp[diagindex], RES2);
				div_rem(y, RES3, RES.value[j][diagindex], RES2);

				for (z = 0; z < RES.columns; z++) {
					TMP.assign(RES.value[j][z]);
					TMP1.assign(RES.value[diagindex][z]);

					mult_mod(TMP2, TMP, x, mod);
					mult_mod(TMP3, TMP1, y, mod);
					sub_mod(RES.value[j][z], TMP2, TMP3, mod);

					mult_mod(TMP2, TMP, RES0, mod);
					mult_mod(TMP3, TMP1, RES1, mod);
					add_mod(RES.value[diagindex][z], TMP2, TMP3, mod);
				}
			}

		// value[diagindex][diagindex] | value[i][j] ???
		TMP = Atmp[diagindex];
		for (j = diagindex+1; j < RES.rows; j++)
			for (z = diagindex+1; z < RES.columns; z++) {
				if (RES.value[j][z] % TMP != 0) {
					if (j != diagindex)
						for (l = diagindex; l < RES.columns; l++)
							add_mod(Atmp[l], Atmp[l], RES.value[j][l], mod);
					j = RES.rows;
					z = RES.columns;
				}
			}

		for (z = diagindex+1; z < RES.columns && Atmp[z].is_zero(); z++);
		if (z != RES.columns)
			diagindex--;
	}

	// Step 29 - 32
	bigint D = mod;
	for (j = 0; j < RES.rows; j++) {
		Atmp = RES.value[j];
		Atmp[j].assign(xgcd(RES0, RES1, Atmp[j], D));
		div_rem(D, TMP, D, Atmp[j]);
	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::snfmod_cohen (MATRIX_TYPE &RES,
							  const bigint & mod) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snfmod_cohen(const bigint &)", DVALUE + 8);

	register lidia_size_t diagindex, j, z, l;

	bigint RES0, RES1, RES2, RES3; // 0 = lggT, 1 = rggt, 2 = ggt
	bigint x, y;
	bigint TMP, TMP1, TMP2, TMP3;
	bigint *Atmp, *Atmp1 = NULL;
	bigint D = mod;

	// Step 1
	for (diagindex = 0; diagindex < RES.rows; diagindex++) {
		Atmp = RES.value[diagindex];

		if (Atmp[diagindex].is_zero())
			Atmp[diagindex].assign(mod);

		// Step 2 - 4
		for (j = diagindex+1; j < RES.columns; j++)
			if (!Atmp[j].is_zero()) {
				// columns reduction
				RES2 = xgcd(RES0, RES1, Atmp[j], Atmp[diagindex]);
				div_rem(x, RES3, Atmp[diagindex], RES2);
				div_rem(y, RES3, Atmp[j], RES2);

				for (z = 0; z < RES.rows; z++) {
					Atmp1 = RES.value[z];
					TMP.assign(Atmp1[j]);
					TMP1.assign(Atmp1[diagindex]);

					mult_mod(TMP2, TMP, x, mod);
					mult_mod(TMP3, TMP1, y, mod);
					sub_mod(Atmp1[j], TMP2, TMP3, mod);

					mult_mod(TMP2, TMP, RES0, mod);
					mult_mod(TMP3, TMP1, RES1, mod);
					add_mod(Atmp1[diagindex], TMP2, TMP3, mod);
				}
			}

		// Step 5 - 7
		for (j = diagindex+1; j < RES.rows; j++)
			if (!RES.value[j][diagindex].is_zero()) {
				// row reduction
				RES2 = xgcd(RES0, RES1, RES.value[j][diagindex], Atmp[diagindex]);
				div_rem(x, RES3, Atmp[diagindex], RES2);
				div_rem(y, RES3, RES.value[j][diagindex], RES2);

				for (z = 0; z < RES.columns; z++) {
					TMP.assign(RES.value[j][z]);
					TMP1.assign(RES.value[diagindex][z]);

					mult_mod(TMP2, TMP, x, mod);
					mult_mod(TMP3, TMP1, y, mod);
					sub_mod(RES.value[j][z], TMP2, TMP3, mod);

					mult_mod(TMP2, TMP, RES0, mod);
					mult_mod(TMP3, TMP1, RES1, mod);
					add_mod(RES.value[diagindex][z], TMP2, TMP3, mod);
				}
			}

		// Step 8, 9
		TMP = Atmp[diagindex];

		for (j = diagindex+1; j < RES.rows; j++)
			for (z = diagindex+1; z < RES.columns; z++) {
				if (RES.value[j][z] % TMP != 0) {
					if (j != diagindex)
						for (l = diagindex; l < RES.columns; l++)
							add_mod(Atmp[l], Atmp[l], RES.value[j][l], mod);
					j = RES.rows;
					z = RES.columns;
				}
			}

		for (z = diagindex+1; z < RES.columns && Atmp[z].is_zero(); z++);
		if (z != RES.columns)
			diagindex--;
		else {
			// Step 10
			Atmp[diagindex] = xgcd(TMP1, TMP2, TMP, D);
			div_rem(D, TMP1, D, Atmp[diagindex]);
		}
	}
}



//
// END: Linear algebra
// PART 2
//

template< class MATRIX_TYPE >
inline void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::gauss (MATRIX_TYPE &RES) const
{
	debug_handler(DMESSAGE, "in member - function gauss()");

	MATRIX_TYPE TR(RES.columns, RES.columns);
	bigint *REM = NULL;
	register lidia_size_t startr = 0, startc = 0, i;

	for (startc = RES.columns - 1, startr = RES.rows - 1; startr >= 0 && startc >= 0; startr--, startc--) {

		bigint *ZU = new bigint[RES.columns];
		RES.get_row(ZU, startr);
		for (i = startc + 1; i < RES.columns; ZU[i].assign_zero(), i++);
		REM = TR.mgcd1(ZU, RES.columns);
		delete[] REM;
		delete[] ZU;

		TR = TR.trans();
		LiDIA::multiply(RES, RES, TR);
		for (i = 0; RES.value[startr][i].is_zero() && i <= startc; i++);
		if (i > startc)
			return;
		RES.swap_columns(startc, i);

	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::basis_completion (MATRIX_TYPE &M, bigint *x,
							      lidia_size_t n) const
{
	bigint *y = new bigint[n+1];

	lidia_size_t i, j;
	for (i = 1; i <= n; i++) {
		y[i] = x[i-1];
	}
	y[0] = -1;

	if (M.rows != n)
		M.set_no_of_rows(n);
	if (M.columns != n+1)
		M.set_no_of_columns(n+1);

	// diag
	for (i = 0; i < n; i++) {
		M.sto(i, i+1, 1);
		M.sto(i, 0, x[i]);
	}

	bigint a, b, g, c, d, gold = y[1], TMP;

	for (i = 1; i < n - 1; i++) {
		if (y[i] == 1) {
			for (j = i; j < n; j++)
				M.swap_columns(j, j+1);
			M.columns--;
			return;
		}

		if (y[i + 1] == 0) {
			continue;
		}
		g = xgcd(a, b, y[i + 1], gold);
		b.negate();
		c = gold / g;
		d = y[i+1] /g;

		for (j = 0; j < n; j++) {

			TMP = M.value[j][i];
			M.sto(j, i, TMP * a + M(j, i+1) * b);
			M.sto(j, i+1, TMP * c + M(j, i+1) * d);
		}
		y[i+1] = g;
		y[i] = 0;

		gold = g;
	}
}



template< class MATRIX_TYPE >
void
sparse_bigint_matrix_kernel< MATRIX_TYPE >::simple_basis_completion (MATRIX_TYPE &M, bigint *x,
								     lidia_size_t n) const
{
	bigint TMP1, TMP2, g;

	if (M.rows != n)
		M.set_no_of_rows(n);
	if (M.columns != n)
		M.set_no_of_columns(n);

	M.diag(1, 0);

	lidia_size_t i, pos;
	for (pos = 0; pos < n && x[pos].is_zero(); pos++);

	//  std::cout << "Position = " << pos << std::endl;

	if (pos + 1 == n) {
		if (x[pos] < 0)
			x[pos].negate();
		return;
	}

	g = xgcd(TMP1, TMP2, x[pos], x[pos + 1]);
	for (i = pos; i < n; i++)
		M.sto(i, pos, x[i]);
	M.sto(pos, pos + 1, -TMP2);
	M.sto(pos+1, pos+1, TMP1);
}



template< class MATRIX_TYPE >
lidia_size_t
sparse_bigint_matrix_kernel< MATRIX_TYPE >::cond_matrix (MATRIX_TYPE &M,
							 bigint *v,
							 lidia_size_t n) const
{
	if (M.rows != n)
		M.set_no_of_rows(n);
	if (M.columns != n)
		M.set_no_of_columns(n);

	lidia_size_t i, pos;

	M.diag(1, 0);

	bigint *b = new bigint[n];
	for (i = 0; i < n; i++)
		b[i] = v[i];

	for (pos = 0; pos < n && b[pos].is_zero(); pos++);

	bigint g, N = b[pos], a, qa, ra, qb, rb, t = 0;
	if (pos + 1 < n)
		a = b[pos + 1];

	for (i = pos + 2; i < n && !g.is_one(); i++)
		if (!b[i].is_zero())
			if (!a.is_zero()) {
				g = gcd(b[i], a);
				div_rem(qa, ra, a/g, N);
				div_rem(qb, rb, b[i]/g, N);

				for (t = 0; gcd((ra + t*rb), N) != 1; t++);

				a = a + t * b[i]; // NEW

				M.sto(pos + 1, i, t);
			}
			else {
				a = a + b[i]; // NEW
				M.sto(pos + 1, i, bigint(1));
			}

	delete[] b;

	return pos;
}



//
// mgcd2
//

template< class MATRIX_TYPE >
bigint *
sparse_bigint_matrix_kernel< MATRIX_TYPE >::mgcd2 (MATRIX_TYPE &RES,
						   const bigint * aconst,
						   lidia_size_t n) const
{
	//
	// DESCRIPTION: RES = T.mgcd2(a, n);
	// =  > RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	// =  > RES[0] = gcd(a[0], ..., a[n-1])
	// =  > T*a = RES
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements, MIN assignments
	// PAPER: Hermite normal form computation for integer matrices, Havas
	// VERSION: 1.8
	//

	debug_handler("multiple_gcd", "in member - function "
		      "mgcd2(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, bound, SW;
	bigint MIN, TMP, q, r, *Ttmp1, *Ttmp2 = NULL;

	if (RES.columns != n)
		RES.set_no_of_columns(n);
	if (RES.rows != n)
		RES.set_no_of_rows(n);
	RES.diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd2 :: "
		       "Error in memory allocation (a)");

	for (i = 0; i < n; i++)
		a[i].assign(aconst[i]);

	// init
	for (index = 0; index < n && a[index].is_zero(); index++);

	if (index == n) {
		delete[] a;
		return new bigint[n];
	}
	else
		bound = index;

	do {
		MIN.assign(a[index]);

		// Pivot search: MINIMUM
		for (i = bound; i < n; i++)
			if ((abs(MIN) > abs(a[i])) && !a[i].is_zero()) {
				MIN.assign(a[i]);
				index = i;
			}

		// all elements
		SW = 0;
		Ttmp2 = RES.value[index];
		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW = 1;
				Ttmp1 = RES.value[i];
				div_rem(q, r, a[i], MIN);
				a[i].assign(r);
				for (j = 0; j < n; j++) {
					LiDIA::multiply(TMP, q, Ttmp2[j]);
					LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
				}
			}
	}
	while (SW == 1);

	Ttmp2 = RES.value[index];

	// gcd < 0 ?
	if (a[index] < 0) {
		a[index].negate();
		for (i = 0; i < n; i++)
			Ttmp2[i].negate();
	}

	if (index != 0)
		a[0].assign(a[index]);
	for (i = 1; i <= n; i++)
		a[i].assign(Ttmp2[i - 1]);

	return a;
}



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_CC_GUARD_
