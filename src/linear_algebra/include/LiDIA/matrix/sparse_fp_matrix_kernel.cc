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


#ifndef LIDIA_SPARSE_FP_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_SPARSE_FP_MATRIX_KERNEL_CC_GUARD_


#include	"LiDIA/modular_operations.inl"
#ifndef LIDIA_SPARSE_FP_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_fp_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// column step form
//

template< class T, class MATRIX_TYPE >
int
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::STF (MATRIX_TYPE &A, const T & mod) const
{
	//
	// Input: **value = values of matrix
	//              r = number of rows
	//              c = number of columns
	//            mod = modulus for Fp - class
	// Task: ex = STF(value, r, c, mod);
	// =  > matrix (value, r, c) in column step form
	// =  > ex = (-1)^Number_of_column_exchanges
	// Version: bigint 1.9
	//



	register lidia_size_t index = 0, i = 0, j = 0;

	lidia_size_t startr = A.rows;
	lidia_size_t startc = A.columns;

	T TMP, TMP1, TMP2;

	// Step 1 - 4
	register int exchange = 1;

	// Step 5 - 8
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		// Step 9 - 13
		for (index = startc; index >= 0 && S_base_modul.member(A, startr, index) == A.Zero; index--);

		// Step 14 - 26
		if (index != -1) {
			if (index != startc) {
				exchange = -exchange;

				// exchange column j0 with column j
				S_base_modul.swap_columns(A, startc, index);
			}

			inv_mod(TMP, S_base_modul.member(A, startr, startc), mod);

			// Step 27 - 29
			for (j = startc - 1; j >= 0; j--) {
				if (S_base_modul.member(A, startr, j) != A.Zero) {
					// Step 30 - 40
					mult_mod(TMP1, S_base_modul.member(A, startr, j), TMP, mod);

					for (i = 0; i <= startr; i++) {
						mult_mod(TMP2, S_base_modul.member(A, i, startc), TMP1, mod);
						sub_mod(TMP2, S_base_modul.member(A, i, j), TMP2, mod);
						S_base_modul.sto(A, i, j, TMP2);
					}
				}
			}
		}
		else {
			startc++;

		}
	}
	return exchange;
}



//
// extended column step form
//

template< class T, class MATRIX_TYPE >
const T *
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::STF_extended (MATRIX_TYPE &A, const T & mod) const
{
	//
	// Input: **value = values of matrix
	//              r = number of rows
	//              c = number of columns
	//            mod = modulus for Fp - class
	// Task: ex = STF_extended(A, mod);
	// =  > matrix A in column step form
	// =  > ex = (-1)^Number_of_column_exchanges
	// Version: bigint 1.9
	//

	register lidia_size_t index = 0, i = 0, j = 0;
	lidia_size_t startr = A.rows;
	lidia_size_t startc = A.columns;

	lidia_size_t line = A.columns - A.rows;

	T *RES = new T[line + 1];
	memory_handler(RES, "STF_extended", "Error in memory allocation (RES)");

	for (i = 0; i <= line; i++)
		RES[i] = 1;

	T TMP, TMP1, TMP2;
	bool SW = true;

	// Step 1 - 4
	register int exchange = 1;

	// Step 5 - 8
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		// Step 9 - 13
		for (index = startc; index > line && S_base_modul.member(A, startr, index) == A.Zero; index--);

		// Step 14 - 26
		if (index > line) {
			if (index != startc) {
				exchange = -exchange;

				// exchange column j0 with column j
				S_base_modul.swap_columns(A, startc, index);
			}

			for (i = 0; i <= line; i++)
				LiDIA::mult_mod(RES[i], RES[i], S_base_modul.member(A, startr, startc), mod);

			inv_mod(TMP, S_base_modul.member(A, startr, startc), mod);

			// Step 27 - 29
			for (j = startc - 1; j >= 0; j--) {
				if (S_base_modul.member(A, startr, j) != A.Zero) {
					// Step 30 - 40
					mult_mod(TMP1, S_base_modul.member(A, startr, j), TMP, mod);

					for (i = 0; i < startr; i++) {
						mult_mod(TMP2, S_base_modul.member(A, i, startc), TMP1, mod);
						sub_mod(TMP2, S_base_modul.member(A, i, j), TMP2, mod);
						S_base_modul.sto(A, i, j, TMP2);
					}
				}
			}
		}
		else {
			if (SW) {
				for (i = 0; i <= line; i++)
					LiDIA::mult_mod(RES[i], RES[i], S_base_modul.member(A, startr, i), mod);
				exchange = -exchange;
				SW = false;
			}
			else {
				for (i = 0; i <= line; i++)
					RES[i] = 0;
			}
			startc++;
		}
	}

	if (exchange == -1)
		for (i = 0; i <= line; i++)
			RES[i] = -RES[i];
	return RES;
}



//
// rank
//

template< class T, class MATRIX_TYPE >
lidia_size_t
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::rank (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: **Avalue = values of matrix
	//               r = number of rows
	//               c = number of columns
	//            Gmod = modulus for Fp - class
	// Task: rank_intern(Avalue, r, c, Gmod) = rank of matrix (Avalue, r, c)
	// Version: 1.9
	//

	debug_handler_l("bigint_matrix",
			"rank(MATRIX_TYPE &A, const T &)", LDBL_MATRIX);

	register lidia_size_t i, j, No = 0;

	// Step 1, 2
	STF(A, mod);

	// Step 3 - 12
	for (j = 0; j < A.columns; j++) {
		for (i = A.rows - 1; i >= 0 && S_base_modul.member(A, i, j) == A.Zero; i--);
		if (i == -1)
			No++;
	}

	// Step 13 - 24
	return (A.columns - No);
}



//
// rank and linearly independent rows or columns
//

template< class T, class MATRIX_TYPE >
lidia_size_t *
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::lininr (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: **Avalue = values of matrix
	//               r = number of rows
	//               c = number of columns
	//            Gmod = modulus for Fp - class
	// Task: IND = lininr_intern(Avalue, r, c, Gmod);
	// =  > IND[0] = Rank of matrix (Avalue, r, c)
	// =  > IND[1], ..., IND[IND[0]], such that
	//                 row(IND[1]), ..., row(IND[IND[0]])
	//                 of matrix (Avalue, r, c) are linearly independent.
	// Version: 1.9
	//

	debug_handler("bigint_matrix",
		      "lininr(MATRIX_TYPE &A, const T &)");

	register lidia_size_t i, j;
	lidia_size_t *l = new lidia_size_t[A.columns];
	for (i = 0; i < A.columns; l[i] = 0, i++);

	STF(A, mod);

	// Step 3 - 12
	for (j = 0; j < A.columns; j++) {
		for (i = A.rows - 1; i >= 0 && S_base_modul.member(A, i, j) == A.Zero; i--);
		l[j] = i;
	}

	// Step 13 - 24
	for (i = 0; i < A.columns && l[i] == -1; i++); // i = number of zero-columns

	lidia_size_t TMP = A.columns - i; // rank
	lidia_size_t *IND = new lidia_size_t[TMP + 1];
	memory_handler(IND, "bigint_matrix", "lininr_intern - "
		       "Error in memory allocation (IND)");
	IND[0] = TMP; // rank
	for (j = 0; j < TMP; j++)
		IND[TMP - j] = l[j + i];
	delete[] l;
	return IND;
}



template< class T, class MATRIX_TYPE >
lidia_size_t *
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::lininc (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: **Avalue = values of matrix
	//               r = number of rows
	//               c = number of columns
	//            Gmod = modulus for Fp - class
	// Task: IND = lininc_intern(Avalue, r, c, Gmod);
	// =  > IND[0] = Rank of matrix (Avalue, r, c)
	// =  > IND[1], ..., IND[IND[0]], such that
	//                 column(IND[1]), ..., column(IND[IND[0]])
	//                 of matrix (Avalue, r, c) are linearly independent.
	// Version: 1.9
	//

	register lidia_size_t i, j;
	lidia_size_t *l = new lidia_size_t[A.columns + 1];
	memory_handler(l, "bigint_matrix", "lininc_intern :: "
		       "Error in memory allocation (l)");

	// Step 1, 2
	STF(A, mod);

	// Step 3 - 12
	for (j = 0; j < A.columns; j++) {
		for (i = A.rows - 1; i >= 0 && S_base_modul.member(A, i, j) == A.Zero; i--);
		l[j] = i;
	}

	// Step 13 - 24
	for (i = 0; i < A.columns && l[i] == -1; i++); // i = number of zero-columns

	lidia_size_t TMP = A.columns - i;
	lidia_size_t *IND = new lidia_size_t[TMP + 1];
	memory_handler(IND, "bigint_matrix", "lininc_intern - "
		       "Error in memory allocation (IND)");
	IND[0] = TMP; // rank
	for (j = 0; j < TMP; j++)
		IND[TMP - j] = l[j + i];
	delete[] l;
	return IND;
}



//
// adjoint matrix
//

template< class T, class MATRIX_TYPE >
void
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::adj (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// Task: adj_intern(value, r, mod);
	// =  > adjoint matrix (value, r, r)
	// Version: 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "adj(MATRIX_TYPE &, cont T &)");

	register lidia_size_t i, j, z;
	T TMP, TMP1, TMP2;
	T *tmp, *tmp1, *Btmp, *Btmp1;
	lidia_size_t exchange = 1;
	T DET = 1;

	// Step 1, 2
	T **Bvalue = new T *[A.rows];
	memory_handler(Bvalue, "bigint_matrix", "adj_intern - T part :: "
		       "Error in memory allocation (Bvalue)");
	for (i = 0; i < A.rows; i++) {
		Btmp = new T[A.rows];
		tmp = A.value[i];
		memory_handler(Btmp, "bigint_matrix", "adj_intern - T part :: "
			       "Error in memory allocation (Btmp)");
		for (j = 0; j < A.rows; j++) {
			Btmp[j] = tmp[j];
			if (i == j)
				tmp[j] = 1;
			else
				tmp[j] = 0;
		}
		Bvalue[i] = Btmp;
	}

	// Step 3 - 5
	for (i = A.rows - 1; i >= 0; i--) {
		Btmp = Bvalue[i];

		// Step 6 - 11
		for (j = i; j >= 0 && Btmp[j] == 0; j--);

		// Step 12 - 19
		if (j != i) {
			exchange = -exchange;
			for (z = 0; z < A.rows; z++) {
				Btmp1 = Bvalue[z];
				tmp1 = A.value[z];

				// A.swap_columns(i, j);
				TMP = Btmp1[j];
				Btmp1[j] = Btmp1[i];
				Btmp1[i] = TMP;

				// B.swap_columns(i, j);
				TMP = tmp1[i];
				tmp1[i] = tmp1[j];
				tmp1[j] = TMP;
			}
		}

		inv_mod(TMP1, Btmp[i], mod);

		// Step 20 - 32
		for (j = 0; j < A.rows; j++) {
			if (j != i) {
				mult_mod(TMP2, Btmp[j], TMP1, mod);
				for (z = 0; z < i; z++) {
					Btmp1 = Bvalue[z];

					mult_mod(TMP, Btmp1[i], TMP2, mod);
					sub_mod(Btmp1[j], Btmp1[j], TMP, mod);
				}

				for (z = 0; z < A.rows; z++) {
					tmp1 = A.value[z];

					mult_mod(TMP, tmp1[i], TMP2, mod);
					sub_mod(tmp1[j], tmp1[j], TMP, mod);
				}
			}
		}
		mult_mod(DET, DET, Btmp[i], mod);
		for (z = 0; z < i; z++) {
			Btmp1 = Bvalue[z];
			mult_mod(Btmp1[i], Btmp1[i], TMP1, mod);
		}
		for (z = 0; z < A.rows; z++) {
			tmp1 = A.value[z];
			mult_mod(tmp1[i], tmp1[i], TMP1, mod);
		}
	}

	// Step 33 - 43
	if (exchange < 0)
		DET = -DET;
	for (j = 0; j < A.rows; j++) {
		tmp = A.value[j];
		for (i = 0; i < A.rows; i++)
			mult_mod(tmp[i], tmp[i], DET, mod);
	}

	for (j = 0; j < A.rows; j++)
		delete[] Bvalue[j];
	delete[] Bvalue;
}



//
// determinant
//

template< class T, class MATRIX_TYPE >
const T
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::det (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// Task: det_intern(value, r, mod) = determinant of matrix (value, r, r);
	// Version: T 1.9
	//

	debug_handler("bigint_matrix", "det(MATRIX_TYPE &, T)");

	register lidia_size_t i, j, z;
	T TMP, TMP1, TMP2;

	// Step 1 - 4
	lidia_size_t ex = 1;
	T ret = 1;

	// Step 5 - 8
	for (i = 0; i < A.rows; i++) {

		// Step 9 - 13
		for (j = i; j < A.rows && S_base_modul.member(A, j, i) == 0; j++);

		// Step 14 - 26
		if (j == A.rows)
			return 0;
		if (j != i) {
			ex = -ex;
			S_base_modul.swap_rows(A, j, i);
		}

		// Step 27 - 29
		inv_mod(TMP1, S_base_modul.member(A, i, i), mod);

		for (j = i + 1; j < A.rows; j++) {
			mult_mod(TMP2, S_base_modul.member(A, j, i), TMP1, mod);
			for (z = i + 1; z < A.rows; z++) {
				mult_mod(TMP, S_base_modul.member(A, i, z), TMP2, mod);
				sub_mod(TMP, S_base_modul.member(A, j, z), TMP, mod);
				S_base_modul.sto(A, j, z, TMP);
			}
		}
		mult_mod(ret, ret, S_base_modul.member(A, i, i), mod);
	}
	if (ex < 0)
		ret = -ret;
	return ret;
}



//
// Hessenberg form
//

template< class T, class MATRIX_TYPE >
void
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::HBF (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: //value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// Task: HBF_intern(value, r, mod);
	// =  > matrix (value, r, r) in Hessenberg form
	// Version: 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "HBF_intern(T **, lidia_size_t, T)");

	// Step 1, 2
	lidia_size_t i, j, z;
	T TMP, TMP1, TMP2;
	T *tmp;

	// Step 3 - 11
	for (i = A.rows - 1; i >= 1; i--) {
		for (j = i - 1; j >= 0 && A.value[i][j] == 0; j--);
		if (j != -1) {

			// Step 12, 13
			if (j != i - 1) {

				// Step 14 - 18
				// exchange columns i-1 and j
				for (z = 0; z < A.rows; z++) {
					TMP = A.value[z][i - 1];
					A.value[z][i - 1] = A.value[z][j];
					A.value[z][j] = TMP;
				}

				// Step 19 - 24
				// exchange rows i-1 and j
				tmp = A.value[i - 1];
				A.value[i - 1] = A.value[j];
				A.value[j] = tmp;
			}
			tmp = A.value[i];

			// Step 25 - 41
			inv_mod(TMP2, tmp[i - 1], mod);
			for (j = i - 2; j >= 0; j--) {
				mult_mod(TMP1, tmp[j], TMP2, mod);
				for (z = 0; z < A.rows; z++) {
					mult_mod(TMP, A.value[z][i - 1], TMP1, mod);
					sub_mod(A.value[z][j], A.value[z][j], TMP, mod);
				}
				for (z = 0; z < A.rows; z++) {
					mult_mod(TMP, A.value[j][z], TMP1, mod);
					add_mod(A.value[i - 1][z], A.value[i - 1][z], TMP, mod);
				}
			}
		}
	}
}



//
// characteristic polynomial
//

template< class T, class MATRIX_TYPE >
T *
sparse_fp_matrix_kernel< T, MATRIX_TYPE >::charpoly (MATRIX_TYPE &A, const T &mod) const
{
	//
	// Input: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// Task: RES = charpoly_intern(value, r, mod);
	// =  > RES[0], ..., RES[r] are the coefficients of
	//                 the characteristic polynomial of matrix (value, r, r)
	// Version:  1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "charpoly_intern(T **, lidia_size_t, T)");

	lidia_size_t i, j, z;
	T TMP;
	T *tmp;
	lidia_size_t sign;

	// Step 1 - 5
	HBF(A, mod);

	T *K = new T[A.rows]; // size = r
	memory_handler(K, "bigint_matrix", "charpoly_intern - Version T ::"
		       "Error in memory allocation (K)");
	for (i = 0; i < A.rows; i++)
		K[i] = 1;

	// Step 6 - 8
	T **P = new T *[A.rows + 1];
	memory_handler(P, "bigint_matrix", "charpoly_intern - Version T :: "
		       "Error in memory allocation (P)");
	for (i = 0; i < A.rows + 1; i++) {
		tmp = new T[A.rows + 1];
		memory_handler(tmp, "bigint_matrix", "charpoly_intern - Version T :: "
			       "Error in memory allocation (tmp)");
		for (j = 0; j < A.rows + 1; j++)
			tmp[j] = 0;
		P[i] = tmp;
	}
	P[0][0] = 1;

	// Step 9 - 11
	for (z = 1; z <= A.rows; z++) {

		// Step 12 - 16
		for (j = 1; j < z; j++)
			mult_mod(K[j - 1], K[j - 1], A.value[z - 1][z - 2], mod);

		// Step 17 - 23
		P[z][z] = mod - P[z - 1][z - 1];
		for (i = 1; i < z; i++) {
			mult_mod(TMP, A.value[z - 1][z - 1], P[i][z - 1], mod);
			sub_mod(P[i][z], TMP, P[i - 1][z - 1], mod);
		}
		mult_mod(P[0][z], A.value[z - 1][z - 1], P[0][z - 1], mod);

		// Step 24 - 34
		sign = 1;
		for (j = z - 1; j >= 1; j--) {
			sign = -sign;
			for (i = 0; i <= j - 1; i++) {
				mult_mod(TMP, sign, P[i][j - 1], mod);
				mult_mod(TMP, TMP, A.value[j - 1][z - 1], mod);
				mult_mod(TMP, TMP, K[j - 1], mod);
				add_mod(P[i][z], P[i][z], TMP, mod);
			}
		}
	}

	// Step 35 - 40
	T *RES = new T[A.rows + 1];
	memory_handler(RES, "bigint_matrix", "charpoly_intern - Version T :: "
		       "Error in memory allocation (RES)");
	for (i = 0; i <= A.rows; i++) {
		RES[i] = P[i][A.rows];
		delete[] P[i];
	}
	delete[] P;
	delete[] K;
	return RES;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_FP_MATRIX_KERNEL_CC_GUARD_
