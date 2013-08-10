// -*- C++ -*-
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
//	$Id: modular_functions.inl,v 2.1 2001/05/15 12:07:16 hamdy Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef __LIDIA_MODULAR_FUNCTIONS_INL
#define __LIDIA_MODULAR_FUNCTIONS_INL


#include	<LiDIA/modular_operations.inl>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// column step form
//

inline int
STF_intern(bigint ** value, lidia_size_t r, lidia_size_t c, const bigint & mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows
	//              c = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: ex = STF_intern(value, r, c, mod);
	// =  > matrix (value, r, c) in column step form
	// =  > ex = (-1)^Number_of_column_exchanges
	// VERSION: bigint 1.9
	//

	debug_handler_l("bigint_matrix", "in inline - function "
			"STF_intern(bigint **, lidia_size_t, lidia_size_t, const bigint &)", LDBL_MATRIX);

	register lidia_size_t i, j, z, j0 = c - 1;
	bigint TMP, TMP1, TMP2;
	bigint *tmp;

	// Step 1 - 4
	register int exchange = 1;

	// Step 5 - 8
	for (i = r - 1; i >= 0; i--) {
		tmp = value[i];

		// Step 9 - 13
		for (j = j0; j >= 0 && tmp[j].is_zero(); j--);

		// Step 14 - 26
		if (j != -1) {
			if (j != j0) {
				exchange = -exchange;

				// exchange column j0 with column j
				for (z = 0; z <= i; z++)
					swap(value[z][j0], value[z][j]);
			}
			inv_mod(TMP, tmp[j0], mod);

			// Step 27 - 29
			for (j = j0 - 1; j >= 0; j--) {
				if (tmp[j] != 0) {
					// Step 30 - 40
					mult_mod(TMP1, tmp[j], TMP, mod);

					for (z = 0; z <= i; z++) {
						mult_mod(TMP2, value[z][j0], TMP1, mod);
						sub_mod(value[z][j], value[z][j], TMP2, mod);
					}
				}
			}

			// Step 41 - 48
			j0--;
		}
	}
	return exchange;
}



inline int
STF_intern(long ** value, lidia_size_t r, lidia_size_t c, long mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows
	//              c = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: ex = STF_intern(value, r, c, mod);
	// =  > matrix (value, r, c) in column step form
	// =  > ex = (-1)^Number_of_column_exchanges
	// VERSION: long 1.9
	//

	debug_handler_l("bigint_matrix", "in member - function "
			"STF_intern(long **, lidia_size_t, lidia_size_t, long)", LDBL_MATRIX);

	register lidia_size_t i, j, z, j0 = c - 1;
	long TMP, TMP1, TMP2;
	long *tmp;

	// Step 1 - 4
	register int exchange = 1;

	// Step 5 - 8
	for (i = r - 1; i >= 0; i--) {
		tmp = value[i];

		// Step 9 - 13
		for (j = j0; j >= 0 && tmp[j] == 0; j--);

		// Step 14 - 26
		if (j != -1) {
			if (j != j0) {
				exchange = -exchange;
				for (z = 0; z <= i; z++) {
					tmp = value[z];

					// exchange column j0 with column j
					TMP = tmp[j0];
					tmp[j0] = tmp[j];
					tmp[j] = TMP;
				}
			}
			inv_mod(TMP, tmp[j0], mod);

			// Step 27 - 29
			for (j = j0 - 1; j >= 0; j--) {
				if (tmp[j] != 0) {
					// Step 30 - 40
					mult_mod(TMP1, tmp[j], TMP, mod);

					for (z = 0; z <= i; z++) {
						tmp = value[z];
						mult_mod(TMP2, tmp[j0], TMP1, mod);
						sub_mod(tmp[j], tmp[j], TMP2, mod);
					}
				}
			}

			// Step 41 - 48
			j0--;
		}
	}
	return exchange;
}



//
// rank
//

inline lidia_size_t
rank_intern(bigint ** Avalue, lidia_size_t r, lidia_size_t c, const bigint & Gmod)
{
	//
	// INPUT: **Avalue = values of matrix
	//               r = number of rows
	//               c = number of columns
	//            Gmod = modulus for Fp - class
	// DESCRIPTION: rank_intern(Avalue, r, c, Gmod) = rank of matrix (Avalue, r, c)
	// VERSION: 1.9
	//

	debug_handler_l("bigint_matrix", "in inline - function "
			"rank_intern(bigint **, lidia_size_t, lidia_size_t, const bigint &)", LDBL_MATRIX);

	register lidia_size_t i, j, No = 0;
	const bigint *Atmp;

	if (Gmod.bit_length() > bigint::bits_per_digit()) {
		// bigint part
		bigint **value = new bigint *[r];
		memory_handler(value, "bigint_matrix", "rank_intern - bigint part ::"
			       "Error in memory allocation (value)");
		bigint *tmp;

		for (i = 0; i < r; i++) {
			tmp = new bigint[c];
			Atmp = Avalue[i];
			memory_handler(tmp, "bigint_matrix", "rank_intern - bigint part :: "
				       "Error in memory allocation (tmp)");
			for (j = 0; j < c; j++)
				best_remainder(tmp[j], Atmp[j], Gmod);
			value[i] = tmp;
		}

		// Step 1, 2
		STF_intern(value, r, c, Gmod);

		// Step 3 - 12
		for (j = 0; j < c; j++) {
			for (i = r - 1; i >= 0 && value[i][j].is_zero(); i--);
			if (i == -1)
				No++;
		}

		for (j = 0; j < r; j++)
			delete[] value[j];
		delete[] value;
	}
	else {
		// long part
		long mod;
		Gmod.longify(mod);

		long **value = new long *[r];
		memory_handler(value, "bigint_matrix", "rank_intern - long part :: "
			       "Error in memory allocation (value)");
		long *tmp;

		for (i = 0; i < r; i++) {
			tmp = new long[c];
			Atmp = Avalue[i];
			memory_handler(tmp, "bigint_matrix", "rank_intern - long part :: "
				       "Error in memory allocation (tmp)");
			for (j = 0; j < c; j++)
				best_remainder(tmp[j], Atmp[j], mod);
			value[i] = tmp;
		}

		// Step 1, 2
		STF_intern(value, r, c, mod);

		// Step 3 - 12
		for (j = 0; j < c; j++) {
			for (i = r - 1; i >= 0 && value[i][j] == 0; i--);
			if (i == -1)
				No++;
		}
		for (j = 0; j < r; j++)
			delete[] value[j];
		delete[] value;
	}

	// Step 13 - 24
	return (c - No);
}



//
// rank and linearly independent rows or columns
//

inline lidia_size_t *
lininr_intern(bigint ** Avalue, lidia_size_t r, lidia_size_t c, const bigint & Gmod)
{
	//
	// INPUT: **Avalue = values of matrix
	//               r = number of rows
	//               c = number of columns
	//            Gmod = modulus for Fp - class
	// DESCRIPTION: IND = lininr_intern(Avalue, r, c, Gmod);
	// =  > IND[0] = Rank of matrix (Avalue, r, c)
	// =  > IND[1], ..., IND[IND[0]], such that
	//                 raw(IND[1]), ..., raw(IND[IND[0]])
	//                 of matrix (Avalue, r, c) are linearly independent.
	// VERSION: 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "lininr_intern(bigint **, lidia_size_t, lidia_size_t, const bigint &)");

	register lidia_size_t i, j;
	const bigint *Atmp;

	long *l = new long[c + 1];
	memory_handler(l, "bigint_matrix", "lininr_intern :: "
		       "Error in memory allocation (l)");

	if (Gmod.bit_length() > bigint::bits_per_digit()) {
		// bigint part
		bigint **value = new bigint *[r];
		memory_handler(value, "bigint_matrix", "lininr_intern - bigint part ::"
			       "Error in memory allocation (value)");
		bigint *tmp;

		for (i = 0; i < r; i++) {
			tmp = new bigint[c];
			Atmp = Avalue[i];
			memory_handler(tmp, "bigint_matrix", "lininr_intern - bigint part :: "
				       "Error in memory allocation (tmp)");
			for (j = 0; j < c; j++)
				best_remainder(tmp[j], Atmp[j], Gmod);
			value[i] = tmp;
		}

		// Step 1, 2
		STF_intern(value, r, c, Gmod);

		// Step 3 - 12
		for (j = 0; j < c; j++) {
			for (i = r - 1; i >= 0 && value[i][j].is_zero(); i--);
			l[j] = i;
		}

		for (j = 0; j < r; j++)
			delete[] value[j];
		delete[] value;
	}
	else {
		// long part
		long mod;
		Gmod.longify(mod);

		long **value = new long *[r];
		memory_handler(value, "bigint_matrix", "lininr_intern - long part :: "
			       "Error in memory allocation (value)");
		long *tmp;

		for (i = 0; i < r; i++) {
			tmp = new long[c];
			Atmp = Avalue[i];
			memory_handler(tmp, "bigint_matrix", "lininr_intern - long part :: "
				       "Error in memory allocation (tmp)");
			for (j = 0; j < c; j++)
				best_remainder(tmp[j], Atmp[j], mod);
			value[i] = tmp;
		}

		// Step 1, 2
		STF_intern(value, r, c, mod);

		// Step 3 - 12
		for (j = 0; j < c; j++) {
			for (i = r - 1; i >= 0 && value[i][j] == 0; i--);
			l[j] = i;
		}
		for (j = 0; j < r; j++)
			delete[] value[j];
		delete[] value;
	}

	// Step 13 - 24
	for (i = 0; i < c && l[i] == -1; i++); // i = number of zero-columns

	lidia_size_t TMP = c - i; // rank
	lidia_size_t *IND = new lidia_size_t[TMP + 1];
	memory_handler(IND, "bigint_matrix", "lininr_intern - "
		       "Error in memory allocation (IND)");
	IND[0] = TMP; // rank
	for (j = 0; j < TMP; j++)
		IND[TMP - j] = l[j + i];
	delete[] l;
	return IND;
}



inline lidia_size_t *
lininc_intern(bigint ** Avalue, lidia_size_t r, lidia_size_t c, const bigint & Gmod)
{
	//
	// INPUT: **Avalue = values of matrix
	//               r = number of rows
	//               c = number of columns
	//            Gmod = modulus for Fp - class
	// DESCRIPTION: IND = lininc_intern(Avalue, r, c, Gmod);
	// =  > IND[0] = Rank of matrix (Avalue, r, c)
	// =  > IND[1], ..., IND[IND[0]], such that
	//                 column(IND[1]), ..., column(IND[IND[0]])
	//                 of matrix (Avalue, r, c) are linearly independent.
	// VERSION: 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "lininc_intern(bigint **, lidia_size_t, lidia_size_t, const bigint &)");

	register lidia_size_t i, j;
	lidia_size_t *l = new lidia_size_t[r + 1];
	memory_handler(l, "bigint_matrix", "lininc_intern :: "
		       "Error in memory allocation (l)");
	if (Gmod.bit_length() > bigint::bits_per_digit()) {
		// bigint part
		bigint **value = new bigint *[c];
		memory_handler(value, "bigint_matrix", "lininc_intern - bigint part ::"
			       "Error in memory allocation (value)");
		bigint *tmp;

		for (i = 0; i < c; i++) {
			tmp = new bigint[r];
			memory_handler(tmp, "bigint_matrix", "lininc_intern - bigint part :: "
				       "Error in memory allocation (tmp)");
			for (j = 0; j < r; j++)
				best_remainder(tmp[j], Avalue[j][i], Gmod);
			value[i] = tmp;
		}

		// Step 1, 2
		STF_intern(value, c, r, Gmod);

		// Step 3 - 12
		for (j = 0; j < r; j++) {
			for (i = c - 1; i >= 0 && value[i][j].is_zero(); i--);
			l[j] = i;
		}

		for (j = 0; j < c; j++)
			delete[] value[j];
		delete[] value;
	}
	else {
		// long part
		long mod;
		Gmod.longify(mod);

		long **value = new long *[c];
		memory_handler(value, "bigint_matrix", "lininc_intern - long part :: "
			       "Error in memory allocation (value)");
		long *tmp;

		for (i = 0; i < c; i++) {
			tmp = new long[r];
			memory_handler(tmp, "bigint_matrix", "lininc_intern - long part :: "
				       "Error in memory allocation (tmp)");
			for (j = 0; j < r; j++)
				best_remainder(tmp[j], Avalue[j][i], mod);
			value[i] = tmp;
		}

		// Step 1, 2
		STF_intern(value, c, r, mod);

		// Step 3 - 12
		for (j = 0; j < r; j++) {
			for (i = c - 1; i >= 0 && value[i][j] == 0; i--);
			l[j] = i;
		}
		for (j = 0; j < c; j++)
			delete[] value[j];
		delete[] value;
	}

	// Step 13 - 24
	for (i = 0; i < r && l[i] == -1; i++); // i = number of zero-columns

	lidia_size_t TMP = r - i;
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

inline void
adj_intern(bigint ** value, lidia_size_t r, const bigint & mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: adj_intern(value, r, mod);
	// =  > adjoint matrix (value, r, r)
	// VERSION: bigint 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "adj_intern(bigint **, lidia_size_t, const bigint &)");

	register lidia_size_t i, j, z;
	bigint TMP, TMP1, TMP2;
	bigint *tmp, *tmp1, *Btmp, *Btmp1;
	lidia_size_t exchange = 1;
	bigint DET = 1;

	// Step 1, 2
	bigint **Bvalue = new bigint *[r];
	memory_handler(Bvalue, "bigint_matrix", "adj_intern - bigint part :: "
		       "Error in memory allocation (Bvalue)");
	for (i = 0; i < r; i++) {
		Btmp = new bigint[r];
		tmp = value[i];
		memory_handler(Btmp, "bigint_matrix", "adj_intern - bigint part :: "
			       "Error in memory allocation (Btmp)");
		for (j = 0; j < r; j++) {
			Btmp[j].assign(tmp[j]);
			if (i == j)
				tmp[j].assign_one();
			else
				tmp[j].assign_zero();
		}
		Bvalue[i] = Btmp;
	}

	// Step 3 - 5
	for (i = r - 1; i >= 0; i--) {
		Btmp = Bvalue[i];

		// Step 6 - 11
		for (j = i; j >= 0 && Btmp[j].is_zero(); j--);

		// Step 12 - 19
		if (j != i) {
			exchange = -exchange;
			for (z = 0; z < r; z++) {
				Btmp1 = Bvalue[z];
				tmp1 = value[z];

				// A.swap_columns(i, j);
				TMP.assign(Btmp1[j]);
				Btmp1[j].assign(Btmp1[i]);
				Btmp1[i].assign(TMP);

				// B.swap_columns(i, j);
				TMP.assign(tmp1[i]);
				tmp1[i].assign(tmp1[j]);
				tmp1[j].assign(TMP);
			}
		}
		inv_mod(TMP1, Btmp[i], mod);

		// Step 20 - 32
		for (j = 0; j < r; j++) {
			if (j != i) {
				mult_mod(TMP2, Btmp[j], TMP1, mod);
				for (z = 0; z < i; z++) {
					Btmp1 = Bvalue[z];

					mult_mod(TMP, Btmp1[i], TMP2, mod);
					sub_mod(Btmp1[j], Btmp1[j], TMP, mod);
				}
				for (z = 0; z < r; z++) {
					tmp1 = value[z];

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
		for (z = 0; z < r; z++) {
			tmp1 = value[z];
			mult_mod(tmp1[i], tmp1[i], TMP1, mod);
		}
	}

	// Step 33 - 43
	if (exchange < 0)
		DET.negate();
	for (j = 0; j < r; j++) {
		tmp = value[j];
		for (i = 0; i < r; i++)
			mult_mod(tmp[i], tmp[i], DET, mod);
	}
	for (j = 0; j < r; j++)
		delete[] Bvalue[j];
	delete[] Bvalue;
}



inline void
adj_intern(long **value, lidia_size_t r, long mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: adj_intern(value, r, mod);
	// =  > adjoint matrix (value, r, r)
	// VERSION: long 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "adj_intern(long **, lidia_size_t, long)");

	register lidia_size_t i, j, z;
	long TMP, TMP1, TMP2;
	long *tmp, *tmp1, *Btmp, *Btmp1;
	lidia_size_t exchange = 1;
	long DET = 1;

	// Step 1, 2
	long **Bvalue = new long *[r];
	memory_handler(Bvalue, "bigint_matrix", "adj_intern - long part :: "
		       "Error in memory allocation (Bvalue)");
	for (i = 0; i < r; i++) {
		Btmp = new long[r];
		tmp = value[i];
		memory_handler(Btmp, "bigint_matrix", "adj_intern - long part :: "
			       "Error in memory allocation (Btmp)");
		for (j = 0; j < r; j++) {
			Btmp[j] = tmp[j];
			if (i == j)
				tmp[j] = 1;
			else
				tmp[j] = 0;
		}
		Bvalue[i] = Btmp;
	}

	// Step 3 - 5
	for (i = r - 1; i >= 0; i--) {
		Btmp = Bvalue[i];

		// Step 6 - 11
		for (j = i; j >= 0 && Btmp[j] == 0; j--);

		// Step 12 - 19
		if (j != i) {
			exchange = -exchange;
			for (z = 0; z < r; z++) {
				Btmp1 = Bvalue[z];
				tmp1 = value[z];

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
		for (j = 0; j < r; j++) {
			if (j != i) {
				mult_mod(TMP2, Btmp[j], TMP1, mod);
				for (z = 0; z < i; z++) {
					Btmp1 = Bvalue[z];

					mult_mod(TMP, Btmp1[i], TMP2, mod);
					sub_mod(Btmp1[j], Btmp1[j], TMP, mod);
				}

				for (z = 0; z < r; z++) {
					tmp1 = value[z];

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
		for (z = 0; z < r; z++) {
			tmp1 = value[z];
			mult_mod(tmp1[i], tmp1[i], TMP1, mod);
		}
	}

	// Step 33 - 43
	if (exchange < 0)
		DET = -DET;
	for (j = 0; j < r; j++) {
		tmp = value[j];
		for (i = 0; i < r; i++)
			mult_mod(tmp[i], tmp[i], DET, mod);
	}

	for (j = 0; j < r; j++)
		delete[] Bvalue[j];
	delete[] Bvalue;
}



//
// determinant
//

inline void
det_intern(bigint & ret, bigint ** value, lidia_size_t r, const bigint & mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: det_intern(ret, value, r, mod);
	// =  > ret = determinant of matrix (value, r, r)
	// VERSION: bigint 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "det_intern(bigint &, bigint **, lidia_size_t, const bigint &)");

	register lidia_size_t i, j, z;
	bigint TMP, TMP1, TMP2;
	bigint *tmp, *tmp1;

	// Step 1 - 4
	lidia_size_t ex = 1;
	ret.assign_one();

	// Step 5 - 8
	for (i = 0; i < r; i++) {

		// Step 9 - 13
		for (j = i; j < r && value[j][i].is_zero(); j++);

		// Step 14 - 26
		if (j == r) {
			ret.assign_zero();
			return;
		}
		if (j != i) {
			ex = -ex;
			tmp1 = value[j];
			value[j] = value[i];
			value[i] = tmp1;
		}
		tmp = value[i];

		// Step 27 - 29
		inv_mod(TMP1, tmp[i], mod);
		for (j = i + 1; j < r; j++) {

			// Step 30 - 40
			tmp1 = value[j];
			mult_mod(TMP2, tmp1[i], TMP1, mod);
			for (z = i + 1; z < r; z++) {
				mult_mod(TMP, tmp[z], TMP2, mod);
				sub_mod(tmp1[z], tmp1[z], TMP, mod);
			}
		}
		mult_mod(ret, ret, tmp[i], mod);
	}
	if (ex < 0)
		ret.negate();
}



inline lidia_size_t
det_intern(long **value, lidia_size_t r, long mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: det_intern(value, r, mod) = determinant of matrix (value, r, r);
	// VERSION: long 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "det_intern(long **, lidia_size_t, long)");

	register lidia_size_t i, j, z;
	long TMP, TMP1, TMP2;
	long *tmp, *tmp1;

	// Step 1 - 4
	lidia_size_t ex = 1;
	long ret = 1;

	// Step 5 - 8
	for (i = 0; i < r; i++) {

		// Step 9 - 13
		for (j = i; j < r && value[j][i] == 0; j++);

		// Step 14 - 26
		if (j == r)
			return 0;
		if (j != i) {
			ex = -ex;
			tmp1 = value[j];
			value[j] = value[i];
			value[i] = tmp1;
		}
		tmp = value[i];

		// Step 27 - 29
		inv_mod(TMP1, tmp[i], mod);

		for (j = i + 1; j < r; j++) {

			// Step 30 - 40
			tmp1 = value[j];
			mult_mod(TMP2, tmp1[i], TMP1, mod);
			for (z = i + 1; z < r; z++) {
				mult_mod(TMP, tmp[z], TMP2, mod);
				sub_mod(tmp1[z], tmp1[z], TMP, mod);
			}
		}
		mult_mod(ret, ret, tmp[i], mod);
	}
	if (ex < 0)
		ret = -ret;
	return ret;
}



//
// Hessenberg form
//

inline void
HBF_intern(bigint ** value, lidia_size_t r, const bigint & mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: HBF_intern(value, r, mod);
	// =  > matrix (value, r, r) in Hessenberg form
	// VERSION: bigint 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "HBF_intern(bigint **, lidia_size_t, const bigint &)");

	// Step 1, 2
	register lidia_size_t i, j, z;
	bigint TMP, TMP1, TMP2;
	bigint *tmp;

	// Step 3 - 11
	for (i = r - 1; i >= 1; i--) {
		for (j = i - 1; j >= 0 && value[i][j].is_zero(); j--);
		if (j != -1) {

			// Step 12, 13
			if (j != i - 1) {

				// Step 14 - 18
				// exchange columns i-1 and j
				for (z = 0; z < r; z++)
					swap(value[z][i-1], value[z][j]);

				// Step 19 - 24
				// exchange rows i-1 and j
				tmp = value[i - 1];
				value[i - 1] = value[j];
				value[j] = tmp;
			}
			tmp = value[i];

			// Step 25 - 41
			inv_mod(TMP2, tmp[i - 1], mod);
			for (j = i - 2; j >= 0; j--) {
				mult_mod(TMP1, tmp[j], TMP2, mod);
				for (z = 0; z < r; z++) {
					mult_mod(TMP, value[z][i - 1], TMP1, mod);
					sub_mod(value[z][j], value[z][j], TMP, mod);
				}
				for (z = 0; z < r; z++) {
					mult_mod(TMP, value[j][z], TMP1, mod);
					add_mod(value[i - 1][z], value[i - 1][z], TMP, mod);
				}
			}
		}
	}
}



inline void
HBF_intern(long **value, lidia_size_t r, long mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: HBF_intern(value, r, mod);
	// =  > matrix (value, r, r) in Hessenberg form
	// VERSION: long 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "HBF_intern(long **, lidia_size_t, long)");

	// Step 1, 2
	lidia_size_t i, j, z;
	long TMP, TMP1, TMP2;
	long *tmp;

	// Step 3 - 11
	for (i = r - 1; i >= 1; i--) {
		for (j = i - 1; j >= 0 && value[i][j] == 0; j--);
		if (j != -1) {

			// Step 12, 13
			if (j != i - 1) {

				// Step 14 - 18
				// exchange columns i-1 and j
				for (z = 0; z < r; z++) {
					TMP = value[z][i - 1];
					value[z][i - 1] = value[z][j];
					value[z][j] = TMP;
				}

				// Step 19 - 24
				// exchange rows i-1 and j
				tmp = value[i - 1];
				value[i - 1] = value[j];
				value[j] = tmp;
			}
			tmp = value[i];

			// Step 25 - 41
			inv_mod(TMP2, tmp[i - 1], mod);
			for (j = i - 2; j >= 0; j--) {
				mult_mod(TMP1, tmp[j], TMP2, mod);
				for (z = 0; z < r; z++) {
					mult_mod(TMP, value[z][i - 1], TMP1, mod);
					sub_mod(value[z][j], value[z][j], TMP, mod);
				}
				for (z = 0; z < r; z++) {
					mult_mod(TMP, value[j][z], TMP1, mod);
					add_mod(value[i - 1][z], value[i - 1][z], TMP, mod);
				}
			}
		}
	}
}



//
// characteristic polynomial
//

inline bigint *
charpoly_intern(bigint ** value, lidia_size_t r, const bigint & mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: RES = charpoly_intern(value, r, mod);
	// =  > RES[0], ..., RES[r] are the coefficients of the
	//                 characteristic polynomial of matrix (value, r, r)
	// VERSION: bigint 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "charpoly_intern(bigint **, lidia_size_t, const bigint &)");

	register lidia_size_t i, j, z;
	bigint TMP;
	lidia_size_t sign;

	// Step 1 - 5
	HBF_intern(value, r, mod);

	bigint *K = new bigint[r]; // size = c
	memory_handler(K, "bigint_matrix", "charpoly_intern - Version bigint :: "
		       "Error in memory allocation (K)");

	for (i = 0; i < r; i++)
		K[i].assign_one();

	// Step 6 - 8
	bigint *tmp;
	bigint **P = new bigint *[r+1];
	memory_handler(P, "bigint_matrix", "charpoly_intern :: "
		       "Error in memory allocation (P)");
	for (i = 0; i < r+1; i++) {
		tmp = new bigint[r+1];
		memory_handler(tmp, "bigint_matrix", "charpoly_intern :: "
			       "Error in memory allocation (tmp)");
		for (j = 0; j < r+1; j++)
			tmp[j].assign_zero();
		P[i] = tmp;
	}

	P[0][0].assign_one();

	// Step 9 - 11
	for (z = 1; z <= r; z++) {

		// Step 12 - 16
		for (j = 1; j <= z - 1; j++)
			mult_mod(K[j - 1], K[j - 1], value[z - 1][z - 2], mod);

		// Step 17 - 23
		subtract(P[z][z], mod, P[z - 1][z - 1]);
		for (i = 1; i <= z - 1; i++) {
			mult_mod(TMP, value[z - 1][z - 1], P[i][z - 1], mod);
			sub_mod(P[i][z], TMP, P[i - 1][z - 1], mod);
		}
		mult_mod(P[0][z], value[z - 1][z - 1], P[0][z - 1], mod);

		// Step 24 - 34
		sign = 1;
		for (j = z - 1; j >= 1; j--) {
			sign = -sign;
			for (i = 0; i <= j - 1; i++) {
				mult_mod(TMP, sign, P[i][j - 1], mod);
				mult_mod(TMP, TMP, value[j - 1][z - 1], mod);
				mult_mod(TMP, TMP, K[j - 1], mod);
				add_mod(P[i][z], P[i][z], TMP, mod);
			}
		}
	}

	// Step 35 - 40
	bigint *RES = new bigint[r + 1];
	memory_handler(RES, "bigint_matrix", "charpoly_intern - Version bigint :: "
		       "Error in memory allocation (RES)");
	for (i = 0; i < r + 1; i++)
		RES[i].assign(P[i][r]);
	delete[] K;

	for (i = 0; i < r+1; i++)
		delete[] P[i];
	delete[] P;

	return RES;
}



inline long *
charpoly_intern(long **value, lidia_size_t r, long mod)
{
	//
	// INPUT: **value = values of matrix
	//              r = number of rows = number of columns
	//            mod = modulus for Fp - class
	// DESCRIPTION: RES = charpoly_intern(value, r, mod);
	// =  > RES[0], ..., RES[r] are the coefficients of
	//                 the characteristic polynomial of matrix (value, r, r)
	// VERSION: long 1.9
	//

	debug_handler("bigint_matrix", "in inline - function "
		      "charpoly_intern(long **, lidia_size_t, long)");
	lidia_size_t i, j, z;
	long TMP;
	long *tmp;
	lidia_size_t sign;

	// Step 1 - 5
	HBF_intern(value, r, mod);
	long *K = new long[r]; // size = r
	memory_handler(K, "bigint_matrix", "charpoly_intern - Version long ::"
		       "Error in memory allocation (K)");
	for (i = 0; i < r; i++)
		K[i] = 1;

	// Step 6 - 8
	long **P = new long *[r + 1];
	memory_handler(P, "bigint_matrix", "charpoly_intern - Version long :: "
		       "Error in memory allocation (P)");
	for (i = 0; i < r + 1; i++) {
		tmp = new long[r + 1];
		memory_handler(tmp, "bigint_matrix", "charpoly_intern - Version long :: "
			       "Error in memory allocation (tmp)");
		for (j = 0; j < r + 1; j++)
			tmp[j] = 0;
		P[i] = tmp;
	}
	P[0][0] = 1;

	// Step 9 - 11
	for (z = 1; z <= r; z++) {

		// Step 12 - 16
		for (j = 1; j < z; j++)
			mult_mod(K[j - 1], K[j - 1], value[z - 1][z - 2], mod);

		// Step 17 - 23
		P[z][z] = mod - P[z - 1][z - 1];
		for (i = 1; i < z; i++) {
			mult_mod(TMP, value[z - 1][z - 1], P[i][z - 1], mod);
			sub_mod(P[i][z], TMP, P[i - 1][z - 1], mod);
		}
		mult_mod(P[0][z], value[z - 1][z - 1], P[0][z - 1], mod);

		// Step 24 - 34
		sign = 1;
		for (j = z - 1; j >= 1; j--) {
			sign = -sign;
			for (i = 0; i <= j - 1; i++) {
				mult_mod(TMP, sign, P[i][j - 1], mod);
				mult_mod(TMP, TMP, value[j - 1][z - 1], mod);
				mult_mod(TMP, TMP, K[j - 1], mod);
				add_mod(P[i][z], P[i][z], TMP, mod);
			}
		}
	}

	// Step 35 - 40
	long *RES = new long[r + 1];
	memory_handler(RES, "bigint_matrix", "charpoly_intern - Version long :: "
		       "Error in memory allocation (RES)");
	for (i = 0; i <= r; i++)
		RES[i] = P[i][r];
	for (i = 0; i <= r; i++)
		delete[] P[i];
	delete[] P;
	delete[] K;
	return RES;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// __LIDIA_MODULAR_FUNCTIONS_INL
