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


#ifndef LIDIA_BIGINT_MATRIX_ALGORITHMS_CC_GUARD_
#define LIDIA_BIGINT_MATRIX_ALGORITHMS_CC_GUARD_


#ifndef LIDIA_INFO_H_GUARD_
# include	"LiDIA/info.h"
#endif
#ifndef LIDIA_DENSE_FP_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_fp_matrix_kernel.h"
#endif
#ifndef LIDIA_MODULAR_ARITHMETIC_H_GUARD_
# include	"LiDIA/matrix/modular_arithmetic.h"
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

#define DVALUE LDBL_MATRIX                   // Debug value
#define DMESSAGE "bigint_matrix_algorithms"  // Debug message
#define EMESSAGE matrix_error_msg            // Error message

//
// divide
//

template< class REP1, class REP2, class REP3 >
inline void bigint_matrix_algorithms< REP1, REP2, REP3 >::
divide(matrix< bigint > &RES, const matrix< bigint > &A,
       const bigint &k) const
{
	register lidia_size_t j, i;
	bigint TMP;

	for (j = 0; j < A.rows; j++)
		for (i = 0; i < A.columns; i++) {
			LiDIA::divide(TMP, rep_modul2.member(A, j, i), k);
			rep_modul1.sto(RES, j, i, TMP);
		}
}



template< class REP1, class REP2, class REP3 >
inline void bigint_matrix_algorithms< REP1, REP2, REP3 >::
compwise_divide(matrix< bigint > &RES, const matrix< bigint > &A,
		const matrix< bigint > &B) const
{
	register lidia_size_t j, i;
	bigint TMP;

	for (j = 0; j < RES.rows; j++)
		for (i = 0; i < RES.columns; i++) {
			LiDIA::divide(TMP, rep_modul2.member(A, j, i), rep_modul3.member(B, j, i));
			rep_modul1.sto(RES, j, i, TMP);
		}
}



//
// remainder
//

template< class REP1, class REP2, class REP3 >
void bigint_matrix_algorithms< REP1, REP2, REP3 >::
remainder(matrix< bigint > &RES, const matrix< bigint > & M,
	  const bigint & mod) const
{
	RES.set_representation(M.bitfield.get_representation());

	if (RES.bitfield.get_representation() == matrix_flags::dense_representation) {
		lidia_size_t i, j;
		bigint *REStmp, *Mtmp;

		for (i = 0; i < M.rows; i++) {
			REStmp = RES.value[i];
			Mtmp = M.value[i];
			for (j = 0; j < M.columns; j++)
				LiDIA::best_remainder(REStmp[j], Mtmp[j], mod);
		}
	}
	else {
		lidia_size_t i, j;
		bigint *REStmp;
		lidia_size_t *REStmp2;
		bigint *Mtmp;

		for (i = 0; i < M.rows; i++) {
			Mtmp = M.value[i];

			if (RES.allocated[i] < M.value_counter[i]) {
				delete[] RES.index[i];
				delete[] RES.value[i];

				REStmp = new bigint[M.value_counter[i]];
				REStmp2 = new lidia_size_t[M.value_counter[i]];
				RES.value_counter[i] = M.value_counter[i];
				RES.allocated[i] = M.value_counter[i];
			}
			else {
				REStmp = RES.value[i];
				REStmp2 = RES.index[i];
				RES.value_counter[i] = M.value_counter[i];
			}

			for (j = 0; j < M.value_counter[i]; j++) {
				LiDIA::best_remainder(REStmp[j], Mtmp[j], mod);
				REStmp2[j] = M.index[i][j];
			}
			RES.value[i] = REStmp;
			RES.index[i] = REStmp2;
		}
	}
}



template< class REP1, class REP2, class REP3 >
inline void bigint_matrix_algorithms< REP1, REP2, REP3 >::
trans_remainder(matrix< bigint > &RES, const matrix< bigint > & M,
		const bigint & mod) const
{

	lidia_size_t i, j;
	if (RES.bitfield.get_representation() == matrix_flags::dense_representation) {
		bigint *REStmp;

		for (i = 0; i < M.columns; i++) {
			REStmp = RES.value[i];
			for (j = 0; j < M.rows; j++)
				LiDIA::best_remainder(REStmp[j], M.value[j][i], mod);
		}
	}
	else {
		bigint TMP;
		for (i = 0; i < M.columns; i++) {
			for (j = 0; j < M.rows; j++) {
				LiDIA::best_remainder(TMP, M(j, i), mod);
				RES.sto(i, j, TMP);
			}
		}
	}
}



template< class REP1, class REP2, class REP3 >
void bigint_matrix_algorithms< REP1, REP2, REP3 >::
remainder(matrix< long > &RES, const matrix< bigint > & M, long mod) const
{
	RES.set_representation(M.bitfield.get_representation());

	if (RES.bitfield.get_representation() == matrix_flags::dense_representation) {
		lidia_size_t i, j;
		long *REStmp;
		bigint *Mtmp;

		for (i = 0; i < M.rows; i++) {
			REStmp = RES.value[i];
			Mtmp = M.value[i];
			for (j = 0; j < M.columns; j++)
				LiDIA::best_remainder(REStmp[j], Mtmp[j], mod);
		}
	}
	else {
		lidia_size_t i, j;
		long *REStmp;
		lidia_size_t *REStmp2;
		bigint *Mtmp;

		for (i = 0; i < M.rows; i++) {
			Mtmp = M.value[i];

			if (RES.allocated[i] < M.value_counter[i]) {
				delete[] RES.index[i];
				delete[] RES.value[i];

				REStmp = new long[M.value_counter[i]];
				REStmp2 = new lidia_size_t[M.value_counter[i]];
				RES.value_counter[i] = M.value_counter[i];
				RES.allocated[i] = M.value_counter[i];
			}
			else {
				REStmp = RES.value[i];
				REStmp2 = RES.index[i];
				RES.value_counter[i] = M.value_counter[i];
			}

			for (j = 0; j < M.value_counter[i]; j++) {
				LiDIA::best_remainder(REStmp[j], Mtmp[j], mod);
				REStmp2[j] = M.index[i][j];
			}
			RES.value[i] = REStmp;
			RES.index[i] = REStmp2;
		}
	}
}



template< class REP1, class REP2, class REP3 >
inline void bigint_matrix_algorithms< REP1, REP2, REP3 >::
trans_remainder(matrix< long > &RES, const matrix< bigint > & M,
		long mod) const
{
	lidia_size_t i, j;
	if (RES.bitfield.get_representation() == matrix_flags::dense_representation) {
		long *REStmp;

		for (i = 0; i < M.columns; i++) {
			REStmp = RES.value[i];
			for (j = 0; j < M.rows; j++)
				LiDIA::best_remainder(REStmp[j], M.value[j][i], mod);
		}
	}
	else {
		long TMP;
		for (i = 0; i < M.columns; i++) {
			for (j = 0; j < M.rows; j++) {
				LiDIA::best_remainder(TMP, M(j, i), mod);
				RES.sto(i, j, TMP);
			}
		}
	}
}



//
// Kernel
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
kernel1(matrix< bigint > &RES, const matrix< bigint > & A, const bigint &H) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"kernel1(const matrix< bigint > &)", DVALUE + 8);

	const modular_arithmetic< DRMK < bigint >,
		dense_fp_matrix_kernel< long, MR < long > >,
		dense_fp_matrix_kernel< bigint, MR < bigint > > > Dm_bigint_modul;

	bigint *ZBAtmp, *Atmp;
	register long i, j;
	lidia_size_t c = A.columns;

	// Step 1
	lidia_size_t *linuz = Dm_bigint_modul.lininr1(A, H);
	lidia_size_t r = linuz[0];
	if (r == c) {
		matrix< bigint > RES1(c, 1);
		RES.assign(RES1);
		delete[] linuz;
		return;
	}

	// Step 2
	matrix< bigint > ZBA(r, c);
	for (i = 1; i <= r; i++) {
		ZBAtmp = ZBA.value[i - 1];
		Atmp = A.value[linuz[r - i + 1]];
		for (j = 0; j < c; j++)
			ZBAtmp[j].assign(Atmp[j]);
	}
	delete[] linuz;

	// Step 3
	matrix< bigint > TRANS(c, c);
	ZBA.hnf(TRANS);

	// Step 4
	matrix< bigint > PART2(c, r);
	if (RES.rows != c)
		RES.set_no_of_rows(c);
	if (RES.columns != c - r)
		RES.set_no_of_columns(c - r);

	//dmodul2.split_h(TRANS, RES, PART2);
	dmodul.insert_at(RES, 0, 0, TRANS, 0, 0, RES.rows, RES.columns);
	dmodul.insert_at(PART2, 0, 0, TRANS, 0, TRANS.columns - PART2.columns, PART2.rows, PART2.columns);
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
inline void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
kernel2(matrix< bigint > &RES, const matrix< bigint > & A) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"kernel2(const matrix< bigint > &)", DVALUE + 8);

	register lidia_size_t i;
	matrix< bigint > B = A;
	B.hnf_havas(RES);

	for (i = 0; i < A.columns && B.is_column_zero(i); i++);

	if (i == 0) {
		matrix< bigint > C(A.columns, 1);
		RES.assign(C);
	}
	else
		RES.set_no_of_columns(i);
}



//
// regular InvImage
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
reginvimage1(matrix< bigint > &RES, const matrix< bigint > & A,
	     const matrix< bigint > & B, const bigint &H) const
{
	//
	// Task: C.reginvimage1(A, B);
	// =  > A * C.column(j) = g(j)*B.column(j), j = 0, ..., B.columns
	// =  > g(j) minimal
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "reginvimage1(const matrix< bigint > &, const matrix< bigint > &",
			DVALUE + 8);

	register long i, j;
	bigint TMP, TMP1;

	// Step 1
	bigint DET;
	A.det(DET, H);
	if (DET == 0) {
		precondition_error_handler(DET, "det(A)", "det(A) != 0",
				    "void matrix< bigint >::"
				    "reginvimage1(const matrix< bigint > & A, const matrix< bigint > & B)",
				    DMESSAGE, EMESSAGE[11]);
		return;
	}
	matrix< bigint > ADJ;

	const modular_arithmetic< DRMK < bigint >,
		dense_fp_matrix_kernel< long, MR < long > >,
		dense_fp_matrix_kernel< bigint, MR < bigint > > > Dm_bigint_modul;

	Dm_bigint_modul.adj1(ADJ, A, H, DET);

  // Step 2
	matrix< bigint > PROD = ADJ * B;
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
reginvimage2(matrix< bigint > &RES, const matrix< bigint > & A,
	     const matrix< bigint > & B, const bigint &H) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"reginvimage2(const matrix< bigint > &, const matrix< bigint > &", DVALUE + 8);

	register lidia_size_t i, j, len, oldlen;
	bigint TMP, TMP1;

	// Step 1
	bigint DET;
	A.det(DET, H);
	if (DET == 0) {
		precondition_error_handler(DET, "det(A)", "det(A) != 0",
				    "void matrix< bigint >::"
				    "reginvimage2(const matrix< bigint > & A, const matrix< bigint > & B)",
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

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
image1(matrix< bigint > &RES, const matrix< bigint > & A, const bigint &H) const
{
	bigint *ZBAtmp, *Atmp;
	register long i, j;

	// Step 1
	const modular_arithmetic< DRMK < bigint >,
		dense_fp_matrix_kernel< long, MR < long > >,
		dense_fp_matrix_kernel< bigint, MR < bigint > > > Dm_bigint_modul;

	lidia_size_t *v = Dm_bigint_modul.lininr1(A, H);
	lidia_size_t RANG = v[0];

	// Step 2, 3
	matrix< bigint > ZBA(RANG, A.columns);
	for (i = 1; i <= RANG; i++) {
		ZBAtmp = ZBA.value[i - 1];
		Atmp = A.value[v[RANG - i + 1]];
		for (j = 0; j < A.columns; j++)
			ZBAtmp[j].assign(Atmp[j]);
	}
	delete[] v;

	// Step 4
	matrix< bigint > TRANS(A.columns, A.columns);
	ZBA.hnf(TRANS);

	// Step 5
	if (RES.rows != A.rows)
		RES.set_no_of_rows(A.rows);
	if (RES.columns != RANG)
		RES.set_no_of_columns(RANG);

	if (A.columns == RANG)
		LiDIA::multiply(RES, A, TRANS);
	else {
		matrix< bigint > M(A.rows, A.columns);
		LiDIA::multiply(M, A, TRANS);
		matrix< bigint > PART1(RANG, A.columns - RANG);
		//dmodul2.split_h(M, PART1, RES);
		dmodul.insert_at(PART1, 0, 0, M, 0, 0, PART1.rows, PART1.columns);
		dmodul.insert_at(RES, 0, 0, M, 0, M.columns - RES.columns, RES.rows, RES.columns);
	}
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
inline void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
image2(matrix< bigint > &RES, const matrix< bigint > & A) const
{
	register lidia_size_t i;
	RES.assign(A);
	RES.hnf_havas();

	for (i = 0; i < RES.columns && RES.is_column_zero(i); i++);

	if (i != 0)
		if (i == RES.columns)
			RES.set_no_of_columns(1);
		else {
			matrix< bigint > M(RES);
			matrix< bigint > PART1(RES.rows, i);
			RES.set_no_of_columns(RES.columns-i);
			//dmodul2.split_h(M, PART1, RES);
			dmodul.insert_at(PART1, 0, 0, M, 0, 0, PART1.rows, PART1.columns);
			dmodul.insert_at(RES, 0, 0, M, 0, M.columns - RES.columns, RES.rows, RES.columns);
		}
}



//
// InvImage
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
invimage(matrix< bigint > &RES, const matrix< bigint > & B, const bigint * b) const
{
	if (b == NULL)
		precondition_error_handler(PRT, "b", "b != NULL",
				    "void matrix< bigint >::"
				    "invimage(const matrix< bigint > & B, const bigint * b)",
				    DMESSAGE, EMESSAGE[1]);

	register long i;
	bigint *tmp;

	// Step 1
	matrix< bigint > A = B;
	A.set_no_of_columns(B.columns + 1);
	for (i = 0; i < B.rows; i++)
		A.value[i][B.columns].assign(-b[i]);

	kernel1(RES, A, LiDIA::hadamard(A));

	// Step 2
	if (RES.is_column_zero(0) || RES.is_row_zero(B.columns)) {
		matrix< bigint > C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	tmp = new bigint[RES.columns];
	RES.get_row(tmp, B.columns);
	bigint *g = LiDIA::mgcd2(tmp, RES.columns);
	delete[] tmp;
	if (g[0] > 1) {
		matrix< bigint > C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	// Step 3, 4 
	bigint *x = (RES) * &(g[1]);
	delete[] g;

	// Step 5
	kernel1(RES, B, LiDIA::hadamard(B));
	RES.set_no_of_columns(RES.columns + 1);
	for (i = 0; i < RES.rows; i++)
		RES.value[i][RES.columns-1].assign(x[i]);
	delete[] x;
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
invimage(matrix< bigint > &RES, const matrix< bigint > & B, const math_vector< bigint > &b) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"invimage(const matrix< bigint > &, const math_vector< bigint > &)", DVALUE + 8);

	register long i;
	// Step 1
	matrix< bigint > A = B;
	A.set_no_of_columns(B.columns + 1);

	if (b.size() != B.rows)
		precondition_error_handler(b.size(), "b.size", "b.size == B.rows",
				    B.rows, "B.rows", "b.size == B.rows",
				    "void matrix< bigint >::"
				    "invimage(const matrix< bigint > & B, const math_vector< bigint > &b)",
				    DMESSAGE, EMESSAGE[1]);

	bigint *tmp = b.get_data_address();
	for (i = 0; i < B.rows; i++)
		A.value[i][B.columns].assign(-tmp[i]);
	kernel1(RES, A, LiDIA::hadamard(A));

	// Step 2
	if (RES.is_column_zero(0) || RES.is_row_zero(B.columns)) {
		matrix< bigint > C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	tmp = new bigint[RES.columns];
	RES.get_row(tmp, B.columns);
	bigint *g = LiDIA::mgcd2(tmp, RES.columns);

	delete[] tmp;
	if (g[0] > 1) {
		matrix< bigint > C(B.rows, 1);
		RES.set_no_of_columns(1);
		RES.assign(C);
		return;
	}

	// Step 3, 4 
	bigint *x = (RES) * &(g[1]);
	delete[] g;

	// Step 5
	kernel1(RES, B, LiDIA::hadamard(B));
	RES.set_no_of_columns(RES.columns + 1);
	for (i = 0; i < RES.rows; i++)
		RES.value[i][RES.columns-1].assign(x[i]);
	delete[] x;
}



//
// Smith normal form
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_hartley(matrix< bigint > &RES) const
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_hartley(matrix< bigint > &RES, matrix< bigint > & T1, matrix< bigint > & T2) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_hartley(matrix< bigint > &, matrix< bigint > &)", DVALUE + 8);

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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_simple(matrix< bigint > &RES) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_simple()", DVALUE + 8);
	bigint PIVOT, TMP1, TMP2;
	bigint *tmp = NULL, *deltmp;

	matrix< bigint > TR1(RES.rows, RES.rows);
	matrix< bigint > TR2(RES.columns, RES.columns);

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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_simple(matrix< bigint > &RES, matrix< bigint > & T1, matrix< bigint > & T2) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_simple(matrix< bigint > &, matrix< bigint > &)", DVALUE + 8);
	bigint PIVOT, TMP1, TMP2;
	bigint *tmp = NULL, *deltmp;

	matrix< bigint > TR1 = T1;
	matrix< bigint > TR2 = T2;
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_havas(matrix< bigint > &RES) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_havas()", DVALUE + 8);

	register lidia_size_t i, j, z, index;
	bigint PIVOT;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, SW, TEILBARKEIT;
	bigint TMP1, TMP2;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {
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
				do
				{
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
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
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_havas(matrix< bigint > &RES, matrix< bigint > & T1, matrix< bigint > & T2) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_havas(matrix< bigint > &, matrix< bigint > &)", DVALUE + 8);

	register lidia_size_t i, j, z, index;
	bigint PIVOT;
	bigint *tmp = NULL;

	register lidia_size_t startr, startc, xpivot, ypivot, SW, TEILBARKEIT;
	bigint TMP1, TMP2;

	for (startc = 0, startr = 0; startr < RES.rows && startc < RES.columns; startr++, startc++) {
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
				do
				{
					SW = 0;
					for (i = 0; i < RES.columns; i++)
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
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
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_mult(matrix< bigint > &RES, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_mult(matrix< bigint > &, matrix< bigint > &)", DVALUE + 8);
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
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
					LiDIA::multiply(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if (((abs_compare(PIVOT, RES.value[i][j]) > 0) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
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
				do
				{
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
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_mult(matrix< bigint > &RES, matrix< bigint > & T1, matrix< bigint > & T2, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_mult(matrix< bigint > &, matrix< bigint > &)", DVALUE + 8);

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
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
					LiDIA::multiply(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if (((abs_compare(PIVOT, RES.value[i][j]) > 0) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
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
				do
				{
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
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_add(matrix< bigint > &RES, long art) const
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
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
					LiDIA::add(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if (((abs_compare(PIVOT, RES.value[i][j]) > 0) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
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
				do
				{
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
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_add(matrix< bigint > &RES, matrix< bigint > & T1, matrix< bigint > & T2, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_add(matrix< bigint > &, matrix< bigint > &, long)", DVALUE + 8);

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
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
					LiDIA::add(TMP1, ROW, COLUMN);
					if (TMP1 < NORM) {
						NORM.assign(TMP1);
						PIVOT.assign(abs(RES.value[i][j]));
						xpivot = i;
						ypivot = j;
					}
				}

				if (((abs_compare(PIVOT, RES.value[i][j]) > 0) && !RES.value[i][j].is_zero()) || PIVOT.is_zero()) {
					PIVOT.assign(abs(RES.value[i][j]));
					RES.row_norm(ROW, i, art);
					RES.column_norm(COLUMN, j, art);
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
				do
				{
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
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_new(matrix< bigint > &RES, long art) const
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
			RES.row_norm(RO[i], i, art);
		for (i = 0; i < RES.columns; i++)
			RES.column_norm(CO[i], i, art);

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
					if (abs_compare(PIVOT, RES.value[i][j]) > 0) {
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
				do
				{
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
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc] , RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snf_new(matrix< bigint > &RES, matrix< bigint > & T1, matrix< bigint > & T2, long art) const
{
	debug_handler_l(DMESSAGE, "in member - function "
			"snf_new(matrix< bigint > &, matrix< bigint > &, long)", DVALUE + 8);
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
			RES.row_norm(RO[i], i, art);
		for (i = 0; i < RES.columns; i++)
			RES.column_norm(CO[i], i, art);

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
					if (abs_compare(PIVOT, RES.value[i][j]) > 0) {
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
				do
				{
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
						if ((abs_compare(tmp[index], tmp[i]) > 0) && !tmp[i].is_zero())
							index = i;
				}
				while (SW == 1);

				for (i = 0; RES.value[startr][i].is_zero(); i++);
				RES.swap_columns(startc, i);
				T2.swap_columns(startc, i);

				// mgcd2(column(startc), rows, TR1);
				index = startr;
				do
				{
					SW = 0;
					for (i = 0; i < RES.rows; i++)
						if ((abs_compare(RES.value[index][startc], RES.value[i][startc]) > 0) && !RES.value[i][startc].is_zero())
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snfmod_dkt(matrix< bigint > &RES, const bigint &mod) const
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
		Atmp[j].assign(gcd(Atmp[j], D));
		div_rem(D, TMP, D, Atmp[j]);
	}
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
snfmod_cohen(matrix< bigint > &RES, const bigint & mod) const
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

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
inline void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
gauss(matrix< bigint > &RES) const
{
	debug_handler(DMESSAGE, "in member - function gauss()");

	matrix< bigint > TR(RES.columns, RES.columns);
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
basis_completion(matrix< bigint > &M, bigint *x, lidia_size_t n) const
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
		M.value[i][i+1].assign_one();
		M.value[i][0].assign(x[i]);
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
			M.value[j][i].assign(TMP * a + M.value[j][i+1] * b);
			M.value[j][i+1].assign(TMP * c + M.value[j][i+1] * d);
		}
		y[i+1] = g;
		y[i] = 0;

		gold = g;
	}
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
cond_matrix(matrix< bigint > &M, bigint *v, lidia_size_t n) const
{
	if (M.rows != n)
		M.set_no_of_rows(n);
	if (M.columns != n)
		M.set_no_of_columns(n);

	M.diag(1, 0);

	lidia_size_t i;

	bigint *b = new bigint[n];
	for (i = 0; i < n; i++)
		b[i] = v[i];

	bigint g, N = b[0], a = b[1], qa, ra, qb, rb, t = 0;

	for (i = 2; i < n && g != 1; i++) {
		if (b[i] != 0) {
			g = gcd(b[i], a);

			div_rem(qa, ra, a/g, N);
			div_rem(qb, rb, b[i]/g, N);

			for (t = 0; gcd((ra + t*rb), N) != 1; t++);

			a = a + t * b[i]; // NEW

			M.sto(1, i, t);
		}
	}
}



//
// mgcd2
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
bigint *modular_bigint_matrix_algorithms< REP, SINGLE_MODUL, MULTI_MODUL >::
mgcd2(matrix< bigint > &RES, const bigint * aconst, lidia_size_t n) const
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

	do
	{
		MIN.assign(a[index]);

		// Pivot search: MINIMUM
		for (i = bound; i < n; i++)
			if ((abs_compare(MIN, a[i]) > 0) && !a[i].is_zero()) {
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



#endif	// LIDIA_BIGINT_MATRIX_ALGORITHMS_CC_GUARD_
