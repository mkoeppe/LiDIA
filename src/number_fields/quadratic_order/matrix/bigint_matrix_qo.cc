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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/timer.h"
#include	"LiDIA/error.h"
#include	"LiDIA/info.h"
#include	"LiDIA/bigint_matrix.h"
#include	"LiDIA/base_power_product.h"
#include	"LiDIA/matrix/dense_bigint_matrix_modules.h"
#include	"LiDIA/matrix/sparse_bigint_matrix_modules.h"
#include	"LiDIA/matrix/hnf_conf.h"
#include	"LiDIA/matrix/hnf_kernel.h"
#include	"LiDIA/matrix/crt_and_prime_handling.h"
#include	"LiDIA/matrix/dense_fp_matrix_kernel.h"
#include	"LiDIA/matrix/sparse_fp_matrix_kernel.h"
#include	"LiDIA/matrix/sparse_fp_matrix_algorithms.h"
#include	"LiDIA/matrix/bigint_matrix_algorithms.h"
#include	"LiDIA/matrix/modular_arithmetic.h"
#include	"LiDIA/bigint_lattice.h"

#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// debug defines / error defines
//

extern const char *PRT;
extern const char *matrix_error_msg[];

#define DVALUE LDBL_MATRIX           // Debug value
#define DMESSAGE "bigint_matrix"     // Debug message
#define EMESSAGE matrix_error_msg    // Error message

//
// debug level
//
//   0 : remainder
//   1 : divide
//   2 : norms and bounds
//   3 : randomize
//   4 : Linear algebra PART 1
//   5 : Linear algebra PART 2
//

#define DRMKex DRMK< bigint >
#define SRMKex SRMK< bigint >


//////////////////////////////////////
// BEGIN: quadratic order hnf stuff //
//////////////////////////////////////

//
// ADDED BY MJJ
//

void matrix< bigint >::
adj(file_adjoint & RES, const bigint &H, const bigint & DET)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj(const matrix< bigint > &)", DVALUE + 4);

	if (bitfield.get_representation() == matrix_flags::sparse_representation)
		std::cout << "Adjoint with files not implemented for sparse matrices!" << std::endl;

	else {
		register lidia_size_t i, z1, z2;

		matrix< bigint > B(rows, columns);
		register bigint *Atmp;
		long *tmp1;

		dense_fp_matrix_kernel< bigint, MR< bigint > > bigint_modul;
		dense_fp_matrix_kernel< long, MR< long > > long_modul;

		const modular_arithmetic< DRMK< bigint >,
			dense_fp_matrix_kernel< long, MR< long > >,
			dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;


		register bigint *PRIM = get_primes(bigint(2*H), abs(DET), true);
		long n;
		PRIM[0].longify(n);

		lidia_qo_info_handler(std::cout << "(adj) " << n << " primes required!" << std::endl;);

		// Step 4
		bigint MODULUS, MOD;
		long Modlong;
		memory_handler(chininput, DMESSAGE, "adj :: "
			       "Error in memory allocation (chininput)");

		for (i = 1; i <= n; i++) {
			debug_handler_c(DMESSAGE, "in function "
					"adj(const matrix< bigint > &)", DVALUE + 4,
					std::cout << "In Iteration " << i << std::endl;);

			MOD.assign(PRIM[i]);

			if (MOD.bit_length() > bigint::bits_per_digit()) {
				lidia_qo_xinfo_handler("adjoint [bigint]", i, n);

				remainder(B, *this, MOD);
				bigint_modul.adj(B, MOD);
			}
			else {
				lidia_qo_xinfo_handler("adjoint [long]", i, n);

				MOD.longify(Modlong);
				matrix< long > Blong(rows, columns);
				remainder(Blong, *this, Modlong);

				long_modul.adj(Blong, Modlong);

				for (z1 = 0; z1 < rows; z1++) {
					tmp1 = Blong.value[z1];
					Atmp = B.value[z1];
					for (z2 = 0; z2 < columns; z2++)
						Atmp[z2].assign(tmp1[z2]);
				}
			}

			RES.add_new_prime(B, MOD);

			if (RES.test_adjoint(*this, DET)) {
				lidia_qo_info_handler(std::cout << "\n(adj) Only " << i << " primes required!\n" << std::flush;);
				break;
			}
		}

		if (qo_special) {
			if (i == n+1)  --i;
			std::cout << " " << i << std::flush;
		}

		delete[] PRIM;
	}
}



void matrix< bigint >::
adj(file_adjoint & RES, const bigint &H, const bigint & DET, int num_mod)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj(const matrix< bigint > &)", DVALUE + 4);

	if (bitfield.get_representation() == matrix_flags::sparse_representation)
		std::cout << "Adjoint with files not implemented for sparse matrices!" << std::endl;

	else {
		register lidia_size_t i, z1, z2;

		matrix< bigint > B(rows, columns);
		register bigint *Atmp;
		long *tmp1;

		dense_fp_matrix_kernel< bigint, MR< bigint > > bigint_modul;
		dense_fp_matrix_kernel< long, MR< long > > long_modul;

		const modular_arithmetic< DRMK< bigint >,
			dense_fp_matrix_kernel< long, MR< long > >,
			dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;


		register bigint *PRIM = get_primes(bigint(2*H), abs(DET), true);
		long n;
		PRIM[0].longify(n);

		lidia_qo_info_handler(std::cout << "(adj) " << num_mod << " primes required!";
				      std::cout << std::endl;);

		// Step 4
		bigint MODULUS, MOD;
		long Modlong;
		memory_handler(chininput, DMESSAGE, "adj :: "
			       "Error in memory allocation (chininput)");

		bool done = false;
		for (i = 1; (i <= num_mod) || !done; i++) {
			debug_handler_c(DMESSAGE, "in function "
					"adj(const matrix< bigint > &)", DVALUE + 4,
					std::cout << "In Iteration " << i << std::endl;);

			MOD.assign(PRIM[i]);

			if (MOD.bit_length() > bigint::bits_per_digit()) {
				lidia_qo_xinfo_handler("adjoint [bigint]", i, n);

				remainder(B, *this, MOD);
				bigint_modul.adj(B, MOD);
			}
			else {
				lidia_qo_xinfo_handler("adjoint [long]", i, n);

				MOD.longify(Modlong);
				matrix< long > Blong(rows, columns);
				remainder(Blong, *this, Modlong);

				long_modul.adj(Blong, Modlong);

				for (z1 = 0; z1 < rows; z1++) {
					tmp1 = Blong.value[z1];
					Atmp = B.value[z1];
					for (z2 = 0; z2 < columns; z2++)
						Atmp[z2].assign(tmp1[z2]);
				}
			}

			RES.add_new_prime(B, MOD);

			if (i >= num_mod)
				done = RES.test_adjoint_diag(*this, DET);

			if ((i >= num_mod) && done) {
				lidia_qo_info_handler(std::cout << "\n(adj) Only " << i << " primes required!" << std::endl;);
			}
		}

		if (qo_special) {
			std::cout << " " << i-1 << std::flush;
		}

		delete[] PRIM;
	}
}



void
matrix< bigint >::size_red_jacobs()
{
	bigint_lattice L(*this);

	L.lll(0.99);
	*this = L;

#if 0
	register lidia_size_t i, j, k;
	bigint q, u, temp;
	math_vector< bigint > D, v1, v2;
	base_vector< int > f;
	matrix< bigint > L;

	D.set_mode(EXPAND);
	f.set_mode(EXPAND);
	L.resize(columns, columns);

	D[0].assign(1);
	v1 = get_column_vector(0);
	LiDIA::multiply(temp, v1, v1);
	if (temp.is_zero()) {
		D[1].assign_one();
		f[0] = 0;
	}
	else {
		D[1].assign(temp);
		f[0] = 1;
	}


	for (k = 1; k < columns; ++k) {
		// incremental Gram-Schmidt
		for (j = 0; j <= k; ++j) {
			if (!f[j] && (j < k))
				L.sto(k, j, bigint(0));
			else {
				v1 = get_column_vector(k);
				v2 = get_column_vector(j);
				LiDIA::multiply(u, v1, v2);

				for (i = 0; i < j; ++i) {
					if (f[i]) {
						LiDIA::multiply(u, u, D[i+1]);
						LiDIA::multiply(temp, L.member(k, i), L.member(j, i));
						LiDIA::subtract(u, u, temp);
						LiDIA::divide(u, u, D[i]);
					}
				}

				if (j < k)
					L.sto(k, j, u);
				else {
					if (!u.is_zero()) {
						D[k+1] = u;
						f[k] = 1;
					}
					else {
						D[k+1] = D[k];
						f[k] = 0;
					}
				}
			}
		}


		// size reduction
		for (i = k-1; i >= 0; --i) {
			if (f[i]) {
				shift_left(temp, L.member(k, i), 1);
				if (abs(temp) > D[i+1]) {
					LiDIA::add(temp, temp, D[i+1]);
					LiDIA::divide(q, temp, (D[i+1] << 1));
					if (L.member(k, i).is_lt_zero())
						dec(q);

					if (!q.is_zero()) {
						for (j = 0; j < rows; ++j) {
							LiDIA::subtract(temp, member(j, k), q*member(j, i));
							sto(j, k, temp);
						}

						LiDIA::subtract(temp, L.member(k, i), q*D[i+1]);
						L.sto(k, i, temp);
						for (j = 0; j < i; ++j) {
							LiDIA::subtract(temp, L.member(k, j), q*L.member(i, j));
							L.sto(k, j, temp);
						}
					}
				}
			}
		}

	}

	for (i = 0; i < columns-1; ++i) {
		if (f[i]) {
			std::cout << "u[ " << i+1 << ", " << i << " ] = ";
			std::cout << abs(bigfloat(L.member(i+1, i)) / bigfloat(D[i+1])) << std::endl;
		}
		else {
			std::cout << "col[ " << i << " ] is dependent!" << std::endl;
			for (j = 0; j < rows; ++j)
				sto(j, i, bigint(0));
		}
	}

	if (!f[columns-1]) {
		std::cout << "col[ " << columns-1 << " ] is dependent!" << std::endl;
		for (j = 0; j < rows; ++j)
			sto(j, columns-1, bigint(0));
	}
#endif
}



void
matrix< bigint >::size_red_jacobs(trans_matrix & TR)
{
	bigint_lattice L(*this);
	matrix< bigint > tran;

	tran.resize(columns, columns);
	tran.diag(1, 0);

	L.lll(tran, 0.99);
	*this = L;

	TR.store_matrix(tran);

#if 0
	register lidia_size_t i, j, k;
	bigint q, u, temp;
	math_vector< bigint > D, v1, v2;
	base_vector< int > f;
	matrix< bigint > L, tran;

	D.set_mode(EXPAND);
	f.set_mode(EXPAND);
	L.resize(columns, columns);
	tran.resize(columns, columns);
	tran.diag(1, 0);

	D[0].assign(1);
	v1 = get_column_vector(0);
	LiDIA::multiply(temp, v1, v1);
	if (temp.is_zero()) {
		D[1].assign_one();
		f[0] = 0;
	}
	else {
		D[1].assign(temp);
		f[0] = 1;
	}


	for (k = 1; k < columns; ++k) {
		// incremental Gram-Schmidt
		for (j = 0; j <= k; ++j) {
			if (!f[j] && (j < k))
				L.sto(k, j, bigint(0));
			else {
				v1 = get_column_vector(k);
				v2 = get_column_vector(j);
				LiDIA::multiply(u, v1, v2);

				for (i = 0; i < j; ++i) {
					if (f[i]) {
						LiDIA::multiply(u, u, D[i+1]);
						LiDIA::multiply(temp, L.member(k, i), L.member(j, i));
						LiDIA::subtract(u, u, temp);
						LiDIA::divide(u, u, D[i]);
					}
				}

				if (j < k)
					L.sto(k, j, u);
				else {
					if (!u.is_zero()) {
						D[k+1] = u;
						f[k] = 1;
					}
					else {
						D[k+1] = D[k];
						f[k] = 0;
					}
				}
			}
		}


		// size reduction
		for (i = k-1; i >= 0; --i) {
			if (f[i]) {
				shift_left(temp, L.member(k, i), 1);
				if (abs(temp) > D[i+1]) {
					LiDIA::add(temp, temp, D[i+1]);
					LiDIA::divide(q, temp, (D[i+1] << 1));
					if (L.member(k, i).is_lt_zero())
						dec(q);

					if (!q.is_zero()) {
						for (j = 0; j < rows; ++j) {
							LiDIA::subtract(temp, member(j, k), q*member(j, i));
							sto(j, k, temp);
						}

						for (j = 0; j < columns; ++j) {
							LiDIA::subtract(temp, tran.member(j, k), q*tran.member(j, i));
							tran.sto(j, k, temp);
						}

						LiDIA::subtract(temp, L.member(k, i), q*D[i+1]);
						L.sto(k, i, temp);
						for (j = 0; j < i; ++j) {
							LiDIA::subtract(temp, L.member(k, j), q*L.member(i, j));
							L.sto(k, j, temp);
						}
					}
				}
			}
		}

	}

	TR.store_matrix(tran);

	for (i = 0; i < columns-1; ++i) {
		if (f[i]) {
			std::cout << "u[ " << i+1 << ", " << i << " ] = ";
			std::cout << abs(bigfloat(L.member(i+1, i)) / bigfloat(D[i+1])) << std::endl;
		}
		else {
			std::cout << "col[ " << i << " ] is dependent!" << std::endl;
			for (j = 0; j < rows; ++j)
				sto(j, i, bigint(0));
		}
	}

	if (!f[columns-1]) {
		std::cout << "col[ " << columns-1 << " ] is dependent!" << std::endl;
		for (j = 0; j < rows; ++j)
			sto(j, columns-1, bigint(0));
	}
#endif
}



void
pre_reduction(matrix< long > & A,
              lidia_size_t actual_row,
              lidia_size_t actual_column,
              base_vector< lidia_size_t > & row_vec)
{
	register lidia_size_t i, j;
	lidia_size_t row, col;
	bool single, ones;

	// initialize vector of row indices
	row_vec.set_capacity(actual_row);
	row_vec.set_size(actual_row);
	for (i = 0; i < actual_row; ++i)
		row_vec[i] = i;

	// change orientation
	A.set_orientation(matrix_flags::row_oriented);

	// find all rows with only a single +1 or -1
	for (i = actual_row-1; i >= 0; --i) {
		single = true;
		col = -1;
		for (j = 0; j < actual_column && single; ++j) {
			if (std::abs(A.member(i, j)) > 1)
				single = false;
			else if (std::abs(A.member(i, j)) == 1) {
				if (col < 0)
					col = j;
				else
					single = false;
			}
		}

		if (single && (col >= 0)) {
			A.swap_rows(actual_row-1, i);
			LiDIA::swap(row_vec[actual_row-1], row_vec[i]);

			--actual_row;
		}
	}


	// find all rows with only 0, 1, or -1 entries
	row = actual_row-1;
	for (i = row; i >= 0; --i) {
		ones = true;
		for (j = 0; j < actual_column && ones; ++j)
			ones = (std::abs(A.member(i, j)) <= 1);

		if (ones) {
			A.swap_rows(row, i);
			LiDIA::swap(row_vec[row-1], row_vec[i]);
			--row;
		}
	}
}



void
pre_reduction(matrix< long > & A,
              trans_matrix & TR,
              lidia_size_t actual_row,
              lidia_size_t actual_column,
              base_vector< lidia_size_t > & row_vec)
{
	register lidia_size_t i, j;
	lidia_size_t row, col;
	bool single, ones;
	matrix< bigint > tran;

	tran.resize(actual_column, actual_column);
	tran.diag(1, 0);

	// initialize vector of row indices
	row_vec.set_capacity(actual_row);
	row_vec.set_size(actual_row);
	for (i = 0; i < actual_row; ++i)
		row_vec[i] = i;

	// change orientation
	A.set_orientation(matrix_flags::row_oriented);


	// find all rows with only a single +1 or -1
	for (i = actual_row-1; i >= 0; --i) {
		single = true;
		col = -1;
		for (j = 0; j < actual_column && single; ++j) {
			if (std::abs(A.member(i, j)) > 1)
				single = false;
			else if (std::abs(A.member(i, j)) == 1) {
				if (col < 0)
					col = j;
				else
					single = false;
			}
		}

		if (single && (col >= 0)) {
			A.swap_rows(actual_row-1, i);
			LiDIA::swap(row_vec[actual_row-1], row_vec[i]);

			tran.swap_rows(actual_row-1, i);

			--actual_row;
		}
	}


	// find all rows with only 0, 1, or -1 entries
	row = actual_row-1;
	for (i = row; i >= 0; --i) {
		ones = true;
		for (j = 0; j < actual_column && ones; ++j)
			ones = (std::abs(A.member(i, j)) <= 1);

		if (ones) {
			A.swap_rows(row, i);
			LiDIA::swap(row_vec[row], row_vec[i]);

			tran.swap_rows(row, i);

			--row;
		}
	}

	TR.store_matrix(tran);
}



void
matrix< bigint >::post_reduction(base_vector< lidia_size_t > & row_vec,
				  lidia_size_t *RET)
{
	register lidia_size_t i, j;
	bool found;

	if (RET) {
		for (i = 1; i <= RET[0]; ++i)
			RET[i] = row_vec[RET[i]];
	}

	// reswap rows back to original order
	for (i = 0; i < rows; ++i) {
		if (row_vec[i] != i) {
			found = false;
			for (j = i+1; j < rows && !found; ++j) {
				if (row_vec[j] == i) {
					found = true;
					swap_rows(i, j);
				}
			}
		}
	}
}



lidia_size_t *
matrix< bigint >::hnf_jacobs0(trans_matrix & TR, bigint & H, int num_same)
{
	//
	// modul definitions
	//

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2;

	const modular_arithmetic< DRMK< bigint >, dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	//
	// Variables
	//

	register lidia_size_t i, j, k;
	int num_mod;
	bigint DET;
	lidia_size_t *RET = NULL;
	lidia_size_t r, *linu, oldcols;



	// hnf with dkt, transformation matrix with adjoint
	matrix< bigint > tran;
	tran.set_storage_mode(matrix_flags::dense_representation);
	tran.set_orientation(matrix_flags::row_oriented);
	tran.resize(columns, columns);
	tran.diag(1, 0);


	// move linear independant columns to front
	linu = lininc(H);
	r = linu[0];
	if (r != rows) {
		lidia_qo_info_handler(
			std::cout << "Submatrix is singular - modular reduction not performed";
			std::cout << std::endl;);

		lidia_size_t *linr = lininr(H);
		RET = new lidia_size_t[rows];
		j = 0;
		RET[0] = 0;
		for (i = 0; i < linr[0]; ++i) {
			while (linr[linr[0]-i] != j) {
				++RET[0];
				RET[RET[0]] = j;
				++j;
			}
			++j;
		}

		delete [] linr;
		delete [] linu;

		return RET;
	}


	for (i = 0; i < r; i++) {
		swap_columns(i, linu[r-i]);
		tran.swap_columns(i, linu[r-i]);
	}
	delete [] linu;

	TR.store_matrix(tran);
	tran.reset();

	// use hnfmod_mueller on full rank part
	oldcols = columns;
	columns = r;


	// compute determinant
	matrix< bigint > OLD(*this);

	num_mod = Dm_bigint_modul.det(*this, DET, H, num_same);

	num_mod -= 4;
	lidia_qo_info_handler(std::cout << "Det = " << DET << std::endl;);

	if (qo_special)
		std::cout << " 0 " << std::flush;

	// compute hnf (with mod_dkt)
	lidia_qo_info_handler(std::cout << "Computing modular hnf (dkt)..." << std::endl;);
	D_bigint_modul.hnfmod_dkt_part(*this, abs(DET));

	normalization_kernel< bigint, RODMM< bigint >, RODMM< bigint > > normalize_dmodul3;

	normalize_dmodul3.normalizeHybrid_ChouCollins(*this, 0, columns - rows);

	// compute adj
	tran.resize(columns, columns);
	Dm_bigint_modul.adj2(tran, OLD, H, DET, num_mod);
	OLD.reset();

	// compute transformation matrix
	bigint val, akj;
	math_vector< bigint > curr_row;
	curr_row.set_mode(EXPAND);
	for (i = 0; i < rows; ++i) {
		tran.get_row_vector(curr_row, i);

		for (j = 0; j < columns; ++j) {
			val.assign_zero();
			for (k = 0; k < rows; ++k) {
				akj = member(k, j);
				if (!akj.is_zero())
					LiDIA::add(val, val, akj*curr_row[k]);
			}
			tran.sto(i, j, val/DET);
		}
	}

	TR.store_matrix(tran);


	if (qo_special) {
		std::cout << std::endl;
	}



	// use havas on rectangular matrix with first r columns in hnf
	columns = oldcols;
	tran.set_storage_mode(matrix_flags::dense_representation);
	tran.set_orientation(matrix_flags::row_oriented);
	tran.resize(columns, columns);
	tran.diag(1, 0);

	bigint BOUND2;
	BOUND2.assign_zero();
	lidia_size_t actual_row = rows, actual_column = columns;
	lidia_qo_info_handler(std::cout << "Computing stf..." << std::endl;);
	hnf_dmodul2.stf(*this, tran, actual_row, actual_column, BOUND2);

	TR.store_matrix(tran);
	tran.reset();

	// reduce
	tran.set_storage_mode(matrix_flags::dense_representation);
	tran.set_orientation(matrix_flags::row_oriented);
	tran.resize(columns, columns);
	tran.diag(1, 0);

	lidia_qo_info_handler(std::cout << "\nReducing (dense)..." << std::endl;);

	normalize_dmodul3.normalizeHybrid_ChouCollins(*this, TR, tran, columns - rows);

	TR.store_matrix(tran);

	return RET;
}



lidia_size_t *
matrix< bigint >::hnf_jacobs1(trans_matrix & TR, bigint & H, int num_same,
			       bigint & DET, timer & t)
{
	//
	// modul definitions
	//

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2;

	const modular_arithmetic< DRMK< bigint >, dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	//
	// Variables
	//

	register lidia_size_t i, j;
	int num_mod;
	bigint DET2, G;
	lidia_size_t *RET = NULL;
	matrix< bigint > tran;


	// make sure matrix has full row-rank!
	lidia_size_t *linr = lininr(H);
	if (linr[0] < rows) {
		lidia_qo_info_handler(
			std::cout << "Submatrix is singular - modular reduction not performed";
			std::cout << std::endl;);

		RET = new lidia_size_t[rows];
		j = 0;
		RET[0] = 0;
		for (i = 0; i < linr[0]; ++i) {
			while (linr[linr[0]-i] != j) {
				++RET[0];
				RET[RET[0]] = j;
				++j;
			}
			++j;
		}

		delete [] linr;

		return RET;
	}



	// compute multiple of lattice determinant
	lidia_size_t *linu = lininc(H);
	lidia_size_t r = linu[0];
	for (i = 0; i < r; i++)
		swap_columns(i, linu[r-i]);

	lidia_size_t old_columns = columns;
	columns = linu[0];

	num_mod = Dm_bigint_modul.det(*this, DET, H, num_same);
	num_mod -= 4;

	columns = old_columns;
	for (i = r-1; i >= 0; i--)
		swap_columns(i, linu[r-i]);
	for (i = 0; i < columns/2; i++)
		swap_columns(i, columns-i-1);

	lidia_size_t *linu2 = lininc(H);
	for (i = 0; i < r; i++)
		swap_columns(i, linu2[r-i]);
	columns = linu2[0];

	Dm_bigint_modul.det(*this, DET2, H, num_same);

	columns = old_columns;
	for (i = r-1; i >= 0; i--)
		swap_columns(i, linu2[r-i]);
	for (i = 0; i < columns/2; i++)
		swap_columns(i, columns-i-1);

	delete[] linu2;
	DET2 = gcd(DET, DET2);

	if (DET2.is_lt_zero())
		DET2.negate();
	lidia_qo_info_handler(std::cout << "Det mult = " << DET2 << std::endl;);

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << "determinant: " << t << std::endl << std::endl;
			      t.cont_timer(););

	// compute copy of *this
	matrix< bigint > A(*this);


	// compute hnf and transformation matrix
	lidia_qo_info_handler(std::cout << "Computing modular hnf..." << std::endl;);
	D_bigint_modul.hnfmod_dkt_part(*this, DET2);

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << "hnfmod_dkt: " << t << std::endl << std::endl;
			      t.cont_timer(););

	normalization_kernel< bigint, RODMM< bigint >, RODMM< bigint > > normalize_dmodul3;
	normalize_dmodul3.normalizeHybrid_ChouCollins(*this, 0, columns - rows);


	// compute regular form of A
	if (linu[0] != rows) {
		std::cout << "ERROR:  not full rank submatrix!!!  PANIC!\n" << std::flush;
		exit(1);
	}

	tran.resize(columns, columns);
	tran.diag(1, 0);
	for (i = 0; i < rows; ++i) {
		A.swap_columns(i, linu[rows-i]);
		tran.swap_columns(i, linu[rows-i]);
	}
	delete [] linu;


	TR.store_matrix(tran);
	tran.reset();


	// store dependent columns in *this (temporarily)
	for (i = 0; i < rows; ++i)
		for (j = rows; j < columns; ++j)
			sto(i, j-rows, A.member(i, j));


	// compute adjoint
	A.set_no_of_columns(rows);
	lidia_qo_info_handler(std::cout << "DET(A) = " << DET << std::endl;);
	if (TR.get_adjoint_mode()) {
		// use files for adjoint
		TR.ADJ.init(rows, rows, static_cast<lidia_size_t>(std::ceil(rows/25.0)));
		A.adj(TR.ADJ, A.hadamard(), DET, num_mod);
		TR.store_adjoint();

		lidia_qo_info_handler(t.stop_timer();
				      std::cout << "adjoint: " << t << std::endl;
				      t.cont_timer(););
	}
	else {
		tran.resize(rows, rows);
		Dm_bigint_modul.adj2(tran, A, A.hadamard(), DET, num_mod);

		lidia_qo_info_handler(t.stop_timer();
				      std::cout << "adjoint: " << t << std::endl;
				      t.cont_timer(););

		TR.store_matrix(tran);
		tran.reset();
	}
	A.reset();

	// compute extended hnf
	tran = *this;
	tran.resize(columns, columns);
	for (i = rows; i < columns; ++i)
		tran.sto(i, i-rows, -DET);

	TR.store_matrix(tran, DET);

	// remove dependent columns of A from *this
	for (i = 0; i < rows; ++i)
		for (j = 0; j < columns-rows; ++j)
			sto(i, j, bigint(0));


	if (qo_special) {
		std::cout << std::endl;
	}

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg2(const matrix< long > &B, int no, bool & do_mod)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;

	havas_kernel< long, COSMM< long >, COSMM< bigint >,
		havas_best_remainder_ext< long, COSMM< long >, COSMM< bigint > > > hnf_smodul2le;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	//
	// Variables
	//

	timer t;
	lidia_qo_info_handler(t.start_timer(););

	lidia_size_t actual_column;
	lidia_size_t actual_row = B.rows;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	register lidia_size_t i, j;

	int BOUND_1;
	bigint BOUND_2, DET, H;
	BOUND_1 = 2147483647;
	power(BOUND_2, bigint(BOUND_1), 2);

	//
	// main algorithm
	//

	// Change representation of matrix B
	set_representation(matrix_flags::sparse_representation);

	// Change dimensions of member matrix
	resize(B.rows, B.columns);

	if (qo_special) {
		std::cout << B.rows << std::flush;
	}

	do_mod = false;

	//
	// Stage 1: non-modular long
	//
	{
		// Create a copy of the original matrix
		matrix< long > A;
		A.set_representation(matrix_flags::sparse_representation);
		A.set_zero_element(0);
		A.resize(B.rows, B.columns);
		for (i = 0; i < B.rows; i++) {
			lidia_size_t len = B.value_counter[i];
			A.value[i] = new long[len];
			A.index[i] = new lidia_size_t[len];
			for (j = 0; j < len; j++) {
				A.value[i][j] = B.value[i][j];
				A.index[i][j] = B.index[i][j];
			}
			A.value_counter[i] = len;
			A.allocated[i] = len;
		}

		// change orientation
		A.set_orientation(matrix_flags::column_oriented);

		// Compute the Hadamard Bound
		modul3.hadamard(A, H);
		lidia_qo_info_handler(t.stop_timer();
				      std::cout << std::endl << "Hadamard Bound = " << decimal_length(H) << " decimal digits" << std::endl;
				      std::cout << "Time: " << t << std::endl;
				      std::cout << "long reduction (row " << actual_row << ")..." << std::endl;
				      t.cont_timer(););

		// STF computation
		if (hnf_smodul2le.stf(A, actual_row, actual_column, BOUND_1)) {
			set_orientation(matrix_flags::row_oriented);

			// copy data
			for (i = 0; i < columns_orig; i++)
				for (j = 0; j < A.value_counter[i]; j++)
					sto(A.index[i][j], i, bigint(A.value[i][j]));

			lidia_qo_info_handler(std::cout << "Reducing...\n" << std::flush;);
			set_orientation(matrix_flags::column_oriented);

			normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;
			RET = normalize_smodul2.normalizeHybrid_ChouCollins(*this, 0, actual_column);

			set_orientation(matrix_flags::row_oriented);

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "total: " << t << std::endl;);

			if (qo_special) {
				std::cout << " 0 0 0 0 " << std::endl;
			}

			return RET;
		}

		// copying matrix
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < columns_orig; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < A.value_counter[i]; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}
	}

	//
	// Stage 2: non-modular bigint
	//

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "stf (long): " << t << std::endl;
			      std::cout << "\nbigint reduction (row " << actual_row << ")..." << std::endl;
			      t.cont_timer(););

	if (qo_special) {
		std::cout << " " << actual_row << std::flush;
	}


	{
		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
		}

		bool OK = hnf_dmodul2e.stf(A3, actual_row, actual_column, BOUND_2);
		lidia_size_t A3rows = A3.rows, A3cols = A3.columns;

		if (!OK && (0)) {
			// use LLL size reduction

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "stf (bigint): " << t << std::endl;
					      std::cout << "\nsize_red (row " << actual_row << ")...";
					      std::cout << std::endl;
					      t.cont_timer(););

			A3.rows = actual_row;
			A3.columns = actual_column;

			A3.size_red_jacobs();

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "size_reduce (bigint): " << t;
					      std::cout << std::endl;
					      t.cont_timer(););

			square(BOUND_2, BOUND_2);
			OK = hnf_dmodul2e.stf(A3, actual_row, actual_column, BOUND_2);

			if (!OK) {
				j = 0;
				for (i = 0; i < A3.columns; ++i) {
					if (!A3.is_column_zero(i)) {
						A3.swap_columns(i, j);
						++j;
					}
				}

				actual_column = j;
			}
		}

		if (!OK) {
			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "stf (bigint): " << t << std::endl;
					      std::cout << "\nhnfmod (row " << actual_row << ")...";
					      std::cout << std::endl;
					      t.cont_timer(););

			if (qo_special) {
				std::cout << " " << actual_row << std::flush;
			}

			A3.rows = actual_row;
			A3.columns = actual_column;

			lidia_size_t *linr = Dm_bigint_modul.latticedet5(A3, DET, H, 5);
			if (!linr) {
				D_bigint_modul.hnfmod_dkt_part(A3, DET);

				if (qo_special) {
					std::cout << std::endl;
				}
			}
			else {
				do_mod = true;

				lidia_qo_info_handler(
					std::cout << "Submatrix is singular - modular reduction not performed";
					std::cout << std::endl;);

				RET = new lidia_size_t[A3.rows];
				j = 0;
				RET[0] = 0;
				for (i = 0; i < linr[0]; ++i) {
					while (linr[linr[0]-i] != j) {
						++RET[0];
						RET[RET[0]] = j;
						++j;
					}
					++j;
				}

				delete [] linr;
			}

			A3.rows = A3rows;
			A3.columns = A3cols;
		}

		else {
			A3.rows = A3rows;
			A3.columns = A3cols;

			if (qo_special) {
				std::cout << " 0 0 0" << std::endl;
			}
		}



		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}
			value_counter[i] = k;
		}

		A3.reset();
	}

	// reduce
	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;
	if (!RET)
		RET = normalize_smodul2.normalizeHybrid_ChouCollins(*this, 0, B.columns - B.rows);


	set_orientation(matrix_flags::row_oriented);
	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg2(const matrix< long > &B, trans_matrix &TR, int no, bigint & DET,
        int strat, bool & do_mod)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;
	COSMM< bigint > modul4;

	havas_kernel< long, COSMM< long >, COSMM< bigint >,
		havas_best_remainder_ext< long, COSMM< long >, COSMM< bigint > > > hnf_smodul2le;

	havas_kernel< long, COSMM< long >, RODMM< bigint >,
		havas_best_remainder_ext< long, COSMM< long >, RODMM< bigint > > > hnf_smodul2le_d;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;

	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;


	//
	// variables
	//

	matrix< bigint > tran;

	timer t;
	lidia_qo_info_handler(t.start_timer(););

	lidia_size_t actual_column;
	lidia_size_t actual_row = B.rows;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	register lidia_size_t i, j;

	int BOUND_1;
	bigint BOUND_2, H;
	BOUND_1 = 2147483647;
	power(BOUND_2, bigint(BOUND_1), 2);


	//
	// main algorithm
	//

	// Change representation of matrix B
	set_representation(matrix_flags::sparse_representation);

	// Change dimensions of member matrix
	resize(B.rows, B.columns);

	if (qo_special) {
		std::cout << B.rows << std::flush;
	}

	DET = 1;
	do_mod = false;


	//
	// Stage 1: non-modular long
	//
	{
		// Create a copy of the original matrix
		matrix< long > A;
		A.set_representation(matrix_flags::sparse_representation);
		A.set_zero_element(0);
		A.resize(B.rows, B.columns);
		for (i = 0; i < B.rows; i++) {
			lidia_size_t len = B.value_counter[i];
			A.value[i] = new long[len];
			A.index[i] = new lidia_size_t[len];
			for (j = 0; j < len; j++) {
				A.value[i][j] = B.value[i][j];
				A.index[i][j] = B.index[i][j];
			}
			A.value_counter[i] = len;
			A.allocated[i] = len;
		}

		// change orientation
		A.set_orientation(matrix_flags::column_oriented);

		// Compute the Hadamard Bound
		modul3.hadamard(A, H);
		lidia_qo_info_handler(t.stop_timer();
				      std::cout << std::endl << "Hadamard Bound = " << decimal_length(H)
				      << " decimal digits" << std::endl;
				      std::cout << "Time: " << t << std::endl;
				      std::cout << "long reduction (row " << actual_row << ")..." << std::endl;
				      t.cont_timer(););

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::column_oriented);
		tran.resize(B.columns, B.columns);
		tran.diag(1, 0);


		// STF computation
		int SW1 = hnf_smodul2le.hnf_Z2(A, TR, tran, actual_row, actual_column, B.rows,
					       BOUND_1);

		if (actual_row < B.rows)
			TR.store_matrix(tran);
		tran.reset();

		if (SW1 < 0) {
			// more long reduction with bigint transformation matrices

			tran.set_storage_mode(matrix_flags::dense_representation);
			tran.set_orientation(matrix_flags::row_oriented);
			tran.resize(actual_column, actual_column);
			tran.diag(1, 0);

			A.rows = actual_row;
			A.columns = actual_column;

			SW1 = hnf_smodul2le_d.hnf_Z2(A, TR, tran, actual_row, actual_column, B.rows,
						     BOUND_1);

			A.rows = B.rows;
			A.columns = B.columns;

			TR.store_matrix(tran);
			tran.reset();
		}

		if (SW1 > 0) {
			set_orientation(matrix_flags::row_oriented);

			// copy data
			for (i = 0; i < columns_orig; i++)
				for (j = 0; j < A.value_counter[i]; j++)
					sto(A.index[i][j], i, bigint(A.value[i][j]));

			lidia_qo_info_handler(std::cout << "Reducing...\n" << std::flush;);

			set_orientation(matrix_flags::column_oriented);

			tran.set_storage_mode(matrix_flags::sparse_representation);
			tran.set_orientation(matrix_flags::column_oriented);
			tran.resize(B.columns, B.columns);
			tran.diag(1, 0);

//        RET = normalize_smodul2.normalizeHybrid_ChouCollins(*this, tran, 0, actual_column);
			RET = normalize_smodul2.normalize_ChouCollins(*this, TR, tran, 0, B.columns - B.rows);

			TR.store_matrix(tran);

			set_orientation(matrix_flags::row_oriented);

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "total: " << t << std::endl;);

			if (qo_special) {
				std::cout << " 0 0 0 0 0" << std::endl;
			}

			return RET;
		}

		// copying matrix
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < columns_orig; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < A.value_counter[i]; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}
	}

	//
	// Stage 2: non-modular bigint
	//

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "stf (long): " << t << std::endl;
			      std::cout << "\nbigint reduction (row " << actual_row << ")..." << std::endl;
			      t.cont_timer(););

	if (qo_special) {
		std::cout << " " << actual_row << std::flush;
	}

	{
		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
		}

		A3.set_orientation(matrix_flags::row_oriented);

		tran.set_storage_mode(matrix_flags::dense_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.resize(actual_column, actual_column);
		tran.diag(1, 0);

		bool OK = hnf_dmodul2e.hnf_Z2(A3, TR, tran, actual_row, actual_column,
					      B.rows, BOUND_2);
		TR.store_matrix(tran);
		tran.reset();

		if (!OK && (0)) {
			// use LLL size reduction

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "stf (bigint): " << t << std::endl;
					      std::cout << "\nsize_red (row " << actual_row << ")...";
					      std::cout << std::endl;
					      t.cont_timer(););

			// insert finished columns into result
			lidia_size_t nonzero;
			for (i = actual_column; i < A3.columns; ++i) {
				nonzero = 0;
				for (j = 0; j < A3.rows; ++j)
					if (A3.value[j][i] != A3.Zero)
						++nonzero;

				value[i] = new bigint[nonzero+5];
				index[i] = new lidia_size_t[nonzero+5];
				allocated[i] = nonzero+5;

				lidia_size_t k = 0;
				for (j = 0; j < A3.rows; j++) {
					if (A3.value[j][i] != A3.Zero) {
						value[i][k] = A3.value[j][i];
						index[i][k] = j;
						k++;
					}
				}

				value_counter[i] = k;
			}

			A3.resize(actual_row, actual_column);

			A3.size_red_jacobs(TR);

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "size_reduce (bigint): " << t;
					      std::cout << std::endl;
					      t.cont_timer(););

			square(BOUND_2, BOUND_2);

			tran.set_storage_mode(matrix_flags::dense_representation);
			tran.set_orientation(matrix_flags::row_oriented);
			tran.resize(actual_column, actual_column);
			tran.diag(1, 0);

			OK = hnf_dmodul2e.hnf_Z2(A3, TR, tran, actual_row, actual_column,
						 B.rows, BOUND_2);

			TR.store_matrix(tran);

			if (!OK) {
				tran.set_storage_mode(matrix_flags::dense_representation);
				tran.set_orientation(matrix_flags::row_oriented);
				tran.resize(actual_column, actual_column);
				tran.diag(1, 0);

				j = 0;
				for (i = 0; i < A3.columns; ++i) {
					if (!A3.is_column_zero(i)) {
						A3.swap_columns(i, j);
						tran.swap_columns(i, j);
						++j;
					}
				}

				actual_column = j;
				TR.store_matrix(tran);
			}
		}

		if (!OK) {
			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "stf (bigint): " << t << std::endl;
					      std::cout << "\nhnfmod (row " << actual_row << ")...";
					      std::cout << std::endl;
					      t.cont_timer(););

			// insert finished columns into result
			lidia_size_t nonzero;
			for (i = actual_column; i < A3.columns; ++i) {
				nonzero = 0;
				for (j = 0; j < A3.rows; ++j)
					if (A3.value[j][i] != A3.Zero)
						++nonzero;

				value[i] = new bigint[nonzero+5];
				index[i] = new lidia_size_t[nonzero+5];
				allocated[i] = nonzero+5;

				lidia_size_t k = 0;
				for (j = 0; j < A3.rows; j++) {
					if (A3.value[j][i] != A3.Zero) {
						value[i][k] = A3.value[j][i];
						index[i][k] = j;
						k++;
					}
				}

				value_counter[i] = k;
			}

			A3.resize(actual_row, actual_column);

			if (qo_special) {
				std::cout << " " << actual_row << std::flush;
			}


			switch (strat) {
			case 0:  // dkt on regular square, havas to eliminate rest
				RET = A3.hnf_jacobs0(TR, H, no);
				break;

			default:  // hnfmod_dkt, transformation computed with adjoint
				RET = A3.hnf_jacobs1(TR, H, no, DET, t);
				break;
			}

			if (RET)
				do_mod = true;

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "hnfmod: " << t << std::endl;
					      t.cont_timer(););
		}

		else {
			if (qo_special) {
				std::cout << " 0 0 0 0" << std::endl;
			}
		}




		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}

			value_counter[i] = k;
		}

		A3.reset();
	}

	set_orientation(matrix_flags::column_oriented);

	// reduce
	lidia_qo_info_handler(std::cout << "reducing..." << std::endl;);

	if (!RET) {
		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::column_oriented);
		tran.resize(B.columns, B.columns);
		tran.diag(1, 0);

//    RET = normalize_smodul2.normalizeHybrid_ChouCollins(*this, TR, tran, B.columns - B.rows);
		RET = normalize_smodul2.normalize_ChouCollins(*this, TR, tran, 0, B.columns - B.rows);

		TR.store_matrix(tran);
		tran.reset();
	}

	set_orientation(matrix_flags::row_oriented);

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg2(const matrix< long > &B, int no)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;

	havas_kernel< long, COSMM< long >, COSMM< bigint >,
		nf_conf3e< long, COSMM< long >, COSMM< bigint > > > hnf_smodul2le;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;


	//
	// Variables
	//

	timer t;
	lidia_qo_info_handler(t.start_timer(););

	lidia_size_t actual_row = B.rows;
	lidia_size_t actual_column = B.columns;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	register lidia_size_t i, j;

	int BOUND_1;
	bigint BOUND_2, H, DET;
	BOUND_1 = 2147483647;

	power(BOUND_2, bigint(BOUND_1), 2);

	//
	// main algorithm
	//

	// Change representation of matrix B
	set_representation(matrix_flags::sparse_representation);

	// Change dimensions of member matrix
	resize(B.rows, B.columns);

	if (qo_special) {
		std::cout << B.rows << std::flush;
	}

	//
	// Stage 1: non-modular long
	//
	{
		// Create a copy of the original matrix
		matrix< long > A = B;
		A.set_zero_element(0);

		// change orientation
		A.set_orientation(matrix_flags::column_oriented);

		// Compute the Hadamard Bound
		modul3.hadamard(A, H);
		lidia_qo_info_handler(t.stop_timer();
				      std::cout << std::endl << "Hadamard Bound = " << decimal_length(H) <<
				      " decimal digits" << std::endl;
				      std::cout << "Time: " << t << std::endl;
				      std::cout << "long reduction (row " << actual_row << ")..."
				      << std::endl;
				      t.cont_timer(););

		// STF computation
		if (hnf_smodul2le.stf(A, actual_row, actual_column, BOUND_1)) {
			std::cout << "Version 1" << A << std::endl;
			set_orientation(matrix_flags::row_oriented);

			std::cout << "Version 2" << A << std::endl;
			// copy data
			for (i = 0; i < columns_orig; i++)
				for (j = 0; j < A.value_counter[i]; j++)
					sto(A.index[i][j], i, bigint(A.value[i][j]));

			std::cout << "Version 3" << *this << std::endl;
			lidia_qo_info_handler(std::cout << "Reducing...\n" << std::flush;);
			set_orientation(matrix_flags::column_oriented);

			std::cout << "Version 4" << *this << std::endl;
			std::cout << "Hier" << std::endl;
			RET = normalize_smodul2.normalize_ChouCollins(*this, 0, actual_column);

			std::cout << "Version 5" << *this << std::endl;

			set_orientation(matrix_flags::row_oriented);

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "total: " << t << std::endl;);

			if (qo_special) {
				std::cout << " 0 0 0 0 " << std::endl;
			}

			return RET;
		}

		// copying matrix
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < columns_orig; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < A.value_counter[i]; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}
	}

	//
	// Stage 2: non-modular bigint
	//

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "stf (long): " << t << std::endl;
			      std::cout << "\nbigint reduction (row " << actual_row << ")..." << std::endl;
			      t.cont_timer(););

	if (qo_special) {
		std::cout << " " << actual_row << std::flush;
	}

	{
		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
		}

		if (!hnf_dmodul2e.stf(A3, actual_row, actual_column, BOUND_2)) {
			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "stf (bigint): " << t << std::endl;
					      std::cout << "\nhnfmod (row " << actual_row << ")..." << std::endl;
					      t.cont_timer(););

			if (qo_special) {
				std::cout << " " << actual_row << std::flush;
			}

			lidia_size_t A3rows = A3.rows, A3cols = A3.columns;
			A3.rows = actual_row;
			A3.columns = actual_column;

			Dm_bigint_modul.latticedet4(A3, DET, H, no);
			D_bigint_modul.hnfmod_dkt_part(A3, DET);

			A3.rows = A3rows;
			A3.columns = A3cols;

			if (qo_special) {
				std::cout << std::endl;
			}
		}
		else {
			if (qo_special) {
				std::cout << " 0 0 0" << std::endl;
			}
		}

		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}
			value_counter[i] = k;
		}

		A3.reset();
	}

	// reduce
	std::cout << "Hier2" << std::endl;
	RET = normalize_smodul2.normalize_ChouCollins(*this, 0, B.columns - B.rows);

	set_orientation(matrix_flags::row_oriented);
	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg2(const matrix< long > &B, trans_matrix &TR, int no, bigint & DET,
        int strat)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;
	COSMM< bigint > modul4;

	havas_kernel< long, COSMM< long >,
		COSMM< bigint >,
		nf_conf3e< long, COSMM< long >,
		COSMM< bigint > > > hnf_smodul3le;

	havas_kernel< long, COSMM< long >,
		RODMM< bigint >,
		nf_conf3e< long, COSMM< long >,
		RODMM< bigint > > > hnf_smodul3le_d;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul3;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul3;

	havas_kernel< bigint, RODMM< bigint >,
		RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >,
		RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul3;

	normalization_kernel< bigint, RODMM< bigint >, RODMM< bigint > > normalize_dmodul3;


	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;


	//
	// variables
	//

	matrix< bigint > tran;

	timer t;
	lidia_qo_info_handler(t.start_timer(););

	lidia_size_t actual_row, actual_column;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	register lidia_size_t i, j;

	int BOUND_1;
	bigint BOUND_2, H;
	BOUND_1 = 2147483647;
	power(BOUND_2, bigint(BOUND_1), 2);


	//
	// main algorithm
	//

	// Change representation of matrix B
	set_representation(matrix_flags::sparse_representation);

	// Change dimensions of member matrix
	resize(B.rows, B.columns);

	if (qo_special) {
		std::cout << B.rows << std::flush;
	}

	DET = 1;


	//
	// Stage 1: non-modular long
	//
	{
		// Create a copy of the original matrix
		matrix< long > A = B;
		A.set_zero_element(0);

		// change orientation
		A.set_orientation(matrix_flags::column_oriented);

		// Compute the Hadamard Bound
		modul3.hadamard(A, H);
		lidia_qo_info_handler(t.stop_timer();
				      std::cout << std::endl << "Hadamard Bound = " << decimal_length(H)
				      << " decimal digits" << std::endl;
				      std::cout << "Time: " << t << std::endl;
				      std::cout << "long reduction (row " << actual_row << ")..." << std::endl;
				      t.cont_timer(););

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::column_oriented);
		tran.resize(B.columns, B.columns);
		tran.diag(1, 0);


		// STF computation
		int SW1 = hnf_smodul3le.hnf_Z2(A, TR, tran, actual_row, actual_column, B.rows,
					       BOUND_1);

		if (actual_row < B.rows)
			TR.store_matrix(tran);
		tran.reset();

		if (SW1 < 0) {
			// more long reduction with bigint transformation matrices

			tran.set_storage_mode(matrix_flags::dense_representation);
			tran.set_orientation(matrix_flags::row_oriented);
			tran.resize(actual_column, actual_column);
			tran.diag(1, 0);

			A.rows = actual_row;
			A.columns = actual_column;

			SW1 = hnf_smodul3le_d.hnf_Z2(A, TR, tran, actual_row, actual_column, B.rows,
						     BOUND_1);

			A.rows = B.rows;
			A.columns = B.columns;

			TR.store_matrix(tran);
			tran.reset();
		}

		if (SW1 > 0) {
			set_orientation(matrix_flags::row_oriented);

			// copy data
			for (i = 0; i < columns_orig; i++)
				for (j = 0; j < A.value_counter[i]; j++)
					sto(A.index[i][j], i, bigint(A.value[i][j]));

			lidia_qo_info_handler(std::cout << "Reducing...\n" << std::flush;);

			set_orientation(matrix_flags::column_oriented);

			tran.set_storage_mode(matrix_flags::sparse_representation);
			tran.set_orientation(matrix_flags::column_oriented);
			tran.resize(B.columns, B.columns);
			tran.diag(1, 0);

			RET = normalize_smodul3.normalize_ChouCollins(*this, tran, 0, actual_column);

			set_orientation(matrix_flags::row_oriented);

			TR.store_matrix(tran);

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "total: " << t << std::endl;);

			if (qo_special) {
				std::cout << " 0 0 0 0 0" << std::endl;
			}

			return RET;
		}

		// copying matrix
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < columns_orig; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < A.value_counter[i]; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}
	}

	//
	// Stage 2: non-modular bigint
	//

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "stf (long): " << t << std::endl;
			      std::cout << "\nbigint reduction (row " << actual_row << ")..." << std::endl;
			      t.cont_timer(););

	if (qo_special) {
		std::cout << " " << actual_row << std::flush;
	}

	{
		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
		}

		A3.set_orientation(matrix_flags::row_oriented);

		tran.set_storage_mode(matrix_flags::dense_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.resize(actual_column, actual_column);
		tran.diag(1, 0);

		bool SW = hnf_dmodul2e.hnf_Z2(A3, TR, tran, actual_row, actual_column,
					      B.rows, BOUND_2);
		TR.store_matrix(tran);
		tran.reset();

		if (!SW) {
			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "stf (bigint): " << t << std::endl;
					      std::cout << "\nhnfmod (row " << actual_row << ")..." << std::endl;
					      t.cont_timer(););

			// insert finished columns into result
			lidia_size_t nonzero;
			for (i = actual_column; i < A3.columns; ++i) {
				nonzero = 0;
				for (j = 0; j < A3.rows; ++j)
					if (A3.value[j][i] != A3.Zero)
						++nonzero;

				value[i] = new bigint[nonzero+5];
				index[i] = new lidia_size_t[nonzero+5];
				allocated[i] = nonzero+5;

				lidia_size_t k = 0;
				for (j = 0; j < A3.rows; j++) {
					if (A3.value[j][i] != A3.Zero) {
						value[i][k] = A3.value[j][i];
						index[i][k] = j;
						k++;
					}
				}

				value_counter[i] = k;
			}
			A3.resize(actual_row, actual_column);


			tran.set_storage_mode(matrix_flags::dense_representation);
			tran.set_orientation(matrix_flags::row_oriented);
			tran.resize(A3.columns, A3.columns);
			tran.diag(1, 0);

			if (qo_special) {
				std::cout << " " << actual_row << std::flush;
			}

			lidia_size_t i, r, *linu, oldcols;

			switch (strat) {
			case 0:  // dkt on regular square, havas to eliminate rest
				// move linear independant columns to front
				linu = A3.lininc(H);
				r = linu[0];
				for (i = 0; i < r; i++) {
					A3.swap_columns(i, linu[r-i]);
					tran.swap_columns(i, linu[r-i]);
				}
				delete [] linu;

				TR.store_matrix(tran);
				tran.reset();

				// use hnfmod_mueller on full rank part
				oldcols = A3.columns;
				A3.columns = r;

				A3.hnf_jacobs0(TR, H, no);

				if (qo_special) {
					std::cout << std::endl;
				}

				// use havas on rectangular matrix with first r columns in hnf
				A3.columns = oldcols;
				tran.set_storage_mode(matrix_flags::dense_representation);
				tran.set_orientation(matrix_flags::row_oriented);
				tran.resize(A3.columns, A3.columns);
				tran.diag(1, 0);

				lidia_qo_info_handler(std::cout << "Computing stf..." << std::endl;);
				hnf_dmodul3.stf(A3, tran, actual_row, actual_column, BOUND_2);

				TR.store_matrix(tran);
				tran.reset();

				// reduce
				tran.set_storage_mode(matrix_flags::dense_representation);
				tran.set_orientation(matrix_flags::row_oriented);
				tran.resize(A3.columns, A3.columns);
				tran.diag(1, 0);

				lidia_qo_info_handler(std::cout << "\nReducing (dense)..." << std::endl;);
				normalize_dmodul3.normalize_ChouCollins(A3, TR, tran, 0, A3.columns - A3.rows);

				TR.store_matrix(tran);
				tran.reset();

				break;

			default:  // hnfmod_dkt, transformation computed with adjoint
				A3.hnf_jacobs1(TR, H, no, DET, t);

				if (qo_special) {
					std::cout << std::endl;
				}

				break;
			}

			lidia_qo_info_handler(t.stop_timer();
					      std::cout << std::endl << "hnfmod: " << t << std::endl;
					      t.cont_timer(););
		}
		else {
			if (qo_special) {
				std::cout << " 0 0 0 0" << std::endl;
			}
		}

		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}

			value_counter[i] = k;
		}

		A3.reset();
	}

	set_orientation(matrix_flags::column_oriented);

	// reduce
	lidia_qo_info_handler(std::cout << "reducing..." << std::endl;);

	if (strat == 3)
		RET = normalize_smodul3.normalize_ChouCollins(*this, 0, B.columns - B.rows);
	else {
		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::column_oriented);
		tran.resize(B.columns, B.columns);
		tran.diag(1, 0);

		RET = normalize_smodul3.normalize_ChouCollins(*this, TR, tran, 0, B.columns - B.rows);

		TR.store_matrix(tran);
		tran.reset();
	}

	set_orientation(matrix_flags::row_oriented);

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg3(const matrix< bigint > &B)
{
	//
	// modul definitions
	//

	COSMM< bigint > modul3;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;


	//
	// variables
	//

	timer t;
	lidia_qo_info_handler(t.start_timer(););


	lidia_size_t  actual_row, actual_column;
	lidia_size_t rows_orig = actual_row = B.rows;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	register lidia_size_t i, j;

	bigint BOUND2;
	matrix< bigint > A;



	// Create a copy of the original matrix
	A.set_representation(matrix_flags::sparse_representation);
	A = B;

	// change orientation
	A.set_orientation(matrix_flags::column_oriented);

	// STF computation
	hnf_smodul2.stf(A, actual_row, actual_column, BOUND2);

	// copying matrix
	set_representation(matrix_flags::sparse_representation);
	kill();
	resize(rows_orig, columns_orig);

	for (i = 0; i < columns_orig; ++i)
		for (j = 0; j < A.value_counter[i]; ++j)
			sto(A.index[i][j], i, bigint(A.value[i][j]));

	// reduce
	lidia_qo_info_handler(std::cout << "Reducing...\n" << std::flush;);
	set_orientation(matrix_flags::column_oriented);
	RET = normalize_smodul2.normalize_ChouCollins(*this, 0, columns - rows);

	set_orientation(matrix_flags::row_oriented);
	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg3_mod(const matrix< bigint > &B, bool & do_mod)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder_ext< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2e;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	//
	// variables
	//

	timer t;
	lidia_qo_info_handler(t.start_timer(););


	lidia_size_t  actual_row, actual_column;
	lidia_size_t rows_orig = actual_row = B.rows;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	bigint DET;
	register lidia_size_t i, j;

	bigint BOUND2;
	matrix< bigint > A;

	square(BOUND2, bigint(2147483647));
	do_mod = false;


	// Create a copy of the original matrix
	A.set_representation(matrix_flags::sparse_representation);
	A = B;

	// change orientation
	A.set_orientation(matrix_flags::column_oriented);

	bigint H = A.hadamard();

	// STF computation
	if (!hnf_smodul2e.stf(A, actual_row, actual_column, BOUND2)) {
		// copying matrix
		set_representation(matrix_flags::sparse_representation);
		kill();
		resize(rows_orig, columns_orig);

		// copying matrix
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < columns_orig; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < A.value_counter[i]; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}


		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
			allocated[i] = 0;
		}

		lidia_qo_info_handler(std::cout << "\nhnfmod (row " << actual_row;
				      std::cout << ")..." << std::endl;);

		lidia_size_t *linr = Dm_bigint_modul.latticedet5(A3, DET, H, 5);
		if (!linr) {
			D_bigint_modul.hnfmod_dkt_part(A3, DET);

			if (qo_special) {
				std::cout << std::endl;
			}
		}
		else {
			do_mod = true;

			lidia_qo_info_handler(
				std::cout << "Submatrix is singular - modular reduction not performed";
				std::cout << std::endl;);

			RET = new lidia_size_t[A3.rows];
			j = 0;
			RET[0] = 0;
			for (i = 0; i < linr[0]; ++i) {
				while (linr[linr[0]-i] != j) {
					++RET[0];
					RET[RET[0]] = j;
					++j;
				}
				++j;
			}
		}

		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}
			value_counter[i] = k;
		}

		A3.reset();
	}
	else {
		// copying matrix
		set_representation(matrix_flags::sparse_representation);
		kill();
		resize(rows_orig, columns_orig);

		for (i = 0; i < columns_orig; ++i)
			for (j = 0; j < A.value_counter[i]; ++j)
				sto(A.index[i][j], i, bigint(A.value[i][j]));
	}

	// reduce

	if (!RET)
		RET = normalize_smodul2.normalize_ChouCollins(*this, 0, B.columns - B.rows);

	set_orientation(matrix_flags::row_oriented);
	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg3(const matrix< bigint > &B, trans_matrix & TR)
{
	//
	// modul definitions
	//

	COSMM< bigint > modul3;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul3;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul3;


	//
	// variables
	//

	timer t;
	lidia_qo_info_handler(t.start_timer(););

	matrix< bigint > tran;


	lidia_size_t  actual_row, actual_column;
	lidia_size_t rows_orig = actual_row = B.rows;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	register lidia_size_t i, j;

	bigint BOUND2;
	matrix< bigint > A;



	// Create a copy of the original matrix
	A.set_representation(matrix_flags::sparse_representation);
	A = B;

	// change orientation
	A.set_orientation(matrix_flags::column_oriented);

	tran.set_storage_mode(matrix_flags::sparse_representation);
	tran.set_orientation(matrix_flags::column_oriented);
	tran.resize(B.columns, B.columns);
	tran.diag(1, 0);


	// STF computation
	hnf_smodul3.stf(A, tran, actual_row, actual_column, BOUND2);

	TR.store_matrix(tran);
	tran.reset();

	// copying matrix
	set_representation(matrix_flags::sparse_representation);
	kill();
	resize(rows_orig, columns_orig);

	for (i = 0; i < columns_orig; ++i)
		for (j = 0; j < A.value_counter[i]; ++j)
			sto(A.index[i][j], i, bigint(A.value[i][j]));

	// reduce
	set_orientation(matrix_flags::column_oriented);

	tran.set_storage_mode(matrix_flags::sparse_representation);
	tran.set_orientation(matrix_flags::column_oriented);
	tran.resize(B.columns, B.columns);
	tran.diag(1, 0);

	if (columns >= rows) {
		RET = normalize_smodul3.normalize_ChouCollins(*this, tran, 0, columns - rows);

		TR.store_matrix(tran);
		tran.reset();
	}

	set_orientation(matrix_flags::row_oriented);
	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



lidia_size_t * matrix< bigint >::
hnf_cg3_mod(const matrix< bigint > &B, trans_matrix &TR, int strat,
            bool &do_mod)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;
	COSMM< bigint > modul4;

	havas_kernel< bigint, COSMM< bigint >,
		COSMM< bigint >,
		havas_best_remainder_ext< bigint, COSMM< bigint >,
		COSMM< bigint > > > hnf_smodul2e;

	havas_kernel< long, COSMM< long >,
		RODMM< bigint >,
		nf_conf3e< long, COSMM< long >,
		RODMM< bigint > > > hnf_smodul3le_d;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul3;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul3;

	havas_kernel< bigint, RODMM< bigint >,
		RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >,
		RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul3;

	normalization_kernel< bigint, RODMM< bigint >, RODMM< bigint > > normalize_dmodul3;


	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;


	//
	// variables
	//

	matrix< bigint > tran;

	timer t;
	lidia_qo_info_handler(t.start_timer(););


	lidia_size_t  actual_row, actual_column;
	lidia_size_t rows_orig = actual_row = B.rows;
	lidia_size_t columns_orig = actual_column = B.columns;
	lidia_size_t *RET = NULL;
	bigint DET;
	register lidia_size_t i, j;

	bigint BOUND2;
	matrix< bigint > A;

	square(BOUND2, bigint(2147483647));
	do_mod = false;


	// Create a copy of the original matrix
	A.set_representation(matrix_flags::sparse_representation);
	A = B;

	// change orientation
	A.set_orientation(matrix_flags::column_oriented);

	bigint H = A.hadamard();

	tran.set_storage_mode(matrix_flags::sparse_representation);
	tran.set_orientation(matrix_flags::column_oriented);
	tran.resize(B.columns, B.columns);
	tran.diag(1, 0);

	// STF computation
	if (!hnf_smodul2e.stf(A, tran, actual_row, actual_column, BOUND2)) {
		TR.store_matrix(tran);

		// copying matrix
		set_representation(matrix_flags::sparse_representation);
		kill();
		resize(rows_orig, columns_orig);

		// copying matrix
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < columns_orig; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < A.value_counter[i]; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}


		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
			allocated[i] = 0;
		}

		lidia_qo_info_handler(std::cout << "\nhnfmod (row " << actual_row;
				      std::cout << ")..." << std::endl;);

		A3.set_orientation(matrix_flags::row_oriented);

		switch (strat) {
		case 0:  // dkt on regular square, havas to eliminate rest
			if (A3.rows <= 800)
				RET = A3.hnf_jacobs0(TR, H, 5);
			else
				RET = A3.hnf_jacobs1(TR, H, 5, DET, t);
			break;

		default:  // hnfmod_dkt, transformation computed with adjoint
			RET = A3.hnf_jacobs1(TR, H, 5, DET, t);
			break;
		}

		if (RET)
			do_mod = true;

		lidia_qo_info_handler(t.stop_timer();
				      std::cout << std::endl << "hnfmod: " << t << std::endl;
				      t.cont_timer(););


		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}

			value_counter[i] = k;
		}

		A3.reset();
	}
	else {
		TR.store_matrix(tran);

		// copying matrix
		set_representation(matrix_flags::sparse_representation);
		kill();
		resize(rows_orig, columns_orig);

		for (i = 0; i < columns_orig; ++i)
			for (j = 0; j < A.value_counter[i]; ++j)
				sto(A.index[i][j], i, bigint(A.value[i][j]));
	}

	set_orientation(matrix_flags::column_oriented);

	// reduce
	lidia_qo_info_handler(std::cout << "reducing..." << std::endl;);

	if (!RET) {
		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::column_oriented);
		tran.resize(B.columns, B.columns);
		tran.diag(1, 0);

		RET = normalize_smodul3.normalize_ChouCollins(*this, TR, tran, 0, B.columns - B.rows);

		TR.store_matrix(tran);
		tran.reset();
	}

	set_orientation(matrix_flags::row_oriented);

	lidia_qo_info_handler(t.stop_timer();
			      std::cout << std::endl << "total: " << t << std::endl;);

	return RET;
}



//
// solve_hnf
//
// Task:
//      solves the linear system, where A is in HNF.
//

bool matrix< bigint >::
solve_hnf(math_vector< bigint > & b, math_vector< bigint > & x)
{
	debug_handler("bigint_matrix", "my_solve");

	register lidia_size_t n, i, j;
	bigint temp, temp2, q, r;
	bool is_sol;

	is_sol = true;
	n = get_no_of_columns();
	x.set_mode(EXPAND);
	x.set_size(n);
	for (i = n-1; i >= 0; --i) {
		if (member(i, i).is_one())
			x[i] = b[i];
		else {
			temp.assign_zero();
			for (j = i+1; j < n; ++j) {
				LiDIA::multiply(temp2, member(i, j), x.member(j));
				LiDIA::add(temp, temp, temp2);
			}
			LiDIA::subtract(temp, b.member(i), temp);
			LiDIA::div_rem(q, r, temp, member(i, i));
			x[i].assign(q);
			if (!r.is_zero())
				is_sol = false;
		}

		if (!is_sol)  break;
	}

	return is_sol;
}



#undef DRMKex
#undef SRMKex

#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
