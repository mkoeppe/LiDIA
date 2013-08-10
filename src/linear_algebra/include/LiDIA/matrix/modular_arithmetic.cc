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


#ifndef LIDIA_MODULAR_ARITHMETIC_CC_GUARD_
#define LIDIA_MODULAR_ARITHMETIC_CC_GUARD_


#ifndef LIDIA_INFO_H_GUARD_
# include	"LiDIA/info.h"
#endif
#include        "LiDIA/matrix/crt_and_prime_handling.h"


#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// debug defines / error defines
//


extern const char *matrix_error_msg[];

#define DMESSAGE "modular_arithmetic"  // Debug message
#define EMESSAGE matrix_error_msg      // Error message



//
// Chinese remaindering theorem for matrices
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::chinrest (matrix< bigint > &RES,
								const matrix< bigint > *v,
								const bigint *prim) const
{
	register lidia_size_t i, j, l;

	long len;
	prim[0].longify(len);
	bigint M, X, mod;
	bigint TMP, TMP0, TMP1, TMP2;
	bigint dummy;

	register bigint *e = new bigint[len];
	memory_handler(e, DMESSAGE, "chinrest :: "
		       "Error in memory allocation (e)");

	lidia_size_t r = v[0].rows;
	lidia_size_t c = v[0].columns;

	// new dimensions
	if (RES.rows != r)
		this->rep_modul.set_no_of_rows(RES, r);
	if (RES.columns != c)
		this->rep_modul.set_no_of_columns(RES, c);

	register bigint *m = new bigint[len];
	memory_handler(m, DMESSAGE, "chinrest :: "
		       "Error in memory allocation (m)");

	// precomputation
	bigint L = 1;
	for (i = 0; i < len; i++)
		LiDIA::multiply(L, L, prim[i + 1]);

	for (i = 0; i < len; i++) {
		mod.assign(prim[i + 1]);
		LiDIA::divide(M, L, mod);
		best_remainder(TMP, M, mod);
		dummy = xgcd_right(TMP1, mod, TMP);
		LiDIA::multiply(e[i], TMP1, M);
	}


	// crt for all elements
	for (i = 0; i < r; i++)
		for (j = 0; j < c; j++) {
			X.assign_zero();

			for (l = 0; l < len; l++)
				m[l].assign(this->rep_modul.member(v[l], i, j));

			for (l = 0; l < len; l++) {
				LiDIA::multiply(TMP, e[l], m[l]);
				LiDIA::add(X, X, TMP);
			}

			best_remainder(X, X, L);
			LiDIA::subtract(TMP0, L, bigint(1));
			shift_left(TMP1, X, 1);
			shift_left(TMP2, L, 1);
			if (!(TMP1 >= -TMP0 && TMP1 <= TMP0))
				if (TMP1 - TMP2 <= TMP0 || TMP1 - TMP2 >= -TMP0)
					LiDIA::subtract(X, X, L);
				else
					lidia_error_handler("void matrix< bigint >::"
							    "chinrest(const matrix< bigint > * v, "
							    "const bigint * prim)", DMESSAGE, EMESSAGE[8]);

			this->rep_modul.sto(RES, i, j, X);
		}
	delete[] e;
	delete[] m;
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
inline void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::chinese_remainder (matrix< bigint > &RES,
									 const bigint & m1,
									 const matrix< bigint > &b,
									 const bigint & m2) const
{
	bigint mm1(abs(m1)), mm2(abs(m2)), u, d, t;
	lidia_size_t i, j;

	d = xgcd_left (u, mm1, mm2); // determine lcm(m1, m2)

	LiDIA::divide(mm2, mm2, d);

	for (i = 0; i < RES.rows; ++i)
		for (j = 0; j < RES.columns; ++j) {
			best_remainder(t, (this->rep_modul.member(b, i, j) - this->rep_modul.member(RES, i, j))/d*u, mm2);
			best_remainder(t, this->rep_modul.member(RES, i, j) + t*mm1, mm1*mm2);

			this->rep_modul.sto(RES, i, j, t);
		}
}



//
// rank
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
lidia_size_t
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::rank (const matrix< bigint > &M,
							    const bigint &H) const
{
	// read primelist from file
	long Number_of_primes;
	register bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	PRIM[0].longify(Number_of_primes);

	// Compute the rank for all primes
	lidia_size_t RANG = 0;

	lidia_size_t i, No;
	bigint Gmod;
	for (i = 1; i <= Number_of_primes; i++) {
		Gmod = PRIM[i];
		if (Gmod.bit_length() > bigint::bits_per_digit()) {
			// bigint part
			matrix< bigint > A(M.rows, M.columns, M.bitfield);
			remainder(A, M, Gmod);

			No = this->bigint_modul.rank(A, Gmod);
		}
		else {
			// long part
			long mod;
			Gmod.longify(mod);

			matrix< long > A(M.rows, M.columns, M.bitfield);
			A.set_zero_element(0);
			remainder(A, M, mod);

			No = this->long_modul.rank(A, mod);
		}

		if (RANG < No)
			RANG = No;

		// shortcut full rank
		if (RANG == M.rows || RANG == M.columns)
			return RANG;
	}
	return RANG;
}



//
// rank and linearly independent rows
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
lidia_size_t *
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::lininr1 (const matrix< bigint > &M,
							       const bigint &H) const
{
	register lidia_size_t i, j;
	bigint Gmod;

	// read primelist from file
	long Number_of_primes;
	register bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	PRIM[0].longify(Number_of_primes);

	// compute the rank and the linear independent rows for all primes
	lidia_size_t *res_vector = new lidia_size_t[M.rows + 1];
	memory_handler(res_vector, DMESSAGE, "lininr :: "
		       "Error in memory allocation (res_vector)");

	for (i = 0; i <= M.rows; res_vector[i] = 0, i++);

	lidia_size_t *v = NULL;
	for (i = 1; i <= Number_of_primes; i++) {
		Gmod = PRIM[i];
		if (Gmod.bit_length() > bigint::bits_per_digit()) {
			// bigint part
			matrix< bigint > A(M.rows, M.columns, M.bitfield);
			remainder(A, M, Gmod);

			v = this->bigint_modul.lininr(A, Gmod);
		}
		else {
			// long part
			long mod;
			Gmod.longify(mod);

			matrix< long > A(M.rows, M.columns, M.bitfield);
			A.set_zero_element(0);
			remainder(A, M, mod);

			v = this->long_modul.lininr(A, mod);
		}

		for (j = 0; j <= v[0]; j++)
			if (res_vector[j] < v[j])
				res_vector[j] = v[j];
		delete[] v;

		// full dimension
		if (res_vector[0] == M.rows || res_vector[0] == M.columns)
			return res_vector;
	}
	return res_vector;
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
lidia_size_t *
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::lininr2 (const matrix< bigint > &M,
							       const bigint &H) const
{
	register lidia_size_t i;
	bigint Gmod;

	// read primelist form file
	bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	long Number_of_primes;
	PRIM[0].longify(Number_of_primes);

	// init
	lidia_size_t *res_vector = NULL;
	lidia_size_t *v = NULL, RANG = 0;

	for (i = 1; i <= Number_of_primes; i++) {
		Gmod = PRIM[i];
		if (Gmod.bit_length() > bigint::bits_per_digit()) {
			// bigint part
			matrix< bigint > A(M.rows, M.columns, M.bitfield);
			remainder(A, M, Gmod);

			v = this->bigint_modul.lininr(A, Gmod);
		}
		else {
			// long part
			long mod;
			Gmod.longify(mod);

			matrix< long > A(M.rows, M.columns, M.bitfield);
			A.set_zero_element(0);
			remainder(A, M, mod);

			v = this->long_modul.lininr(A, mod);
		}

		if (RANG < v[0]) {
			RANG = v[0];
			if (res_vector != NULL)
				delete[] res_vector;
			res_vector = v;
		}

		// shortcut full rank
		if (RANG == M.rows || RANG == M.columns)
			return res_vector;
	}
	return res_vector;
}



//
// rank linearly independent columns
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
lidia_size_t *
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::lininc1 (const matrix< bigint > &M,
							       const bigint &H) const
{
	const bigint_matrix_algorithms< REP, REP, REP > modul;

	register lidia_size_t i, j;
	bigint Gmod;

	// read primlist from file
	register bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	long Number_of_primes;
	PRIM[0].longify(Number_of_primes);

	//  compute the rank and the linear independent columns for all primes
	lidia_size_t *res_vector = new lidia_size_t[M.columns + 1];
	memory_handler(res_vector, DMESSAGE, "lininc :: "
		       "Error in memory allocation (res_vector)");

	for (i = 0; i <= M.columns; res_vector[i] = 0, i++);

	lidia_size_t *v = NULL;

	for (i = 1; i <= Number_of_primes; i++) {
		Gmod = PRIM[i];
		if (Gmod.bit_length() > bigint::bits_per_digit()) {
			// bigint part
			matrix< bigint > A(M.columns, M.rows, M.bitfield);
			modul.trans_remainder(A, M, Gmod);

			v = this->bigint_modul.lininc(A, Gmod);
		}
		else {
			// long part
			long mod;
			Gmod.longify(mod);

			matrix< long > A(M.columns, M.rows, M.bitfield);
			A.set_zero_element(0);
			modul.trans_remainder(A, M, mod);

			v = this->long_modul.lininc(A, mod);
		}

		for (j = 0; j <= v[0]; j++)
			if (res_vector[j] < v[j])
				res_vector[j] = v[j];
		delete[] v;

		// full dimension
		if (res_vector[0] == M.rows || res_vector[0] == M.columns)
			return res_vector;
	}
	return res_vector;
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
lidia_size_t *
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::lininc2 (const matrix< bigint > &M,
							       const bigint &H) const
{
	const bigint_matrix_algorithms< REP, REP, REP > modul;

	register lidia_size_t i;
	bigint Gmod;

	// read primelist from file
	register bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	long Number_of_primes;
	PRIM[0].longify(Number_of_primes);

	// Init
	lidia_size_t *res_vector = NULL;
	lidia_size_t *v = NULL, RANG = 0;

	for (i = 1; i <= Number_of_primes; i++) {
		Gmod = PRIM[i];
		if (Gmod.bit_length() > bigint::bits_per_digit()) {
			// bigint part
			matrix< bigint > A(M.columns, M.rows, M.bitfield);
			modul.trans_remainder(A, M, Gmod);

			v = this->bigint_modul.lininc(A, Gmod);
		}
		else {
			// long part
			long mod;
			Gmod.longify(mod);

			matrix< long > A(M.columns, M.rows, M.bitfield);
			A.set_zero_element(0);
			modul.trans_remainder(A, M, mod);

			v = this->long_modul.lininc(A, mod);
		}

		if (RANG < v[0]) {
			RANG = v[0];
			if (res_vector != NULL)
				delete[] res_vector;
			res_vector = v;
		}

		// full dimension
		if (RANG == M.rows || RANG == M.columns)
			return res_vector;
	}
	return res_vector;
}



//
// adjoint matrix
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::adj1 (matrix< bigint > &RES,
							    const matrix< bigint > & A,
							    const bigint &H,
							    const bigint &DET) const
{
	register lidia_size_t i, z1, z2;

	// read primelist from file
	register bigint *PRIM = get_primes(bigint(2)*H, DET, true);
	long n;
	PRIM[0].longify(n);

	// compute adj for all primes
	bigint MOD;
	long Modlong;
	matrix< bigint > *chininput = new matrix< bigint > [n];
	memory_handler(chininput, DMESSAGE, "adj :: "
		       "Error in memory allocation (chininput)");

	for (i = 1; i <= n; i++) {
		debug_handler_c(DMESSAGE, "adj(const matrix< bigint > &)", DVALUE + 7,
				std::cout << "In Iteration " << i << std::endl;);

		MOD.assign(PRIM[i]);

		if (MOD.bit_length() > bigint::bits_per_digit()) {
			matrix< bigint > B(A.rows, A.columns, A.bitfield);
			remainder(B, A, MOD);

			this->bigint_modul.adj(B, MOD);

			chininput[i - 1] = B;
		}
		else {
			MOD.longify(Modlong);
			matrix< long > B(A.rows, A.columns, A.bitfield);
			remainder(B, A, Modlong);

			this->long_modul.adj(B, Modlong);

			chininput[i - 1].set_no_of_rows(A.rows);
			chininput[i - 1].set_no_of_columns(A.columns);

			for (z1 = 0; z1 < A.rows; z1++)
				for (z2 = 0; z2 < A.columns; z2++)
					this->rep_modul.sto(chininput[i-1], z1, z2, B(z1, z2));
		}
	}

	// CRT
	this->chinrest(RES, chininput, PRIM);

	delete[] PRIM;
	delete[] chininput;
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::adj2 (matrix< bigint > &RES,
							    const matrix< bigint > & A,
							    const bigint &H,
							    const bigint & DET) const
{
	register lidia_size_t i, ii, jj, kk, z1, z2;
	bigint val, val2;
	bool done;

	register bigint *PRIM = get_primes(bigint(2)*H, abs(DET), true);
	long n;
	PRIM[0].longify(n);

	lidia_qo_info_handler(std::cout << "(adj) " << n << " primes required!" << std::endl;);

	// Step 4
	bigint MODULUS, MOD;
	long Modlong;
	memory_handler(chininput, DMESSAGE, "adj :: "
		       "Error in memory allocation (chininput)");

	matrix< bigint > B(A.rows, A.columns, A.bitfield);

	for (i = 1; i <= n; i++) {
		debug_handler_c(DMESSAGE, "in function "
				"adj(const matrix< bigint > &)", DVALUE + 7,
				std::cout << "In Iteration " << i << std::endl;);

		MOD.assign(PRIM[i]);

		if (MOD.bit_length() > bigint::bits_per_digit()) {
			lidia_qo_xinfo_handler("adjoint [bigint]", i, n);

			remainder(B, A, MOD);

			this->bigint_modul.adj(B, MOD);
		}
		else {
			lidia_qo_xinfo_handler("adjoint [long]", i, n);

			MOD.longify(Modlong);
			matrix< long > Blong(A.rows, A.columns, A.bitfield);
			remainder(Blong, A, Modlong);

			this->long_modul.adj(Blong, Modlong);

			for (z1 = 0; z1 < A.rows; z1++)
				for (z2 = 0; z2 < A.columns; z2++)
					this->rep_modul.sto(B, z1, z2, Blong(z1, z2));
		}

		if (i == 1) {
			RES = B;
			MODULUS = MOD;
		}
		else {
			this->chinese_remainder(RES, MODULUS, B, MOD);
			MODULUS *= MOD;
		}

		// perform test (A*RES = det I)
		done = true;
		for (ii = 0; (ii < A.rows) && done; ++ii)
			for (jj = 0; (jj < A.columns) && done; ++jj) {
				val.assign_zero();
				for (kk = 0; kk < A.columns; ++kk)
					add(val, val, A.member(ii, kk)* RES.member(kk, jj));
				if (ii == jj)
					done = (val == DET);
				else
					done = val.is_zero();
			}

		if (done) {
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



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::adj2 (matrix< bigint > &RES,
							    const matrix< bigint > & A,
							    const bigint &H,
							    const bigint & DET,
							    int num_mod) const
{
	register lidia_size_t i, z1, z2, ii, kk;
	bigint val, val2;

	register bigint *PRIM = get_primes(bigint(2)*H, abs(DET), true);
	long n;
	PRIM[0].longify(n);

	lidia_qo_info_handler(std::cout << "(adj) " << num_mod << " primes required!";
			      std::cout << std::endl;);

	// Step 4
	bigint MODULUS, MOD;
	long Modlong;
	memory_handler(chininput, DMESSAGE, "adj :: "
		       "Error in memory allocation (chininput)");

	matrix< bigint > B(A.rows, A.columns, A.bitfield);

	bool done = false;
	for (i = 1; (i <= num_mod) || !done; i++) {
		debug_handler_c(DMESSAGE, "in function "
				"adj(const matrix< bigint > &)", DVALUE + 7,
				std::cout << "In Iteration " << i << std::endl;);

		MOD.assign(PRIM[i]);

		if (MOD.bit_length() > bigint::bits_per_digit()) {
			lidia_qo_xinfo_handler("adjoint [bigint]", i, n);

			remainder(B, A, MOD);

			this->bigint_modul.adj(B, MOD);
		}
		else {
			lidia_qo_xinfo_handler("adjoint [long]", i, n);

			MOD.longify(Modlong);
			matrix< long > Blong(A.rows, A.columns, A.bitfield);
			remainder(Blong, A, Modlong);

			this->long_modul.adj(Blong, Modlong);

			for (z1 = 0; z1 < A.rows; z1++)
				for (z2 = 0; z2 < A.columns; z2++)
					this->rep_modul.sto(B, z1, z2, Blong(z1, z2));
		}

		if (i == 1) {
			RES = B;
			MODULUS = MOD;
		}
		else {
			this->chinese_remainder(RES, MODULUS, B, MOD);
			MODULUS *= MOD;
		}


		if (i >= num_mod) {
			// perform test (A*RES = det I), only check diagonals
			done = true;
			for (ii = 0; (ii < A.rows) && done; ++ii) {
				val.assign_zero();
				for (kk = 0; kk < A.columns; ++kk)
					add(val, val, A.member(ii, kk)* RES.member(kk, ii));
				done = (val == DET);

				if (done && (ii > 0)) {
					val.assign_zero();
					for (kk = 0; kk < A.columns; ++kk)
						add(val, val, A.member(ii, kk)* RES.member(kk, 0));
					done = (val.is_zero());
				}
			}
		}

		if ((i >= num_mod) && done) {
			lidia_qo_info_handler(std::cout << "\n(adj) Only " << i << " primes required!" << std::endl;);
		}
	}

	if (qo_special) {
		std::cout << " " << i-1 << std::flush;
	}

	delete[] PRIM;
}



//
// lattice determinant
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet1 (const matrix< bigint > &RES,
								   bigint & DET,
								   const bigint &H) const
{
	register lidia_size_t i, j;

	//  Step 1
	lidia_size_t *linuz = lininr1(RES, H);
	lidia_size_t r = linuz[0];

	// Step 2
	matrix< bigint > C(r, RES.columns, RES.bitfield);
	for (i = 0; i < r; i++)
		for (j = 0; j < RES.columns; j++)
			this->rep_modul.sto(C, i, j, this->rep_modul.member(RES, linuz[i+1], j));

	delete[] linuz;

	// Step 3
	if (r == RES.columns)
		det(C, DET, H);
	else {
		linuz = lininc1(C, H);

		matrix< bigint > D(r, r, RES.bitfield);
		for (i = 0; i < r; i++)
			for (j = 0; j < r; j++)
				this->rep_modul.sto(D, j, i, this->rep_modul.member(C, j, linuz[i + 1]));
		det(D, DET, H);
		delete[] linuz;
	}

	if (DET.is_lt_zero())
		DET.negate();
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet2 (const matrix< bigint > &RES,
								   bigint & DET,
								   const bigint &H) const
{
	register lidia_size_t i, j;
	bigint TMP, TMP1;

	// Step 1
	lidia_size_t *linu = lininr1(RES, H);
	lidia_size_t r = linu[0];

	// Step 2
	matrix< bigint > C(r, RES.columns, RES.bitfield);
	matrix< bigint > C1(r, RES.columns, RES.bitfield);

	for (i = 0; i < r; i++)
		for (j = 0; j < RES.columns; j++) {
			this->rep_modul.sto(C, i, j, this->rep_modul.member(RES, linu[i+1], j));
			this->rep_modul.sto(C1, i, RES.columns - 1 - j, this->rep_modul.member(RES, linu[i+1], j));
		}
	delete[] linu;

	// Step 3
	if (r == RES.columns)
		C.det(DET, H);
	else {
		linu = lininc1(C, H);
		matrix< bigint > D(r, r, RES.bitfield);

		for (i = 0; i < r; i++)
			for (j = 0; j < r; j++)
				this->rep_modul.sto(D, j, i, this->rep_modul.member(C, j, linu[i + 1]));
		D.det(TMP, H);
		delete[] linu;

		linu = lininc1(C1, H);
		for (i = 0; i < r; i++)
			for (j = 0; j < r; j++)
				this->rep_modul.sto(D, j, i, this->rep_modul.member(C1, j, linu[i + 1]));
		D.det(TMP1, H);
		delete[] linu;
		DET = gcd(TMP, TMP1);
	}
	if (DET.is_lt_zero())
		DET.negate();
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet2 (const matrix< bigint > &RES,
								   bigint & DET,
								   bigint &H,
								   int num_same) const
{
	register lidia_size_t i, j;
	bigint TMP, TMP1, *tmp, *tmp1, *tmp2;

	// Step 1
	lidia_size_t *linu = lininr1(RES, H);
	lidia_size_t r = linu[0];

	// Step 2
	matrix< bigint > C(r, RES.columns);
	matrix< bigint > C1(r, RES.columns);

	for (i = 0; i < r; i++) {
		tmp = C.value[i];
		tmp1 = RES.value[linu[i+1]];
		tmp2 = C1.value[i];
		for (j = 0; j < RES.columns; j++) {
			tmp[j].assign(tmp1[j]);
			tmp2[RES.columns - 1 - j].assign(tmp1[j]);
		}
	}
	delete[] linu;

	// Step 3
	if (r == RES.columns)
		C.det(DET, H, num_same);
	else {
		linu = lininc1(C, H);
		matrix< bigint > D(r, r);

		for (i = 0; i < r; i++)
			for (j = 0; j < r; j++)
				D.value[j][i].assign(C.value[j][linu[i + 1]]);
		D.det(TMP, H, num_same);
		delete[] linu;

		linu = lininc1(C1, H);
		for (i = 0; i < r; i++)
			for (j = 0; j < r; j++)
				D.value[j][i].assign(C1.value[j][linu[i + 1]]);
		D.det(TMP1, H, num_same);
		delete[] linu;

		DET = gcd(TMP, TMP1);
	}
	if (DET.is_lt_zero())
		DET.negate();
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet3 (const matrix< bigint > &RES,
								   bigint & DET,
								   const bigint &H) const
{
	register lidia_size_t i, j;
	register bigint *tmp, *tmp1;

	// Step 1
	lidia_size_t *linuz = lininr1(RES, H);
	lidia_size_t r = linuz[0];

	// Step 2
	if (r == RES.columns) {
		matrix< bigint > C(r, RES.columns);
		for (i = 0; i < r; i++) {
			tmp = C.value[i];
			tmp1 = RES.value[linuz[i+1]];
			for (j = 0; j < RES.columns; j++)
				tmp[j].assign(tmp1[j]);
		}
		C.det(DET, H);
	}
	else {
		lidia_size_t *linuz1 = lininc1(RES, H);

		matrix< bigint > D(r, r);
		for (i = 0; i < r; i++)
			for (j = 0; j < r; j++)
				D.value[j][i].assign(RES.value[linuz[j+1]][linuz1[i + 1]]);
		D.det(DET, H);
		delete[] linuz1;
	}
	delete[] linuz;
	if (DET.is_lt_zero())
		DET.negate();
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet4 (const matrix< bigint > &RES,
								   bigint & DET,
								   bigint &H,
								   int num_same) const
{
	register lidia_size_t i, j;
	bigint TMP, TMP1, *tmp, *tmp1, *tmp2;

	// Step 1
	lidia_size_t *linu2, *linu = lininr2(RES, H);
	lidia_size_t r = linu[0];

	matrix< bigint > C(r, RES.columns);
	for (i = 0; i < r; i++) {
		tmp = C.value[i];
		tmp1 = RES.value[linu[i+1]];
		for (j = 0; j < RES.columns; j++)
			tmp[j].assign(tmp1[j]);
	}

	linu2 = lininc1(C, H);
	for (i = 0; i < r; i++)
		for (j = 0; j < r; j++)
			C.value[j][i].assign(C.value[j][linu2[r-i]]);
	C.resize(r, r);
	delete[] linu2;
	linu2 = NULL;
	det(C, TMP, H, num_same);


	// second submatrix
	C.resize(r, RES.columns);
	for (i = 0; i < r; i++) {
		tmp1 = RES.value[linu[i+1]];
		tmp2 = C.value[i];
		for (j = 0; j < RES.columns; j++)
			tmp2[RES.columns - 1 - j].assign(tmp1[j]);
	}
	delete[] linu;

	linu2 = lininc1(C, H);
	for (i = 0; i < r; i++)
		for (j = 0; j < r; j++)
			C.value[j][i].assign(C.value[j][linu2[r-i]]);
	C.resize(r, r);
	delete[] linu2;
	det(C, TMP1, H, num_same);
	DET = gcd(TMP, TMP1);

	if (DET.is_lt_zero())
		DET.negate();
	lidia_info_handler(std::cout << "Det mult = " << DET << std::endl;);
}



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
lidia_size_t *
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet5 (matrix< bigint > &RES,
								   bigint & DET,
								   bigint &H,
								   int num_same) const
{
	register lidia_size_t i;
	bigint TMP, TMP1;

	lidia_size_t *linu1 = lininr1(RES, H);
	if (linu1[0] < RES.rows)
		return linu1;

	if (RES.rows == RES.columns) {
		det(RES, DET, H, num_same);
		if (qo_special) {
			std::cout << " 0" << std::flush;
		}
		lidia_qo_info_handler(std::cout << "Det mult = " << DET << std::endl;);

		return NULL;
	}

	for (i = 0; i < linu1[0]; ++i)
		RES.swap_rows(i, linu1[linu1[0] - i]);

	// Step 1
	lidia_size_t *linu = lininc1(RES, H);
	lidia_size_t r = linu[0];
	for (i = 0; i < r; i++)
		RES.swap_columns(i, linu[r-i]);

	lidia_size_t old_columns = RES.columns;
	RES.columns = linu[0];
	lidia_size_t old_rows = RES.rows;
	RES.rows = linu1[0];

	det(RES, TMP, H, num_same);

	RES.columns = old_columns;
	RES.rows = old_rows;
	for (i = r-1; i >= 0; i--)
		RES.swap_columns(i, linu[r-i]);
	for (i = 0; i < RES.columns/2; i++)
		RES.swap_columns(i, RES.columns-i-1);
	delete [] linu;

	linu = lininc1(RES, H);
	for (i = 0; i < r; i++)
		RES.swap_columns(i, linu[r-i]);
	RES.columns = linu[0];
	RES.rows = linu1[0];

	det(RES, TMP1, H, num_same);

	RES.rows = old_rows;
	RES.columns = old_columns;
	for (i = r-1; i >= 0; i--)
		RES.swap_columns(i, linu[r-i]);
	for (i = 0; i < RES.columns/2; i++)
		RES.swap_columns(i, RES.columns-i-1);

	delete[] linu;
	DET = gcd(TMP, TMP1);

	if (DET.is_lt_zero())
		DET.negate();
	lidia_qo_info_handler(std::cout << "Det mult = " << DET << std::endl;);

	for (i = linu1[0]-1; i >= 0; --i)
		RES.swap_rows(i, linu1[linu1[0]-i]);
	delete [] linu1;

	linu1 = NULL;
	return linu1;
}



//
// lattice determinant special
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::latticedet_special (const matrix< bigint > &C,
									  bigint & DET,
									  const bigint &H,
									  int num_same) const
{
	matrix< bigint > B = C;

	// Step 1, 2
	lidia_size_t num = 0;
	bigint *PRIM = get_primes(bigint(2) * H, bigint(1));

	lidia_size_t *linu = lininc1(B, H);
	lidia_size_t r = linu[0];
	lidia_size_t i;
	for (i = 0; i < r; i++)
		B.swap_columns(B.columns - B.rows + i, linu[r-i]);
	delete[] linu;

	long n, Modlong;
	PRIM[0].longify(n);
	lidia_info_handler(std::cout << n << " primes required for modulo lattice determinant computation.\n" << std::flush;);

	bigint MODULUS = 1;

	lidia_size_t diff = B.columns - B.rows;
	bigint *RES = new bigint[diff + 1];

	const bigint *bigintRES = NULL;
	bigint Gmod, DETOLD;
	const long *longRES = NULL;

	matrix< bigint > Abigint(B.rows, B.columns, B.bitfield);
	matrix< long > Along(B.rows, B.columns, B.bitfield);
	Along.set_zero_element(0);

	DET = 0;

	// Step 3
	lidia_size_t j;
	for (i = 1; i <= n; i++) {
		Gmod.assign(PRIM[i]);
		if (Gmod.bit_length() >= bigint::bits_per_digit()) {
			// bigint part
			lidia_qo_xinfo_handler("lattice determinante [bigint]", i, n);

			remainder(Abigint, B, Gmod);

			bigintRES = this->bigint_modul.STF_extended(Abigint, Gmod);

			if (i == 1) {
				for (j = 0; j <= diff; j++)
					RES[j] = bigintRES[j];
				MODULUS = Gmod;
			}
			else {
				for (j = 0; j <= diff; j++)
					RES[j] = LiDIA::chinese_remainder(RES[j], MODULUS, bigintRES[j], Gmod);
				MODULUS *= Gmod;
				for (j = 0; j <= diff; j++)
					best_remainder(RES[j], RES[j], MODULUS);

			}
			delete[] bigintRES;
		}
		else {
			// long part
			lidia_qo_xinfo_handler("lattice determinante [long]", i, n);

			Gmod.longify(Modlong);
			remainder(Along, B, Modlong);

			longRES = this->long_modul.STF_extended(Along, Modlong);
			Gmod = Modlong;

			if (i == 1) {
				for (j = 0; j <= diff; j++)
					RES[j] = longRES[j];
				MODULUS = Gmod;
			}
			else {
				for (j = 0; j <= diff; j++)
					RES[j] = LiDIA::chinese_remainder(RES[j], MODULUS, longRES[j], Gmod);
				MODULUS *= Gmod;
				for (j = 0; j <= diff; j++)
					best_remainder(RES[j], RES[j], MODULUS);
			}
			delete[] longRES;
		}

		DET = RES[diff];
		if (DETOLD == DET) {
			num++;
			if (num == num_same) {
				lidia_info_handler(std::cout << "\nOnly " << i << " primes required!\n" << std::flush;);
				break;
			}
		}
		else {
			num = 0;
			DETOLD = DET;
		}
	}
	DET = RES[0];
	for (j = 1; j <= diff; j++)
		DET = gcd(DET, RES[j]);
	delete[] RES;
}



//
// determinant
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::det (const matrix< bigint > &B,
							   bigint & DET,
							   const bigint &H) const
{
	// read primes from file
	bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	long n, Modlong;
	PRIM[0].longify(n);

	bigint Gmod;

	// init
	bigint *chininput = new bigint[n];
	memory_handler(chininput, DMESSAGE, "det :: "
		       "Error in memory allocation (chininput)");

	matrix< bigint > Abigint(B.rows, B.columns, B.bitfield);
	matrix< long > Along(B.rows, B.columns, B.bitfield);
	Along.set_zero_element(0);

	for (lidia_size_t i = 1; i <= n; i++) {
		Gmod.assign(PRIM[i]);
		if (Gmod.bit_length() >= bigint::bits_per_digit()) {
			// bigint part
			lidia_qo_xinfo_handler("determinante [bigint]", i, n);

			remainder(Abigint, B, Gmod);

			chininput[i - 1] = this->bigint_modul.det(Abigint, Gmod);
		}
		else {
			// long part
			lidia_qo_xinfo_handler("determinante [long]", i, n);

			Gmod.longify(Modlong);
			remainder(Along, B, Modlong);

			chininput[i - 1] = this->long_modul.det(Along, Modlong);
		}
		best_remainder(chininput[i - 1], chininput[i - 1], Gmod);
	}

	// CRT
	LiDIA::chinrest(DET, (const bigint *)chininput, (const bigint *)PRIM);
	delete[] chininput;
}



//
// determinant (with iterative chinese remaindering)
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
int
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::det (const matrix< bigint > &B,
							   bigint & DET,
							   const bigint &H,
							   int num_same) const
{
	// read primes from file
	lidia_size_t i, num = 0;
	bigint *PRIM = get_primes(bigint(2) * H, bigint(1));
	long n, Modlong;
	PRIM[0].longify(n);
	lidia_qo_info_handler(std::cout << n << " primes required for modulo determinant computation.\n"
			      << std::flush;);

	bigint MODULUS = 1;
	bigint RES, Gmod, DETOLD;

	matrix< bigint > Abigint(B.rows, B.columns, B.bitfield);
	matrix< long > Along(B.rows, B.columns, B.bitfield);
	Along.set_zero_element(0);

	// compute determinante for all finite fields
	for (i = 1; i <= n; i++) {
		Gmod.assign(PRIM[i]);
		if (Gmod.bit_length() >= bigint::bits_per_digit()) {
			// bigint part
			lidia_qo_xinfo_handler("determinante [bigint]", i, n);

			remainder(Abigint, B, Gmod);

			RES = this->bigint_modul.det(Abigint, Gmod);
		}
		else {
			// long part
			lidia_qo_xinfo_handler("determinante [long]", i, n);

			Gmod.longify(Modlong);
			remainder(Along, B, Modlong);

			RES = this->long_modul.det(Along, Modlong);
			Gmod = Modlong;
		}

		if (i == 1) {
			DET = RES;
			MODULUS = Gmod;
		}
		else {
			DET = LiDIA::chinese_remainder(DET, MODULUS, RES, Gmod);
			MODULUS *= Gmod;
			best_remainder(DET, DET, MODULUS);
		}

		if (DETOLD == DET) {
			num++;
			if (num == num_same) {
				lidia_qo_info_handler(std::cout << "\nOnly " << i << " primes required!\n" << std::flush;);
				break;
			}
		}
		else {
			num = 0;
			DETOLD = DET;
		}

	}

	if (qo_special) {
		if (i == n+1)  --i;
		std::cout << " " << i << std::flush;
	}

	return i;
}



//
// characteristic polynomial
//

template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
void
modular_arithmetic< REP, SINGLE_MODUL, MULTI_MODUL >::charpoly (const matrix< bigint > &A,
								bigint *RES,
								const bigint &H) const
{
	register lidia_size_t i, j;
	bigint TMP = H, TMP1;
	long len;

	// Computing bound
	j = static_cast<lidia_size_t>(A.columns) / 2;

	TMP1.assign_one();
	for (i = A.columns; i > A.columns - j; i--)
		LiDIA::multiply(TMP, TMP, (i*i));
	for (i = j; i > 0; i--)
		LiDIA::multiply(TMP1, TMP1, (i*i));
	LiDIA::divide(TMP, TMP, TMP1);

	if (TMP.is_zero())
		RES[A.columns].assign((A.columns % 2 == 0) ? 1 : -1);

	// read primes from file
	shift_left(TMP, TMP, 1);
	bigint *PRIM = get_primes(TMP, bigint(1));
	PRIM[0].longify(len);

	// Init
	matrix< bigint > U(A.columns + 1, static_cast<int>(len));

	bigint MOD;
	long Modlong;

	// computing charpoly
	for (i = 1; i <= len; i++) {
		MOD.assign(PRIM[i]);
		if (MOD.bit_length() > bigint::bits_per_digit()) {
			matrix< bigint > B(A.rows, A.columns);
			remainder(B, A, MOD);

			bigint *zwbigint = this->bigint_modul.charpoly(B, MOD);

			for (j = 0; j < A.columns + 1; j++)
				U.value[j][i - 1].assign(zwbigint[j]);
			delete[] zwbigint;
		}
		else {
			MOD.longify(Modlong);
			matrix< long > B(A.rows, A.columns);
			B.set_zero_element(0);
			remainder(B, A, Modlong);

			long *zwlong = this->long_modul.charpoly(B, Modlong);
			for (j = 0; j < A.columns + 1; j++)
				U.value[j][i - 1].assign(zwlong[j]);
			delete[] zwlong;
		}

	}

	// CRT
	for (i = 0; i < A.columns + 1; i++)
		LiDIA::chinrest(RES[i], U.value[i], PRIM);
}



#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_MODULAR_ARITHMETIC_CC_GUARD_
