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
#include	"LiDIA/bigint_matrix.h"
#include	"LiDIA/modular_operations.inl"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bigint * matrix< bigint >::
mgcd(const bigint * a, lidia_size_t n)
{
	//
	// Task: RES = mgcd(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	// ALGORITHM: Bradley
	// Version: 1.8
	//

	debug_handler_l("multiple_gcd", "in function "
			"mgcd(const bigint *, lidia_size_t)", LDBL_MATRIX);

	bigint *RES = NULL;
	register lidia_size_t i;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		bigint TMP, TMP1, TMP2;

		bigint *tmp = new bigint[n];
		memory_handler(tmp, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (tmp)");
		bigint *y = new bigint[n];
		memory_handler(y, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (y)");
		bigint *z = new bigint[n];
		memory_handler(z, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (z)");

		tmp[0].assign(a[0]);
		for (i = 1; i < n; i++)
			tmp[i] = xgcd(y[i], z[i], tmp[i - 1], a[i]);

		// Step 9 - 19
		bigint *x = new bigint[n];
		memory_handler(x, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (x)");
		bigint *y1 = new bigint[n];
		memory_handler(y1, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (y1)");

		x[n-1].assign(z[n-1]);
		y1[n-1].assign(y[n-1]);

		bigint *G = new bigint[n];
		memory_handler(G, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (G)");
		bigint *S = new bigint[n];
		memory_handler(S, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (S)");
		for (i = n - 2; i >= 1; i--) {
			div_rem(G[i], TMP, tmp[i - 1], tmp[i]);
			if (G[i].is_zero())
				S[i].assign_zero();
			else
				nearest(S[i], y1[i + 1] * z[i], G[i]);

			// y1[i] = y1[i+1]*y[i] + S[i]*(a[i]/tmp[i]);
			div_rem(TMP1, TMP, a[i], tmp[i]);
			LiDIA::multiply(TMP1, S[i], TMP1);
			LiDIA::multiply(TMP2, y1[i + 1], y[i]);
			LiDIA::add(y1[i], TMP2, TMP1);

			// x[i] = y1[i+1]*z[i] - S[i]*G[i];
			LiDIA::multiply(TMP1, S[i], G[i]);
			LiDIA::multiply(TMP2, y1[i + 1], z[i]);
			LiDIA::subtract(x[i], TMP2, TMP1);
		}
		x[0].assign(y1[1]);

		// Step 20,21
		RES = new bigint[n + 1];
		memory_handler(RES, "matrix_bigint", "mgcd :: "
			       "Error in memory allocation (RES)");
		RES[0].assign(tmp[n - 1]);
		for (i = 1; i <= n; i++)
			RES[i].assign(x[i - 1]);

		delete[] tmp;
		delete[] y;
		delete[] z;
		delete[] x;
		delete[] y1;
		delete[] G;
		delete[] S;
	}
	if (n == 1) {
		RES = new bigint [2];
		RES[0].assign(a[0]);
		RES[1].assign_one();
	}
	if (RES[0].is_lt_zero())
		for (i = 0; i <= n; i++)
			RES[i].negate();
	return RES;
}



bigint * matrix< bigint >::
mgcd_new(const bigint * a, lidia_size_t n)
{
	//
	// Task: RES = mgcd(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	// ALGORITHM: Bradley
	// Version: 1.8
	//

	debug_handler_l("multiple_gcd", "in function "
			"mgcd(const bigint *, lidia_size_t)", LDBL_MATRIX);

	bigint *RES = new bigint[n+1];
	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		matrix< bigint > Tr(n, n);
		Tr.diag(1, 0);

		bigint y, z, g, TMP, TMP1, t1, t2;
		bigint g_old = a[0];
		for (i = 1; i < n; i++) {
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

		// Balaciertes Rueckwaertseinsetzen
		bigint q, r;
		lidia_size_t l;
		for (i = n - 1, j = n - 2; i >= 0 && j >= 0; i--, j--) {
			div_rem(q, r, Tr.value[i][n-1], Tr.value[i][j]);
			if (abs(q) > 0)
				for (l = 0; l < n; l++)
					Tr.value[l][n-1] = Tr.value[l][n-1] - q*Tr.value[l][j];
		}
		//std::cout << Tr << std::endl;
		for (j = 0; j < n; j++) {
			TMP = 0;
			for (i = 0; i < n; i++)
				TMP += a[i] * Tr.value[i][j];
			//std::cout << TMP << " " << std::flush;
		}
		//std::cout << std::endl;

		for (j = 0; j < n; j++)
			RES[j+1] = Tr.value[j][n-1];
		RES[0] = g;
	}
	return RES;
}



bigint * matrix< bigint >::
mgcd_new2(const bigint * a, lidia_size_t n)
{
	//
	// Task: RES = mgcd(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	// ALGORITHM: Bradley
	// Version: 1.8
	//

	debug_handler_l("multiple_gcd", "in function "
			"mgcd(const bigint *, lidia_size_t)", LDBL_MATRIX);

	bigint *RES = new bigint[n+1];
	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		matrix< bigint > Tr(n, n);
		Tr.diag(1, 0);

		bigint y, z, g, TMP, TMP1, t1, t2;
		bigint g_old = a[0];
		for (i = 1; i < n; i++) {
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

		// Balaciertes Rueckwaertseinsetzen
		bigint q, r;
		lidia_size_t l, k;
		for (i = n - 1, j = n - 2; i >= 0 && j >= 0; i--, j--)
			for (k = j + 1; k < n; k++) {
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
		//std::cout << Tr << std::endl;
		for (j = 0; j < n; j++) {
			TMP = 0;
			for (i = 0; i < n; i++)
				TMP += a[i] * Tr.value[i][j];
			//std::cout << TMP << " " << std::flush;
		}
		//std::cout << std::endl;

		for (j = 0; j < n; j++)
			RES[j+1] = Tr.value[j][n-1];
		RES[0] = g;
	}
	return RES;
}



bigint * matrix< bigint >::
mgcd_new3(const bigint * b, lidia_size_t n)
{
	//
	// Task: RES = mgcd(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	// ALGORITHM: Bradley
	// Version: 1.8
	//

	debug_handler_l("multiple_gcd", "in function "
			"mgcd(const bigint *, lidia_size_t)", LDBL_MATRIX);
	bigint *a = new bigint[n];
	for (lidia_size_t h = 0; h < n; h++)
		a[h] = b[h];

	bigint *RES = new bigint[n+1];
	register lidia_size_t i, j;

	if (n <= 0)
		lidia_error_handler("multiple_gcd", "mgcd :: Error in parameter !!");

	if (n > 1) {
		// Step 1 - 8
		matrix< bigint > Tr(n, n);
		Tr.diag(1, 0);

		bigint y, z, g, TMP, TMP1, t1, t2;

		lidia_size_t step = 1, step2 = 2;
		while (step < n) {
			for (i = 0; i < n; i += step2) {
				if (i+step < n) {
					g = xgcd(y, z, a[i], a[i+step]);
					TMP = -a[i+step]/g;
					TMP1 = a[i]/g;
					for (j = 0; j < n; j++) {
						t1 = Tr.value[j][i];
						t2 = Tr.value[j][i + step];
						Tr.value[j][i] = y*t1+z*t2;
						Tr.value[j][i + step] = TMP*t1+TMP1*t2;
					}
					t1 = a[i];
					t2 = a[i+step];
					a[i] = y*t1+z*t2;
					a[i+step] = TMP*t1+TMP1*t2;
				}
			}
			step = step2;
			step2 = 2*step2;
		}
		Tr.swap_columns(0, n-1);
		//std::cout << Tr << std::endl;
		for (j = 0; j < n; j++) {
			TMP = 0;
			for (i = 0; i < n; i++)
				TMP += b[i] * Tr.value[i][j];
			//std::cout << TMP << " " << std::flush;
		}
		//std::cout << std::endl;

		for (j = 0; j < n; j++)
			RES[j+1] = Tr.value[j][n-1];
		RES[0] = g;
	}
	delete [] a;
	return RES;
}


bigint * matrix< bigint >::
mgcd_tree(const bigint * a, lidia_size_t n)
{
	//
	// Task: RES = mgcd(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	// ALGORITHM: Bradley
	// Version: 1.8
	//

	debug_handler_l("multiple_gcd", "in function "
			"mgcd(const bigint *, lidia_size_t)", LDBL_MATRIX);
	bigint *g = new bigint[2*n+1];
	bigint *y = new bigint[2*n+1];
	bigint *z = new bigint[2*n+1];
	lidia_size_t i;

	for (i = 0; i < n; i++)
		g[n+i] = a[i];

	// Bottom-Up Phase
	for (i = n-1; i > 0; i--)
		g[i] = xgcd(y[2*i], y[2*i+1], g[2*i], g[2*i+1]);

	// Top-Down Phase
	z[0] = g[1]; // ggt
	z[2] = y[2];
	z[3] = y[3];
	bigint u, u1, w, w1, k, r;
	for (i = 2; i < n; i++) {
		//g[i] = g[i]*y[i];
		u = g[2*i]/g[i]; u1 = g[2*i+1]/g[i];

		w = y[2*i]*y[i]; w1 = y[2*i+1]*y[i];

		if (u1 == 0)
			k = 0;
		else
			nearest(k, w, u1);
		y[2*i+1] = w1 + (k*u);
		y[2*i] = w - (k*u1);
	}
	y[n-1] = g[1];

	delete[] z;
	delete[] g;

	return &y[n-1];
}


//
// mgcd2
//

void
mgcd2(bigint & RES, const bigint * aconst, lidia_size_t n)
{
	//
	// Task: mgcd2(Res,a,n);
	//              => RES = gcd(a[0],...,a[n-1])
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements, MIN assignments
	// PAPER: Hermite normal form computation for integer matrices, Havas
	// Version: 1.8
	//

	debug_handler("multiple_gcd", "in member - function "
		      "mgcd2(bigint &, const bigint *, lidia_size_t)");

	register lidia_size_t i, index, SW, bound;
	bigint MIN, TMP, q, r;

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd2 :: "
		       "Error in memory allocation (a)");

	for (i = 0; i < n; i++)
		a[i].assign(aconst[i]);

	// init
	for (index = 0; index < n && a[index].is_zero(); index++);

	if (index == n) {
		RES.assign_zero();
		return;
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

		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW = 1;
				div_rem(q, r, a[i], MIN);
				a[i].assign(r);
			}
	} while (SW == 1);

	// gcd < 0 ?
	if (a[index] < 0)
		LiDIA::negate(RES, a[index]);
	else
		RES.assign(a[index]);
}



bigint *
mgcd2(const bigint * aconst, lidia_size_t n)
{
	//
	// Task: RES = mgcd2(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements, MIN assignments
	// PAPER: Hermite normal form computation for integer matrices, Havas
	// Version: 1.8
	//

	debug_handler("multiple_gcd", "in function "
		      "mgcd2(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, SW, bound;
	bigint MIN, TMP, q, r, *Ttmp1, *Ttmp2 = NULL;

	matrix< bigint > T(n, n);
	T.diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd2 :: "
		       "Error in memory allocation (a)");

	for (i = 0; i < n; i++)
		a[i].assign(aconst[i]);


	for (index = 0; index < n && a[index].is_zero(); index++);

	if (index == n) {
		delete[] a;
		return new bigint[n];
	}
	else
		bound = index;

	do {
		MIN.assign(a[index]);


		for (i = bound; i < n; i++)
			if ((abs(MIN) > abs(a[i])) && !a[i].is_zero()) {
				MIN.assign(a[i]);
				index = i;
			}


		SW = 0;

		Ttmp2 = T.value[index];
		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW = 1;
				Ttmp1 = T.value[i];
				div_rem(q, r, a[i], MIN);
				a[i].assign(r);
				for (j = 0; j < n; j++) {
					LiDIA::multiply(TMP, q, Ttmp2[j]);
					LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
				}
			}
	} while (SW == 1);

	Ttmp2 = T.value[index];


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



#if 0
bigint *matrix < bigint >::
mgcd2(const bigint * aconst, lidia_size_t n)
{
	//
	// Task: RES = T.mgcd2(a,n);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	//              => T*a = RES
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements, MIN assignments
	// PAPER: Hermite normal form computation for integer matrices, Havas
	// Version: 1.8
	//

	debug_handler("multiple_gcd", "in member - function "
		      "mgcd2(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, bound, SW;
	bigint MIN, TMP, q, r, *Ttmp1, *Ttmp2 = NULL;

	if (columns != n)
		set_no_of_columns(n);
	if (rows != n)
		set_no_of_rows(n);
	diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd2 :: "
		       "Error in memory allocation (a)");

	for (i = 0; i < n; i++)
		a[i].assign(aconst[i]);

	// init
	for (index=0; index<n && a[index].is_zero();index++);

	if (index==n) {
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
		SW=0;
		Ttmp2 = value[index];
		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW=1;
				Ttmp1 = value[i];
				div_rem(q, r, a[i], MIN);
				a[i].assign(r);
				for (j = 0; j < n; j++) {
					LiDIA::multiply(TMP, q, Ttmp2[j]);
					LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
				}
			}
	} while (SW == 1);

	Ttmp2 = value[index];

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



#endif



//
// mgcd3
//

bigint *matrix< bigint >::
mgcd3(const bigint * aconst, lidia_size_t n)
{
	//
	// Task: RES = mgcd3(a,n,A);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	//              => T*a = RES
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements,
	//               index assignments
	// PAPER: Hermite normal form computation for lidia_size_teger matrices, Havas
	// Version: 1.8
	//

	debug_handler("multiple_gcd", "in member - function "
		      "mgcd3(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, bound, SW;
	bigint MIN, TMP, q, r, *Ttmp2 = NULL, *Ttmp1;

	if (columns != n)
		set_no_of_columns(n);
	if (rows != n)
		set_no_of_rows(n);
	diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd3 :: "
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
		// Pivot search: MINIMUM
		for (i = bound; i < n; i++)
			if ((abs(a[index]) > abs(a[i])) && !a[i].is_zero())
				index = i;

		MIN.assign(a[index]);

		// all elements
		SW = 0;
		Ttmp2 = value[index];

		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW = 1;
				Ttmp1 = value[i];
				div_rem(q, r, a[i], MIN);
				a[i].assign(r);
				for (j = 0; j < n; j++) {
					LiDIA::multiply(TMP, q, Ttmp2[j]);
					LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
				}
			}
	} while (SW == 1);

	Ttmp2 = value[index];

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



//
// mgcd4
//

bigint *matrix< bigint >::
mgcd4(const bigint * aconst, lidia_size_t n)
{
	//
	// Task: RES = mgcd4(a,n,A);
	//              => RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	//              => RES[0] = gcd(a[0],...,a[n-1])
	//              => T*a = RES
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements,
	//               MIN assignments, best remainder
	// PAPER: Hermite normal form computation for lidia_size_teger matrices,
	//        Havas + best remainder
	// Version: 1.8
	//

	debug_handler("multiple_gcd", "in member - function "
		      "mgcd4(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, bound, SW;
	bigint MIN, TMP, q, r, *Ttmp1, *Ttmp2 = NULL;

	if (columns != n)
		set_no_of_columns(n);
	if (rows != n)
		set_no_of_rows(n);
	diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "in mgcd4 :: Error in memory allocation (a) !!");

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
		Ttmp2 = value[index];
		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW = 1;

				Ttmp1 = value[i];
				pos_div_rem(q, r, a[i], MIN);
				if (r * bigint(2) > abs(MIN)) {
					if (!MIN.is_lt_zero()) {
						inc(q);
						LiDIA::subtract(r, r, MIN);
					}
					else {
						dec(q);
						LiDIA::add(r, r, MIN);
					}
				}

				a[i].assign(r);
				for (j = 0; j < n; j++) {
					LiDIA::multiply(TMP, q, Ttmp2[j]);
					LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
				}
			}
	} while (SW == 1);

	Ttmp2 = value[index];

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



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
