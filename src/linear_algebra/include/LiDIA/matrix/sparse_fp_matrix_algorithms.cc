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


#ifndef LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_CC_GUARD_
#define LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_CC_GUARD_


#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
# include	"LiDIA/random_generator.h"
#endif
#include	"LiDIA/modular_operations.inl"
#ifndef LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_H_GUARD_
# include	"LiDIA/matrix/sparse_fp_matrix_algorithms.h"
#endif

#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
# include	"LiDIA/bigint_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// multiply
//

math_vector< bigint >
operator * (const base_power_product< ring_matrix < bigint >, lidia_size_t > &A,
	    const math_vector< bigint > &v)
{
	math_vector< bigint > w = v;
	lidia_size_t len = A.get_no_of_components();
	for (lidia_size_t i = 0; i < len; i++)
		w = A.get_base(i) * w;
	return w;
}



math_vector< long >
operator * (const base_power_product< ring_matrix < long >, lidia_size_t > &A,
	    const math_vector< long > &v)
{
	math_vector< long > w = v;
	lidia_size_t len = A.get_no_of_components();
	for (lidia_size_t i = 0; i < len; i++)
		w = A.get_base(i) * w;
	return w;
}



template< class T, class MATRIX_TYPE >
inline const T
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::multiply_special(const math_vector< T > &w,
								const math_vector< T > &v,
								const T &mod) const
{
	T TMP, RES = 0;
	for (register lidia_size_t i = 0; i < v.size(); i++) {
		mult_mod(TMP, v[i], w[i], mod);
		add_mod(RES, RES, TMP, mod);
	}
	return RES;
}



template< class T, class MATRIX_TYPE >
void
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::multiply_special (math_vector< T > &w,
								 const base_power_product< ring_matrix < T >, lidia_size_t > &B,
								 const T &mod) const
{
	lidia_size_t *itmp;
	T *vtmp;

	T TMP, RES;
	for (lidia_size_t l = B.get_no_of_components() - 1; l >= 0; l--) {
		math_vector< T > v((B.get_base(l)).rows, (B.get_base(l)).rows);
		lidia_size_t i = 0;
		for (i = 0; i < (B.get_base(l)).rows; i++) {
			itmp = (B.get_base(l)).index[i];
			vtmp = (B.get_base(l)).value[i];

			RES = 0;
			for (lidia_size_t j = 0; j < (B.get_base(l)).value_counter[i]; j++) {
				if (vtmp[i] == 1)
					add_mod(RES, RES, w[itmp[j]], mod);
				else {
					mult_mod(TMP, vtmp[j], w[itmp[j]], mod);
					add_mod(RES, RES, TMP, mod);
				}
			}
			v[i] = RES;
		}
		w = v;
	}
}



template< class T, class MATRIX_TYPE >
void
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::multiply_special (math_vector< T > &w,
								 const MATRIX_TYPE &B,
								 const T &mod) const
{
	lidia_size_t *itmp;
	T *vtmp;
	T TMP, RES;

	math_vector< T > v(B.rows, B.rows);
	lidia_size_t i = 0;
	for (i = 0; i < B.rows; i++) {
		itmp = B.index[i];
		vtmp = B.value[i];

		RES = 0;
		for (lidia_size_t j = 0; j < B.value_counter[i]; j++) {
			if (vtmp[i] == 1)
				add_mod(RES, RES, w[itmp[j]], mod);
			else if (vtmp[i] == -1)
				sub_mod(RES, RES, w[itmp[j]], mod);
			else {
				mult_mod(TMP, vtmp[j], w[itmp[j]], mod);
				add_mod(RES, RES, TMP, mod);
			}
		}
		v[i] = RES;
	}
	w = v;
}



template< class T, class MATRIX_TYPE >
void
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::multiply_special(math_vector< T > &v,
								const math_vector< T > &w,
								const MATRIX_TYPE &B,
								const T &mod) const
{
	if (v.get_size() < w.get_size()) {
		v.set_capacity(w.get_size());
		v.set_size(w.get_size());
	}

	lidia_size_t *itmp;
	T *vtmp;
	T TMP, RES;
	lidia_size_t i = 0;
	for (i = 0; i < B.rows; i++) {
		itmp = B.index[i];
		vtmp = B.value[i];

		RES = 0;
		for (lidia_size_t j = 0; j < B.value_counter[i]; j++) {
			if (vtmp[j] == 1)
				add_mod(RES, RES, w[itmp[j]], mod);
			else {
				mult_mod(TMP, vtmp[j], w[itmp[j]], mod);
				add_mod(RES, RES, TMP, mod);
			}
		}
		v[i] = RES;
	}
}



//
// berlekamp massay algorithm
//

template< class T, class MATRIX_TYPE >
T *
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::berlekamp_massay(T *s, lidia_size_t n, const T &mod) const
{
	// Memory allocation
	T D, Dold, TMP2;
	T *C = new T[n + 1];
	T *B = new T[n + 1];
	T *tmp = new T[n + 1];
	lidia_size_t m, N, i, L;

	for (i = 0; i <= n; i++) {
		C[i] = 0;
		B[i] = 0;
		tmp[i] = 0;
	}

	// Initialization
	C[0] = 1;
	L = 0;
	m = -1;
	B[0] = 1;
	N = 0;
	Dold = 1;

	// Iteration
	while (N < n) {
		//Compute the next discrepancy
		best_remainder(D, s[N], mod);
		for (i = 1; i <= L; i++) {
			//D += (C[i]*s[N-i]); D %= mod;
			mult_mod(TMP2, C[i], s[N-i], mod);
			add_mod(D, D, TMP2, mod);
		}

		if (D != 0) {
			for (i = 0; i <= L; i++)
				tmp[i] = C[i];

			for (i = 0; i <= L; i++) {
				inv_mod(TMP2, Dold, mod);
				mult_mod(TMP2, D, TMP2, mod);
				mult_mod(TMP2, TMP2, B[i], mod);
				sub_mod(C[i + N - m], C[i + N - m], TMP2, mod); // -= (D*TMP2)*B[i]; C[i + N - m] %= mod;
			}
			if (L*2 <= N+1) {
				L = N + 1 - L;
				m = N;
				for (i = 0; i <= L; i++)
					B[i] = tmp[i];
				Dold = D;
			}
		}
		N++;
	}
	delete[] B;
	delete[] tmp;

	T *RES = new T[L+2];
	RES[0] = L;

	for (i = 0; i <= L; i++)
		RES[i+1] = C[L-i];
	delete[] C;
	return RES;
}



//
// column step form
//

template< class T, class MATRIX_TYPE >
inline int
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::STF(MATRIX_TYPE &, const T &) const
{
	return 0;
}



//
// rank
//

template< class T, class MATRIX_TYPE >
inline lidia_size_t
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::rank(MATRIX_TYPE &A, const T &) const
{
	return A.columns;
}



//
// ran kand linearly independent rows or columns
//

template< class T, class MATRIX_TYPE >
inline lidia_size_t *
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::lininr(MATRIX_TYPE &, const T &) const
{
	return NULL;
}



template< class T, class MATRIX_TYPE >
inline lidia_size_t *
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::lininc(MATRIX_TYPE &, const T &) const
{
	return NULL;
}



//
// adjoint matrix
//

template< class T, class MATRIX_TYPE >
inline void
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::adj(MATRIX_TYPE &, const T &) const
{

}



//
// determinant
//

template< class T, class MATRIX_TYPE >
inline const T
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::det(MATRIX_TYPE &A, const T &mod) const
{
	T *data = charpoly(A, mod);
	T ret;
	if (A.rows % 2 == 1)
		ret = -data[0];
	else
		ret = data[0];
	delete[] data;
	best_remainder(ret, ret, mod);
	return ret;
}



template< class T, class MATRIX_TYPE >
inline const T
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::det (const base_power_product< ring_matrix < T >, lidia_size_t > &A,
						    const T &mod) const
{
	T *data = charpoly(A, mod);
	T ret;
	if (A.get_base(0).rows % 2 == 1)
		ret = -data[0];
	else
		ret = data[0];
	delete[] data;
	return ret;
}



//
// Hessenberg form
//

template< class T, class MATRIX_TYPE >
inline void
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::HBF (MATRIX_TYPE &,
						    const T &) const
{
}



//
// characteristic polynomial
//

template< class T, class MATRIX_TYPE >
T *
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::charpoly (MATRIX_TYPE &A,
							 const T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t n = A.rows;

	random_generator gen;
	long TMPlong;
	Fp_polynomial res;
	res.set_modulus(mod);
	math_vector< T > u(n, n), b(n, n);
	bool SW = false;
	math_vector< T > w(n, n);
	T *s = new T[n+n];

	lidia_size_t i;
	do {
		for (i = 0; i < n; i++) {
			gen >> TMPlong;
			best_remainder(u[i], TMPlong, mod);
			gen >> TMPlong;
			best_remainder(b[i], TMPlong, mod);
		}

		w = b;
		//s[0] = (w*u) % mod;
		s[0] = multiply_special(w, u, mod);


		for (i = 1; i < n+n; i++) {
			//w = (A * w) % mod;
			multiply_special(w, A, mod);
			//s[i] = (u * w) % mod;
			s[i] = multiply_special(u, w, mod);
		}

		T *c = berlekamp_massay(s, n+n, mod);


		bigint LEN = c[0];
		lidia_size_t len;
		LEN.sizetify(len);

		polynomial< bigint > va(&(c[1]), len);

		Fp_polynomial vb(va, mod);
		delete[] c;

		if (res.degree() <= 0)
			res = vb;
		else {
			res = (res * vb)/gcd(res, vb);

			if (res == vb)
				SW = true;
		}
		//std::cout << "grad = " << res.degree() << std::endl;
	}
	while (res.degree() != n && SW == false);
	delete[] s;

	return Fp_polynomial_convert(res, mod);
}



template< class T, class MATRIX_TYPE >
T *
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::charpoly (const base_power_product< ring_matrix < T >, lidia_size_t > &A,
							 const T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t n = A.get_base(0).rows;
	random_generator gen;
	long TMPlong;
	Fp_polynomial res;
	res.set_modulus(mod);
	bool SW = false;
	math_vector< T > u(n, n), b(n, n);

	math_vector< T > w(n, n);
	T *s = new T[n+n];

	lidia_size_t i;
	do {
		for (i = 0; i < n; i++) {
			gen >> TMPlong;
			best_remainder(u[i], TMPlong, mod);
			gen >> TMPlong;
			best_remainder(b[i], TMPlong, mod);
		}

		w = b;
		//s[0] = (w * u);
		s[0] = multiply_special(w, u, mod);

		for (i = 1; i < n+n; i++) {
			//w = A * w;
			multiply_special(w, A, mod);
			//s[i] = (u * w) % mod;
			s[i] = multiply_special(u, w, mod);
		}

		T *c = berlekamp_massay(s, n+n, mod);
		bigint LEN = c[0];
		lidia_size_t len;
		LEN.sizetify(len);

		polynomial< bigint > va(&(c[1]), len);
		Fp_polynomial vb(va, mod);

		delete[] c;

		if (res.degree() <= 0)
			res = vb;
		else {
			res = (res * vb)/gcd(res, vb);
			if (res == vb)
				SW = true;
		}
	}
	while (res.degree() != n && SW == false);

	delete[] s;

	return Fp_polynomial_convert(res, mod);
}



template< class T, class MATRIX_TYPE >
bool
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::lanczos (const MATRIX_TYPE &A,
							math_vector< T > &x,
							const math_vector< T > &b,
							T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t i, n = A.columns;

	math_vector< T > *w = new math_vector< T > [n+n];
	math_vector< T > *v = new math_vector< T > [n+n];

	for (i = 0; i < n+n; i++) {
		w[i].set_capacity(n);
		w[i].set_size(n);
		v[i].set_capacity(n);
		v[i].set_size(n);
	}

	T *wv = new T[n+n];

	// Init
	w[0] = b % mod;
	multiply_special(v[1], w[0], A, mod);

	T TMP1, TMP2, TMP3, TMP4, TMP = multiply_special(w[0], v[1], mod);
	inv_mod(wv[0], TMP, mod);

	TMP1 = multiply_special(v[1], v[1], mod);

	//w[1] = (v[1] - (TMP1*wv[0])*w[0]) % mod;
	mult_mod(TMP3, TMP1, wv[0], mod);
	for (i = 0; i < n; i++) {
		mult_mod(TMP2, TMP3, w[0][i], mod);
		sub_mod(w[1][i], v[1][i], TMP2, mod);
	}

	//x = ((w[0]*b)*wv[0]*w[0]) % mod;
	TMP2 = multiply_special(w[0], b, mod);
	mult_mod(TMP2, TMP2, wv[0], mod);
	for (i = 0; i < n; i++)
		mult_mod(x[i], TMP2, w[0][i], mod);

	i = 1;
	lidia_size_t j;

	while (TMP != 0) {
		//std::cout << "Lanczos Iteration: " << i << std::endl;
		multiply_special(v[i+1], w[i], A, mod);

		TMP = multiply_special(w[i], v[i+1], mod);
		if (TMP != 0) {
			inv_mod(wv[i], TMP, mod);

			TMP1 = multiply_special(v[i+1], v[i+1], mod);
			mult_mod(TMP1, TMP1, wv[i], mod);
			TMP2 = multiply_special(w[i], v[i+1], mod);
			mult_mod(TMP2, TMP2, wv[i-1], mod);

			//w[i+1] = (v[i+1] - TMP1*w[i] - TMP2*w[i-1]) % mod;
			for (j = 0; j < n; j++) {
				mult_mod(TMP3, TMP1, w[i][j], mod);
				mult_mod(TMP4, TMP2, w[i-1][j], mod);
				add_mod(TMP3, TMP3, TMP4, mod);
				sub_mod(w[i+1][j], v[i+1][j], TMP3, mod);
			}
			// Update of solution
			//x = (x+((w[i]*b)*wv[i])*w[i]) % mod;
			TMP2 = multiply_special(w[i], b, mod);
			mult_mod(TMP2, TMP2, wv[i], mod);
			for (j = 0; j < n; j++) {
				mult_mod(TMP3, TMP2, w[i][j], mod);
				add_mod(x[j], x[j], TMP3, mod);
			}
			i++;
		}
	}
	delete[] w;
	delete[] v;
	delete[] wv;
	return true;
}



template< class T, class MATRIX_TYPE >
T
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::lanczos_ZmZ (const MATRIX_TYPE &A,
							    math_vector< T > &x,
							    const math_vector< T > &b,
							    T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t i, n = A.columns;
	T factor_1, factor_2, factor_3, factor_1_2;

	math_vector< T > *w = new math_vector< T > [n+n];
	math_vector< T > *v = new math_vector< T > [n+n];

	for (i = 0; i < n+n; i++) {
		w[i].set_capacity(n);
		w[i].set_size(n);
		v[i].set_capacity(n);
		v[i].set_size(n);
	}

	// Init
	w[0] = b % mod;

	multiply_special(v[1], w[0], A, mod);

	T TMP1, TMP2, TMP3, TMP4;
	T TMP = multiply_special(w[0], v[1], mod);

	// FAKTOR_1, FACTOR_3
	factor_1 = TMP;
	factor_3 = factor_1;

	TMP1 = multiply_special(v[1], v[1], mod);

	//w[1]
	for (i = 0; i < n; i++) {
		mult_mod(TMP2, TMP1, w[0][i], mod);
		mult_mod(TMP3, TMP, v[1][i], mod);
		sub_mod(w[1][i], TMP3, TMP2, mod);
	}

	//FACTOR_2
	mult_mod(factor_2, factor_1, factor_1, mod);

	//x
	TMP2 = multiply_special(w[0], b, mod);
	for (i = 0; i < n; i++)
		mult_mod(x[i], TMP2, w[0][i], mod);

	i = 1;
	lidia_size_t j;

	while (TMP != 0) {
		multiply_special(v[i+1], w[i], A, mod);

		TMP = multiply_special(w[i], v[i+1], mod);

		if (TMP != 0) {

			TMP1 = multiply_special(v[i+1], v[i+1], mod);
			TMP2 = multiply_special(w[i], v[i+1], mod);

			//FACTOR_1
			factor_1 = TMP2;
			mult_mod(factor_3, factor_3, factor_2, mod);
			mult_mod(factor_1_2, factor_2, factor_1, mod);

			mult_mod(TMP1, TMP1, factor_2, mod);
			mult_mod(TMP2, TMP2, factor_1, mod);

			//w[i+1]
			for (j = 0; j < n; j++) {
				mult_mod(TMP3, TMP1, w[i][j], mod);
				mult_mod(TMP4, TMP2, w[i-1][j], mod);
				add_mod(TMP3, TMP3, TMP4, mod);
				mult_mod(TMP4, factor_1_2, v[i+1][j], mod);
				sub_mod(w[i+1][j], TMP4, TMP3, mod);
			}

			// Update of solution
			TMP2 = multiply_special(w[i], b, mod);
			for (j = 0; j < n; j++) {
				mult_mod(TMP3, TMP2*factor_3, w[i][j], mod);
				mult_mod(TMP4, factor_1_2, x[j], mod);
				add_mod(x[j], TMP4, TMP3, mod);
			}

			//FACTOR_2
			mult_mod(factor_3, factor_3, factor_1, mod);
			mult_mod(factor_2, factor_2, factor_1, mod);
			mult_mod(factor_2, factor_2, factor_1, mod);

			i++;
		}
	}

	delete[] w;
	delete[] v;
	return factor_3;
}



//
// wiedemann algorithm
//

template< class T, class MATRIX_TYPE >
bool
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::wiedemann (const ring_matrix< T > &A,
							  math_vector< T > &x,
							  const math_vector< T > &b,
							  T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t n = A.columns;

	// vector u
	math_vector< T > u(n, n);
	long TMPlong;
	random_generator gen;
	lidia_size_t i;

	for (i = 0; i < n; i++) {
		gen >> TMPlong;
		u[i] = TMPlong % mod;
	}

	math_vector< T > *w = new math_vector< T > [n+n];
	T *s = new T[n+n];

	w[0] = b;
	s[0] = (w[0]*u);

	for (i = 1; i < n+n; i++) {
		w[i] = A * w[i-1] % mod;
		s[i] = (u* w[i]) % mod;
	}

	T *c = berlekamp_massay(s, n+n, mod);

	for (i = 0; i < x.size(); i++)
		x[i] = 0;
	for (i = 2; i <= c[0]+1; i++) {
		x += c[i]*w[i-2];
	}

	T TMP;
	inv_mod(TMP, c[1], mod);
	x = x * (-TMP) % mod;
	return true;
}



template< class T, class MATRIX_TYPE >
bool
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::wiedemann (const base_power_product< ring_matrix < T >, lidia_size_t > &A,
							  math_vector< T > &x,
							  const math_vector< T > &b,
							  T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t n = A.get_base(0).columns;

	// vector u
	math_vector< T > u(n, n);
	long TMPlong;
	random_generator gen;
	lidia_size_t i;

	for (i = 0; i < n; i++) {
		gen >> TMPlong;
		u[i] = TMPlong % mod;
	}

	math_vector< T > *w = new math_vector< T > [n+n];
	T *s = new T[n+n];

	w[0] = b;
	s[0] = (w[0]*u);

	for (i = 1; i < n+n; i++) {
		w[i] = A * w[i-1] % mod;
		s[i] = (u* w[i]) % mod;
	}

	T *c = berlekamp_massay(s, n+n, mod);

	for (i = 0; i < x.size(); i++)
		x[i] = 0;

	for (i = 2; i <= c[0]+1; i++) {
		x += c[i]*w[i-2];
	}

	T TMP;
	inv_mod(TMP, c[1], mod);
	x = x * (-TMP) % mod;
	return true;
}



template< class T, class MATRIX_TYPE >
bool
sparse_fp_matrix_algorithms< T, MATRIX_TYPE >::conjugate_gradient (const ring_matrix< T > &A,
								   math_vector< T > &x,
								   const math_vector< T > &b,
								   T &mod) const
{
	// A symmetrische n x n Matrix
	lidia_size_t n = A.columns;

	// vector u
	math_vector< T > u(n, n), p(n, n), r(n, n);
	long TMPlong;
	random_generator gen;
	for (register lidia_size_t i = 0; i < n; i++) {
		gen >> TMPlong;
		u[i] = TMPlong % mod;
	}

	T TMP, TMP1, a, beta;
	// Init
	p = (b - A*u) % mod;
	r = p;

	TMP = (r*r);
	while (TMP != 0) {
		// PART 1
		TMP1 = (p * (A * p)) % mod;

		inv_mod(TMP1, TMP1, mod);
		a = TMP*TMP1;

		u = (u + a * p) % mod;
		r = (r - a*A*p) % mod;

		// PART 2
		inv_mod(TMP1, TMP, mod);
		TMP = (r*r);
		beta = (TMP*TMP1) % mod;

		p = (r + beta*p) % mod;

	}
	x = u;
	return true;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_FP_MATRIX_ALGORITHMS_CC_GUARD_
