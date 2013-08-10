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
//	Author	: Thomas Pfahler (TPf)
//                based on code written by Victor Shoup
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/Fp_poly_modulus.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void null_space(lidia_size_t & r, base_vector< sdigit > & D,
		Fp_polynomial** & M, lidia_size_t n, const galois_field& K)
{
	debug_handler("factoring.c", "null_space(lidia_size_t&, base_vector< sdigit > &, Fp_polynomial**&, lidia_size_t, galois_field&)");

	//input : (n x n)-matrix M, modulus p, vector D
	//output : "triangularized" "compressed" matrix M with same nullspace
	//                      as input-matrix M, empty rows are deleted (?)
	//              r = dimension of nullspace of M
	//              D[i] == -1 iff row i of M is zero
	//              D[i] == m iff M[m][m]!=0 and M[m][j]==0 for all 0<=j<m

	lidia_size_t k, l;
	lidia_size_t i, j;
	lidia_size_t pos;
	Fp_polynomial t1, t2;
	Fp_polynomial *x, *y;
	bool verbose = single_factor< gf_polynomial >::verbose();

	Fp_polynomial minus_one;
	minus_one.set_modulus(K.characteristic());
	minus_one.assign_one();
	minus_one.negate();

	Fp_poly_modulus F(K.irred_polynomial());

	if (D.capacity() < n)
		D.set_capacity(n);
	D.set_size(n);

	for (j = 0; j < n; j++)
		D[j] = -1;


	r = 0;
	l = 0;
	for (k = 0; k < n; k++) {
		if (verbose && k % 10 == 0) std::cerr << "+";

		pos = -1;
		for (i = l; i < n; i++) {
			if (pos == -1 && !M[i][k].is_zero())
				pos = i;
		}

		if (pos != -1) {
			//swap(M[pos], M[l]) :
			x = M[pos];
			M[pos] = M[l];
			M[l] = x;

			// make M[l, k] == -1, and make row l reduced
			invert_mod(t1, M[l][k], K.irred_polynomial());
			t1.negate(); //t1 = -1/M[l][k]
			M[l][k].assign(minus_one);
			for (j = k + 1; j < n; j++) {
				remainder(t2, M[l][j], F);
				multiply(M[l][j], t2, t1, F);
			}

			for (i = l + 1; i < n; i++) {
				// M[i] = M[i] + M[l]*M[i,k]

				t1.assign(M[i][k]); // this is already reduced

				if (!t1.is_zero()) {
					x = &M[i][k+1];
					y = &M[l][k+1];

					for (j = k + 1; j < n; j++, x++, y++) {
						// *x = *x + (*y)*t1

						multiply(t2, *y, t1);
						add(*x, *x, t2);
					}
				}

				M[i][k].assign_zero(); //this has been eliminated

			}

			D[k] = l; // variable k is defined by row l

			l++;

		}
		else {
			r++;
		}

		//reduce next column of remaining rows
		if (k != n - 1)
			for (j = l; j < n; j++)
				remainder(M[j][k + 1], M[j][k + 1], F);
	}
	//one could delete rows i with D[i] == -1, couldn't one ?
}


void build_matrix(Fp_polynomial** & M, lidia_size_t n,
		  const gf_polynomial & g, const gf_poly_modulus & F)
{
	debug_handler("gf_polynomial", "build_matrix(Fp_polynomial**&, lidia_size_t, gf_polynomial&, gf_poly_modulus&, lidia_size_t)");

	//input : polynomial g, deg(g) < n; gf_poly_modulus F, deg(F) < n
	//output: (n x n)-matrix M, where column i is equivalent to g^i mod F   (i = 0..n-1)

	if (g.degree() >= n || F.modulus().degree() > n)
		lidia_error_handler("gf_polynomial", "build_matrix(...)::degree of polynomil too big");

	lidia_size_t i, j, m;
	bool verbose = single_factor< gf_polynomial >::verbose();

	M = new Fp_polynomial *[n];
	if (!M) lidia_error_handler("gf_polynomial",
				    "build_matrix(...)::out of memory");

	for (i = 0; i < n; i++) {
		M[i] = new Fp_polynomial[n];
		if (!M[i]) lidia_error_handler("gf_polynomial",
					       "build_matrix(...)::out of memory");
	}

	const galois_field &K = F.modulus().get_field();
	Fp_polynomial ZERO;
	ZERO.set_modulus(K.characteristic());

	gf_polynomial h;
	h.assign_one(K);

	for (j = 0; j < n; j++) {
		if (verbose && j % 10 == 0) std::cerr << "+";

		m = h.degree();
		for (i = 0; i < n; i++) {
			if (i <= m) M[i][j].assign(h[i].polynomial_rep());
			else (M[i][j]).assign(ZERO);
		}

		if (j < n - 1)
			multiply(h, h, g, F);

	}

	//we might move this into the function null_space
	ZERO.assign_one(); //I know, not the best name for this variable...
	for (i = 0; i < n; i++)
		subtract(M[i][i], M[i][i], ZERO);
}



void random_basis_elt2(gf_polynomial & g, const base_vector< sdigit > &D, Fp_polynomial** M, const galois_field &K)
{
	debug_handler("gf_polynomial", "random_basis_elt2(gf_polynomial&, base_vector< sdigit > &, Fp_polynomial**, lidia_size_t, galois_field&)");

	//input : vector D, Matrix M, modulus p as output of function null_space
	//output: g = random basis element of ker M

	Fp_polynomial t1, t2;
	t1.set_modulus(K.characteristic());
	lidia_size_t n = D.size();
	lidia_size_t i, j, s;

	g.ffield = K;
	gf_polynomial::build_frame(K);
	g.set_degree(n - 1);

	for (j = n - 1; j >= 0; j--) {
		if (D[j] == -1)  g[j].randomize();
		else {
			i = static_cast<lidia_size_t>(D[j]);

			// g[j] = sum_{s=j+1}^{n-1} v[s]*M[i,s]

			t1.assign_zero();
			for (s = j + 1; s < n; s++) {
				multiply(t2, g[s].polynomial_rep(), M[i][s]);
				add(t1, t1, t2);
			}
			g[j].set_polynomial_rep(t1);
		}
	}

	g.remove_leading_zeros();
	gf_polynomial::delete_frame();
}



void
sf_berlekamp_work(factorization< gf_polynomial > &factors,
		  const gf_polynomial &X_to_the_power_of_q, const gf_poly_modulus &F)
	// Assumes f is square-free, monic, deg(f) > 0 and f(0) != 0 .
	// returns list of factors of f.
{
	debug_handler("gf_polynomial", "sf_berlekamp(base_vector< gf_polynomial > &, gf_polynomial&, gf_polynomial&)");

	const gf_polynomial &f = F.modulus();

	if (f.const_term().is_zero())
		lidia_error_handler("gf_polynomial", "sf_berlekamp::const_term == 0");

	factors.kill();

	bool verbose = single_factor< gf_polynomial >::verbose();
	my_timer t;

	lidia_size_t n = f.degree();
	const galois_field &K = f.get_field();

	gf_polynomial g, h;

	base_vector< sdigit > D;
	lidia_size_t r;
	Fp_polynomial** M;

	if (verbose) t.start("building matrix...");
	build_matrix(M, n, X_to_the_power_of_q, F);
	if (verbose) t.stop();

	if (verbose) t.start("diagonalizing...");
	null_space(r, D, M, n, K);
	if (verbose) t.stop();

	if (verbose) std::cerr << "number of factors = " << r << "\n";

	if (r == 1) {
		append_irred_factor(factors, f);
		for (r = 0; r < n; r++) delete[] M[r];
		delete[] M;
		return;
	}

	const bigint &p = K.characteristic();
	const bigint &q = K.number_of_elements();

	int method;
	//if (ext_method<0)
	//    method = -ext_method;
	//else

	if (q > 2*r)
		method = 3;
	else
		method = 1;

	//    if (f.degree() < f.modulus())  method = 0;
	//    else                           method = 1;

	lidia_size_t i, j;
	switch(method) {
	case(3) :
	{
		if (verbose) t.start("computing splitting polynomial...");
		do {
			if (verbose) std::cerr << "+";
			random_basis_elt2(g, D, M, K);
		} while (!checked_min_poly(h, g, r, F));
		if (verbose) t.stop();

		for (i = 0; i < n; i++) delete[] M[i];
		delete[] M;

		base_vector< gf_element > roots;

		if (verbose) t.start("finding roots of splitting polynomial...");
		roots.assign(find_roots(h));
		if (verbose) t.stop();

		if (verbose) t.start("finding factors...");
		find_irred_factors(factors, f, g, roots);
		if (verbose) t.stop();
	}; break;
	case(2) :
	{
		if (!is_int(q) || p != q) lidia_error_handler("gf_polynomial",
							      "sf_berlekamp(...)::field too large");

		if (verbose) t.start("computing gcd's ...");

		gf_element ONE;
		ONE.assign_one(K);
		//we could use a simple method to get ALL elements in a finite field:
		//use all polynomials of degree < pb.field().degree()
		//but we think this method is terribly slow even for small fields

		gf_polynomial *container;
		container = new gf_polynomial[r];
		if (!container) lidia_error_handler("factoring.c",
						    "sf_berlekamp(...)::out of memory");
		container[0].assign(f);
		lidia_size_t len = 1;

		while (len < r) {
			if (verbose) std::cerr << "+";

			random_basis_elt2(g, D, M, K);
			for (i = 0; (i < q) && (len < r); i++) {
				for (j = 0; (j < len) && (len < r); j++) {
					gcd(h, container[j], g);
					if (!h.is_one() && h != container[j]) {
						container[len] = h;
						len++;
						divide(container[j], container[j], h);
					}
				}
				add(g, g, ONE);
			}
		}
		for (i = 0; i < n; i++) delete[] M[i];
		delete[] M;

		for (i = 0; i < r; i++)
			append_irred_factor(factors, container[i]);
		delete[] container;

		if (verbose) t.stop();
	}; break;
	case(1) :
	{
		if (verbose) t.start("direct powering method...");

		gf_element c(K), ONE(K);
		ONE.assign_one();

		bigint q2(q);
		q2.divide_by_2(); //q2 = (q-1)/2

		gf_polynomial tmp, hh;
		gf_polynomial *container = new gf_polynomial[r];
		if (!container) lidia_error_handler("factoring.c",
						    "sf_berlekamp(...)::out of memory");
		container[0].assign(f);
		lidia_size_t len = 1;

		do {
			if (verbose) std::cerr << "+";
			random_basis_elt2(g, D, M, K);
			if (p != 2) {
				power(hh, g, q2, F);
				add(h, hh, ONE);
			}
			else
				trace_map(h, g, K.degree(), F, X_to_the_power_of_q);


			for (i = 0; (i < len) && (len < r); i++) {
				gcd(tmp, h, container[i]);
				if (!tmp.is_one() && tmp != container[i]) {
					divide(container[len], container[i], tmp);
					container[i] = tmp;
					len++;
				}
				if (len < r) {
					gcd(tmp, hh, container[i]);
					if (!tmp.is_one() && tmp != container[i]) {
						divide(container[len], container[i], tmp);
						container[i] = tmp;
						len++;
					}
				}
			}
		} while (len < r);

		for (i = 0; i < n; i++) delete[] M[i];
		delete[] M;

		for (i = 0; i < r; i++)
			append_irred_factor(factors, container[i]);
		delete[] container;

		if (verbose) t.stop();
	}; break;
	default : lidia_error_handler("factoring.c",
				      "sf_berlekamp::wrong index for method");
	};

	if (verbose) {
		std::cerr << "degrees:";
		for (i = 0; i < factors.no_of_prime_components(); i++)
			std::cerr << " " << factors.prime_base(i).base().degree();
		std::cerr << "\n";
	}
}



void
sf_berlekamp(factorization< gf_polynomial > &factors, const gf_polynomial & f)
{
	bool verbose = single_factor< gf_polynomial >::verbose();
	my_timer t;
	gf_polynomial g;
	gf_poly_modulus F(f);
	if (verbose) t.start("computing X^q...");
	power_x(g, f.get_field().number_of_elements(), F);
	if (verbose) t.stop();
	sf_berlekamp_work(factors, g, F);
}


void
berlekamp(factorization< gf_polynomial > &factors, const gf_polynomial & f)
	// f must be monic
	// returns a list of factors, with multiplicities.
{
	debug_handler("gf_polynomial", "berlekamp(factorization< gf_polynomial > &, gf_polynomial&)");

	if (f.is_zero())
		lidia_error_handler("gf_polynomial", "berlekamp(...)::input is zero polynomial");

	bool verbose = single_factor< gf_polynomial >::verbose();

	my_timer t;
	factorization< gf_polynomial > sfd, x;
	gf_polynomial ff;
	gf_element lc = f.lead_coeff();
	gf_element lc_inv = inverse(lc);
	multiply(ff, f, lc_inv);

	if (verbose) t.start("square-free decomposition...");
	square_free_decomp(sfd, ff);
	if (verbose) t.stop();

	lidia_size_t i;
	factors.kill();
	gf_polynomial tmp, X(f.get_field());
	X.assign_x();

	for (i = 0; i < sfd.no_of_composite_components(); i++) {
		if (verbose) {
			std::cerr << "factoring multiplicity " << sfd.composite_exponent(i);
			std::cerr << ", deg = " << sfd.composite_base(i).base().degree() << "\n";
		}

		// check for trivial factor 'X'  
		if (!sfd.composite_base(i).base().const_term().is_zero())
			sf_berlekamp(x, sfd.composite_base(i).base());
		else {
			divide(tmp, sfd.composite_base(i).base(), X);
			if (tmp.degree() > 0)
				sf_berlekamp(x, tmp);
			else
				x.kill();

			// append trivial factor  
			append_irred_factor(factors, X, sfd.composite_exponent(i));
		}
		x.power(sfd.composite_exponent(i));

		//append factors
		multiply(factors, factors, x);
	}
	ff.set_degree(0);
	ff[0].assign(lc);
	factors.append(ff);
}



factorization< gf_polynomial > berlekamp(const gf_polynomial& f)
{
	factorization< gf_polynomial > F;
	berlekamp(F, f);
	return F;
}



factorization< gf_polynomial >
single_factor< gf_polynomial >::berlekamp() const
{
	factorization< gf_polynomial > F;
	LiDIA::berlekamp(F, rep);
	return F;
}



factorization< gf_polynomial >
single_factor< gf_polynomial >::sf_berlekamp() const
{
	factorization< gf_polynomial > F;
	LiDIA::sf_berlekamp(F, rep);
	return F;
}



factorization< gf_polynomial > sf_berlekamp(const gf_polynomial& f)
{
	factorization< gf_polynomial > F;
	sf_berlekamp(F, f);
	return F;
}



bool
det_irred_test(const gf_polynomial &f)
	//irreducibility test according to Berlekamp
{
	debug_handler("gf_polynomial", "det_irred_test(...)");

	if (f.degree() < 1)
		return false;
	if (f.degree() == 1)
		return true;
	if (f.const_term().is_zero())
		return false;

	gf_polynomial g;
	derivative(g, f);
	gcd(g, f, g);
	if (!g.is_one())
		return false;

	lidia_size_t i, n = f.degree();
	const galois_field &K = f.get_field();
	const bigint &q = K.number_of_elements();
	gf_polynomial X_to_the_power_of_q;
	gf_poly_modulus F(f);
	power_x(X_to_the_power_of_q, q, F);

	base_vector< sdigit > D;
	lidia_size_t r;
	Fp_polynomial** M;
	build_matrix(M, n, X_to_the_power_of_q, F);
	null_space(r, D, M, n, K);

	for (i = 0; i < n; i++) delete[] M[i];
	delete[] M;

	return (r == 1);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
