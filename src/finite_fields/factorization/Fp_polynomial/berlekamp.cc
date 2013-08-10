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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/Fp_poly_multiplier.h"

#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/udigit.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if 0
inline
long mul_mod_sp2(long a, long b, long p, double bpinv)
{
	register long q = static_cast<long>(static_cast<double>(a) * bpinv);
	register long res = a*b - q*p;

	res += (res >> ((32) -1)) & p;
	res -= p;
	res += (res >> ((32) -1)) & p;

	return res;
}
#endif



//***********************************************************************
//
//			class berl_matrix
//	a class for the matrix step in the Berlekamp algorithm
//	which distinguishes between small and large moduli
//
//***********************************************************************

class berl_matrix
{
	int flag; // flag = 0 for moduli that fit into an udigit
    				// flag = 1 otherwise
	udigit p_small;
	bigint p_big;
	int n;
	udigit ** M_small;
	bigint ** M_big;
	lidia_size_t * D; //a vector for row information

	lidia_size_t r; // dimension of the kernel

	void build(const Fp_polynomial &g, const Fp_poly_modulus &F);

	void nullspace_small();
	void nullspace_big();

	berl_matrix(); // disable

public:
	berl_matrix(const Fp_polynomial &g, lidia_size_t N, const Fp_poly_modulus &F);
	~berl_matrix();
	void kill();

	void nullspace(); // computes the kernel

	lidia_size_t dim_nullspace() const
	{
		return r;
	}		// returns the dimension of the kernel
    				// result is only defined if nullspace()
				// was called before

	Fp_polynomial random_basis_elt() const;
    				// returns a random element in the kernel
				// result is only defined if nullspace()
    				// was called before
};



berl_matrix::berl_matrix(const Fp_polynomial &g, lidia_size_t N,
			 const Fp_poly_modulus &F) :
	n(N), M_small(0), M_big(0), D(0), r(-1)
{
	debug_handler("berl_matrix", "constructor");

	const Fp_polynomial &f = F.modulus();
	g.comp_modulus(f, "berl_matrix::constructor");
	if (g.degree() >= n || f.degree() > n) {
		lidia_error_handler("berl_matrix", "constructor :: "
				    "degree of polynomial too large");
		return; // LC
	}

	p_big = f.modulus();
	long tmp;
	if (p_big.longify(tmp))
		flag = 1;
	else {
		flag = 0;
		p_small = tmp; //p_small has sucessfully been converted into a long
	}

	D = new lidia_size_t[n];
	memory_handler(D, "berl_matrix", "constructor :: "
		       "Error in memory allocation");

	build(g, F);
}



berl_matrix::~berl_matrix()
{
	debug_handler("berl_matrix", "destructor");
	kill();
}



void berl_matrix::kill()
{
	debug_handler("berl_matrix", "kill()");
	lidia_size_t i;
	if (M_small)
		for (i = 0; i < n; i++)
			delete[] M_small[i];
	if (M_big)
		for (i = 0; i < n; i++)
			delete[] M_big[i];
	delete[] M_small;
	delete[] M_big;
	delete[] D;
	M_small = 0;
	M_big = 0;
	D = 0;
}




void berl_matrix::build(const Fp_polynomial &g, const Fp_poly_modulus &F)
	//assumes that 'flag' and 'n' are already properly set
{
	debug_handler("berl_matrix", "build(Fp_polynomial, lidia_size_t)");

	lidia_size_t i, j, m;
	bool verbose = single_factor< Fp_polynomial >::verbose();
	Fp_polynomial h;

	//allocate space for the matrix
	if (flag == 0) {
		M_small = new udigit *[n];
		memory_handler(M_small, "berl_matrix", "build :: "
			       "Error in memory allocation (M_small)");
		for (i = 0; i < n; i++) {
			M_small[i] = new udigit[n];
			memory_handler(M_small[i], "berl_matrix", "build :: "
				       "Error in memory allocation (M_small[i])");
		}
	}
	else {
		M_big = new bigint *[n];
		memory_handler(M_big, "berl_matrix", "build :: "
			       "Error in memory allocation (M_big)");
		for (i = 0; i < n; i++) {
			M_big[i] = new bigint[n];
			memory_handler(M_big[i], "berl_matrix", "build :: "
				       "Error in memory allocation (M_big[i])");
		}
	}

	Fp_poly_multiplier G(g, F);
	h.set_modulus(g);
	h.assign_one();
	bigint tmp;
	long tmp2;

	// construct the Berlekamp matrix
	// row i "is equal" to g^i mod f
	for (j = 0; j < n; j++) {
		if (verbose && j % 10 == 0) std::cerr << "+";

		m = h.degree();
		for (i = 0; i < n; i++) {
			if (flag == 0) {
				if (i <= m) {
					tmp.assign(h[i]);

					tmp.longify(tmp2);
					M_small[i][j] = tmp2;
				}
				else        M_small[i][j] = 0;
			}
			else {
				if (i <= m) M_big[i][j].assign(h[i]);
				else        M_big[i][j].assign_zero();
			}

		}

		if (j < n - 1)
			multiply(h, h, G, F); //Fp_poly_multiplier
	}

	//decrease diagonal entries by one
	for (i = 0; i < n; i++) {
		if (flag == 0)
			M_small[i][i] =
				subtract_mod(M_small[i][i], static_cast<udigit>(1), p_small);
		else
			SubMod(M_big[i][i], M_big[i][i], bigint(1), p_big);
	}
}


void berl_matrix::nullspace()
	//input : (n x n)-matrix M, modulus p, vector D
	//output : "triangularized" "compressed" matrix
	//output : "triangularized" "compressed" matrix M with same nullspace
	//                      as imput-matrix M, empty rows are deleted (?)
	//              r = dimension of nullspace of M
	//              D[i] == -1 iff row i of M is zero
	//              D[i] == m iff M[m][m]!=0 and M[m][j]==0 for all 0<=j<m
{
	debug_handler("berl_matrix", "nullspace()");
	lidia_size_t j;
	for (j = 0; j < n; j++)
		D[j] = -1;
	//std::cout<<"NULLSPACE: "<<flag<<std::endl;
	if (flag == 0)
		nullspace_small();
	else
		nullspace_big();
}



void berl_matrix::nullspace_small()
{
	debug_handler("berl_matrix", "nullspace_small()");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	lidia_size_t i, j, k, l, pos;
	udigit t1, t2, p = p_small;
	udigit *x, *y, **M = M_small;
	double bpinv, pinv = 1.0/static_cast<double>(p);

	r = 0;
	l = 0;
	for (k = 0; k < n; k++) {
		if (verbose && k % 10 == 0) std::cerr << "+";
		pos = -1;
		for (i = l; i < n; i++)
			if (pos == -1 && M[i][k] != 0) pos = i;

		if (pos != -1) {
			// swap rows 'pos' and 'l'
			x = M[pos];
			M[pos] = M[l];
			M[l] = x;

			// multiply row 'l' with t1 = -1/M[l][k] mod p
			t1 = invert_mod(M[l][k], p);
			t1 = negate_mod(t1, p);
			M[l][k] = p - 1;
			bpinv = static_cast<double>(t1)*pinv;
			for (j = k + 1; j < n; j++)
//		M[l][j] = mul_mod_sp2(M[l][j], t1, p, bpinv);
				M[l][j] = multiply_mod(M[l][j], t1, p);

			// add multiples of row 'l' to the other rows
			for (i = l + 1; i < n; i++) {
				// M[i] = M[i] + M[l]*M[i,k]
				t1 = M[i][k];
				bpinv = static_cast<double>(t1)*pinv;
				if (t1 != 0) {
					M[i][k] = 0;
					x = &M[i][k+1];
					y = &M[l][k+1];
					for (j = k + 1; j < n; j++, x++, y++) {
//			t2 = mul_mod_sp2(*y, t1, p, bpinv);
						t2 = multiply_mod(*y, t1, p);
						*x = add_mod(*x, t2, p);
					}
				}
			}
			D[k] = l; // variable k is defined by row l
			l++;
		}
		else {
			r++; // finally, r will be the dimension of the kernel
		}
	}
}



// compare berl_matrix::nullspace_small() - it's nearly identical
void berl_matrix::nullspace_big()
{
	debug_handler("berl_matrix", "nullspace_big()");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	lidia_size_t i, j, k, l, pos;
	bigint t1, t2, p = p_big;
	bigint *x, *y, **M = M_big;

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
			// swap rows 'pos' and 'l'
			x = M[pos];
			M[pos] = M[l];
			M[l] = x;

			// multiply row 'l' with t1 = -1/M[l][k] mod p
			InvMod(t1, M[l][k], p);
			NegateMod(t1, t1, p);
			M[l][k] = -1;
			for (j = k + 1; j < n; j++)
				MulMod(M[l][j], M[l][j], t1, p);

			// add multiples of row 'l' to the other rows
			for (i = l + 1; i < n; i++) {
				// M[i] = M[i] + M[l]*M[i,k]
				t1.assign(M[i][k]); // this is already reduced
				if (!t1.is_zero()) {
					M[i][k].assign_zero();
					x = &M[i][k+1];
					y = &M[l][k+1];
					for (j = k + 1; j < n; j++, x++, y++) {
						multiply(t2, *y, t1);
						add(*x, *x, t2);
					}
				}
			}
			D[k] = l; // variable k is defined by row l
			l++;
		}
		else {
			r++;
		}

		//reduce next column of remaining rows
		if ((flag != 0) && (k != n - 1))
			for (j = l; j < n; j++)
				remainder(M[j][k + 1], M[j][k + 1], p);
	}
	//one could delete rows i with D[i] == -1, couldn't one ?
}



Fp_polynomial berl_matrix::random_basis_elt() const
	//output: g = random basis element of ker M
{
	debug_handler("berl_matrix", "random_basis_elt()");

	bigint t1, t2;
	lidia_size_t i, j, s;
	bigint p = p_big;
	Fp_polynomial g;
	g.set_modulus(p);
	for (j = n - 1; j >= 0; j--) {
		if (D[j] == -1)
			g.set_coefficient(randomize(p), j);
		else {
			i = D[j];

			// g[j] = sum_{s=j+1}^{n-1} g[s]*M[i,s]
			t1.assign_zero();
			for (s = j + 1; s < n; s++) {
				if (flag == 0) multiply(t2, g[s], M_small[i][s]);
				else           multiply(t2, g[s], M_big[i][s]);
				add(t1, t1, t2);
			}
			g.set_coefficient(t1, j); //reduction mod p is done automatically
		}
	}
	//leading zeros are removed by Fp_polynomial::set_coefficient
	return g;
}



//***********************************************************************
//
//	several methods for performing the last step of the
//	Berlekamp algorithm (i.e. really finding the factors)
//
//***********************************************************************

//
// Algorithm:	Compute gcd(v, g - c)
//		 for all already found factors v of f,
//		 for all c \in F_p (!!!) and
//		 for g \in ker(M) (i.e. g^p == g mod f)
//		until we have found all r irred. factors.
//		Actually, we'll succeed if we use (for g) all elements of a
//		basis of the kernel, but using random elements of the kernel is
//		also possible.
//
//		This step is of course only feasable for small p.
//

inline
static
bool berl_gcd_method(factorization< Fp_polynomial > &factors,
		     const berl_matrix &M, const Fp_poly_modulus &F)
{
	debug_handler("berlekamp.c", "berl_gcd_method(factorization< Fp_polynomial > &, berl_matrix&, Fp_poly_modulus&)");

	const Fp_polynomial &f = F.modulus();
	const bigint &p = f.modulus();
	if (p > 100) {
		lidia_error_handler("Fp_polynomial", "berl_gcd_method::"
				    "modulus too large");
		return false; // LC
	}

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;
	if (verbose) t.start("computing gcd's ...");

	lidia_size_t r = M.dim_nullspace();
	Fp_polynomial *container = new Fp_polynomial[r];
	memory_handler(container, "Fp_polynomial", "berl_gcd_method::"
		       "Error in memory allocation");

	container[0].assign(f);
	lidia_size_t i, j, len = 1, count = 0;
	Fp_polynomial g, h;

	do {
		if (verbose) std::cerr << "+";
		count++;

		g = M.random_basis_elt();

		// compute gcd(g+a, f_i) for all a \in \F_p and all already found
		// factors f_i
		for (i = 0; (i < p) && (len < r); i++) {
			for (j = 0; (j < len) && (len < r); j++) {
				gcd(h, container[j], g);
				if (!h.is_one() && h != container[j]) {
					container[len].assign(h);
					len++;
					divide(container[j], container[j], h);
				}
			}
			add(g, g, bigint(1));
		}
	} while (len < r && count < 100);

	if (len == r) {
		for (i = 0; i < r; i++)
			append_irred_factor(factors, container[i]);
	}
	else
		if (verbose) std::cerr << "failed." << std::endl;

	delete[] container;
	if (verbose) t.stop();
	return (len == r);
}



//
// Algorithm:	The idea is:
// 		If g is a random elt. of F_p[x]/(f), then
//		f = gcd(f, g)*gcd(f, g^{(p-1)/2} - 1)*gcd(f, g^{(p-1)/2} + 1),
//		and this factorization is nontrivial with a probability > 1/2.
//

inline
static
bool berl_dirpower_method(factorization< Fp_polynomial > &factors,
			  const berl_matrix &M, const Fp_poly_modulus &F)
{
	debug_handler("berlekamp.c", "berl_dirpower_method(factorization< Fp_polynomial > &, berl_matrix&, Fp_poly_modulus&)");

	const Fp_polynomial &f = F.modulus();
	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;
	if (verbose) t.start("direct powering method...");

	lidia_size_t r = M.dim_nullspace();
	Fp_polynomial *container = new Fp_polynomial[r];
	memory_handler(container, "Fp_polynomial", "berl_dirpower_method::"
		       "Error in memory allocation");

	container[0].assign(f);
	lidia_size_t i, len = 1, count = 0;;

	bigint p2;
	shift_right(p2, f.modulus(), 1); //p2 = (p-1)/2
	Fp_polynomial g, h, tmp;

	do {
		if (verbose) std::cerr << "+";
		count++;

		g = M.random_basis_elt(); // g = rand. elt. of F_p[x]/(f)
//	bigint c = randomize(p);
//	add(g, g, c);
		power(h, g, p2, F);
		add(h, h, bigint(1)); // h = g^{(p-1)/2} + 1

		for (i = 0; (i < len) && (len < r); i++) {
			Fp_polynomial &f_i = container[i];

			gcd(tmp, h, f_i);
			if (!tmp.is_one() && tmp != f_i) {
				divide(container[len], f_i, tmp);
				f_i.assign(tmp);
				len++;
			}
			if (len < r) {
				gcd(tmp, g, f_i);
				if (!tmp.is_one() && tmp != f_i) {
					divide(container[len], f_i, tmp);
					f_i.assign(tmp);
					len++;
				}
			}
		}
	} while (len < r && count < 1000);

	if (len == r) {
		for (i = 0; i < r; i++)
			append_irred_factor(factors, container[i]);
	}
	else
		if (verbose) std::cerr << "failed." << std::endl;

	if (verbose) t.stop();
	delete[] container;
	return (len == r);
}



//
// Algorithm:	Given a polynomial g \in ker(M) (i.e. g^p == g mod f),
//		we're interested in those c \in F_p, for which
//		gcd(f, g - c) is nontrivial.
//		If we compute a minimum polynomial h of g, then
//		h(f) == 0 mod f, and the roots of h are exactly those c
//		we are looking for
//

inline
static
bool berl_minpoly_method(factorization< Fp_polynomial > &factors,
			 berl_matrix &M, const Fp_poly_modulus &F)
{
	debug_handler("berlekamp.c", "berl_minpoly_method(factorization< Fp_polynomial > &, berl_matrix&, Fp_poly_modulus&)");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;
	if (verbose) t.start("computing splitting polynomial...");

	lidia_size_t r = M.dim_nullspace();
	lidia_size_t count = 0;
	const Fp_polynomial &f = F.modulus();
	Fp_polynomial g, h;
	base_vector< bigint > roots;

	do {
		if (verbose) std::cerr << "+";
		count++;

		g = M.random_basis_elt();
		min_poly(h, g, r, F);
	} while (h.degree() < r && count < 100);

	if (verbose) t.stop();

	if (count < 100) {
		M.kill();
		if (verbose) t.start("finding roots of splitting polynomial...");
		rec_find_roots(roots, h);
		if (verbose) t.stop();

		if (verbose) t.start("finding factors...");
		find_irred_factors(factors, f, g, roots);
		if (verbose) t.stop();
	}
	else
		if (verbose) std::cerr << "failed." << std::endl;

	return (count < 100);
}



//***********************************************************************
//
//			the main routine
//
//***********************************************************************

//
// Task:	Computes factorization of f = F.modulus()
//
// Conditions:	f is square-free, monic, deg(f) > 0 and f(0) != 0.
//		x_to_the_p = x^p mod f
//
// Algorithm:	Berlekamp
//

void sf_berlekamp_work(factorization< Fp_polynomial > &factors,
		       const Fp_polynomial &x_to_the_p, const Fp_poly_modulus &F)
{
	debug_handler("Fp_polynomial", "sf_berlekamp_work(base_vector< Fp_polynomial > &, Fp_polynomial&, Fp_poly_modulus&)");

	const Fp_polynomial &f = F.modulus();

	if (f.const_term().is_zero() || f.degree() <= 0) {
		lidia_error_handler("Fp_polynomial", "sf_berlekamp::input polynomial"
				    "must be monic with non-zero");
		return; // LC
	}

	factors.kill();

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;

	const bigint & p = f.modulus();
	lidia_size_t n = f.degree();

	Fp_polynomial g, h;

	lidia_size_t r;

	if (verbose) t.start("building matrix...");
	berl_matrix M(x_to_the_p, n, F); // compute Berlekamp-matrix
	if (verbose) t.stop();

	if (verbose) t.start("diagonalizing...");
	M.nullspace(); // compute nullspace of it
	if (verbose) t.stop();

	r = M.dim_nullspace();
	if (verbose) std::cerr << "number of factors = " << r << "\n";

	if (r == 1) {
		append_irred_factor(factors, f);
		return;
	}

	bool success = false;

//    if (p < 5)
//	success = berl_gcd_method(factors, M, F);
//    else

	if (p > 2*r)
		success = berl_minpoly_method(factors, M, F);

	if (!success)
		success = berl_dirpower_method(factors, M, F);

	if (!success) {
		lidia_error_handler_c("Fp_polynomial", "sf_berlekamp::couldn't factor "
				      "polynomial", std::cout << f << std::endl;);
		return; // LC
	}

	if (verbose) {
		std::cerr << "degrees:";
		lidia_size_t i;
		for (i = 0; i < factors.no_of_prime_components(); i++)
			std::cerr << " " << factors.prime_base(i).base().degree();
		std::cerr << "\n";
	}
}



//***********************************************************************
//
//			"interface"
//
//***********************************************************************

void
sf_berlekamp(factorization< Fp_polynomial > &factors, const Fp_polynomial & f)
{
	debug_handler("Fp_polynomial", "sf_berlekamp(factorization< Fp_polynomial > &, Fp_polynomial&)");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;
	Fp_polynomial g;
	Fp_poly_modulus F(f);

	if (verbose) t.start("computing x^p...");
	power_x(g, f.modulus(), F);
	if (verbose) t.stop();

	sf_berlekamp_work(factors, g, F);
}



void
berlekamp(factorization< Fp_polynomial > &factors, const Fp_polynomial & f)
{
	debug_handler("Fp_polynomial", "berlekamp(factorization< Fp_polynomial > &, Fp_polynomial&)");

	factor_generic(factors, f, &LiDIA::sf_berlekamp);
}



factorization< Fp_polynomial > berlekamp(const Fp_polynomial& f)
{
	debug_handler("Fp_polynomial", "berlekamp(Fp_polynomial&)");

	factorization< Fp_polynomial > F;
	berlekamp(F, f);
	return F;
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::berlekamp() const
{
	debug_handler("single_factor< Fp_polynomial >", "berlekamp()");

	factorization< Fp_polynomial > F;
	LiDIA::berlekamp(F, rep);
	return F;
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::sf_berlekamp() const
{
	debug_handler("single_factor< Fp_polynomial >", "sf_berlekamp()");

	factorization< Fp_polynomial > F;
	LiDIA::sf_berlekamp(F, rep);
	return F;
}



factorization< Fp_polynomial > sf_berlekamp(const Fp_polynomial& f)
{
	debug_handler("Fp_polynomial", "sf_berlekamp(Fp_polynomial&)");

	factorization< Fp_polynomial > F;
	sf_berlekamp(F, f);
	return F;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
