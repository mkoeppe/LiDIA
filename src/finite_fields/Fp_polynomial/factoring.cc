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
#include	"LiDIA/udigit.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//***********************************************************************
//
//		Algorithms for finding roots
//
//***********************************************************************

//
// Task:	appends the list of roots to x
//
// Conditions:	f is monic, and has deg(f) distinct roots.
//

void rec_find_roots(base_vector< bigint > & x, const Fp_polynomial& f)
{
	debug_handler("factoring.c", "find_roots(base_vector< bigint > &, Fp_polynomial&)");

	if (f.degree() == 0)
		return;

	const bigint &p = f.modulus();

	if (f.degree() == 1) {
		lidia_size_t k = x.size();

		if (x.capacity() < k+1)
			x.set_capacity(k+1);
		x.set_size(k+1);

		NegateMod(x[k], f.const_term(), p);
		return;
	}

	Fp_polynomial h;

	bigint r;
	bigint p1(p);
	p1.divide_by_2(); // p1 = (p-1)/2


	{
		Fp_poly_modulus F(f);
		do
		{
			r.assign(randomize(p));
			power_x_plus_a(h, r, p1, F);
			add(h, h, -1);
			gcd(h, h, f); // split f
		} while (h.degree() <= 0 || h.degree() == f.degree());
	}

	rec_find_roots(x, h); // recursively
	divide(h, f, h);
	rec_find_roots(x, h);
}



//
// Task:	Returns the list of roots of f (without multiplicities)
//
// Conditions:	always: f mustn't be the zero polynomial
//		if (flag == 0): no assumptions [default]
//		if (flag != 0): f has deg(f) distinct roots.
//

base_vector< bigint > find_roots(const Fp_polynomial& f, int flag)
{
	debug_handler("factoring.c", "find_roots(Fp_polynomial&)");

	if (f.is_zero()) {
		lidia_error_handler("factoring.c",
				    "find_roots(Fp_polynomial&)::input polynomial is zero");
		return base_vector< bigint > (); // LC
	}

	Fp_polynomial g(f);
	g.make_monic();

	base_vector< bigint > roots;

	if (flag != 0) {
		rec_find_roots(roots, g); // f has deg(f) distinct roots
	}
	else {
		Fp_polynomial x_to_the_p, x;
		x.set_modulus(g);
		x.assign_x();
		Fp_poly_modulus F(g);

		// this is the costly step we want to avoid if flag != 0
		power_x(x_to_the_p, g.modulus(), g);

		gcd(g, g, x_to_the_p - x); // split off linear factors

		factorization< Fp_polynomial > sqfr;
		square_free_decomp(sqfr, g); // decompose acc. to multiplicities

		// any component of sqfr is a product
		// of distinct linear factors
		lidia_size_t i;
		for (i = 0; i < sqfr.no_of_composite_components(); i++)
			rec_find_roots(roots, sqfr.composite_base(i).base());
	}
	return roots;
}



//
// Task:        Finds a root of ff.
//
// Conditions:  f splits into distinct linear factors.
//

bigint find_root(const Fp_polynomial& ff)
{
	debug_handler("factoring.cc", "find_root(Fp_polynomial&)");

	Fp_poly_modulus F;
	Fp_polynomial h, f(ff);

	bigint r;

	bigint minus_one, root;
	minus_one.assign_one();
	minus_one.negate();

	const bigint &p = ff.modulus();

	bigint p1(p);
	p1.divide_by_2(); // p1 = (p-1)/2


	while (f.degree() > 1) {
		F.build(f);
		r.assign(randomize(p));
		power_x_plus_a(h, r, p1, F);
		add(h, h, -1);
		gcd(h, h, f); // split f
		if (h.degree() > 0 && h.degree() < f.degree()) {
			if (h.degree() > f.degree()/2)
				divide(f, f, h);
			else
				f.assign(h);
		}
	}

	NegateMod(root, f.const_term(), p);
	return root;
}

//
// Task:        Finds an irreducible factor of degree d of ff.
//
// Conditions:  ff splits into distinct irreducible factors of degree d.
//
Fp_polynomial
find_factor( const Fp_polynomial & ff, lidia_size_t d )
{
    debug_handler( "factoring.cc", "find_factor( Fp_polynomial& , degree )" );

    Fp_poly_modulus F;
    Fp_polynomial h, f(ff);

    bigint r;

    bigint minus_one, root;
    minus_one.assign_one();
    minus_one.negate();

    const bigint &p = ff.modulus();

    bigint p1;
	power( p1, p, d );
    p1.divide_by_2();				// p1 = (p^d-1)/2
        
    while( f.degree() > d )
    {
		F.build( f );
		r.assign( randomize( p ) );
		power_x_plus_a( h, r, p1, F );
		add(h, h, -1);
		gcd(h, h, f);				// split f
		if (h.degree() > 0 && h.degree() < f.degree())
		{
			if ( h.degree() > f.degree()/2 )
				divide(f, f, h);
			else
				f.assign( h );
		}
    }
    
    return( f );
}

//
// Task:        Finds a factor of ff of degree at most d-1 of ff.
//
// Conditions:  ff splits into distinct linear factors.
//
Fp_polynomial
find_factor_degree_d( const Fp_polynomial & ff, lidia_size_t d )
{
    debug_handler( "factoring.cc", "find_root( Fp_polynomial& )" );

    Fp_poly_modulus F;
    Fp_polynomial h, f(ff);

    bigint r;

    bigint minus_one, root;
    minus_one.assign_one();
    minus_one.negate();

    const bigint &p = ff.modulus();

    bigint p1(p);
    p1.divide_by_2();				// p1 = (p-1)/2
    
    while( f.degree() >= 1 )
    {
		F.build(f);
		r.assign(randomize(p));
		power_x_plus_a(h, r, p1, F);
		add(h, h, -1);
		gcd(h, h, f);	// split f
		if (h.degree() > 0 && h.degree() < f.degree())
		{
			if (h.degree() > f.degree()/2)
				divide(f, f, h);
			else
				f.assign( h );
		}
		if( f.degree() < d )
			return( f );
    }
    return f;
}



//***********************************************************************
//
//		    trace_map, power_compose
//
//***********************************************************************

//
// Literature:	J. von~zur~Gathen and V. Shoup:
//		Computing Frobenius Maps and factoring polynomials,
//		Computational Complexity, 2:187-224, 1992
//		Algorithm 5.2.
//
// Task:	computes w = a+a^q+...+^{q^{d-1}} mod f    (f = F.modulus())
//
// Conditions:  d >= 0,
//		b = x^q mod f,
//		q a power of p
//
// Algorithm:	trace map
//		Space allocation can be controlled via ComposeBound
//
//	The main idea is the following:
//	Let	w_m = sum_{0<=i<m} a^{q^i}
//	and 	z_m = x^{q^m}.
//
//	Then	w_{2m} = w_m + w_m(z_m)
//	and	z_{2m} = z_m(z_m)
//

void trace_map(Fp_polynomial& w, const Fp_polynomial& a, lidia_size_t d,
	       const Fp_poly_modulus& F, const Fp_polynomial& b)
{
	debug_handler("factoring.c", "trace_map(Fp_polynomial&, Fp_polynomial&, lidia_size_t, Fp_poly_modulus&, Fp_polynomial&)");

	a.comp_modulus(b, "trace_map");
	a.comp_modulus(F.modulus(), "trace_map");

	Fp_polynomial y, z, t;

	z.assign(b);
	y.assign(a);
	w.set_modulus(b); // assigns zero

	// Let d = sum_{0<=i<=l} d_i*2^i.
	// Set i=0.

	while (d) {
		if (d == 1) {
			if (w.is_zero())
				w.assign(y); // w = w(z) + y = y
			else {
				compose(w, w, z, F);
				add(w, w, y); // w = w(z) + y
			}
			// we don't need to compute new values for y and z
			// because we're done
		}
		else {
			if ((d & 1) == 0) {
				compose2(z, t, z, y, z, F); // z = z(z)
				add(y, t, y); // y = y(z) + y
			}
			else {
				if (w.is_zero()) {
					w.assign(y); // w = w(z) + y = y
					compose2(z, t, z, y, z, F); // z = z(z)
					add(y, t, y); // y = y(z) + y
				}
				else {
					compose3(z, t, w, z, y, w, z, F); // z = z(z)
					add(w, w, y); // w = w(z) + y = y
					add(y, t, y); // y = y(z) + y
				}
			}
		}
		// Now, we have checked bit i of (the original) d.
		//	z = x^{q^{2^i}}
		//	y = sum_{1<=j<=2^i} a^{q^j}
		//	w = sum_{0<=j<=D_i} a^{q^j}, where D_i = sum_{0<=j<=i} d_j*2^j
		// Set i = i+1

		d = d >> 1;
	}
}



//
// Literature:	J. von~zur~Gathen and V. Shoup:
//		Computing Frobenius Maps and factoring polynomials,
//		Computational Complexity, 2:187-224, 1992
//		Algorithm 5.2.
//
// Task:	computes y = x^{q^d} mod f	(f = F.modulus())
//
// Conditions:  d >= 0,
//		b = x^q mod f,
//		q a power of p
//
// Algorithm:	uses modular composition
//		Space allocation can be controlled via ComposeBound
//

void power_compose(Fp_polynomial& y, const Fp_polynomial& b, lidia_size_t d, const Fp_poly_modulus& F)
{
	debug_handler("factoring.c", "power_compose(Fp_polynomial&, Fp_polynomial&, lidia_size_t, Fp_poly_modulus&)");

	b.comp_modulus(F.modulus(), "power_compose");

	Fp_polynomial z;
	z.set_max_degree(F.modulus().degree() - 1);
	lidia_size_t sw;

	z.assign(b);

	y.set_modulus(b);
	y.assign_x();

	// Let d = sum_{0<=i<=l} d_i*2^i.
	// Set i=0.

	while (d) {
		sw = 0;

		if (d > 1) sw = 2;
		if (d & 1) {
			if (y.is_x())				// very simple compos.:
				y.assign(z); // y = y(z)
			else
				sw = sw | 1;
		}

		switch (sw) {
		case 0:		// this should never happen (only poss. if d = 0)
			break;

		case 1:		// this will be our last step
			compose(y, y, z, F); // y = y(z)
			break;

		case 2:
			compose(z, z, z, F); // z = z(z)
			break;

		case 3:
			compose2(y, z, y, z, z, F); // y = y(z), z = z(z)
			break;
		}
		// Now, we have checked bit i of (the original) d.
		//	z = x^{q^{2^i}}
		//	y = x^{q^{D_i}}, where D_i = sum_{0<=j<=i} d_j*2^j
		// Set i = i+1.

		d = d >> 1;
	}
}



//***********************************************************************
//
//		Algorithms for irreducibility testing:
//	    prob_irred_test, iter_irred_test, det_irred_test
//
//***********************************************************************

//
// Literature:	Shoup, J. Symbolic Comp. 17:371-391, 1994
//
// Task:	Performs a fast, probabilistic irreduciblity test.
//		The test can err only if f is reducible, and the
//		error probability is bounded by p^{-iter}.
//		(p = f.modulus())
//
// Conditions:	Works for any p.
//

bool prob_irred_test(const Fp_polynomial& f, lidia_size_t iter)
{
	debug_handler("factoring.c", "prob_irred_test(Fp_polynomial&, lidia_size_t)");

	lidia_size_t n = f.degree();

	if (n <= 0)
		return false;
	if (n == 1)
		return true;

	const bigint& p = f.modulus();

	Fp_poly_modulus F(f);

	Fp_polynomial b, r, s;

	power_x(b, p, F);

	lidia_size_t i;
	for (i = 0; i < iter; i++) {
		randomize(r, p, n-1); // r = random element of F_p[X]/(f)
		trace_map(s, r, n, F, b); // s = "Trace"(r)

		if (s.degree() > 0)		// if f is irreducible, then F_p[X]/(f)
			return false; // is a field, and s lies in F_p
		// <==> deg(s) <= 0
	}

	if (p >= n)
		return true;

	lidia_size_t pp; // p is now < n
	p.sizetify(pp); // n is of type lidia_size_t

	if (n % pp != 0)
		return true; // f is definitely not a pth power

	power_compose(s, b, n/pp, F);
	return !s.is_x();
}



static
bool irred_base_case(const Fp_polynomial& h, lidia_size_t q, int a,
		     const Fp_poly_modulus& F)
{
	debug_handler("factoring.c", "irred_base_case(Fp_polynomial&, lidia_size_t, int, Fp_poly_modulus&)");

	lidia_size_t e;
	Fp_polynomial x, s, d;

	e = static_cast<lidia_size_t>(power(static_cast<udigit>(q), static_cast<unsigned int>(a-1)));
	power_compose(s, h, e, F);

	x.set_modulus(h);
	x.assign_x();
	subtract(s, s, x);
	gcd(d, F.modulus(), s);
	return d.is_one();
}



static
int compute_split(int lo, int hi, const fac_vec& fvec)
{
	debug_handler("factoring.c", "compute_split(int, int, fac_vec&)");

	int mid, i;
	double total, sum;

	total = 0;
	for (i = lo; i <= hi; i++)
		total = total + fvec[i].len;

	mid = lo-1;
	sum = 0;
	while (sum < total/2) {
		mid++;
		sum = sum + fvec[mid].len;
	}

	if (mid == hi || (mid != lo && 2*sum > total + fvec[mid].len))
		mid--;

	return mid;
}



static
void tandem_power_compose(Fp_polynomial& y1, Fp_polynomial& y2,
			  const Fp_polynomial& h, lidia_size_t q1, lidia_size_t q2,
			  const Fp_poly_modulus& F)
{
	debug_handler("factoring.c", "tandem_power_compose(Fp_polynomial&, Fp_polynomial&, Fp_polynomial& h, lidia_size_t, lidia_size_t, Fp_poly_modulus&)");

	Fp_polynomial z;
	z.set_max_degree(F.modulus().degree() - 1);
	lidia_size_t sw;

	z.assign(h);

	y1.set_modulus(h);
	y2.set_modulus(h);
	y1.assign_x();
	y2.assign_x();

	while (q1 || q2) {
		sw = 0;

		if (q1 > 1 || q2 > 1)
			sw = 4;

		if (q1 & 1) {
			if (y1.is_x())
				y1.assign(z);
			else
				sw = sw | 2;
		}

		if (q2 & 1) {
			if (y2.is_x())
				y2.assign(z);
			else
				sw = sw | 1;
		}

		switch (sw) {
		case 0:
			break;

		case 1:
			compose(y2, y2, z, F);
			break;

		case 2:
			compose(y1, y1, z, F);
			break;

		case 3:
			compose2(y1, y2, y1, y2, z, F);
			break;

		case 4:
			compose(z, z, z, F);
			break;

		case 5:
			compose2(z, y2, z, y2, z, F);
			break;

		case 6:
			compose2(z, y1, z, y1, z, F);
			break;

		case 7:
			compose3(z, y1, y2, z, y1, y2, z, F);
			break;
		}

		q1 = q1 >> 1;
		q2 = q2 >> 1;
	}
}



static
bool rec_irred_test(int lo, int hi, const Fp_polynomial& h,
		    const Fp_poly_modulus& F, const fac_vec& fvec)
{
	debug_handler("factoring.c", "rec_irred_test(int, int, Fp_polynomial&, Fp_poly_modulus&, fac_vec&)");

	lidia_size_t q1, q2;
	int mid;
	Fp_polynomial h1, h2;

	if (h.is_x()) return false;

	if (lo == hi) {
		return irred_base_case(h, fvec[lo].q, fvec[lo].a, F);
	}

	mid = compute_split(lo, hi, fvec);

	q1 = fvec.prod(lo, mid);
	q2 = fvec.prod(mid+1, hi);

	tandem_power_compose(h1, h2, h, q1, q2, F);
	return 	rec_irred_test(lo, mid, h2, F, fvec)
		&& rec_irred_test(mid+1, hi, h1, F, fvec);
}



//
// Literature:	Shoup, J. Symbolic Comp. 17:371-391, 1994
//
// Task:	Performs a recursive deterministic irreducibility test;
//		fast in the worst-case (when input is irreducible).
//

bool det_irred_test(const Fp_polynomial& f)
{
	debug_handler("factoring.c", "det_irred_test(Fp_polynomial&)");

	if (f.degree() <= 0)
		return false;
	if (f.degree() == 1)
		return true;

	Fp_poly_modulus F(f);

	Fp_polynomial h;
	const bigint &p = f.modulus();

	power_x(h, p, F);

	Fp_polynomial s;
	lidia_size_t n = f.degree();

	power_compose(s, h, n, F);
	if (!s.is_x())
		return false;

	fac_vec fvec(n);

	int NumFactors = fvec.number_of_factors();

	return rec_irred_test(0, NumFactors-1, h, F, fvec);
}



//
// Literature:	Shoup, J. Symbolic Comp. 17:371-391, 1994
//
// Task:	Performs an iterative deterministic irreducibility test,
//		based on DDF.
//		Fast on average (when f has a small factor).
//

bool iter_irred_test(const Fp_polynomial& f)
{
	debug_handler("factoring.c", "iter_irred_test(Fp_polynomial&)");

	const int MaxLimit = 8;

	if (f.degree() <= 0)
		return false;
	if (f.degree() == 1)
		return true;

	Fp_poly_modulus F(f);

	Fp_polynomial h;
	const bigint &p = f.modulus();

	power_x(h, p, F);

	lidia_size_t CompTableSize = 2*square_root(f.degree());

	poly_argument H;

	H.build(h, F, CompTableSize);

	lidia_size_t i, d, limit, old_n;
	Fp_polynomial g, x, t, prod;

	x.set_modulus(f);
	x.assign_x();

	i = 0;
	g = h;
	d = 1;
	limit = 1;

	prod.set_modulus(f);
	prod.assign_one();


	while (2*d <= f.degree()) {
		old_n = f.degree();
		subtract(t, g, x);
		multiply(prod, prod, t, F); //Fp_poly_modulus
		i++;
		if (i == limit) {
			gcd(t, f, prod);
			if (!t.is_one())
				return false;

			prod.assign_one();
			limit = comparator< lidia_size_t >::min(MaxLimit, 2 * limit);
			i = 0;
		}

		d = d + 1;
		if (2*d <= f.degree()) {
			H.compose(g, g, F);
		}
	}

	if (i > 0) {
		gcd(t, f, prod);
		if (!t.is_one())
			return false;
	}

	return true;
}



//***********************************************************************
//
//	    Algorithms for constructing irreducible polynomials
//		    build_random_irred, build_irred
//
//***********************************************************************

static
void multiply_by_x_plus_y(base_vector< Fp_polynomial > & h,
			  const Fp_polynomial& f, const Fp_polynomial& g)
// h represents the bivariate polynomial h[0] + h[1]*y + ... + h[n-1]*y^k,
// where the h[i]'s are polynomials in x, each of degree < deg(f),
// and k < deg(g).
// h is replaced by the bivariate polynomial h*(x+y) (mod f(x), g(y)).

{
	debug_handler("factoring.c", "multiply_by_x_plus_y(base_vector< Fp_polynomial > &, Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t n = g.degree();
	lidia_size_t k = h.size()-1;
	lidia_size_t i;

	if (k < 0)
		return;

	if (k < n-1) {
		if (h.capacity() < k+2)
			h.set_capacity(k+2);
		h.set_size(k+2);

		h[k+1].assign(h[k]);
		for (i = k; i >= 1; i--) {
			multiply_by_x_mod(h[i], h[i], f);
			add(h[i], h[i], h[i-1]);
		}
		multiply_by_x_mod(h[0], h[0], f);
	}
	else {
		Fp_polynomial b, t;

		b.assign(h[n-1]);
		for (i = n-1; i >= 1; i--) {
			multiply_by_scalar(t, b, g[i]);
			multiply_by_x_mod(h[i], h[i], f);
			add(h[i], h[i], h[i-1]);
			subtract(h[i], h[i], t);
		}
		multiply_by_scalar(t, b, g[0]);
		multiply_by_x_mod(h[0], h[0], f);
		subtract(h[0], h[0], t);
	}

	// normalize

	k = h.size()-1;
	while (k >= 0 && h[k].is_zero())
		k--;
	if (h.capacity() < k+1)
		h.set_capacity(k+1);
	h.set_size(k+1);

}



static
void irred_combine(Fp_polynomial& x, const Fp_polynomial& f,
		   const Fp_polynomial& g)
{
	debug_handler("factoring.c", "irred_combine(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t df = f.degree();
	lidia_size_t dg = g.degree();
	lidia_size_t m = df*dg;

	if (df < dg) {
		irred_combine(x, g, f);
		return;
	}

	// f.degree() >= deg(g)...not necessary, but maybe a little more
	// time & space efficient

	Fp_polynomial a, b;
	bigint t;
	const bigint &p = g.modulus();

	base_vector< Fp_polynomial > h(dg, dg);

	lidia_size_t i;
	for (i = 0; i < dg; i++)
		h[i].set_max_degree(df-1);

	if (h.capacity() == 0)
		h.set_capacity(1);
	h.set_size(1);
	h[0].set_modulus(f);
	h[0].assign_one();

	a.set_modulus(f);
	b.set_modulus(f);

	poly_matrix M(p);

	a.set_max_degree(2*m-1);
	for (i = 0; i < 2*m; i++) {
		a.set_coefficient(h[0].const_term(), 2*m-1-i);
		if (i < 2*m-1)
			multiply_by_x_plus_y(h, f, g);
	}

	b.set_coefficient(2*m);

	M.half_gcd(b, a, m+1);

	// make monic
	InvMod(t, (M(1, 1)).lead_coeff(), p);
	multiply_by_scalar(x, M(1, 1), t);
}



static
void build_prime_power_irred(Fp_polynomial& f, const bigint & p,
			     lidia_size_t q, int e)
{
	debug_handler("factoring.c", "build_prime_power_irred(Fp_polynomial&, bigint&, lidia_size_t, int)");

	lidia_size_t n = static_cast<lidia_size_t>(power(static_cast<unsigned int>(q), static_cast<unsigned int>(e)));

	do {
		randomize(f, p, n-1);
		f.set_coefficient(n);
	} while (!iter_irred_test(f));
}



//
// Literature:	Shoup, J. Symbolic Comp. 17:371-391, 1994 (?)
//
// Task:	Build a monic irreducible poly of degree n (over Z/pZ).
//

void build_irred(Fp_polynomial& f, const bigint& p, lidia_size_t n)
{
	debug_handler("factoring.c", "build_irred(Fp_polynomial&, bigint&, lidia_size_t)");

	if (n <= 0) {
		lidia_error_handler("factoring.c", "build_irred(Fp_polynomial&, "
				    "bigint&, lidia_size_t)::n must be positive");
		return; // LC
	}

	if (n == 1) {
		f.set_modulus(p);
		f.assign_x();
		return;
	}

	fac_vec fvec(n);

	int i, NumFactors = fvec.number_of_factors();

	f.set_max_degree(n);


	build_prime_power_irred(f, p, fvec[0].q, fvec[0].a);

	Fp_polynomial g;
	g.set_max_degree(n);

	for (i = 1; i < NumFactors; i++) {
		build_prime_power_irred(g, p, fvec[i].q, fvec[i].a);
		irred_combine(f, f, g);
	}
}



//
// Literature:	Shoup, J. Symbolic Comp. 17:371-391, 1994 (?)
//
// Task:	constructs a random monic irreducible polynomial f of degree
//		g.degree()
//
// Conditions:	g is a monic irreducible polynomial.
//

void build_random_irred(Fp_polynomial& f, const Fp_polynomial& g)
{
	debug_handler("factoring.c", "build_random_irred(Fp_polynomial&, Fp_polynomial&)");

	const bigint &p = g.modulus();
	Fp_poly_modulus G(g);
	Fp_polynomial h, ff;

	do
	{
		randomize(h, p, g.degree()-1);
		irred_poly(ff, h, g.degree(), G);
	} while (ff.degree() < g.degree());

	f.assign(ff);
}



//***********************************************************************
//
//  Algorithms which are useful in counting points on elliptic curves:
//		compute_degree, prob_compute_degree
//
//***********************************************************************

static
lidia_size_t base_case(const Fp_polynomial& h, lidia_size_t q, int a,
		       const Fp_poly_modulus& F)
{
	debug_handler("factoring.c", "base_case(Fp_polynomial&, lidia_size_t, int, Fp_poly_modulus&)");

	lidia_size_t b;
	Fp_polynomial lh;
	lh.set_max_degree(F.modulus().degree() - 1);

	lh.assign(h);
	b = 0;
	while (b < a-1 && !lh.is_x()) {
		b++;
		power_compose(lh, lh, q, F);
	}

	if (!lh.is_x())
		b++;

	return b;
}



static
void rec_compute_degree(int lo, int hi, const Fp_polynomial& h,
			const Fp_poly_modulus& F, fac_vec& fvec)
{
	debug_handler("factoring.c", "rec_compute_degree(int, int, Fp_polynomial&, Fp_poly_modulus&, fac_vec&)");

	int mid;
	lidia_size_t q1, q2;
	Fp_polynomial h1, h2;

	if (h.is_x()) {
		fvec.clear(lo, hi);
		return;
	}

	if (lo == hi) {
		fvec[lo].b = base_case(h, fvec[lo].q, fvec[lo].a, F);
		return;
	}

	mid = compute_split(lo, hi, fvec);

	q1 = fvec.prod(lo, mid);
	q2 = fvec.prod(mid+1, hi);

	tandem_power_compose(h1, h2, h, q1, q2, F);
	rec_compute_degree(lo, mid, h2, F, fvec);
	rec_compute_degree(mid+1, hi, h1, F, fvec);
}



//
// Task:	The common degree of the irreducible factors of f is computed.
//		This routine is useful in counting points on elliptic curves.
//
// Conditions:	f = F.modulus() is assumed to be an "equal degree" polynomial,
//		h = x^p mod f,
//              d is multiple of common degree, if d == -1, such information
//              not known
//

lidia_size_t compute_degree(const Fp_polynomial& h, const Fp_poly_modulus& F,
                            lidia_size_t d)
{
  debug_handler("factoring.c", "compute_degree(Fp_polynomial&, Fp_poly_modulus&)");
  
  h.comp_modulus(F.modulus(), "compute_degree");
  
  if (h.is_x())
    return 1;
  
  lidia_size_t res;
  
  if (d == -1)
    d = F.modulus().degree();
  
  fac_vec fvec(d);
  
  int i, NumFactors = fvec.number_of_factors();
  
  rec_compute_degree(0, NumFactors-1, h, F, fvec);
  
  res = 1;
  
  for (i = 0; i < NumFactors; i++)
    {
      res = res * 
	static_cast<lidia_size_t>(power(static_cast<udigit>(fvec[i].q),
					static_cast<unsigned int>(fvec[i].b)));
    }
  
  return res;
}



//
// Task:	Same as above, but uses a slightly faster probabilistic
//		algorithm.
//		The return value may be 0 or may be too big, but for large p
//		(relative to n), this happens with very low probability.
//
// Conditions:	Same as above.
//

lidia_size_t prob_compute_degree(const Fp_polynomial& h, const Fp_poly_modulus& F)
{
	debug_handler("factoring.c", "prob_compute_degree(Fp_polynomial&, Fp_poly_modulus&)");

	h.comp_modulus(F.modulus(), "prob_compute_degree");

	if (h.is_x())
		return 1;

	lidia_size_t n = F.modulus().degree();
	const bigint &p = h.modulus();

	Fp_polynomial P1, P2, P3;

	randomize(P1, p, n-1);
	trace_map(P2, P1, n, F, h);
	min_poly(P3, P2, n/2, F);

	lidia_size_t r = P3.degree();

	if (r <= 0 || n % r != 0)
		return 0;
	else
		return n/r;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
