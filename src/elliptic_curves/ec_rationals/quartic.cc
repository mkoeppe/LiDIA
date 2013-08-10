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
//	Author	: Nigel Smart
//                Adaption of John Cremona's code
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/quartic.h"
#include	"LiDIA/elliptic_curves/ec_arith.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



quartic::quartic()
{
	roots = new bigcomplex[4];
}



quartic::quartic(const bigint& qa, const bigint& qb, const bigint& qc,
                 const bigint& qd, const bigint& qe)
	: a(qa),
	  b(qb),
	  c(qc),
	  d(qd),
	  e(qe)
{
	ii = bigint(12)*a*e - bigint(3)*b*d + c*c;
	jj = bigint(72)*a*c*e + bigint(9)*b*c*d - bigint(27)*a*d*d - bigint(27)*b*b*e - bigint(2)*c*c*c;
	power(disc, ii, 3);
	disc = bigint(4)*disc-jj*jj;

	roots = solve_real_quartic(bigfloat(a), bigfloat(b), bigfloat(c),
				   bigfloat(d), bigfloat(e));

	long i, nrr = 0;
	for (i = 0; i < 4; i++) {
		if (roots[i].is_approx_real()) {
			nrr++;
		}
	}
	type = nrr+1;
}



quartic::quartic(const quartic& q)
	: a(q.a),
	  b(q.b),
	  c(q.c),
	  d(q.d),
	  e(q.e),
	  ii(q.ii),
	  jj(q.jj),
	  disc(q.disc),
	  type(q.type)
{
	roots = new bigcomplex[4];
	for (int i = 0; i < 4; i++)
		roots[i] = q.roots[i];
}



void quartic::assign(const bigint& qa, const bigint& qb, const bigint& qc,
                     const bigint& qd, const bigint& qe)
{
	a = qa;
	b = qb;
	c = qc;
	d = qd;
	e = qe;
	ii = bigint(12)*a*e - bigint(3)*b*d + c*c;
	jj = bigint(72)*a*c*e + bigint(9)*b*c*d - bigint(27)*a*d*d - bigint(27)*b*b*e - bigint(2)*c*c*c;
	power(disc, ii, 3);
	disc = 4*disc-jj*jj;

	delete[] roots;
	roots = solve_real_quartic(bigfloat(a), bigfloat(b), bigfloat(c),
				   bigfloat(d), bigfloat(e));

	long i, nrr = 0;
	for (i = 0; i < 4; i++) {
		if (roots[i].is_approx_real()) {
			nrr++;
		}
	}
	type = nrr+1;
}



quartic & quartic::operator = (const quartic& q)
{
	a = q.a;
	b = q.b;
	c = q.c;
	d = q.d;
	e = q.e;
	for (int i = 0; i < 4; i++)
		roots[i] = q.roots[i];
	type = q.type;
	ii = q.ii;
	jj = q.jj;
	disc = q.disc;

	return *this;
}



void quartic::assign(const quartic& q)
{
	a = q.a;
	b = q.b;
	c = q.c;
	d = q.d;
	e = q.e;
	for (int i = 0; i < 4; i++)
		roots[i] = q.roots[i];
	type = q.type;
	ii = q.ii;
	jj = q.jj;
	disc = q.disc;
}



int quartic::trivial() const // Checks for a rational root
{
	int i, found;
	bigint num;
	bigfloat realroot;
	bigint ac = a*c, a2d = a*a*d, a3e = a*a*a*e;
	bigfloat ra = bigfloat(a);

	for (i = 0, found = 0; i < 4 && ! found; i++) {
		if (roots[i].is_approx_real()) {
			realroot = roots[i].real();
			(ra*realroot).bigintify(num);
			found = (((((num+b)*num+ac)*num+a2d)*num+a3e) == 0);
		}
	}
	return(found);
}



// find the number of roots of aX^4 + bX^3 + cX^2 + dX + e = 0 (mod p)
// except 4 is returned as 3 so result is 0,1,2 or 3.

// This last bit isn't exactly very LiDIA like,
// may need to change later
long quartic::no_roots_mod(long p) const
{
	long ap; best_remainder(ap, a, p);
	long bp; best_remainder(bp, b, p);
	long cp; best_remainder(cp, c, p);
	long dp; best_remainder(dp, d, p);
	long ep; best_remainder(ep, e, p);

	long nroots = (ap == 0); // must count infinity as a root!
	for (long i = 0; (i < p) && (nroots < 3); i++) {
		// The following line could overflow.....
		long temp = ((((ap*i+bp)*i + cp)*i + dp)*i + ep);
		if ((temp%p) == 0) {
			nroots++;
		}
	}
	if (nroots == 4)
		return 3;
	return nroots;
}



std::ostream& operator << (std::ostream& s, const quartic& q)
{
	s << "(" << q.a << ", " << q.b << ", " << q.c << ", " << q.d << ", " << q.e << ")" << std::flush;
	return s;
}



//
// BSD - CHECK FOR LOCAL SOLUBILITY (NOT QP-SOL HERE ONLY ZP-SOL)
//

// tests if aa is a p-adic square
int psquare(const bigint& aa, const bigint& p)
{
	if (aa.is_zero())
		return 1;
	long v = valuation(p, aa);
	if ((v%2) != 0)
		return 0;
	bigint a(aa);
	while (v--)
		a /= p;
	if (p == 2) {
		return ((a-bigint(1))%8).is_zero();
	}
	else {
		return jacobi(a, p) == 1;
	}
}



// returns -1 for insoluble, 0 for undecided, +1 for soluble --- odd p
int quartic::lemma6(const bigint& p, int nu, const bigint& x)
{
	bigint gx = (((a*x+b)*x+c)*x+d)*x+e;

	if (psquare(gx, p))
		return +1;

	bigint gdashx = ((4*a*x+3*b)*x+2*c)*x+d;
	long lambda = valuation(p, gx);
	long mu = valuation(p, gdashx);

	if ((lambda-mu >= nu) && (nu > mu))
		return +1;
	if ((lambda >= 2*nu) && (mu >= nu))
		return 0;
	return -1;
}



// returns -1 for insoluble, 0 for undecided, +1 for soluble --- p=2
int quartic::lemma7(const bigint& p, int nu, const bigint& x)
{
	bigint gx = (((a*x+b)*x+c)*x+d)*x+e;
	if (psquare(gx, p)) {
	    return +1;
	}
	// gx != 0, because psquare(0, p) == true.
	long lambda = valuation(p, gx);  

	// mu == -1 iff valuation(p, gdashx) == infinity 
	bigint gdashx = ((4*a*x+3*b)*x+2*c)*x+d;
	long mu = gdashx.is_zero() ? -1L : valuation(p, gdashx);

	bigint oddgx = gx;
	if (oddgx.is_zero())
		oddgx = 1;
	else {
		while (oddgx.is_even())
			oddgx /= 2;
	}

	int odd4 = ((oddgx-1)%4).is_zero();

	if (mu > -1) {
	    // mu = valuation(p, gdashx) \in N
	    if ((lambda-mu >= nu) && (nu > mu))
		return +1;
	    if ((nu > mu) && (lambda == mu+nu-1) && (lambda%2) == 0)
		return +1;
	    if ((nu > mu) && (lambda == mu+nu-2) && (lambda%2) == 0 && odd4)
		return +1;
	    if ((mu >= nu) && (lambda >= 2*nu))
		return 0;
	    if ((mu >= nu) && (lambda == 2*nu-2) && odd4)
		return 0;
	}
	else {
	    // mu = valuation(p, gdashx) = +infinity
	    if (lambda >= 2*nu)
		return 0;
	    if ((lambda == 2*nu-2) && odd4)
		return 0;
	}
	return -1;
}



// Checks for solublility in Zp with x=x0 (mod p^nu)
// Fully recursive (depth-first) version
int quartic::zp_soluble(const bigint& p, const bigint& x0, long nu)
{
	long result = (p == 2) ? lemma7(p, nu, x0) : lemma6(p, nu, x0);
	if (result == +1)
		return 1;
	if (result == -1)
		return 0;
//else result==0, so refine to look modulo p^(nu+1):
	bigint x = x0, pnu;
	power(pnu, p, nu);
	for (bigint i = 0; i < p; i++, x += pnu) {
		if (zp_soluble(p, x, nu+1))
			return 1;
	}
	return 0;
}



int quartic::qp_soluble(const bigint& p)
{
	if (p > 10) {
		return samir_qp_soluble(p);
	}
	else {
		if (zp_soluble(p, 0, 0))
			return 1;
		else
		{
			quartic q(e, d, c, b, a);
			return q.zp_soluble(p, 0, 1);
		}
	}
} // end of qpsoluble


int quartic::locally_soluble(const sort_vector< bigint > & plist)
{
	if (a.is_lt_zero() && type == 1) {
		// Not soluble at infinity
		return 0;
	}
	for (long i = 0; i < plist.size(); i++) {
		if (!qp_soluble(plist[i])) {
			return 0;
		}
	}
	return 1;
}



// The next two functions implement Siksek's local solubility
// test for "large" primes p.
// This runs in conjectured polynomial time and was explained
// in the paper on higher descents by
//	Merriman, Siksek and Smart

int samir_local_sol(const bigint& p, bigint *c);



int quartic::samir_qp_soluble(const bigint& p)
{
	bigint *cc = new bigint[5];

	cc[0] = a;
	cc[1] = b;
	cc[2] = c;
	cc[3] = d;
	cc[4] = e;

	int res = samir_local_sol(p, cc);

	if (!res) {
		swap(cc[0], cc[4]);
		swap(cc[1], cc[3]);
		res = samir_local_sol(p, cc);
	}
	delete[] cc;
	return res;
}



factorization< Fp_polynomial > fact_c(const bigint& p, bigint *c)
{
	Fp_polynomial f;

	f.set_modulus(p);
	for (long i = 0; i < 5; i++) {
		f.set_coefficient(c[i], i);
	}
	factorization< Fp_polynomial > fact_f = factor(f);
	//berlekamp(fact_f,f);
	return fact_f;
}



int samir_local_sol(const bigint& p, bigint *c)
{
	bigint r[2], t;
	long e, fl, i;
	bigint p2 = p*p;
	int xxx = 1;

	for (i = 0; i < 5 && xxx; i++) {
		xxx = (c[i]%p).is_zero();
	}

	single_factor< Fp_polynomial > term;
	Fp_polynomial F;

	if (xxx) {
		// STEP II  (s is zero on entry)
		xxx = 1;
		for (i = 0; i < 5 && xxx; i++) {
			xxx = (c[i]%p2).is_zero();
		}
		bigint *dd = new bigint[5];
		if (xxx) {
			for (i = 0; i < 5; i++) {
				dd[i] = c[i]/p2;
			}
			int res = samir_local_sol(p, dd);
			delete[] dd;
			return res;
		}
		for (i = 0; i < 5; i++) {
			dd[i] = c[i]/p;
		}
		factorization< Fp_polynomial > fact_f = fact_c(p, dd);
		// Is there a non-repeated root
		for (i = 0; i < fact_f.no_of_components(); i++) {
			e = fact_f.prime_exponent(i);
			F = fact_f.prime_base(i).base();
			if (F.degree() == 1 && e == 1) {
				delete[] dd;
				return 1;
			}
		}
		// Go thru each repeated root and make the
		// required transformation
		bigint *d = new bigint[5];
		fl = 0;
		for (i = 0; i < fact_f.no_of_components() && !fl; i++) {
			term = fact_f.prime_base(i);
			F = term.base();
			e = fact_f.prime_exponent(i);
			if (F.degree() == 1 && e != 1) {
				bigint r = -F[0];

				d[4] = dd[4]*p2*p;
				d[3] = p2*(dd[3]+4*dd[4]*r);
				d[2] = p*(dd[2]+6*dd[4]*r*r+3*dd[3]*r);
				d[1] = (2*dd[2]*r+4*dd[4]*r*r*r+3*dd[3]*r*r+dd[1]);
				d[0] = ((((dd[4]*r+dd[3])*r+dd[2])*r+dd[1])*r+dd[0])/p;
				fl = samir_local_sol(p, d);
			}
		}
		delete[] d;
		delete [] dd;
		return fl;
	}
	// STEP I
	bigint un;
	xxx = 1;
	for (i = 4; i >= 0 && xxx; i--) {
		un = c[i];
		xxx = (c[i]%p).is_zero();
	}
	// If leading non-zero term is a square return 1
	if (jacobi(un, p) == 1) {
		return 1;
	}
	// If f is a constant mod p and constant not a square return 0
	if (i == -1) {
		return 0;
	}
	// Factorize f
	factorization< Fp_polynomial > fact_f = fact_c(p, c);
	long nc = fact_f.no_of_components();
	// Check if of the form un g^2
	for (i = 0; i < nc; i++) {
		if (fact_f.prime_exponent(i)%2 != 0) {
			return 1;
		}
	}
	// Compute roots mod p
	long num_r = 0;
	for (i = 0; i < nc; i++) {
		term = fact_f.prime_base(i);
		F = term.base();
		if (F.degree() == 1) {
			r[num_r] = -F[0];
			num_r = num_r+1;
		}
	}
	if (num_r == 0) {
		return 0;
	}

	// Now need to determine g
	if (nc == 1) {
		term = fact_f.prime_base(0);
		F = term.base();
		if (fact_f.prime_exponent(0)%4 == 0) {
			F = F*F;
		}
	}
	else {
		term = fact_f.prime_base(0);
		F = term.base();
		term = fact_f.prime_base(1);
		F = F*term.base();
	}

	bigint te, g0, g1, g2;
	g2 = F[0];
	g1 = F[1];
	g0 = 0;
	if (F.degree() == 2) {
		g0 = 1;
	}
	// Now determine h
	bigint h0, h1, h2, h3, h4;
	h4 = (c[4]-un*g0*g0)/p;
	h3 = (c[3]-2*un*g1*g0)/p;
	h2 = (c[2]-2*un*g2*g0-un*g1*g1)/p;
	h1 = (c[1]-2*un*g1*g2)/p;
	h0 = (c[0]-un*g2*g2)/p;

	// For each root which is also a root of h
	// transform the equations and call again
	bigint *d = new bigint[5];
	fl = 0;
	for (i = 0; i < num_r && !fl; i++) {
		te = (((h4*r[i]+h3)*r[i]+h2)*r[i])%p;
		te = ((te+h1)*r[i]+h0)%p;
		if (te.is_zero()) {
			d[4] = c[4]*p2;
			d[3] = p*(c[3]+4*c[4]*r[i]);
			d[2] = (c[2]+6*c[4]*r[i]*r[i]+3*c[3]*r[i]);
			d[1] = (2*c[2]*r[i]+4*c[4]*r[i]*r[i]*r[i]+3*c[3]*r[i]*r[i]+c[1])/p;
			d[0] = ((((c[4]*r[i]+c[3])*r[i]+c[2])*r[i]+c[1])*r[i]+c[0])/(p2);
			fl = samir_local_sol(p, d);
		}
	}

	delete[] d;
	return fl;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
