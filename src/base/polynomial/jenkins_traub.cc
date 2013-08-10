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
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/polynomial.h"
#include	"LiDIA/bigcomplex.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// *********************************************************************
// *  Jenkins-Traub algorithm for computing the zeros of a polynomial  *
// *********************************************************************
// * Algorithm 419 for polynomials with complex coefficients,
// *		Comm. ACM 15 (1972) 97-99
// * Algorithm 493 for polynomials with real coefficients,
// *		ACM TOMS 1 (1975) 178-189
// * Both algorithms are available as Fortran code via
// * http://www.mirror.ac.uk/sites/netlib.bell-labs.com/netlib/toms/
// *
// * The code below is a direct translation of this Fortran code, so
// * vector indices run from 1 to nn (instead of starting from zero);
// * besides, polynomial coefficients are stored in order of decreasing
// * powers. Oh, and don't tell me that the code looks ugly. F77 is
// * quite an archaic language, and I couldn't be asked to reinvent the
// * wheel just for beauty's sake.
// * TPf, Feb 2000
// *
// * TODO:
// * - the original algorithms suggest scaling the polynomial
// *   coefficients under some circumstances (here skipped completely)
// * - the error bound computations in both algorithms correspond to
// *   Fortran double precision; not much thinking has been done for
// *   the adaption to bigfloats/bigcomplex
// *********************************************************************



// *********************************************************************
// * Algorithm 419
// *********************************************************************

// evaluates a polynomial p at s by the Horner recurrence
// placing the partial sums in q and returning the computed value
static bigcomplex
polyev(int nn, bigcomplex &s, const bigcomplex *p,
       bigcomplex *q)
{
	if (s.imag().is_approx_zero())
		s.assign(s.real());
	if (s.real().is_approx_zero())
		s.assign(bigcomplex(0, s.imag()));

	bigcomplex pv = p[1];
	q[1] = p[1];
	for (int i = 2; i <= nn; i++) {
		pv = pv*s + p[i];
		if (pv.imag().is_approx_zero())
			pv.assign(pv.real());
		if (pv.real().is_approx_zero())
			pv.assign(bigcomplex(0, pv.imag()));
		q[i] = pv;
	}
	return pv;
}



// bounds on the error in evaluating the polynomial by the Horner recurrence
// q	the partial sums
// ms	modulus of the point
// mp	modulus of the polynomial value
static bigfloat
errev(int nn, const bigcomplex *q, const bigfloat &ms,
      const bigfloat &mp, const bigfloat &eta)
{
	bigfloat are = eta;
	bigfloat mre = 2.0*std::sqrt(2.0)*eta;
// are, mre	error bounds on complex addition and multiplication

	bigfloat e = abs(q[1])*mre / (are+mre);
	for (int i = 1; i <= nn; i++)
		e = e*ms + abs(q[i]);
	return e * (are+mre) - mp*mre;
}



// returns a lower bound on the moduli of the zeros of a polynomial p
static bigfloat
cauchy(int nn, bigcomplex *p)
{
	bigfloat x, xm, f, dx, df;
	int i, n;
	bigfloat *q = new bigfloat[nn+1];

	bigfloat *pt = new bigfloat[nn+1]; // pt is the modulus of the coefficients
	for (i = 1; i <= nn; i++)
		pt[i] = abs(p[i]);

	pt[nn] = -pt[nn];

	// compute upper estimate of bound
	n = nn-1;
	x = exp((log(-pt[nn]) - log(pt[1])) / n);
	if (!pt[n].is_approx_zero()) {
		// if Newton step at origin is better, use it
		xm = -pt[nn]/pt[n];
		if (xm < x)
			x = xm;
	}

	// chop the interval (0,x) until f <= 0
	for (;;) {
		xm = x*0.1;
		f = pt[1];
		for (i = 2; i <= nn; i++)
			f = f*xm + pt[i];
		if (f > 0.0)
			x = xm;
		else
			break;
	}
	dx = x;

	// do Newton iteration until x converges to two decimal places
//    while (abs(dx/x) > 0.005)
	while (norm(dx/x) > 0.000025) {
		q[1] = pt[1];
		for (i = 2; i <= nn; i++)
			q[i] = q[i-1]*x + pt[i];
		f = q[nn];
		df = q[1];
		for (i = 2; i <= n; i++)
			df = df*x + q[i];
		dx = f/df;
		x = x-dx;
	}

	delete[] q;
	delete[] pt;
	return x;
}



//-----------------------------------------------------------------------------
class algorithm_419
{
	bigcomplex	s, //
		t, //
		pv; // p(s)
	bigcomplex	*p, // coeff vector in order of DECREASING powers
		*h, //
		*qp, // partial sums of pol. p (see polyev)
		*qh, // partial sums of pol. h
		*sh; // "save h"
	bigfloat	eta; // machine precision
	int		nn, // degree+1
		old_mode; // store rounding mode


	void noshft(int l1);
	bool fxshft(int l2, bigcomplex &z);
	bool vrshft(int l3, bigcomplex &z);
	bool calct();
	void nexth(bool bl);

	void init();
	void cleanup();

public:
  // The default constructor does the job and we don't want to 
  // manually initialize all members to make cppcheck happy...
  // algorithm_419() {};
	~algorithm_419() {};
	bool cpoly(const base_polynomial< bigcomplex > &, base_vector< bigcomplex > &);
};



void
algorithm_419::init()
{
	// store rounding mode
	old_mode = bigfloat::get_mode();
	bigfloat::set_mode(MP_RND);

	// init machine constants
	// compute eta such that 1+eta != 1, but 1+eta/2 is rounded towards 1
	double tt = bigfloat::get_precision()*std::log(10.0)/std::log(bigint::radix());
	long t = static_cast<long>(tt);
	if (t < tt)
		t++;
	t = (t + 3)*bigint::bits_per_digit();
	eta = 1;
	eta >>= (t-1);

	// allocate "global" arrays
	p = new bigcomplex[5*(nn+1)];
	h = p +   (nn+1);
	qp = p + 2*(nn+1);
	qh = p + 3*(nn+1);
	sh = p + 4*(nn+1);
//    p  = new bigcomplex[nn+1];
//    h  = new bigcomplex[nn+1];
//    qp = new bigcomplex[nn+1];
//    qh = new bigcomplex[nn+1];
//    sh = new bigcomplex[nn+1];
}



void
algorithm_419::cleanup()
{
	delete[] p;
//    delete[] h;
//    delete[] qp;
//    delete[] qh;
//    delete[] sh;

	// restore rounding mode
	bigfloat::set_mode(old_mode);
}



//-----------------------------------------------------------------------------

// finds the zeros of a complex polynomial
// op		polynomial
// zero		output vector of the zeros
// return false only if leading coefficient is zero or if cpoly has found
// fewer than degree zeros
bool
algorithm_419::cpoly(const base_polynomial< bigcomplex > &op,
		     base_vector< bigcomplex > &zero)
{
	bigfloat xx, yy, sinr, cosr;
#if 0
	xx = .70710678;
	yy = -xx;
	cosr = -.069756474; // was originally -.060756474, probably a typo
	sinr = .99756405;
#else
	sqrt(xx, 2);
	xx.invert();
	yy = -xx;
	sinr = sin((94*Pi())/180);
	cosr = cos((94*Pi())/180);
#endif

	int i;
	int degree = op.degree();
	nn = degree + 1;
	bool conv;

	// algorithm fails if leading coefficient is zero
	if (op[degree].is_zero())
		return false;

	zero.set_capacity(degree);
	zero.set_size(degree);

	init();

	// copy polynomial coefficients
	for (i = 0; i <= degree; i++)
		p[degree-i+1] = op[i];

	// remove zeros at the origin if any
	while (p[nn].is_zero()) {
		zero[degree-nn+1] = 0;
		nn--;
	}


	// TPf: we skip scaling the polynomial


	// start the algorithm for one zero
 l40:
	if (nn <= 2)		// calculate the final zero and return
	{
		zero[degree-1] = -p[2]/p[1];
		cleanup();
		return true;
	}

	// calculate a lower bound on the modulus of the zeros
	bigfloat bnd = cauchy(nn, p);


	bigcomplex z;
	for (int cnt1 = 1; cnt1 <= 2; i++)	// outer loop to contol 2 major passes
		// with different sequences of shifts
	{
		// first calculation, no shift
		noshft(5);

		for (int cnt2 = 1; cnt2 <= 9; cnt2++)	// inner loop to select a shift
		{
			// shift is chosen with modulus bnd and amplitude rotated by
			// 94 degrees from the previous shift
			bigfloat xxx = cosr*xx - sinr*yy;
			yy = sinr*xx + cosr*yy;
			xx = xxx;
			s = bigcomplex(bnd*xx, bnd*yy);

			// second stage calculation, fixed shift
			conv = fxshft(10*cnt2, z);

			if (conv) {
				// the second stage jumps directly to the third stage iteration
				// if successful the zero is stored and the polynomial deflated
				zero[degree-nn+1] = z;
				nn--;
				for (int i = 1; i <= nn; i++)
					p[i] = qp[i];

				goto l40;
			}
			// if the iteration is unsuccessful another shift is chosen
		}
		// if 9 shifts fail, the outer loop is repeated with another sequence
		// of shifts
	}
	// the zerofinder has failed on two major passes.
	// return empty handed
	cleanup();
	return false;
}



// computes the derivative polynomial as the initial h
// polynomial and computes l1 no-shift h polynomials
void
algorithm_419::noshft(int l1)
{
	bigfloat xni, t1, t2;
	int i, j, jj, n, nm1;

	n = nn-1;
	nm1 = n - 1;
	for (i = 1; i <= n; i++) {
		xni = nn - i;
		h[i] = xni*p[i] / n;
	}
	for (jj = 1; jj <= l1; jj++) {
//	if (abs(h[n]) > eta*10*abs(p[n]))
		if (norm(h[n]) > eta*eta*100*norm(p[n])) {
			t = -p[nn]/h[n];
			for (i = 1; i <= nm1; i++) {
				j = nn - i;
				h[j] = t*h[j-1] + p[j];
			}
			h[1] = p[1];
		}
		else {
			// if the constant term is essentially zero, shift h coefficients
			for (i = 1; i <= nm1; i++) {
				j = nn - i;
				h[j] = h[j-1];
			}
			h[1] = 0;
		}
	}
}



// computes l2 fixed-shift h polynomials and tests for convergence
// initiates a variable-shift iteration and returns with the approximate
// zero if successful
// l2	limit of fixed shift steps
// z	approximate zero if conv is true
// conv	boolean indicating convergence of stage 3 iteration
bool
algorithm_419::fxshft(int l2, bigcomplex &z)
{
	bigcomplex ot, svs;
	bool test, passed, bl;
	int i, j, n;

	bool conv;

	n = nn-1;

	// evaluate p at s
	pv = polyev(nn, s, p, qp);
	test = true;
	passed = false;

	// calculate first t = -p(s)/h(s)
	bl = calct();

	for (j = 1; j <= l2; j++)		// main loop for one second stage step
	{
		ot = t;

		// compute next h polynomial and new t
		nexth(bl);
		bl = calct();
		z = s+t;

		// test for convergence unless stage 3 has failed once or this is the
		// last h polynomial
		if (!(bl || !test || j == l2)) {
//	    if (abs(t-ot) < 0.5*abs(z))
			if (norm(t-ot) < 0.25*norm(z)) {
				if (passed) {
					// the weak convergence test has been passed twice, start
					// the third stage iteration, after saving the current h
					// polynomial and shift
					for (i = 1; i <= n; i++)
						sh[i] = h[i];
					svs = s;

					conv = vrshft(10, z);
					if (conv)
						return true;

					// the iteration failed to convergence, turn off testing
					// and restore h, s, pv and t
					test = false;
					for (i = 1; i <= n; i++)
						h[i] = sh[i];
					s = svs;

					pv = polyev(nn, s, p, qp);
					bl = calct();
					continue;
				}
				passed = true;
			}
			else
				passed = false;
		}
	}

	// attempt an iteration with final h polynomial from second stage
	return vrshft(10, z);
}



// carries out the third stage iteration
// l3	limit of steps in stage 3
// z	on entry contains the initial iterate, if the iteration converges it
//	contains the final iterate on exit
// conv	true if iteration converges
bool
algorithm_419::vrshft(int l3, bigcomplex &z)
{
	bigfloat mp, ms, omp, relstp, r1, r2, tp;
	bool b, bl;
	int i, j;

	bool conv = false;

	b = false;
	s = z;

	for (i = 1; i <= l3; i++)		// main loop for stage 3
	{
		// evaluate p at s and test for convergence
		pv = polyev(nn, s, p, qp);
		mp = abs(pv);
		ms = abs(s);
		if (mp <= 20.0*errev(nn, qp, ms, mp, eta)) {
			// polynomial value is smaller in value than a bound on the error
			// in evaluating p, terminate the iteration
			conv = true;
			z = s;
			break;
		}
		if (i != 1) {
			if (!(b || mp< omp || relstp >= 0.05)) {
				// iteration has stalled. probably a cluster of zeros. do 5
				// fixed shift steps into the cluster to force one zero to
				// dominate
				tp = relstp;
				b = true;
				if (relstp < eta)
					tp = eta;
				r1 = sqrt(tp);
				s = s * bigcomplex(1.0 + r1, r1);
				pv = polyev(nn, s, p, qp);
				for (j = 1; j <= 5; j++) {
					bl = calct();
					nexth(bl);
				}
				omp.flag() = PlusInf;
				goto l20;
			}
			if (mp*0.1 > omp)
			{	// exit if polynomial value increases significantly
				break;
			}
		}
		omp = mp;
	l20:
		bl = calct();
		nexth(bl);
		bl = calct();
		if (!bl) {
			relstp = abs(t/s);
			s += t;
		}
	}
	return conv;
}



// computes t = -p(s)/h(s)
// return true if h(s) is essentially zero
bool
algorithm_419::calct()
{
	bigcomplex hv;
	int n = nn-1;

	// evaluate h at s
	hv = polyev(n, s, h, qh);
//    bool bl = abs(hv) <= eta*10.0*abs(h[n]);
	bool bl = norm(hv) <= eta*eta*100.0*norm(h[n]);
	if (!bl)
		t = -pv/hv;
	else
		t = 0;
	return bl;
}



// calculates the next shifted h polynomial
// if bl == true then h(s) is essentially zero
void
algorithm_419::nexth(bool bl)
{
	int j, n;

	n = nn-1;
	if (!bl) {
		for (j = 2; j <= n; j++)
			h[j] = t*qh[j-1] + qp[j];
		h[1] = qp[1];
	}
	else {
		// if h(s) is zero replace h with qh
		for (j = 2; j <= n; j++)
			h[j] = qh[j-1];
		h[1] = 0;
	}
}



// *********************************************************************
// * Algorithm 419
// *********************************************************************


// divides p by the quadratic  1,u,v  placing the quotient in q and the
// remainder in  a,b
static void quadsd(int nn, const bigfloat &u, const bigfloat &v,
		   const bigfloat *p, bigfloat *q, bigfloat &a, bigfloat &b)
{
	int i;
	bigfloat c;
	b = p[1];
	q[1] = b;
	a = p[2] -u*b;
	q[2] = a;
	for (i = 3; i <= nn; i++) {
		c = p[i] - u*a - v*b;
		q[i] = c;
		b = a;
		a = c;
	}
}



// calculate the zeros of the quadratic a*z^2 + b1*z + c.
// the quadratic formula, modified to avoid overflow, is used to find the
// larger zero if the zeros are real and both zeros are complex.
// the smaller real zero is found directly from the product of the zeros c/a.
static void quadr(const bigfloat &a, const bigfloat &b1, const bigfloat &c,
		  bigcomplex &s, bigcomplex &l)
{
	bigfloat b, d, e, dabs, dsqrt;
	if (a.is_approx_zero()) {
		if (!b1.is_approx_zero())
			s = -c/b1;
		else
			s.assign_zero();
		l.assign_zero();
		return;
	}
	if (c.is_approx_zero()) {
		s.assign_zero();
		l = -b1/a;
		return;
	}

	// compute discriminant avoiding overflow
	b = b1 >> 1;
	if (abs(b) >= abs(c)) {
		e = 1.0 - (a/b)*(c/b);
		d = sqrt(abs(e))*abs(b);
	}
	else {
		e = a;
		if (c < 0.0)
			e.negate();
		e = b*(b/abs(c)) - e;
		d = sqrt(abs(e))*sqrt(abs(c));
	}
	if (e >= 0.0) {
		// real zeros
		if (b.is_gt_zero())
			d.negate();
		l = (-b+d)/a;
		if (!l.is_approx_zero())
			s = c/(-b+d);
		else {
			l.assign_zero(); // l is_approx_zero anyway
			s.assign_zero();
		}
	}
	else {
		// complex conjugate zeros
		s = bigcomplex(-b/a, abs(d/a));
		conj(l, s);
	}
}



//-----------------------------------------------------------------------------
class algorithm_493
{
	bigfloat	*p, *qp, *k, *qk, *svk;
	bigfloat	u, v, a, b, c, d, a1, a2, a3, a6, a7, e, f, g, h;
	bigcomplex	sz, lz, s;
	bigfloat	eta, are, mre;
	int		n, nn;
	int		old_mode;

	void newest(int type, bigfloat &uu, bigfloat &vv);
	void nextk(int &type);
	void calcsc(int &type);
	void realit(bigfloat &sss, int &nz, int &iflag);
	void quadit(bigfloat &uu, bigfloat &vv, int &nz);
	void fxshfr(int l2, int& nz);

	void init();
	void cleanup();

public:
	algorithm_493() {};
	~algorithm_493() {};
	bool rpoly(const base_polynomial< bigfloat > &, base_vector< bigcomplex > &);
	bool rpoly2(const base_polynomial< bigfloat > &, base_vector< bigcomplex > &);
};

//-----------------------------------------------------------------------------



void
algorithm_493::init()
{
	// store rounding mode
	old_mode = bigfloat::get_mode();
	bigfloat::set_mode(MP_RND);

	// init machine constants
	// compute eta such that 1+eta != 1, but 1+eta/2 is rounded towards 1
	double tt = bigfloat::get_precision()*std::log(10.0)/std::log(bigint::radix());
	long t = static_cast<long>(tt);
	if (t < tt)
		t++;
	t = (t + 3)*bigint::bits_per_digit();
	eta = 1;
	eta >>= (t-1);

	// are and mre refer to the unit error in + and * respectively. they are
	// assumed to be the same as eta.
	are = eta;
	mre = eta;

	// allocate "global" arrays
	p = new bigfloat[5*(nn+1)];
	qp = p +   (nn+1);
	k = p + 2*(nn+1);
	qk = p + 3*(nn+1);
	svk = p + 4*(nn+1);
}



void
algorithm_493::cleanup()
{
	delete[] p;

	// restore rounding mode
	bigfloat::set_mode(old_mode);
}

//-----------------------------------------------------------------------------



// finds the zeros of a real polynomial
// op      - vector of coefficients in order of decreasing powers
// degree  - integer degree of polynomial
// zero    - output vector of the zeros
// returns false only of leading coefficient is zero or if rpoly has found
// fewer than degree zeros.
bool
algorithm_493::rpoly(const base_polynomial< bigfloat > &op, base_vector< bigcomplex > & zero)
{
	// algorithm fails if leading coefficient is zero
	if (op.lead_coeff().is_approx_zero())
		return false;

	int degree = op.degree();

	bigfloat *pt = new bigfloat[degree+2];
	bigfloat *temp = new bigfloat[degree+2];
	bigfloat t, aa, bb, cc, xx, yy, cosr, sinr, xxx, x, bnd, xm, ff, df, dx;
	int cnt, nz, i, j, jj, nm1;
	bool zerok;

	// initialization of constants for shift rotation
#if 0
	xx = .70710678;
	yy = -xx;
	cosr = -.069756474;
	sinr = .99756405;
#else
	sqrt(xx, 2);
	xx.invert();
	yy = -xx;
	sinr = sin((94*Pi())/180);
	cosr = cos((94*Pi())/180);
#endif

	zero.set_capacity(degree);
	zero.set_size(degree);
	n = degree;
	nn = n+1;

	init();

	// remove zeros at the origin if any
	while (op[degree-n].is_approx_zero()) {
		zero[degree-n] = 0;
		nn--;
		n--;
	}

	// copy polynomial coefficients
	for (i = 1; i <= nn; i++)
		p[i] = op[degree-i+1];

 l40:// start the algorithm for one zero
	if (n <= 2)
	{   // calculate the final zero or pair of zeros
		if (n == 1)
			zero[degree-1] = -p[2]/p[1];
		else if (n == 2)
			quadr(p[1], p[2], p[3], zero[degree-2], zero[degree-1]);
		cleanup();
		delete[] pt;
		delete[] temp;
		return true;
	}

	// ...
	// TPf: skip scaling
	// ...

	// compute lower bound on moduli of zeros
	for (i = 1; i <= nn; i++)
		pt[i] = abs(p[i]);
	// compute upper estimate of bound
	x = exp((log(pt[nn])-log(pt[1]))/n);
	if (!pt[n].is_approx_zero()) {
		// if Newton step at the origin is better, use it.
		xm = pt[nn]/pt[n];
		if (xm < x)
			x = xm;
	}
	pt[nn].negate();
	// chop the interval (0,x) until ff < 0
	for (;;) {
		xm = x*.1;
		ff = pt[1];
		for (i = 2; i <= nn; i++)
			ff = ff*xm + pt[i];
		if (ff < 0)
			break;
		x = xm;
	}
	dx = x;
	// do Newton iteration until x converges to two decimal places
	while (abs(dx/x) >= .005) {
		ff = pt[1];
		df = ff;
		for (i = 2; i <= n; i++) {
			ff = ff*x + pt[i];
			df = df*x + ff;
		}
		ff = ff*x + pt[nn];
		dx = ff/df;
		x = x-dx;
	}

	bnd = x;
	// compute the derivative as the initial k polynomial
	// and do 5 steps with no shift
	nm1 = n-1;
	for (i = 2; i <= n; i++)
		k[i] = (nn-i)*p[i]/n;
	k[1] = p[1];
	aa = p[nn];
	bb = p[n];
	zerok = k[n].is_approx_zero();
	for (jj = 1; jj <= 5; jj++) {
		cc = k[n];
		if (!zerok)
		{   // use scaled form of recurrence if value of k at 0 is nonzero
			t = -aa/cc;
			for (i = 1; i <= nm1; i++) {
				j = nn - i;
				k[j] = t*k[j-1] + p[j];
			}
			k[1] = p[1];
			zerok = (abs(k[n]) <= abs(bb)*eta*10); // || k[n].is_approx_zero();
		}
		else
		{   // use unscaled form of recurrence
			for (i = 1; i <= nm1; i++) {
				j = nn - i;
				k[j] = k[j-1];
			}
			k[1] = 0;
			zerok = k[n].is_approx_zero();
		}
	}
	// save k for restarts with new shifts
	for (i = 1; i <= n; i++)
		temp[i] = k[i];
	// loop to select the quadratic corresponding to each new shift
	for (cnt = 1; cnt <= 20; cnt++) {
		// quadratic corresponds to a double shift to a non-real point and its
		// complex conjugate. the point has modulus bnd and amplitude rotated
		// by 94 degrees from the previous shift
		xxx = cosr*xx - sinr*yy;
		yy = sinr*xx + cosr*yy;
		xx = xxx;
		s = bigcomplex(bnd*xx, bnd*yy);
		u = -2*s.real();
		v = bnd;
		// second stage calculation, fixed quadratic
		fxshfr(20*cnt, nz);
		if (nz != 0) {
			// the second stage jumps directly to one of the third stage
			// iterations and returns here if successful. deflate the
			// polynomial, store the zero or zeros abd return to the main
			// algorithm
			j = degree - n;
			zero[j] = sz;
			nn -= nz;
			n = nn-1;
			for (i = 1; i <= nn; i++)
				p[i] = qp[i];
			if (nz == 1)
				goto l40;
			zero[j+1] = lz;
			goto l40;
		}
		// if the iteration is unsuccessful another quadratic is chosen after
		// restoring k
		for (i = 1; i <= n; i++)
			k[i] = temp[i];
	}
	// return with failure if no convergence with 20 shifts
	delete[] pt;
	delete[] temp;
	return false;
}



// computes up to l2 fixed shift k-polynomials, testing for convergence in the
// linear or quadratic case. initiates one of the variable shift iterations and
// returns with the number of zeros found.
// l2	- limit of fixed shift steps
// nz	- number of zeros found
void
algorithm_493::fxshfr(int l2, int &nz)
{
	bigfloat betas, betav, oss, ovv, ss, vv, ts, tv, ots, otv, tvv, tss;
	bigfloat svu, svv, ui, vi, S;
	int type, i, j, iflag;
	bool vpass, spass, vtry, stry;

	nz = 0;
	betav = .25;
	betas = .25;
	oss = s.real();
	ovv = v;
	// evaluate polynomial by synthetic division
	quadsd(nn, u, v, p, qp, a, b);
	calcsc(type);
	for (j = 1; j <= l2; j++) {
		// calculate next k polynomial and estimate v
		nextk(type);
		calcsc(type);
		newest(type, ui, vi);
		vv = vi;
		ss = 0;
		// estimate s
		if (!k[n].is_approx_zero())
			ss = -p[nn]/k[n];
		tv = 1;
		ts = 1;
		if (j != 1 && type != 3) {
			// compute relative measures of convergence of s and v seqences
			if (!vv.is_approx_zero())
				tv = abs((vv-ovv)/vv);
			if (!ss.is_approx_zero())
				ts = abs((ss-oss)/ss);
			// if decreasing, multiply two most recent convergence measures
			tvv = 1;
			if (tv < otv)
				tvv = tv*otv;
			tss = 1;
			if (ts < ots)
				tss = ts*ots;
			// compare with convergence criteria
			vpass = tvv < betav;
			spass = tss < betas;
			if (vpass || spass) {
				// at least one sequence has passed the convergence test. store
				// variables before iterating
				svu = u;
				svv = v;
				for (i = 1; i <= n; i++)
					svk[i] = k[i];
				S = ss;
				// choose iteration according to the fastest converging sequence
				vtry = false;
				stry = false;
				if (spass && ((!vpass) || tss < tvv))
					goto l40;
			l20:
				quadit(ui, vi, nz);
				if (nz > 0)
					return;

				// quadratic iteration has failed. flag that it has been tried
				// and decrease the convergence criterion
				vtry = true;
				betav *= .25;
				// try linear iteration if it has not been tried and the s
				// sequence is converging
				if (!stry && spass) {
					for (i = 1; i <= n; i++)
						k[i] = svk[i];
				l40:
					realit(S, nz, iflag);
					if (nz > 0)
						return;

					// linear iteration has failed. flag that it has been tried
					// and decrease the convergence criterion
					stry = true;
					betas *= .25;
					// if linear iteration signals an almost double real zero
					// attempt quadratic iteration
					if (iflag != 0) {
						ui = -(S+S);
						vi = S*S;
						goto l20;
					}
				}
				// restore variables
				u = svu;
				v = svv;
				for (i = 1; i <= n; i++)
					k[i] = svk[i];
				// try quadratic iteration if it has not been tried and the v
				// sequence is converging
				if (vpass && (!vtry))
					goto l20;
				// recompute qp and scalar values to continue the second stage
				quadsd(nn, u, v, p, qp, a, b);
				calcsc(type);
			}
		}
		ovv = vv;
		oss = ss;
		otv = tv;
		ots = ts;
	}
}



// variable-shift k-polynomial iteration for a quadratic factor converges only
// if the zeros are equimodular or nearly so.
// uu,vv  - coefficients of starting quadratic
// nz     - number of zero found
void
algorithm_493::quadit(bigfloat &uu, bigfloat &vv, int &nz)
{
	bigfloat ui, vi, ms, mp, omp, ee, relstp, t, zm; int type, i, j;
	bool tried;

	nz = 0;
	tried = false;
	u = uu;
	v = vv;
	j = 0;

	// main loop
	for (;;) {
		quadr(1.0, u, v, sz, lz);
		// return if roots of the quadratic are real and not close to multiple
		// or nearly equal and of opposite sign
		if (abs(abs(sz.real()) - abs(lz.real())) > .01*abs(lz.real()))
			return;

		// evaluate polynomial by quadratic synthetic division
		quadsd(nn, u, v, p, qp, a, b);
		mp = abs(a - sz.real()*b) + abs(sz.imag()*b);
		// compute a rigorous bound on the rounding error in evaluating p
		zm = sqrt(abs(v));
		ee = 2*abs(qp[1]);
		t = -sz.real()*b;
		for (i = 2; i <= n; i++)
			ee = ee*zm + abs(qp[i]);
		ee = ee*zm + abs(a+t);
		ee = (5*mre+4*are)*ee - (5*mre+2*are)*(abs(a+t)+abs(b)*zm)
			+ 2*are*abs(t);
		// iteration has converged sufficiently if the polynomial value is less
		// than 20 times this bound
		if (mp <= 20*ee) {
			nz = 2;
			return;
		}

		j++;
		// stop iteration after 20 steps
		if (j > 20)
			return;
		if (j >= 2 && !(relstp > .01 || mp < omp || tried)) {
			// a cluster appeares to be stalling the convergence.
			// five fixed shift steps are taken with a u,v close to the
			// cluster
			if (relstp < eta)
				relstp = eta;
			relstp = sqrt(relstp);
			u -= u*relstp;
			v += v*relstp;
			quadsd(nn, u, v, p, qp, a, b);
			for (i = 1; i <= 5; i++) {
				calcsc(type);
				nextk(type);
			}
			tried = true;
			j = 0;
		}
		omp = mp;
		// calculate next k polynomial and new u and v
		calcsc(type);
		nextk(type);
		calcsc(type);
		newest(type, ui, vi);
		// if vi is zero the iteration is not converging
		if (vi.is_approx_zero())
			return;
		relstp = abs((vi-v)/vi);
		u = ui;
		v = vi;
	}
}



// variable-shift h polynomial iteration for a real zero
// sss	  - starting iterate
// nz	  - number of zero found
// iflag  - flag to indicate a pair of zeros near real axis
void
algorithm_493::realit(bigfloat &sss, int &nz, int &iflag)
{
	bigfloat pv, kv, t, S;
	bigfloat ms, mp, omp, ee;
	int i, j, nm1;
	nm1 = n-1;
	nz = 0;
	S = sss;
	iflag = 0;
	j = 0;

	// main loop
	for (;;) {
		pv = p[1]; // evaluate p at s
		qp[1] = pv;
		for (i = 2; i <= nn; i++) {
			pv = pv*S + p[i];
			qp[i] = pv;
		}
		mp = abs(pv); // compute a rigorous bound on the error in
		ms = abs(S); // evaluating p
		ee = (mre/(are+mre))*abs(qp[1]);
		for (i = 2; i <= nn; i++)
			ee = ee*ms + abs(qp[i]);
		// iteration has converged sufficiently if the polynomial value is less
		// than 20 times this bound
		if (mp <= 20*((are+mre)*ee-mre*mp)) {
			nz = 1;
			sz = S;
			return;
		}

		j++;
		if (j > 10)		// stop iteration after 10 steps
			return;
		if (j >= 2 && !(abs(t) > .001*abs(s-t) || mp < omp))
		{   // a cluster of zeroes near the real axis has been encountered
			// return with iflag set to initiate a quadratic iteration
			iflag = 1;
			sss = S;
			// return if the polynomial value has increased significantly
			return;
		}

		omp = mp; // compute t, the next polynomial,
		kv = k[1]; // and the new iterate
		qk[1] = kv;
		for (i = 2; i <= n; i++) {
			kv = kv*S + k[i];
			qk[i] = kv;
		}
		if (abs(kv) >= abs(k[n])*10*eta)
		{   // use the scaled form of the recurrence if the value of k at s
			// is nonzero
			t = -pv/kv;
			k[1] = qp[1];
			for (i = 2; i <= n; i++)
				k[i] = t*qk[i-1] + qp[i];
		}
		else
		{   // use unscaled form
			k[1].assign_zero();
			for (i = 2; i <= n; i++)
				k[i] = qk[i-1];
		}

		kv = k[1];
		for (i = 2; i <= n; i++)
			kv = kv*S + k[i];
		if (abs(kv) > abs(k[n])*10*eta)
			t = -pv/kv;
		else
			t.assign_zero();
		S += t;
	}
}



// this routine calculates scalar quantities used to compute the next k
// polynomial and new estimates of the quadratic coefficients.
// return 'type' - integer variable set here indicating how the calculations
// are normalized to avoid overflow
void
algorithm_493::calcsc(int &type)
{
	// synthetic division of k by the quadratic  1,u,v
	quadsd(n, u, v, k, qk, c, d);
	if ((abs(c) <= abs(k[n])*100*eta) && (abs(d) <= abs(k[n-1])*100*eta)) {
		type = 3; // type = 3 indicates the quadratic is almost a factor of k
		return;
	}
	if (abs(d) >= abs(c)) {
		type = 2; // type = 2 indicates that all formulas are divided by d
		e = a/d;
		f = c/d;
		g = u*b;
		h = v*b;
		a3 = (a+g)*e + h*(b/d);
		a1 = b*f - a;
		a7 = (f+u)*a + h;
	}
	else {
		type = 1; // type = 1 indicates that all formulas are divided by c
		e = a/c;
		f = d/c;
		g = u*e;
		h = v*b;
		a3 = a*e + (h/c+g)*b;
		a1 = b - a*(d/c);
		a7 = a + g*d + h*f;
	}
}



// computes the next k polynomials using scalars computed in calcsc
void
algorithm_493::nextk(int &type)
{
	bigfloat temp;
	int i;

	if (type == 3)
	{	// use unscaled form of the recurrence if type is 3
		k[1].assign_zero();
		k[2].assign_zero();
		for (i = 3; i <= n; i++)
			k[i] = qk[i-2];
	}

	if (type == 1)
		temp = b;
	else
		temp = a;
	if (abs(a1) <= abs(temp)*eta*10)
	{	// if a1 is nearly zero the use a special form of the recurrence
		k[1] = 0;
		k[2] = -a7*qp[1];
		for (i = 3; i <= n; i++)
			k[i] = a3*qk[i-2] - a7*qp[i-1];
	}
	else
	{	// use scaled form of the recurrence
		a7 = a7/a1;
		a3 = a3/a1;
		k[1] = qp[1];
		k[2] = qp[2] - a7*qp[1];
		for (i = 3; i <= n; i++)
			k[i] = a3*qk[i-2] - a7*qp[i-1] + qp[i];
	}
}



// compute new estimates of the quadratic coefficients using the scalars
// computed in calcsc
void
algorithm_493::newest(int type, bigfloat &uu, bigfloat &vv)
{
	// use formulas appropriate to setting of type
	if (type != 3) {
		bigfloat a4, a5, b1, b2, c1, c2, c3, c4, temp;
		if (type != 2) {
			a4 = a + u*b + h*f;
			a5 = c + (u+v*f)*d;
		}
		else {
			a4 = (a+g)*f + h;
			a5 = (f+u)*c + v*d;
		}
		// evaluate new quadratic coefficients
		b1 = -k[n]/p[nn];
		b2 = -(k[n-1]+b1*p[n])/p[nn];
		c1 = v*b2*a1;
		c2 = b1*a7;
		c3 = b1*b1*a3;
		c4 = c1 - c2 - c3;
		temp = a5 + b1*a4 - c4;
		if (!temp.is_approx_zero()) {
			uu = u - (u*(c3+c2) + v*(b1*a1+b2*a7))/temp;
			vv = v*(1. + c4/temp);
			return;
		}
	}
	// if type=3 the quadratic is zeroed
	uu.assign_zero();
	vv.assign_zero();
}



// *********************************************************************
// * the LiDIA interface
// *********************************************************************

// if algorithm 419/493 fails
static void
last_resort(const base_polynomial< bigcomplex > &p,
	    base_vector< bigcomplex > &r)
{
	bigcomplex *z = new bigcomplex[p.degree()];
	roots(p, z);
	r = base_vector< bigcomplex > (z, p.degree());
	delete[] z;
}



base_vector< bigcomplex >
roots(const base_polynomial< bigcomplex > &p)
{
	base_vector< bigcomplex > r;

	if (p.lead_coeff().is_zero()) {
		lidia_error_handler("jenkins_traub.c",
				    "roots(base_polynomial< bigcomplex > )::"
				    "leading coefficient may not be zero");
	}
	else {
		bool real = true;
		int i, n = p.degree();
		for (i = 0; (i <= n) && real; i++)
			real = p[i].imag().is_zero();

		if (real) {
			base_polynomial< bigfloat > P;
			P.set_degree(n);
			for (i = 0; i <= n; i++)
				P[i] = p[i].real();
			r = roots(P);
		}
		else {
			algorithm_419 algorithm_419_instance;
			bool ok = algorithm_419_instance.cpoly(p, r);

			if (!ok) {
				lidia_warning_handler("algorithm_419",
						      "could not compute zeros");
				last_resort(p, r);
			}
		}
	}
	return r;
}



base_vector< bigcomplex >
roots(const base_polynomial< bigfloat > &p)
{
	base_vector< bigcomplex > r;

	if (p.lead_coeff().is_zero()) {
		lidia_error_handler("jenkins_traub.c",
				    "roots(base_polynomial< bigfloat > )::"
				    "leading coefficient may not be zero");
	}
	else {
		algorithm_493 algorithm_493_instance;
		bool ok = algorithm_493_instance.rpoly(p, r);

		if (!ok) {
			lidia_warning_handler("algorithm_493", "could not compute zeros");
			base_polynomial< bigcomplex > P = polynomial< bigfloat > (p);
			last_resort(P, r);
		}
	}
	return r;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
