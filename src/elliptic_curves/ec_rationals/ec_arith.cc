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
//	Author	: Nigel Smart (NS)
//                Adaption of John Cremona's code
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/elliptic_curves/ec_arith.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void convert(bigint& a, const bigrational& b)
{
	if (!b.denominator().is_one()) {
		lidia_error_handler("convert", "Cannot convert bigrational to bigint");
	}
	a.assign(b.numerator());
}



base_vector< bigint >
int_roots_cubic(const bigint& a,
		const bigint& b,
		const bigint& c)
{
	bigcomplex za(a), zb(b), zc(c);
	bigcomplex* croots = solve_cubic(za, zb, zc);
	base_vector< bigint > iroots(3, 0);
	int niroots = 0;
	bigint x, cx;
	for (int i = 0; i < 3; i++) {
		croots[i].real().bigintify(x);
		if (x.is_zero()) {
			if (c.is_zero()) {
				int newone = 1;
				for (int j = 0; (j < niroots) && newone; j++)
					newone = (x != iroots[j]);
				if (newone)
					iroots[niroots++] = x;
			}
		}
		else {
			cx = c/x;
			if (x*cx == c)
				if (((x+a)*x+b+cx).is_zero()) {
					//x is a root, so add to list if not there already:
					int newone = 1;
					for (int j = 0; (j < niroots) && newone; j++)
						newone = (x != iroots[j]);
					if (newone)
						iroots[niroots++] = x;
				}
		}
	}
	delete[] croots;
	return iroots;
}



bigcomplex
discriminant(const bigcomplex& b,
	     const bigcomplex& c,
	     const bigcomplex& d)
{
	bigcomplex bb = b*b, cc = c*c, bc = b*c;
	return bigfloat(27.0)*d*d-bc*bc+bigfloat(4.0)*bb*b*d-bigfloat(18.0)*bc*d+bigfloat(4.0)*c*cc;
}



bigcomplex*
solve_cubic(const bigcomplex& c1,
	    const bigcomplex& c2,
	    const bigcomplex& c3)
{
	bigfloat three(3.0), two(2.0), one(1.0);
	bigfloat third = one/three;
	bigcomplex w = bigfloat(0.5)*bigcomplex(-1, sqrt(three));
	bigcomplex disc = discriminant(c1, c2, c3);
	bigcomplex p3 = three*c2 - c1*c1;
	bigcomplex mc1over3 = -c1*third;
	bigcomplex *roots = new bigcomplex[3];

	if (abs(disc).is_approx_zero()) {
		if (abs(p3).is_approx_zero()) {
			// triple root
			roots[0] = roots[1] = roots[2] = mc1over3;
		}
		else {
			// double root
			roots[0] = roots[1] = (c1*c2 -  bigfloat(9.0)*c3)/(p3+p3);
			roots[2] = -(roots[0] + roots[0]+c1);
		}
	}
	else {
		// distinct roots
		bigcomplex q = (((mc1over3+c1)*mc1over3 +c2)*mc1over3 +c3);
		// = F(mc1over3);
		if (abs(p3).is_approx_zero()) {
			// pure cubic
			roots[0] = power(-q, third);
			roots[1] = w*roots[0];
			roots[2] = w*roots[1];
			roots[0] += mc1over3;
			roots[1] += mc1over3;
			roots[2] += mc1over3;
		}
		else {
			bigcomplex d = bigfloat(729.0)*q*q+ bigfloat(4.0)*p3*p3*p3;
			bigcomplex t1cubed = bigfloat(0.5)*(sqrt(d)- bigfloat(27.0)*q);
			bigcomplex t1 = power(t1cubed, third);
			bigcomplex t2 = t1*w;
			bigcomplex t3 = t2*w;
			roots[0] = (-c1+t1-p3/t1)* third;
			roots[1] = (-c1+t2-p3/t2)* third;
			roots[2] = (-c1+t3-p3/t3)* third;
		}
	}
	long i, iter, niter = 2; // number of iterations in Newton refinement of roots
	for (i = 0; i < 3; i++) {
		bigcomplex z = roots[i], fz, fdashz;
		for (iter = 0; iter < niter; iter++) {
			fz = ((z+c1)*z+c2)*z+c3;
			fdashz = (3*z+2*c1)*z+c2;
			if (!fdashz.is_approx_zero()) {
				z -= fz/fdashz;
			}
		}
		roots[i] = z;
	}
	return roots;
}



bigcomplex*
solve_quartic(const bigcomplex& a,
	      const bigcomplex& b,
	      const bigcomplex& c,
	      const bigcomplex& d)
{
	bigcomplex p, q, r, aa, e, f1, f2;
	bigcomplex a4 = a/4;
	bigcomplex* roots = new bigcomplex[4];
	if (d.is_approx_zero()) {
		roots[0] = 0;
		bigcomplex* cuberoots = solve_cubic(a, b, c);
		roots[1] = cuberoots[0];
		roots[2] = cuberoots[1];
		roots[3] = cuberoots[2];
		delete[] cuberoots;
	}
	else {
		p = b -  3*a*a / 8;
		q = ((a / 2) * (a*a4 - b)) + c;
		r = ((256*d) - (64*a*c) + (16*a*a*b) - (3*a*a*a*a)) / 256;
		if (q.is_approx_zero()) {
			bigcomplex s = sqrt(p*p - 4*r);
			roots[0] = sqrt((-p + s) / 2) - a4;
			roots[1] = -sqrt((-p + s) / 2) - a4;
			roots[2] = sqrt((-p - s) / 2) - a4;
			roots[3] = -sqrt((-p - s) / 2) - a4;
		}
		else {
			bigcomplex* aaroots = solve_cubic(-p/2, -r, ((p*r)/2-(q*q)/8));
			aa = aaroots[0];
			if (aa.is_approx_zero()) aa = 0;
			e = sqrt(-p + 2*aa);
			f1 = (aa + q / (2*e));
			f2 = (aa - q / (2*e));
			bigcomplex s1 = sqrt(e*e -  4*f1);
			bigcomplex s2 = sqrt(e*e -  4*f2);
			roots[0] = ((e + s1) /  2) - a4;
			roots[1] = ((e - s1) /  2) - a4;
			roots[2] = ((-e + s2) /  2) - a4;
			roots[3] = ((-e - s2) /  2) - a4;
			delete[] aaroots;
		}
	}
	return roots;
}



// puts in decreasing order
void
orderreal(bigfloat& e1, bigfloat& e2, bigfloat& e3)
{
	if (e1 < e3)
		swap(e1, e3);
	if (e1 < e2)
		swap(e1, e2);
	else if (e2 < e3)
		swap(e2, e3);
}



bigcomplex*
solve_real_quartic(const bigfloat& a,
		   const bigfloat& b,
		   const bigfloat& c,
		   const bigfloat& d,
		   const bigfloat& e)
{
	bigfloat ii = 12*a*e - 3*b*d + c*c;
	bigfloat jj = (72*a*e + 9*b*d - 2*c*c) * c - 27*(a*d*d + b*b*e);
	bigfloat disc = 4*ii*ii*ii-jj*jj;
	bigfloat  H = 8*a*c - 3*b*b;
	bigfloat  R = b*b*b + 8*d*a*a - 4*a*b*c;
	bigfloat  Q = H*H-16*a*a*ii; // = 3*Q in JC's standard notation
	int type, nrr;

	if (disc < 0) {
		type = 3;
		nrr = 2;
	}       // 2 real roots
	else {
		if ((H< 0) && (Q > 0)) {
			type = 2;
			nrr = 4;
		}   // 4 real roots
		else {
			type = 1;
			nrr = 0;
		}   // 0 real roots
	}
	bigcomplex c1(0), c2(-3*ii), c3(jj);
	bigcomplex* cphi = solve_cubic(c1, c2, c3);
	bigcomplex* roots = new bigcomplex[4];
	bigfloat a4 = 4*a;
	bigfloat one(1);
	bigfloat oneover4a = one/a4;
	bigfloat athird = one/bigfloat(3.0);

	if (type < 3) {
		// all the phi are real;  order them so that a*phi[i] decreases
		bigfloat phi1 = real(cphi[0]);
		bigfloat phi2 = real(cphi[1]);
		bigfloat phi3 = real(cphi[2]);

		if (a > 0)
			orderreal(phi1, phi2, phi3);
		else
			orderreal(phi3, phi2, phi1);

		if (type == 2) {
			// all roots are real
			bigfloat r1 = sqrt((a4*phi1-H)*athird);
			bigfloat r2 = sqrt((a4*phi2-H)*athird);
			bigfloat r3 = sqrt((a4*phi3-H)*athird);
			if (R < 0)  r3 = -r3;
			roots[0] = bigcomplex((r1 + r2 - r3 - b) * oneover4a);
			roots[1] = bigcomplex((r1 - r2 + r3 - b) * oneover4a);
			roots[2] = bigcomplex((-r1 + r2 + r3 - b) * oneover4a);
			roots[3] = bigcomplex((-r1 - r2 - r3 - b) * oneover4a);
			// Those are all real and in descending order of size
		}
		else {
			// no roots are real
			bigfloat r1 = sqrt((a4*phi1-H)/3);
			bigfloat ir2 = sqrt(-((a4*phi2-H)/3));
			bigfloat ir3 = sqrt(-((a4*phi3-H)/3));
			if (R > 0)  r1 = -r1;
			roots[0] = bigcomplex(r1-b, ir2 - ir3) * oneover4a;
			roots[1] = conj(roots[0]);
			roots[2] = bigcomplex(-r1-b, ir2 + ir3) * oneover4a;
			roots[3] = conj(roots[2]);
		}
	}
	else {
		// disc < 0
		bigfloat realphi; // will hold the real root, which will be cphi[2]
		if (cphi[1].is_approx_real()) {
			realphi = real(cphi[1]);
			cphi[1] = cphi[2];
			cphi[2] = realphi;
		}
		else
			if (cphi[2].is_approx_real()) {
				realphi = real(cphi[2]);
			}
			else {
				realphi = real(cphi[0]);
				cphi[0] = cphi[2];
				cphi[2] = realphi;
			}
		bigcomplex r1 = sqrt((a4*cphi[0]-H)/3);
		bigfloat r3 = sqrt((a4*realphi-H)/3);
		if (R < 0)  r3 = -r3;
		roots[0] = bigcomplex(r3 - b, 2*imag(r1)) * oneover4a;
		roots[1] = conj(roots[0]);
		roots[2] = bigcomplex((2*real(r1) - r3 - b)) * oneover4a;
		roots[3] = bigcomplex((-2*real(r1) - r3 - b)) * oneover4a;
		// roots[2] and roots[3] are real
	}
	delete[] cphi;
	return roots;
}



base_vector< bigint >
square_divs(const bigint& N)
{
	rational_factorization plist(N);
	plist.factor();

	int np = plist.no_of_comp();
	int i, e = 0, nu = 1, nd = nu;

	int* elist = new int[np];
	bigint p;
	for (i = 0; i < np; i++) {
		p = plist.base(i);
		elist[i] = e = valuation(p, N)/2;
		nd *= (1+e);
	}
	base_vector< bigint > dlist(nd, 0);
	dlist[0] = 1;
	nd = nu;
	for (i = 0; i < np; i++) {
		p = plist.base(i);
		e = elist[i];
		for (int j = 0; j < e; j++)
			for (int k = 0; k < nd; k++)
				dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
		nd *= (e+1);
	}
	delete[] elist;
	return dlist;
}



bigcomplex
cagm(const bigcomplex& a,
     const bigcomplex& b)
{
	bigcomplex x = a, y = b, z;
	bigfloat theta, two(2), pi = Pi();
	while (!(abs(x-y)).is_approx_zero()) {
		z = (x+y)/two;
		y = sqrt(x*y);
		theta = two*arg(y/z);
		if (theta > pi || theta <= -pi)
			y = -y;
		x = z;
	}
	return x;
}



bigcomplex
normalize(bigcomplex& w1, bigcomplex& w2)
{
	bigcomplex tau = w1/w2, w3;
	if (tau.imag() < 0) {
		w1 = -w1;
		tau = -tau;
	}
	w1 = w1-w2*round(tau.real());
	tau = w1/w2;

	for (int i = 1; i < 50 && (abs(tau) < 1); i++) {
		// {Just to stop infinite loop due to rounding}
		w3 = -w1;
		w1 = w2;
		w2 = w3;
		tau = w1/w2;
		w1 = w1-w2*round(tau.real());
		tau = w1/w2;
	}
	return tau;
}



void
getc4c6(const bigcomplex& w1,
	const bigcomplex& w2,
	bigcomplex& c4,
	bigcomplex &c6)
{
	bigcomplex tau = w1/w2;
	bigfloat zero(0), one(1), two(2);
	bigfloat pi(Pi());
	bigfloat x = two*pi*tau.real();
	bigfloat y = two*pi*tau.imag();
	bigcomplex q = exp(-y) *  bigcomplex(cos(x), sin(x));
	bigcomplex f = two*pi/w2;
	bigcomplex f2 = f*f;
	bigcomplex f4 = f2*f2;

	bigcomplex term = bigcomplex(one);
	bigcomplex qpower = bigcomplex(one);
	bigcomplex sum4 = bigcomplex(zero);
	bigcomplex sum6 = bigcomplex(zero);

	bigfloat n, n2;

	for (n = 1; !term.is_approx_zero(); n += 1) {
		n2 = n*n;
		qpower *= q;
		term = n*n2*qpower/(one-qpower);
		sum4 += term;
		term *= n2;
		sum6 += term;
	}
	c4 = (one +  bigfloat(240.0)*sum4)*f4;
	c6 = (one -  bigfloat(504.0)*sum6)*f4*f2;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
