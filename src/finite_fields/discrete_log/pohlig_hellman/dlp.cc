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
//	Author	: Damian Weber (DW)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/dlp.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



static char *ord_suffix(const bigint &k)
{
	static char suffix[3];

	bigint ord = k % 10;
	int i;
	ord.intify(i);

	switch(i % 10) {
	case 1 : if (k % 100 != bigint(11)) strcpy(suffix, "st");
	else        strcpy(suffix, "th");
	break;
	case 2 : if (k % 100 != bigint(12)) strcpy(suffix, "nd");
	else        strcpy(suffix, "th");
	break;
	case 3 : if (k % 100 != bigint(13)) strcpy(suffix, "rd");
	else        strcpy(suffix, "th");
	break;
	default: if (k.is_gt_zero()) strcpy(suffix, "th");
	else  suffix[0] = 0;
	break;
	}
	return suffix;
}



inline bigint chinese_remainder_dlp(const bigint &x1, const bigint &p1,
				    const bigint &x2, const bigint &p2)
{
	bigint m, x0, u, v, up1x2, vp2x1, sum;

	multiply(m, p1, p2);
	xgcd(u, v, p1, p2);
	multiply(up1x2, p1, u);
	multiply(up1x2, up1x2, x2);
	multiply(vp2x1, p2, v);
	multiply(vp2x1, vp2x1, x1);
	add(sum, up1x2, vp2x1);
	remainder(x0, sum, m);
	if (x0.is_lt_zero())
		add(x0, x0, m);
	return x0;
}



bigint search_generator(const bigint &p,
			const bigint &g0,
			rational_factorization &f)
	// given the factorization of p-1 this function determines
	// the least generator of Fp greater or equal g0
{
	if (!is_prime(p, 7))
		lidia_error_handler("search_generator(p, g0, f)", "p is not a prime number.");

	bigint y1, g, q;
	int np = f.no_of_comp();
	int j, no_gen;

	g = g0;
	while (true) {
		no_gen = 0;
		for (j = 0; j <= np-1; j++) {
			q = f.base(j);
			power_mod(y1, g, (p-bigint(1))/q, p);
			if (y1.is_one()) no_gen++;
		}
		if (!no_gen) break;
		g += 1;
	}
	return g;
}



bigint search_generator(const bigint &p, const bigint &g0)
	// this function determines the least generator of Fp
	// greater or equal g0
{
	if (!is_prime(p, 7))
		lidia_error_handler("search_generator(p, g0)", "p is not a prime number.");
	rational_factorization f(p - bigint(1));
	f.factor();
	return search_generator(p, g0, f);
}



bigint search_generator(const bigint &p)
	// this function determines the least generator of Fp
{
	if (!is_prime(p, 7))
		lidia_error_handler("search_generator(p)", "p is not a prime number.");
	return search_generator(p, bigint(2));
}



void recursion(bigint &x, bigint &a, bigint &b,
	       const bigint &p, const bigint &q,
	       const bigint &a0, const bigint &g0)

{
	static unsigned long rem;

	rem = x.least_significant_digit() & 7;

	if (rem <= 2)		// set S1
	{
		multiply(x, x, a0);
		remainder(x, x, p);
		inc(a);
		return;
	}

	if (rem <= 5)		// set S2
	{
		multiply(x, x, g0);
		remainder(x, x, p);
		inc(b);
		return;
	}

	square(x, x); // set S3
	remainder(x, x, p);
	a.multiply_by_2();
	b.multiply_by_2();
	if (a.compare(q) >= 0)
		subtract(a, a, q);
	if (b.compare(q) >= 0)
		subtract(b, b, q);

}



bigint pohlig_hellman_shanks(const bigint &g,
			     const bigint &a,
			     const bigint &p,
			     const bigint &q)
{
	static bigint m, g0, a0, l;
	static bigint ai, bi, xi, a2i, b2i, x2i;
	static bigint u, v, s, t;
	static bigint d, k, gd, gk, i;

	m = (p-bigint(1))/q;
	power_mod(a0, a, m, p);
	power_mod(g0, g, m, p);

	xi = 1;
	ai = 0;
	bi = 0;
	recursion(xi, ai, bi, p, q, a0, g0);

	x2i = xi; a2i = ai, b2i = bi;
	//    std::cout << " " << x2i;

	recursion(x2i, a2i, b2i, p, q, a0, g0);

	//    while (x2i.compare(xi)!=0)

	while (x2i != xi) {
		recursion(x2i, a2i, b2i, p, q, a0, g0);
		//          std::cout << " " << x2i;
		recursion(x2i, a2i, b2i, p, q, a0, g0);
		//          std::cout << " " << x2i;
		//	  j+=2;
		//	  if (j%10==0) std::cout << "\n";
		recursion(xi, ai, bi, p, q, a0, g0);
	}

	//    std::cout << "\n";

	s = (ai-a2i) % q;
	if (s.is_lt_zero()) s += q;

	t = (b2i-bi) % q;
	if (t.is_lt_zero()) t += q;

	//    d=xgcd(u,v,s,q);

	d = xgcd_left(u, s, q);
	if (u.is_lt_zero()) u += q;

	if (d.is_one())
		l = (u*t) % q;
	else {
		k = ((u*t) % q)/d;
		if (k.is_le_zero()) k = k+q;
		power_mod(gd, g0, q/d, p);
		power_mod(gk, g0, k, p);

		i = 0;

		while (gk != a0) {
			i += 1;
			gk = (gk*gd) % p;

		}
		l = (k+i*q/d) % q;
	}

	return l;
}



bigint discrete_log(const bigint &a,
		    const bigint &b,
		    const bigint &p,
		    const rational_factorization &f)
{
	int j;
	int np = f.no_of_comp();
	bigint q, tq, xc, mc;

	timer ti;
	long  td;

	std::cout << "intermediate results" << "\n";

	for (j = 0; j <= np-1; j++) {
		power(q, f.base(j), f.exponent(j));
		ti.start_timer();

		tq = pohlig_hellman_shanks(a, b, p, q);
		ti.stop_timer();
		td = ti.real_time()/100;
		//
		// ------------------------------------------
		//  output mod q^j
		// ------------------------------------------
		//
		std::cout << tq << " mod " << f.base(j);
		if (f.exponent(j) > 1)
			std::cout << "^" << f.exponent(j);
		std::cout << "    (" << td/60 << ":" << td%60 << ")" << "\n";
		// ------------------------------------------

		if (j) {
			xc = chinese_remainder_dlp(tq, q, xc, mc);
			mc *= q;
		}
		else {
			xc = tq;
			mc = q;
		}
	}

	return xc;
}



// compute the discrete log of b to the base a mod p
// return positive residue mod (p-1)
//   or   -1 : p isn't a prime number
//        -2 : b doesn't lie in the subgroup generated by a

bigint  dl(const bigint &a,
	   const bigint &b,
	   const bigint &p, int verbose)
{
	bigint q, tq;
	bigint g;
	bigint up, ua, xa, da, xc, x, mc;
	bigint y1, y2;
	long td;
	int j, k, np;
	int no_gen;
	timer tr, t1, t2;

	rational_factorization f;

	tr.start_timer();

	if (!is_prime(p, 8)) {
		std::cerr << "dl(): " << p << " isn't a prime number." << std::endl;
		return bigint(-1);
	}

	if (verbose) {
		std::cout << "The DLP is " << a << "^x = " << b << " mod " << p << std::endl;
		std::cout << "factoring p-1 ..." << std::endl;
	}

	f.assign(p-bigint(1));
	f.factor();

	if (verbose)
		std::cout << p-bigint(1) << " = ";

	np = f.no_of_comp();

	if (verbose)
		for (j = 0; j <= np-1; j++) {
			if (j) std::cout << " * ";
			std::cout << f.base(j);
			if (f.exponent(j) > 1)
				std::cout << "^" << f.exponent(j);
		}

	if (verbose)
		std::cout << std::endl << std::endl;

	no_gen = 0;

	for (j = 0; j <= np-1; j++) {
		for (k = 1; k <= f.exponent(j); k++) {
			power(q, f.base(j), k);
			power_mod(y1, a, (p-bigint(1))/q, p);
			if (y1.is_one()) {
				no_gen++;
				power_mod(y2, b, (p-bigint(1))/q, p);
				if (!y2.is_one()) {
					if (verbose) {
						std::cout << "The subgroup of "
							  << q << ord_suffix(q) << " powers" << std::endl;
						std::cout << "mod " << p << std::endl
							  << "contains " << a
							  << " but not " << b << "." << std::endl;
						std::cout << "So there doesn't exist a solution."
							  << std::endl;
					}
					return bigint(-2);
				}
			}
		}
	}

	mc = p-bigint(1);

	if (no_gen) {
		g = search_generator(p, bigint(2), f);
		if (verbose)
			std::cout << "using generator " << g << std::endl;
	}
	else   g = a;

	if (verbose && g != a)
		std::cout << "solving " << g << "^y = " << b << std::endl;

	t1.start_timer();
	xc = discrete_log(g, b, p, f);
	t1.stop_timer();
	if (verbose && g != a) {
		std::cout << std::endl << "y = " << xc << " mod " << mc;
		td = t1.real_time()/100;
		std::cout << "    (" << td/60 << " min " << td%60 << " sec)"
			  << std::endl << std::endl;
	}

	if (g == a) {
		xa = 1;
	}
	else {
		if (verbose)
			std::cout << "solving " << g << "^y = " << a << std::endl;
		t2.start_timer();
		xa = discrete_log(g, a, p, f);
		t2.stop_timer();
		if (verbose)
			std::cout << std::endl << "y = " << xa << " mod " << mc;
		td = t2.real_time()/100;
		if (verbose)
			std::cout << "    (" << td/60 << " min " << td%60 << " sec)"
				  << std::endl << std::endl;
	}

	if (!xa.is_one()) {
		da = gcd(xa, p-bigint(1));
		xgcd(ua, up, xa/da, (p-bigint(1))/da);

		x = (ua*xc/da) % ((p-bigint(1))/da);
		if (x.is_lt_zero()) x = x+p-bigint(1);
	}
	else x = xc;

	if (verbose) {
		std::cout << std::endl;
		std::cout << "Solution:" << std::endl;
		std::cout << "x = " << x << " mod " << mc << std::endl;
		std::cout << "-------------------------------------" << std::endl;
	}

	tr.stop_timer();

	td = tr.real_time()/100;
	if (verbose) {
		std::cout << "    CPU time : " << td/60 << " min " << td%60
			  << " sec." << std::endl;
		std::cout << "=====================================" << std::endl;
	}

	return x;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
