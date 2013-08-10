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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/xbigfloat.h"
#include	<cassert>
#include        <cstring>


#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define LIDIA_XBIGFLOAT_QUIET        0
#define LIDIA_XBIGFLOAT_INTERACTIVE  1
#define LIDIA_XBIGFLOAT_VERBOSE      2

//
//
//

void debug_print(const xbigfloat & a)
{
	std::cout << a.get_mantissa() << " * 2^" << a.get_exponent() << std::endl;
}



//
// exact operations: +, -, *
//

void identity_test(const xbigfloat &a,
		   const xbigfloat &b,
		   const xbigfloat &c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);

	bigint    m = a.get_mantissa();
	bigint    ba, bb, q1, q2;
	xbigfloat xa, xb, xq;

	ba = randomize(m);
	bb = randomize(m);

	multiply (xa, ba, a);
	multiply (xb, bb, a);

	divide (q1, ba, bb);
	divide (q2, xa, xb);
	assert (q1 == q2);

	divide (q1, bb, ba);
	divide (q2, xb, xa);
	assert (q1 == q2);
}



void accum_test(xbigfloat& a, xbigfloat& b, xbigfloat& c)
{
	xbigfloat x = a;
	x *= b;
	assert(x == (a * b));
	x += c;
	assert(x == ((a * b) + c));
	x -= a;
	assert(x == (((a * b) + c) - a));
}



//
// operations returning relative approximations:
//
// k >= 1, c >= 0
// |x/op(f,g) -1| < 2^{-k}
// |y/op(f,g) -1| < 2^{-k-c}
//
// => |x-y| < 2^{-k+1} (2^{b(x)} + 2^{b(y)-c})
//
// requires g != 0
//

void relative_approximation_test(const xbigfloat & f,
				 const xbigfloat & g,
				 long k,
				 long c)
{
	xbigfloat x, y, z, v;
	xbigfloat gpos;

	if (k < 1)
		lidia_error_handler ("xbigfloat_appl::relative_approximation_test",
				     "k < 1");
	if (k+c < 1)
		lidia_error_handler ("xbigfloat_appl::relative_approximation_test",
				     "k+c < 1");
	if (g.is_zero())
		lidia_error_handler ("xbigfloat_appl::relative_approximation_test",
				     "g must not be equal to zero.");

	gpos = g;
	if (gpos.get_sign() < 0)
		gpos.negate();


	//
	// divide
	//
	divide (x, f, g, k);
	divide (y, f, g, k+c);
	assert(xbigfloat::check_relative_error(x, y, k, c));

	//
	// sqrt
	//
	sqrt (x, gpos, k);
	sqrt (y, gpos, k+c);
	assert(xbigfloat::check_relative_error(x, y, k, c));

	//
	// exp
	//
	sqrt (x, f, k);
	sqrt (y, f, k+c);
	assert(xbigfloat::check_relative_error(x, y, k, c));
}



//
// operations returning absolute approximations:
//
// |x - op(f,g)| < 2^{-k}
// |y - op(f,g)| < 2^{-k-c}
//
// => |x-y| < 2^{-k} + 2^{-k-c}
//
// Requires g != 0
//

void absolute_approximation_test(const xbigfloat & f,
				 const xbigfloat & g,
				 long k,
				 long c)
{
	xbigfloat x, y, z, v;
	xbigfloat gpos;

	(void)f;

	if (g.is_zero())
		lidia_error_handler ("xbigfloat_appl::absolute_approximation_test",
				     "g must not be equal to zero.");

	gpos = g;
	if (gpos.get_sign() < 0)
		gpos.negate();

	//
	// log
	//
	log(x, gpos, k);
	log(y, gpos, k+c);
	assert(xbigfloat::check_absolute_error(x, y, k, c));
}



//
// Test with constant values
//

void constants_test()
{
	xbigfloat x, v;
	random_generator rg;
	long k;

	//
	// exp
	//

	// |x/exp(0)-1| < 2^{-k}
	rg >> k;
	k %= 100;
	if (k < 0) k += 100;
	exp(x, 0, k);
	dec(x);
	x.absolute_value();

	v.assign_one();
	v >>= k;
	assert(x < v);
}



int lidia_application_runs_quiet(int argc, char *argv[])
{
	if (argc != 2)
		return LIDIA_XBIGFLOAT_VERBOSE;
	else {
		if (!std::strcmp(argv[1], "--quiet"))
			return LIDIA_XBIGFLOAT_QUIET;
		else if (!std::strcmp(argv[1], "--interactive"))
			return LIDIA_XBIGFLOAT_INTERACTIVE;
		else
			return LIDIA_XBIGFLOAT_VERBOSE;
	}
}



int main_LiDIA(int argc, char** argv)
{
	xbigfloat xa, xb, xc;
	long k, c;
	int i, noftests;
	int info, argv_check;

	argv_check = lidia_application_runs_quiet(argc, argv);

	random_generator rg;
	bigint n;
	long e;

	if (argv_check == LIDIA_XBIGFLOAT_QUIET)
		info = 0;
	else
		info = 1;

	// Test with constant values
	constants_test();

	if (argv_check == LIDIA_XBIGFLOAT_INTERACTIVE) {
		std::cout << "Please, enter number of tests : ";
		std::cin >> noftests;
	}
	else {
		// randomly chosen long and bigint
		rg >> e;
		n = randomize(e);
		noftests = 5;
	}

	for (i = 0; i < noftests; i++) {
		if (argv_check == LIDIA_XBIGFLOAT_INTERACTIVE) {
			std::cout << "Please, enter xa : ";
			std::cin >> xa;

			std::cout << "Please, enter xb : ";
			std::cin >> xb;

			std::cout << "Please, enter xc : ";
			std::cin >> xc;
		}
		else {
			n = randomize(n*n*e);
			e = e % 1000;
			//n = randomize(100) + 2;
			//e = e % 10;

			if (info) {
				std::cout << "n = " << n << std::endl;
				std::cout << "e = " << e << std::endl;
			}

			// randomly chosen xbigfloat
			xa.randomize(n, e);

			xb.randomize(n, e);
			if (xb.is_zero())
				inc(xb);

			xc.randomize(n, e);
		}

		if (info) {
			std::cout << "xa = " << xa << std::endl;
			std::cout << "xb = " << xb << std::endl;
			std::cout << "xc = " << xc << std::endl;
		}

		if (argv_check == LIDIA_XBIGFLOAT_INTERACTIVE) {
			std::cout << "Please, enter k : ";
			std::cin >> k;

			std::cout << "Please, enter c : ";
			std::cin >> c;
		}
		else {
			// randomly chosen precision k and c
			rg >> k;
			rg >> c;

			k %= 500;
			if (k <= 0)
				k += 500;
			c %= 100;
			if (c <= 0)
				c += 100;
		}

		if (info) {
			std::cout << "k = " << k << std::endl;
			std::cout << "c = " << c << std::endl;
		}

		if (info)	std::cout << "identity_test ... ";
		identity_test(xa, xb, xc);
		if (info) std::cout << "passed" << std::endl;

		if (info)	std::cout << "accum_test ... ";
		accum_test(xa, xb, xc);
		if (info) std::cout << "passed" << std::endl;

		if (info)	std::cout << "relative_approximation_test ... ";
		relative_approximation_test(xa, xb, k, c);
		if (info) std::cout << "passed" << std::endl;

		if (info)	std::cout << "absolute_approximation_test ... ";
		absolute_approximation_test(xa, xb, k, c);
		if (info) std::cout << "passed" << std::endl;
	}

	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
