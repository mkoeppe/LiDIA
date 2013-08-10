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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf2n_rational_function.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>

#include        <cstdlib>


#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



bool is_quiet = false;



void identitytest(const gf2n_rational_function & a,
		  const gf2n_rational_function & b,
                  const gf2n_rational_function & c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));

	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-b)) == -(a / b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a + (b + c)) == ((a + c) + b));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b * c)) == ((c * a) * b));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * b) / b) == a);
	assert(((a * b) / a) == b);
	assert((a / (a / b)) == b);

	gf2n_rational_function d;

	negate(d, a); assert(-(d) == a);
	add(d, b, a); assert((a + b) == d);
	subtract(d, a , b); assert((a + (-b)) == d);
	multiply(d, a, b); assert(d == (b * a));
	d.assign(b); d.negate(); multiply(d, a, d); assert(d == -(a * b));
	divide(d, a, b); assert((a / (-b)) == -d);
	add(d, b, c); multiply(d, a, d); assert(d == ((a * b) + (a * c)));
	d = b; d.invert(); assert(((a / b) / a) == d);
	d = b; d.invert();
	d.assign(b.numerator(), b.denominator()); assert(b == d);
}

void identity_mod_test (const gf2n_rational_function & aa,
                        const gf2n_rational_function & bb,
			const gf2n_polynomial & f)
{
	gf2n_poly_modulus ff;
	gf2n_polynomial g;
	gf2n_rational_function a(aa), b(bb);


	a.reduce(f); b.reduce(f);
	ff.build(f);

	gf2n_rational_function d, e, dd;

	add_mod(d, b, a, f); add(e, b, a, ff);
	assert(e == d);
	add(dd, b, a);
	assert(equal_mod(dd, d, f));

	subtract_mod(d, a , b, f);
	subtract(e, a , b, ff);
	assert(e == d);
	subtract(dd, a, b);
	assert(equal_mod(dd, d, f));

	d.assign(a);
	subtract(a, a, b, ff);
	add(a, a, b, ff);
	assert(equal(a, d, ff));

	multiply_mod(d, a , b, f);
	multiply(e, a , b, ff);
	assert(e == d);
	multiply(dd, a, b);
	assert(equal_mod(dd, d, f));

	square_mod(d, a, f);
	square(e, a, ff);
	assert(e == d);
	square(dd, a);
	assert(equal_mod(dd, d, f));

	assert(convert_mod_status(g, a, f));
	assert(equal_mod(g, a, f));
	assert(equal(g, a, ff));
}

void utiltest(const gf2n_rational_function & a)
{
	gf2n_rational_function b, c;

	square(b, a);
	multiply(c, a, a);
	assert(b == c);

	gf2n_rational_function x(a), y(a);

	x.assign_one(); y.assign_one();

	for (int i = 0; i < 10; ++i) {
		x *= a;
		y = a * y;
		assert(y == x);
	}
	x.assign_one();
	assert(x.is_one());
}


void accumtest(const gf2n_rational_function & a,
	       const gf2n_rational_function & b,
               const gf2n_rational_function & c)
{
	gf2n_rational_function x = a;
	x *= b;
	assert(x == (b * a));
	x += c;
	assert(x == (c+(b * a)));
	x -= a;
	assert(x == (-a + c + (b * a)));
	x /= b;
	assert(x == (((b *a) -a + c) / b));

	x.assign(a);
	assert(x == a);
	multiply(x, x, b);
	assert(x == (b * a));
	add(x, x, c);
	assert(x == (c+(b * a)));
	subtract(x, x, a);
	assert(x == (-a + c + (b * a)));
	divide(x, x, b);
	assert(x == (((b *a) -a + c) / b));
}

void iotest()
{
	gf2n_rational_function result(5);

	std::cout << "\nIO Test ";
	std::cout << "\n--------\n\nThe input format for a ";
	std::cout << "gf2n_rational_function is ";
	std::cout << "\ngiven as [c_n ... c_0] / [d_m ... d_0]\n\n";
	std::cout << "Please input a gf2n_rational_function now : ";
	std::cin >> result;
	std::cout << "\nInput in raw format  :  "; result.print();
	std::cout << "\nInput in nice format :  " << result << "\n";
}


void all_test()
{
	gf2n_rational_function a, b;

	a.randomize(15, 10);
	b.randomize(20, 3);

	utiltest(a);
	identitytest(a, b, b);
	identitytest(a, a, b);

	gf2n_polynomial f;

	f.randomize(20);
	f.make_monic();
	identity_mod_test(a, b, f);

	accumtest(a, a, b);
	accumtest(b, a, b);

	utiltest(a);
	utiltest(b);
}


int main_LiDIA(int argc, char** argv)
{
	int j, n;

	if (argc == 2 && strcmp(argv[1], "--quiet") == 0)
		is_quiet = true;

	std::cout << "\n Extension degree n: ";
	std::cin >> n;
	if (n < 5) {
		std::cout << "\n\n Please use reasonable large degrees, abort ...\n\n";
		std::exit(0);
	}
	gf2n_init(n);

	if (!is_quiet) {
		// iotest();
		std::cout << "\nNow several redundancy tests are started ... \n" << std::flush;
	}

	for (j = 0; j < 10; j++) {
		all_test();
		if (!is_quiet)
			std::cout << "   -->OK " << std::flush;
	}
	if (!is_quiet)
		std::cout << "\n\nSuccessful end of Test\n";
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
