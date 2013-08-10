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


#include	"LiDIA/Fp_rational_function.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



bool is_quiet = false;



void identitytest(const Fp_rational_function & a, const Fp_rational_function & b,
                  const Fp_rational_function & c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));

	assert(a == (2*a - a));
	assert((a + (-b)) == (b - a+2*a - 2*b));
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

	Fp_rational_function d(a.modulus());

	negate(d, a); assert(-(d) == a);
	add(d, b, a); assert((a + b) == d);
	subtract(d, a , b); assert((a + (-b)) == d);
	negate(d, b); assert((a + d) == (b - a+2*a - 2*b));
	d.assign(a); d.multiply_by_2(); assert((a + a) == d);
	multiply(d, a, b); assert(d == (b * a));
	d.assign(b); d.negate(); multiply(d, a, d); assert(d == -(a * b));
	divide(d, a, b); assert((a / (-b)) == -d);
	add(d, b, c); multiply(d, a, d); assert(d == ((a * b) + (a * c)));
	d = b; d.invert(); assert(((a / b) / a) == d);
	d = b; d.invert();
	d.assign(b.numerator(), b.denominator()); assert(b == d);
}



void identity_mod_test (const Fp_rational_function & aa,
                        const Fp_rational_function & bb, const Fp_polynomial & f)
{
	Fp_poly_modulus ff;
	Fp_polynomial g;
	Fp_rational_function a(aa), b(bb);


	a.reduce(f); b.reduce(f);
	ff.build(f);

	Fp_rational_function d(a.modulus()), e(a.modulus());
	Fp_rational_function dd(a.modulus());

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
	assert(equal(g, a, f));
	assert(equal(g, a, ff));
}



void utiltest(const Fp_rational_function & a)
{
	Fp_rational_function b, c;

	square(b, a);
	multiply(c, a, a);
	assert(b == c);

	Fp_rational_function x(a), y(a);

	x.assign_one(); y.assign_one();

	for (int i = 0; i < 10; ++i) {
		x *= a;
		y = a * y;
		assert(y == x);
	}
	x.assign_one();
	assert(x.is_one());
}



void accumtest(const Fp_rational_function & a, const Fp_rational_function & b,
               const Fp_rational_function & c)
{
	Fp_rational_function x = a;
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
	Fp_rational_function result(5);

	std::cout << "\nIO Test ";
	std::cout << "\n--------\n\nThe input format for a ";
	std::cout << "Fp_rational_function is either";
	std::cout << "\ngiven as [c_n ... c_0] / [d_m ... d_0] mod p or";
	std::cout << "\nc_n*x^n + ... + c_0 / d_m*x^m + ... d_0 mod p\n\n";
	std::cout << "Please input a Fp_rational_function now : ";
	std::cin >> result;
	std::cout << "\nInput in raw format  :  "; result.print();
	std::cout << "\nInput in nice format :  " << result << "\n";
}



void all_test(const bigint & ii)
{
	Fp_rational_function a, b;
	bigint p, aa, bb;

	p = next_prime(ii);

	if (!is_quiet)
		std::cout << "\nTests with prime p = " << p << std::flush;

	a.set_modulus(p);
	b.set_modulus(p);
	a.randomize(15, 10);
	b.randomize(20, 3);

	utiltest(a);
	identitytest(a, b, b);
	identitytest(a, a, b);

	Fp_polynomial f;

	build_irred(f, a.modulus(), a.degree_numerator()+1);
	identity_mod_test(a, b, f);

	accumtest(a, a, b);
	accumtest(b, a, b);

	utiltest(a);
	utiltest(b);
}



int main_LiDIA(int argc, char** argv)
{
	int j;
	bigint i = 100000;

	if (argc == 2 && strcmp(argv[1], "--quiet") == 0)
		is_quiet = true;


	if (!is_quiet) {
		iotest();
		std::cout << "\nNow several redundancy tests are started ... \n" << std::flush;
	}

	for (j = 0; j < 10; j++) {
		all_test(randomize(i));
		multiply(i, i, 100);
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
