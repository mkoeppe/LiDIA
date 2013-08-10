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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/gf_element.h"

#include	<cassert>
#include	<ctype.h>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void identitytest(const point< gf_element > & a, const point< gf_element > & b,
		  const point< gf_element > & c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a+b) == (b-a+(2*a)));
	assert((a + (-b)) == ((b - a + (2*a)) - (2*b)));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a + (b + c)) == ((a + c) + b));

	point< gf_element > d(a), e(a);

	assert(d == a);
	assert((d != a) == false);

	negate(d, a);
	assert(d.is_negative_of(a));
	assert(d == (-a));
	add(d, b, a);
	assert((a + b) == d);
	subtract(d, a , b);
	assert((a + (-b)) == d);
	negate(d, b);
	assert((a + d) == (b - a + 2*a - 2*b));
	multiply_by_2(d, a);
	assert((a + a) == d);

	subtract(d, a, a);
	assert(d.is_zero());
	add(d, a, a);
	multiply_by_2(e, a);
	assert(d == e);
}



void utiltest(const point< gf_element > & a)
{
	point< gf_element > b(a), c(a);
	point< gf_element > x(a), y(a);


	x.assign_zero();

	int i;

	for (i = 0; i < 10; ++i) {
		multiply(b, i, a);
		assert(b.on_curve());
		assert(b == x);
		x = x + a;
	}
	x.assign_zero();
	assert(x.is_zero());

	negate(c, a);

	for (i = 0; i > -10; --i) {
		multiply(b, i, a);
		assert(b == x);
		x = x + c;
	}
}



void scalar_multiplication_test(const point< gf_element > & a)
{
	bigint m1, m2;
	point< gf_element > P1(a), P2(a), P(a);

	m1.randomize(a.get_x().characteristic());
	m2.randomize(a.get_x().characteristic());

	multiply(P1, m1, a);
	multiply(P2, m2, a);
	add(P1, P1, P2);

	multiply(P, m1 + m2, a);
	assert(P == P1);
}



int main_LiDIA(int argc, char** argv)
{
	bigint p;
	int tests, i, n;

	std::cout << "\nELLIPTIC CURVES OVER FINITE FIELD\n\n";

	std::cout << "\n\nInput  characteristic : "; std::cin >> p;
	p = next_prime(p-1);
	std::cout << "--> Choosing next prime " << p;
	std::cout << "\n\nExtension degree n : "; std::cin >> n;
	std::cout << "\n";

	char model;
	do
	{
		std::cout << "Affine or projective model ? Enter 'a' or 'p' : ";
		std::cin >> model;
                model = tolower(model);
	}
	while (model != 'a' && model != 'p');

	galois_field F(p, n);
	gf_element a1, a2, a3, a4, a6;

	a1.assign_zero(F);
	a2.assign_zero(F);
	a3.assign_zero(F);
	a4.assign_zero(F);
	a6.assign_zero(F);

	elliptic_curve< gf_element > e;

	std::cout << "\n\nInput a1 : "; std::cin >> a1;
	std::cout << "\nInput a2 : "; std::cin >> a2;
	std::cout << "\nInput a3 : "; std::cin >> a3;
	std::cout << "\nInput a4 : "; std::cin >> a4;
	std::cout << "\nInput a6 : "; std::cin >> a6;
	std::cout << "\n# Tests : "; std::cin >> tests;

	std::cout << "Initialization ... " << std::flush;
	if (model == 'a')
		e.set_coefficients(a1, a2, a3, a4, a6);
	else
		e.set_coefficients(a1, a2, a3, a4, a6, elliptic_curve_flags::PROJECTIVE);

	point< gf_element > P(e);
	point< gf_element > Q(e), H(e);

	std::cout << "done." << std::endl;
	std::cout << "Starting tests " << std::flush;

	for (i = 1; i <= tests; i++) {
		P = e.random_point();
		Q = e.random_point();
		H = e.random_point();

		assert(P.on_curve());
		assert(Q.on_curve());
		assert(H.on_curve());

		identitytest(P, Q, H);
		utiltest(P);
		if (i % 20 == 1)
			scalar_multiplication_test(P);
		std::cout << "*" << std::flush;
	}
	std::cout << "\n\n";

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
