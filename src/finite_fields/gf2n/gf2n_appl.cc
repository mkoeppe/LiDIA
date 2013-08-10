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


#include	"LiDIA/gf2n.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>




#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void identitytest(const gf2n & a, gf2n & b, const gf2n & c)
{
	assert((a + b) == (b + a));
	assert((a * b) == (b * a));
	assert((a - b) == (a + b));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a * b) / b) == a);
	assert(((a / b) * b) == a);
}



void assigntest(unsigned long tests)
{
	gf2n a, b, c, d;
	bigint x, q;
	unsigned long xx;

	a.assign_zero(); b = gf2n(0); assert(a == b);
	b = gf2n(bigint(0)); assert(a == b);
	assert(a.is_zero() == true);


	a.assign_one(); b = gf2n(1); assert(a == b);
	b = gf2n(bigint(1)); assert(a == b);
	assert(a.is_one() == true);

	shift_left(q, bigint(1), gf2n::get_absolute_degree());

	for (unsigned int i = 0; i < tests; i++) {
		x = randomize(q);
		xx = x.least_significant_digit();

		a.assign(x); b = gf2n(x); assert(a == b);
		c = b; assert(a == c);
		c.assign(b); assert(a == c);

		a.assign(xx); b = gf2n(xx); assert(a == b);
		c = b; assert(a == c);
		c.assign(b); assert(a == c);
	}

	a.randomize(); c = a;
	b.randomize(); d = b;

	swap(a, b);
	assert(a == d);
	assert(b == c);
}



void accumtest(const gf2n & a, const gf2n & b)
{
	gf2n c, d;

	c = a + b; add(d, a, b); assert(c == d);
	c = a - b; subtract(d, a, b); assert(c == d);
	c = a * b; multiply(d, a, b); assert(c == d);
	c = a * a; square(d, a); assert(c == d);
	c = a / b; divide(d, a, b); assert(c == d);

	c.assign(b); c.invert(); assert(c == inverse(b));
	invert(c, b); assert((c*b) == gf2n(1));

	c.assign(a); c += b; add(d, a, b); assert(c == d);
	c.assign(a); c -= b; subtract(d, a, b); assert(c == d);
	c.assign(a); c *= b; multiply(d, a, b); assert(c == d);
	c.assign(a); c /= b; divide(d, a, b); assert(c == d);
}



void utiltest(const gf2n & a)
{
	gf2n b;
	bigint bord, aord, q;

	square(b, a);
	sqrt(b, b);
	assert(b == a);
	assert(b == sqrt(a*a));

	aord = compute_order(a);
	power(b, a, aord);
	assert(b == gf2n(1));
	invert(b, a);
	bord = compute_order(b);
	assert(bord == aord);

	shift_left(q, bigint(1), gf2n::get_absolute_degree());
	b = get_generator();
	bord = compute_order(b);
	assert(bord == q-bigint(1));

	unsigned int i;
	random_generator rg;

	i = a.relative_degree();
	assert(gf2n::get_absolute_degree() % i == 0);

	rg >> i;
	i = i % gf2n::get_absolute_degree();

	do {
		i = (i+1) % gf2n::get_absolute_degree();
		if (i == 0)
			i = 1;
	}
	while (gf2n::get_absolute_degree() % i != 0);

	b = get_generator(i);
	assert(i == b.relative_degree());

	b.randomize(i);
	assert(i = b.relative_degree());
}



void iotest()
{
	gf2n result;
	unsigned int i;

	for (i = 0; i < 5; i++) {
		std::cout << "\n\nEnter a gf2n-element : ";
		std::cin >> result;

		std::cout << "\noutput in predefined output mode = " << result;
		gf2nIO::setbase(gf2nIO::Hex);
		std::cout << "\noutput in 'Hex' mode = " << result;
		gf2nIO::setbase(gf2nIO::Dec);
		std::cout << "\noutput in 'Dec' mode = " << result;
		gf2nIO::noprefix();
		std::cout << "\noutput after noprefix = " << result;

		identitytest(result, result, result);
		accumtest(result, result);
	}
	std::cout << "\n\n";
}



int main_LiDIA(int argc, char** argv)
{
	unsigned int  deg, tests, i, j;
	int quiet = 0;

	if (argc == 2 && strcmp(argv[1], "--quiet") == 0)
		quiet = 1;

	gf2n_init(10);
	gf2n a, b, c;

	if (!quiet) {
		std::cout << "\nTEST PROGRAM FOR SOME GF2N FUNCTIONS";
		std::cout << "\n====================================";
		std::cout << "\nInput Lower Bound for Absolute Extension Degree : ";
		std::cin >> j;
		std::cout << "\nInput Upper Bound for Absolute Extension Degree : ";
		std::cin >> deg;
		std::cout << "\nNumber of Random Tests : "; std::cin >> tests;
		if (j > deg)
			j = deg;
	}
	else {
		j = 2;
		deg = 100;
		tests = 10;
	}

	for (; j <= deg; j++) {
		if (!quiet) {
			std::cout << "\n------------------------------------------\n";
			std::cout << "\nTest for GF(2^" << j << ") .... \n" << std::flush;
		}

		gf2n_init(j);
		a.re_initialize(); b.re_initialize(); c.re_initialize();

		for (i = 0; i < tests; i++) {
			a.randomize(); c.randomize();
			do {
				b.randomize();
			} while (b.is_zero());

			identitytest(a, b, c);
		}
		if (!quiet)
			std::cout << "\n\nIdentitytests finished\n" << std::flush;

		assigntest(tests);
		if (!quiet)
			std::cout << "Assigntests finished\n" << std::flush;

		for (i = 0; i < tests; i++) {
			a.randomize();
			do {
				b.randomize();
			} while (b.is_zero());
			accumtest(a, b);
		}
		if (!quiet)
			std::cout << "Accumtest finished\n" << std::flush;

		tests /= 100;
		if (tests < 10)
			tests = 10;

		for (i = 0; i < tests; i++) {
			do {
				a.randomize();
			} while (a.is_zero());
			utiltest(a);
		}
		if (!quiet)
			std::cout << "Utiltest finished" << std::flush;

		std::cout << "\nSolve_quadratic tests (with all subfields) ... ";
		gf2n::initialize_table_for_solve_quadratic();
		unsigned int fdeg = gf2n::extension_degree();

		for (i = 0; i < tests; i++) {
			for (unsigned int k = 2; k <= deg; k++)
				if (fdeg % k == 0) {
					c.randomize(k);
					b.randomize(k);

					add(a, c, b); // coefficient a1
					multiply(b, c, b); // coefficient a0

					assert(solve_quadratic(c, a, b));
					assert((c*c + a*c + b).is_zero());

					assert(solve_quadratic_with_table(c, a, b));
					assert((c*c + a*c + b).is_zero());
				}
		}
		if (!quiet)
			std::cout << " -->finished\n" << std::flush;

		std::cout << "\n\n -->Tests for degree " << j << " finished without errors\n\n";
		gf2n::delete_table_for_solve_quadratic();
	}
	std::cout << "\n\n";
	return (0);
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
