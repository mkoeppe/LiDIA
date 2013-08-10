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


#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/quadratic_number_standard.h"

#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void identitytest(quadratic_number_standard & a,
		  quadratic_number_standard & b,
		  quadratic_number_standard & c)
{
	// needs c != 0

	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-c)) == -(a / c));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b + c)) == ((a * b) + (a * c)));

	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * c) / c) == a);

	quadratic_number_standard minus_one, d;

	minus_one.assign_one(a.get_order());
	minus_one.negate();

	d = a;
	d.negate();

	assert(a*minus_one == d);
}



void accumtest(quadratic_number_standard & a,
	       quadratic_number_standard & b,
	       quadratic_number_standard & c)
{
	quadratic_number_standard x = a;
	x *= b;
	assert(x == (a * b));
	x += c;
	assert(x == ((a * b) + c));
	x -= a;
	assert(x == (((a * b) + c) - a));
}



void misctest (quadratic_number_standard & a)
{
	bigrational na, nb;
	quadratic_number_standard sa;

	// a * sigma(a) = norm(a)
	a.norm(na);
	conjugate(sa, a);
	assert ((a*sa) == na);

	// is_positive, is_negative
	//
	quadratic_number_standard minus_a;

	minus_a = a;
	minus_a.negate();

	if (a.is_positive()) {
		assert(minus_a.is_negative());
	}
	else {
		assert(minus_a.is_positive());
	}

	minus_a.absolute_value();
	assert(minus_a.is_positive());

	// verify Ln approximation
	assert(a.check_Ln_correctness(20));
}



void qnappl_random_order(quadratic_order & O, lidia_size_t bits)
{
	bigint Delta, r;
	bigint S;

	S = (bigint(1) << bits);
	do {
		Delta.randomize(S);

		if (Delta.is_negative())
			Delta.negate();

		if (Delta < 5)
			Delta = 5;
	} while (!is_quadratic_discriminant(Delta));

	O.assign(Delta);
}



int main_LiDIA(int argc, char** argv)
{
	lidia_size_t i, no_of_tests;

	// should we give output ?
	bool quiet;

	if (argc == 2)
		if (!strcmp(argv[1], "--quiet"))
			quiet = true;
		else
			quiet = false;
	else
		quiet = false;

	// should we be interactive
	bool interactive;

	if (quiet) {
		interactive = false;
		no_of_tests = 5;
	}
	else {
		// read the number of tests
		// or set it to 5
		std::cout << "Please, enter the number of tests (0 = default) : ";
		std::cin >> no_of_tests;

		if (no_of_tests <= 0)
			no_of_tests = 5;

		std::cout << "Doing " << no_of_tests << " tests." << std::endl;

		std::cout << "Do you want to enter all parameters (0 = no / 1 = yes) ? ";
		std::cin >> i;
		if (i == 0)
			interactive = false;
		else
			interactive = true;
	}

	// alpha, beta with empty constructor;
	// still uninitialized
	quadratic_number_standard alpha;
	quadratic_number_standard beta;

	for (i = 0; i < no_of_tests; i++) {
		if (!quiet) {
			std::cout << std::endl;
			std::cout << "===== Test No. " << i+1 << " =====" << std::endl;
		}

		// read the order
		quadratic_order  O;

		if (interactive) {
			std::cout << "Please, enter a quadratic order = ";
			std::cin >> O;
			std::cout << std::endl;
		}
		else
			qnappl_random_order(O, (i+1)*6);

		if (!quiet)
			std::cout << "quadratic order O = " << O << std::endl;

		// set the global order to O
		quadratic_order::qo_l().set_last(&O);

		// we initialize alpha, beta with order O
		alpha.assign_order(O);
		beta.assign_order(O);

		// the empty constructor now should
		// use the globally set order O
		quadratic_number_standard gamma;

		if (interactive) {
			std::cout << "quadratic number alpha = ";
			std::cin >> alpha;

			std::cout << "quadratic number beta = ";
			std::cin >> beta;

			std::cout << "non-zero quadratic number gamma = ";
			std::cin >> gamma;
		}
		else {
			alpha.randomize();
			beta.randomize();
			do { gamma.randomize(); } while (gamma.is_zero());
		}

		if (!quiet) {
			std::cout << "alpha = " << alpha << std::endl;
			std::cout << "beta = " << beta << std::endl;
			std::cout << "gamma = " << gamma << std::endl;
		}

		if (!quiet)
			std::cout << "identitytest ... " << std::flush;

		identitytest (alpha, beta, gamma);

		if (!quiet) {
			std::cout << "passed" << std::endl;
			std::cout << "accumtest ... " << std::flush;
		}

		accumtest (alpha, beta, gamma);

		if (!quiet) {
			std::cout << "passed" << std::endl;
			std::cout << "misctest ... " << std::flush;
		}

		misctest	(alpha);

		if (!quiet) {
			std::cout << "passed" << std::endl;
		}
	}

	//std::cout << "alpha = " << alpha << "\n";
	//std::cout << "b(A) = " << alpha.b_value() << "\n";
	//std::cout << "b(A) = " << alpha.b_value(sq_Delta) << "\n";
	//std::cout.flush();

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
