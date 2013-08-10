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
#include	"LiDIA/quadratic_ideal.h"
#include	"LiDIA/quadratic_number_standard.h"
#include	"LiDIA/quadratic_number_power_product.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/xbigfloat.h"

#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void qi_appl_assert(bool expr,
		    const quadratic_order & O,
		    const quadratic_ideal & I,
		    const quadratic_ideal & J,
		    const quadratic_ideal & K,
		    const quadratic_number_standard & g,
		    const quadratic_number_power_product & h,
		    int   line)
{
	if (!expr) {
		std::cout << "ERROR in line " << line << "." << std::endl;
		std::cout << "O = " << O << std::endl;
		std::cout << "I = " << I << std::endl;
		std::cout << "J = " << J << std::endl;
		std::cout << "K = " << K << std::endl;
		std::cout << "g = " << g << std::endl;
		std::cout << "h = " << h << std::endl;
		std::cout << "evaluated h = " << h.evaluate() << std::endl;
		exit(1);
	}
}



void qiappl_get_random_order(quadratic_order & O,
			     lidia_size_t bits,
			     int real)
{
	bigint Delta, r;
	bigint S;

	S = (bigint(1) << bits);
	do {
		Delta.randomize(S);

		if (real) {
			if (Delta.is_negative())
				Delta.negate();
			if (Delta < 5)
				Delta = 5;
		}
		else {
			if (Delta.is_positive())
				Delta.negate();
			if (Delta > -3)
				Delta = -3;
		}
	} while (!is_quadratic_discriminant(Delta));

	O.assign(Delta);
}



int main (int argc, char *argv[])
{
	quadratic_order O;
	quadratic_ideal I, J, K;
	quadratic_number_standard      a_std, h;
	quadratic_number_power_product a_pp;

	lidia_size_t i, j, no_of_tests, default_no_of_tests;
	bigint Delta, lower, upper, p, start_p;
	random_generator rg;
	long exp, k;
	xbigfloat l, t;
	lidia_size_t bit_multiplier;
	bool use_rho_for_cycle;

	// Should we talk to the user ?
	//
	bool quiet;

	if (argc == 2)
		if (!strcmp(argv[1], "--quiet"))
			quiet = true;
		else
			quiet = false;
	else
		quiet = false;

	// Should we be interactive ?
	//
	bool interactive;
	default_no_of_tests = 6;

	if (quiet) {
		interactive = false;
		no_of_tests = default_no_of_tests;
	}
	else {
		// read the number of tests
		// or set it to default_no_of_tests
		std::cout << "Please, enter the number of tests (0 = default) : ";
		std::cin >> no_of_tests;

		if (no_of_tests <= 0)
			no_of_tests = default_no_of_tests;

		std::cout << "Doing " << no_of_tests << " tests." << std::endl;

		std::cout << "Do you want to enter all parameters (0 = no / 1 = yes) ? ";
		std::cin >> i;
		if (i == 0)
			interactive = false;
		else
			interactive = true;
	}


	// determine the prime bounds
	//
	if (interactive) {
		std::cout << "For each order a prime ideal over p is tested." << std::endl;

		std::cout << "Please, enter lower prime bound: ";
		std::cin >> lower;

		std::cout << "Please, enter upper prime bound: ";
		std::cin >> upper;

		std::cout << "Read lower = " << lower << std::endl;
		std::cout << "Read upper = " << upper << std::endl;
	}
	else {
		lower = 3;
		upper = 200;
	}

	bit_multiplier = 6;
	use_rho_for_cycle = true;

	// Run the tests
	//
	for (i = 0; i < no_of_tests; i++) {
		if (!quiet) {
			std::cout << std::endl;
			std::cout << "===== Test No. " << i+1 << " =====" << std::endl;
		}

		// read the order
		quadratic_order  O;

		//std::cout << "Please, enter bit multiplier : ";
		//std::cin  >> bit_multiplier;

		if (interactive) {
			std::cout << "Please, enter a quadratic order = ";
			std::cin >> O;
			std::cout << std::endl;
		}
		else
			//qiappl_get_random_order(O, (i+1)*bit_multiplier, i&1 );
			qiappl_get_random_order(O, (i+1)*bit_multiplier, 1);

		if (!quiet)
			std::cout << "quadratic order O = " << O << std::endl;

		start_p = next_prime(lower);

		for (p = start_p; p < upper; p = next_prime(p)) {
			//
			// generate ideal I
			//

			if (!interactive) {
				while (!generate_prime_ideal(I, p, O))
					p = next_prime(p);

				if (!quiet) {
					std::cout << " ----- Test for prime " << p << " ----- " << std::endl;
					std::cout << std::endl << "Generating ideal ... " << std::flush;
				}
				qi_appl_assert(I.is_normalized() == true, O, I, J, K, a_std, a_pp, __LINE__);

				//
				// raise I to power exp
				//
				rg >> exp;
				exp %= 7;
				if (exp < 0)
					exp = -exp;
				if (exp <= 1)
					exp = 2;

				power(I, I, exp);

				if (!quiet)
					std::cout << "DONE" << std::endl;
				qi_appl_assert(I.is_normalized() == true, O, I, J, K, a_std, a_pp, __LINE__);
			}
			else {
				I.assign_one(O);
				std::cout << "Please, enter ideal I : ";
				std::cin >> I;
				std::cout << "Read ideal I = " << I << std::endl;
			}

			//
			// Testing local_close in real case
			//

			if (O.discriminant() > 0) {
				if (!interactive) {
					rg >> exp;
					exp %= 100;
					if (exp < 0)
						exp = -exp;
				}
				else {
					std::cout << std::endl;
					std::cout << "Please, enter the number of rho() steps : ";
					std::cin >> exp;
					std::cout << std::endl;
				}

				if (!quiet) {
					std::cout << "Testing local_close " << std::flush;
					if (use_rho_for_cycle)
						std::cout << "(using " << exp << " rho steps) ... " << std::flush;
					else
						std::cout << "(using " << exp << " inverse_rho steps) ... " << std::flush;
				}

				// Determine minimum
				//
				J.assign(I);
				J.reduce(h);
				a_pp.assign(h);

				if (use_rho_for_cycle)
					for (j = 0; j < exp; j++) {
						J.rho(h);
						a_pp.multiply(a_pp, h);
					}
				else
					for (j = 0; j < exp; j++) {
						J.inverse_rho(h);
						a_pp.multiply(a_pp, h);
					}

				// Find that minimum with local_close
				//
				k = b_value(O.discriminant())+3;
				a_pp.get_absolute_Ln_approximation(t, k);

				J.assign(I);
				J.local_close(a_std, l, t, k);

				// Verify the result of local_close
				//
				qi_appl_assert(a_std == a_pp, O, I, J, K, a_std, a_pp, __LINE__);

				if (!quiet)
					std::cout << "DONE" << std::endl;
			}

			//
			// Testing order_close in real case
			//

			if (O.discriminant() > 0 && p == start_p) {
				if (!interactive) {
					rg >> exp;
					exp %= 20;
					if (exp < 0)
						exp = -exp;
				}
				else {
					std::cout << std::endl;
					std::cout << "Please, enter the number of rho() steps : ";
					std::cin >> exp;
					std::cout << std::endl;
				}

				if (!quiet) {
					std::cout << "Testing order_close ";
					if (use_rho_for_cycle)
						std::cout << "(using " << exp << " rho steps) ... " << std::flush;
					else
						std::cout << "(using " << exp << " inverse_rho steps) ... " << std::flush;
				}

				// Determine minimum
				//
				J.assign(O);
				a_std.assign_one(O);

				if (use_rho_for_cycle)
					for (j = 0; j < exp; j++) {
						J.rho(h);
						a_std *= h;
					}
				else
					for (j = 0; j < exp; j++) {
						J.inverse_rho(h);
						a_std *= h;
					}

				// Find that minimum with order_close
				//
				k = b_value(O.discriminant())+3;
				a_std.get_absolute_Ln_approximation(t, k);

				//std::cout << "t = "; t.print_as_bigfloat(); std::cout << std::endl;
				//std::cout << "k = " << k << std::endl;

				J.assign(O);
				J.order_close(a_pp, l, t, k);

				// Verify the result of order_close
				//
				qi_appl_assert(a_std == a_pp, O, I, J, K, a_std, a_pp, __LINE__);

				if (!quiet)
					std::cout << "DONE" << std::endl;
			}

			//
			// Testing close in real case
			//

			if (O.discriminant() > 0) {
				if (!interactive) {
					rg >> exp;
					exp %= 20;
					if (exp < 0)
						exp = -exp;
				}
				else {
					std::cout << std::endl;
					std::cout << "Please, enter the number of rho() steps : ";
					std::cin >> exp;
					std::cout << std::endl;
				}

				if (!quiet) {
					std::cout << "Testing close ";
					if (use_rho_for_cycle)
						std::cout << "(using " << exp << " rho steps) ... " << std::flush;
					else
						std::cout << "(using " << exp << " inverse_rho steps) ... " << std::flush;
				}

				// Determine minimum
				//
				J = I;
				J.reduce(a_std);

				if (use_rho_for_cycle)
					for (j = 0; j < exp; j++) {
						J.rho(h);
						a_std *= h;
					}
				else
					for (j = 0; j < exp; j++) {
						J.inverse_rho(h);
						a_std *= h;
					}

				// Find that minimum with close
				//
				k = b_value(O.discriminant())+3;
				a_std.get_absolute_Ln_approximation(t, k);

				J.assign(I);
				J.close(a_pp, l, t, k);

				// Verify the result of close
				//
				qi_appl_assert(a_std == a_pp, O, I, J, K, a_std, a_pp, __LINE__);

				if (!quiet)
					std::cout << "DONE" << std::endl;
			}

			if (use_rho_for_cycle == true)
				use_rho_for_cycle = false;
			else
				use_rho_for_cycle = true;

		} // end for (p
	} // end for (i

	if (!quiet)
		std::cout << std::endl << "All tests passed :-)" << std::endl << std::endl;

}
