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


#include	"LiDIA/quadratic_number_power_product.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/quadratic_number_standard.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/xbigfloat.h"

#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



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


int main_LiDIA(int argc, char** argv)
{

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
	lidia_size_t no_of_tests;
	lidia_size_t default_no_of_tests = 30;
	lidia_size_t all_params;

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
		std::cin >> all_params;
		if (all_params == 0)
			interactive = false;
		else
			interactive = true;
	}


	// Run the tests
	//
	lidia_size_t i, j;
	quadratic_order O;
	quadratic_number_standard q1, q2;
	base_vector< quadratic_number_standard > v;
	base_vector< bigint > exp;
	quadratic_number_power_product_basis basis_quadratic_number;
	quadratic_number_power_product p1, p2, p3;
	xbigfloat l1, l2, l3, l4, x, y;

	random_generator rg;
	lidia_size_t sizePPB;
	long k1, k2;

	bigint max_exp = 5;
	lidia_size_t max_sizePPB = 5;
	lidia_size_t bit_multiplier = 6;

	for (i = 0; i < no_of_tests; i++) {
		if (!quiet) {
			std::cout << std::endl;
			std::cout << "===== Test No. " << i+1 << "===== " << std::endl;
		}

		// read the order
		//
		if (interactive) {
			std::cout << "Please, enter a quadratic order = ";
			std::cin >> O;
			std::cout << std::endl;
		}
		else
			qiappl_get_random_order(O, (i+1)*bit_multiplier, 1);

		if (!quiet)
			std::cout << "quadratic order O = " << O << std::endl;

		quadratic_order::qo_l().set_last(&O);


		// generate power product
		//
		if (interactive) {
			std::cout << "Please, enter a power product of quadratic numbers ";
			std::cin >> p1;
		}
		else {
			rg >> sizePPB;
			sizePPB %= max_sizePPB;
			if (sizePPB < 0)
				sizePPB = -sizePPB;
			if (sizePPB <= 1)
				sizePPB = 2;

			v.set_capacity(sizePPB);
			exp.set_capacity(sizePPB);

			for (j = 0; j < sizePPB; j++) {
				v[j].randomize();
				exp[j].randomize(max_exp);
			}

			basis_quadratic_number.set_basis(v);
			p1.set_basis(basis_quadratic_number);
			p1.set_exponents(exp);
			p1.square(p1);
		}

		if (!quiet)
			std::cout << "Power product p1 = " << p1 << std::endl;


		//
		// Test construction and operator ==
		//
		v = p1.get_basis().get_basis();
		exp = p1.get_exponents();

		assert(v.get_size() == exp.get_size());
		sizePPB = v.get_size();

		// Divide exponents by 2
		// and square basis (should equal to p1).
		for (j = 0; j < sizePPB; j++) {
			square(v[j], v[j]);
			exp[j].divide_by_2();
		}

		basis_quadratic_number.set_basis(v);
		p2.set_basis(basis_quadratic_number);
		p2.set_exponents(exp);

		if (!quiet)
			std::cout << "Power product p2 = " << p2 << std::endl;

		assert(p1 == p2);
		assert(!(p1 != p2));

		//
		// Test evaluation
		//
		q1 = p1.evaluate();
		q2 = p2.evaluate();
		assert(q1 == q2);

		//
		// Operator ==, != again
		//
		assert(p1 == q1);
		assert(q1 == p1);
		assert(p2 == q2);
		assert(q2 == p2);

		assert(!(p1 != q1));
		assert(!(q1 != p1));
		assert(!(p2 != q2));
		assert(!(q2 != p2));

		//
		// Test assignment
		//
		p3.assign(p1);
		assert(p3 == p1);

		p3.assign(q2);
		assert(p3 == q2);

		//
		// Test swap
		//
		swap(p1, p2);
		assert(p1 == p2);

		//
		// Test Ln approximation
		//
		rg >> k1;
		k1 %= 100;

		if (!quiet) {
			std::cout << "Determing absolute " << k1 << " approx. l1, l3 to Ln(p1)";
			std::cout << std::endl;
		}
		p1.get_absolute_Ln_approximation(l1, k1);
		q1.get_absolute_Ln_approximation(l3, k1);

		rg >> k2;
		k2 %= 100;

		if (!quiet) {
			std::cout << "Determing absolute " << k2 << " approx. l2, l4 to Ln(p2)";
			std::cout << std::endl;
		}
		p2.get_absolute_Ln_approximation(l2, k2);
		q2.get_absolute_Ln_approximation(l4, k2);

		if (!quiet) {
			std::cout << "l1 = " << l1 << std::endl;
			std::cout << "l2 = " << l2 << std::endl;
			std::cout << "l3 = " << l3 << std::endl;
			std::cout << "l4 = " << l4 << std::endl;
		}

		if (k1 <= k2)
			assert(xbigfloat::check_absolute_error(l1, l2, k1, k2-k1)
			       == true);
		else
			assert(xbigfloat::check_absolute_error(l2, l1, k2, k1-k2)
			       == true);

		assert(xbigfloat::check_absolute_error(l1, l3, k1, 0) == true);
		assert(xbigfloat::check_absolute_error(l2, l4, k2, 0) == true);
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
