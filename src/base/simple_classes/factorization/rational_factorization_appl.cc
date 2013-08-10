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


#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/timer.h"

#include        <cstdlib>

#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	rational_factorization f, g, h;
	long quiet = 0;
	long i, num_of_tests;
	bigint n, n2, upper_bound;
	timer t;

	if (argc == 2 && strcmp(argv[1], "--quiet") == 0)
		quiet = 1;

	t.set_print_mode(HMS_MODE);
	power(upper_bound, bigint(10), 20);

	if (quiet)
		num_of_tests = 6;
	else {
		std::cout << "\nTest Program for Rational Factorization";
		std::cout << "\n---------------------------------------\n\n";
		std::cout << "\n\nInput Number of Tests : "; std::cin >> num_of_tests;
	}

	for (i = 0; i < num_of_tests; i++) {
		if (quiet) {
			n = randomize(upper_bound);
			upper_bound *= 10000;
		}
		else {
			std::cout << "\n\nPlease enter number to factor : "; std::cin >> n;
		}

		f.assign(n);
		if (!quiet)
			f.verbose(1);
		g = f;

		t.start_timer();
		f.factor();
		t.stop_timer();

		if (!quiet) {
			std::cout << "\n\nFactorization of n : " << f;
			std::cout << "\nFactorization took time : " << t << std::flush;
		}
		n2 = evaluate_to_bigint(f);
		if (n != n2 && f != g) {
			std::cerr << "\n\n Error found for n = " << n << "\n\n" << std::flush;
			std::exit(1);
		}
	}
	if (!quiet)
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
