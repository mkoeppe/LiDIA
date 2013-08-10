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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/factorization.h"


#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	std::cout << "This program computes the factorization of a polynomial over\n";
	std::cout << "a finite prime field using an algorithm of Shoup and v.z.Gathen.\n";
	std::cout << "Options :\n";
	std::cout << " -v   print additional information (such as timinigs)\n";
	std::cout << " -cn  set poly_arg_bound to n\n";
	std::cout << "      poly_arg_bound controls space allocation for tables that\n";
	std::cout << "      are used in various routines.\n";
	std::cout << "      If n > 0, space is only allocated for about n numbers mod p.\n";
	std::cout << "      If n <= 0, space is allocated as to maximize speed.\n";
	std::cout << " -bm  set DDF_GCD_BLOCKING_FACTOR to m\n";
	std::cout << "      In the DDF routine, gcds are computed in blocks of this size.\n\n";

	lidia_size_t verbose = 0;

	argc--;
	argv++;
	single_factor< Fp_polynomial >::set_verbose_mode(0);

	while (argc > 0) {
		if (strcmp(argv[0], "-v") == 0) {
			verbose = 1;
			single_factor< Fp_polynomial >::set_verbose_mode(1);
		}
		else if (argv[0][0] == '-' && argv[0][1] == 'c')
			poly_argument::set_poly_arg_bound(atoi(argv[0] + 2));
		//see Fp_polynomial_util.h, class poly_argument
		else if (argv[0][0] == '-' && argv[0][1] == 'b')
			DDF_GCD_BLOCKING_FACTOR = atoi(argv[0]+2);
		//see DDF.cc
		else {
			std::cerr << "unknown option: " << argv[0] << "\n";
			return 1;
		}

		argc--;
		argv++;
	}

	std::cout << "Please enter a monic polynomial (example: x^54 - 3*x^2 + 1 mod 5) :\n";
	Fp_polynomial f;
	std::cin >> f;
	const bigint & p = f.modulus();

	if (verbose) {
		std::cerr << "can_zass:  degree of polynomial = " << f.degree()
			  << ", bitlength of modulus = " << p.bit_length()
			  << ", \n poly_arg_bound = " << poly_argument::get_poly_arg_bound()
			  << ", gcd blocking = " << DDF_GCD_BLOCKING_FACTOR
			  << "\n";
	}

	factorization< Fp_polynomial > factors;

	timer t;
	t.start_timer();

	can_zass(factors, f);

	if (verbose) {
		t.stop_timer();
		std::cout << "total time: " << t.user_time() << std::endl;
	}

	std::cout << "factorization :\n";
	std::cout << factors << "\n";

	Fp_polynomial ff;

	ff = factors.value();
	if (f != ff) {
		std::cerr << "Incorrect factorization!!\n";
		std::cout << "Product of factors :";
		std::cout << ff << std::endl;
		std::cout << "should be :" << std::endl;
		std::cout << f << std::endl;
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
