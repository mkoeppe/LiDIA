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
#include	"LiDIA/timer.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	std::cout << "This program computes the factorization of a polynomial over\n";
	std::cout << "a finite prime field using the Berlekamp algorithm.\n";
	std::cout << "If the flag '-v' is used, additional information (such as\n";
	std::cout << "timinigs) are printed during the factorization process.\n\n";


	lidia_size_t verbose = 0;

	if (argc > 1 && !strcmp(argv[1], "-v")) {
		verbose = 1;
		single_factor< Fp_polynomial >::set_verbose_mode(1);
	}
	else
		single_factor< Fp_polynomial >::set_verbose_mode(0);


	std::cout << "Please enter a polynomial (example: "
	    "\"17x^54 - 3x^2 + 1 101\" means "
	    "\"17x^54 - 3x^2 + 1 mod 101\") :\n";
	Fp_polynomial f;
	std::cin >> f;

	factorization< Fp_polynomial > factors;

	timer t;
	t.start_timer();

	berlekamp(factors, f);

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
