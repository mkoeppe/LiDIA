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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/factorization.h"

#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char **argv)
{
	std::cout << "This program computes the factorization of a polynomial over\n";
	std::cout << "a finite prime field using different algorithms.\n";
	std::cout << "Options :\n";
	std::cout << " -v   print additional information (such as timinigs)\n";

	lidia_size_t verbose = 0;

	argc--;
	argv++;
	single_factor< Fp_polynomial >::set_verbose_mode(0);

	while (argc > 0) {
		if (strcmp(argv[0], "-v") == 0) {
			verbose = 1;
			single_factor< Fp_polynomial >::set_verbose_mode(1);
		}
		else {
			std::cout << "unknown option: " << argv[0] << "\n";
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
		std::cout << "degree of polynomial = " << f.degree()
			  << ", bitlength of modulus = " << p.bit_length()
			  << "\n\n";
	}

	factorization< Fp_polynomial > Berl, Canz, Fact;

	timer t;
	t.start_timer();
	berlekamp(Berl, f);
	t.stop_timer();
	std::cout << "Berlekamp : " << t.user_time() << std::endl;

	t.start_timer();
	can_zass(Canz, f);
	t.stop_timer();
	std::cout << "Cantor/Zassenhaus (Shoup/v.z.Gathen) : " << t.user_time() << std::endl;

	t.start_timer();
	factor(Fact, f);
	t.stop_timer();
	std::cout << "factor : " << t.user_time() << std::endl;

	Fp_polynomial ff;
	if (f.degree() > 0)
		ff = Berl.value();
	else
		ff = f;

	if (f != ff || Berl != Canz || Berl != Fact) {
		std::cout << "Incorrect factorization!!\n";
		std::cout << "Polynomial was: \n";
		std::cout << f << std::endl;
		std::cout << "Please report this bug.\n";
	}
	else {
		std::cout << "\nfactorization :\n";
		std::cout << Berl << "\n";
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
