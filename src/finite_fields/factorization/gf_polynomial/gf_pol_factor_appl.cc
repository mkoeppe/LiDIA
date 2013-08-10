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
//	Author	: Victor Shoup, Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/timer.h"

#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char **argv)
{
	std::cout << "This program computes the factorization of a polynomial over\n";
	std::cout << "a finite field using different algorithms.\n";
	std::cout << "Options :\n";
	std::cout << " -v   print additional information (such as timinigs)\n";

	lidia_size_t verbose = 0;

	argc--;
	argv++;
	single_factor< gf_polynomial >::set_verbose_mode(0);

	while (argc > 0) {
		if (strcmp(argv[0], "-v") == 0) {
			verbose = 1;
			single_factor< gf_polynomial >::set_verbose_mode(1);
		}
		else {
			std::cout << "unknown option: " << argv[0] << "\n";
			return (-1);
		}
		argc--;
		argv++;
	}

	std::cout << "Please enter the finite field in form of \"[p, n]\", where p is the\n";
	std::cout << "characteristic and n the degree of the finite field : ";
	galois_field GF;
	std::cin >> GF;

	lidia_size_t d = 0;

	while (d < 1) {
		std::cout << "Please enter the degree of the polynomial: ";
		std::cin >> d;
	}
	std::cout << "generating a random monic polynomial..." << std::flush;

	gf_polynomial f;

	f = randomize(GF, d);
	f[d].assign_one(GF);
	std::cout << " done.\n";
	std::cout << "f = " << f << std::endl << std::endl;

	factorization< gf_polynomial > Berl, Canz, Fact;

	timer t;
	std::cout << "Berlekamp algorithm:" << std::endl;
	t.start_timer();
	berlekamp(Berl, f);
	t.stop_timer();
	std::cout << "  Total Berlekamp : " << t.user_time() << std::endl;

	std::cout << "Cantor/Zassenhaus algorithm:" << std::endl;
	t.start_timer();
	can_zass(Canz, f);
	t.stop_timer();
	std::cout << "  Total Cantor/Zassenhaus : " << t.user_time() << std::endl;

	std::cout << "Factor algorithm:" << std::endl;
	t.start_timer();
	factor(Fact, f);
	t.stop_timer();
	std::cout << "  Total Factor : " << t.user_time() << std::endl;

	if (Berl != Canz || Berl != Fact || Berl.value() != f) {
		std::cout << "Incorrect factorization!!\n";
		std::cout << "Polynomial was: \n";
		std::cout << f << std::endl;
		std::cout << "Please report this bug.\n";
		std::cout << Berl << std::endl << Canz << std::endl << Fact << std::endl;
	}
	else {
		std::cout << "\nfactorization :\n";
		std::cout << Berl << "\n";
	}

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
