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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/alg_number.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	timer T;
	int i; // MM, for SunPro


	std::cout << "\nPlease enter an order\n";
	order O;
	std::cin >> O;
	order O2 = O.maximize();

	std::cout << "Maximized order is " << O2 << std::flush;

	factorization< alg_ideal > f, g;
	bigrational n;


	while (true) {
		std::cout << "\nPlease enter the number to factor "; std::cin >> n;

		std::cout << "Assigning " << n << std::endl;
		f.assign(alg_ideal(alg_number(n.numerator(), O2))/n.denominator());
		std::cout << "f is " << f << std::endl;

		g = f;

		single_factor< alg_ideal >::set_verbose_mode(true);

		T.start_timer();
		f.factor_all_components();
		T.stop_timer();

		//  std::cout << "Length of the factorization : " << f.no_of_components() << "\n";
		std::cout << "\n\tFactorization :\n" << f;

		if (f.is_prime_factorization())
			std::cout << "  is prime factorization";
		else
			std::cout << "  is no prime factorization";

		std::cout << " and was computed in time " << T << ".\n" << std::flush;

		std::cout << "\n\n Original factorization was " << g << ".\n";

		std::cout << "\n\n Checking the factorization: ";

		alg_ideal M = O2;
		alg_ideal N;

		for (i = 0; i < f.no_of_prime_components(); i++) {
			power(N, f.prime_base(i).prime_base(), f.prime_exponent(i));
			M *= N;
		}
		if (!f.is_prime_factorization())
			for (i = 0; i < f.no_of_composite_components(); i++) {
				power(N, f.composite_base(i).base(), f.composite_exponent(i));
				M *= N;
			}
		std::cout << " f evaluates to " << M;
		if ((f.no_of_prime_components() == 0 &&
		     f.no_of_composite_components() == 0 &&
		     g.no_of_composite_components() == 0) ||
		    M == g.composite_base(0).base())
			// Note that g.composite_base(0) does
			// not exist, if we factored the order!
			std::cout << " which is equal to the input.\n" << std::endl;
		else
			std::cout << "\b which differs from the input.\n" << std::endl;

		std::cout << "Ramification indizes and degrees of inertia for the primes:\n";
		std::cout << "===========================================================\n\n";
		for (i = 0; i < f.no_of_prime_components(); i++) {
			std::cout << f.prime_base(i).prime_base() << " has ram.ind. ";
			std::cout << f.prime_base(i).prime_base().ramification_index();
			std::cout << " and deg.inert. ";
			std::cout << f.prime_base(i).prime_base().degree_of_inertia() << std::endl << std::flush;
		}
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
