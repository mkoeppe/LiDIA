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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/eco_gf2n.h"

#include        <cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint q;
	int n;

	//
	// Read characteristic and curve parameters
	//
	std::cout << "\n\nECO GF(2^n) Point Counting Program (LiDIA Version)";
	std::cout << "\n====================================================\n\n";

	std::cout << "\nExtension degree n: ";
	std::cin >> n;

	if (n < 5) {
		std::cerr << "\n\nPlease use reasonable large degrees, aborting ...\n\n";
		std::exit(1);
	}

	shift_left(q, bigint(1), n);
	gf2n_init(n);
	gf2n a6;

	std::cout << "\nField has " << q << " (" << decimal_length(q) << ") elements";
	std::cout << std::flush;

	std::cout << "\n\nElliptic curve :  Y^2 + X * Y = X^3 + a_6\n";
	std::cout << "\nInput coefficient a6 : ";

	std::cin >> a6;

	eco_gf2n ec(a6);
	ec.set_info_mode(1);
	int strat;

	std::cout << "\n\nGeneral strategies: \n";
	std::cout << "\n(1) Splitting degree computation always";
	std::cout << "\n(2) Splitting degree computation only for Elkies primes";
	std::cout << "\n(3) Splitting degree computation only for Atkin primes";
	std::cout << "\n(4) Splitting degree computation never\n";
	std::cout << "\nPlease input index of strategy : "; std::cin >> strat;

	switch (strat) {
	case 1:
		ec.set_strategy(eco_gf2n::COMPUTE_SP_DEGREE);
		break;
	case 2:
		ec.set_strategy(eco_gf2n::COMPUTE_SP_DEGREE_IF_ELKIES);
		break;
	case 3:
		ec.set_strategy(eco_gf2n::COMPUTE_SP_DEGREE_IF_ATKIN);
		break;
	case 4:
		ec.set_strategy(eco_gf2n::DONT_COMPUTE_SP_DEGREE);
		break;
	default:
		std::cerr << "\n\nNo implemented strategy, aborting ...\n\n";
		std::exit(1);
	}

	std::cout << "\n\nStrategies for finding the eigenvalue: \n";
	std::cout << "\n(1) Division polynomials with test of X-coordinate";
	std::cout << "\n(2) Division polynomials with test of Y-coordinate";
	std::cout << "\n(3) Rational functions with test of X-coordinate";
	std::cout << "\n(4) Rational functions with test of Y-coordinate";
	std::cout << "\n(5) Rational functions with test of X-coordinate and table";
	std::cout << "\n(6) Rational functions with test of Y-coordinate and table";
	std::cout << "\n(7) Babystep Giantstep"; 
	std::cout << "\n(8) Funny Babystep Giantstep";
	std::cout << "\n(9) Funny Babystep Giantstep with computation of sign";
	std::cout << "\n\nPlease input index of eigenvalue search strategy : "; std::cin >> strat;

	switch (strat) {
	case 1:
		ec.set_ev_search_strategy(eco_gf2n::EV_DIVISION_POLYNOMIAL |
					  eco_gf2n::TEST_X_COORDINATE);
		break;
	case 2:
		ec.set_ev_search_strategy(eco_gf2n::EV_DIVISION_POLYNOMIAL |
					  eco_gf2n::TEST_Y_COORDINATE);
		break;
	case 3:
		ec.set_ev_search_strategy(eco_gf2n::EV_RATIONAL_FUNCTION |
					  eco_gf2n::TEST_X_COORDINATE);
		break;
	case 4:
		ec.set_ev_search_strategy(eco_gf2n::EV_RATIONAL_FUNCTION |
					  eco_gf2n::TEST_Y_COORDINATE);
		break;
	case 5:
		ec.set_ev_search_strategy(eco_gf2n::EV_RATIONAL_FUNCTION_TABLE |
					  eco_gf2n::TEST_X_COORDINATE);
		break;
	case 6:
		ec.set_ev_search_strategy(eco_gf2n::EV_RATIONAL_FUNCTION_TABLE |
					  eco_gf2n::TEST_Y_COORDINATE);
		break;
	case 7:
		ec.set_ev_search_strategy(eco_gf2n::EV_BABYSTEP_GIANTSTEP);
		break;
	case 8:
		ec.set_ev_search_strategy(eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP);
		break;
	case 9:
		ec.set_ev_search_strategy(eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN);
		break;
	default:  std::cerr << "\n\nNo implemented strategy, aborting ...\n\n";
		std::exit(1);
	}


	//
	// compute the group order
	//

	bigint order = ec.compute_group_order();
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
