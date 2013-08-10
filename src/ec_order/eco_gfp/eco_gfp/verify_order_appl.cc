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
//	Author	: Markus Maurer
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/eco_prime.h"
#include	"LiDIA/trace_list.h"
#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/galois_field.h"
#include	"LiDIA/gf_element.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint p, ec_order;
	unsigned int n;

	//
	// Read characteristic, curve parameters, and order
	//

	std::cout << "\n\nECO Verify Order Program (LiDIA Version)";
	std::cout << "\n==========================================\n\n";

	std::cout << "\nCharacteristic p : ";
	std::cin >> p;

	if (!is_prime(p)) {
		std::cout << "Sorry, p not prime. Abort." << std::endl;
		return 1;
	}
	std::cout << "\np = " << p;
	std::cout << " (" << decimal_length(p) << ")" << std::flush;

	bigmod::set_modulus(p);
	galois_field theField(p);

	gf_element a4(theField);
	gf_element a6(theField);

	std::cout << "\n\nCoefficient a4 : ";
	std::cin >> a4;

	std::cout << "\nCoefficient a6 : ";
	std::cin >> a6;

	std::cout << "\nGroup order : ";
	std::cin >> ec_order;

	std::cout << "Number of tests n = ";
	std::cin >> n;

	//
	// verify the group order
	//
	elliptic_curve< gf_element > E(a4, a6);

	if (E.probabilistic_test_of_group_order(ec_order, n)) {
		std::cout << "\n\nProbabilistic correctness proof agrees !!";
		std::cout << "\n\nCharacteristic : "; std::cout << p;
		std::cout << "\nCoefficient a4 : "; std::cout << a4;
		std::cout << "\nCoefficient a6 : "; std::cout << a6;
		std::cout << "\n\nGroup order is " << ec_order << "." << std::endl;
		return 0;
	}
	else {
		std::cout << "\n\nProbabilistic correctness rejects candidate !!"
			  << "\n" << std::endl;
		return 1;
	}
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
