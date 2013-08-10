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
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/multi_bigmod.h"
#include	"LiDIA/gf2n.h"
#include	"LiDIA/galois_field.h"
#include	"LiDIA/gf_element.h"
#include	"LiDIA/elliptic_curve.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	// Read parameters.
	//
	std::cout << "\n\nECO_GF2N Verify Strong Curve Program (LiDIA Version)";
	std::cout << "\n====================================================\n\n";

	lidia_size_t n;
	std::cout << "\nExtension degree n: ";
	std::cin >> n;

	galois_field F(2, n);
	gf2n_init(n);

	gf_element a2(F);
	std::cout << "Coefficient a2 = ";
	std::cin >> a2;

	gf_element a6(F);
	std::cout << "Coefficient a6 = ";
	std::cin >> a6;

	bigint ec_order;
	std::cout << "Group order : ";
	std::cin >> ec_order;

	lidia_size_t nofProbTests;
	std::cout << "Prob. proof, number of tests = ";
	std::cin  >> nofProbTests;

	bigint cofactor;
	std::cout << "Maximal cofactor = ";
	std::cin  >> cofactor;

	int movExponent;
	std::cout << "Maximal degree of finite field (embedding) = ";
	std::cin  >> movExponent;

	std::cout << std::endl;


  // Verify the group order.
  //
	std::cout << nofProbTests << " probabilistic checks of group order ... " << std::flush;

	elliptic_curve<gf_element> E(a2, a6);

	if (E.probabilistic_test_of_group_order(ec_order, nofProbTests))
		std::cout << "PASSED." << std::endl;
	else {
		std::cout << "FAILED." << std::endl;
		return 1;
	}


	// Determine prime factorization of group order.
	//
	std::cout << "Determine prime factorization of group order ... " << std::flush;

	rational_factorization rf;
	rf.assign(ec_order);
	rf.trialdiv();
	if (decimal_length(cofactor) > 6)
		rf.ecm(decimal_length(cofactor));
	else
		rf.ecm(15);

	if (rf.is_prime_factorization())
		std::cout << "found." << std::endl;
	else {
		std::cout << "not found." << std::endl;
		return 1;
	}


	// Check size of cofactor.
	//
	std::cout << "Checking for small cofactor ... " << std::flush;

	if (ec_order/rf.base(rf.no_of_comp()-1) <= cofactor)
		std::cout << "YES (" << ec_order/rf.base(rf.no_of_comp()-1) << ")." << std::endl;
	else {
		std::cout << "NO." << std::endl;
		return 1;
	}


	// Check embedding in finite field.
	//
	std::cout << "Test whether point order divides q^k-1 ... ";

	bigint pointOrder = rf.base(rf.no_of_comp()-1);

	bigint q;
	shift_left(q, bigint(1), n);

	multi_bigmod p_mod(q, pointOrder);
	multi_bigmod ppower_mod(1, pointOrder);

	int i;

	for (i = 1; i <= movExponent; i++) {
		multiply(ppower_mod, ppower_mod, p_mod);
		if (ppower_mod.is_one()) {
			std::cout << "YES (k = " << i << ").";
			return 1;
		}
	}

	std::cout << "NO (All k <= " << movExponent << " tested)." << std::endl;


	// Resumee.
	//
	std::cout << std::endl;
	std::cout << "All tests passed. Strong curve." << std::endl;
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
