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
//	Author	:Markus Maurer (MM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/wep_rat_function.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint p;
	int i;
	bool quiet = false;

	if (argc == 2 && strcmp(argv[1], "--quiet") == 0)
		quiet = true;

	if (!quiet) {
		std::cout << "\n\nTest Program for wep_rat_functions"
			 << "\n==================================\n"
			 << "\nCharacteristic : ";
		std::cin >> p;
	}
	else power(p, 10, 8);

	p = next_prime(p-1);
	if (!quiet)
		std::cout << "\nChoose next prime, p = " << p << " (" << decimal_length(p) << ")" << std::flush;
	bigmod::set_modulus(p);
	bigmod a4, a6;
	int deg_modulus;

	if (!quiet) {
		std::cout << "\n\nCoefficient a4 : "; std::cin >> a4;
		std::cout << "\nCoefficient a6 : "; std::cin >> a6;
		std::cout << "\nDegree of Modulus : "; std::cin >> deg_modulus;
	}
	else {
		a4.randomize();
		a6.randomize();
		deg_modulus = 50;
	}

	Fp_polynomial m;
	build_irred(m, p, deg_modulus);
	std::cout << "\nThe modulus chosen is " << m << "." << std::endl;
	Fp_poly_modulus pm(m);

	wep_rat_function::initialize(a4, a6, pm);
	wep_rat_function P, mP, m2P, m3P;

	P.assign_xy(); mP.assign_zero(); m2P.assign_zero(); m3P.assign_zero();
	assert(P.on_curve());
	assert(mP.on_curve());


	if (!quiet)
		std::cout << "\nAddition Tests ... " << std::flush;

	for (i = 1; i <= 10; i++) {
		add(mP, mP, P);
		assert(mP.on_curve());
	}
	if (!quiet)
		std::cout << "  DONE.  " << std::flush;

	if (!quiet)
		std::cout << "\nSubtraction Tests ... " << std::flush;

	mP.assign_zero(); m2P.assign_zero();

	for (i = 1; i <= 10; i++) {
		subtract(mP, mP, P);
		add(m2P, m2P, P);
		assert(mP.on_curve());
		negate(m3P, m2P);
		assert(equal(m3P.get_x(), mP.get_x(), pm));
		assert(equal(m3P.get_y(), mP.get_y(), pm));
	}
	if (!quiet)
		std::cout << "  DONE.  " << std::flush;

	if (!quiet)
		std::cout << "\nDoubling Tests ... " << std::flush;

	assert(P.on_curve());

	multiply_by_2(mP, P);
	assert(mP.on_curve());

	for (i = 1; i < 10; i++) {
		multiply_by_2(mP, mP);
		assert(mP.on_curve());
	}
	if (!quiet)
		std::cout << "  DONE.  " << std::flush;

	if (!quiet)
		std::cout << "\nMultiplication Tests ... " << std::flush;


	mP.assign_zero();
	for (i = 1; i <= 50; i++) {
		add(mP, mP, P);
		multiply(m2P, i, P);
		if(!equal(m2P.get_x(), mP.get_x(), pm) ||
		   !equal(m2P.get_y(), mP.get_y(), pm)) {
		    std::cerr << "mP = " << mP << "\n"
			      << "m2P = " << m2P << "\n"
			      << "i = " << i << "\n";
		    lidia_error_handler("wep_rat_function_appl",
					"successive add and multiply differ");
		}
	}
	if (!quiet)
		std::cout << "  DONE.\n\n" << std::flush;
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
