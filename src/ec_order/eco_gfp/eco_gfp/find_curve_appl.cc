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
//	Author	: Andrea Rau (AR)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/EC_domain_parameters_P1363.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/point.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	lidia_size_t bitsize_factor = 0;
	lidia_size_t percentage = 0;
	int info = 10;

	EC_domain_parameters_P1363 gc;


	// Ask for info level.
	//
	std::cout << " 2 : plenty of information" << std::endl;
	std::cout << " 1 : some information" << std::endl;
	std::cout << " 0 : tries and time only" << std::endl;
	std::cout << "Please, enter info level = ";
	std::cin >> info;

	if (info < 0) {
	    info = 0;
	}
	else if (info > 2) {
	    info = 2;
	}
	std::cout << "Info level is " << info << std::endl;

	// bit size
	//
	std::cout << std::endl;
	std::cout << "The approx. group order of E will be #E = r*cofactor, ";
	std::cout << "r prime." << std::endl;
	std::cout << "Please, input the number of bits for r ";
	std::cout << "(" << gc.default_bitsize() << " recommended) = ";

	do {
		std::cin >> bitsize_factor;
	} while (bitsize_factor < 0);

	std::cout << "Bitsize is " << bitsize_factor << std::endl;


	// percentage
	//
	std::cout << std::endl;
	std::cout << "Please enter a percentage x. The number of bits for ";
	std::cout << "the cofactor" << std::endl;
	std::cout << "will be maximal x percent of the number of bits of r.";
	std::cout << std::endl;

	do {
		std::cout << "percentage (0-100), ";
		std::cout << "(" << gc.default_percentage() << " recommended) = ";
		std::cin >> percentage;
	} while (percentage < 0 || percentage > 100);

	std::cout << "Percentage is " << percentage << std::endl;


	//
	// Generate parameters.
	//
	int GFP = 1;

	if (percentage == 0 && bitsize_factor == 0)
		gc.generate_parameters(GFP, info);
	else if (percentage == 0)
		gc.generate_parameters(GFP, bitsize_factor, info);
	else
		gc.generate_parameters(GFP, bitsize_factor, percentage, info);


	//
	// Output parameters.
	//
	std::cout << "q : " << gc.get_q() << std::endl;
	std::cout << "a : " << gc.get_a() << std::endl;
	std::cout << "b : " << gc.get_b() << std::endl;
	std::cout << "Prime divisor r of #E : " << gc.get_r() << std::endl;
	std::cout << "Cofactor k : " << gc.get_k() << std::endl;
	std::cout << "Point of order r : " << gc.get_G() << std::endl;

	//
	// Verify parameters.
	//
	// k*G     != 0,
	// r*(k*G) == 0,
	// r*G     == 0.
	//
	point< gf_element > G, P;

	G = gc.get_G();
	multiply (P, gc.get_k(), G);
	assert(!P.is_zero());

	multiply (P, gc.get_r(), P);
	assert(P.is_zero());

	multiply (P, gc.get_r(), G);
	assert (P.is_zero());

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
