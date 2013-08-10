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
//	Author	: Nigel Smart (NS), John Cremona (JC)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/point_bigint.h"
#include	"LiDIA/curve_isomorphism.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	elliptic_curve< bigrational > Er1(1, 2, 3, 4, 5);
	elliptic_curve< bigrational > Er2(6, 0, 24, 16, 320);
	elliptic_curve< bigint > Ei(1, 2, 3, 4, 5);
	std::cout << "E1: " << Er1 << "\nE2: " << Er2 << "\nE3: " << Ei << std::endl << std::endl;

	curve_isomorphism< bigrational, bigint > iso1(Er1, Ei);
	curve_isomorphism< bigrational, bigint > iso2(Er2, Ei);
	curve_isomorphism< bigrational, bigrational > iso3(Er1, Er2);

	std::cout << "isomorphism E1->E3: " << iso1 << std::endl;
	std::cout << "isomorphism E2->E3: " << iso2 << std::endl;
	std::cout << "isomorphism E1->E2: " << iso3 << std::endl << std::endl;

	point< bigint > Pi(2, 3, Ei);

	std::cout << Pi << " is a point on E3" << std::endl;

	point< bigrational > Pr1 = iso1.inverse(Pi);
	point< bigrational > Pr2 = iso2.inverse(Pi);

	std::cout << Pr1 << " is its preimage on E1: " << Pr1.on_curve() << std::endl;
	std::cout << Pr2 << " is its preimage on E2: " << Pr2.on_curve() << std::endl << std::endl;

	point< bigint > twoPi = Pi+Pi;
	point< bigrational > twoPr1 = Pr1+Pr1;
	point< bigrational > twoPr2 = Pr2+Pr2;

	std::cout << "Test One (the following should be the same point on E3):\n";
	std::cout << twoPi << std::endl;
	std::cout << iso1.map(twoPr1) << std::endl;
	std::cout << iso2.map(twoPr2) << std::endl;

	std::cout << "Test Two (the following should be the same point on E1):\n";
	std::cout << twoPr1 << std::endl;
	std::cout << iso1.inverse(twoPi) << std::endl;
	std::cout << iso3.inverse(twoPr2) << std::endl;

	std::cout << "Test Three (the following should be the same point on E2):\n";
	std::cout << twoPr2 << std::endl;
	std::cout << iso2.inverse(twoPi) << std::endl;
	std::cout << iso3.map(twoPr1) << std::endl;

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
