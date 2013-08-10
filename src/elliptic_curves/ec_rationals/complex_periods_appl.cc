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
#include	"LiDIA/point_bigint.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	long lidia_precision = 40;

	bigfloat::set_precision(lidia_precision);


	elliptic_curve< bigint > E1;

	std::cout << "Enter a curve: ";
	std::cin >> E1;
	std::cout << "Minimal model E1 = " << E1 << std::endl;

	complex_periods cp = E1.get_periods();
	std::cout << "Periods of " << E1 << " are:\n" << cp << std::endl;
	elliptic_curve< bigint > E2 = trans_to_curve(cp);
	std::cout << "Curve reconstructed from those periods = " << E2 << std::endl;

	if ((E1.get_a1() == E2.get_a1()) && (E1.get_a2() == E2.get_a2()) &&
	    (E1.get_a3() == E2.get_a3()) && (E1.get_a4() == E2.get_a4()) &&
	    (E1.get_a6() == E2.get_a6()))
		std::cout << "Curves are the same!\n";
	else
		std::cout << "Error: curves are different!\n";

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
