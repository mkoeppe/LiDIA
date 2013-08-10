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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint D, Dmod4, x, y, p;

	std::cout << "Tries to find a (x, y), such that" << std::endl;
	std::cout << "x^2 + |D| * y^2 = 4 * p, where" << std::endl;
	std::cout << "D = 0, 1 mod 4, D < 0, |D| < 4*p, p odd" << std::endl;
	std::cout << std::endl;

	std::cout << "Please, enter odd p = ";
	std::cin >> p;
	std::cout << "Please, enter D = 0, 1 mod 4, D = ";
	std::cin >> D;

	if (!p.is_odd()) {
		std::cout << "p must be odd" << std::endl;
		return 1;
	}

	if (D >= 0) {
		std::cout << "D must be negative" << std::endl;
		return 1;
	}

	Dmod4 = D % 4;  // Note: sign(D % 4) == sign(D)

	if (Dmod4 != 0 && Dmod4 != -3) {
		std::cout << "D must be = 0, 1 mod 4" << std::endl;
		return 1;
	}

	if (cornacchia(x, y, D, p)) {
		std::cout << "Found x = " << x << std::endl;
		std::cout << "Found y = " << y << std::endl;
		std::cout << "Verification ... ";
		assert(x*x - D*y*y == 4*p);
		std::cout << "correct." << std::endl;
	}
	else
		std::cout << "No solution found!" << std::endl;

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
