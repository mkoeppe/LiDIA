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
//	Author	: Roland Dreier (RD) (dreier@math.berkeley.edu)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint_polynomial.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv) {
	polynomial< bigint > f, g;

	std::cout << "Enter f : ";
	std::cin >> f;

	std::cout << "Enter g : ";
	std::cin >> g;

	std::cout << "resultant(f, g) = " << resultant(f, g) << "\n";
	std::cout << "discriminant(f) = " << discriminant(f) << "\n";

	return(0);
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
