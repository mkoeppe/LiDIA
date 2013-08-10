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
//	Author	: Andreas M"uller (AM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint n;
	int a;

	std::cout << "\nEnter the number to test for primality: ";
	std::cin >> n;

	std::cout << "\nHow many tests shall I do? ";
	std::cin >> a;

	if (is_prime(n, a))
		std::cout << "\nis_prime(" << n << ", " << a << ") = true\n";
	else
		std::cout << "\nis_prime(" << n << ", " << a << ") = false\n";

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
