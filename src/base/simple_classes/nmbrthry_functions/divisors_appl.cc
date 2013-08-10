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


#include	"LiDIA/nmbrthry_functions.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint n;

	std::cout << "Please enter a bigint n = ";
	std::cin >> n;

	rational_factorization f(n);
	f.factor();

	std::cout << "All positive divisors of " << n << ": " << divisors(n);
	std::cout << "\n";

	std::cout << "All positive and negative divisors of " << n << ": ";
	std::cout << all_divisors(n);
	std::cout << "\n";

	std::cout << "All positive square-free divisors of " << n << ": ";
	std::cout << square_free_divisors(n);
	std::cout << "\n";

	std::cout << "All positive and negative square-free divisors of " << n << ": ";
	std::cout << all_square_free_divisors(n);
	std::cout << "\n";

	std::cout << "All positive divisors of " << n << " whose square divides " << n << " : ";
	std::cout << square_divides_n_divisors(n);
	std::cout << "\n";

	std::cout << "All positive and negative divisors of " << n << " whose square divides " << n << " : ";
	std::cout << all_square_divides_n_divisors(n);
	std::cout << "\n";

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
