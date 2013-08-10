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
//	Author	: Oliver Morsch (OM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	int i = 1;
	bigint a, p;
	bigint r;
	timer T;

	std::cout << "Computing Square Roots r of a mod p, \n";
	std::cout << "where p is a prime and 0 <= a < p.\n";
	std::cout.flush ();

	while (i <= 5) {
		std::cout << "\n\n";
		std::cout << "RESSOL: please enter a: ";
		std::cin >> a;
		std::cout << "        please enter p: ";
		std::cin >> p;
		i++;

		if (jacobi(a, p) != 1)
			std::cout << "\n\n a is quadratic non residue mod p";
		else {
			T.start_timer();
			ressol_p (r, a, p);
			T.stop_timer();

			std::cout << "        square root of " << a << " mod " << p;
			std::cout << " is " << r << "\n";
			std::cout << "        Time: " << T;
		}
	}

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
