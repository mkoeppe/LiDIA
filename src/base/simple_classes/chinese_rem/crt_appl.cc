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
//	$Id: crt_appl.cc,v 2.3 2004/06/15 10:19:26 lidiaadm Exp $
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/crt.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	int    i , n;
	bigint big, B;
	sdigit small_value;

	crt_table T;
	crt       C;

	std::cout << " Bound for absolute value : ";
	std::cin >> B;

	T.init (B);
	C.init (T);
	n = C.number_of_primes();

	std::cout << " " << n << " congruences needed \n\n";

	for (i = 0; i < n; i++) {
		std::cout << " Congruence modulo ";
		std::cout << C.get_prime(i) << " : ";
		std::cin >> small_value;
		C.combine (small_value, i);
	}
	C.get_result (big);
	std::cout << "\n Solution : " << big << std::endl;
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
