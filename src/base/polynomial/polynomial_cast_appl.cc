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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/polynomial.h"
#include	"LiDIA/bigcomplex.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	polynomial< bigint > A, B;

	polynomial< bigcomplex > C;

	std::cout << "\n input a polynomial< bigint > : " << std::flush;
	std::cin >> A;
	std::cout << std::endl;
	std::cout << " input another polynomial< bigint > : " << std::flush;
	std::cin >> B;

	std::cout << "\n Number of real roots of " << A << " is "
		  << no_of_real_roots(A) << std::flush;
	std::cout << "\n gcd is " << gcd(A, B) << std::flush;
	std::cout << "\n gcd over bigcomplex is " << gcd(polynomial< bigcomplex > (A),
							 polynomial< bigcomplex > (B));
	std::cout << std::flush;

	bigcomplex * zeros = new bigcomplex[A.degree()+1];

	roots (A, zeros);
	std::cout << "\n Zeros of " << A << "are :" << std::endl;
	for (long int i = 0; i < A.degree(); i++)
		std::cout << "\n\t" << zeros[i];
	std::cout << "\n input a polynomial< bigcomplex > : " << std::flush;
	std::cin >> C;
	std::cout << "(" << A << ") * (" << C << ") = " << base_polynomial< bigcomplex > (A) * C << std::endl;

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
