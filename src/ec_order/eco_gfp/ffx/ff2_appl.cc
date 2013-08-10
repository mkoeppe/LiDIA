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
//	Author	:Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/ff2.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	udigit p;
	udigit exp;
	ff2 a, b, c;
	base_vector< ff2 > v;

	std::cout << " Characteristic of GF(l^2) : "; std::cin >> p; std::cout << "\n" << std::endl;
	ff2::init_field(p);

	std::cout << " Input a : "; std::cin >> a; std::cout << std::endl;
	std::cout << " Input b : "; std::cin >> b; std::cout << std::endl;

	std::cout << " Read a = " << a << std::endl;
	std::cout << " Read b = " << b << std::endl;

	add      (c, a, b); std::cout << " a + b = " << c << std::endl;
	subtract (c, a, b); std::cout << " a - b = " << c << std::endl;
	multiply (c, a, b); std::cout << " a * b = " << c << std::endl;
	divide   (c, a, b); std::cout << " a / b = " << c << std::endl;
	square   (c, a); std::cout << " a ^ 2 = " << c << std::endl;
	//  sqrt	   (b, c);
	negate   (c, a); std::cout << " - a   = " << c << std::endl;
	invert   (c, a); std::cout << " 1 / a = " << c << std::endl;
	power    (c, a, 3); std::cout << " a ^ 3 = " << c << std::endl;


	std::cout << " order of element to be determined : "; std::cin >> exp; std::cout << "\n" << std::endl;

	nearly_all_of_order (exp, v);
	std::cout << " nearly all elements of order " << exp << " = " << v << std::endl;

	std::cout << " multiplicative order of " << c << " = " << c.multiplicative_order() << std::endl;

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
