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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/mv_term.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	int n;

	bigint v1, v2;

	std::cout << "\nmv_term tests";
	std::cout << "\n=============\n\n";
	std::cout << "Input Extension degree of Field : ";
	std::cin >> n;
	gf2n_init(n);

	mv_term a, b, c;
	gf2n bb;

	std::cout << "\nInput element in field : ";
	std::cin >> bb;

	std::cout << "\nInput Variables in Bigint form : ";
	std::cin >> v1;

	a.assign(bb, v1);

	std::cout << "\n\nYou entered this term a : " << a << std::endl;

	std::cout << "\n b = " << b;


	square(b, a);
	std::cout << "Square(a) = " << b << std::endl;
	c = a + b;
	std::cout << "\na + square(a) = " << c << std::endl;

	std::cout << "\n\nPlease enter new variables : ";
	std::cin >> v2;

	a.assign_var(v2);
	std::cout << "\n We test with a =  " << a;
	std::cout << "\n and          b =  " << b << std::endl;

#if 0
	c = a * b;
	std::cout << a << " * " << b << " = " << c << std::endl;
	if (a == b)
		std::cout << a << " == " << b << std::endl;
	if (a != b)
		std::cout << a << " != " << b << std::endl;
#endif

	std::cout << "\nInput index of variable from a that will be substituted by b : ";
	std::cin >> n;
	substitute(a, b, n);
	std::cout << "\nResult of substitution : " << a << std::endl;

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
