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



#include	"LiDIA/ff1.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA (int argc, char** argv)
{
	udigit p, exp;
	udigit a, b;
	ff1 A, B, C;
	base_vector< ff1 > v;
	random_generator rg;

	std::cout << "\n Test Program for class ff1 \n\n";
	std::cout << " Input Characteristic (has to be an odd prime) : "; std::cin >> p;
	ff1::set_characteristic (p);
	bigmod::set_modulus(bigint(p));
	std::cout << "\n\n Characteristic was set to " << ff1::characteristic() << std::endl;

	bigmod aa, bb, cc;
	long e;

	for (a = 1; a < p; a++)
		for (b = 1; b < p; b++) {
			aa.assign(bigint(a)); bb.assign(bigint(b));
			A.assign(a); B.assign(b);

			add (C, A, B); add(cc, aa, bb);
			std::cerr << "C == " << C << "\n";
			std::cerr << "C.as_udigit() == "
				  << C.as_udigit() << "\n";
			std::cerr << "bigint(C.as_udigit()) == "
				  << bigint(C.as_udigit()) << "\n";
			std::cerr << "cc == " << cc << "\n";
			std::cerr << "cc == bigint(C.as_udigit()) == "
				  << (cc == bigint(C.as_udigit())) << "\n";
			assert(cc == bigint(C.as_udigit()));

			subtract (C, A, B); subtract (cc, aa, bb);
			assert(cc == bigint(C.as_udigit()));

			multiply (C, A, B); multiply (cc, aa, bb);
			assert(cc == bigint(C.as_udigit()));

			divide (C, A, B); divide (cc, aa, bb);
			assert(cc == bigint(C.as_udigit()));

			square(C, A); square   (cc, aa);
			assert(cc == bigint(C.as_udigit()));

			negate(C, A); negate   (cc, aa);
			assert(cc == bigint(C.as_udigit()));

			invert(C, A); invert   (cc, aa);
			assert(cc == bigint(C.as_udigit()));

			rg >> e;
			e %= p;
			power  (C, A, e); power  (cc, aa, e);
			assert(cc == bigint(C.as_udigit()));
		}

	std::cout << " Input order of elements to be determined : ";
	std::cin >> exp;
	if((p-1) % exp != 0) {
	    std::cout << " Order has to be a divisor of p-1 = "
		      << (p-1) << "!" << std::endl;
	}
	else {
	    nearly_all_of_order (exp, v);
	    std::cout << " Nearly all elements of order " << exp
		      << " (w/o inverses) " << v << std::endl;
	}
	std::cout << " Multiplicative order of " << C
		  << " = " << C.multiplicative_order() << std::endl;

	std::cout << " Computing square roots, please enter some element : ";
	std::cin >> C;

	if (!C.is_square()) {
		std::cout << " " << C << " is not a square." << std::endl;
	}
	else {
		sqrt(A, C);
		square(B, A);
		assert(B == C);
		std::cout << " The square root of " << C << " is " << A << "."
			  << std::endl;
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
