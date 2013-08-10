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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	std::cout << "Simple demonstration of the class polynomial< gf_element > \n";
	std::cout << "(an implementation of polynomials over finite fields)\n\n";
	std::cout << "Please enter the finite field in form of \"[p, n]\", where p is the\n";
	std::cout << "characteristic and n the degree of the finite field : ";
	galois_field A;
	std::cin >> A;

	lidia_size_t d = 0;

	while (d < 1) {
		std::cout << "Please enter the degree of the polynomials: ";
		std::cin >> d;
	}

	bigint p = A.characteristic();
	std::cout << "multiplication of two polynomials modulo another polynomial\n";
	std::cout << "finite field:\n";
	std::cout << "  - characteristic " << p << " (" << p.bit_length() << " bits)\n";
//    std::cout<<"  - generating polynomial "<<A.irred_polynomial()<<std::endl;
	std::cout << "  - degree " << A.degree() << " over prime field\n";
	std::cout << "degree of the polynomials = " << d << "\n\n";

	std::cout << "generating random polynomials..." << std::flush;
	polynomial< gf_element > f, g, h, r1, r2, r3;
	gf_element one;
	one.assign_one(A);
	f = randomize(A, d);
	f[d] = one;

	timer t;
	gf_poly_modulus F(f);
	g = randomize(A, d-1);
	h = randomize(A, d-1);
	std::cout << " done.\n\n";

	std::cout << "classical algorithm... \t\t\t" << std::flush;
	t.start_timer();
	plain_multiply(r1, g, h);
	plain_remainder(r1, r1, f);
	t.stop_timer();
	std::cout << t.user_time() << " hsec\n";

	std::cout << "Kronecker substitution algorithm... \t" << std::flush;
	t.start_timer();
	multiply(r2, g, h);
	remainder(r2, r2, f);
	t.stop_timer();
	std::cout << t.user_time() << " hsec\n";

	std::cout << "multiplication with gf_poly_modulus... \t" << std::flush;
	t.start_timer();
	multiply(r3, g, h, F);
	t.stop_timer();
	std::cout << t.user_time() << " hsec\n";

	std::cout << std::endl;

	// compare the results...

	if (r1 != r2)
		std::cout << "r1 != r2!!\n";
	else if (r1 != r3)
		std::cout << "r1 != r3!!\n";
	else
		std::cout << "results are OK.\n";

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
