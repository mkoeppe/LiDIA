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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf2n_polynomial.h"
#include	"LiDIA/timer.h"

#include        <cstdlib>

#undef PRINT



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	unsigned int degree;
	int low, high;
	timer t1, t2;

	t1.set_print_mode();
	t2.set_print_mode();

	std::cout << "\nPlease input extension degree of field : "; std::cin >> degree;

	gf2n_init(degree);

	std::cout << "\nInput lower bound deg(f) = "; std::cin >> low;
	std::cout << "\nInput upper bound deg(f) = "; std::cin >> high;

	gf2n_polynomial f(1), g(1), h1(100), h2(100);
	gf2n_polynomial a0, a1, b0, b1;


	for (int i = low; i <= high; i++) {
		std::cout << "\n\nTests for deg " << i << std::flush;
		int j = i-1;
		for (int t = 0; t < 10; t++) {
			f.randomize(i);

			t1.cont_timer();
			plain_inv(h1, f, j);
			t1.stop_timer();
			t2.cont_timer();
			newton_inv(h2, f, j);
			t2.stop_timer();

			if (h1 != h2) {
				std::cout << "\ncorrect h1 = " << h1;
				std::cout << "\nwrong   h2 = " << h2;
				std::cout << "\nERROR";
				std::exit(1);
			}
		}
		std::cout << " successful: \ntime(plain) = " << t1;
		std::cout << "\ntime(newton) = " << t2 << std::flush;
		t1.start_timer(); t1.stop_timer();
		t2.start_timer(); t2.stop_timer();
	}

	std::cout << std::endl;
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
