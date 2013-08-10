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


#include	"LiDIA/polynomial.h"
#include	"LiDIA/bigcomplex.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	long prec;
	std::cout << "Enter the working precision : ";
	std::cin >> prec;
	bigfloat::set_precision(prec);

	polynomial< bigcomplex > f;
	std::cout << "Enter a polynomial :\n";
	std::cin >> f;

	timer t;

	std::cout << "The roots of f(x) : = " << f << "\nare :\n";
	t.start_timer();
	base_vector< bigcomplex > z = roots(f);
	t.stop_timer();
	for (int i = 0; i < z.size(); i++) {
		std::cout << "%" << i + 1 << ": " << z[i] << std::endl;
		std::cout << "  eval f(%" << i + 1 << ") : " << f(z[i]) << std::endl << std::endl;
	}
	std::cout << "The computation took " << t.user_time() << " hs.\n";

#if 0
	std::cout << "Alternative computation with the algorithm from Cohen's book :\n";
	bigcomplex *zz = new bigcomplex[f.degree()];
	t.start_timer();
	roots(f, zz);
	t.stop_timer();
	std::cout << "The roots are \n" << base_vector< bigcomplex > (zz, f.degree())
		  << std::endl;
	std::cout << "The computation took " << t.user_time() << " hs.\n";
	delete[] zz;
#endif
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
