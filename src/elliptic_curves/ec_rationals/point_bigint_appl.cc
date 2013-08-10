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
//	Author	: Nigel Smart, John Cremona
//                Adaption of John Cremona's code; some code came
//                originally from Oisin McGuiness
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/point_bigint.h"
#include	"LiDIA/quartic.h"

#include        <cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	long lidia_precision = 40;
	bigfloat::set_precision(lidia_precision);

	elliptic_curve< bigint > ER;
	std::cout << "Enter a curve: ";
	std::cin >> ER;

	std::cout << ER << std::endl;

	point< bigint > P(ER), Q(ER), H(ER);
	std::cout << "Enter point on this curve: [x, y, z] or (x, y) \n" << std::flush;
	std::cin >> P;
	if (P.on_curve())
		std::cout << " The point " << P << " is on the curve\n\n";
	else
		lidia_error_handler("ec", "point not on curve");

	bigfloat hP = P.get_height(), hQ;
	std::cout << "P = " << P << ", h(P) = " << hP;
	long oP = P.get_order();
	if (oP > 0) std::cout << ", order = " << oP << std::endl;
	else     std::cout << ", infinte order" << std::endl;

	long i;
	for (i = 1; i <= 10; i++) {
		add(Q, Q, P);
		if (!Q.on_curve()) {
			std::cerr << "\n Q = " << i << " * P = " << Q << " is not on curve any more\n";
			std::exit(1);
		}
		multiply(H, i, P);

		if (Q != H) {
			std::cerr << "\n Q = " << i << " * P = " << Q;
			std::cerr << "\n H = " << H << " is different\n";
			std::exit(1);
		}
		hQ = Q.get_height();
		std::cout << i << "  " << Q << "\n     " << hQ << "\n     ";
		std::cout << i*i*hP;
		std::cout << "\n     " << Q.get_naive_height();
		std::cout << std::endl;
	}
	std::cout << "\n found no errors\n\n";

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
