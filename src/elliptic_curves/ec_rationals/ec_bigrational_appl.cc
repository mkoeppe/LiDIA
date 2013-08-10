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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/bigrational.h"

#include        <cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	std::cout << "\n ELLIPTIC CURVES OVER RATIONALS\n\n";

	bigrational a1, a2, a3, a4, a6, x, y;

	std::cout << "\n Input of curve:\n\n";
	std::cout << "\n Input a1 : "; std::cin >> a1;
	std::cout << "\n Input a2 : "; std::cin >> a2;
	std::cout << "\n Input a3 : "; std::cin >> a3;
	std::cout << "\n Input a4 : "; std::cin >> a4;
	std::cout << "\n Input a6 : "; std::cin >> a6;
	elliptic_curve< bigrational > e(a4, a6);

	std::cout << "\n We work with elliptic curve " << e;
	std::cout << "\n\n Input point P:";
	std::cout << "\n x(P) = "; std::cin >> x;
	std::cout << "\n y(P) = "; std::cin >> y;
	point< bigrational > P(x, y, e);
	point< bigrational > Q (e), H (e);

	if (P.on_curve())
		std::cout << "\n Given P = " << P << " is point on curve\n\n";
	else
		lidia_error_handler("", "point not on curve");

	unsigned long i, na;

	std::cout << "\n Number of additions: "; std::cin >> na; std::cout << "\n";

	for (i = 1; i <= na; i++) {
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
		std::cout << "\n " << i << " * P = " << Q;
	}
	std::cout << "\n\n Found no errors during " << na << " additions\n\n";

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
