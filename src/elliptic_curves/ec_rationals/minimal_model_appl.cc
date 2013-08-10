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
//	Author	: Nigel Smart (NS)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/base_vector.h"
#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/point_bigint.h"
#include	"LiDIA/curve_isomorphism.h"
#include	"LiDIA/minimal_model.h"
#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	long lidia_precision = 40;
	bigfloat::set_precision(lidia_precision);


	elliptic_curve< bigrational > ER;
	elliptic_curve< bigint > Emin;
	curve_isomorphism< bigrational, bigint > iso;

	std::cout << "Enter a curve: ";
	std::cin >> ER;
	std::cout << "Input curve E0 = " << ER << std::endl;

	minimal_model(Emin, ER, iso);

	std::cout << "Minimal model E1 = " << Emin << std::endl;
	std::cout << "Isomorphism E0->E1: [u, r, s, t] = " << iso << std::endl;
	std::cout << "Unimodular?...";
	if (iso.is_unimodular()) std::cout << "yes"; else std::cout << "no";
	std::cout << std::endl;
	std::cout << "Identity?...";
	if (iso.is_identity()) std::cout << "yes"; else std::cout << "no";
	std::cout << std::endl;

	std::cout << "Long output:\n";
	Emin.output_long(std::cout);
	std::cout << std::endl;
	std::cout << "Torsion points on the minimal model: \n";
	Emin.output_torsion(std::cout);
	std::cout << std::endl;

	std::cout << "Torsion points on the input model: \n";
	base_vector< point < bigint > > tor = Emin.get_torsion();
	long it, nt = tor.get_size();
	for (it = 0; it < nt; it++) {
		point< bigrational > P = iso.inverse(tor[it]);
		std::cout << P;
		if (!P.on_curve()) std::cout << " -- error! not on curve!";
		std::cout << std::endl;
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
