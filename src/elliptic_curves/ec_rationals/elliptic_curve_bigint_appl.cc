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


#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/point_bigint.h"
#include	"LiDIA/quartic.h"



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

	std::cout << "Default output:  \t"; std::cout << ER << std::endl;
	std::cout << "output_short():  \t"; ER.output_short(std::cout); std::cout << std::endl;
	std::cout << "output_long():   \t"; ER.output_long(std::cout); std::cout << std::endl;
	std::cout << "output_pretty(): \t"; ER.output_pretty(std::cout); std::cout << std::endl;
	std::cout << "output_tex():    \t"; ER.output_tex(std::cout); std::cout << std::endl;

	ER.set_output_mode(1);
	std::cout << "Default output after set_output_mode(1): \t"; std::cout << ER << std::endl;
	ER.set_output_mode(2);
	std::cout << "Default output after set_output_mode(2): \t"; std::cout << ER << std::endl;
	ER.set_output_mode(3);
	std::cout << "Default output after set_output_mode(3): \t"; std::cout << ER << std::endl;
	ER.set_output_mode(0);
	std::cout << "Default output after set_output_mode(0): \t"; std::cout << ER << std::endl;

	std::cout << std::endl;
	ER.output_torsion(std::cout);
	std::cout << std::endl;

	std::cout << "Periods \n";
	std::cout << "Omega_1 = " << ER.get_omega_1() << std::endl;
	std::cout << "Omega_2 = " << ER.get_omega_2() << std::endl;

	std::cout << "Height Constant\n";
	std::cout << ER.get_height_constant() << std::endl;

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
