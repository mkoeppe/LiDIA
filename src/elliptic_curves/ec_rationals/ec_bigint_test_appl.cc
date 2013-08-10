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
	long lidia_precision = 20;
	bigfloat::set_precision(lidia_precision);

	elliptic_curve< bigint > ER;
	bigint prog_N, file_N;
	long   prog_T, file_T;

	long ncurves = 0, nbadN = 0, nbadT = 0;
	std::cout << "Starting test, output . after every 100 curves\n";

	while (std::cin >> ER >> file_N >> file_T, std::cin.good()) {
		ncurves++;
		if ((ncurves%100) == 0)
			std::cout << "." << std::flush;
		prog_N = ER.get_conductor();
//ER.output_short(std::cout); std::cout<<std::endl;
//ER.output_torsion(std::cout);
		prog_T = ER.get_no_tors();
		if (file_N != prog_N) {
			nbadN++;
			std::cout << "Curve " << ER << " wrong conductor (";
			std::cout << prog_N << " computed, should be " << file_N << ")" << std::endl;
		}
		if (file_T != prog_T) {
			nbadT++;
			std::cout << "Curve " << ER << " wrong torsion order (";
			std::cout << prog_T << " computed, should be " << file_T << ")" << std::endl;
		}
	}
	std::cout << ncurves << " curves tested\n";
	if (nbadN) std::cout << nbadN << " wrong conductors\n";
	else std::cout << "All conductors correct\n";
	if (nbadT) std::cout << nbadT << " wrong torsion orders\n";
	else std::cout << "All torsion orders correct\n";

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
