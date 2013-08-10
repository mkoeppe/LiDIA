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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/quadratic_order.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint D;
	quadratic_order QO;
	int vb, strat;

	if (argc > 1)
		strat = atoi(argv[1]);
	else
		strat = -1;

	std::cout << "Enter verbosity level: ";
	std::cin >> vb;
	quadratic_order::verbose(vb);

	std::cout << "Enter verification level: ";
	std::cin >> vb;

	std::cout << "Enter a quadratic discriminant: ";
	std::cin >> D;
	std::cout << "\n";
	while (D != 0) {
		if (QO.assign(D)) {
			switch (strat) {
			case 0:  QO.class_group_BJT();
				break;

			case 1:  QO.class_group_shanks();
				break;

			case 2:  QO.class_group_randexp();
				break;

			case 3:  QO.class_group_mpqs();
				break;

			case 4:  QO.class_group_siqs();
				break;

			case 5:  QO.class_group_lp();
				break;

			case 6:  QO.regulator_BJT();
				break;

			case 7:  QO.class_group_siqs(atoi(argv[2]));
				break;

			case 8:  QO.class_group_siqs(atoi(argv[2]), 0, true);
				break;

			default: QO.class_group();
			}

			QO.generators();
			QO.Lfunction();
			QO.factor_discriminant();

			bigfloat::set_precision(QO.prec);

			std::cout << QO;
			std::cout << "   ULI = " << QO.ULI() << "\n";
			std::cout << "   LLI = " << QO.LLI() << "\n";

			if ((vb) && QO.is_real() && !QO.verify_regulator())
				std::cout << "ERROR:  regulator incorrect!!!" << std::endl;

			if (!QO.verify_class_group(vb))
				std::cout << "ERROR:  class group incorrect!!!" << std::endl;
		}

		std::cout << "\nEnter a quadratic discriminant: ";
		std::cin >> D;
		std::cout << "\n";
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
