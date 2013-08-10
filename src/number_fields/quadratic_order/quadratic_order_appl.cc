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
	quadratic_order QO1, QO2(-161651), QO3(14401), QO4, QO5, QO;
	bigint D;
	int i;

	// test whether an uninitialized order is still valid
	if (QO5.is_zero())
		std::cout << "Uninitialized quadratic order is zero? true\n" << std::flush;
	else
		std::cout << "Uninitialized quadratic order is zero? false\n" << std::flush;

	// test assign function
	string_to_bigint("307807807807807525", D);
	if (!QO4.assign(D)) {
		std::cout << "ERROR:  can't assign " << D << " to a quadratic order!\n" << std::flush;
		return 1;
	}

	string_to_bigint("-400000000000000004", D);
	if (!QO5.assign(10000))
		QO5.assign(D);


	// perform tests on 4 selected discriminants
	for (i = 0; i < 4; ++i) {
		quadratic_order::verification(0);

		switch (i) {
		case 0:  QO = QO2;
			break;
		case 1:  QO = QO3;
			break;
		case 2:  QO = QO4;
			break;
		default:
			QO = QO5;
		}

		std::cout << "\n=======================================\n" << std::flush;
		std::cout << QO << "\n" << std::flush;
		if (QO.is_imaginary())
			std::cout << "imaginary? true\n" << std::flush;
		else
			std::cout << "imaginary? false\n" << std::flush;
		if (QO.is_real())
			std::cout << "real? true\n" << std::flush;
		else
			std::cout << "real? false\n" << std::flush;
		std::cout << "disriminant = " << QO.discriminant() << "\n\n" << std::flush;

		// compute the class group
		std::cout << "computing class group...\n" << std::flush;
		QO.class_group();
		std::cout << "\n" << QO << "\n" << std::flush;

		// test L(s,chi)-related functions
		std::cout << "testing LD, ULI, LLI, etc...\n" << std::flush;
		std::cout << "L(3) (p < 10000, truncated prod) = " << QO.estimate_L(3, 10000) << "\n" << std::flush;
		std::cout << "L(2) (p < 10000, truncated prod) = " << QO.estimate_L(2, 10000) << "\n" << std::flush;
		std::cout << "L(1) (p < 10000, truncated prod) = " << QO.estimate_L(1, 10000) << "\n" << std::flush;
		std::cout << "L(1) (p < 10000, Bach) = " << QO.estimate_L1(5000) << "\n" << std::flush;
		std::cout << "approx. error = " << QO.estimate_L1_error(5000) << "\n" << std::flush;
		std::cout << "L(1) = " << QO.Lfunction() << "\n" << std::flush;
		std::cout << "LLI = " << QO.LLI() << "\n" << std::flush;
		std::cout << "ULI = " << QO.ULI() << "\n" << std::flush;
		std::cout << "LD(1) = " << QO.LDfunction() << "\n" << std::flush;
		std::cout << "LLI_D = " << QO.LLI_D() << "\n" << std::flush;
		std::cout << "ULI_D = " << QO.ULI_D() << "\n" << std::flush;
		std::cout << "C(D) = " << QO.Cfunction() << "\n\n" << std::flush;

		// test factor discriminant function
		std::cout << "testing factor_discriminant...\n" << std::flush;
		std::cout << QO.discriminant() << " = " << QO.factor_discriminant() << "\n\n" << std::flush;

		// test factor_h function
		std::cout << "testing factor_h...\n" << std::flush;
		std::cout << QO.class_number() << " = " << QO.factor_h() << "\n\n" << std::flush;

		// test class group related functions
		std::cout << "CL = " << QO.class_group() << "\n" << std::flush;
		std::cout << "generators = " << QO.generators() << "\n" << std::flush;
		std::cout << "exponent = " << QO.exponent() << "\n" << std::flush;
		std::cout << "2-rank = " << QO.p_rank(bigint(2)) << "\n" << std::flush;
		std::cout << "3-rank = " << QO.p_rank(bigint(3)) << "\n\n" << std::flush;

		// test maximize functions
		std::cout << "testing maximize..." << "\n" << std::flush;
		std::cout << "QO.conductor() = " << QO.conductor() << "\n" << std::flush;
		if (QO.is_maximal())
			std::cout << "maximal? true\n" << std::flush;
		else
			std::cout << "maximal? false\n" << std::flush;
		QO.maximize(QO1);
		std::cout << "maximal order: ";
		std::cout << QO1 << "\n" << std::flush;

		// test comparisions
		if (QO.is_equal(QO1))
			std::cout << "equal to max ord? true\n" << std::flush;
		else
			std::cout << "equal to max ord? false\n" << std::flush;
		if (QO == QO1)
			std::cout << "equal (operator)? true\n" << std::flush;
		else
			std::cout << "equal (operator)? false\n" << std::flush;
		if (QO < QO1)
			std::cout << "proper subset of max ord? true\n" << std::flush;
		else
			std::cout << "proper subset of max ord? false\n" << std::flush;
		if (QO <= QO1)
			std::cout << "subset of max ord? true\n\n" << std::flush;
		else
			std::cout << "subset of max ord? false\n\n" << std::flush;

		// test verification functions
		quadratic_order::verification(2);
		if (QO.is_real()) {
			std::cout << "Verifying regulator:\n" << std::flush;
			if (QO.verify_regulator())
				std::cout << "Regulator is correct!\n\n" << std::flush;
			else {
				std::cout << "ERROR:  regulator incorrect!\n\n" << std::flush;
				return 1;
			}
		}

		std::cout << "Verifying class group:\n" << std::flush;
		if (QO.verify_class_group())
			std::cout << "Class group is correct!\n\n" << std::flush;
		else {
			std::cout << "ERROR:  class group incorrect!\n\n" << std::flush;
			return 1;
		}
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
