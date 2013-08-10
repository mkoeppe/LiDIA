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
//	Author	: Thomas Pfahler(TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



typedef Fp_polynomial elem_type;
//typedef long elem_type;
//typedef bigint elem_type;



int main_LiDIA(int argc, char** argv)
{
	int err_num = 0;

	std::cout << "\n\n This (very simple) program runs several tests to verify\n";
	std::cout << " the correctness of the template class factorization< Fp_polynomial > .\n";
	std::cout << " If an error occurs, a corresponding message\n";
	std::cout << " is displayed. At the end, the program reports\n";
	std::cout << " the number of detected errors.\n\n" << std::endl;


	factorization< elem_type > F1;

	std::cout << " Please enter a factorization with at least two entries\n";
	std::cout << " (format: [ [base1, exp1], [base2, exp2], ... [basen, expn] ] :\n\n";
	std::cin >> F1;
	std::cout << "\n\n" << F1 << "\n" << std::endl;

	std::cout << "Please, check whether the output for a is correct." << std::endl;
	std::cout << "==================================================\n" << std::endl;




	factorization< elem_type > F2(F1), F3;
	if (F1 != F2) {
		std::cout << "factorization(factorization) or operator != failed\n" << std::endl;
		err_num++;
	}

	F2.kill();
	if (F2 != F3) {
		std::cout << "factorization() or kill() or operator != failed\n" << std::endl;
		err_num++;
	}

	F2 = F1;
	if (F1 != F2) {
		std::cout << "operator = or operator != failed\n" << std::endl;
		err_num++;
	}



	single_factor< elem_type > s1;
	if (F1.no_of_components() > 0) {
		s1 = F1.composite_base(0);
		F2.kill();
		F2 = s1;
		factorization< elem_type > F4(s1);
		if (F2 != F4) {
			std::cout << "factorization(single_factor) or operator != failed\n" << std::endl;
			err_num++;
		}
	}

	elem_type t1;
	t1 = F1.unit();
	single_factor< elem_type > s2(t1), s3;
	s3 = t1;
	if (s2 != s3) {
		std::cout << "single_factor(single_factor) or operator = failed\n" << std::endl;
		err_num++;
	}


	std::cout << "\nNumber of errors : " << err_num << std::endl;
	if (err_num == 0)
		std::cout << "OK.\n";

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
