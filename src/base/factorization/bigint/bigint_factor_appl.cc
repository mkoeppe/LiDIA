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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/rational_factorization.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	factorization< bigint > res;
	single_factor< bigint > g;
	rational_factorization rf;
	bigint n;
	int choice;
	timer t;

	t.set_print_mode(HMS_MODE);
	single_factor< bigint >::set_verbose_flag(1);

	std::cout << "\nPlease enter a number "; std::cin >> n;
	std::cout << "\n Enter factorization method: \n\n (1) Trial Division";
	std::cout << "\n (2) Pollard p-1 \n (3) Williams p+1";
	std::cout << "\n (4) Pollard Rho \n (5) Fermat's Method";
	std::cout << "\n (6) ECM\n (7) MPQS\n (8) Built-in Strategy\n";

	std::cout << "\n Please enter your Choice : "; std::cin >> choice; std::cout << "\n";
	t.start_timer();

	switch (choice) {
	case 1:
		res = TrialDiv(n);
		break;
	case 2:
		res = PollardPminus1 (n, 15);
		break;
	case 3:
		res = WilliamsPplus1 (n, 15);
		break;
	case 4:
		res = PollardRho(n, 10);
		break;
	case 5:
		res = Fermat(n);
		break;
	case 6:
		res = ECM(n);
		break;
	case 7:
		res = MPQS(n);
		break;
	case 8:
		res = sf_factor(n);
		break;
	default:
		std::cout << "\n\n Sorry, bad choice, abort. \n";
		return 1;
	}
	t.stop_timer();

	std::cout << "\nResult is : " << res;
	if (res.is_prime_factorization())
		std::cout << "  is prime factorization\n";
	else
		std::cout << "  is no prime factorization\n";

	std::cout << "\nUsed Time for this Factorization : " << t << "\n";

	t.start_timer();
	rf = n;
	rf.factor();
	t.stop_timer();

	std::cout << "\n\nRational factorization gives result : " << rf;
	std::cout << "\nNeeded Time : " << t << "\n\n";

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
