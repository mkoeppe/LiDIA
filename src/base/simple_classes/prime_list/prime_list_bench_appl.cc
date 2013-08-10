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


#include	"LiDIA/prime_list.h"
#include	"LiDIA/timer.h"
#include	<cmath>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA (int, char**)
{
	int i;
	PRIME_LIST_NUMBER lb, ub;
	prime_list *pl;
	timer t;
	char *name = NULL, mode = '\0';
	PRIME_LIST_COUNTER mem, su;
	lidia_size_t np;

	std::cout << "Please enter lower bound: " << std::flush;
	std::cin >> lb;
	std::cout << "Please enter upper bound: " << std::flush;
	std::cin >> ub;

	std::cout << "\n";
	std::cout << "*******************************************************************\n";
	std::cout << "**                Benchmark-Test for prime_list                  **\n";
	std::cout << "*******************************************************************\n";
	std::cout << "\n";
	std::cout << "Mode  Algorithm                         Memory* (bytes)  Time (sec)\n";
	std::cout << "-------------------------------------------------------------------\n";
	for (i = 0; i < 5; i++) {
		su = static_cast<PRIME_LIST_COUNTER>(std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(ub)));
		mem = su;
		if (static_cast<PRIME_LIST_NUMBER>(mem) > lb)
			mem = lb;
		mem += (ub - lb);
		switch (i) {
		case 0:
			name = "Sieve of Erathostenes";
			mode = 'E';
			mem *= sizeof(PRIME_LIST_SIEVE);
			break;
		case 1:
			name = "6k+-1 Sieve";
			mode = 'K';
			mem /= 3;
			mem *= sizeof(PRIME_LIST_SIEVE);
			break;
		case 2:
			name = "Sieve of Erathostenes (Bit-Level)";
			mode = 'B';
			mem /= 8;
			break;
		case 3:
			name = "6k+-1 Bit Sieve";
			mode = '6';
			mem /= 24;
			break;
		case 4:
			name = "Interval Sieve";
			mode = 'I';
			mem = ub - lb;
			if (su > mem)
				mem = su;
			if (mem > 1000000)
				mem = 1000000;
			mem *= sizeof(PRIME_LIST_SIEVE);
			mem += static_cast<PRIME_LIST_COUNTER>(static_cast<PRIME_LIST_FLOAT_NUMBER>(su) /
							       std::log(static_cast<PRIME_LIST_FLOAT_NUMBER>(su)))
				* sizeof(PRIME_LIST_COUNTER);
			break;
		}
		std::cout << "'" << mode << "'   ";
		std::cout.setf(std::ios::left, std::ios::adjustfield);
		std::cout.width(35);
		std::cout << name;
		std::cout.setf(std::ios::right, std::ios::adjustfield);
		std::cout.width(12);
		std::cout << mem << std::flush;
		t.start_timer();
		pl = new prime_list(lb, ub, mode);
		t.stop_timer();
		np = pl->get_number_of_primes();
		delete pl;
		std::cout.width(12);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout.precision(2);
		std::cout << static_cast<float>(t.user_time()) / 100 << "\n";
	}
	std::cout << "\n";
	std::cout << "Number of primes in interval: " << np << "\n";
	std::cout << "\n";
	std::cout << "* approximate amount of memory required for calculation" << std::endl;

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
