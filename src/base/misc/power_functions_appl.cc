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
//	Author	: Markus Maurer (MM), Stefan Neis (SN)
//		  Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/power_functions.h"
#include	"LiDIA/timer.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{

	std::cout << "Testing repeated squaring left_to_right and right_to_left" << std::endl;
	std::cout << "and sliding window, window size 1-5.\n" << std::endl << std::endl;

	std::cout << "Please enter the magnitude of the modulus (the modulus will" << std::endl;
	std::cout << "be randomly chosen and will be less than 10^< your input > ):" << std::endl;

	bigint tmp;
	std::cin >> tmp;
	lidia_power(tmp, bigint(10), tmp);
	tmp = randomize (tmp);
	bigmod::set_modulus (tmp);

	bigmod powl, powr, pows[5];
	bigmod a;
	unsigned long exp;
	bigint bigint_exp;

	random_generator rg;

	a.randomize();
	//std::cin >> a;
	//std::cout << a << std::endl;

	// exp randomly chosen
	std::cout << "Please enter the magnitude of the exponent: " << std::endl;

	std::cin >> bigint_exp;
	lidia_power(bigint_exp, bigint(10), bigint_exp);
	bigint_exp = randomize (bigint_exp);

	//  rg >> exp;

	//exp %= 100000;
	//std::cin >> exp;
	//std::cout << exp << std::endl;
	remainder(exp, bigint_exp, 100000);

	std::cout << "Base: \t" << a << std::endl;
	std::cout << "Exponent:\t" << exp << std::endl;
	std::cout << "Big exponent:\t" << bigint_exp << std::endl;
	std::cout << "Modulus:\t" << tmp << std::endl;
	std::cout << "\n";
	std::cout << "This may take a while ...\n";
	std::cout << "\n";
	std::cout.flush();


	// bigint_exp = exp;
	timer T1, T2, TS1, TS2, TS3, TS4, TS5;

	std::cout << "Testing l2r" << std::endl;

	T1.start_timer();
	lidia_power_right_to_left(powr, a, bigint_exp);
	T1.stop_timer();

	std::cout << "Testing r2l" << std::endl;

	T2.start_timer();
	lidia_power_left_to_right(powl, a, bigint_exp);
	T2.stop_timer();

	std::cout << "Testing pws1" << std::endl;

	TS1.start_timer();
	lidia_power_sliding_window(pows[0], a, bigint_exp, 1);
	TS1.stop_timer();

	std::cout << "Testing pws2" << std::endl;

	TS2.start_timer();
	lidia_power_sliding_window(pows[1], a, bigint_exp, 2);
	TS2.stop_timer();

	std::cout << "Testing pws3" << std::endl;

	TS3.start_timer();
	lidia_power_sliding_window(pows[2], a, bigint_exp, 3);
	TS3.stop_timer();

	std::cout << "Testing pws4" << std::endl;

	TS4.start_timer();
	lidia_power_sliding_window(pows[3], a, bigint_exp, 4);
	TS4.stop_timer();

	std::cout << "Testing pws5" << std::endl;

	TS5.start_timer();
	lidia_power_sliding_window(pows[4], a, bigint_exp, 5);
	TS5.stop_timer();

	std::cout << "power_right_to_left(bigint): " << T1 << std::endl;
	std::cout << "power_left_to_right(bigint): " << T2 << std::endl;
	std::cout << "power_sliding_window, k = 1: " << TS1 << std::endl;
	std::cout << "power_sliding_window, k = 2: " << TS2 << std::endl;
	std::cout << "power_sliding_window, k = 3: " << TS3 << std::endl;
	std::cout << "power_sliding_window, k = 4: " << TS4 << std::endl;
	std::cout << "power_sliding_window, k = 5: " << TS5 << std::endl;
	//std::cout << "powl == " << powl << std::endl;
	//std::cout << "powr == " << powr << std::endl;

	// verify result
	assert(powl == powr);
	assert(powl == pows[0]);
	for (int i = 0; i < 4; i++)
		assert(pows[i] == pows[i+1]);

	T1.start_timer();
	lidia_power_right_to_left(powr, a, exp);
	T1.stop_timer();

	// verify result
	if (bigint_exp == bigint(exp))
		assert(powl == powr);

	T2.start_timer();
	lidia_power_left_to_right(powl, a, exp);
	T2.stop_timer();

	std::cout << "power_right_to_left(long): " << T1 << std::endl;
	std::cout << "power_left_to_right(long): " << T2 << std::endl;
	//std::cout << "powl == " << powl << std::endl;
	//std::cout << "powr == " << powr << std::endl;

	// verify result
	assert(powl == powr);

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
