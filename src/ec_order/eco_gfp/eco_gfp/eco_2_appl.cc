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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/eco_prime.h"
#include	"LiDIA/trace_list.h"
#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	lidia_size_t digits, run;
	lidia_size_t step_size;
	lidia_size_t lower, upper;
	lidia_size_t runs_per_digit;
	int info;

	std::cout << "Test program for point counting routines." << std::endl;
	std::cout << "This programs randomly generates primes p, where the number of" << std::endl;
	std::cout << "decimal digits (dd) of p runs from lower to upper in steps of stepsize." << std::endl;
	std::cout << "For each number of dd, the program generates two primes p, " << std::endl;
	std::cout << "counts the number of points on a random E, and " << std::endl;
	std::cout << "verifies the correctness of the result with a probablistic test." << std::endl;
	std::cout << std::endl;

	std::cout << "Please, enter the lower bound for the number of dd of p (0 = default 10): ";
	std::cin  >> lower;
	std::cout << "Please, enter the upper bound for the number of dd of p (0 = default 30): ";
	std::cin  >> upper;
	std::cout << "Please, enter the stepsize (0 = default 5): ";
	std::cin  >> step_size;
	std::cout << "Please, enter the runs per digit (0 = default 3): ";
	std::cin  >> runs_per_digit;
	std::cout << "Please, set info mode: (0 = only parameters and group order, 1 = verbose): ";
	std::cin  >> info;

	if (lower <= 0)
		lower = 10;
	if (upper <= 0)
		upper = 30;
	if (step_size <= 0)
		step_size = 5;
	if (runs_per_digit <= 0)
		runs_per_digit = 3;
	if (info == 1)
		info = 100;
	else if (info != 0)
		info = 0;

	std::cout << std::endl;
	std::cout << "Doing tests for p with number of dd from " << lower << " to " << upper;
	std::cout << " with stepsize " << step_size << " and " << runs_per_digit << " tests per digit. " << std::endl;

	if (info == 0)
		std::cout << "Outputting only parameters and group order." << std::endl;
	else
		std::cout << "Verbose output." << std::endl;

	std::cout << std::endl << "This may take a while .... " << std::endl;

	timer t;
	t.set_print_mode (HMS_MODE);

	bigint p;
	bigint start, order;

	for (digits = lower; digits < upper; digits += step_size) {
		for (run = 0; run < runs_per_digit; run++) {
			if (run == 0) {
				power(start, 10, digits-1);
				p = next_prime(start);
			}
//			else
//				p = next_prime(p+2);

			std::cout <<"\nCharacteristic : "; std::cout << p;
			std::cout << " (" << decimal_length(p) << ")" << std::flush;
			bigmod::set_modulus(p);

			galois_field theField(p);
			gf_element a4(theField), a6(theField);

			do {
				if (run == 0)
					a4.assign_zero();
				else
					a4.randomize();

				if (run == 1)
					a6.assign_zero();
				else
					a6.randomize();
			} while ((4*a4*a4*a4 + 27*a6*a6) == 0); // Choose new curve

			std::cout << "\nCoefficient a4 : "; std::cout << a4 << std::endl;
			std::cout << "Coefficient a6 : "; 	std::cout << a6 << std::endl;

			t.start_timer();
			order = compute_group_order(a4, a6, info);
			t.stop_timer();

			std::cout << "Group order = " << order << std::endl; 
			t.print();
			std::cout << std::endl;
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
