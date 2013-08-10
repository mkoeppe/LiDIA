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



#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/trace_list.h"
#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/timer.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	int n;
	bigint q, co_factor, found_co_factor, found_co_factor_tw, res;
	timer t, t2; t.set_print_mode (HMS_MODE); t2.set_print_mode(HMS_MODE);
	int tries = 1, percentage;
	bool good_curve = false, stop_computation;
	rational_factorization rf;
	int both_abort = 0, no_abort = 0;
	int info;
	bool no_mod_equations;

	char printCurrentTime[1024];
	sprintf(printCurrentTime, "/bin/date");

	std::cout << "\nCryptographic Strong Elliptic Curve Finder (LiDIA version)";
	std::cout << "\n==========================================================";

	std::cout << "\n\nThis program can be used to find a cryptographically strong";
	std::cout << " elliptic curve over a field of characteristic 2. These curves ";
	std::cout << " can always be chosen in the form" << std::endl;
	std::cout << "    Y^2 + X Y = X^3 + a2 X^2 + a6.\n\n" << std::flush;

	if (argc == 2) {
		if (!strcmp(argv[1], "--quiet"))
			info = 0;
		else
			info = 1;

		if (!strcmp(argv[1], "--help")) {
			std::cout << "usage: find_curve_gf2n_appl [--quiet]" << std::endl;
			exit(0);
		}
	}
	else
		info = 1;

	//
	// Determine characteristic
	//
	std::cout << "\nPlease input extension degree n : "; std::cin >> n;

	if (n < 5 || n > 1000) {
		std::cout << "\nPlease use reasonable degrees, aborting ... " << std::endl;
		exit(1);
	}

	shift_left(q, bigint(1), n);
	std::cout << "\nUsing field GF(2^" << n << ") with " << q << " elements (";
	std::cout << decimal_length(q) << ")" << std::endl;

	gf2n_init(n);

	if (n > 300) {
		std::cout << "\n\nInitializing table for solving quadratic equations. This";
		std::cout << "\nmay take some time ... " << std::flush;
		gf2n::initialize_table_for_solve_quadratic();
		std::cout << " done.\n" << std::flush;
	}

	//
	// Determine upper bound
	//

	do {
		std::cout << std::endl;
		std::cout << "Please input a percentage x. The number of decimal digits in the ";
		std::cout << "upper bound for the composite factor will be approximately ";
		std::cout << "x percent of the number of decimal digits of p. Note";
		std::cout << " that all elliptic curves of the special form Y^2 + X Y = ";
		std::cout << "X^3 + a2 X^2 + a6 are always divisible by 4 for a2 = 0 and ";
		std::cout << "2 otherwise. "<<std::endl;
		std::cout << "Percentage (0-100) = ";
		std::cin >> percentage;
	} while (percentage < 0 || percentage > 100);

	if (percentage == 100)
		co_factor = q;
	else
		if (percentage == 0)
			co_factor = 1;
		else
			power(co_factor, bigint(10), static_cast<int>((static_cast<double>(decimal_length(q)) /
								       100.0) * static_cast<double>(percentage)));

	if (co_factor > q)
		co_factor = q;

	if (co_factor < 4)
		co_factor = 4;

	std::cout << "Upper bound for the composite factor is " << co_factor;
	std::cout << " (" << decimal_length(co_factor) << ")" << std::endl;

	//
	// initialize computation
	//



	gf2n a6, a2;

	a2.assign_zero();

	elliptic_curve<gf_element> e;
	trace_list tl;
	trace_mod tm;
	udigit ll;
	long qq;
	eco_gf2n ep;
	bigint g;

	ep.set_info_mode(info);
	ep.set_schoof_bound(11);
	ep.set_atkin_bound(70);
	ep.set_power_bound(100);
	ep.set_strategy(eco_gf2n::COMPUTE_SP_DEGREE_IF_ATKIN);

	trace_list::set_info_mode(info);

	t.start_timer();

	//
	// Choose curves until good one is found.
	// Use early abort strategy on curve and twist.
	//
	do {
		// Choose next curve
		//
		std::cout << " \n\n---------------------------------------------------------------";
		a6.randomize();
		std::cout << "\nChoose random E :   a6 = " << a6 << std::endl;
		std::system(printCurrentTime);

		// initialize for new curve
		//
		tries ++;
		found_co_factor = 1;
		found_co_factor_tw = 1;
		stop_computation = false;
		no_mod_equations = false;
		e = ep.generateCurve(a6);

		ep.set_curve(a6);

		t2.start_timer();
		ll = 3;

		tl.clear();
		tl.set_curve(e);

		// trace mod 2
		//
		if (info) {
			std::cout << "\n\n-------------------------------------------------------\n";
			std::cout << "Working on prime l = 2";
		}
		ep.compute_mod_2_power();
		tm.set_vector(ep.get_prime(), ep.get_relation());
		tl.append(tm);

		qq = remainder(q, ep.get_prime());
		g = gcd(q + 1 - (ep.get_relation())[0], ep.get_prime());
		if (! g.is_one())
			multiply(found_co_factor, found_co_factor, g);
		g = gcd(q + 1 + (ep.get_relation())[0], ep.get_prime());
		if (! g.is_one())
			multiply(found_co_factor_tw, found_co_factor_tw, g);

		// for each prime ll compute trace mod ll
		// and check for divisor of group order
		//

		bool enough_relations = false;

		do {
			if (found_co_factor > co_factor && found_co_factor_tw > co_factor) {
				t2.stop_timer();
				std::cout << "\nTest No. " << tries-1 << " (early abort for E "
					"and E^tw) took time " << t2;
				found_co_factor = 1;
				found_co_factor_tw = 1;
				both_abort ++;
				stop_computation = true;
			}
			else {
				if (info) {
					std::cout << "\n\n-------------------------------------------------------\n";
					std::cout << "Working on prime l = " << ll;
				}

				if (ep.set_prime(ll)) {
					no_mod_equations = false;
					ep.compute_splitting_type();

					if (ep.is_elkies()) {
						ep.compute_trace_elkies();

						qq = remainder(q, ep.get_prime());
						g = gcd(q + 1 - static_cast<long>((ep.get_relation())[0]),
							ep.get_prime());
						if (! g.is_one())
							multiply(found_co_factor, found_co_factor, g);
						g = gcd(q + 1 + static_cast<long>((ep.get_relation())[0]),
							ep.get_prime());
						if (! g.is_one())
							multiply(found_co_factor_tw, found_co_factor_tw, g);
					}
					else
						ep.compute_trace_atkin();

					tm.set_vector(ep.get_prime(), ep.get_relation());
					enough_relations = tl.append(tm);
				}
				else if (ll >= 1020) {
					std::cout << "\nNot enough modular equations. Aborting." << std::endl;
					stop_computation = true;
				}

				ll = next_prime(ll+1);
			}
		} while (!enough_relations && !stop_computation);

		// Exit, if there are no modular equations
		if (no_mod_equations) {
			std::cout << "\nNot enough modular equations." << std::endl;
			std::cout << "\nExiting program." << std::endl;
			exit(1);
		}

		// Early abort
		//
		if (stop_computation)
			// found factor
			continue;

		// Determine group order
		//
		res = tl.bg_search_for_order();
		t2.stop_timer();
		std::cout << "\nTest No. " << tries-1 << " took time " << t2 << std::flush;

		// Test for good curve
		//
		if (!(co_factor.is_one() && !is_prime(res)))
			if (found_co_factor <= co_factor) {
				rf.assign(res);
				rf.verbose(info);
				rf.trialdiv();
				if (decimal_length(co_factor) > 6)
					rf.ecm(decimal_length(co_factor));
				else
					rf.ecm(15);
				std::cout << "\n#E = " << rf << std::flush;

				if (res/rf.base(rf.no_of_comp()-1) <= co_factor) {
					std::cout << "\nE is a good curve !!\n" << std::endl;
					good_curve = true;
				}
			}
			else {
				std::cout << "\nearly abort for E" << std::flush;
			}

		// Test for good twist
		//
		if (!good_curve) {
			if (found_co_factor_tw <= co_factor) {
				res = q + 1 + (q+1-res);
				rf.assign(res);
				rf.verbose(info);
				rf.trialdiv();
				if (decimal_length(co_factor) > 6)
					rf.ecm(decimal_length(co_factor));
				else
					rf.ecm(15);
				std::cout << "\n#E^tw = " << rf << std::flush;

				if (!(co_factor.is_one() && !is_prime(res)))
					if (res/rf.base(rf.no_of_comp()-1) <= co_factor) {
						std::cout << "\n#E^tw is a good curve !!\n" << std::flush;

						do {
							a2.randomize();
						} while (a2.trace() != 1);

						std::cout << "\nE^tw :  Y^2 + X Y = X^3 + " << a2 << " X^2 + ";
						std::cout << a6 << "\n";
						good_curve = true;
					}
			}
			else {
				std::cout << "\nearly abort for E^tw" << std::flush;
			}
		}

		no_abort ++;
	} while (!good_curve);
	t.stop_timer();

	std::cout << "\n\n========================================================";


	galois_field theField(2, gf2n::get_absolute_degree());
	gf_element One(theField);
	gf_element Zero(theField);
	One.assign_one();
	Zero.assign_zero();
	gf_element gfElementA2 = eco_gf2n::convertToGFElement(a2);
	gf_element gfElementA6 = eco_gf2n::convertToGFElement(a6);

	e.set_coefficients(One, gfElementA2, Zero, Zero, gfElementA6);
	point<gf_element> P(e), H(e);

	int i;
	for (i = 0; i < 5; i++) {
		P = e.random_point();
		multiply(H, res, P);

		if (!H.is_zero()) {
			std::cout << "\n\nProbabilistic correctness rejects candidate !!" << std::flush;
			std::cout << "\n\n"; exit(1);
		}
	}

	std::cout << "\n\nProbabilistic correctness proof agrees !!";
	std::cout << "\n\nField GF(2^" << n << ") with " << q << " elements (";
	std::cout << decimal_length(q) << ")" << std::flush;
	std::cout << "\nCoefficient a2 : "; std::cout << a2;
	std::cout << "\nCoefficient a6 : "; std::cout << a6;
	std::cout << std::endl << std::endl;
	std::cout << "Group order is " << res << std::endl;
	std::cout << "             = " << rf;
	if (rf.is_prime_factorization())
		std::cout << " (Prime Factorization)" << std::endl;
	else {
		std::cout << " (No Prime Factorization)" << std::endl;

		std::cout << "\nWARNING: This is no prime factorization, but no 'reasonable small'" << std::endl;
		std::cout << "composite factor could be found with ECM. Therefore the curve" << std::endl;
		std::cout << "is good with very high probability. You should however use" << std::endl;
		std::cout << "a factorization program to compute the prime factorization" << std::endl;
		std::cout << "and double check the 'goodness' of the curve." << std::endl;
	}

	// Test that order does not divide (p^B)-1 for B small
	// (Reduction to finite fields attack).
	//

	int log2 = q.bit_length() - 1;
	int stop = static_cast<int>(std::floor(2000.0 / log2));

	std::cout << "\n\nTest whether order divides q^k-1 for some 1 <= k <= ";
	std::cout << stop << " ... " << std::flush;

	multi_bigmod p_mod(q, res);
	multi_bigmod ppower_mod(1, res);

	for (i = 1; i <= stop; i++) {
		multiply(ppower_mod, ppower_mod, p_mod);
		if (ppower_mod.is_one()) {
			std::cout << "\n\nYES (k = " << i << ")";
			std::cout << "\nIf you want to avoid the MOV attack, you should";
			std::cout << "\nstart the program again and compute a new curve.\n";
			i = stop+2;
		}
	}
	if (i == stop+1)
		std::cout << "\n\nNO ==> MOV attack is not feasable for given curve.\n";


	std::cout << "\n\nComplete Test took " << tries-1 << " tries and time " << t << "\n";
	std::cout << "\nearly abort for both curves was used in " << both_abort << " tests ";
	std::cout << "\nfull group computation was used in " << no_abort << " tests \n\n";
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
