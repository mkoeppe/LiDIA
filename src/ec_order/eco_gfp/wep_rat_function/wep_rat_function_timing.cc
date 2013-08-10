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
//	Author	:Markus Maurer (MM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/wep_rat_function.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main(int argc, char* argv[])
{
	bigint p;
	timer t; t.set_print_mode(HMS_MODE);
	int tests;

	std::cout << "\n\nTiming Program for wep_rat_functions";
	std::cout << "\n====================================\n";
	std::cout << "\nCharacteristic : "; std::cin >> p;

	p = next_prime(p-1);
	std::cout << "\nChoose next prime, p = " << p << " (" << decimal_length(p) << ")" << std::flush;
	bigmod::set_modulus(p);
	bigmod a4, a6;
	int deg_modulus;

	std::cout << "\n\nCoefficient a4 : "; std::cin >> a4;
	std::cout << "\nCoefficient a6 : "; std::cin >> a6;
	std::cout << "\nDegree of Modulus : "; std::cin >> deg_modulus;
	std::cout << "\nNumber of Tests : "; std::cin >> tests;

	Fp_polynomial m;
	m.set_modulus(p);
	m.randomize(deg_modulus);
	m.make_monic();
	Fp_poly_modulus pm(m);

	wep_rat_function::initialize(a4, a6, pm);
	wep_rat_function P, mP, m2P;

	P.assign_xy();

	random_generator rg;
	long i;

	rg >> i; i = i % 100; multiply(mP, i, P);
	rg >> i; i = i % 100; multiply(m2P, i, P);

	t.start_timer();
	for (i = 0; i < tests; i++)
		multiply_by_2(mP, mP);
	t.stop_timer();

	std::cout << "\n" << tests << " doublings take time " << t << std::flush;

	t.start_timer();
	for (i = 0; i < tests; i++)
		add(mP, mP, m2P);
	t.stop_timer();

	std::cout << "\n" << tests << " additions take time " << t << "\n\n" << std::flush;
}
