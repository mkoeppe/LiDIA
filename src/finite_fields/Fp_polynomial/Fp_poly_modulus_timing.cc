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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main(int, char **)
{
	lidia_size_t n, d, i, tests;
	bigint p;
	timer t;
	t.set_print_mode(HMS_MODE);

	std::cout << "\n Timing Program for Fp_polynomial and Fp_poly_modulus";
	std::cout << "\n ====================================================\n\n";

	std::cout << "\n Input prime : ";
	std::cin >> p;

	p = next_prime(p);
	std::cout << "\n Choosing next prime : p = " << p << " (" << decimal_length(p) << ")";

	std::cout << "\n\n Input polynomial degree : ";
	std::cin >> d;

	std::cout << "\n Input number of tests : ";
	std::cin >> tests;

	Fp_polynomial f;
	Fp_poly_modulus fpm;
	Fp_polynomial erg1, erg2, x, xq;

	std::cout << "\n\n Arithmetical operations first :";
	std::cout << "\n -------------------------------\n\n" << std::flush;

	randomize(f, p, d);
	f.make_monic();
	fpm.build(f);

	randomize(erg1, p, d-1);
	randomize(erg2, p, d-1);

	t.start_timer();
	for (i = 1; i <= tests; i++)
		multiply(x, erg1, erg2);
	t.stop_timer();
	std::cout << tests << " Multiplications without reduction take time " << t << "\n" << std::flush;

	t.start_timer();
	for (i = 1; i <= tests; i++)
		square(x, erg1);
	t.stop_timer();
	std::cout << tests << " Squarings without reduction take time " << t << "\n" << std::flush;

	t.start_timer();
	for (i = 1; i <= tests; i++)
		multiply(x, erg1, erg2, fpm);
	t.stop_timer();
	std::cout << tests << " Multiplications with reduction take time " << t << "\n" << std::flush;

	t.start_timer();
	for (i = 1; i <= tests; i++)
		square(x, erg1, fpm);
	t.stop_timer();
	std::cout << tests << " Squarings with reduction take time " << t << "\n" << std::flush;

	std::cout << "\n\n High level operations for factoring polynomials : ";
	std::cout << "\n -------------------------------------------------" << std::flush;

#if 0
	for (i = 1; i <= tests; i++) {
		randomize(f, p, d);
		f.make_monic();
		fpm.build(f);

		x.assign_x();
		t.start_timer();
		power_x(xq, p, fpm);
		t.stop_timer();
		std::cout << "\n\nTest " << i << ":  Power_x() needs time " << t << std::flush;

		erg1.assign_x();
		t.start_timer();
		Xq(erg1, fpm);
		t.stop_timer();
		std::cout << "\nTest " << i << ":  Xq needs time " << t << std::flush;

		t.start_timer();
		compute_degree(xq, fpm, d/2);
		t.stop_timer();
		std::cout << "\nTest " << i << ":  compute_degree(" << d/2 << ") needs time " << t << std::flush;

		t.start_timer();
		compute_degree(xq, fpm, d);
		t.stop_timer();
		std::cout << "\nTest " << i << ":  compute_degree(" << d << ") needs time " << t << std::flush;
	}
#endif
	std::cout << "\n\n";

	return 0;
}
