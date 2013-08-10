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


#include	"LiDIA/gf2n_poly_modulus.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main()
{
	lidia_size_t n, d, i, tests;
	bigint q;
	timer t;
	t.set_print_mode(HMS_MODE);

	std::cout << "\n Timing Program for gf2n_polynomial and gf2n_poly_modulus ";
	std::cout << "\n========================================================\n\n";

	std::cout << "\n Input field extension degree over GF(2):  ";
	std::cin >> n;

	gf2n_init(n);
	shift_left(q, bigint(1), n);

	gf2nIO::noprefix();
	int low, high;

	std::cout << "\nInput lower bound deg(f) = "; std::cin >> low;
	std::cout << "\nInput upper bound deg(f) = "; std::cin >> high;


	gf2n_polynomial f;
	gf2n_poly_modulus fpm;
	gf2n_polynomial erg1[100], erg2[100], erg3, x, xq, erg4;

	std::cout << "\n\n Arithmetical operations first :";
	std::cout << "\n -------------------------------\n" << std::flush;
	//  std::cout<<"\n Input number of arithmetical tests : ";
	//std::cin >> tests;

	tests = 100;


	for (d = low; d <= high; d++) {
		for (i = 0; i < 100; i++) {
			erg1[i].randomize(d-1);
			erg2[i].randomize(d-1);
		}

		t.start_timer();
		for (i = 1; i <= tests; i++)
			multiply(erg3, erg1[i-1], erg2[i-1]);
		t.stop_timer();

		std::cout << "\nMult (" << d << ") = " << t;
		continue;

		t.start_timer();
		for (i = 1; i <= tests; i++)
			remainder(erg3, erg1[i-1], f);
		t.stop_timer();
		std::cout << "\nremainder = " << t;
	}

#if 0
	t.start_timer();
	for (i = 1; i <= tests; i++)
		plain_rem(erg3, erg1[i-1], f);
	t.stop_timer();
	std::cout << "\n " << tests << " plain reduction take time " << t << std::flush;

	t.start_timer();
	for (i = 1; i <= tests; i++)
		kara_rem(erg3, erg1[i-1], f);
	t.stop_timer();
	std::cout << "\n " << tests << " kara reduction take time " << t << std::flush;

#endif

#if 0
	std::cout << "\n\n High level operations for factoring polynomials : ";
	std::cout << "\n -------------------------------------------------" << std::flush;
	std::cout << "\n Input number of high-level tests : ";
	std::cin >> tests;

	for (i = 1; i <= tests; i++) {
		f.randomize(d);
		f.make_monic();
		fpm.build(f);

		x.assign_x();
		t.start_timer();
		power(xq, x, q, fpm);
		t.stop_timer();
		std::cout << "\n\nTest " << i << ":  Power(q) needs time " << t << std::flush;

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






