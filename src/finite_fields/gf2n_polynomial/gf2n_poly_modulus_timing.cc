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

	std::cout << "\n Input polynomial degree : ";
	std::cin >> d;

	gf2n_polynomial f;
	gf2n_poly_modulus fpm;

	std::cout << "\n\n Arithmetical operations first :";
	std::cout << "\n -------------------------------\n" << std::flush;
	std::cout << "\n Input number of arithmetical tests : ";
	std::cin >> tests;

	gf2n_polynomial erg1[tests], erg2[tests], erg3, x, xq;

	f.randomize(d);
	f.make_monic();
	fpm.build(f);

	gf2n_polynomial * table;

	for (i = 0; i < tests; i++) {
		erg1[i].randomize(d-1);
		erg2[i].randomize(d-1);
	}

	table = new gf2n_polynomial[d];
	table[0].assign(gf2n_polynomial(d) % f);

	for (i = 1; i < d; i++) {
		shift_left(table[i], table[i-1], 1);
		remainder(table[i], table[i], f);
	}

	std::cout << "\n Initialisation done, starting tests ... " << std::flush;

	t.start_timer();
	for (i = 0; i < tests; i++)
		multiply(x, erg1[i], erg2[i], fpm);
	t.stop_timer();
	std::cout << "\n " << tests << " Multiplications with reduction take time " << t << std::flush;

	t.start_timer();
	for (i = 0; i < tests; i++)
		square(x, erg1[i], fpm);
	t.stop_timer();
	std::cout << "\n " << tests << " Squarings with reduction take time " << t << std::flush;

	for (i = 0; i < tests; i++)
		multiply(erg1[i], erg1[i], erg2[i]);

	t.start_timer();
	for (i = 0; i < tests; i++)
		remainder(x, erg1[i], f);
	t.stop_timer();
	std::cout << "\n " << tests << " Reduction take time " << t << std::flush;

	int j, dd;
	gf2n e;

	t.start_timer();
	for (i = 0; i < tests; i++) {
		dd = erg1[i].degree();
		x.assign_zero();
		for (j = d; j < dd; j++) {
			erg1[i].get_coefficient(e, j);
			if (!e.is_zero())
				add(x, x, e * table[j-d]);
		}
	}
	t.stop_timer();
	std::cout << "\n " << tests << " Reduction with table take time " << t << std::flush;

	delete[] table;

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

		t.start_timer();
		Xq(erg1[0], fpm);
		t.stop_timer();
		std::cout << "\nTest " << i << ":  Xq needs time " << t << std::flush;
		xq.assign(erg1[0]);

		table = new gf2n_polynomial[d+1];
		t.start_timer();
		table = power_table(erg1[0], fpm, fpm.modulus().degree());
		compose(erg1[0], erg1[0], xq, table);
		t.stop_timer();
		delete[] table;
		std::cout << "\nTest " << i << ":  compose_table needs time " << t << std::flush;

		t.start_timer();
		compose(erg1[0], erg1[0], xq, fpm);
		t.stop_timer();
		std::cout << "\nTest " << i << ":  compose_polymod needs time " << t << std::flush;
	}

	std::cout << "\n\n";
	return 0;
}
