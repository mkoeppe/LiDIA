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



#include	"LiDIA/mv_poly.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main(int, char**)
{
	int n, k, i, j;
	int anz_terms;
	int anz;

	std::cout << "\nTiming program for mv_poly : \n";
	std::cout <<   "============================\n\n";
	std::cout << "\nPlease input degree n of field GF(2^n) : ";
	std::cin >> n;
	gf2n_init(n);

	bigint bound;
	gf2n r;

	std::cout << "\nPlease input maximal number of variables : "; std::cin >> k;

	shift_left(bound, bigint(1), k);
	dec(bound);

	std::cout << "\nPlease input maximal number of terms : "; std::cin >> anz_terms;
	std::cout << "\nPlease enter number of tests : "; std::cin >> anz;
	std::cout << "\n\nChoose random polynomials ... " << std::flush;

	mv_poly p1[anz], p2[anz], p3, p4;
	bigint xx;

	for (j = 0; j < anz; j++)
		for (i = 0; i < anz_terms; i++) {
			r.randomize();
			xx.randomize(bound);
			p1[j] += mv_term(r, xx);

			r.randomize();
			xx.randomize(bound);
			p2[j] += mv_term(r, xx);
		}

	std::cout << "\nStarting multiplication tests ... " << std::flush;
	timer t; t.set_print_mode();

	t.start_timer();
	for (i = 0; i < anz; i++)
		p3 = p1[i] * p2[i];
	t.stop_timer();

	std::cout << "\n" << anz << " mults need time " << t << std::endl;


	std::cout << "\n==> DONE.\n\n";
}
