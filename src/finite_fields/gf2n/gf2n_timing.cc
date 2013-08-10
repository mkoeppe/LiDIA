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


#include	"LiDIA/gf2n.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



const unsigned int tests = 10000;
const unsigned int jmax = 10;



int main(int, char **)
{
	unsigned int  deg;
	register unsigned int i, j;
	timer t;

	t.set_print_mode(HMS_MODE);

	std::cout << "\nInput Extension Degree: "; std::cin >> deg;

	gf2n_init(deg);
	gf2n a;
	gf2n *tv;

	tv = new gf2n[tests+1];

	for (i = 0; i <= tests; i++)
		tv[i].randomize();

	t.start_timer();
	for (j = 0; j < jmax; j++) {
		for (i = 0; i < tests; i++)
			add(a, tv[i], tv[i+1]);
	}
	t.stop_timer();
	std::cout << "\n " << jmax*tests << " additions take time         " << t << std::flush;

	t.start_timer();
	for (j = 0; j < jmax; j++) {
		for (i = 0; i < tests; i++)
			square(a, tv[i]);
	}
	t.stop_timer();
	std::cout << "\n " << jmax*tests << " squarings take time          " << t << std::flush;

	t.start_timer();
	for (j = 0; j < jmax; j++) {
		for (i = 0; i < tests; i++)
			multiply(a, tv[i], tv[i+1]);
	}
	t.stop_timer();
	std::cout << "\n " << jmax*tests << " multiplications take time   " << t << std::flush;

	t.start_timer();
	for (j = 0; j < jmax; j++) {
		for (i = 0; i < tests; i++)
			invert(a, tv[i]);
	}
	t.stop_timer();
	std::cout << "\n " << jmax*tests << " inversions take time        " << t << std::flush;
	std::cout << "\n\n";

        delete[] tv;
	return 0;
}
