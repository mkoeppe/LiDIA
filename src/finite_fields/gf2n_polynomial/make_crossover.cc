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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf2n_poly_modulus.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define WDH 70
#define STARTVALUE 5
//STARTVALUE must be > 3

gf2n_polynomial test_array1[WDH];
gf2n_polynomial test_array2[WDH];


void init_test_arrays(lidia_size_t n1, lidia_size_t n2)
{
	lidia_size_t i;
	for (i = 0; i < WDH; i++) {
		test_array1[i].randomize(n1);
		test_array2[i].randomize(n2);
	}
}

//------------------------------------------------------------------

void
test_inv(long &t_p, long &t_f, int n)
{
	timer t;
	int k;

	for (k = 0; k < WDH; k++)
		do {
			test_array1[k].randomize(n - 1);
		} while (test_array1[k].const_term().is_zero());

	gf2n_polynomial c;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_inv(c, test_array1[k], n);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		newton_inv(c, test_array1[k], n);
	t.stop_timer();
	t_f = t.user_time();
}

//------------------------------------------------------------------
void
test_red(long &t_p, long &t_f, int n)
{
	timer t;
	int k;

	init_test_arrays(2*n-2, n);

	gf2n_polynomial c;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_rem(c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		kara_rem(c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_f = t.user_time();
}


void
test_square_red(long &t_p, long &t_f, int n)
{
	timer t;
	int k;
	gf2n_poly_modulus fpm;
	gf2n_polynomial f;

	f.randomize(n);
	f.make_monic();
	fpm.build(f);

	for (k = 0; k < WDH; k++) {
		test_array1[k].randomize(n-1);
		square(test_array1[k], test_array1[k]);
	}

	gf2n_polynomial c;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_rem(c, test_array1[k], fpm);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		kara_rem(c, test_array1[k], fpm);
	t.stop_timer();
	t_f = t.user_time();
}



int
crossover(void (*FUNCTION)(long&, long&, int))
{
	long l1, l2;
	long plain_time[20], new_time[20];
	int counter = -1;
	int first = 1 << (STARTVALUE);
	int n;

	std::cout << "\tdegree\t\tplain \tnew" << std::endl;
	std::cout << "\t-----------------------------------" << std::endl;
	n = 1 << (STARTVALUE-1);

	do {
		counter++;
		n = 2*n;
		FUNCTION(l1, l2, n);
		std::cout << "\t" << n << "\t\t" << l1 << "\t" << l2 << std::endl;
		plain_time[counter] = l1;
		new_time[counter] = l2;
	} while ((static_cast<double>(l1)/static_cast<double>(l2) < 1.1) || (l1 + l2 < 100));

	bool measured = true;
	int interval = n/2;
	while (interval >= 2) {
		if (l1 < l2) {
			n += interval;
			measured = false;
		}
		else
			n -= interval;
		interval = interval/2;

		if (n < first)
			measured = false;

		if (measured) {
			counter--;
			l1 = plain_time[counter];
			l2 = new_time[counter];
			std::cout << "\t" << n << "\t\t*" << l1 << "\t*" << l2 << std::endl;
		}
		else {
			FUNCTION(l1, l2, n);
			std::cout << "\t" << n << "\t\t" << l1 << "\t" << l2 << std::endl;
		}
	}
	std::cout << "Set crossover to " << n << std::endl;
	return n;
}


int main()
{
	std::cout << "This program tries to find some crossover points." << std::endl;
	std::cout << "This may take a while." << std::endl << std::endl;

	int n = 521; // this will be flexibel later
	int y_inv;
	int y_red;
	int y_square_red;

	gf2n_init(n);

	std::cout << "\ninvert : \n\n" << std::flush;
	y_inv = crossover(test_inv);

	std::cout << "\nreduce : \n\n" << std::flush;
	y_red = crossover(test_red);

	std::cout << "\nreduce after squaring : \n\n" << std::flush;
	y_square_red = crossover(test_square_red);

	std::cout << "Writing file...";
	std::ofstream s;
	s.open("crossover.tbl");
	if (!s) {
		lidia_error_handler("gf2n_polynomial", "crossover::error while "
				    "writing the file \"crossover.tbl\"");
		return 0;
	}

	s << "// LiDIA - a library for computational number theory\n";
	s << "//   Copyright (c) 1994 - 2000 by the LiDIA Group\n";
	s << "//\n";
	s << "// This file was generated automatically by 'make crossover'.\n";
	s << "// For default values, copy 'crossover.tbl.default' to\n";
	s << "// 'crossovers.tbl', type  'touch crossover.tbl' and\n";
	s << "// 'make'.\n";
	s << "//\n";
	s << "\n";

	s << "const int CROSSOVER_PLAININV_NEWTON = " << y_inv << "; ";
	s << "const int CROSSOVER_REDUCE = " << y_red << "; ";
	s << "const int CROSSOVER_REDUCE_SQUARE = " << y_square_red << "; ";
	s << "\n";
	s.close();

	std::cout << " done." << std::endl;
	return 0;
}
