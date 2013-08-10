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


#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/base/interface_lib.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define STARTVALUE 4
//STARTVALUE must be > 3

const int max_prime_size = 640; // for primes bigger 2^640 the last computed
const int upper_bound_gcd = 160; // value is taken (smaller for gcd)

int WDH = 50;
Fp_polynomial test_array1[50];
Fp_polynomial test_array2[50];


void init_test_arrays(const bigint &p, lidia_size_t n1, lidia_size_t n2)
{
	lidia_size_t i;
	for (i = 0; i < WDH; i++) {
		randomize(test_array1[i], p, n1);
		randomize(test_array2[i], p, n2);
	}
}



void
test_mul(long &t_p, long &t_f, int n, const bigint &p)
{
	init_test_arrays(p, n, n);
	Fp_polynomial c;
	int k;
	timer t;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_mul(c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		fft_mul(c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_f = t.user_time();
}



void
test_rem(long &t_p, long &t_f, int n, const bigint &p)
{
	init_test_arrays(p, 2*n, n);
	timer t;
	int k;
	Fp_polynomial d;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_rem(d, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		fft_rem(d, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_f = t.user_time();
}



void
test_div(long &t_p, long &t_f, int n, const bigint &p)
{
	init_test_arrays(p, 2*n, n);
	timer t;
	int k;
	Fp_polynomial c, d;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_div_rem(d, c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		fft_div_rem(d, c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_f = t.user_time();
}



void
test_inv(long &t_p, long &t_f, int n, const bigint &p)
{
	timer t;
	int k;

	for (k = 0; k < WDH; k++)
		do
			randomize(test_array1[k], p, n);
		while (test_array1[k][0] == 0);
	Fp_polynomial c;
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



void
test_gcd(long &t_p, long &t_f, int n, const bigint &p)
{
	init_test_arrays(p, n, n);
	Fp_polynomial c;
	int k;
	timer t;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_gcd(c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		gcd(c, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_f = t.user_time();
}



void
test_xgcd(long &t_p, long &t_f, int n, const bigint &p)
{
	init_test_arrays(p, n, n);
	Fp_polynomial c, s, tt;
	int k;
	timer t;
	t.start_timer();
	for (k = 0; k < WDH; k++)
		plain_xgcd(c, s, tt, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_p = t.user_time();

	t.start_timer();
	for (k = 0; k < WDH; k++)
		xgcd(c, s, tt, test_array1[k], test_array2[k]);
	t.stop_timer();
	t_f = t.user_time();
}



bigint seek_prime_with_n_bits(int n)
{
	if (n < 2) n = 2;
	bigint t(1);
	t = (t << n) - 1; //the largest n-bit number

	if (n > max_prime_size)
		return bigint(-1);
	else
		return previous_prime(t);
}



int
crossover(void (*FUNCTION)(long&, long&, int, const bigint&), const bigint &p,
          const char * NAME)
	//compute crossover-point for a given FUNCTION and a fixed modulus p
{
	long l1, l2;
	long plain_time[20], new_time[20];
	int counter = -1;
	std::cout << NAME << " -- prime : " << p << " (" << p.bit_length() << " bits)" << std::endl;

	int n;

	std::cout << "\tdegree\t\tplain\tnew" << std::endl;
	std::cout << "\t-----------------------------------" << std::endl;
	n = 1 << (STARTVALUE-1);

	do {
		counter++;
		n = 2*n;
		FUNCTION(l1, l2, n-1, p);
		std::cout << "\t" << n-1 << "\t\t" << l1 << "\t" << l2 << std::endl;
		plain_time[counter] = l1;
		new_time[counter] = l2;
	} while ((static_cast<double>(l1)/static_cast<double>(l2) < 1.01) || (l1 + l2 < 100));

	bool measured = true;
	int interval = n/2, h;
	while (interval >= 2) {
		if (l1 + l2 > 100 && !measured)  // if sufficiently large: stop if
			// very close
		{
			h = l1 - l2;
			if (h < 0)
				h = - h;
			if (h <= (0.005 * comparator< lidia_size_t >::max(l1, l2)))
				break;
		}

		if (static_cast<double>(l1) / static_cast<double>(l2) < 1.01) {
			n += interval;
			measured = false;
		}
		else
			n -= interval;
		interval = interval/2;
		if (measured) {
			counter--;
			l1 = plain_time[counter];
			l2 = new_time[counter];
			std::cout << "\t" << n-1 << "\t\t*" << l1 << "\t*" << l2 << std::endl;
		}
		else {
			FUNCTION(l1, l2, n, p);
			std::cout << "\t" << n << "\t\t" << l1 << "\t" << l2 << std::endl;
		}
	}
	std::cout << "Set crossover to " << n << std::endl;
	return n;
}



void do_work(int crov[CROV_NUM_VALUES],
	     void (*FUNCTION)(long&, long&, int, const bigint&),
	     const char* NAME, const bigint prime_table[CROV_NUM_VALUES])
{
	int i;
	std::cout << NAME << " : " << std::endl;
	for (i = 0; i < CROV_NUM_VALUES; i++)
		if (prime_table[i] == -1)
			crov[i] = crov[i-1];
		else
			crov[i] = crossover(FUNCTION, prime_table[i], NAME);
}



int main()
{
	std::cout << "This program tries to find some crossover points." << std::endl;
	std::cout << "This may take a while (up to several hours)." << std::endl << std::endl;

	int x[CROV_NUM_VALUES] = {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120};
	//bitlengths of moduli
	int y_mul[CROV_NUM_VALUES];
	int y_div[CROV_NUM_VALUES];
	int y_rem[CROV_NUM_VALUES];
	int y_inv[CROV_NUM_VALUES];
	int y_gcd[CROV_NUM_VALUES];
	int y_xgcd[CROV_NUM_VALUES];
	int halfgcd = 16;
	int log2_newton = 5;

	int i;
	bigint prime_table[CROV_NUM_VALUES];
	std::cout << "Generating primes" << std::endl;
	for (i = 0; i < CROV_NUM_VALUES; i++) {
		prime_table[i] = seek_prime_with_n_bits(x[i]);
		std::cout << "+" << std::flush;
	}
	std::cout << " done." << std::endl << std::endl;

	do_work(y_mul, test_mul, "multiply", prime_table);
	for (i = 0; i < CROV_NUM_VALUES; i++)
		std::cout << "y_mul[" << i << "] = " << y_mul[i] << std::endl;

	do_work(y_rem, test_rem, "remainder", prime_table);
	for (i = 0; i < CROV_NUM_VALUES; i++)
		std::cout << "y_rem[" << i << "] = " << y_rem[i] << std::endl;

	do_work(y_inv, test_inv, "invert", prime_table);
	for (i = 0; i < CROV_NUM_VALUES; i++)
		std::cout << "y_inv[" << i << "] = " << y_inv[i] << std::endl;

	do_work(y_div, test_div, "divide", prime_table);
	for (i = 0; i < CROV_NUM_VALUES; i++)
		std::cout << "y_div[" << i << "] = " << y_div[i] << std::endl;


	for (i = 0; i < CROV_NUM_VALUES; i++)
		if (prime_table[i].bit_length() > upper_bound_gcd)
			prime_table[i].assign(-1);

	WDH = WDH / 4;

	do_work(y_gcd, test_gcd, "gcd", prime_table);
	for (i = 0; i < CROV_NUM_VALUES; i++)
		std::cout << "y_gcd[" << i << "] = " << y_gcd[i] << std::endl;

	do_work(y_xgcd, test_xgcd, "xgcd", prime_table);
	for (i = 0; i < CROV_NUM_VALUES; i++)
		std::cout << "y_xgcd[" << i << "] = " << y_xgcd[i] << std::endl;

	std::cout << "Writing file...";

	std::ofstream s;
	s.open("Fp_pol_crossover.h");
	if (!s) {
		lidia_error_handler("Fp_polynomial", "crossover::error while writing"
				    "the file \"crossover.tbl\"");
		return -1;
	}

	s << "// LiDIA - a library for computational number theory\n";
	s << "//   Copyright (c) 1994 - 1999 by the LiDIA Group\n";
	s << "//\n";
	s << "// This file was generated automatically by 'make optimize'.\n";
	s << "// For default values, copy 'crossover.tbl.default' to\n";
	s << "// 'crossover.tbl', type  'touch crossover.tbl' and\n";
	s << "// 'make'.\n";
	s << "//\n";
	s << "\n";

	s << "static int x_val[CROV_NUM_VALUES] = { " << x[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << x[i];
	s << " }; \n";

	s << "static int fftmul_val[CROV_NUM_VALUES] = { " << y_mul[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << y_mul[i];
	s << " }; \n";

	s << "static int fftdiv_val[CROV_NUM_VALUES] = { " << y_div[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << y_div[i];
	s << " }; \n";

	s << "static int fftrem_val[CROV_NUM_VALUES] = { " << y_div[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << y_rem[i];
	s << " }; \n";

	s << "static int inv_val[CROV_NUM_VALUES] = { " << y_inv[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << y_inv[i];
	s << " }; \n";

	s << "static int gcd_val[CROV_NUM_VALUES] = { " << y_inv[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << y_gcd[i];
	s << " }; \n";

	s << "static int xgcd_val[CROV_NUM_VALUES] = { " << y_inv[0];
	for (i = 1; i < CROV_NUM_VALUES; i++)
		s << ", " << y_xgcd[i];
	s << " }; \n";

	s << "\n";
	s << "int halfgcd_val = " << halfgcd << "; \n";
	s << "int log2_newton_val = " << log2_newton << "; \n";
	s << std::endl;
	s.close();

	std::cout << " done." << std::endl;
	return 0;
}
