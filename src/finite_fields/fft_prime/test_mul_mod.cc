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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/timer.h"
#include	"LiDIA/udigit.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



// long mul_mod_sp(long a, long b, long q, double qinv)
//    returns a*b mod q
//    qinv must be  1/double(q)
//
// long mul_mod_sp2(long a, long b, long q, double bqinv)
//    returns a*b mod q
//    bqinv must be  double(b)/double(q)
//
// Since their implementation is based on double-arithmetic,
// long-arguments must not be greater than 2^26 !!


inline
long mul_mod_sp(long a, long b, long p, double pinv)
{
	double ab = static_cast<double>(a) * static_cast<double>(b);
	register long q = static_cast<long>(ab * pinv);
	register long res = static_cast<long>(ab - (static_cast<double>(q) * static_cast<double>(p)));

	if (res >= p)
		res -= p;
	else if (res < 0)
		res += p;
	return
		res;
}



inline
long mul_mod_sp2(long a, long b, long p, double bpinv)
{
	double ab = static_cast<double>(a) * static_cast<double>(b);
	register long q = static_cast<long>(a* bpinv);
	register long res = static_cast<long>(ab - (static_cast<double>(q) * static_cast<double>(p)));

	if (res >= p)
		res -= p;
	else if (res < 0)
		res += p;
	return
		res;
}



const int loop = 1000;
const int tests = 10000;



int main ()
{
	int j, i;
	//udigit *a, *b;
	udigit p, c;
	timer t;
	double speed1, speed2;
	random_generator rg;

	std::cout << "Testing the relative speed of udigits in comparison to "
		"doubles ..." << std::flush;

	t.set_print_mode();
	p = previous_prime(max_udigit() / 2);

        std::vector<udigit> a(tests);
        std::vector<udigit> b(tests);

	for (i = 0; i < tests; i++) {
		rg >> a[i];
		a[i] %= p;
		rg >> b[i];
		b[i] %= p;
	}

	t.start_timer();
	for (j = 0; j < loop; j++)
		for (i = 0; i < tests; i++)
			c = multiply_mod(a[i], b[i], p);
	t.stop_timer();
	speed1 = t.user_time();
	std::cout << "\nmul_mod : " << speed1 << std::flush;


	p = previous_prime(1 << 26 - 1);
	for (i = 0; i < tests; i++) {
		a[i] = a[i] % p;
		b[i] = b[i] % p;
	}

	double qinv = 1.0 / static_cast<double>(p);

	t.start_timer();
	for (j = 0; j < loop; j++)
		for (i = 0; i < tests; i++)
			c = mul_mod_sp(a[i], b[i], p, qinv);
	t.stop_timer();

	speed2 = t.user_time() * (1.0 * (bits_per_udigit()-1)) / 26.0;
	std::cout << "\nmul_mod_double : " << speed2 << std::flush;

	std::ofstream s;
	s.open("fft_mul_mod.inl");
	if (!s) {
		lidia_error_handler("fft_prime", "error while writing"
				    "the file \"fft_mul_mod.inl\"");
		return -1;
	}

	if (speed1 < 1.01 * speed2) {
		std::cout << "\n ==  > Choosing udigit's as basis\n";
		s << "#undef LIDIA_MUL_MOD_SPECIAL\n";
	}
	else {
		std::cout << "\n ==  > Choosing doubles as basis\n";
		s << "#define LIDIA_MUL_MOD_SPECIAL 1 \n";
	}
	s.close();
	return 0;
}
