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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/udigit.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



lidia_size_t next_power_of_two(lidia_size_t m)
{
	debug_handler("tools.c", "next_power_of_two (lidia_size_t)");
	lidia_size_t k, n;
	n = 1;
	k = 0;
	while (n < m) {
		n = n << 1;
		k++;
	}
	return k;
}



lidia_size_t square_root(lidia_size_t n)
{
	debug_handler("tools.c", "square_root (lidia_size_t)");
	if (n <= 0)
		return 0;
	else if (n == 1)
		return 1;
	else
		return static_cast<lidia_size_t>(std::sqrt(static_cast<double>(n)));
}



//*************************************************************************
//				my_timer
//*************************************************************************

my_timer::my_timer() :
	msg(0)
{
	debug_handler("my_timer", "my_timer()");
	t.set_print_mode(HMS_MODE);
}



my_timer::~my_timer()
{
	debug_handler("my_timer", "~my_timer()");
	delete[] msg;
}



void my_timer::start(const char* message)
{
	debug_handler("my_timer", "start(char*)");
	delete[] msg;
	msg = new char[strlen(message)+1];
	memory_handler(msg, "my_timer",
		       "start(char*)::Error in memory allocation");
	strcpy(msg, message);
	t.start_timer();
}



void my_timer::stop()
{
	debug_handler("my_timer", "stop()");
	if (msg == 0)
		return;
	t.stop_timer();
	std::cerr << msg << "   ";
	t.print(std::cerr);
	std::cerr << std::endl;
	delete[] msg;
	msg = 0;
}



//*************************************************************************
//			    int_factor/fac_vec
//*************************************************************************


void swap(int_factor& i1, int_factor& i2)
{
	debug_handler("int_factor", "swap (int_factor&, int_factor&)");

	lidia_size_t t_lidia_size_t;
	int t_int;
	double t_double;

	t_lidia_size_t = i1.q;
	i1.q = i2.q;
	i2.q = t_lidia_size_t;

	t_int = i1.a;
	i1.a = i2.a;
	i2.a = t_int;

	t_int = i1.b;
	i1.b = i2.b;
	i2.b = t_int;

	t_double = i1.len;
	i1.len = i2.len;
	i2.len = t_double;
}



fac_vec::fac_vec(lidia_size_t n) :
	fvec(0),
	num_factors(0)
{
	debug_handler("fac_vec", "fac_vec(lidia_size_t)");

	compute_factorization(n);
}



void fac_vec::compute_factorization(lidia_size_t n)
{
	debug_handler("fac_vec", "compute_factorization (lidia_size_t)");
	lidia_size_t q;

	delete[] fvec;
	lidia_size_t m = next_power_of_two(n);
	fvec = new int_factor[m];
	memory_handler(fvec, "fac_vec", "compute_factorization(lidia_size_t)::"
		       "Error in memory allocation");

	num_factors = 0;
	q = 2;

	// compute factorization of n  
	while (n != 1) {
		if (n%q == 0) {
			fvec[num_factors].q = q;
			n = n/q;
			fvec[num_factors].a = 1;
			while (n%q == 0) {
				n = n/q;
				(fvec[num_factors].a)++;
			}

			num_factors++;
		}

		q++;
	}

	// set 'len'  
	lidia_size_t i;
	for (i = 0; i < num_factors; i++)
		fvec[i].len = fvec[i].a * std::log(static_cast<double>(fvec[i].q));

	sort();
}



void fac_vec::sort()
{
	debug_handler("fac_vec", "sort (void)");
	int i, j;

	for (i = 1; i <= num_factors - 1; i++)
		for (j = 0; j <= num_factors - i - 1; j++)
			if (fvec[j].len > fvec[j+1].len)
				swap(fvec[j], fvec[j+1]);
}



void fac_vec::clear(int lo, int hi)
{
	debug_handler("fac_vec", "assign_zero_to_b(int, int)");
	int i;
	for (i = lo; i <= hi; i++)
		fvec[i].b = 0;
}



lidia_size_t fac_vec::prod(int lo, int hi) const
{
	debug_handler("fac_vec", "prod(int, int)");

	lidia_size_t i, res;

	res = 1;
	for (i = lo; i <= hi; i++)
		res = res * power(static_cast<unsigned int>(fvec[i].q), static_cast<unsigned int>(fvec[i].a));

	return res;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
