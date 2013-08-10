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
//	Author	: Andreas Mueller (AM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"
#include	<cmath>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const unsigned int ecm_primes::BIT = 0x80000000;



void ecm_primes::
initprimes(unsigned int us, unsigned int os, unsigned int sp)
{
	unsigned int h, j, k, step, i;

	max = os;
	space = sp;

	p_max = static_cast<unsigned int>(std::sqrt(static_cast<double>(max))) >> 1;

	h = static_cast<unsigned int>(500 + (p_max >> 3) + (p_max << 2));

	if (space< ((max >> 4)+500) && space < h) {
		lidia_error_handler("ecm_primes", "initprimes::not enough memory for primetables");
		return;
	}

	if (space< ((max >> 4)+500)) {
		// sieving on several intervals
		dim_ar = (space - h) >> 2;
		flag = 1;
	}
	else {
		// sieving on one interval
		dim_ar = max >> 6;
		flag = 0;
	}

	last_prime = dim_ar << 6;

	table = new unsigned int[dim_ar];

	if (flag) {
		h = (p_max >> 5)+1;
		over = new unsigned int [h << 5];
		table2 = new unsigned int [h];
	}

	maske35 = 0x5b4cb699; // patterns for the primes up to 31
	maske7 = 0xfdfbf7ef;
	maske11 = 0xfdffbff7;
	maske13 = 0x7ffbffdf;
	maske17 = 0xffdfffef;
	maske19 = 0xfffeffff;
	maske23 = 0xdfffffbf;
	maske29 = 0xffefffff;
	maske31 = 0xfffdffff;

	table[0] = 0x76d32d26;

	sieve_array(1);

	for (i = 17; i < p_max; i++) {
		// sieving primes between 31 and sqrt(max)
		j = 2 * i + 3;
		if (table[j >> 6] & (BIT >> ((j & 0x0000003f) >> 1))) {
			step = j << 1;
			for (k = j*j; k < last_prime; k += step) {
				pl = (k & 0x0000003f) >> 1;
				feld = k >> 6;
				table[feld] = table[feld] & ~(BIT >> pl);
			}
			if (flag)
				over[i] = k - last_prime; // overhead to the next interval
		}
	}

	if (flag) {
		for (i = 0; i < h; i++)
			table2[i] = table[i]; // saving primes up to sqrt(max)

		if (us > last_prime) {
			maske35 = 0x96d32da6;
			maske7 = 0xefdfbf7e;
			maske11 = 0xfbff7fef;
			maske13 = 0xfdffefff;
			maske17 = 0xff7fffbf;
			maske19 = 0xffbffff7;
			maske23 = 0xffefffff;
			maske29 = 0xfffdffff;
			maske31 = 0xfffeffff;

			k = us >> 6;
			j = k % 15;

			for (i = 0; i < j; i++) {
				maske35 = maske35 << 2;
				maske35 = ((maske35 & 0xc0000000) >> 30) | maske35;
			}

			j = k % 7;
			for (i = 0; i < j; i++) {
				maske7 = maske7 << 4;
				maske7 = ((maske7 & 0xf0000000) >> 28) | maske7;
			}

			j = k % 11;
			for (i = 0; i < j; i++) {
				maske11 = maske11 << 10;
				maske11 = ((maske11 & 0xffc00000) >> 22) | maske11;
			}

			j = k % 13;
			for (i = 0; i < j; i++) {
				maske13 = maske13 << 6;
				maske13 = ((maske13 & 0xfc000000) >> 26) | maske13;
			}

			j = k % 17;
			for (i = 0; i < j; i++) {
				maske17 = maske17 << 15;
				maske17 = ((maske17 & 0xfffe0000) >> 17) | maske17;
			}

			j = k % 19;
			for (i = 0; i < j; i++) {
				maske19 = maske19 << 13;
				maske19 = ((maske19 & 0xfff80000) >> 19) | maske19;
			}

			j = k % 23;
			for (i = 0; i < j; i++) {
				maske23 = maske23 << 9;
				maske23 = ((maske23 & 0xff800000) >> 23) | maske23;
			}

			j = k % 29;
			for (i = 0; i < j; i++) {
				maske29 = maske29 << 3;
				maske29 = ((maske29 & 0xe0000000) >> 29) | maske29;
			}

			j = k % 31;
			for (i = 0; i < j; i++) {
				maske31 = maske31 << 1;
				maske31 = ((maske31 & 0x80000000) >> 31) | maske31;
			}

			k = us & ~0x0000003f;

			for (i = 17; i < p_max; i++) {
				j = 2 * i + 3;
				if (table2[j >> 6] & (BIT >> ((j & 0x0000003f) >> 1))) {
					h = ((k / j) + 1) * j - k;
					if (!(h&1))   h += j;
					over[i] = h;
				}
			}
		}
	}

	feld = us >> 6; // pointer must be set
	feld_abs = feld;
	pl = (us & 0x03f) >> 1;
	maske_i = BIT >> pl;
	return;
}



unsigned int ecm_primes::
getprimes()
{
 st:
	while ((feld < dim_ar) && ((feld_abs << 6) < max)) {
		while (pl < 32)         // loop: find 1-bit in unsigned int
		{
			pl++;
			if (table[feld] & maske_i) {
				maske_i >>= 1;
				return((pl << 1)+(feld_abs << 6)-1);
			}
			maske_i >>= 1;
		}
		pl = 0;
		feld_abs++;
		feld++;
		maske_i = 0x80000000;
	}

	if (flag && (feld_abs< (max >> 6)))      // sieve next interval
	{
		sieve_array(0);

		unsigned int j, k = 0, step, i;

		for (i = 17; i < p_max; i++) {
			j = 2 * i + 3;
			if (table2[j >> 6] & (BIT >> ((j & 0x03f) >> 1))) {
				step = j << 1;

				for (k = over[i]; k < last_prime; k += step) {
					pl = (k & 0x0000003f) >> 1;
					feld = k >> 6;
					table[feld] = table[feld] & ~(BIT >> pl);
				}
			}

			over[i] = k - last_prime;
		}

		feld = 0;
		pl = 0;
		maske_i = 0x80000000;

		goto st;
	}

	warning_handler("ecm_primes", "getprimes::no more primes in primetables");
	return 1;
}



void ecm_primes::
resetprimes(unsigned int us)            // set the pointer to a new prime
{
	if (us > max) {
		lidia_error_handler("ecm_primes", "resetprimes::lower bound > maxprime");
		return;
	}

	if (flag) {
		delete [] table;
		delete [] table2;
		delete [] over;

		initprimes(us, max, space);
	}
	else {
		feld_abs = us >> 6;
		feld = feld_abs;
		pl = (us & 0x0000003f) >> 1;
		maske_i = BIT >> pl;
	}

	return;
}



void ecm_primes::
killprimes()
{

	delete [] table;

	if (flag) {
		delete [] table2;
		delete [] over;
	}
}



void ecm_primes::
sieve_array(int start)               // sieving multiples of primes up to 31
{
	unsigned int i = start;

	while (i < dim_ar) {
		table[i++] = maske35 & maske7 & maske11 & maske13 & maske17 & maske19 & maske23 & maske29 & maske31;

		maske35 = maske35 << 2;
		maske35 = ((maske35 & 0xc0000000) >> 30) | maske35;

		maske7 = maske7 << 4;
		maske7 = ((maske7 & 0xf0000000) >> 28) | maske7;

		maske11 = maske11 << 10;
		maske11 = ((maske11 & 0xffc00000) >> 22) | maske11;

		maske13 = maske13 << 6;
		maske13 = ((maske13 & 0xfc000000) >> 26) | maske13;

		maske17 = maske17 << 15;
		maske17 = ((maske17 & 0xfffe0000) >> 17) | maske17;

		maske19 = maske19 << 13;
		maske19 = ((maske19 & 0xfff80000) >> 19) | maske19;

		maske23 = maske23 << 9;
		maske23 = ((maske23 & 0xff800000) >> 23) | maske23;

		maske29 = maske29 << 3;
		maske29 = ((maske29 & 0xe0000000) >> 29) | maske29;

		maske31 = maske31 << 1;
		maske31 = ((maske31 & 0x80000000) >> 31) | maske31;
	}
	return;
}



const bigmod ec_point_W::zero(0UL);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
