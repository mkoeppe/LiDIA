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
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/Fp_poly_multiplier.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void randomize_bitlength(bigint& x, lidia_size_t num_bits)
{
	// x = "random number" with precisely num_bits bits.
	x.assign_zero();

	bigint tmp;
	tmp.assign_one();
	shift_left(tmp, tmp, num_bits); //tmp = 2^num_bits;

	do {
		x = randomize(tmp);
	} while (x.bit_length() != num_bits);

	return;
}



int main_LiDIA(int, char**)
{
	int n, k;

	n = 200;
	k = 200;

	bigint p, t1, t2;

	// p is computed as a random k-bit number.
	// It is not prime, but who cares.

	randomize_bitlength(t1, k);
	t2.assign_one();
	shift_left(t2, t2, k);
	add(p, t1, t2);

	std::cout << "multiplication of two polynomials modulo another polynomial "
		"over Z/pZ\n";
	std::cout << "  degree of the polynomials = " << n << std::endl;
	std::cout << "  modulus p = " << p << " (about " << k << " bits)\n\n";


	Fp_polynomial f, g, h, r1, r2, r3, r4;

	randomize(g, p, n-1); // g = random polynomial of degree < n
	randomize(h, p, n-1); // h = "   "
	randomize(f, p, n-1); // f = "   "
	f.set_coefficient(n); // Sets coefficient of X^n to 1


	// For doing arithmetic mod f quickly, one must pre-compute
	// some information.

	Fp_poly_modulus F(f);
	timer t;

	std::cout << "classical algorithm : \t\t\t\t\t\t" << std::flush;
	t.start_timer();
	plain_mul(r1, g, h); // uses classical arithmetic
	plain_rem(r1, r1, f);
	t.stop_timer();
	std::cout << t.user_time() << " hsec" << std::endl;

	std::cout << "FFT algorithm : \t\t\t\t\t\t" << std::flush;
	t.start_timer();
	multiply_mod(r3, g, h, f); // uses FFT
	t.stop_timer();
	std::cout << t.user_time() << " hsec" << std::endl;

	std::cout << "multiplication with Fp_poly_modulus : \t\t\t\t" << std::flush;
	t.start_timer();
	multiply(r2, g, h, F); // uses Fp_poly_modulus
	t.stop_timer();
	std::cout << t.user_time() << " hsec" << std::endl;

	std::cout << "multiplication with Fp_poly_modulus, Fp_poly_multiplier : \t"
		  << std::flush;
	Fp_poly_multiplier H(h, F); // we don't count precomputations
	t.start_timer();
	multiply(r4, g, H, F); // uses Fp_poly_modulus and Fp_poly_multiplier
	t.stop_timer();
	std::cout << t.user_time() << " hsec" << std::endl;


	// compare the results...

	if (r1 != r2)
		std::cout << "r1 != r2!!\n";
	else if (r1 != r3)
		std::cout << "r1 != r3!!\n";
	else if (r1 != r4)
		std::cout << "r1 != r4!!\n";
	else
		std::cout << "results are OK.\n";

	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
