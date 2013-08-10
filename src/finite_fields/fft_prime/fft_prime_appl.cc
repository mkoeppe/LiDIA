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
//	Author	: Thorsten Rottschaefer (TR)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/fft_prime.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	lidia_size_t i;
	fft_prime fft;

	std::cout << "\ntest 1 : " << std::flush;
	// test fft. result should be:
	// 1 4 10 3 8 7 16 0

	fft.set_prime(17);

	fft_prime_t y[8], c[4], d[4], yy[8];

	//pa=1  + 2x^1  + 3x^2  + 4x^3
	c[0] = 1; c[1] = 2; c[2] = 3; c[3] = 4;

	//pb=1  + 2x^1  + 3x^2  + 4x^3
	d[0] = 1; d[1] = 2; d[2] = 3; d[3] = 4;

	//expected result:
	yy[0] = 1; yy[1] = 4; yy[2] = 10; yy[3] = 3; yy[4] = 8; yy[5] = 7; yy[6] = 16; yy[7] = 0;

	//std::cout << "\n\nsmall test of fft_prime\n\n";
	//std::cout << "  multiplying pa=1 + 2x^1 + 3x^2 + 4x^3\n";
	//std::cout << "          and pb=1 + 2x^1 + 3x^2 + 4x^3\n\n";

	multiply_fft(static_cast<void*>(y), 0, static_cast<void*>(c), 3, static_cast<void*>(d), 3, fft);
	// a and b have degree 3
	// the second parameter "0" means that y,c,d are of type fft_prime_t

	//std::cout << "     result is : ";
	//for (i = 0; i < 7; i++)
	//	std::cout << y[i] << "x^" << i << " + ";

	bool ok = true;
	for (i = 0; i < 8; i++)
		if (y[i] != yy[i]) ok = false;

	if (ok)
		std::cout << "passed.\n";
	else
		std::cout << "ERROR: not passed.\n";


	std::cout << "test 2 : " << std::flush;
	// test fft now for bigger polynom. result should be:
	// 1 2 3 .... 998 999 1000 999 998 .... 3 2 1

	fft.set_prime(12289);
	udigit_mod::set_modulus(12289);

	udigit_mod* x = new udigit_mod[2048];
	udigit_mod* a = new udigit_mod[1000];
	udigit_mod* b = new udigit_mod[1000];

	for (i = 0; i < 1000; i++)
		a[i] = b[i] = 1;

	multiply_fft(static_cast<void*>(x), 1, static_cast<void*>(a), 999, static_cast<void*>(b), 999, fft);
	// a and b have degree 999
	// the second parameter "1" means that y,c,d are of type udigit_mod

	ok = true;
	for (i = 0; i < 1000; i++)
		if (x[i] != i+1) ok = false;
	for (i = 0; i < 1000; i++)
		if (x[1000+i] != 999 - i) ok = false;
	for (i = 2000; i < 2048; i++)
		if (x[i] != 0) ok = false;

	if (ok)
		std::cout << "passed.\n";
	else
		std::cout << "ERROR: not passed.\n";

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
