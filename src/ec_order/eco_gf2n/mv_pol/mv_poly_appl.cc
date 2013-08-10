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
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void identitytest(const mv_poly & a, const mv_poly & b, const mv_poly & c)
{
	mv_poly d;

	d = a;
	assert(d == a);
	assert(-(-a) ==  a);
	assert((a + b) ==  (b + a));
	assert((a + (-b)) ==  (a - b));
	assert((a * b) ==  (b * a));
	assert((a * (-b)) ==  -(a * b));
	assert((a - b) ==  -(b - a));
	assert((a + (b + c)) ==  ((a + b) + c));
	assert((a * (b * c)) ==  ((a * b) * c));
	assert((a * (b + c)) ==  ((a * b) + (a * c)));
	assert(((a - b) + b) ==  a);
	assert(((a + b) - b) ==  a);
}



int main_LiDIA(int argc, char** argv)
{
	int n, k, i;
	int anz_terms;
	int anz;

	std::cout << "\nTest program for mv_poly : \n";
	std::cout <<   "==========================\n\n";
	std::cout << "\nPlease input degree n of field GF(2^n) : ";
	std::cin >> n;
	gf2n_init(n);

	mv_poly p1, p2, p3, p4;
	bigint x, bound;
	gf2n r;

	std::cout << "\nPlease input maximal number of variables : "; std::cin >> k;

	shift_left(bound, bigint(1), k);
	dec(bound);

	std::cout << "\nPlease input maximal number of terms : "; std::cin >> anz_terms;
	std::cout << "\nPlease enter number of tests : "; std::cin >> anz;
	std::cout << "\n\nChoose random polynomials ... " << std::flush;

	for (i = 0; i < anz_terms; i++) {
		r.randomize();
		x.randomize(bound);
		p1 += mv_term(r, x);
	}

	for (i = 0; i < anz_terms; i++) {
		r.randomize();
		x.randomize(bound);
		p2 += mv_term(r, x);

		r.randomize();
		x.randomize(bound);
		p3 += mv_term(r, x);
	}

	std::cout << "\nStarting addition tests ... " << std::flush;

	for (i = 0; i < anz; i++) {
		p3 = p1 + p2;
		p4 = p2 + p1;
		assert(p3 == p4);
		p2 = p3;
		identitytest(p1, p2, p3);
	}

	std::cout << "\nStarting multiplication tests ... " << std::flush;
	gf2n xx;

	for (i = 0; i < anz; i++) {
		p3 = p1 * p2;
		p4 = p2 * p1;
		assert(p3 == p4);

		p1 = p3;
		p3 = p1 * p2;
		p1 *= p2;
		assert(p3 == p1);


		square(p2, p3);
		p4 = p3 * p3;
		assert(p2 == p4);

		p3.assign_zero();
		xx.randomize();
		p3 = mv_term(xx);

		multiply(p1, p3, p2);
		multiply(p2, xx, p2);
		assert(p2 == p1);
	}
	mv_poly::delete_free_list();
	std::cout << "\n==> DONE.\n\n";
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
