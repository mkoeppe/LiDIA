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


#include	"LiDIA/gf_element.h"
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/gf2n.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define LOOP1(a) t1.start_timer(); for (w = 0; w < wdh; w++) {a; } t1.stop_timer();
#define LOOP2(a) t2.start_timer(); for (w = 0; w < wdh; w++) {a; } t2.stop_timer();
#define PRINT_RES(text) std::cout << text << " \t\t" << t1.user_time() << " hs\t\t"\
 << t2.user_time() << " hs\t\t(" << wdh << ")\n";

#define ERROR(a, A) std::cout << "an error occurred : " << a << " != " << A << std::endl;

#define CHECK1(a, A) { if (Fp_polynomial(A.mantissa(), A.modulus()) != a.polynomial_rep()) ERROR(a, A); }
#define CHECK2(a, A) { if (to_gf2n(a) != A) ERROR(a, A); }
#define CHECK3(a, A) { if (A != a.polynomial_rep()) ERROR(a, A); }



gf2n to_gf2n(const gf_element &a)
{
	Fp_polynomial pol = a.polynomial_rep();
	bigint b;
	int i;
	for (i = pol.degree(); i > 0; i--) {
		if (pol[i] == 1)
			b += 1;
		b <<= 1;
	}
	if (pol[0] == 1)
		b += 1;
	return gf2n(b);
}



void case1()
{
	std::cout << "gf_element < ->bigmod" << std::endl
		  << "---------------------" << std::endl;

	timer t1, t2;
	long w, wdh, WDH;

	bigint p;
	std::cout << "Enter prime modulus : ";
	std::cin >> p;
	std::cout << "Enter number of reps : ";
	std::cin >> WDH;

	bigmod::set_modulus(p);
	bigmod A, B, C;

	galois_field K(p);
	gf_element a(K), b(K), c(K);

	//assign random elts.
	A.randomize();
	B.randomize();
	a.assign(A.mantissa());
	b.assign(B.mantissa());

	std::cout << std::endl;
	std::cout << "operation:\t\tgf_element\tbigmod\t\t#loops" << std::endl;

	wdh = WDH;
	LOOP1(add(c, a, b));
	LOOP2(add(C, A, B));
	PRINT_RES("add     ");
	CHECK1(c, C);

	LOOP1(subtract(c, a, b));
	LOOP2(subtract(C, A, B));
	PRINT_RES("subtract");
	CHECK1(c, C);

	wdh = WDH/10;
	LOOP1(multiply(c, a, b));
	LOOP2(multiply(C, A, B));
	PRINT_RES("multiply");
	CHECK1(c, C);

	wdh = WDH/100;
	LOOP1(divide(c, a, b));
	LOOP2(divide(C, A, B));
	PRINT_RES("divide ");
	CHECK1(c, C);

	std::cout << std::endl;
}



//---------------------------------------------------------------------------

void case2()
{
	std::cout << "gf_element < ->gf2n" << std::endl
		  << "---------------------" << std::endl;

	long deg, wdh, w, WDH;
	timer t1, t2;
	std::cout << "Enter degree : ";
	std::cin >> deg;
	std::cout << "Enter number of reps : ";
	std::cin >> WDH;

	gf2n_init(deg);
	gf2n A, B, C;

	galois_field K(2, deg);
	gf_element a(K), b(K), c(K);

	//assign random elts.
	a.randomize();
	b.randomize();
	A = to_gf2n(a);
	B = to_gf2n(b);

	std::cout << std::endl;
	std::cout << "operation:\t\tgf_element\tgf2n\t\t#loops" << std::endl;

	wdh = WDH;
	LOOP1(add(c, a, b));
	LOOP2(add(C, A, B));
	PRINT_RES("add     ");
	CHECK2(c, C);

	LOOP1(subtract(c, a, b));
	LOOP2(subtract(C, A, B));
	PRINT_RES("subtract");
	CHECK2(c, C);

	wdh = WDH/10;
	LOOP1(multiply(c, a, b));
	LOOP2(multiply(C, A, B));
	PRINT_RES("multiply");
	CHECK2(c, C);

	wdh = WDH/100;
	LOOP1(divide(c, a, b));
	LOOP2(divide(C, A, B));
	PRINT_RES("divide ");
	CHECK2(c, C);

	std::cout << std::endl;
}



//---------------------------------------------------------------------------

void case3()
{
	std::cout << "gf_element < ->Fp_polynomial" << std::endl
		  << "----------------------------" << std::endl;

	long w, wdh, WDH;
	timer t1, t2;
	bigint p;
	int deg;
	std::cout << "Enter prime modulus : ";
	std::cin >> p;
	std::cout << "Enter degree : ";
	std::cin >> deg;
	std::cout << "Enter number of reps : ";
	std::cin >> WDH;

	galois_field K(p, deg);
	gf_element a(K), b(K), c(K);

	Fp_polynomial A, B, C, D;
	Fp_poly_modulus F(K.irred_polynomial());

	A.set_modulus(p);
	B.set_modulus(p);
	C.set_modulus(p);
	D.set_modulus(p);

	//assign random elts.
	a.randomize();
	b.randomize();
	A = a.polynomial_rep();
	B = b.polynomial_rep();

	std::cout << std::endl;
	std::cout << "operation:\t\tgf_element\tFp_polynomial\t#loops" << std::endl;

	wdh = WDH;
	LOOP1(add(c, a, b));
	LOOP2(add(C, A, B));
	PRINT_RES("add     ");
	CHECK3(c, C);

	LOOP1(subtract(c, a, b));
	LOOP2(subtract(C, A, B));
	PRINT_RES("subtract");
	CHECK3(c, C);

	wdh = WDH/10;
	LOOP1(multiply(c, a, b));
	LOOP2(multiply(C, A, B, F));
	PRINT_RES("multiply");
	CHECK3(c, C);

	wdh = WDH/100;
	LOOP1(divide(c, a, b));
	LOOP2(invert_mod(D, B, K.irred_polynomial()); multiply(C, A, D, F));
	PRINT_RES("divide ");
	CHECK3(c, C);

	std::cout << std::endl;
}



//---------------------------------------------------------------------------

int main_LiDIA(int argc, char** argv)
{
	int x;
	do {
		std::cout << "gf_element timings" << std::endl
			  << "==================" << std::endl;
		std::cout << "  gf_element < ->bigmod  ...............  (1)" << std::endl;
		std::cout << "  gf_element < ->gf2n  .................  (2)" << std::endl;
		std::cout << "  gf_element < ->Fp_polynomial  ........  (3)" << std::endl;
		std::cout << "  quit  ................................  (0)" << std::endl;
		std::cout << std::endl;
		std::cout << "Enter your choice : ";
		std::cin >> x;
		std::cout << std::endl;

		switch(x) {
		case 1:	case1(); break;
		case 2:	case2(); break;
		case 3:	case3(); break;
		};

	} while (x != 0);

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
