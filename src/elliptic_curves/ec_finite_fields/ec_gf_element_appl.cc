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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"
#include	"LiDIA/galois_field.h"
#include	"LiDIA/gf_element.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/timer.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint p;
	int deg;

	std::cout << "\nELLIPTIC CURVES OVER FINITE FIELD\n\n";
	std::cout << "\nEquation Y^2 + a1 * Y*X + a3 * Y = X^3 + a2 * X^2 + a4 * X + a6\n";

	std::cout << "\nInput charateristic of field : "; std::cin >> p;
	p = next_prime(p-1);
	std::cout << "\nChoosing next prime, p = " << p << " (" << decimal_length(p) << ")\n\n";
	std::cout << "\nInput extension degree: "; std::cin >> deg;

	galois_field F(p, deg);

	gf_element a1, a2, a3, a4, a6, x, y;

	a1.assign_zero(F);
	a2.assign_zero(F);
	a3.assign_zero(F);
	a4.assign_zero(F);
	a6.assign_zero(F);
	x.assign_zero(F);
	y.assign_zero(F);

	std::cout << "\nInput Curve Coefficients:\n";
	std::cout << "\na1 : "; std::cin >> a1;
	std::cout << "\na2 : "; std::cin >> a2;
	std::cout << "\na3 : "; std::cin >> a3;
	std::cout << "\na4 : "; std::cin >> a4;
	std::cout << "\na6 : "; std::cin >> a6;

	elliptic_curve< gf_element > e(a1, a2, a3, a4, a6);
	timer t; t.set_print_mode(HMS_MODE);

	std::cout << "\nDetermine Group Order, might take some time ... " << std::flush;

	t.start_timer();
	bigint o = e.group_order(true);
	t.stop_timer();
	std::cout << "\nGroup Order = " << o << "  (Time " << t << ") " << std::flush;
	std::cout << "\nFactorization of Order = " << factor(o) << std::flush;

	point< gf_element > P(e), Q(e);

	P = e.random_point();
	bigint xx, n1, n2;

	std::cout << "\n\nChoosing Random Point on E : P = " << P << std::flush;
	xx = order_point(P);

	std::cout << "\n\nOrder (P) = " << xx << " = " << factor(xx) << std::flush;
	assert((xx*P).is_zero());

	std::cout << "\n\nDL Algorithm:\nInput multiplier : "; std::cin >> xx;

	Q = xx*P;

	std::cout << "\nsolve DL-problem : x*P = " << Q << "\n";
	t.start_timer();
	n1 = discrete_logarithm(P, Q, true);
	std::cout << "\n\nComputed Result : x = " << n1 << std::flush;
	t.stop_timer();
	std::cout << "  (Time " << t << ")" << std::flush;

	std::cout << "\n\n\nIsomorphism Type of Curve : ";
	t.start_timer();
	e.isomorphism_type(n1, n2, P, Q);
	t.stop_timer();
	std::cout << " [ " << n1 << " , " << n2 << " ]  \n";
	std::cout << "\nGenerator for Z/" << n1 << "Z-subgroup is " << P;
        if (n2 > 1)
	  std::cout << "\nGenerator for Z/" << n2 << "Z-subgroup is " << Q;
        std::cout<<"\nTime needed : "<<t;
	std::cout << "\n\n";

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
