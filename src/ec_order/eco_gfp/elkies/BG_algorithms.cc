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
//	Author	:Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"
#include	"LiDIA/wep_rat_function.h"
#ifdef DEBUG
#include	<cassert>
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//-----------------------------------------------------------------
//  Let M be the matrix whose columns are given as a(x) * x^i mod F
//  for i = 0, ..., F.degree()-1. The function determines the first
//  row of this matrix and returns it as vector v.


void eco_prime::
mult_matrix(base_vector< bigint > & v, const ff_pol & aa,
	    const ff_polmod & F)
{
	int i, d;
	ff_pol a(aa);

	d = F.modulus().degree();
	v[0].assign(a.const_term());

	for (i = 1; i < d; i++) {
		multiply_by_x_mod(a, a, F.modulus());
		v[i].assign(a.const_term());
	}
}



//******************************************************************
// This functions implements the fast search from the paper
// Vmueller and Mmaurer: "Fast ev search". The function inner_product
// is already available in Fp_polynomial.


int eco_prime::search_in_table(const ff_rat & f,
			       ff_rat * table, const ff_polmod & F,
			       int size)
{
	int d = F.modulus().degree();
	base_vector< bigint > v1(d, d), v2(d, d);
	bigint r1, r2;
	ff_element res1, res2;

	mult_matrix(v1, f.numerator(), F);
	mult_matrix(v2, f.denominator(), F);

	for (int i = 0; i < size; i++) {
		inner_product(r1, v1, table[i].denominator());
		inner_product(r2, v2, table[i].numerator());

		res1.assign(r1);
		res2.assign(r2);

#ifdef DEBUG
		if (equal(f, table[i], F)) {
			std::cout << "\nDEBUG: found match for i = " << i << std::flush;
			if (res1 == res2)
				std::cout << " --> Test also " << std::flush;
			else
				std::cout << "\n --> Test NOT " << std::flush;
		}
#endif

		if (res1 == res2) {
			// pseudo_match found
			if (equal(f, table[i], F))
				return i;
		}
	}
	return -1;
}



//---------------------------------------------------------------------
//  the Babystep Giantstep version of schoofpart

bool eco_prime::schoofpart_BG (lidia_size_t & ev, const ff_pol & fC,
			       lidia_size_t* klist)
{
	if (l < lower_bound_for_BG) {
		// don't use BG for such small primes !!
		char strat = ev_strategy;
		bool rc;

		ev_strategy = eco_prime::EV_DIVISION_POLYNOMIAL;
		rc = schoofpart(ev, fC, klist, false);
		ev_strategy = strat;
		return rc;
	}

	lidia_size_t   i, j, d;
	lidia_size_t   number_BG;
	ff_pol  xp, yp; // Y^(q-1) = (X^3+aX+b)^((q-1)/2)  and  X^q mod fC
	ff_polmod ffC;

	ffC.build (fC);
	// compute number of babysteps/giantsteps first

	d = klist[klist[0]]; // the maximal element
	if (l-klist[1] > static_cast<unsigned int>(d))   // take care of the sign !!
		d = l - klist[1];

	number_BG = static_cast<lidia_size_t>(std::sqrt(static_cast<double>(d) / 2.0)) + 1;

	while(2*number_BG*(number_BG + 1) < d)
		number_BG ++;

	if (info) {
		std::cout << "\nTesting " << d;
		if (d == 1) {
			std::cout << " possibility using Babystep Giantstep Algorithm\n(";
			std::cout << number_BG << " Baby- / Giantsteps) :  " << std::flush;
		}
		else {
			std::cout << " possibilities using Babystep Giantstep Algorithm\n(";
			std::cout << number_BG << " Baby-/Giantsteps) :  " << std::flush;
		}
	}

	xp.set_modulus(pn);
	yp.set_modulus(pn);

#ifdef TIMING
	timer t;
	t.set_print_mode();
	t.start_timer();
#endif

	power_x (xp, pn, ffC); // compute frob(P) = (X, Y)^p
	Ytop_f (yp, ffC);

#ifdef TIMING
	t.stop_timer();
	if (info)
		std::cout << "\nComputation of Phi(X, Y) needs time : " << t << std::flush;
#endif

	if (xp.is_x()) {
		// match for k = +/- 1.
		if (yp.is_one()) {
			ev = 1;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is 1" << std::flush;
			return true;
		}
		else {
			ev = -1;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is -1" << std::flush;
			return true;
		}
	}

	// initialize the formal point computations

	wep_rat_function::initialize(eco_prime::A, eco_prime::B, ffC);

	wep_rat_function P, frob_P, kP;
	ff_rat *table;

	frob_P.assign(xp, yp);

#ifdef DEBUG
	wep_rat_function frob(frob_P), H;
#endif

	P.assign_xy();
	kP.assign_xy();

	if (info)
		std::cout << "B " << std::flush;

	// start preparing the babystep table

	table = new ff_rat [number_BG]; // table[i] = (i+1) * P
	table[0].assign(P.get_y());
	if (info)
		std::cout << "B " << std::flush;

	multiply_by_2(kP, kP);
	table[1].assign(kP.get_y());
	if (info)
		std::cout << "B " << std::flush;

	if (equal(frob_P.get_x(), kP.get_x(), ffC)) {
		if (equal(frob_P.get_y(), kP.get_y(), ffC)) {
			ev = 2;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}

		if (equal(frob_P.get_y(), -kP.get_y(), ffC)) {
			ev = -2;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}
	}

        for (i = 3; i <= number_BG; i++)
          {
            add_P_XY(kP, kP);
	    
            if (equal(frob_P.get_x(), kP.get_x(), ffC)) {
              if (equal(frob_P.get_y(), kP.get_y(), ffC)) {
                ev = i;
                if (info)
                  std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
                delete[] table;
                return true;
              }
              if (equal(frob_P.get_y(), -kP.get_y(), ffC)) {
                ev = -i;
                if (info)
                  std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
                delete[] table;
                return true;
              }
            }
           
            table[i-1].assign(kP.get_y());
            if (info)
              std::cout << "B " << std::flush;
          }

        if (info)
          std::cout << "  ";

        subtract(frob_P, frob_P, kP);
        multiply_by_2(kP, kP);
        negate(kP, kP);

	for (j = 0; j <= number_BG; j++) {
		// iteration j:  frob_P = (X, Y)^p - (2*j+1)
		// * number_BG * (X, Y)
		if (info)
			std::cout << "G " << std::flush;

		if (frob_P.is_zero()) {
			ev = (2*j+1) * number_BG;

#ifdef DEBUG
			multiply(H, ev, P);
			assert(frob == H);
#endif

			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}
	
		i = search_in_table(frob_P.get_y(), table, ffC, number_BG);

		if (i != -1) {
			// match found: frob_P = (i+1) * (X, Y)
			ev = (2*j+1) * number_BG + i + 1;

#ifdef DEBUG
			P.assign_xy();
			multiply(H, ev, P);
			assert(frob == H);
#endif

			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}

		i = search_in_table(-(frob_P.get_y()), table, ffC, number_BG);

		if (i != -1) {
			// match found: frob_P = -(i+1) * (X, Y)
			ev = (2*j+1) * number_BG - i - 1;

#ifdef DEBUG
			P.assign_xy();
			multiply(H, ev, P);
			assert(frob == H);
#endif

			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}

		add(frob_P, frob_P, kP);
	}

	delete[] table;
	lidia_error_handler("eco_prime", "schoofpart_BG::found no solution");
	return false;
}



//********************************************************************
//
// Now all the functions for the funny Babystep Giantstep idea are coming.
//


//--------------------------------------------------------------------
// The following function doubles the point with x-coordinate x1 and
// returns only the x-coordinate of the result (using formulas for
// projective coordinates)

void eco_prime::double_point_x_only(ff_rat & x3, const ff_rat & x1,
				    const ff_polmod &F)
{
	ff_pol x(x1.numerator()), z(x1.denominator());
	ff_pol x2, z2, xz, n, d;

	x2.set_modulus(A.modulus());
	z2.set_modulus(x2); xz.set_modulus(x2); n.set_modulus(x2);
	d.set_modulus(x2);

	n.assign_zero(); d.assign_zero();

	square(x2, x, F);
	square(z2, z, F);
	multiply(xz, x, z, F);

	add(d, x2, A.mantissa()*z2);
	multiply(d, d, xz, F);

	add(n, (-2 * A.mantissa())*xz, (-8 * B.mantissa()) * z2);
	multiply(n, n, xz, F);

	square(x2, x2, F);
	square(z2, z2, F);

	add(d, d, B.mantissa()* z2);
	multiply_by_scalar(d, d, 4);
	add(n, n, x2);
	add(n, n, (A*A).mantissa() * z2);
	x3.assign(n, d);
}



//------------------------------------------------------------------
// Since the cost for Babysteps are lower than the cost for Giantsteps
// (double_point_x_only needs 4 squarings and 3 mults --> weight 5,
// Giantsteps with divpol:  7 mults, 1 squarings --> weight 7.5)
// we use a weighted formula for computing costs.

void eco_prime::
find_number_of_FBG_steps(lidia_size_t & baby, lidia_size_t & giant,
			 double & average,
			 lidia_size_t *klist)
{
	int i, j;
	float min;
	int b_tmp, g_tmp;
	lidia_size_t alpha;
	const float weight_B = 5.0;
	const float weight_G = 7.5;

	baby = giant = 0;

	for (i = 1; i <= klist[0]; i++) {
		// compute maximal exponent 2^j*frob
		// if alpha = klist[i], set b_tmp
		min = weight_G * (klist[i] - 1);
		alpha = klist[i];
		b_tmp = 0;
		g_tmp = alpha;

		for (j = 1; j < static_cast<int>(l); j ++) {
			alpha = (alpha << 1) % l;
			alpha = comparator< lidia_size_t >::min(alpha, l-alpha);
			if (weight_B * j + weight_G * alpha  < min) {
				// 2^j * ev = alpha mod l
				min = weight_B * j + weight_G * alpha;
				b_tmp = j;
				g_tmp = alpha;
			}
		}

		if (b_tmp > baby)
			baby = b_tmp;

		if (g_tmp > giant)
			giant = g_tmp;
	}

	average = 0.0;
	giant = 0;

	for (i = 1; i <= klist[0]; i++) {
		alpha = klist[i];
		g_tmp = alpha;

		for (j = 1; j <= baby; j++) {
			alpha = (alpha << 1) % l;
			alpha = comparator< lidia_size_t >::min(alpha, l-alpha);
			if (alpha < g_tmp)
				g_tmp = alpha;
		}
		if (g_tmp > giant)
			giant = g_tmp;
		average += static_cast<double>(weight_B* baby + weight_G * g_tmp);
	}
	average /= static_cast<double>(klist[0]);
	average /= static_cast<double>(weight_G);
}



bool eco_prime::schoofpart_FBG (lidia_size_t & ev, const ff_pol & fC,
				lidia_size_t* klist)
{
	lidia_size_t   i, j;
	lidia_size_t   number_baby, number_giant;
	double average;
	ff_pol  xp_pol; // X^q mod fC
	ff_rat  xp; // and as rational function
	ff_polmod ffC;

	if (l < lower_bound_for_FBG) {
		// don't use FBG for such small primes !!
		char strat = ev_strategy;
		bool rc;

		ev_strategy = eco_prime::EV_DIVISION_POLYNOMIAL;
		rc = schoofpart(ev, fC, klist, false);
		ev_strategy = strat;
		return rc;
	}

	ffC.build (fC);
	find_number_of_FBG_steps(number_baby, number_giant, average, klist);

	if (info) {
		std::cout << "\nTesting " << klist[0];
		if (klist[0] == 1) {
			std::cout << " possibility using Funny Babystep Giantstep Algorithm\n(";
			std::cout << number_baby << " Babysteps, " << number_giant << " Giantsteps, ";
			std::cout << "expected cost : " << average;
			std::cout << " additions) :  " << std::flush;
		}
		else {
			std::cout << " possibilities using Funny Babystep Giantstep Algorithm\n(";
			std::cout << number_baby << " Babysteps, " << number_giant << " Giantsteps, ";
			std::cout << average << " expected additions) :   " << std::flush;
		}
	}

	xp_pol.set_modulus(pn);

#ifdef TIMING
	timer t;
	t.set_print_mode();
	t.start_timer();
#endif

	power_x (xp_pol, pn, ffC);

#ifdef TIMING
	t.stop_timer();
	if (info)
		std::cout << "\nComputing X^q mod f_C(X) needs time : " << t << std::flush;
#endif

	xp.set_modulus(pn);
	xp.assign(xp_pol);

	if (xp_pol.is_x()) {
		// match for k = +/- 1.
		ev = 1;
		if (info)
			std::cout << "\n\nEigenvalue of Frobenius is +/- 1" << std::flush;
		return true;
	}

	PsiPowers* psi_pow;
	ff_pol  sqrYY, right_side; // y^4 = (x^3+ax+b)^2 and y^2
	ff_pol num, den, help;
	ff_rat *table;
	ff_element inv2;
	int top, deleted = -1, kk;

	sqrYY.set_modulus(p);
	right_side.set_modulus(p);
	CurveEqn (right_side, p, A, B, ffC);
	square   (sqrYY, right_side, ffC);
	invert (inv2, ff_element(2));

	init_psi_comp (psi_pow, top, number_giant + 1, ffC);

	while (top <= number_giant + 1) {
		next_psi (psi_pow, top, ffC, sqrYY, inv2);
		if (top & 1)
			kk = ((top+1) >> 1)-3;
		else
			kk = ((top) >> 1)-2;

		while (deleted <= kk) {
			deleted ++;
			psi_pow[deleted].pow3->kill();
		}
	}

	// now I'm testing first whether (X^q, .) = alpha *(X, Y) with
	// alpha even, 2 <= alpha <= number_giant

	num.set_modulus(p);
	den.set_modulus(p);
	help.set_modulus(p);

	for (i = 2; i <= number_giant; i += 2) {
		multiply(den, *(psi_pow[i].pow2), right_side, ffC);
		multiply_by_x_mod(num, den, ffC.modulus());
		multiply(help, *(psi_pow[i-1].pow1), *(psi_pow[i+1].pow1), ffC);
		subtract(num, num, help);

		if (equal(xp, ff_rat(num, den), ffC)) {
			ev = i;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is +/- " << ev << std::flush;
			free_psi(psi_pow, 0, number_giant+4);
			return true;
		}
	}

	// next we determine a table of babysteps as
	// table[i] = x-coord(2^j * (X^q, .))
	
	table = new ff_rat [number_baby+1];
	table[0].assign(xp);

	for (i = 1; i <= number_baby; i++) {
		double_point_x_only(xp, xp, ffC);
		table[i].assign(xp);
		if (info)
			std::cout << "B " << std::flush;
	}

	// now do Giantsteps, but necessary only for odd multipliers !!

	udigit uinv;

	if (info)
		std::cout << "  ";

	for (j = 1; j <= number_giant; j += 2) {
		if (info)
			std::cout << "G G " << std::flush;

		multiply_by_x_mod(num, *(psi_pow[j].pow2), ffC.modulus());
		multiply(help, *(psi_pow[j-1].pow1), *(psi_pow[j+1].pow1), ffC);
		multiply(help, help, right_side, ffC);
		subtract(num, num, help);

		i = search_in_table(ff_rat(num, *(psi_pow[j].pow2)), table,
				    ffC, number_baby+1);

		if (i != -1) {
			// match found
			uinv = power_mod(static_cast<udigit>(2),
						static_cast<udigit>(i),
						static_cast<udigit>(l));
			uinv = invert_mod(uinv, static_cast<udigit>(l));
			ev = static_cast<int>(uinv * j) % l;
			delete[] table;
			free_psi(psi_pow, 0, number_giant+4);
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is +/- " << ev << std::flush;
	
			return true;
		}
	}

	delete[] table;
	free_psi(psi_pow, 0, number_giant+4);
	lidia_error_handler("eco_prime", "FBG::nothing found");
	return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
