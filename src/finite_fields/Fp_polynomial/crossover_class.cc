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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



crossover_class Fp_polynomial::crossovers;


void spline(const int x[], const int y[], int n, double yp1, double ypn,
	    double y2[])
	//see Numerical Recipes in C, ISBN 0-521-43108-5, p. 113ff
{
	int i, k;
	double p, qn, sig, un;
	double *u = new double[n];
	memory_handler(u, "Fp_polynomial", "spline::Error in memory allocation");

	if (yp1 > 0.99e30)
		y2[0] = u[0] = 0.0;
	else {
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0])) * static_cast<double>((y[1]-y[0]) / static_cast<double>(x[1]-x[0]) - yp1);
	}

	for (i = 1; i < n-1; i++) {
		sig = static_cast<double>(x[i]-x[i-1]) / static_cast<double>(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = static_cast<double>(y[i+1]-y[i]) / static_cast<double>(x[i+1]-x[i])
			- static_cast<double>(y[i]-y[i-1]) / static_cast<double>(x[i]-x[i-1]);
		u[i] = (6.0 * u[i] / static_cast<double>(x[i+1]-x[i-1]) - sig * u[i-1]) / p;
	}
	if (ypn > 0.99e30)
		qn = un = 0.0;
	else {
		qn = 0.5;
		un = (3.0 / static_cast<double>(x[n-1]-x[n-2]))
			* (ypn - static_cast<double>(y[n-1]-y[n-2]) / static_cast<double>(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k = n-2; k >= 0; k--)
		y2[k] = y2[k]*y2[k+1]+u[k];
	delete[] u;
}



double splint(const int xa[], const int ya[], const double y2a[], int n,
	      double x)
	//see Numerical Recipes in C, ISBN 0-521-43108-5, p. 113ff
{
	int klo, khi, k;
	double h, b, a, y;

	klo = 0;
	khi = n-1;
	while (khi-klo > 1) {
		k = (khi+klo) >> 1;
		if (xa[k] > x) khi = k;
		else           klo = k;
	}
	h = xa[khi]-xa[klo];
	if (h == 0.0) {
		lidia_error_handler("crossover_class",
				    "splint(...)::identical x-values for interpolation)");
		return 0.0;
	}
	a = static_cast<double>(xa[khi]-x) / h;
	b = static_cast<double>(x-xa[klo]) / h;
	y = a*ya[klo] + b*ya[khi] +
		((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h) / 6.0;
	return y;
}



crossover_class::crossover_class()
{
//
#include	"LiDIA/Fp_pol_crossover.h"
//
	init(x_val, fftmul_val, fftdiv_val, fftrem_val, inv_val, gcd_val, xgcd_val);
	halfgcd = halfgcd_val;
	log2_newton = log2_newton_val;
}



void crossover_class::init(const int x_val[CROV_NUM_VALUES],
			   const int fftmul_val[CROV_NUM_VALUES],
			   const int fftdiv_val[CROV_NUM_VALUES],
			   const int fftrem_val[CROV_NUM_VALUES],
			   const int inv_val[CROV_NUM_VALUES],
			   const int gcd_val[CROV_NUM_VALUES],
			   const int xgcd_val[CROV_NUM_VALUES])
{
	int i;
	for (i = 0; i < CROV_NUM_VALUES; i++) {
		x[i] = x_val[i];
		y_fftmul[i] = fftmul_val[i];
		y_fftdiv[i] = fftdiv_val[i];
		y_fftrem[i] = fftrem_val[i];
		y_inv[i] = inv_val[i];
		y_gcd[i] = gcd_val[i];
	}
	spline(x, y_fftmul, CROV_NUM_VALUES, 0.0, 0.0, y2_fftmul);
	spline(x, y_fftdiv, CROV_NUM_VALUES, 0.0, 0.0, y2_fftdiv);
	spline(x, y_fftrem, CROV_NUM_VALUES, 0.0, 0.0, y2_fftrem);
	spline(x, y_inv, CROV_NUM_VALUES, 0.0, 0.0, y2_inv);
	spline(x, y_gcd, CROV_NUM_VALUES, 0.0, 0.0, y2_gcd);
	spline(x, y_xgcd, CROV_NUM_VALUES, 0.0, 0.0, y2_xgcd);
}



int crossover_class::fftmul_crossover(const bigint &modulus) const
{
	int bl = modulus.bit_length(), r;
	if (bl >= x[CROV_NUM_VALUES - 1])
		r = y_fftmul[CROV_NUM_VALUES - 1];
	else
		r = int(splint(x, y_fftmul, y2_fftmul, CROV_NUM_VALUES, bl));
	return comparator< int >::max(r, 1);
}



int crossover_class::fftdiv_crossover(const bigint &modulus) const
{
	int bl = modulus.bit_length(), r;
	if (bl >= x[CROV_NUM_VALUES - 1])
		r = y_fftdiv[CROV_NUM_VALUES - 1];
	else
		r = int(splint(x, y_fftdiv, y2_fftdiv, CROV_NUM_VALUES, bl));
	return comparator< int >::max(r, 1);
}



int crossover_class::fftrem_crossover(const bigint &modulus) const
{
	int bl = modulus.bit_length(), r;
	if (bl >= x[CROV_NUM_VALUES - 1])
		r = y_fftrem[CROV_NUM_VALUES - 1];
	else
		r = int(splint(x, y_fftrem, y2_fftrem, CROV_NUM_VALUES, bl));
	return comparator< int >::max(r, 1);
}



int crossover_class::inv_crossover(const bigint &modulus) const
{
	int bl = modulus.bit_length(), r;
	if (bl >= x[CROV_NUM_VALUES - 1])
		r = y_inv[CROV_NUM_VALUES - 1];
	else
		r = int(splint(x, y_inv, y2_inv, CROV_NUM_VALUES, bl));
	return comparator< int >::max(r, 1);
}



int crossover_class::halfgcd_crossover(const bigint &modulus) const
{
	(void)modulus;
	return halfgcd;
}



int crossover_class::gcd_crossover(const bigint &modulus) const
{
	int bl = modulus.bit_length(), r;
	if (bl >= x[CROV_NUM_VALUES - 1])
		r = y_gcd[CROV_NUM_VALUES - 1];
	else
		r = int(splint(x, y_gcd, y2_gcd, CROV_NUM_VALUES, bl));
	return comparator< int >::max(r, 1);
}



int crossover_class::xgcd_crossover(const bigint &modulus) const
{
	int bl = modulus.bit_length(), r;
	if (bl >= x[CROV_NUM_VALUES - 1])
		r = y_xgcd[CROV_NUM_VALUES - 1];
	else
		r = int(splint(x, y_xgcd, y2_xgcd, CROV_NUM_VALUES, bl));
	return comparator< int >::max(r, 1);
}



int crossover_class::log2_newton_crossover(const bigint &modulus) const
{
	(void)modulus;
	return log2_newton;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
