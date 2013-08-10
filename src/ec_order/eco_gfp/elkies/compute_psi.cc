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
//	Author	: Frank Lehmann, Markus Maurer, Peter Noss, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
//
//  Compute all division polynomials up to a certain bound
//
//

//----------------------------------------------------------------
// the following function allocates memory for a list of division
// polynomials psi[0], ..., psi[nmax] modulo the ff_polmod f.
// top is set to maximally already computed division polynomial
// triple.
// YY = Y^2 = (X^3+aX+b) mod f
// odd degree div pols are correct, even degree div polynomials are
// without one (Y, Y^2, Y^3), respectively

void eco_prime::init_psi_comp (PsiPowers* & psi_pow,
			       lidia_size_t &  top,
			       lidia_size_t nmax,
			       const ff_polmod & f)
{
	lidia_size_t i;
	ff_element  tmp;
	ff_element  tmp2;

	psi_pow = new PsiPowers[nmax+4];

	memory_handler (psi, "eco_prime::init_psi_comp()",
			"Allocating psi_pow");

	if (nmax >= 2)
		for (i = 5; i <= nmax+3; i++) {
			(psi_pow[i]).pow1 = new ff_pol;
			(psi_pow[i]).pow2 = new ff_pol;
			(psi_pow[i]).pow3 = new ff_pol;

			(psi_pow[i]).pow1->set_modulus(p);
			(psi_pow[i]).pow2->set_modulus(p);
			(psi_pow[i]).pow3->set_modulus(p);
		}

	// psi_0, psi_0^2, psi_0^3

	(psi_pow[0]).pow1 = new ff_pol;
	(psi_pow[0]).pow2 = new ff_pol;
	(psi_pow[0]).pow3 = new ff_pol;

	(psi_pow[0]).pow1->set_modulus(p);
	(psi_pow[0]).pow2->set_modulus(p);
	(psi_pow[0]).pow3->set_modulus(p);

	(*(psi_pow[0].pow1)).assign_zero();
	(*(psi_pow[0].pow2)).assign_zero();
	(*(psi_pow[0].pow3)).assign_zero();

	// psi_1, psi_1^2, psi_1^3

	(psi_pow[1]).pow1 = new ff_pol;
	(psi_pow[1]).pow2 = new ff_pol;
	(psi_pow[1]).pow3 = new ff_pol;

	(psi_pow[1]).pow1->set_modulus(p);
	(psi_pow[1]).pow2->set_modulus(p);
	(psi_pow[1]).pow3->set_modulus(p);

	(*(psi_pow[1].pow1)).assign_one();
	(*(psi_pow[1].pow2)).assign_one();
	(*(psi_pow[1].pow3)).assign_one();

	// psi_2, psi_2^2, psi_2^3 

	psi_pow[2].pow1 = new ff_pol;
	psi_pow[2].pow2 = new ff_pol;
	psi_pow[2].pow3 = new ff_pol;


	psi_pow[2].pow1->set_modulus(p);
	psi_pow[2].pow2->set_modulus(p);
	psi_pow[2].pow3->set_modulus(p);

	tmp = ff_element(2);
	psi_pow[2].pow1->set_coefficient (tmp.mantissa(), 0);

	tmp = ff_element(4);
	psi_pow[2].pow2->set_coefficient (tmp.mantissa(), 0);

	tmp = ff_element(8);
	psi_pow[2].pow3->set_coefficient (tmp.mantissa(), 0);


	// psi_3, psi_3^2, psi_3^3

	psi_pow[3].pow1 = new ff_pol;
	psi_pow[3].pow2 = new ff_pol;
	psi_pow[3].pow3 = new ff_pol;


	psi_pow[3].pow1->set_modulus(p);
	psi_pow[3].pow2->set_modulus(p);
	psi_pow[3].pow3->set_modulus(p);

	tmp = ff_element(3);
	psi_pow[3].pow1->set_coefficient (tmp.mantissa(), 4);
	tmp = ff_element(0);
	psi_pow[3].pow1->set_coefficient (tmp.mantissa(), 3);
	multiply (tmp2, ff_element(6), A);
	psi_pow[3].pow1->set_coefficient (tmp2.mantissa(), 2);
	multiply (tmp2, ff_element(12), B);
	psi_pow[3].pow1->set_coefficient (tmp2.mantissa(), 1);
	square (tmp, A);
	negate (tmp, tmp);
	psi_pow[3].pow1->set_coefficient (tmp.mantissa(), 0);

	remainder (*(psi_pow[3].pow1), *(psi_pow[3].pow1), f);
	square (*(psi_pow[3].pow2), *(psi_pow[3].pow1), f);
	multiply (*(psi_pow[3].pow3), *(psi_pow[3].pow2), *(psi_pow[3].pow1), f);

	// psi_4, psi_4^2, psi_4^3

	psi_pow[4].pow1 = new ff_pol;
	psi_pow[4].pow2 = new ff_pol;
	psi_pow[4].pow3 = new ff_pol;


	psi_pow[4].pow1->set_modulus(p);
	psi_pow[4].pow2->set_modulus(p);
	psi_pow[4].pow3->set_modulus(p);

	tmp =  ff_element(4);
	psi_pow[4].pow1->set_coefficient (tmp.mantissa(), 6);

	tmp =  ff_element(0);
	psi_pow[4].pow1->set_coefficient (tmp.mantissa(), 5);

	multiply (tmp2, ff_element(20), A);
	psi_pow[4].pow1->set_coefficient (tmp2.mantissa(), 4);

	multiply (tmp2, ff_element(80), B);
	psi_pow[4].pow1->set_coefficient (tmp2.mantissa(), 3);

	multiply (tmp, A, B);
	multiply (tmp, tmp, ff_element(16));
	negate (tmp2, tmp);
	psi_pow[4].pow1->set_coefficient (tmp2.mantissa(), 1);

	square(tmp, A);
	multiply (tmp2, tmp, ff_element(20));
	negate (tmp2, tmp2);
	psi_pow[4].pow1->set_coefficient (tmp2.mantissa(), 2);

	multiply    (tmp, tmp, A);
	square (tmp2, B);
	multiply(tmp2, tmp2, ff_element(8));
	add    (tmp, tmp, tmp2);
	multiply(tmp, tmp, ff_element(4));
	negate (tmp2, tmp);
	psi_pow[4].pow1->set_coefficient (tmp2.mantissa(), 0);

	remainder (*(psi_pow[4].pow1), *(psi_pow[4].pow1), f);
	square (*(psi_pow[4].pow2), *(psi_pow[4].pow1), f);
	multiply (*(psi_pow[4].pow3), *(psi_pow[4].pow2), *(psi_pow[4].pow1), f);
	top = 4;
}



//---------------------------------------------------------------
// determine div[top+1] under knowledge of div[0, ..., top]
// sqrYY = Y^4 = (X^3+AX+B)^2, inv = 2^-1.
//


void eco_prime::next_psi (PsiPowers* & psi_pow,
			  lidia_size_t & top,
			  const ff_polmod & f,
			  const ff_pol & sqrYY,
			  const ff_element& inv2)
{
	lidia_size_t  k;

	if (top & 1) {
		// top+1 is even
		top ++;
		k = top >> 1;

		multiply (*(psi_pow[top].pow1), *(psi_pow[k+2].pow1), *(psi_pow[k-1].pow2), f);
		multiply (*(psi_pow[top].pow2), *(psi_pow[k-2].pow1), *(psi_pow[k+1].pow2), f);

		subtract (*(psi_pow[top].pow1), *(psi_pow[top].pow1), *(psi_pow[top].pow2));

		multiply (*(psi_pow[top].pow1), *(psi_pow[k].pow1), *(psi_pow[top].pow1), f);
		multiply (*(psi_pow[top].pow1), *(psi_pow[top].pow1), inv2.mantissa());
	}
	else {
		top++;
		k = top >> 1;

		multiply(*(psi_pow[top].pow1), *(psi_pow[k+2].pow1), *(psi_pow[k].pow3), f);
		multiply(*(psi_pow[top].pow2), *(psi_pow[k-1].pow1), *(psi_pow[k+1].pow3), f);

		if (! (k & 1))
			multiply(*(psi_pow[top].pow1), *(psi_pow[top].pow1), sqrYY, f);
		else
			multiply(*(psi_pow[top].pow2), *(psi_pow[top].pow2), sqrYY, f);

		subtract (*(psi_pow[top].pow1), *(psi_pow[top].pow1), *(psi_pow[top].pow2));
	}
	square(*(psi_pow[top].pow2), *(psi_pow[top].pow1), f);
	multiply(*(psi_pow[top].pow3), *(psi_pow[top].pow2), *(psi_pow[top].pow1), f);
}



//------------------------------------------------------------------
// the following function frees memory of psi_pow[b, ..., e], e >= b.


void eco_prime::free_psi (PsiPowers* & psi_pow, lidia_size_t b, lidia_size_t e)
{
	lidia_size_t  i;

	for (i = b; i <= e; i++) {
		//      if (psi_pow[i].pow1 != NULL)
		delete psi_pow[i].pow1;
		//      if (psi_pow[i].pow2 != NULL)
		delete psi_pow[i].pow2;
		//      if (psi_pow[i].pow3 != NULL)
		delete psi_pow[i].pow3;
	}
}



//
//
//  Compute one specific division polynomial
//
//

//  several functions for computation of one specific division polynomial:
//  build optimal plan, print plan, compute the one specific div-pol.

//-------------------------------------------------------------------
// to_use[i][0] || to_use[i][1] || to_use[i][2] == 1 iff i-th
// division polynomial is used for computation of k-th divpol
// Note: computation is recursive !!


void eco_prime::build_plan (lidia_size_t ** & to_use, lidia_size_t k)
{
	lidia_size_t n;
	to_use[k][0] = 1;

	while (k > 4) {
		n = k >> 1;
		if (k & 1) {
			to_use[n+2][0] ++;
			to_use[ n ][2] ++;
			to_use[n+1][2] ++;
			to_use[n-1][0] ++;
		}
		else {
			to_use[ n ][0] ++;
			to_use[n+2][0] ++;
			to_use[n-1][1] ++;
			to_use[n-2][0] ++;
			to_use[n+1][1] ++;
		}
		do {
			k --;
		} while (to_use[k][0] == 0 && to_use[k][1] == 0 && to_use[k][2] == 0);
	}
}



//----------------------------------------------------------------------
// output the indices of all division polynomials necessary for
// computation of k-th divpol.

void print_plan (lidia_size_t **to_use, lidia_size_t k)
{
	lidia_size_t i;

	for (i = 1; i <= k; i++) {
		std::cout << "{ (" << i << ") : [" << to_use[i][0] << ", " << to_use[i][1];
		std::cout << ", " << to_use[i][2] << "]";

		if (i % 3 == 0)
			std::cout << std::endl;
	}
	std::cout << "\n" << std::endl;
}



//----------------------------------------------------------------------
// compute the k-th divpol modulo f.

void eco_prime::compute_psi (ff_pol &res, lidia_size_t k,
                             const ff_polmod & f)
{
	lidia_size_t  **to_use;
	lidia_size_t i, pos, n, j;
	lidia_size_t size;

	ff_pol  ***psi;
	ff_pol  sqrY;
	ff_element qqq;
	ff_element ttt;
	ff_element inv2;
	ff_element tmp, tmp2;

	ttt = ff_element(0);
	sqrY.set_modulus(p);
	size = comparator< lidia_size_t >::max(k, 4);

	inv2 = ff_element(2);
	invert (inv2, inv2);

	to_use = new lidia_size_t*[size+1];
	memory_handler (to_use, "eco_prime::compute_psi", \
			"Allocating to_use");

	psi = new ff_pol**[size+1];
	memory_handler (psi, "eco_prime::compute_psi", \
			"Allocating psi");

	for (i = 0; i <= size; i ++) {
		to_use[i] = new lidia_size_t[3];
		psi[i] = new ff_pol*[3];
	}

	for (i = 0; i <= size; i++)
		for (j = 0; j < 3; j++)
			to_use[i][j] = 0;

	build_plan (to_use, k);

	CurveEqn (sqrY, p, A, B, f);
	square (sqrY, sqrY, f);

	psi[0][0] = new ff_pol; to_use[0][0] ++;
	psi[0][1] = new ff_pol; to_use[0][1] ++;
	psi[0][2] = new ff_pol; to_use[0][2] ++;

	psi[0][0]->set_modulus (p);
	psi[0][1]->set_modulus (p);
	psi[0][2]->set_modulus (p);

	// psi_1, psi_1^2, psi_1^3 

	psi[1][0] = new ff_pol;
	psi[1][0]->set_modulus (p);
	psi[1][0]->set_coefficient (0);

	if (to_use[1][1] > 0) {
		psi[1][1] = new ff_pol;
		psi[1][1]->set_modulus (p);
		psi[1][1]->set_coefficient (0);
	}

	if(to_use[1][2] > 0) {
		psi[1][2] = new ff_pol;
		psi[1][2]->set_modulus (p);
		psi[1][2]->set_coefficient (0);
	}

	to_use[1][0]++;

	// psi_2, psi_2^2, psi_2^3 

	psi[2][0] = new ff_pol;
	psi[2][0]->set_modulus(p);

	ttt = ff_element(2);
	psi[2][0]->set_coefficient (ttt.mantissa(), 0);

	if (to_use[2][1] > 0) {
		psi[2][1] = new ff_pol;
		psi[2][1]->set_modulus(p);
		ttt = ff_element(4);
		psi[2][1]->set_coefficient (ttt.mantissa(), 0);
	}

	if (to_use[2][2] > 0)   // ?? wo sind die y^i terme ??
	{
		psi[2][2] = new ff_pol;
		psi[2][2]->set_modulus(p);
		ttt = ff_element(8);
		psi[2][2]->set_coefficient(ttt.mantissa(), 0);
	}

	to_use[2][0]++;


	// psi_3, psi_3^2, psi_3^3 

	psi[3][0] = new ff_pol;
	psi[3][0]->set_modulus(p);
	ttt = ff_element(3);
	psi[3][0]->set_coefficient(ttt.mantissa(), 4);
	ttt = ff_element(0);
	psi[3][0]->set_coefficient(ttt.mantissa(), 3);
	multiply(ttt, ff_element(6), A);
	psi[3][0]->set_coefficient(ttt.mantissa(), 2);
	multiply (ttt, ff_element(12), B);
	psi[3][0]->set_coefficient(ttt.mantissa(), 1);
	square(ttt, A);
	negate (ttt, ttt);
	psi[3][0]->set_coefficient(ttt.mantissa(), 0);
	remainder (*psi[3][0], *psi[3][0], f);

	if (to_use[3][1] > 0 || to_use[3][2] > 0) {
		psi[3][1] = new ff_pol;
		square (*psi[3][1], *psi[3][0], f);
	}

	if (to_use[3][2] > 0) {
		psi[3][2] = new ff_pol;
		multiply(*psi[3][2], *psi[3][1], *psi[3][0], f);
		if (to_use[3][1] == 0)
			delete (psi[3][1]);
	}

	to_use[3][0]++;

	// psi_4, psi_4^2, psi_4^3

	psi[4][0] = new ff_pol;
	psi[4][0]->set_modulus(p);
	ttt = ff_element(4);
	psi[4][0]->set_coefficient(ttt.mantissa(), 6);
	ttt = ff_element(0);
	psi[4][0]->set_coefficient(ttt.mantissa(), 5);
	multiply (ttt, ff_element(20), A);
	psi[4][0]->set_coefficient(ttt.mantissa(), 4);
	multiply (ttt, ff_element(80), B);
	psi[4][0]->set_coefficient(ttt.mantissa(), 3);
	multiply (tmp, A, B);
	multiply (ttt, tmp, ff_element(16));
	negate (ttt, ttt);
	psi[4][0]->set_coefficient(ttt.mantissa(), 1);

	square (tmp, A);
	multiply(ttt, tmp, ff_element(20));
	negate (ttt, ttt);
	psi[4][0]->set_coefficient(ttt.mantissa(), 2);

	multiply (tmp, tmp, A);
	multiply (tmp, tmp, ff_element(4)); // tmp = 4.A^3
	square (tmp2, B);
	multiply (tmp2, ff_element(32), tmp2);
	add (ttt, tmp, tmp2);
	negate(ttt, ttt);
	psi[4][0]->set_coefficient(ttt.mantissa(), 0);

	remainder (*psi[4][0], *psi[4][0], f);

	if (to_use[4][1] > 0 || to_use[4][2] > 0) {
		psi[4][1] = new ff_pol;
		square (*psi[4][1], *psi[4][0], f);
	}

	if (to_use[4][2] > 0) {
		psi[4][2] = new ff_pol;
		multiply (*psi[4][2], *psi[4][1], *psi[4][0], f);
		if (to_use[4][1] == 0)
			delete (psi[4][1]);
	}
	to_use[4][0]++;

	//* * * * *  Now we start the recursion  * * * * * * * * * * *

	if (k >= 1 && k <= 4)
		res = *psi[k][0];
	else {
		pos = 4;
		while (k > pos) {
			pos ++;
			while (to_use[pos][0] + to_use[pos][1] + to_use[pos][2] == 0)
				pos ++;

			psi[pos][0] = new ff_pol;
			psi[pos][0]->set_modulus(p);

			// Computation of psi[pos]
			if (pos & 1) {
				n = pos >> 1;
				multiply (*psi[0][0], *psi[n+2][0], *psi[ n ][2], f);
				multiply (*psi[0][1], *psi[n+1][2], *psi[n-1][0], f);

				if (n & 1) {
					multiply (*psi[0][1], *psi[0][1], sqrY, f);
				}
				else {
					multiply (*psi[0][0], *psi[0][0], sqrY, f);
				}

				subtract (*psi[pos][0], *psi[0][0], *psi[0][1]);
				to_use[n+2][0] --;
				to_use[ n ][2] --;
				to_use[n+1][2] --;
				to_use[n-1][0] --;

				if (!to_use[n+2][0])
					delete (psi[n+2][0]);
				if (!to_use[ n ][2])
					delete (psi[ n ][2]);
				if (!to_use[n+1][2])
					delete (psi[n+1][2]);
				if (!to_use[n-1][0])
					delete (psi[n-1][0]);
			}
			else {
				n = pos >> 1;

				multiply (*psi[0][0], *psi[n+2][0], *psi[n-1][1], f);
				multiply (*psi[0][1], *psi[n-2][0], *psi[n+1][1], f);
				subtract  (*psi[pos][0], *psi[0][0], *psi[0][1]);
				multiply  (*psi[pos][0], *psi[pos][0], *psi[n][0], f);

				qqq = inv2;
				multiply (*psi[pos][0], *psi[pos][0], qqq.mantissa());
				to_use[n+2][0] --;
				to_use[n-1][1] --;
				to_use[n-2][0] --;
				to_use[n+1][1] --;
				to_use[ n ][0] --;

				if (!to_use[n+2][0])
					delete (psi[n+2][0]);
				if (!to_use[n-1][1])
					delete (psi[n-1][1]);
				if (!to_use[n-2][0])
					delete (psi[n-2][0]);
				if (!to_use[n+1][1])
					delete (psi[n+1][1]);
				if (!to_use[ n ][0])
					delete (psi[ n ][0]);
			}

			if (to_use[pos][1] > 0 || to_use[pos][2] > 0) {
				psi[pos][1] = new ff_pol;
				square (*psi[pos][1], *psi[pos][0], f);
			}

			if (to_use[pos][2] > 0) {
				psi[pos][2] = new ff_pol;
				multiply (*psi[pos][2], *psi[pos][1], *psi[pos][0], f);

				if (to_use[pos][1] == 0)
					delete (psi[pos][1]);
			}

			if (to_use[pos][0] == 0)
				delete (psi[pos][0]);
		} // end of while

		res = *psi[k][0];
	} // end of else


	for (i = 0; i <= size; i++)
		for (j = 0; j < 3; j++) {
			if (to_use[i][j] != 0)
				delete psi[i][j];
		}

	for (i = 0; i <= size; i++) {
		delete[] psi[i];
		delete[] to_use[i];
	}
	delete[] psi;
	delete[] to_use;
}



//----------------------------------------------------------------------
// compute the k-th divpol.

void eco_prime::compute_psi (ff_pol & res, lidia_size_t k)
{
	lidia_size_t  **to_use;
	lidia_size_t i, pos, n, j;
	lidia_size_t size;

	ff_pol  ***psi;
	ff_pol  sqrY;
	ff_element qqq;
	ff_element ttt;
	ff_element inv2;
	ff_element tmp, tmp2;

	ttt = ff_element(0);
	sqrY.set_modulus(p);
	size = comparator< lidia_size_t >::max(k, 4);

	inv2 = ff_element(2);
	invert (inv2, inv2);

	to_use = new lidia_size_t*[size+1];
	memory_handler (to_use, "eco_prime::compute_psi", \
			"Allocating to_use");

	psi = new ff_pol**[size+1];
	memory_handler (psi, "eco_prime::compute_psi", \
			"Allocating psi");

	for (i = 0; i <= size; i ++) {
		to_use[i] = new lidia_size_t[3];
		psi[i] = new ff_pol*[3];
	}
	for (i = 0; i <= size; i++)
		for (j = 0; j < 3; j++)
			to_use[i][j] = 0;

	build_plan (to_use, k);

	ff_pol f;
	f.set_modulus(p);
	f.set_coefficient(4);

	CurveEqn (sqrY, p, A, B, f);
	square (sqrY, sqrY);

	psi[0][0] = new ff_pol; to_use[0][0] ++;
	psi[0][1] = new ff_pol; to_use[0][1] ++;
	psi[0][2] = new ff_pol; to_use[0][2] ++;

	psi[0][0]->set_modulus (p);
	psi[0][1]->set_modulus (p);
	psi[0][2]->set_modulus (p);

	// psi_1, psi_1^2, psi_1^3 

	psi[1][0] = new ff_pol;
	psi[1][0]->set_modulus (p);
	psi[1][0]->set_coefficient (0);

	if (to_use[1][1] > 0) {
		psi[1][1] = new ff_pol;
		psi[1][1]->set_modulus (p);
		psi[1][1]->set_coefficient (0);
	}

	if(to_use[1][2] > 0) {
		psi[1][2] = new ff_pol;
		psi[1][2]->set_modulus (p);
		psi[1][2]->set_coefficient (0);
	}

	to_use[1][0]++;

	// psi_2, psi_2^2, psi_2^3 

	psi[2][0] = new ff_pol;
	psi[2][0]->set_modulus(p);

	ttt = ff_element(2);
	psi[2][0]->set_coefficient (ttt.mantissa(), 0);

	if (to_use[2][1] > 0) {
		psi[2][1] = new ff_pol;
		psi[2][1]->set_modulus(p);
		ttt = ff_element(4);
		psi[2][1]->set_coefficient (ttt.mantissa(), 0);
	}

	if (to_use[2][2] > 0) {
		// ?? wo sind die y^i terme ??
		psi[2][2] = new ff_pol;
		psi[2][2]->set_modulus(p);
		ttt = ff_element(8);
		psi[2][2]->set_coefficient(ttt.mantissa(), 0);
	}

	to_use[2][0]++;


	// psi_3, psi_3^2, psi_3^3 

	psi[3][0] = new ff_pol;
	psi[3][0]->set_modulus(p);
	ttt = ff_element(3);
	psi[3][0]->set_coefficient(ttt.mantissa(), 4);
	ttt = ff_element(0);
	psi[3][0]->set_coefficient(ttt.mantissa(), 3);
	multiply(ttt, ff_element(6), A);
	psi[3][0]->set_coefficient(ttt.mantissa(), 2);
	multiply (ttt, ff_element(12), B);
	psi[3][0]->set_coefficient(ttt.mantissa(), 1);
	square(ttt, A);
	negate (ttt, ttt);
	psi[3][0]->set_coefficient(ttt.mantissa(), 0);

	if (to_use[3][1] > 0 || to_use[3][2] > 0) {
		psi[3][1] = new ff_pol;
		square (*psi[3][1], *psi[3][0]);
	}

	if (to_use[3][2] > 0) {
		psi[3][2] = new ff_pol;
		multiply(*psi[3][2], *psi[3][1], *psi[3][0]);
		if (to_use[3][1] == 0)
			delete (psi[3][1]);
	}

	to_use[3][0]++;

	// psi_4, psi_4^2, psi_4^3

	psi[4][0] = new ff_pol;
	psi[4][0]->set_modulus(p);
	ttt = ff_element(4);
	psi[4][0]->set_coefficient(ttt.mantissa(), 6);
	ttt = ff_element(0);
	psi[4][0]->set_coefficient(ttt.mantissa(), 5);
	multiply (ttt, ff_element(20), A);
	psi[4][0]->set_coefficient(ttt.mantissa(), 4);
	multiply (ttt, ff_element(80), B);
	psi[4][0]->set_coefficient(ttt.mantissa(), 3);
	multiply (tmp, A, B);
	multiply (ttt, tmp, ff_element(16));
	negate (ttt, ttt);
	psi[4][0]->set_coefficient(ttt.mantissa(), 1);

	square (tmp, A);
	multiply(ttt, tmp, ff_element(20));
	negate (ttt, ttt);
	psi[4][0]->set_coefficient(ttt.mantissa(), 2);

	multiply (tmp, tmp, A);
	multiply (tmp, tmp, ff_element(4)); // tmp = 4.A^3
	square (tmp2, B);
	multiply (tmp2, ff_element(32), tmp2);
	add (ttt, tmp, tmp2);
	negate(ttt, ttt);
	psi[4][0]->set_coefficient(ttt.mantissa(), 0);

	if (to_use[4][1] > 0 || to_use[4][2] > 0) {
		psi[4][1] = new ff_pol;
		square (*psi[4][1], *psi[4][0]);
	}

	if (to_use[4][2] > 0) {
		psi[4][2] = new ff_pol;
		multiply (*psi[4][2], *psi[4][1], *psi[4][0]);
		if (to_use[4][1] == 0)
			delete (psi[4][1]);
	}
	to_use[4][0]++;

	//* * * * *  Now we start the recursion  * * * * * * * * * * *

	if (k >= 1 && k <= 4)
		res = *psi[k][0];
	else {
		pos = 4;
		while (k > pos) {
			pos ++;
			while (to_use[pos][0] + to_use[pos][1] + to_use[pos][2] == 0)
				pos ++;

			psi[pos][0] = new ff_pol;
			psi[pos][0]->set_modulus(p);

			// Computation of psi[pos]
			if (pos & 1) {
				n = pos >> 1;
				multiply (*psi[0][0], *psi[n+2][0], *psi[ n ][2]);
				multiply (*psi[0][1], *psi[n+1][2], *psi[n-1][0]);

				if (n & 1) {
					multiply (*psi[0][1], *psi[0][1], sqrY);
				}
				else {
					multiply (*psi[0][0], *psi[0][0], sqrY);
				}

				subtract (*psi[pos][0], *psi[0][0], *psi[0][1]);
				to_use[n+2][0] --;
				to_use[ n ][2] --;
				to_use[n+1][2] --;
				to_use[n-1][0] --;

				if (!to_use[n+2][0])
					delete (psi[n+2][0]);
				if (!to_use[ n ][2])
					delete (psi[ n ][2]);
				if (!to_use[n+1][2])
					delete (psi[n+1][2]);
				if (!to_use[n-1][0])
					delete (psi[n-1][0]);
			}
			else {
				n = pos >> 1;

				multiply (*psi[0][0], *psi[n+2][0], *psi[n-1][1]);
				multiply (*psi[0][1], *psi[n-2][0], *psi[n+1][1]);
				subtract  (*psi[pos][0], *psi[0][0], *psi[0][1]);
				multiply  (*psi[pos][0], *psi[pos][0], *psi[n][0]);

				qqq = inv2;
				multiply (*psi[pos][0], *psi[pos][0], qqq.mantissa());
				to_use[n+2][0] --;
				to_use[n-1][1] --;
				to_use[n-2][0] --;
				to_use[n+1][1] --;
				to_use[ n ][0] --;

				if (!to_use[n+2][0])
					delete (psi[n+2][0]);
				if (!to_use[n-1][1])
					delete (psi[n-1][1]);
				if (!to_use[n-2][0])
					delete (psi[n-2][0]);
				if (!to_use[n+1][1])
					delete (psi[n+1][1]);
				if (!to_use[ n ][0])
					delete (psi[ n ][0]);
			}

			if (to_use[pos][1] > 0 || to_use[pos][2] > 0) {
				psi[pos][1] = new ff_pol;
				square (*psi[pos][1], *psi[pos][0]);
			}

			if (to_use[pos][2] > 0) {
				psi[pos][2] = new ff_pol;
				multiply (*psi[pos][2], *psi[pos][1], *psi[pos][0]);

				if (to_use[pos][1] == 0)
					delete (psi[pos][1]);
			}

			if (to_use[pos][0] == 0)
				delete (psi[pos][0]);
		} // end of while

		res = *psi[k][0];
	} // end of else

	for (i = 0; i <= size; i++)
		for (j = 0; j < 3; j++) {
			if (to_use[i][j] != 0)
				delete psi[i][j];
		}

	for (i = 0; i <= size; i++) {
		delete[] psi[i];
		delete[] to_use[i];
	}
	delete[] psi;
	delete[] to_use;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
