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



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"



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
// polynomials psi[0], ..., psi[nmax+3] modulo the modulus of F
// top is set to maximally already computed division polynomial
// triple.

void eco_gf2n::
init_psi_comp (PsiPowers* & psi_pow,
	       lidia_size_t & top,
	       lidia_size_t nmax,
	       ff_polmod & F)
{
	lidia_size_t i;

	psi_pow = new PsiPowers[nmax+4];
	memory_handler (psi, "eco_gf2n::init_psi_comp()",
			"Allocating psi_pow");

	if (nmax >= 2)
		for (i = 5; i <= nmax+3; i++) {
			(psi_pow[i]).pow1 = new ff_pol;
			(psi_pow[i]).pow2 = new ff_pol;
			(psi_pow[i]).pow3 = new ff_pol;
		}

	// psi_0, psi_0^2, psi_0^3

	(psi_pow[0]).pow1 = new ff_pol;
	(psi_pow[0]).pow2 = new ff_pol;
	(psi_pow[0]).pow3 = new ff_pol;
	(*(psi_pow[0].pow1)).assign_zero();
	(*(psi_pow[0].pow2)).assign_zero();
	(*(psi_pow[0].pow3)).assign_zero();

	// psi_1, psi_1^2, psi_1^3

	(psi_pow[1]).pow1 = new ff_pol;
	(psi_pow[1]).pow2 = new ff_pol;
	(psi_pow[1]).pow3 = new ff_pol;
	(*(psi_pow[1].pow1)).assign_one();
	(*(psi_pow[1].pow2)).assign_one();
	(*(psi_pow[1].pow3)).assign_one();

	// psi_2, psi_2^2, psi_2^3 

	psi_pow[2].pow1 = new ff_pol;
	psi_pow[2].pow2 = new ff_pol;
	psi_pow[2].pow3 = new ff_pol;
	(*(psi_pow[2].pow1)).assign_x();

	(*(psi_pow[2].pow2)).assign(ff_pol(2));
	remainder(*(psi_pow[2].pow2), *(psi_pow[2].pow2), F.modulus());

	(*(psi_pow[2].pow3)).assign(ff_pol(3));
	remainder(*(psi_pow[2].pow3), *(psi_pow[2].pow3), F.modulus());

	// psi_3, psi_3^2, psi_3^3

	psi_pow[3].pow1 = new ff_pol;
	psi_pow[3].pow2 = new ff_pol;
	psi_pow[3].pow3 = new ff_pol;

	psi_pow[3].pow1->assign(ff_pol(4));
	psi_pow[3].pow1->set_coefficient (gf2n(1), 3);
	psi_pow[3].pow1->set_coefficient (A6, 0);
	remainder(*(psi_pow[3].pow1), *(psi_pow[3].pow1), F.modulus());

	square (*(psi_pow[3].pow2), *(psi_pow[3].pow1), F);
	multiply (*(psi_pow[3].pow3), *(psi_pow[3].pow2), *(psi_pow[3].pow1), F);

	// psi_4, psi_4^2, psi_4^3

	psi_pow[4].pow1 = new ff_pol;
	psi_pow[4].pow2 = new ff_pol;
	psi_pow[4].pow3 = new ff_pol;
	psi_pow[4].pow1->assign(ff_pol(6));
	psi_pow[4].pow1->set_coefficient (A6, 2);
	remainder(*(psi_pow[4].pow1), *(psi_pow[4].pow1), F.modulus());

	square (*(psi_pow[4].pow2), *(psi_pow[4].pow1), F);
	multiply (*(psi_pow[4].pow3), *(psi_pow[4].pow2), *(psi_pow[4].pow1), F);

	top = 4;
}



//---------------------------------------------------------------
// determine and stores div[top+1] under knowledge of div[0, ..., top]
//

void eco_gf2n::next_psi (PsiPowers* & psi_pow,
			 lidia_size_t & top,
			 ff_polmod & F)
{
	lidia_size_t  k;

	if (top & 1)     // (top + 1) is even
	{
		top++;
		k = top >> 1;

		multiply (*(psi_pow[top].pow1), *(psi_pow[k+2].pow1),
			  *(psi_pow[k-1].pow2), F);
		multiply (*(psi_pow[top].pow2), *(psi_pow[k-2].pow1),
			  *(psi_pow[k+1].pow2), F);
		add (*(psi_pow[top].pow1), *(psi_pow[top].pow1), *(psi_pow[top].pow2));

		multiply (*(psi_pow[top].pow1), *(psi_pow[k].pow1),
			  *(psi_pow[top].pow1), F);
		shift_right(*(psi_pow[top].pow1), *(psi_pow[top].pow1)
			    + ((psi_pow[top].pow1)->const_term() /
			       F.modulus().const_term())
			    * F.modulus());
	}
	else {                        // (top + 1) odd
		top++;
		k = top >> 1;

		multiply(*(psi_pow[top].pow1), *(psi_pow[k+2].pow1),
			 *(psi_pow[k].pow3), F);
		multiply(*(psi_pow[top].pow2), *(psi_pow[k-1].pow1),
			 *(psi_pow[k+1].pow3), F);
		add (*(psi_pow[top].pow1), *(psi_pow[top].pow1), *(psi_pow[top].pow2));
	}
	square(*(psi_pow[top].pow2), *(psi_pow[top].pow1), F);
	multiply(*(psi_pow[top].pow3), *(psi_pow[top].pow2),
		 *(psi_pow[top].pow1), F);
}



//------------------------------------------------------------------
// the following function frees memory of psi_pow[b, ..., e], e >= b.

void eco_gf2n::free_psi (PsiPowers* & psi_pow, lidia_size_t b, lidia_size_t e)
{
	lidia_size_t  i;

	for (i = b; i <= e; i++) {
		delete psi_pow[i].pow1;
		delete psi_pow[i].pow2;
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

void eco_gf2n::build_plan (lidia_size_t ** & to_use, lidia_size_t k)
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

void eco_gf2n::print_plan (lidia_size_t **to_use, lidia_size_t k)
{
	lidia_size_t i;

	for (i = 1; i <= k; i++) {
		std::cout << "{ (" << i << ") : [" << to_use[i][0] << ", " << to_use[i][1];
		std::cout << ", " << to_use[i][2] << "]";

		if (i % 3 == 0) std::cout << "\n";
	}
	std::cout << "\n\n";
}



//----------------------------------------------------------------------
// compute the k-th divpol mod F.mod

void eco_gf2n::compute_psi (ff_pol & res, lidia_size_t k, ff_polmod & F)
{
	lidia_size_t  **to_use;
	lidia_size_t i, pos, n, j;
	lidia_size_t size;
	ff_pol  ***psi;

	size = comparator< lidia_size_t >::max(k, 4);

	to_use = new lidia_size_t*[size+1];
	memory_handler (to_use, "eco_gf2n::compute_psi", "Allocating to_use");

	psi = new ff_pol**[size+1];
	memory_handler (psi, "eco_gf2n::compute_psi", "Allocating psi");

	for (i = 0; i <= size; i ++) {
		to_use[i] = new lidia_size_t[3];
		psi[i] = new ff_pol*[3];
	}

	for (i = 0; i <= size; i++)
		for (j = 0; j < 3; j++)
			to_use[i][j] = 0;

	build_plan (to_use, k);

	psi[0][0] = new ff_pol; to_use[0][0] ++;
	psi[0][1] = new ff_pol; to_use[0][1] ++;
	psi[0][2] = new ff_pol; to_use[0][2] ++;

	// psi_1, psi_1^2, psi_1^3 

	psi[1][0] = new ff_pol;
	psi[1][0]->assign_one();

	if (to_use[1][1] > 0) {
		psi[1][1] = new ff_pol;
		psi[1][1]->assign_one();
	}

	if (to_use[1][2] > 0) {
		psi[1][2] = new ff_pol;
		psi[1][2]->assign_one();
	}

	to_use[1][0]++;

	// psi_2, psi_2^2, psi_2^3 

	psi[2][0] = new ff_pol;
	psi[2][0]->assign_x();

	if (to_use[2][1] > 0 || to_use[2][2] > 0) {
		psi[2][1] = new ff_pol;
		square(*psi[2][1], *psi[2][0], F);
	}

	if (to_use[2][2] > 0) {
		psi[2][2] = new ff_pol;
		multiply(*psi[2][2], *psi[2][1], *psi[2][0], F);
		if (to_use[2][1] == 0)
			delete psi[2][1];
	}

	to_use[2][0]++;

	// psi_3, psi_3^2, psi_3^3 

	psi[3][0] = new ff_pol;
	psi[3][0] -> assign(ff_pol(4));
	psi[3][0] -> set_coefficient(gf2n(1), 3);
	psi[3][0] -> set_coefficient(A6, 0);
	remainder(*psi[3][0], *psi[3][0], F.modulus());

	if (to_use[3][1] > 0 || to_use[3][2] > 0) {
		psi[3][1] = new ff_pol;
		square (*psi[3][1], *psi[3][0], F);
	}

	if (to_use[3][2] > 0) {
		psi[3][2] = new ff_pol;
		multiply(*psi[3][2], *psi[3][1], *psi[3][0], F);
		if (to_use[3][1] == 0)
			delete (psi[3][1]);
	}

	to_use[3][0]++;

	// psi_4, psi_4^2, psi_4^3

	psi[4][0] = new ff_pol;
	psi[4][0] -> assign(ff_pol(6));
	psi[4][0] -> set_coefficient(A6, 2);
	remainder(*psi[4][0], *psi[4][0], F.modulus());

	if (to_use[4][1] > 0 || to_use[4][2] > 0) {
		psi[4][1] = new ff_pol;
		square (*psi[4][1], *psi[4][0], F);
	}

	if (to_use[4][2] > 0) {
		psi[4][2] = new ff_pol;
		multiply (*psi[4][2], *psi[4][1], *psi[4][0], F);
		if (to_use[4][1] == 0)
			delete (psi[4][1]);
	}
	to_use[4][0]++;

	// * * * * *  Now we start the recursion  * * * * * * * * * * *

	if (k >= 1 && k <= 4)
		res = *psi[k][0];
	else {
		pos = 4;
		while (k > pos) {
			pos ++;
			while (to_use[pos][0] + to_use[pos][1] + to_use[pos][2] == 0)
				pos ++;

			psi[pos][0] = new ff_pol;

			// Computation of psi[pos]
			if (pos & 1) {
				n = pos >> 1;
				multiply (*psi[0][0], *psi[n+2][0], *psi[ n ][2], F);
				multiply (*psi[0][1], *psi[n+1][2], *psi[n-1][0], F);

				add (*psi[pos][0], *psi[0][0], *psi[0][1]);
				to_use[n+2][0] --;
				to_use[ n ][2] --;
				to_use[n+1][2] --;
				to_use[n-1][0] --;

				if (!to_use[n+2][0]) delete (psi[n+2][0]);
				if (!to_use[ n ][2]) delete (psi[ n ][2]);
				if (!to_use[n+1][2]) delete (psi[n+1][2]);
				if (!to_use[n-1][0]) delete (psi[n-1][0]);
			}
			else {
				n = pos >> 1;

				multiply (*psi[0][0], *psi[n+2][0], *psi[n-1][1], F);
				multiply (*psi[0][1], *psi[n-2][0], *psi[n+1][1], F);
				add  (*psi[pos][0], *psi[0][0], *psi[0][1]);
				multiply  (*psi[pos][0], *psi[pos][0], *psi[n][0], F);
				shift_right(*(psi[pos][0]), *(psi[pos][0])
					    + (psi[pos][0]->const_term() /
					       F.modulus().const_term())
					    * F.modulus());

				to_use[n+2][0] --;
				to_use[n-1][1] --;
				to_use[n-2][0] --;
				to_use[n+1][1] --;
				to_use[ n ][0] --;

				if (!to_use[n+2][0]) delete (psi[n+2][0]);
				if (!to_use[n-1][1]) delete (psi[n-1][1]);
				if (!to_use[n-2][0]) delete (psi[n-2][0]);
				if (!to_use[n+1][1]) delete (psi[n+1][1]);
				if (!to_use[ n ][0]) delete (psi[ n ][0]);
			}

			if (to_use[pos][1] > 0 || to_use[pos][2] > 0) {
				psi[pos][1] = new ff_pol;
				square (*psi[pos][1], *psi[pos][0], F);
			}

			if (to_use[pos][2] > 0) {
				psi[pos][2] = new ff_pol;
				multiply (*psi[pos][2], *psi[pos][1], *psi[pos][0], F);

				if (to_use[pos][1] == 0)
					delete (psi[pos][1]);
			}

			if (to_use[pos][0] == 0) delete (psi[pos][0]);
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

void eco_gf2n::compute_psi (ff_pol & res, lidia_size_t k)
{
	lidia_size_t  **to_use;
	lidia_size_t i, pos, n, j;
	lidia_size_t size;
	ff_pol  ***psi;

	size = comparator< lidia_size_t >::max(k, 4);

	to_use = new lidia_size_t*[size+1];
	memory_handler (to_use, "eco_gf2n::compute_psi", "Allocating to_use");

	psi = new ff_pol**[size+1];
	memory_handler (psi, "eco_gf2n::compute_psi", "Allocating psi");

	for (i = 0; i <= size; i ++) {
		to_use[i] = new lidia_size_t[3];
		psi[i] = new ff_pol*[3];
	}

	for (i = 0; i <= size; i++)
		for (j = 0; j < 3; j++)
			to_use[i][j] = 0;

	build_plan (to_use, k);

	psi[0][0] = new ff_pol; to_use[0][0] ++;
	psi[0][1] = new ff_pol; to_use[0][1] ++;
	psi[0][2] = new ff_pol; to_use[0][2] ++;

	// psi_1, psi_1^2, psi_1^3 

	psi[1][0] = new ff_pol;
	psi[1][0]->assign(ff_pol(0));

	if (to_use[1][1] > 0) {
		psi[1][1] = new ff_pol;
		psi[1][1]->assign(ff_pol(0));
	}

	if(to_use[1][2] > 0) {
		psi[1][2] = new ff_pol;
		psi[1][2]->assign(ff_pol(0));
	}

	to_use[1][0]++;

	// psi_2, psi_2^2, psi_2^3 

	psi[2][0] = new ff_pol;
	psi[2][0]->assign(ff_pol(1));

	if (to_use[2][1] > 0 || to_use[2][2] > 0) {
		psi[2][1] = new ff_pol;
		square(*psi[2][1], *psi[2][0]);
	}

	if (to_use[2][2] > 0) {
		psi[2][2] = new ff_pol;
		multiply(*psi[2][2], *psi[2][1], *psi[2][0]);
		if (to_use[2][1] == 0)
			delete psi[2][1];
	}

	to_use[2][0]++;


	// psi_3, psi_3^2, psi_3^3 

	psi[3][0] = new ff_pol;
	psi[3][0] -> assign(ff_pol(4));
	psi[3][0] -> set_coefficient(gf2n(1), 3);
	psi[3][0] -> set_coefficient(A6, 0);

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
	psi[4][0] -> assign(ff_pol(6));
	psi[4][0] -> set_coefficient(A6, 2);

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

	// * * * *  Now we start the recursion  * * * * * * * * * *

	if (k >= 1 && k <= 4)
		res = *psi[k][0];
	else {
		pos = 4;
		while (k > pos) {
			pos ++;
			while (to_use[pos][0] + to_use[pos][1] + to_use[pos][2] == 0)
				pos ++;

			psi[pos][0] = new ff_pol;

			// Computation of psi[pos]
			if (pos & 1) {
				n = pos >> 1;
				multiply (*psi[0][0], *psi[n+2][0], *psi[ n ][2]);
				multiply (*psi[0][1], *psi[n+1][2], *psi[n-1][0]);

				add (*psi[pos][0], *psi[0][0], *psi[0][1]);
				to_use[n+2][0] --;
				to_use[ n ][2] --;
				to_use[n+1][2] --;
				to_use[n-1][0] --;

				if (!to_use[n+2][0]) delete (psi[n+2][0]);
				if (!to_use[ n ][2]) delete (psi[ n ][2]);
				if (!to_use[n+1][2]) delete (psi[n+1][2]);
				if (!to_use[n-1][0]) delete (psi[n-1][0]);
			}
			else {
				n = pos >> 1;

				multiply (*psi[0][0], *psi[n+2][0], *psi[n-1][1]);
				multiply (*psi[0][1], *psi[n-2][0], *psi[n+1][1]);
				add  (*psi[pos][0], *psi[0][0], *psi[0][1]);
				multiply  (*psi[pos][0], *psi[pos][0], *psi[n][0]);
				shift_right(*psi[pos][0], *psi[pos][0]);

				to_use[n+2][0] --;
				to_use[n-1][1] --;
				to_use[n-2][0] --;
				to_use[n+1][1] --;
				to_use[ n ][0] --;

				if (!to_use[n+2][0]) delete (psi[n+2][0]);
				if (!to_use[n-1][1]) delete (psi[n-1][1]);
				if (!to_use[n-2][0]) delete (psi[n-2][0]);
				if (!to_use[n+1][1]) delete (psi[n+1][1]);
				if (!to_use[ n ][0]) delete (psi[ n ][0]);
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

			if (to_use[pos][0] == 0) delete (psi[pos][0]);
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
