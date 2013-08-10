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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint_lattice.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// using `only` bigints
//
//
// Benne de Weger
//
// * TrD_lll()
// * TrD_lll_trans(T)
// * TrD_lll_rand()
//
// Buchmann Kessler
//
// * Tr_lin_gen_system(T)
//

//
// Benne de Weger version of lll
// result : lll - reduced lattice for parameter a / b
//
void
bigint_lattice::TrD_lll (dense_alg< bigint > & da, lattice_info& li)
{
	debug_handler("bigint_lattice", "TrD_lll(da, li)");

	bigint tempbin0;
	bigint tempbin1;
	bigint tempbin2;
	bigint tempbin3;
	bigint tempbin4;
	bigint tempbin5;
	bigint *tempvectbin0;
	bigint *tempvectbin1;
	bigint *lzdel;
	bigint **lz;
	bigint **tempmatrixbin0;
	p_vector< bigint > vector;


//
// Allocating memory for lattice
//
	lz = new bigint*[2 * da.b.rows];
	memory_handler(lz, "bigint_lattice", "TrD_lll() :: "
		       "not enough memory !");
	lzdel = new bigint[da.b.rows * da.b.columns + da.b.rows * da.b.rows +
			  2 * da.b.rows + da.b.columns + 2];
	memory_handler(lzdel, "bigint_lattice", "TrD_lll() :: "
		       "not enough memory !");
	tempmatrixbin0 = &lz[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		lz[i] = &lzdel[i * da.b.rows];
		tempmatrixbin0[i] = &lzdel[da.b.rows * da.b.rows + i * da.b.columns];
	}

//
// Allocating memory for vectors
//
	tempvectbin0 = &lzdel[da.b.rows * da.b.rows + da.b.rows * da.b.columns];
	tempvectbin1 = &tempvectbin0[da.b.rows + 1];


	li.lll.reduction_steps = 0;
	li.lll.swaps = 0;
	li.lll.correction_steps = 0;
	vector.vectsize = da.b.rows + 1;
	vector.assign_zero(tempvectbin0);
	tempvectbin0[0].assign_one();
	vector.vectsize = da.b.columns;

//
// Compute a bigint approximation of the Gram - Schmidt - Matrix
// Compute my`s by giving a matrix and a vector
//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		vector.assign(tempmatrixbin0[i], da.s.value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			vector.scalprod(lz[j][i], da.s.value[i], tempmatrixbin0[j]);
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				LiDIA::multiply(tempmatrixbin0[i][cx], tempvectbin0[j + 1],
						   tempmatrixbin0[i][cx]);
				LiDIA::multiply(tempvectbin1[cx], lz[j][i],
						   tempmatrixbin0[j][cx]);
				LiDIA::subtract(tempmatrixbin0[i][cx], tempmatrixbin0[i][cx],
						   tempvectbin1[cx]);
			}
			for (lidia_size_t m = 0; m < da.b.columns; m++) {
				LiDIA::divide(tempmatrixbin0[i][m], tempmatrixbin0[i][m],
						 tempvectbin0[j]);
			}
		}
		vector.scalprod(tempvectbin0[i + 1], tempmatrixbin0[i],
				tempmatrixbin0[i]);
		LiDIA::divide(tempvectbin0[i + 1], tempvectbin0[i + 1], tempvectbin0[i]);
	}

	lidia_size_t l, k = 1;

//
// Begin of lll - Algoritm
//
	while (k != da.b.rows) {
		l = k - 1;
		tempbin0.assign(abs(lz[l][k]));
		tempbin0.multiply_by_2();
		if (tempbin0.compare(tempvectbin0[l + 1]) > 0) {
			li.lll.reduction_steps++;
			tempbin1.assign(abs(tempvectbin0[l + 1]));
			LiDIA::add(tempbin5, tempbin0, tempbin1);

			LiDIA::divide(tempbin5, tempbin5, tempbin1);
			tempbin5.divide_by_2();

			if (lz[l][k].sign() != tempvectbin0[l + 1].sign()) {
				tempbin5.negate();
			}
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				LiDIA::multiply(tempvectbin1[cx], tempbin5, da.s.value[l][cx]);
				LiDIA::subtract(da.s.value[k][cx], da.s.value[k][cx],
						   tempvectbin1[cx]);
			}

			for (lidia_size_t j = 0; j < l; j++) {
				LiDIA::multiply(tempbin3, tempbin5, lz[j][l]);
				LiDIA::subtract(lz[j][k], lz[j][k], tempbin3);
			}
			LiDIA::multiply(tempbin3, tempbin5, tempvectbin0[l + 1]);
			LiDIA::subtract(lz[l][k], lz[l][k], tempbin3);
		}

		LiDIA::square(tempbin0, tempvectbin0[k]);
		LiDIA::multiply(tempbin1, tempbin0, da.b.y_nom);
		LiDIA::square(tempbin0, lz[l][k]);
		LiDIA::multiply(tempbin2, tempbin0, da.b.y_denom);
		LiDIA::subtract(tempbin0, tempbin1, tempbin2);
		LiDIA::multiply(tempbin1, tempvectbin0[k - 1], tempvectbin0[k + 1]);
		LiDIA::multiply(tempbin2, tempbin1, da.b.y_denom);

		if (tempbin2.compare(tempbin0) < 0) {
			vector.swap(da.s.value[k - 1], da.s.value[k]);
			++li.lll.swaps;
			for (lidia_size_t j = 0; j <= k - 2; j++) {
				LiDIA::swap(lz[j][k - 1], lz[j][k]);
			}
			for (lidia_size_t i = k + 1; i < da.b.rows; i++) {
				LiDIA::multiply(tempbin5, lz[k - 1][i], lz[k - 1][k]);
				LiDIA::multiply(tempbin4, lz[k][i], tempvectbin0[k - 1]);
				LiDIA::add(tempbin3, tempbin5, tempbin4);
				LiDIA::divide(tempbin5, tempbin3, tempvectbin0[k]);
				LiDIA::multiply(tempbin4, lz[k - 1][i], tempvectbin0[k + 1]);
				LiDIA::multiply(tempbin3, lz[k][i], lz[k - 1][k]);
				LiDIA::subtract(tempbin1, tempbin4, tempbin3);
				LiDIA::divide(tempbin4, tempbin1, tempvectbin0[k]);
				lz[k - 1][i].assign(tempbin5);
				lz[k][i].assign(tempbin4);
			}
			LiDIA::square(tempbin3, lz[k - 1][k]);
			LiDIA::multiply(tempbin5, tempvectbin0[k - 1], tempvectbin0[k + 1]);
			LiDIA::add(tempbin4, tempbin5, tempbin3);
			LiDIA::divide(tempvectbin0[k], tempbin4, tempvectbin0[k]);

			if (k > 1) {
				k--;
			}
		}
		else {
//
// Reduction Step
//
			for (l = k - 2; l >= 0; l--) {
				tempbin3.assign(abs(lz[l][k]));
				tempbin3.multiply_by_2();
				if (tempbin3.compare(tempvectbin0[l + 1]) > 0) {
//                  si = (lz[l][k]).sign() * (tempvectbin0[l + 1]).sign();
// Was passiert, wenn eines der beiden Vorzeichen Null ist ?
// Kann dies auftreten ?
//
					li.lll.reduction_steps++;
					tempbin5.assign(abs(tempvectbin0[l + 1]));
					LiDIA::add(tempbin1, tempbin3, tempbin5);
					LiDIA::divide(tempbin4, tempbin1, tempbin5);
					tempbin4.divide_by_2();
					if ((lz[l][k]).sign() != (tempvectbin0[l + 1]).sign()) {
						tempbin4.negate();
					}

					for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
						LiDIA::multiply(tempvectbin1[cx], tempbin4,
								   da.s.value[l][cx]);
						LiDIA::subtract(da.s.value[k][cx], da.s.value[k][cx],
								   tempvectbin1[cx]);
					}
					for (lidia_size_t j = 0; j < l; j++) {
						LiDIA::multiply(tempbin0, tempbin4, lz[j][l]);
						LiDIA::subtract(lz[j][k], lz[j][k], tempbin0);
					}
					LiDIA::multiply(tempbin0, tempbin4, tempvectbin0[l + 1]);
					LiDIA::subtract(lz[l][k], lz[l][k], tempbin0);
				}
			}
			k++;
		}
	}
	delete[] lzdel;
	delete[] lz;
}



//
// Benne de Weger version of lll
// result : lll - reduced lattice for parameter a / b
//
void
bigint_lattice::TrD_lll_trans (dense_alg< bigint > & da, lattice_info& li)
{
	debug_handler("bigint_lattice", "TrD_lll_trans(da, li)");

	bigint tempbin0;
	bigint tempbin1;
	bigint tempbin2;
	bigint tempbin3;
	bigint tempbin4;
	bigint tempbin5;
	bigint *tempvectbin0;
	bigint *tempvectbin1;
	bigint *lzdel;
	bigint **lz;
	bigint **tempmatrixbin0;
	bigint **T;
	p_vector< bigint > vector;

//
// Allocating memory for matrix
//
	lz = new bigint*[3 * da.b.rows];
	memory_handler(lz, "bigint_lattice", "TrD_lll_trans(da, li) :: "
		       "not enough memory !");
	lzdel = new bigint[2 * da.b.rows * da.b.columns +
			  2 * da.b.rows + da.b.columns + 2];
	memory_handler(lzdel, "bigint_lattice", "TrD_lll_trans(da, li) :: "
		       "not enough memory !");
	T = &lz[da.b.rows];
	tempmatrixbin0 = &lz[da.b.rows * 2];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		lz[i] = &lzdel[i * da.b.columns];
		tempmatrixbin0[i] = &lzdel[(da.b.rows + i) * da.b.columns];
//
// Tr_dense_create() has created identic matrix
//
		T[i] = &da.s.value[i][da.b.columns];
	}

//
// Allocating memory for vectors
//
	tempvectbin0 = &lzdel[2 * da.b.rows * da.b.columns];
	tempvectbin1 = &tempvectbin0[da.b.rows + 1];


	li.lll.reduction_steps = 0;
	li.lll.swaps = 0;
	li.lll.correction_steps = 0;
	vector.vectsize = da.b.columns;
	vector.assign_zero(tempvectbin0);
	tempvectbin0[0].assign_one();


//
// Compute a bigint approximation of the Gram - Schmidt - Matrix
// Compute my`s by giving a matrix and a vector
//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		vector.assign(tempmatrixbin0[i], da.s.value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			vector.scalprod(lz[j][i], da.s.value[i], tempmatrixbin0[j]);
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				LiDIA::multiply(tempmatrixbin0[i][cx], tempvectbin0[j + 1],
						   tempmatrixbin0[i][cx]);
				LiDIA::multiply(tempvectbin1[cx], lz[j][i],
						   tempmatrixbin0[j][cx]);
				LiDIA::subtract(tempmatrixbin0[i][cx], tempmatrixbin0[i][cx],
						   tempvectbin1[cx]);
			}
			for (lidia_size_t m = 0; m < da.b.columns; m++) {
				LiDIA::divide(tempmatrixbin0[i][m], tempmatrixbin0[i][m],
						 tempvectbin0[j]);
			}
		}
		vector.scalprod(tempvectbin0[i + 1], tempmatrixbin0[i],
				tempmatrixbin0[i]);
		LiDIA::divide(tempvectbin0[i + 1], tempvectbin0[i + 1], tempvectbin0[i]);
	}

//
// Begin of lll - Algoritm
//
	lidia_size_t l, k = 1;
	while (k != da.b.rows) {
		l = k - 1;
		tempbin0.assign(abs(lz[l][k]));
		tempbin0.multiply_by_2();
		if (tempbin0.compare(tempvectbin0[l + 1]) > 0) {
			li.lll.reduction_steps++;
			tempbin1.assign(abs(tempvectbin0[l + 1]));
			LiDIA::add(tempbin5, tempbin0, tempbin1);

			LiDIA::divide(tempbin5, tempbin5, tempbin1);
			tempbin5.divide_by_2();

			if (lz[l][k].sign() != tempvectbin0[l + 1].sign()) {
				tempbin5.negate();
			}
			for (lidia_size_t cx = 0; cx < da.b.rows; cx++) {
				LiDIA::multiply(tempvectbin1[0], tempbin5, T[l][cx]);
				LiDIA::subtract(T[k][cx], T[k][cx], tempvectbin1[0]);
			}

			for (lidia_size_t j = 0; j < l; j++) {
				LiDIA::multiply(tempbin3, tempbin5, lz[j][l]);
				LiDIA::subtract(lz[j][k], lz[j][k], tempbin3);
			}
			LiDIA::multiply(tempbin3, tempbin5, tempvectbin0[l + 1]);
			LiDIA::subtract(lz[l][k], lz[l][k], tempbin3);
		}

		LiDIA::square(tempbin0, tempvectbin0[k]);
		LiDIA::multiply(tempbin1, tempbin0, da.b.y_nom);
		LiDIA::square(tempbin0, lz[l][k]);
		LiDIA::multiply(tempbin2, tempbin0, da.b.y_denom);
		LiDIA::subtract(tempbin0, tempbin1, tempbin2);
		LiDIA::multiply(tempbin1, tempvectbin0[k - 1], tempvectbin0[k + 1]);
		LiDIA::multiply(tempbin2, tempbin1, da.b.y_denom);

		if (tempbin2.compare(tempbin0) < 0) {
			li.lll.swaps++;
			vector.swap(T[k - 1], T[k]);
			for (lidia_size_t j = 0; j <= k - 2; j++) {
				LiDIA::swap(lz[j][k - 1], lz[j][k]);
			}
			for (lidia_size_t i = k + 1; i < da.b.rows; i++) {
				LiDIA::multiply(tempbin5, lz[k - 1][i], lz[k - 1][k]);
				LiDIA::multiply(tempbin4, lz[k][i], tempvectbin0[k - 1]);
				LiDIA::add(tempbin3, tempbin5, tempbin4);
				LiDIA::divide(tempbin5, tempbin3, tempvectbin0[k]);
				LiDIA::multiply(tempbin4, lz[k - 1][i], tempvectbin0[k + 1]);
				LiDIA::multiply(tempbin3, lz[k][i], lz[k - 1][k]);
				LiDIA::subtract(tempbin1, tempbin4, tempbin3);
				LiDIA::divide(tempbin4, tempbin1, tempvectbin0[k]);
				lz[k - 1][i].assign(tempbin5);
				lz[k][i].assign(tempbin4);
			}
			LiDIA::square(tempbin3, lz[k - 1][k]);
			LiDIA::multiply(tempbin5, tempvectbin0[k - 1], tempvectbin0[k + 1]);
			LiDIA::add(tempbin4, tempbin5, tempbin3);
			LiDIA::divide(tempvectbin0[k], tempbin4, tempvectbin0[k]);
			if (k > 1) {
				k--;
			}
		}
		else {
//
// Reduction step on transformation lattice
//
			for (l = k - 2; l >= 0; l--) {
				tempbin3.assign(abs(lz[l][k]));
				tempbin3.multiply_by_2();
				if (tempbin3.compare(tempvectbin0[l + 1]) > 0) {
//                  si = (lz[l][k]).sign() * (tempvectbin0[l + 1]).sign();
// Was passiert, wenn eines der beiden Vorzeichen Null ist ?
// Kann dies auftreten ?
//
					li.lll.reduction_steps++;
					tempbin5.assign(abs(tempvectbin0[l + 1]));
					LiDIA::add(tempbin1, tempbin3, tempbin5);
					LiDIA::divide(tempbin4, tempbin1, tempbin5);
					tempbin4.divide_by_2();
					if ((lz[l][k]).sign() != (tempvectbin0[l + 1]).sign()) {
						tempbin4.negate();
					}

					for (lidia_size_t cx = 0; cx < da.b.rows; cx++) {
						LiDIA::multiply(tempvectbin1[0], tempbin4, T[l][cx]);
						LiDIA::subtract(T[k][cx], T[k][cx], tempvectbin1[0]);
					}
					for (lidia_size_t j = 0; j < l; j++) {
						LiDIA::multiply(tempbin0, tempbin4, lz[j][l]);
						LiDIA::subtract(lz[j][k], lz[j][k], tempbin0);
					}
					LiDIA::multiply(tempbin0, tempbin4, tempvectbin0[l + 1]);
					LiDIA::subtract(lz[l][k], lz[l][k], tempbin0);
				}
			}
			k++;
		}
	}
//
// Free allocated storage
//
	delete[] lzdel;
	delete[] lz;
}



//
// Benne de Weger version of lll
// result : lll - reduced lattice for parameter y_nom / y_denom
//
void
bigint_lattice::TrD_lll_rand (dense_alg< bigint > & da, lattice_info& li)
{
	debug_handler("bigint_lattice", "Tr_lll_rand(da, li)");
	random_generator rg;
	bool erfolg;
	lidia_size_t *swaparray;
	lidia_size_t next_swap = -1, count, intervall;
	bigint tempbin0;
	bigint tempbin1;
	bigint tempbin2;
	bigint tempbin3;
	bigint tempbin4;
	bigint tempbin5;
	bigint *tempvectbin0;
	bigint *tempvectbin1;
	bigint *lzdel;
	bigint **lz;
	bigint **tempmatrixbin0;
	p_vector< bigint > vector;

//
//  allocating storage for lattices lz(da.b.rows, da.b.columns),
//                               and tempmatrixbin0(da.b.rows, da.b.columns);
//
	lz = new bigint * [2 * da.b.rows];
	memory_handler(lz, "bigint_lattice", "Tr_lll_rand(da, li) :: "
		       "not enough memory");
	lzdel = new bigint[da.b.rows * da.b.columns * 2 + da.b.rows + da.b.columns + 2];
	memory_handler(lzdel, "bigint_lattice", "Tr_lll_rand(da, li) :: "
		       "not enough memory");
	tempmatrixbin0 = &lz[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		lz[i] = &lzdel[i * da.b.columns];
		tempmatrixbin0[i] = &lzdel[(da.b.rows + i) * da.b.columns];
	}
	tempvectbin0 = &lzdel[2 * da.b.rows * da.b.columns];
	tempvectbin1 = &tempvectbin0[da.b.rows + 1];

	swaparray = new lidia_size_t[da.b.rows];
	memory_handler(swaparray, "bigint_lattice", "Tr_lll_rand(da, li) :: "
		       "not enough memory !");

//
// Clear swaparray
//
	for (lidia_size_t i = 0; i < da.b.rows; swaparray[i++] = 0)
		;
	li.lll.reduction_steps = 0;
	li.lll.swaps = 0;
	li.lll.correction_steps = 0;
	vector.vectsize = da.b.columns;

	vector.assign_zero(tempvectbin0);
	tempvectbin0[0].assign_one();


//
// Compute a bigint approximation of the Gram - Schmidt - Matrix
// Compute my`s by giving a matrix and a vector
//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		vector.assign(tempmatrixbin0[i], da.s.value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			vector.scalprod(lz[j][i], da.s.value[i], tempmatrixbin0[j]);
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				LiDIA::multiply(tempmatrixbin0[i][cx], tempvectbin0[j + 1],
						   tempmatrixbin0[i][cx]);
				LiDIA::multiply(tempvectbin1[cx], lz[j][i],
						   tempmatrixbin0[j][cx]);
				LiDIA::subtract(tempmatrixbin0[i][cx], tempmatrixbin0[i][cx],
						   tempvectbin1[cx]);
			}
			for (lidia_size_t m = 0; m < da.b.columns; m++)
				LiDIA::divide(tempmatrixbin0[i][m], tempmatrixbin0[i][m],
						 tempvectbin0[j]);
		}
		vector.scalprod(tempvectbin0[i + 1], tempmatrixbin0[i],
				tempmatrixbin0[i]);
		LiDIA::divide(tempvectbin0[i + 1], tempvectbin0[i + 1], tempvectbin0[i]);
	}

//
// begin of lll algorithm
//
	lidia_size_t k = 1;
	lidia_size_t l;

	rg >> l;
	l = static_cast<lidia_size_t>(l%(da.b.rows - 2));
	while (l < da.b.rows - 1) {
		k = l + 1;
		tempbin0.assign(abs(lz[l][k]));
		tempbin0.multiply_by_2();
		if (tempbin0.compare(tempvectbin0[l + 1]) > 0) {
			li.lll.reduction_steps++;
			tempbin1.assign(abs(tempvectbin0[l + 1]));
			LiDIA::add(tempbin5, tempbin0, tempbin1);

			LiDIA::divide(tempbin5, tempbin5, tempbin1);
			tempbin5.divide_by_2();

			if (lz[l][k].sign() != tempvectbin0[l + 1].sign()) {
				tempbin5.negate();
			}
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				LiDIA::multiply(tempvectbin1[cx], tempbin5, da.s.value[l][cx]);
				LiDIA::subtract(da.s.value[k][cx], da.s.value[k][cx],
						   tempvectbin1[cx]);
			}

			for (lidia_size_t j = 0; j < l; j++) {
				LiDIA::multiply(tempbin3, tempbin5, lz[j][l]);
				LiDIA::subtract(lz[j][k], lz[j][k], tempbin3);
			}
			LiDIA::multiply(tempbin3, tempbin5, tempvectbin0[l + 1]);
			LiDIA::subtract(lz[l][k], lz[l][k], tempbin3);
		}

//
// set mark after reduction
//
		swaparray[l] = 1;

		LiDIA::square(tempbin0, tempvectbin0[k]);
		LiDIA::multiply(tempbin1, tempbin0, da.b.y_nom);
		LiDIA::square(tempbin0, lz[l][k]);
		LiDIA::multiply(tempbin2, tempbin0, da.b.y_denom);
		LiDIA::subtract(tempbin0, tempbin1, tempbin2);
		LiDIA::multiply(tempbin1, tempvectbin0[k - 1], tempvectbin0[k + 1]);
		LiDIA::multiply(tempbin2, tempbin1, da.b.y_denom);

		if (tempbin2.compare(tempbin0) < 0) {
			li.lll.swaps++;
			vector.swap(da.s.value[k - 1], da.s.value[k]);
			for (lidia_size_t j = 0; j <= k - 2; j++) {
				LiDIA::swap(lz[j][k - 1], lz[j][k]);
			}
			for (lidia_size_t i = k + 1; i < da.b.rows; i++) {
				LiDIA::multiply(tempbin5, lz[k - 1][i], lz[k - 1][k]);
				LiDIA::multiply(tempbin4, lz[k][i], tempvectbin0[k - 1]);
				LiDIA::add(tempbin3, tempbin5, tempbin4);
				LiDIA::divide(tempbin5, tempbin3, tempvectbin0[k]);
				LiDIA::multiply(tempbin4, lz[k - 1][i], tempvectbin0[k + 1]);
				LiDIA::multiply(tempbin3, lz[k][i], lz[k - 1][k]);
				LiDIA::subtract(tempbin1, tempbin4, tempbin3);
				LiDIA::divide(tempbin4, tempbin1, tempvectbin0[k]);
				lz[k - 1][i].assign(tempbin5);
				lz[k][i].assign(tempbin4);
			}
			LiDIA::square(tempbin3, lz[k - 1][k]);
			LiDIA::multiply(tempbin5, tempvectbin0[k - 1], tempvectbin0[k + 1]);
			LiDIA::add(tempbin4, tempbin5, tempbin3);
			LiDIA::divide(tempvectbin0[k], tempbin4, tempvectbin0[k]);
//
// next index
//
			if (l > 0) {
				swaparray[l - 1] = 0;
			}
			if (l < da.b.rows - 2) {
				swaparray[l + 1] = 0;
			}
			if (l == 0) {
				next_swap = 0;
			}
			else {
				if (l == 1) {
					next_swap = -1;
				}
				else {
					if (next_swap >= l - 1) {
						next_swap = l - 2;
					}
				}
			}
		}
		else {
//
// correct next_swap
//
			if (l == next_swap + 1) {
				lidia_size_t cx = l;
				while (swaparray[cx++] == 1) {
					next_swap++;
				}
			}
			for (l = k - 2; l >= 0; l--) {
				tempbin3.assign(abs(lz[l][k]));
				tempbin3.multiply_by_2();
				if (tempbin3.compare(tempvectbin0[l + 1]) > 0) {
//                  si = (lz[l][k]).sign() * (tempvectbin0[l+1]).sign();
// Was passiert, wenn eines der beiden Vorzeichen Null ist ?
// Kann dies auftreten ?
//
					li.lll.reduction_steps++;
					tempbin5.assign(abs(tempvectbin0[l + 1]));
					LiDIA::add(tempbin1, tempbin3, tempbin5);
					LiDIA::divide(tempbin4, tempbin1, tempbin5);
					tempbin4.divide_by_2();
					if ((lz[l][k]).sign() != (tempvectbin0[l + 1]).sign()) {
						tempbin4.negate();
					}

					for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
						LiDIA::multiply(tempvectbin1[cx], tempbin4,
								   da.s.value[l][cx]);
						LiDIA::subtract(da.s.value[k][cx], da.s.value[k][cx],
								   tempvectbin1[cx]);
					}
					for (lidia_size_t j = 0; j < l; j++) {
						LiDIA::multiply(tempbin0, tempbin4, lz[j][l]);
						LiDIA::subtract(lz[j][k], lz[j][k], tempbin0);
					}
					LiDIA::multiply(tempbin0, tempbin4, tempvectbin0[l + 1]);
					LiDIA::subtract(lz[l][k], lz[l][k], tempbin0);
				}
			}
		}
//
// Fetch new index
//
		count = 0;
		erfolg = false;
		intervall = (da.b.rows - 1) - (next_swap + 2);
		while (count < intervall / 2 && (erfolg == false)) {
			rg >> l;
			// Note that STLport cannot resolve abs(int)
			l = ((l < 0 ? -l : l) % intervall) + next_swap + 2;
			erfolg = ((swaparray[l] == 0) ? true : false);
			count++;
		}
		if (erfolg == false) {
			l = next_swap + 1;
		}
	}
//
// Free allocated storage
//
	delete[] lzdel;
	delete[] lz;
	delete[] swaparray;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
