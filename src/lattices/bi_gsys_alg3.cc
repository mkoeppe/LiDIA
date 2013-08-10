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
#include	"LiDIA/lattices/bi_lattice_gensys.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// using `only` bigints
//
//
// Benne de Weger
//
// * Tr_lll()
// * Tr_lll_trans(T)
// * Tr_lll_rand()
//
// Buchmann Kessler
//
// * Tr_lin_gen_system(T)
//

//
// Benne de Weger version of lll
// result : lll - reduced lattice for parameter a / b
//
void bigint_lattice_gensys::Tr_lll()
{
	debug_handler("bigint_lattice_gensys", "Tr_lll()");


	bigint *tempvect0;
	bigint *tempvect1;
	bigint *lzdel;
	bigint **lz;
	bigint **help_matrix;

	//
	// Allocating memory for lattice
	//
	lz = new bigint*[2*rows];
	memory_handler(lz, "bigint_lattice_gensys", "Tr_lll() :: "
		       "not enough memory !");
	lzdel = new bigint[rows*columns+rows*rows+2*rows+columns+2];
	memory_handler(lzdel, "bigint_lattice_gensys", "Tr_lll() :: "
		       "not enough memory !");
	help_matrix = &lz[rows];
	for (lidia_size_t i = 0; i < rows; i++) {
		lz[i] = &lzdel[i*rows];
		help_matrix[i] = &lzdel[rows*rows+i*columns];
	}

	//
	// Allocating memory for vectors
	//
	tempvect0 = &lzdel[rows*rows+rows*columns];
	tempvect1 = &tempvect0[rows+1];


	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;

	bin_assign_zero_bin(tempvect0);

	tempvect0[0].assign_one();

	//
	// Compute a bigint approximation of the Gram - Schmidt - Matrix
	// Compute my`s by giving a matrix and a vector
	//
	for (lidia_size_t i = 0; i < rows; i++) {
		bin_assign_bin(help_matrix[i], value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			bin_scalprod_bin(lz[j][i], value[i], help_matrix[j]);
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(help_matrix[i][cx], tempvect0[j+1], help_matrix[i][cx]);
				LiDIA::multiply(tempvect1[cx], lz[j][i], help_matrix[j][cx]);
				LiDIA::subtract(help_matrix[i][cx], help_matrix[i][cx], tempvect1[cx]);
			}
			for (lidia_size_t m = 0; m < columns; m++)
				LiDIA::divide(help_matrix[i][m], help_matrix[i][m], tempvect0[j]);
		}
		bin_scalquad_bin(tempvect0[i+1], help_matrix[i]);
		LiDIA::divide(tempvect0[i+1], tempvect0[i+1], tempvect0[i]);
	}

	lidia_size_t l, k = 1;

	//
	// Begin of lll - Algoritm
	//
	while (k != rows) {
		l = k-1;
		tempmz0.assign(abs(lz[l][k]));
		tempmz0.multiply_by_2();
		if (tempmz0.compare(tempvect0[l+1]) > 0) {
			reduction_steps++;
			tempmz1.assign(abs(tempvect0[l+1]));
			LiDIA::add(ergmz, tempmz0, tempmz1);

			LiDIA::divide(ergmz, ergmz, tempmz1);
			ergmz.divide_by_2();

			if (lz[l][k].sign() != tempvect0[l+1].sign())
				ergmz.negate();
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(tempvect1[cx], ergmz, value[l][cx]);
				LiDIA::subtract(value[k][cx], value[k][cx], tempvect1[cx]);
			}

			for (lidia_size_t j = 0; j < l; j++) {
				LiDIA::multiply(tempmz3, ergmz, lz[j][l]);
				LiDIA::subtract(lz[j][k], lz[j][k], tempmz3);
			}
			LiDIA::multiply(tempmz3, ergmz, tempvect0[l+1]);
			LiDIA::subtract(lz[l][k], lz[l][k], tempmz3);
		}

		LiDIA::square(tempmz0, tempvect0[k]);
		LiDIA::multiply(tempmz1, tempmz0, y_nom);
		LiDIA::square(tempmz0, lz[l][k]);
		LiDIA::multiply(tempmz2, tempmz0, y_denom);
		LiDIA::subtract(tempmz0, tempmz1, tempmz2);
		LiDIA::multiply(tempmz1, tempvect0[k-1], tempvect0[k+1]);
		LiDIA::multiply(tempmz2, tempmz1, y_denom);

		if (tempmz2.compare(tempmz0) < 0) {
			bin_swap_bin(value[k-1], value[k]);
			++swapping_steps;
			for (lidia_size_t j = 0; j <= k-2; j++)
				LiDIA::swap(lz[j][k-1], lz[j][k]);
			for (lidia_size_t i = k+1; i < rows; i++) {
				LiDIA::multiply(ergmz, lz[k-1][i], lz[k-1][k]);
				LiDIA::multiply(tempmz4, lz[k][i], tempvect0[k-1]);
				LiDIA::add(tempmz3, ergmz, tempmz4);
				LiDIA::divide(ergmz, tempmz3, tempvect0[k]);
				LiDIA::multiply(tempmz4, lz[k-1][i], tempvect0[k+1]);
				LiDIA::multiply(tempmz3, lz[k][i], lz[k-1][k]);
				LiDIA::subtract(tempmz1, tempmz4, tempmz3);
				LiDIA::divide(tempmz4, tempmz1, tempvect0[k]);
				lz[k-1][i].assign(ergmz);
				lz[k][i].assign(tempmz4);
			}
			LiDIA::square(tempmz3, lz[k-1][k]);
			LiDIA::multiply(ergmz, tempvect0[k-1], tempvect0[k+1]);
			LiDIA::add(tempmz4, ergmz, tempmz3);
			LiDIA::divide(tempvect0[k], tempmz4, tempvect0[k]);

			if (k > 1)
				k--;
		}
		else {
			//
			// Reduction Step
			//
			for (l = k-2; l >= 0; l--) {
				tempmz3.assign(abs(lz[l][k]));
				tempmz3.multiply_by_2();
				if (tempmz3.compare(tempvect0[l+1]) > 0) {
					//                  si = (lz[l][k]).sign() * (tempvect0[l+1]).sign();
					// Was passiert, wenn eines der beiden Vorzeichen Null ist ?
					// Kann dies auftreten ?
					//
					reduction_steps++;
					ergmz.assign(abs(tempvect0[l+1]));
					LiDIA::add(tempmz1, tempmz3, ergmz);
					LiDIA::divide(tempmz4, tempmz1, ergmz);
					tempmz4.divide_by_2();
					if ((lz[l][k]).sign() != (tempvect0[l+1]).sign())
						tempmz4.negate();

					for (lidia_size_t cx = 0; cx < columns; cx++) {
						LiDIA::multiply(tempvect1[cx], tempmz4, value[l][cx]);
						LiDIA::subtract(value[k][cx], value[k][cx], tempvect1[cx]);
					}
					for (lidia_size_t j = 0; j < l; j++) {
						LiDIA::multiply(tempmz0, tempmz4, lz[j][l]);
						LiDIA::subtract(lz[j][k], lz[j][k], tempmz0);
					}
					LiDIA::multiply(tempmz0, tempmz4, tempvect0[l+1]);
					LiDIA::subtract(lz[l][k], lz[l][k], tempmz0);
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
void bigint_lattice_gensys::Tr_lll_trans(math_matrix< bigint > & Tr)
{
	debug_handler("bigint_lattice_gensys", "Tr_lll_trans(Tr)");


	bigint *tempvect0;
	bigint *tempvect1;
	bigint *lzdel;
	bigint **lz;
	bigint **help_matrix;
	bigint **TrAddr;
	bigint **T;

	//
	// Allocating memory for matrix
	//
	lz = new bigint*[3*rows];
	memory_handler(lz, "bigint_lattice_gensys", "Tr_lll_trans(Tr) :: "
		       "not enough memory !");
	lzdel = new bigint[2*rows*columns+rows*rows+2*rows+columns+2];
	memory_handler(lzdel, "bigint_lattice_gensys", "Tr_lll_trans(Tr) :: "
		       "not enough memory !");
	T = &lz[rows];
	help_matrix = &lz[rows*2];
	for (lidia_size_t i = 0; i < rows; i++) {
		lz[i] = &lzdel[i*columns];
		help_matrix[i] = &lzdel[(rows+i)*columns];
		T[i] = &lzdel[(2*rows*columns)+i*rows];
		T[i][i].assign_one();
	}

	//
	// Allocating memory for vectors
	//
	tempvect0 = &lzdel[2*rows*columns+rows*rows];
	tempvect1 = &tempvect0[rows+1];


	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;

	bin_assign_zero_bin(tempvect0);

	tempvect0[0].assign_one();


	//
	// Compute a bigint approximation of the Gram - Schmidt - Matrix
	// Compute my`s by giving a matrix and a vector
	//
	for (lidia_size_t i = 0; i < rows; i++) {
		bin_assign_bin(help_matrix[i], value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			bin_scalprod_bin(lz[j][i], value[i], help_matrix[j]);
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(help_matrix[i][cx], tempvect0[j+1], help_matrix[i][cx]);
				LiDIA::multiply(tempvect1[cx], lz[j][i], help_matrix[j][cx]);
				LiDIA::subtract(help_matrix[i][cx], help_matrix[i][cx], tempvect1[cx]);
			}
			for (lidia_size_t m = 0; m < columns; m++)
				LiDIA::divide(help_matrix[i][m], help_matrix[i][m], tempvect0[j]);
		}
		bin_scalquad_bin(tempvect0[i+1], help_matrix[i]);
		LiDIA::divide(tempvect0[i+1], tempvect0[i+1], tempvect0[i]);
	}

	//
	// Begin of lll - Algoritm
	//
	lidia_size_t l, k = 1;
	while (k != rows) {
		l = k-1;
		tempmz0.assign(abs(lz[l][k]));
		tempmz0.multiply_by_2();
		if (tempmz0.compare(tempvect0[l+1]) > 0) {
			reduction_steps++;
			tempmz1.assign(abs(tempvect0[l+1]));
			LiDIA::add(ergmz, tempmz0, tempmz1);

			LiDIA::divide(ergmz, ergmz, tempmz1);
			ergmz.divide_by_2();

			if (lz[l][k].sign() != tempvect0[l+1].sign())
				ergmz.negate();
			for (lidia_size_t cx = 0; cx < rows; cx++) {
				LiDIA::multiply(tempvect1[0], ergmz, T[l][cx]);
				LiDIA::subtract(T[k][cx], T[k][cx], tempvect1[0]);
			}

			for (lidia_size_t j = 0; j < l; j++) {
				LiDIA::multiply(tempmz3, ergmz, lz[j][l]);
				LiDIA::subtract(lz[j][k], lz[j][k], tempmz3);
			}
			LiDIA::multiply(tempmz3, ergmz, tempvect0[l+1]);
			LiDIA::subtract(lz[l][k], lz[l][k], tempmz3);
		}

		LiDIA::square(tempmz0, tempvect0[k]);
		LiDIA::multiply(tempmz1, tempmz0, y_nom);
		LiDIA::square(tempmz0, lz[l][k]);
		LiDIA::multiply(tempmz2, tempmz0, y_denom);
		LiDIA::subtract(tempmz0, tempmz1, tempmz2);
		LiDIA::multiply(tempmz1, tempvect0[k-1], tempvect0[k+1]);
		LiDIA::multiply(tempmz2, tempmz1, y_denom);

		if (tempmz2.compare(tempmz0) < 0) {
			swapping_steps++;
			bin_swap_bin(T[k-1], T[k]);
			for (lidia_size_t j = 0; j <= k-2; j++)
				LiDIA::swap(lz[j][k-1], lz[j][k]);
			for (lidia_size_t i = k+1; i < rows; i++) {
				LiDIA::multiply(ergmz, lz[k-1][i], lz[k-1][k]);
				LiDIA::multiply(tempmz4, lz[k][i], tempvect0[k-1]);
				LiDIA::add(tempmz3, ergmz, tempmz4);
				LiDIA::divide(ergmz, tempmz3, tempvect0[k]);
				LiDIA::multiply(tempmz4, lz[k-1][i], tempvect0[k+1]);
				LiDIA::multiply(tempmz3, lz[k][i], lz[k-1][k]);
				LiDIA::subtract(tempmz1, tempmz4, tempmz3);
				LiDIA::divide(tempmz4, tempmz1, tempvect0[k]);
				lz[k-1][i].assign(ergmz);
				lz[k][i].assign(tempmz4);
			}
			LiDIA::square(tempmz3, lz[k-1][k]);
			LiDIA::multiply(ergmz, tempvect0[k-1], tempvect0[k+1]);
			LiDIA::add(tempmz4, ergmz, tempmz3);
			LiDIA::divide(tempvect0[k], tempmz4, tempvect0[k]);
			if (k > 1)
				k--;
		}
		else {
			//
			// Reduction step on transformation lattice
			//
			for (l = k-2; l >= 0; l--) {
				tempmz3.assign(abs(lz[l][k]));
				tempmz3.multiply_by_2();
				if (tempmz3.compare(tempvect0[l+1]) > 0) {
					//                  si = (lz[l][k]).sign() * (tempvect0[l+1]).sign();
					// Was passiert, wenn eines der beiden Vorzeichen Null ist ?
					// Kann dies auftreten ?
					//
					reduction_steps++;
					ergmz.assign(abs(tempvect0[l+1]));
					LiDIA::add(tempmz1, tempmz3, ergmz);
					LiDIA::divide(tempmz4, tempmz1, ergmz);
					tempmz4.divide_by_2();
					if ((lz[l][k]).sign() != (tempvect0[l+1]).sign())
						tempmz4.negate();

					for (lidia_size_t cx = 0; cx < rows; cx++) {
						LiDIA::multiply(tempvect1[0], tempmz4, T[l][cx]);
						LiDIA::subtract(T[k][cx], T[k][cx], tempvect1[0]);
					}
					for (lidia_size_t j = 0; j < l; j++) {
						LiDIA::multiply(tempmz0, tempmz4, lz[j][l]);
						LiDIA::subtract(lz[j][k], lz[j][k], tempmz0);
					}
					LiDIA::multiply(tempmz0, tempmz4, tempvect0[l+1]);
					LiDIA::subtract(lz[l][k], lz[l][k], tempmz0);
				}
			}
			k++;
		}
	}
	//
	//  Store transformation lattice in math_matrix< bigint >
	//
	Tr.set_no_of_rows(rows);
	Tr.set_no_of_columns(rows);
	TrAddr = Tr.get_data_address();
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < rows; j++)
			if (trans_flag)
				LiDIA::swap(TrAddr[j][i], T[i][j]);
			else
				LiDIA::swap(TrAddr[j][i], T[j][i]);
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
void bigint_lattice_gensys::Tr_lll_rand()
{
	debug_handler("bigint_lattice_gensys", "Tr_lll_rand()");
	bool erfolg;
	lidia_size_t *swaparray;
	lidia_size_t next_swap = -1, count, intervall;
	bigint *tempvect0;
	bigint *tempvect1;
	bigint *lzdel;
	bigint **lz;
	bigint **help_matrix;

	//
	//  allocating storage for  lattices  lz(rows, columns),
	//                               and  help_matrix(rows, columns);
	//
	lz = new bigint*[2*rows];
	memory_handler(lz, "bigint_lattice_gensys", "Tr_lll_rand() :: "
		       "not enough memory");
	lzdel = new bigint[rows*columns*2+rows+columns+2];
	memory_handler(lzdel, "bigint_lattice_gensys", "Tr_lll_rand() :: "
		       "not enough memory");
	help_matrix = &lz[rows];
	for (lidia_size_t i = 0; i < rows; i++) {
		lz[i] = &lzdel[i*columns];
		help_matrix[i] = &lzdel[(rows+i)*columns];
	}
	tempvect0 = &lzdel[2*rows*columns];
	tempvect1 = &tempvect0[rows+1];

	swaparray = new lidia_size_t[rows];
	memory_handler(swaparray, "bigint_lattice_gensys", "Tr_lll_rand() :: "
		       "not enough memory !");

	//
	// Clear swaparray
	//
	for (lidia_size_t i = 0; i < rows; swaparray[i++] = 0);
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;

	bin_assign_zero_bin(tempvect0);

	tempvect0[0].assign_one();


	//
	// Compute a bigint approximation of the Gram - Schmidt - Matrix
	// Compute my`s by giving a matrix and a vector
	//
	for (lidia_size_t i = 0; i < rows; i++) {
		bin_assign_bin(help_matrix[i], value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			bin_scalprod_bin(lz[j][i], value[i], help_matrix[j]);
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(help_matrix[i][cx], tempvect0[j+1], help_matrix[i][cx]);
				LiDIA::multiply(tempvect1[cx], lz[j][i], help_matrix[j][cx]);
				LiDIA::subtract(help_matrix[i][cx], help_matrix[i][cx], tempvect1[cx]);
			}
			for (lidia_size_t m = 0; m < columns; m++)
				LiDIA::divide(help_matrix[i][m], help_matrix[i][m], tempvect0[j]);
		}
		bin_scalquad_bin(tempvect0[i+1], help_matrix[i]);
		LiDIA::divide(tempvect0[i+1], tempvect0[i+1], tempvect0[i]);
	}

	//
	// begin of lll algorithm
	//
	lidia_size_t k = 1;
	lidia_size_t l;
	random_generator rg;

	rg >> l;
	l = static_cast<lidia_size_t>(l%(rows-2));
	while (l < rows-1) {
		k = l+1;
		tempmz0.assign(abs(lz[l][k]));
		tempmz0.multiply_by_2();
		if (tempmz0.compare(tempvect0[l+1]) > 0) {
			reduction_steps++;
			tempmz1.assign(abs(tempvect0[l+1]));
			LiDIA::add(ergmz, tempmz0, tempmz1);

			LiDIA::divide(ergmz, ergmz, tempmz1);
			ergmz.divide_by_2();

			if (lz[l][k].sign() != tempvect0[l+1].sign())
				ergmz.negate();
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(tempvect1[cx], ergmz, value[l][cx]);
				LiDIA::subtract(value[k][cx], value[k][cx], tempvect1[cx]);
			}

			for (lidia_size_t j = 0; j < l; j++) {
				LiDIA::multiply(tempmz3, ergmz, lz[j][l]);
				LiDIA::subtract(lz[j][k], lz[j][k], tempmz3);
			}
			LiDIA::multiply(tempmz3, ergmz, tempvect0[l+1]);
			LiDIA::subtract(lz[l][k], lz[l][k], tempmz3);
		}

		//
		// set mark after reduction
		//
		swaparray[l] = 1;

		LiDIA::square(tempmz0, tempvect0[k]);
		LiDIA::multiply(tempmz1, tempmz0, y_nom);
		LiDIA::square(tempmz0, lz[l][k]);
		LiDIA::multiply(tempmz2, tempmz0, y_denom);
		LiDIA::subtract(tempmz0, tempmz1, tempmz2);
		LiDIA::multiply(tempmz1, tempvect0[k-1], tempvect0[k+1]);
		LiDIA::multiply(tempmz2, tempmz1, y_denom);

		if (tempmz2.compare(tempmz0) < 0) {
			swapping_steps++;
			bin_swap_bin(value[k-1], value[k]);
			for (lidia_size_t j = 0; j <= k-2; j++)
				LiDIA::swap(lz[j][k-1], lz[j][k]);
			for (lidia_size_t i = k+1; i < rows; i++) {
				LiDIA::multiply(ergmz, lz[k-1][i], lz[k-1][k]);
				LiDIA::multiply(tempmz4, lz[k][i], tempvect0[k-1]);
				LiDIA::add(tempmz3, ergmz, tempmz4);
				LiDIA::divide(ergmz, tempmz3, tempvect0[k]);
				LiDIA::multiply(tempmz4, lz[k-1][i], tempvect0[k+1]);
				LiDIA::multiply(tempmz3, lz[k][i], lz[k-1][k]);
				LiDIA::subtract(tempmz1, tempmz4, tempmz3);
				LiDIA::divide(tempmz4, tempmz1, tempvect0[k]);
				lz[k-1][i].assign(ergmz);
				lz[k][i].assign(tempmz4);
			}
			LiDIA::square(tempmz3, lz[k-1][k]);
			LiDIA::multiply(ergmz, tempvect0[k-1], tempvect0[k+1]);
			LiDIA::add(tempmz4, ergmz, tempmz3);
			LiDIA::divide(tempvect0[k], tempmz4, tempvect0[k]);
			//
			// next index
			//
			if (l > 0)
				swaparray[l-1] = 0;
			if (l < rows-2)
				swaparray[l+1] = 0;
			if (l == 0)
				next_swap = 0;
			else
				if (l == 1)
					next_swap = -1;
				else
					if (next_swap >= l-1)
						next_swap = l-2;


		}
		else {
			//
			// correct next_swap
			//
			if (l == next_swap+1) {
				lidia_size_t cx = l;
				while (swaparray[cx++] == 1)
					next_swap++;
			}
			for (l = k-2; l >= 0; l--) {
				tempmz3.assign(abs(lz[l][k]));
				tempmz3.multiply_by_2();
				if (tempmz3.compare(tempvect0[l+1]) > 0) {
					//                  si = (lz[l][k]).sign() * (tempvect0[l+1]).sign();
					// Was passiert, wenn eines der beiden Vorzeichen Null ist ?
					// Kann dies auftreten ?
					//
					reduction_steps++;
					ergmz.assign(abs(tempvect0[l+1]));
					LiDIA::add(tempmz1, tempmz3, ergmz);
					LiDIA::divide(tempmz4, tempmz1, ergmz);
					tempmz4.divide_by_2();
					if ((lz[l][k]).sign() != (tempvect0[l+1]).sign())
						tempmz4.negate();

					for (lidia_size_t cx = 0; cx < columns; cx++) {
						LiDIA::multiply(tempvect1[cx], tempmz4, value[l][cx]);
						LiDIA::subtract(value[k][cx], value[k][cx], tempvect1[cx]);
					}
					for (lidia_size_t j = 0; j < l; j++) {
						LiDIA::multiply(tempmz0, tempmz4, lz[j][l]);
						LiDIA::subtract(lz[j][k], lz[j][k], tempmz0);
					}
					LiDIA::multiply(tempmz0, tempmz4, tempvect0[l+1]);
					LiDIA::subtract(lz[l][k], lz[l][k], tempmz0);
				}
			}
		}
		//
		// Fetch new index
		//
		count = 0;
		erfolg = false;
		intervall = (rows-1)-(next_swap+2);
		while (count < intervall/2 && (erfolg == false)) {
			rg >> l;
			// Note that STLport cannot resolve abs(int)
			l = ((l<0?-l:l)%intervall)+next_swap+2;
			erfolg = ((swaparray[l] == 0)?true:false);
			count++;
		}
		if (erfolg == false)
			l = next_swap+1;
	}
	//
	// Free allocated storage
	//
	delete[] lzdel;
	delete[] lz;
	delete[] swaparray;
}



//
// Buchmann - Kessler version for generating systems
//
void bigint_lattice_gensys::Tr_lin_gen_system(math_matrix< bigint > & T, lidia_size_t& rk)
{
	debug_handler("bigint_lattice_gensys", "Tr_lin_gen_system(T, rk) [1]");

	bigint_lattice_gensys Atilde(rows, rows+columns), // The approximate lattice
		Ahead(rows, columns);

	bigint *help = new bigint[columns];
	bigint *rel = new bigint[rows];
	bigint *temp = new bigint[columns];
	bigint **TAddr;

	bigfloat vor, rechts;
	bigfloat zweipotq, alpha;
	bigfloat norm1, norm2;

	bigint bi_norm1, bi_norm2;
	bigint zwpq;

	sdigit n2 = rows;
	sdigit prec;
	sdigit bit_prec;

	//
	// Compute bigint approximation of lattice
	// to use the schnorr - euchner version of lll
	//
	prec = compute_precision();
	bigfloat::set_precision(prec);
	alpha_compute(alpha);
	zwei_pot_q_compute(zweipotq, n2, alpha);
	zweipotq.bigintify(zwpq);
	for (lidia_size_t i = 0; i < rows; ++i)
		for (lidia_size_t j = 0; j < columns; ++j)
			LiDIA::multiply(Ahead.value[i][j], zwpq, value[i][j]);

	for (lidia_size_t i = 0; i < rows; ++i) {
		Atilde.value[i][i].assign_one(); // 1, wenn i = j, 0 sonst
		for (lidia_size_t j = rows; j < rows+columns; ++j)
			Atilde.value[i][j].assign(Ahead.value[i][j-rows]);
	}

	//
	// Compute needed Precision for approximation bigfloats
	//
	prec = Atilde.compute_read_precision();
	bit_prec = static_cast<sdigit>(static_cast<double>(prec)*(std::log(10.0)/std::log(2.0)));
	bit_prec = ((bit_prec/300)+1)*52;

	//
	// Perform lll
	//
	debug_handler("bigint_lattice_gensys", "Tr_lin_gen_system(T, rk) [2]");
	Atilde.assign_the_rest(*this);
	Atilde.Tr_lll_bfl(bit_prec);
	assign_the_rest(Atilde);
	debug_handler("bigint_lattice_gensys", "Tr_lin_gen_system(T, rk) [3]");

	//
	// Check rank of the lattice
	//
	lidia_size_t l = 0;
	do {
		vectsize = columns;
		bin_assign_zero_bin(help); // Initializes help with the zero - vector
		for (lidia_size_t j = 0; j < rows; ++j) {
			rel[j].assign(Atilde.value[l][j]);
			bin_scalmul_bin(temp, rel[j], Ahead.value[j]);
			bin_add_bin(help, help, temp);
		}
		bin_l2_norm_bin(bi_norm2, help);
		norm2.assign(bi_norm2);
		vectsize = rows;
		bin_l1_norm_bin(bi_norm1, rel);
		norm1.assign(bi_norm1);
		sqrt(norm1, norm1);
		++l;
		LiDIA::divide(vor, n2, rows);
		vor.multiply_by_2();
		sqrt(vor, vor);
		LiDIA::multiply(vor, vor, norm1);
		LiDIA::divide(rechts, bigfloat(std::sqrt(static_cast<double>(n2))), bigfloat(2.0));
		LiDIA::multiply(rechts, rechts, norm1);
	}
	while ((zweipotq.compare(vor) > 0) && (norm2.compare(rechts) <= 0) && (l < rows));
	debug_handler("bigint_lattice_gensys", "Tr_lin_gen_system(T, rk) [4]");
	if (l >= rows)
		warning_handler("bigint_lattice_gensys", "Tr_lin_gen_system(T, rk) :: "
				"lattice of dimension 1");

	rk = rows-l+1; // rows is the dimension of the lattice
	//
	// Store transformation lattice in math_matrix< bigint >
	//
	T.set_no_of_rows(rows);
	T.set_no_of_columns(rows);
	TAddr = T.get_data_address();
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < rows-rk; j++)
			if (trans_flag)
				TAddr[i][j].assign_zero();
			else
				TAddr[j][i].assign_zero();
	for (lidia_size_t i = 0; i < rows; ++i)
		for (lidia_size_t j = rows-rk; j < rows; ++j)
			if (trans_flag)
				LiDIA::swap(TAddr[i][j], Atilde.value[j][i]);
			else
				LiDIA::swap(TAddr[j][i], Atilde.value[j][i]);

	//
	// Free allocated storage
	//
	delete[] rel;
	delete[] help;
	delete[] temp;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
