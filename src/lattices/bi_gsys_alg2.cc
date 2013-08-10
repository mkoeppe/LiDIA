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
// using `bigfloat` approximation
//
// Schnorr - Euchner
// * Tr_lll_bfl()
// * Tr_lll_trans_bfl(T)
// * Tr_lll_bfl_deep_insert() (not yet implemented, due this is not working !!)
// * Tr_lll_bfl_gensys()
// * Tr_lll_trans_bfl_gensys(T)
//
// Modified lll
// * Tr_mlll_bfl(v)
//

//
// Schnorr - Euchner version of lll optimized for bigints usign bigfloats
// result : lll - reduced lattice for parameter y
//
void bigint_lattice_gensys::Tr_lll_bfl(sdigit bit_prec)
{
	debug_handler("bigint_lattice", "Tr_lll_bfl(bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	sdigit prec;
	sdigit old_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	lidia_size_t k;
	bigint *tempvect1;
	bigfloat bound;
	bigfloat invbound;
	bigfloat tempbfl0;
	bigfloat tempbfl1;
	bigfloat tempbfl2;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat conv;
	bigfloat halb(0.5);

	//
	// Compute digit precision and bound for lll - reduction
	// also set precision
	//
	old_prec = bigfloat::get_precision();
	prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	bound = std::exp(std::log(2.0)*bit_prec/2);
	LiDIA::divide(invbound, 1, bound);
	bigfloat::set_precision(prec);

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_bfl(bit_prec) :: "
		       "not enough memory !");
	mydel = new bigfloat[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_bfl(bit_prec) :: "
                       "not enough memory !");
	lattice = &my[rows];
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	tempvect1 = new bigint[columns];
	memory_handler(tempvect1, "bigint_lattice_gensys", "Tr_lll_bfl(bit_prec) :: "
		       "not enough memory !");

	//
	// Setting pointers
	//
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into small bigfloat
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = value[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(ergmz, value[i][j], bit_diff);
				lattice[i][j].assign(ergmz);
				shift_left(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(value[i][j]);
		}

	//      lattice[i][j].assign(value[i][j]);

	k = 1;
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	vectsize = columns;
	bin_assign_zero_bfl(vect);
	bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_bfl(bit_prec) [2]");
			bin_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_bfl(tempbfl0, lattice[k], lattice[j]);
				LiDIA::multiply(tempbfl1, invbound, vect_sq[k]);
				LiDIA::multiply(tempbfl1, tempbfl1, vect_sq[j]);
				//
				// Correction Step I, only  if nessesary |lattice[k]*lattice[j]| to big
				//
				if (abs(tempbfl0).compare(tempbfl1) < 0) {
					//
					// Compute correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					bi_bit_len = tempmz0.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, tempmz0, bit_diff);
						my[j][k].assign(ergmz);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(tempmz0);

					//                my[j][k].assign(tempmz0);
				}
				else
					my[j][k].assign(tempbfl0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempbfl1, my[i][j], my[i][k]);
					LiDIA::multiply(tempbfl1, tempbfl1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempbfl1);
				}

				tempbfl0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempbfl1, my[j][k], tempbfl0);
				LiDIA::subtract(vect[k], vect[k], tempbfl1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_bfl(bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempbfl0, my[j][k]);
						if (abs(tempbfl0).compare(bound) > 0)
							Fc = true;

						tempbfl0.bigintify(tempmz0);
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempbfl1, tempbfl0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempbfl1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempbfl0);
					}
			}
			//
			// Correction Step II
			// needed after Reduction, correct only reduced vector
			// lattice[k] = (bigfloat ) value[k];
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = value[k][j].bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, value[k][j], bit_diff);
						lattice[k][j].assign(ergmz);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(value[k][j]);
				}
			}

			//
			// end of reduction
			//
			if (Fc) {
				if (k > 1)
					--k;
			}
			else {
				Fc = (!repeat_loop) && korr;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("bigint_lattice_gensys", "Tr_lll_bfl(bit_prec) [4]");
		tempbfl0.assign(y_par);
		LiDIA::multiply(tempbfl0, tempbfl0, vect[k-1]);
		LiDIA::square(tempbfl1, my[k-1][k]);
		LiDIA::multiply(tempbfl1, tempbfl1, vect[k-1]);
		LiDIA::add(tempbfl1, tempbfl1, vect[k]);

		//
		// Check lll - conditions for parameter y_par
		//
		if (tempbfl0.compare(tempbfl1) > 0) {
			//
			// ::swap vectors at position k and k-1
			//
			++swapping_steps;
			bin_swap_bin(value[k], value[k-1]);
			bin_swap_bfl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}
	//
	// Free allocated storage
	// and restore precision
	//
	bigfloat::set_precision(old_prec);
	delete[] mydel;
	delete[] my;
	delete[] tempvect1;
}



//
// Schnorr - Euchner version of lll optimized for bigints usign bigfloats
// result : lll - reduced lattice for parameter y
//
void bigint_lattice_gensys::Tr_lll_trans_bfl(math_matrix< bigint > & Tr, sdigit bit_prec)
{
	debug_handler("bigint_lattice", "Tr_lll_trans_bfl(Tr, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	sdigit prec;
	sdigit old_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	lidia_size_t k;
	bigint *tempvect1;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat bound;
	bigfloat invbound;
	bigfloat tempbfl0;
	bigfloat tempbfl1;
	bigfloat tempbfl2;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat conv;
	bigfloat halb(0.5);

	//
	// Computing digit precision, and bound thats needed by the algorithm
	// also the precision is set
	//
	old_prec = bigfloat::get_precision();
	prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	bound = std::exp(std::log(2.0)*bit_prec/2);
	LiDIA::divide(invbound, 1, bound);
	bigfloat::set_precision(prec);

	vectsize = columns;
	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "not enough memory !");
	mydel = new bigfloat[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
                       "not enough memory !");
	T = new bigint*[rows];
	memory_handler(T, "bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "not enough memory ! ");
	//
	// columns >= rows, or no basis
	//
	Tdel = new bigint[rows*rows+columns];
	memory_handler(Tdel, "bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
		T[i] = &Tdel[rows*i];
		T[i][i].assign_one();
	}

	//
	// Allocating Memory for vector
	// (-> Setting pointers)
	//
	tempvect1 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into small bigfloats
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = value[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(ergmz, value[i][j], bit_diff);
				lattice[i][j].assign(ergmz);
				shift_left(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(value[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	vectsize = columns;
	bin_assign_zero_bfl(vect);
	bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [2]");
			bin_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_bfl(tempbfl0, lattice[k], lattice[j]);
				LiDIA::multiply(tempbfl1, invbound, vect_sq[k]);
				LiDIA::multiply(tempbfl1, tempbfl1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]|
				//
				if (abs(tempbfl0).compare(tempbfl1) < 0) {
					//
					// Compute the correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);

					bi_bit_len = tempmz0.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, tempmz0, bit_diff);
						my[j][k].assign(ergmz);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(tempmz0);

					//                my[j][k].assign(tempmz0);
				}
				else
					my[j][k].assign(tempbfl0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempbfl1, my[i][j], my[i][k]);
					LiDIA::multiply(tempbfl1, tempbfl1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempbfl1);
				}

				tempbfl0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempbfl1, my[j][k], tempbfl0);
				LiDIA::subtract(vect[k], vect[k], tempbfl1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempbfl0, my[j][k]);
						if (abs(tempbfl0).compare(bound) > 0)
							Fc = true;

						tempbfl0.bigintify(tempmz0);
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);
						vectsize = rows;
						bin_scalmul_bin(tempvect1, tempmz0, T[j]);
						bin_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;
						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempbfl1, tempbfl0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempbfl1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempbfl0);
					}
			}

			//
			// Correction Step II
			// lattice[k] = (bigfloat ) value[k]
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = value[k][j].bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, value[k][j], bit_diff);
						lattice[k][j].assign(ergmz);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(value[k][j]);
				}
			}

			//
			// end of reduction
			//
			if (Fc) {
				if (k > 1)
					--k;
			}
			else {
				Fc = (!repeat_loop) && korr;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("bigint_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [4]");
		tempbfl0.assign(y_par);
		LiDIA::multiply(tempbfl0, tempbfl0, vect[k-1]);
		LiDIA::square(tempbfl1, my[k-1][k]);
		LiDIA::multiply(tempbfl1, tempbfl1, vect[k-1]);
		LiDIA::add(tempbfl1, tempbfl1, vect[k]);

		//
		// Check lll - condition
		//
		if (tempbfl0.compare(tempbfl1) > 0) {
			//
			// ::swap vectors k and k-1 ann don`t forget the transformation lattice
			//
			++swapping_steps;
			bin_swap_bin(value[k], value[k-1]);
			bin_swap_bin(T[k], T[k-1]);
			bin_swap_bfl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}

	//
	// Copy into math_matrix< bigint >
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
	// and restore precision
	//
	bigfloat::set_precision(old_prec);
	delete[] mydel;
	delete[] Tdel;
	delete[] my;
	delete[] T;
}



void bigint_lattice_gensys::Tr_lll_bfl_gensys(lidia_size_t& rank, sdigit bit_prec)
{
	debug_handler("bigint_lattice", "Tr_lll_bfl_gensys(rank, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	bool break_flag;
	sdigit prec;
	sdigit old_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	bigint *tempvect1;
	bigfloat bound;
	bigfloat invbound;
	bigfloat tempbfl0;
	bigfloat tempbfl1;
	bigfloat tempbfl2;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat conv;
	bigfloat halb(0.5);

	//
	// Compute digit precision and bound for the lll algorithm
	// also the precision is set
	//
	old_prec = bigfloat::get_precision();
	prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	bound = std::exp(std::log(2.0)*bit_prec/2);
	LiDIA::divide(invbound, 1, bound);
	bigfloat::set_precision(prec);
	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) :: "
		       "not enough memory !");
	mydel = new bigfloat[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) :: "
                       "not enough memory !");
	lattice = &my[rows];
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//
	tempvect1 = new bigint[columns];
	memory_handler(tempvect1, "bigint_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) :: "
		       "not enough memory !");

	//
	// Setting pointers
	//
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into small bigfloats
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = value[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(ergmz, value[i][j], bit_diff);
				lattice[i][j].assign(ergmz);
				shift_left(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(value[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	vectsize = columns;
	break_flag = false;
	bin_assign_zero_bfl(vect);
	bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [2]");
			bin_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_bfl(tempbfl0, lattice[k], lattice[j]);
				LiDIA::multiply(tempbfl1, invbound, vect_sq[k]);
				LiDIA::multiply(tempbfl1, tempbfl1, vect_sq[j]);
				//
				// Correction Step I, only if neeeded |lattice[k]*lattice[j]| to big
				//
				if (abs(tempbfl0).compare(tempbfl1) < 0) {
					//
					// Compute the correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					bi_bit_len = tempmz0.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, tempmz0, bit_diff);
						my[j][k].assign(ergmz);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(tempmz0);

					//                my[j][k].assign(tempmz0);
				}
				else
					my[j][k].assign(tempbfl0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempbfl1, my[i][j], my[i][k]);
					LiDIA::multiply(tempbfl1, tempbfl1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempbfl1);
				}

				tempbfl0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempbfl1, my[j][k], tempbfl0);
				LiDIA::subtract(vect[k], vect[k], tempbfl1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempbfl0, my[j][k]);
						if (abs(tempbfl0).compare(bound) > 0)
							Fc = true;

						tempbfl0.bigintify(tempmz0);
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempbfl1, tempbfl0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempbfl1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempbfl0);
					}
			}

			//
			// Correction Step II
			// lattice[k]=(bigfloat )value[k]
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = value[k][j].bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, value[k][j], bit_diff);
						lattice[k][j].assign(ergmz);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(value[k][j]);
				}

				bin_scalquad_bfl(vect[k], lattice[k]);
				sqrt(vect_sq[k], vect[k]);
			}
			//
			// Check if it`s a 0 - vector
			// delete this
			//
			if (vect_sq[k] < halb) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bin_swap_bfl(lattice[j], lattice[j+1]);
					bin_swap_bin(value[j], value[j+1]);
				}
				bin_assign_zero_bin(value[--mom_rows]);
				k = 1;
				bin_assign_zero_bfl(vect);
				bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
				break_flag = true;
				korr = false;
				Fc = false;
				break;
			}



			//
			// end of reduction
			//
			if (Fc) {
				if (k > 1)
					--k;
			}
			else {
				Fc = (!repeat_loop) && korr;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("bigint_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [4]");
		if (!break_flag) {
			tempbfl0.assign(y_par);
			LiDIA::multiply(tempbfl0, tempbfl0, vect[k-1]);
			LiDIA::square(tempbfl1, my[k-1][k]);
			LiDIA::multiply(tempbfl1, tempbfl1, vect[k-1]);
			LiDIA::add(tempbfl1, tempbfl1, vect[k]);

			//
			// Check lll - condition for parameter y_par
			//
			if (tempbfl0.compare(tempbfl1) > 0) {
				++swapping_steps;
				bin_swap_bin(value[k], value[k-1]);
				bin_swap_bfl(lattice[k], lattice[k-1]);
				tempbfl2.assign(vect_sq[k]);
				vect_sq[k].assign(vect_sq[k-1]);
				vect_sq[k-1].assign(tempbfl2);

				if (k > 1)
					--k;
				if (k == 1) {
					bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);
					sqrt(vect_sq[0], vect[0]);
				}
			}
			else
				++k;
		}
		break_flag = false;
	}
	rank = mom_rows;
	//
	// Free allocated storage
	// and restore precision
	//
	bigfloat::set_precision(old_prec);
	delete[] mydel;
	delete[] my;
	delete[] tempvect1;
}



void bigint_lattice_gensys::Tr_lll_trans_bfl_gensys(math_matrix< bigint > & Tr, lidia_size_t& rank, sdigit bit_prec)
{
	debug_handler("bigint_lattice", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	bool break_flag;
	sdigit prec;
	sdigit old_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	bigint *tempvect1;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat bound;
	bigfloat invbound;
	bigfloat tempbfl0;
	bigfloat tempbfl1;
	bigfloat tempbfl2;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat conv;
	bigfloat halb(0.5);

	//
	// Compute digit precision and bound for the lll algorithm
	// also the precision is set
	//
	old_prec = bigfloat::get_precision();
	prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	bound = std::exp(std::log(2.0)*bit_prec/2);
	LiDIA::divide(invbound, 1, bound);
	bigfloat::set_precision(prec);
	vectsize = columns;


	my = new bigfloat*[2*rows];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "not enough memory !");
	mydel = new bigfloat[(rows+2)*(rows+columns)-columns];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(rank, bit_prec) :: "
                       "not enough memory !");
	T = new bigint*[rows];
	memory_handler(T, "bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "not enough memory ! ");
	Tdel = new bigint[rows*rows+((columns > rows)?columns:rows)];
	memory_handler(Tdel, "bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
		T[i] = &Tdel[rows*i];
		T[i][i].assign_one();
	}

	//
	// Allocating Memory for vector
	// (-> Setting pointers)
	//
	tempvect1 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];


	//
	// Start
	//
	// Copy into small bigfloats
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = value[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(ergmz, value[i][j], bit_diff);
				lattice[i][j].assign(ergmz);
				shift_left(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(value[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	vectsize = columns;
	break_flag = false;
	bin_assign_zero_bfl(vect);
	bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [2]");
			bin_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_bfl(tempbfl0, lattice[k], lattice[j]);
				LiDIA::multiply(tempbfl1, invbound, vect_sq[k]);
				LiDIA::multiply(tempbfl1, tempbfl1, vect_sq[j]);
				//
				// Correction Step I, only if neeeded |lattice[k]*lattice[j]| to big
				//
				if (abs(tempbfl0).compare(tempbfl1) < 0) {
					//
					// Compute the correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					bi_bit_len = tempmz0.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, tempmz0, bit_diff);
						my[j][k].assign(ergmz);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(tempmz0);

					//                my[j][k].assign(tempmz0);
				}
				else
					my[j][k].assign(tempbfl0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempbfl1, my[i][j], my[i][k]);
					LiDIA::multiply(tempbfl1, tempbfl1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempbfl1);
				}

				tempbfl0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempbfl1, my[j][k], tempbfl0);
				LiDIA::subtract(vect[k], vect[k], tempbfl1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempbfl0, my[j][k]);
						if (abs(tempbfl0).compare(bound) > 0)
							Fc = true;

						tempbfl0.bigintify(tempmz0);
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);
						vectsize = rows;
						bin_scalmul_bin(tempvect1, tempmz0, T[j]);
						bin_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempbfl1, tempbfl0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempbfl1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempbfl0);
					}
			}

			//
			// Correction Step II
			// lattice[k]=(bigfloat )value[k]
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = value[k][j].bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(ergmz, value[k][j], bit_diff);
						lattice[k][j].assign(ergmz);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(value[k][j]);
				}

				bin_scalquad_bfl(vect[k], lattice[k]);
				sqrt(vect_sq[k], vect[k]);
			}
			//
			// Check if it`s a 0 - vector
			// delete this
			//
			if (vect_sq[k] < halb) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bin_swap_bfl(lattice[j], lattice[j+1]);
					bin_swap_bin(T[j], T[j+1]);
					bin_swap_bin(value[j], value[j+1]);
				}
				bin_assign_zero_bin(value[--mom_rows]);
				vectsize = rows;
				bin_assign_zero_bin(T[mom_rows]);
				vectsize = columns;
				k = 1;
				bin_assign_zero_bfl(vect);
				bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
				break_flag = true;
				korr = false;
				Fc = false;
				break;
			}

			//
			// end of reduction
			//
			if (Fc) {
				if (k > 1)
					--k;
			}
			else {
				Fc = (!repeat_loop) && korr;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("bigint_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [4]");
		if (!break_flag) {
			tempbfl0.assign(y_par);
			LiDIA::multiply(tempbfl0, tempbfl0, vect[k-1]);
			LiDIA::square(tempbfl1, my[k-1][k]);
			LiDIA::multiply(tempbfl1, tempbfl1, vect[k-1]);
			LiDIA::add(tempbfl1, tempbfl1, vect[k]);

			//
			// Check lll - condition for parameter y_par
			//
			if (tempbfl0.compare(tempbfl1) > 0) {
				++swapping_steps;
				bin_swap_bin(value[k], value[k-1]);
				bin_swap_bin(T[k], T[k-1]);
				bin_swap_bfl(lattice[k], lattice[k-1]);
				LiDIA::swap(vect_sq[k], vect_sq[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bin_scalprod_bfl(vect[0], lattice[0], lattice[0]);
					sqrt(vect_sq[0], vect[0]);
				}
			}
			else
				++k;
		}
		break_flag = false;
	}

	//
	// Copy into math_matrix< bigint >
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
	rank = mom_rows;
	//
	// Free allocated storage
	// and restore precision
	//
	bigfloat::set_precision(old_prec);
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}



//
// the modified lll-algorithm
//
void bigint_lattice_gensys::Tr_mlll_bfl(bigint* v)
{
	debug_handler("bigint_lattice_gensys", "Tr_mlll_bfl(v)");
	lidia_size_t j, k = columns, l, m;
	sdigit prec;
	sdigit old_prec;
	bigint *tempvect0;
	bigint *tempvect1;
	bigfloat *tempvect2;
	bigfloat *tempvect3;
	bigfloat *B;
	bigfloat halb(0.5);
	bigfloat dreiviertel(0.75);
	bigfloat tempbflz0, tempbflz1, Mu, Bz;
	bigfloat *mudel;
	bigfloat **mu, **bs;
	bigint *Trdel;
	bigint **Tr;

	//
	// Allocating needed storage
	//
	// Lattices
	//

	Tr = new bigint*[rows];
	memory_handler(T, "bigint_lattice_gensys", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");
	Trdel = new bigint[rows*rows+rows+columns];
	memory_handler(Trdel, "bigint_lattice_gensys", "Tr_mlll_bfl(v) :: "
                       "not enough memory !");

	mu = new bigfloat*[2*rows];
	memory_handler(mu, "bigint_lattice_gensys", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");
	mudel = new bigfloat[rows*((rows+columns)+3)];
	memory_handler(mudel, "bigint_lattice_gensys", "Tr_mlll_bfl(v) :: "
                       "not enough memory !");
	bs = &mu[rows];

	//
	// Restore and set precision
	//
	old_prec = bigfloat::get_precision();
	prec = compute_read_precision();
	halb.set_precision(rows*prec);
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	for (lidia_size_t i = 0; i < rows; i++) {
		Tr[i] = &Trdel[rows*i];
		Tr[i][i].assign_one();
		mu[i] = &mudel[rows*i];
		bs[i] = &mudel[rows*rows+i*columns];
	}


	//
	// Vectors
	// Setting pointers
	//
	tempvect2 = &mudel[(rows+columns)*rows];
	tempvect3 = &tempvect2[rows];
	B = &tempvect2[2*rows];

	tempvect0 = &Trdel[rows*rows];
	tempvect1 = &tempvect0[rows];

	//
	// Start Timer
	//
	vectsize = columns;


	for (lidia_size_t i = 0; i < k+1; i++) {
		for (lidia_size_t h = 0; h < i; h++) {
			mu[h][i].assign_zero();
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				vectbflz.assign(value[i][cx]);
				LiDIA::multiply(vectbflz, vectbflz, bs[h][cx]);
				LiDIA::add(mu[h][i], mu[h][i], vectbflz);
			}
			//
			//          for mixed datatypes
			//          bin_scalprod_bfl(mu[h][i], bs[h], value[i]);
			//
			LiDIA::divide(mu[h][i], mu[h][i], B[h]);
		}
		bin_assign_zero_bfl(tempvect3);
		for (lidia_size_t h = 0; h < i; h++) {
			bin_scalmul_bfl(tempvect2, mu[h][i], bs[h]);
			bin_add_bfl(tempvect3, tempvect3, tempvect2);
		}
		for (lidia_size_t cx = 0; cx < columns; cx++) {
			vectbflz.assign(value[i][cx]);
			LiDIA::subtract(bs[i][cx], vectbflz, tempvect3[cx]);
		}
		//
		//      Mixed datatypes
		//      bin_subtract_bfl(bs[i], value[i], tempvect3);
		//
		bin_scalprod_bfl(B[i], bs[i], bs[i]);
	}

	m = 1;
	l = 0;

	while (true) {
		//
		// Reduction Step
		//
		if (abs(mu[l][m]).compare(halb) > 0) {
			reduction_steps++;
			round(tempbflz0, mu[l][m]);
			tempbflz0.bigintify(ergmz);
			bin_scalmul_bin(tempvect0, ergmz, value[l]);
			bin_subtract_bin(value[m], value[m], tempvect0);
			vectsize = rows;
			bin_scalmul_bin(tempvect0, ergmz, Tr[l]);
			bin_subtract_bin(Tr[m], Tr[m], tempvect0);
			vectsize = columns;

			tempbflz1.assign(ergmz);
			LiDIA::subtract(mu[l][m], mu[l][m], tempbflz1);
			for (lidia_size_t h = 0; h < l; h++) {
				LiDIA::multiply(tempbflz0, tempbflz1, mu[h][l]);
				LiDIA::subtract(mu[h][m], mu[h][m], tempbflz0);
			}
		}


		j = 0;
		while ((j < k) && (value[m][j].is_zero())) {
			j++;
		}

		//
		// Exiting - condition
		//
		if (j == k) {
			for (lidia_size_t i = m; i < k; i++)
				bin_assign_bin(value[i], value[i+1]);
			bin_assign_zero_bin(value[rows-1]);
			vectsize = rows;
			bin_assign_bin(v, Tr[m]);
			vectsize = columns;
			//
			// Free allocated storage
			// and restore precision
			//
			bigfloat::set_precision(old_prec);
			delete[] Trdel;
			delete[] Tr;
			delete[] mudel;
			delete[] mu;
			return;
		}


		if (l >= m-1) {
			LiDIA::square(tempbflz1, mu[m-1][m]);
			LiDIA::subtract(tempbflz0, dreiviertel, tempbflz1);
			LiDIA::multiply(tempbflz0, tempbflz0, B[m-1]);

			if (B[m].compare(tempbflz0) < 0) {
				Mu.assign(mu[m-1][m]);
				LiDIA::square(tempbflz0, Mu);
				LiDIA::multiply(tempbflz0, tempbflz0, B[m-1]);
				LiDIA::add(Bz, B[m], tempbflz0);
				if (!Bz.is_approx_zero()) {
					LiDIA::divide(tempbflz0, B[m-1], Bz);
					LiDIA::multiply(mu[m-1][m], Mu, tempbflz0);
					LiDIA::multiply(B[m], B[m], tempbflz0);
					for (lidia_size_t i = m+1; i < k+1; i++) {
						LiDIA::multiply(tempbflz0, Mu, mu[m][i]);
						LiDIA::subtract(tempbflz1, mu[m-1][i], tempbflz0);
						LiDIA::multiply(tempbflz0, tempbflz1, mu[m-1][m]);
						LiDIA::add(mu[m-1][i], mu[m][i], tempbflz0);
						mu[m][i].assign(tempbflz1);
					}
				}
				B[m-1].assign(Bz);
				swapping_steps++;
				bin_swap_bin(value[m-1], value[m]);
				bin_swap_bin(Tr[m-1], Tr[m]);
				for (j = 0; j <= m-2; j++) {
					tempbflz1.assign(mu[j][m-1]);
					mu[j][m-1].assign(mu[j][m]);
					mu[j][m].assign(tempbflz1);
				}
				if (m > 1)
					m--;
				l = m-1;
				continue;
			}
		}
		l--;
		if (l < 0) {
			m++;
			l = m-1;
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
