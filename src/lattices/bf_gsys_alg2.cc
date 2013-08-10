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
#include	"LiDIA/lattices/bf_lattice_gensys.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// using small bigfloats for approximation
//
//
// Schnorr - Euchner
// * Tr_lll_bfl()
// * Tr_lll_var_bfl();
// * Tr_lll_trans_bfl(T)
// * Tr_lll_trans_var_bfl(T);
// * Tr_lll_deep_insert_bfl()  (not yet implemented)
// * Tr_lll_bfl_gensys()
// * Tr_lll_var_bfl_gensys()
// * Tr_lll_trans_bfl_gensys(T)
// * Tr_lll_trans_var_bfl_gensys(T)
//
// Modified lll
// * Tr_mlll_bfl(v)
//


//
// Schnorr - Euchner version of lll optimized for bigfloats using small bigfloats
// result : lll - reduced lattice for parameter y
//

//
// only high precision if needed
//
void bigfloat_lattice_gensys::Tr_lll_bfl(sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	sdigit old_prec;
	sdigit short_prec;
	sdigit r_prec;
	sdigit n_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	sdigit x_exponent;
	bigint x_factor;
	bigfloat *tempvect1;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat conv;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	r_prec = compute_read_precision();
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= columns; // columns >= rows, basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	//
	// think again
	//
	mydel = new bigfloat[rows*(rows+columns+2)+columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//
	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) :: "
		       "not enough memory !");

	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];


	//
	// Start
	//
	// Copy into floating point
	//
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			x_factor = value[i][j].mantissa();
			x_exponent = value[i][j].exponent();
			bi_bit_len = x_factor.bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(x_factor, x_factor, bit_diff);
				lattice[i][j].assign(x_factor);
				shift_right(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(x_factor);
			shift_left(lattice[i][j], lattice[i][j], x_exponent);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bigfloat::set_precision(n_prec);
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					x_factor = tempmz0.mantissa();
					x_exponent = tempmz0.exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						my[j][k].assign(x_factor);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(x_factor);
					shift_left(my[j][k], my[j][k], x_exponent);
					bigfloat::set_precision(short_prec);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						bigfloat::set_precision(n_prec);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);
						bigfloat::set_precision(short_prec);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);

					}
			}

			//
			// Correction Step II
			// lattice=(bigfloat )value
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					x_factor = value[k][j].mantissa();
					x_exponent = value[k][j].exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						lattice[k][j].assign(x_factor);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(x_factor);
					shift_left(lattice[k][j], lattice[k][j], x_exponent);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl(bit_prec) [4]");
		//
		// Check lll - condition for parameter y_par
		//
		tempmz0.assign(y_par);
		LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
		LiDIA::square(tempmz1, my[k-1][k]);
		LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
		LiDIA::add(tempmz1, tempmz1, vect[k]);

		if (tempmz0.compare(tempmz1) > 0) {
			//
			// ::swap vectors k and k-1
			//
			++swapping_steps;
			bfl_swap_bfl(value[k], value[k-1]);
			bfl_swap_bfl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}

	//
	// Free allocated storage and
	// restore precision
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}



//
// Schnorr - Euchner version of lll optimized for bigints usign bigfloats
// result : lll - reduced lattice for parameter y
//

//
// only high precision if needed
//
void bigfloat_lattice_gensys::Tr_lll_trans_bfl(math_matrix< bigint > & Tr, sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	sdigit old_prec;
	sdigit short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	sdigit r_prec;
	sdigit n_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	sdigit x_exponent;
	bigint x_factor;
	bigfloat *tempvect1;
	bigint *tempvect2;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat *vect;
	bigfloat *vect1;
	bigfloat *vect_sq;
	bigint conv;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);
	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= columns; // columns >= rows, basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	//
	// think again
	//
	mydel = new bigfloat[(rows+2)*(rows+columns)-columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
                       "not enough memory !");
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "Not enough memory !");

	Tdel = new bigint[rows*rows+rows];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "Not enough memory !");

	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
		T[i] = &Tdel[i*rows];
		T[i][i].assign_one();
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//
	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) :: "
		       "not enough memory !");

	tempvect2 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];
	vect1 = &vect[2*rows];


	//
	// Start
	//
	// Copy into floating point
	//
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			x_factor = value[i][j].mantissa();
			x_exponent = value[i][j].exponent();
			bi_bit_len = x_factor.bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(x_factor, x_factor, bit_diff);
				lattice[i][j].assign(x_factor);
				shift_right(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(x_factor);
			shift_left(lattice[i][j], lattice[i][j], x_exponent);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bigfloat::set_precision(n_prec);
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					x_factor = tempmz0.mantissa();
					x_exponent = tempmz0.exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						my[j][k].assign(x_factor);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(x_factor);
					shift_left(my[j][k], my[j][k], x_exponent);
					bigfloat::set_precision(short_prec);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						bigfloat::set_precision(n_prec);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);
						vectsize = rows;
						tempmz0.bigintify(conv);
						bfl_scalmul_bin(tempvect2, conv, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect2);
						vectsize = columns;
						bigfloat::set_precision(short_prec);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);
					}
			}

			//
			// Correction Step II
			// lattice[k]=(bigfloat )value[k]
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					x_factor = value[k][j].mantissa();
					x_exponent = value[k][j].exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						lattice[k][j].assign(x_factor);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(x_factor);
					shift_left(lattice[k][j], lattice[k][j], x_exponent);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl(Tr, bit_prec) [4]");
		//
		// Check lll - condition for parameter y_par
		//
		tempmz0.assign(y_par);
		LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
		LiDIA::square(tempmz1, my[k-1][k]);
		LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
		LiDIA::add(tempmz1, tempmz1, vect[k]);

		if (tempmz0.compare(tempmz1) > 0) {
			//
			// ::swap vectors k and k-1
			//
			++swapping_steps;
			bfl_swap_bfl(value[k], value[k-1]);
			bfl_swap_bfl(lattice[k], lattice[k-1]);
			bfl_swap_bin(T[k], T[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}
	//
	// Copy into math_matrix< bigint > Tr
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
	// and restore old precision
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}



//
// Schnorr - Euchner version of lll optimized for bigfloats using small bigfloats
// result : lll - reduced lattice for parameter y
//

//
// only high precision if needed
//
void bigfloat_lattice_gensys::Tr_lll_bfl_gensys(lidia_size_t& rank, sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool repeat_loop;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	sdigit old_prec;
	sdigit short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	sdigit r_prec;
	sdigit n_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	sdigit x_exponent;
	bigint x_factor;
	bigfloat *tempvect1;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect1;
	bigfloat *vect_sq;
	bigfloat conv;
	bigfloat u_zero;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);
	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	log(u_zero, 10);
	LiDIA::multiply(u_zero, u_zero, r_prec);
	exp(u_zero, u_zero);
	LiDIA::divide(u_zero, 1, u_zero);
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= columns; // columns >= rows, basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	//
	// think again
	//
	mydel = new bigfloat[(rows+2)*(rows+columns)-columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//
	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) :: "
		       "not enough memory !");

	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];
	vect1 = &vect[2*rows];


	//
	// Start
	//
	// Copy into floating point
	//
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			x_factor = value[i][j].mantissa();
			x_exponent = value[i][j].exponent();
			bi_bit_len = x_factor.bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(x_factor, x_factor, bit_diff);
				lattice[i][j].assign(x_factor);
				shift_right(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(x_factor);
			shift_left(lattice[i][j], lattice[i][j], x_exponent);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bigfloat::set_precision(n_prec);
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					x_factor = tempmz0.mantissa();
					x_exponent = tempmz0.exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						my[j][k].assign(x_factor);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(x_factor);
					shift_left(my[j][k], my[j][k], x_exponent);
					bigfloat::set_precision(short_prec);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						bigfloat::set_precision(n_prec);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);
						bigfloat::set_precision(short_prec);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);

					}
			}

			//
			// Correction Step II
			// lattice=(bigfloat )value
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					x_factor = value[k][j].mantissa();
					x_exponent = value[k][j].exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						lattice[k][j].assign(x_factor);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(x_factor);
					shift_left(lattice[k][j], lattice[k][j], x_exponent);
				}
				bfl_scalquad_bfl(vect[k], lattice[k]);
				sqrt(vect_sq[k], vect_sq[k]);
			}
			//
			// Check if it`s 0 - vector
			// delete it
			//
			if (abs(vect_sq[k]) < u_zero) {
				for (lidia_size_t j = 0; j < mom_rows-1; j++) {
					bfl_swap_bfl(lattice[j], lattice[j+1]);
					bfl_swap_bfl(value[j], value[j+1]);
				}
				bfl_assign_zero_bfl(value[--mom_rows]);
				bfl_assign_zero_bfl(lattice[mom_rows]);
				k = 1;
				bfl_assign_zero_bfl(vect);
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
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

		if (!break_flag) {
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_bfl_gensys(rank, bit_prec) [4]");
			//
			// Check lll - condition for parameter y_par
			//
			tempmz0.assign(y_par);
			LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
			LiDIA::square(tempmz1, my[k-1][k]);
			LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
			LiDIA::add(tempmz1, tempmz1, vect[k]);

			if (tempmz0.compare(tempmz1) > 0) {
				//
				// ::swap vectors k and k-1
				//
				++swapping_steps;
				bfl_swap_bfl(value[k], value[k-1]);
				bfl_swap_bfl(lattice[k], lattice[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
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
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}



//
// Schnorr - Euchner version of lll optimized for bigints usign bigfloats
// result : lll - reduced lattice for parameter y
//

//
// only high precision if needed
//
void bigfloat_lattice_gensys::Tr_lll_trans_bfl_gensys(math_matrix< bigint > & Tr, lidia_size_t& rank,
						      sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool repeat_loop;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	sdigit old_prec;
	sdigit short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	sdigit r_prec;
	sdigit n_prec;
	sdigit bi_bit_len;
	sdigit bit_diff;
	sdigit x_exponent;
	bigint x_factor;
	bigfloat *tempvect1;
	bigint *tempvect2;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat u_zero;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigint conv;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);
	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	log(u_zero, 10);
	LiDIA::multiply(u_zero, u_zero, r_prec);
	exp(u_zero, u_zero);
	LiDIA::divide(u_zero, 1, u_zero);
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= columns; // columns >= rows, basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	//
	// think again
	//
	mydel = new bigfloat[(rows+2)*(rows+columns)-columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
                       "not enough memory !");
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "Not enough memory !");

	Tdel = new bigint[rows*rows+rows];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "Not enough memory !");

	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
		T[i] = &Tdel[i*rows];
		T[i][i].assign_one();
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//
	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) :: "
		       "not enough memory !");

	tempvect2 = &T[0][rows*rows];
	vect = &my[0][(rows+columns)*rows];
	vect_sq = &vect[rows];


	//
	// Start
	//
	// Copy into floating point
	//
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			x_factor = value[i][j].mantissa();
			x_exponent = value[i][j].exponent();
			bi_bit_len = x_factor.bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(x_factor, x_factor, bit_diff);
				lattice[i][j].assign(x_factor);
				shift_right(lattice[i][j], lattice[i][j], bit_diff);
			}
			else
				lattice[i][j].assign(x_factor);
			shift_left(lattice[i][j], lattice[i][j], x_exponent);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bigfloat::set_precision(n_prec);
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					x_factor = tempmz0.mantissa();
					x_exponent = tempmz0.exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						my[j][k].assign(x_factor);
						shift_left(my[j][k], my[j][k], bit_diff);
					}
					else
						my[j][k].assign(x_factor);
					shift_left(my[j][k], my[j][k], x_exponent);
					bigfloat::set_precision(short_prec);

				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						bigfloat::set_precision(n_prec);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect2, conv, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect2);
						vectsize = columns;
						bigfloat::set_precision(short_prec);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);
					}
			}

			//
			// Correction Step II
			// lattice[k]=(bigfloat )value[k]
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					x_factor = value[k][j].mantissa();
					x_exponent = value[k][j].exponent();
					bi_bit_len = x_factor.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(x_factor, x_factor, bit_diff);
						lattice[k][j].assign(x_factor);
						shift_left(lattice[k][j], lattice[k][j], bit_diff);
					}
					else
						lattice[k][j].assign(x_factor);
					shift_left(lattice[k][j], lattice[k][j], x_exponent);
				}

				bfl_scalquad_bfl(vect[k], lattice[k]);
				sqrt(vect_sq[k], vect_sq[k]);
			}
			//
			// Check if it`s 0 - vector
			// delete it
			//
			if (abs(vect_sq[k]) < u_zero) {
				for (lidia_size_t j = 0; j < mom_rows-1; j++) {
					bfl_swap_bfl(lattice[j], lattice[j+1]);
					bfl_swap_bfl(value[j], value[j+1]);
					bfl_swap_bin(T[j], T[j+1]);
				}
				bfl_assign_zero_bfl(value[--mom_rows]);
				bfl_assign_zero_bin(T[mom_rows]);
				k = 1;
				bfl_assign_zero_bfl(vect);
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
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

		if (!break_flag) {
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_bfl_gensys(Tr, rank, bit_prec) [4]");
			//
			// Check lll - condition for parameter y_par
			//
			tempmz0.assign(y_par);
			LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
			LiDIA::square(tempmz1, my[k-1][k]);
			LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
			LiDIA::add(tempmz1, tempmz1, vect[k]);

			if (tempmz0.compare(tempmz1) > 0) {
				//
				// ::swap vectors k and k-1
				//
				++swapping_steps;
				bfl_swap_bfl(value[k], value[k-1]);
				bfl_swap_bfl(lattice[k], lattice[k-1]);
				bfl_swap_bin(T[k], T[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
					sqrt(vect_sq[0], vect[0]);
				}
			}
			else
				++k;
		}
		break_flag = false;
	}
	//
	// Copy into math_matrix< bigint > Tr
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
	// and restore old precsion
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}



//
// Buchmann Kessler version of lll for generating systems
// result : transformation lattice
//
void bigfloat_lattice_gensys::Tr_lin_gen_system(math_matrix< bigint > & T, lidia_size_t& rk)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lin_gen_system(T, rk) [1]");

	bigint_lattice_gensys Atilde(rows, rows+columns), // The approximate lattice
		Ahead(rows, columns);
	bigint *help = new bigint[columns];
	bigint *rel = new bigint[rows];
	bigint *temp = new bigint[columns];

	bigfloat vor, rechts;
	bigfloat zweipotq, alpha;
	bigfloat norm1, norm2;

	bigint bi_norm1, bi_norm2;
	bigint zero;
	bigint **TAddr;

	sdigit n2 = rows;
	sdigit old_prec;
	sdigit prec;
	sdigit bit_prec;


	zero.assign_zero();
	old_prec = bigfloat::get_precision();
	prec = compute_precision();
	bigfloat::set_precision(prec);
	alpha_compute(alpha);
	zwei_pot_q_compute(zweipotq, n2, alpha);


	for (lidia_size_t i = 0; i < rows; ++i)
		for (lidia_size_t j = 0; j < columns; ++j) {
			LiDIA::multiply(tempmz0, zweipotq, value[i][j]);
			tempmz0.bigintify(Ahead.value[i][j]);
		}
	for (lidia_size_t i = 0; i < rows; ++i) {
		Atilde.value[i][i].assign_one(); // 1, wenn i = j, 0 sonst
		for (lidia_size_t j = rows; j < (rows+columns); ++j)
			Atilde.value[i][j].assign(Ahead.value[i][j-rows]);
	}

	prec = Atilde.compute_read_precision();
	bit_prec = static_cast<sdigit>(static_cast<double>(prec)*(std::log(10.0)/std::log(2.0)));
	bit_prec = ((bit_prec/300)+1)*52;

	debug_handler("bigfloat_lattice_gensys", "Tr_lin_gen_system(T, rk) [2]");
	Atilde.y_par = y_par;
	Atilde.Tr_lll_bfl(bit_prec);
	debug_handler("bigfloat_lattice_gensys", "Tr_lin_gen_system(T, rk) [3]");
	reduction_steps = Atilde.reduction_steps;
	swapping_steps = Atilde.swapping_steps;
	correction_steps = Atilde.correction_steps;


	lidia_size_t l = 0;
	do {
		vectsize = columns;
		bfl_assign_zero_bin(help); // Initializes help with the zero - vector
		for (lidia_size_t j = 0; j < rows; ++j) {
			rel[j].assign(Atilde.value[l][j]);
			bfl_scalmul_bin(temp, rel[j], Ahead.value[j]);
			bfl_add_bin(help, help, temp);
		}
		bfl_l2_norm_bin(bi_norm2, help);
		norm2.assign(bi_norm2);
		vectsize = rows;
		bfl_l1_norm_bin(bi_norm1, rel);
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
	debug_handler("bigfloat_lattice_gensys", "Tr_lin_gen_system(T, rk) [4]");
	if (l >= rows)
		warning_handler("bigfloat_lattice_gensys", "Tr_lin_gen_system(T, rk) :: "
				"lattice of dimension 1");

	rk = rows-l+1; // rows is the dimension of the lattice
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
	// and restore old precsion
	//
	bigfloat::set_precision(old_prec);
	delete[] help;
	delete[] rel;
	delete[] temp;
}



//
// the modified lll-algorithm
//
void bigfloat_lattice_gensys::Tr_mlll_bfl(bigint* v)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_mlll_bfl(v)");
	lidia_size_t j, k = columns, l, m;
	sdigit old_prec;
	sdigit prec;
	bigint *tempvect0;
	bigfloat *tempvect1;
	bigfloat *tempvect2;
	bigfloat *tempvect3;
	bigfloat *B;
	bigfloat halb(0.5);
	bigfloat dreiviertel(0.75);
	bigfloat Mu, Bz;
	bigint tempbin0;
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
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_mll_bfl(v) :: "
		       "not enough memory !");
	mu = new bigfloat*[2*rows];
	memory_handler(mu, "bigfloat_lattice_gensys", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");

	bs = &mu[rows];
	mudel = new bigfloat[rows*(rows+columns)+3*rows+columns];
	memory_handler(mudel, "bigfloat_lattice_gensys", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");
	Trdel = new bigint[rows*rows+rows];
	memory_handler(Trdel, "bigfloat_lattice_gensys", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");


	old_prec = bigfloat::get_precision();
	prec = compute_read_precision();
	bigfloat::set_precision(rows*prec);
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	for (lidia_size_t i = 0; i < rows; i++) {
		Tr[i] = &Trdel[i*rows];
		Tr[i][i].assign_one();
		mu[i] = &mudel[i*rows];
		bs[i] = &mudel[rows*rows+i*columns];
	}


	//
	// Vectors
	//
	tempvect2 = &mudel[rows*(rows+columns)];
	tempvect3 = &tempvect2[rows];
	B = &tempvect2[2*rows];

	tempvect0 = &Trdel[rows*rows];
	tempvect1 = &tempvect2[3*rows];

	//
	// Start Timer
	//
	vectsize = columns;


	for (lidia_size_t i = 0; i < k+1; i++) {
		for (lidia_size_t h = 0; h < i; h++) {
			mu[h][i].assign_zero();
			bfl_scalprod_bfl(mu[h][i], bs[h], value[i]);

			LiDIA::divide(mu[h][i], mu[h][i], B[h]);
		}
		bfl_assign_zero_bfl(tempvect3);
		for (lidia_size_t h = 0; h < i; h++) {
			bfl_scalmul_bfl(tempvect2, mu[h][i], bs[h]);
			bfl_add_bfl(tempvect3, tempvect3, tempvect2);
		}
		bfl_subtract_bfl(bs[i], value[i], tempvect3);
		bfl_scalprod_bfl(B[i], bs[i], bs[i]);
	}

	m = 1;
	l = 0;

	while (true) {
		//
		// Reduction Step
		//
		if (abs(mu[l][m]).compare(halb) > 0) {
			reduction_steps++;
			round(tempmz0, mu[l][m]);
			tempmz0.bigintify(tempbin0);
			bfl_scalmul_bfl(tempvect1, tempmz0, value[l]);
			bfl_subtract_bfl(value[m], value[m], tempvect1);
			vectsize = rows;
			bfl_scalmul_bin(tempvect0, tempbin0, Tr[l]);
			bfl_subtract_bin(Tr[m], Tr[m], tempvect0);
			vectsize = columns;

			LiDIA::subtract(mu[l][m], mu[l][m], tempmz0);
			for (lidia_size_t h = 0; h < l; h++) {
				LiDIA::multiply(tempmz1, tempmz0, mu[h][l]);
				LiDIA::subtract(mu[h][m], mu[h][m], tempmz1);
			}
		}


		j = 0;
		while ((j < k) && (value[m][j].is_approx_zero())) {
			j++;
		}

		//
		// Exiting - condition
		//
		if (j == k) {
			for (lidia_size_t i = m; i < k; i++)
				bfl_assign_bfl(value[i], value[i+1]);
			bfl_assign_zero_bfl(value[rows-1]);
			vectsize = rows;
			bfl_assign_bin(v, Tr[m]);
			vectsize = columns;
			bigfloat::set_precision(old_prec);
			delete[] mudel;
			delete[] mu;
			delete[] Trdel;
			delete[] Tr;
			return;
		}


		if (l >= m-1) {
			LiDIA::square(tempmz1, mu[m-1][m]);
			LiDIA::subtract(tempmz0, dreiviertel, tempmz1);
			LiDIA::multiply(tempmz0, tempmz0, B[m-1]);

			if (B[m].compare(tempmz0) < 0) {
				Mu.assign(mu[m-1][m]);
				LiDIA::square(tempmz0, Mu);
				LiDIA::multiply(tempmz0, tempmz0, B[m-1]);
				LiDIA::add(Bz, B[m], tempmz0);
				if (!Bz.is_approx_zero()) {
					LiDIA::divide(tempmz0, B[m-1], Bz);
					LiDIA::multiply(mu[m-1][m], Mu, tempmz0);
					LiDIA::multiply(B[m], B[m], tempmz0);
					for (lidia_size_t i = m+1; i < k+1; i++) {
						LiDIA::multiply(tempmz0, Mu, mu[m][i]);
						LiDIA::subtract(tempmz1, mu[m-1][i], tempmz0);
						LiDIA::multiply(tempmz0, tempmz1, mu[m-1][m]);
						LiDIA::add(mu[m-1][i], mu[m][i], tempmz0);
						mu[m][i].assign(tempmz1);
					}
				}
				B[m-1].assign(Bz);
				swapping_steps++;
				bfl_swap_bfl(value[m-1], value[m]);
				bfl_swap_bin(Tr[m-1], Tr[m]);
				for (j = 0; j <= m-2; j++) {
					tempmz1.assign(mu[j][m-1]);
					mu[j][m-1].assign(mu[j][m]);
					mu[j][m].assign(tempmz1);
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



void bigfloat_lattice_gensys::Tr_lll_var_bfl(math_matrix< bigint > & Bi,
					     sdigit bit_len, sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, bit_len, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	sdigit short_prec;
	lidia_size_t k;
	bigint tempbiz;
	bigint *tempvect1;
	bigint **BiAddr;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	BiAddr = Bi.get_data_address();

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	//
	// think again
	//
	mydel = new bigfloat[rows*(rows+columns+2)+columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, Factor, bit_prec) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//
	tempvect1 = new bigint[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, Factor, bit_prec) :: "
		       "not enough memory !");

	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];


	//
	// Start
	//
	// Copy into floating point
	//
	old_prec = bigfloat::get_precision();
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(lattice[i][j], tempmz0, bit_len-bit_diff);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, Factor, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(my[j][k], tempmz0, 2*bit_len-bit_diff);
					//                ::divide(my[j][k], tempbiz, SqrFactor);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;

			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, Factor, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);

					}
			}

			//
			// Correction Step II
			// lattice=(bigfloat )value
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = BiAddr[k][j].bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(tempbiz, BiAddr[k][j], bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(BiAddr[k][j]);
					}
					shift_right(lattice[k][j], tempmz0, bit_len-bit_diff);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl(Bi, Factor, bit_prec) [4]");
		//
		// Check lll - condition for parameter y_par
		//
		tempmz0.assign(y_par);
		LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
		LiDIA::square(tempmz1, my[k-1][k]);
		LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
		LiDIA::add(tempmz1, tempmz1, vect[k]);

		if (tempmz0.compare(tempmz1) > 0) {
			//
			// ::swap vectors k and k-1
			//
			++swapping_steps;
			bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
			bfl_swap_bfl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}
	bigfloat::set_precision(old_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > old_bit_prec) {
				bit_diff = bi_bit_len-old_bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(value[i][j], tempmz0, bit_len-bit_diff);
		}

	//
	// Free allocated storage
	// and restore old precision
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}



void bigfloat_lattice_gensys::Tr_lll_trans_var_bfl(math_matrix< bigint > & Bi, math_matrix< bigint > & Tr,
						   sdigit bit_len, sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_len, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	sdigit short_prec;
	lidia_size_t k;
	bigint tempbiz;
	bigint *tempvect1;
	bigint **BiAddr;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	BiAddr = Bi.get_data_address();

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_prec) :: "
		       "Not enough memory !");


	//
	// think again
	//
	mydel = new bigfloat[rows*(rows+columns+2)+columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_prec) :: "
                       "not enough memory !");
	Tdel = new bigint[rows*rows+columns];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_prec) :: "
		       "Not enough memory !");

	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		T[i] = &Tdel[i*rows];
		T[i][i].assign_one();
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//

	tempvect1 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	old_prec = bigfloat::get_precision();
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(lattice[i][j], tempmz0, bit_len-bit_diff);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Factor, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(my[j][k], tempmz0, 2*bit_len-bit_diff);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect1, tempbiz, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);

					}
			}

			//
			// Correction Step II
			// lattice=(bigfloat )value
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = BiAddr[k][j].bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(tempbiz, BiAddr[k][j], bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(BiAddr[k][j]);
					}
					shift_right(lattice[k][j], tempmz0, bit_len-bit_diff);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl(Bi, Tr, Factor, bit_prec) [4]");
		//
		// Check lll - condition for parameter y_par
		//
		tempmz0.assign(y_par);
		LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
		LiDIA::square(tempmz1, my[k-1][k]);
		LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
		LiDIA::add(tempmz1, tempmz1, vect[k]);

		if (tempmz0.compare(tempmz1) > 0) {
			//
			// ::swap vectors k and k-1
			//
			++swapping_steps;
			bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
			bfl_swap_bfl(lattice[k], lattice[k-1]);
			bfl_swap_bin(T[k], T[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}
	//
	// restore precision and
	// generate bigfloats
	//
	bigfloat::set_precision(old_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > old_bit_prec) {
				bit_diff = bi_bit_len-old_bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(value[i][j], tempmz0, bit_len-bit_diff);
		}

	//
	// Copy into math_matrix< bigint > Tr
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
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}



void bigfloat_lattice_gensys::Tr_lll_var_bfl_gensys(math_matrix< bigint > & Bi,
						    sdigit bit_len,
						    lidia_size_t& rank, sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_len, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool vector_is_zero;
	bool repeat_loop;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	sdigit short_prec;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	bigint tempbiz;
	bigint *tempvect1;
	bigint **BiAddr;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	BiAddr = Bi.get_data_address();

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	//
	// think again
	//
	mydel = new bigfloat[rows*(rows+columns+2)+columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_prec) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//
	tempvect1 = new bigint[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_prec) :: "
		       "not enough memory !");

	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];


	//
	// Start
	//
	// Copy into floating point
	//
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(lattice[i][j], tempmz0, bit_len-bit_diff);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(my[j][k], tempmz0, 2*bit_len-bit_diff);
					//                ::divide(my[j][k], tempbiz, SqrFactor);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);

					}
			}

			//
			// Correction Step II
			// lattice=(bigfloat )value
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					if (!(BiAddr[k][j].is_zero())) {
						bi_bit_len = BiAddr[k][j].bit_length();
						if (bi_bit_len > bit_prec) {
							bit_diff = bi_bit_len-bit_prec;
							shift_right(tempbiz, BiAddr[k][j], bit_diff);
							tempmz0.assign(tempbiz);
						}
						else {
							bit_diff = 0;
							tempmz0.assign(BiAddr[k][j]);
						}
						shift_right(lattice[k][j], tempmz0, bit_len-bit_diff);
					}
					else
						lattice[k][j].assign_zero();
			}

			//
			// Check if it`s 0 - vector
			// delete then
			//
			vector_is_zero = true;
			for  (lidia_size_t j = 0; j < columns; j++)
				if (!(BiAddr[k][j].is_zero())) {
					vector_is_zero = false;
					break;
				}

			if (vector_is_zero) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bfl_swap_bfl(lattice[j], lattice[j+1]);
					bfl_swap_bin(BiAddr[j], BiAddr[j+1]);
				}
				bfl_assign_zero_bin(BiAddr[--mom_rows]);
				bfl_assign_zero_bfl(lattice[mom_rows]);
				k = 1;
				bfl_assign_zero_bfl(vect);
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
				vector_is_zero = false;
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

		if (!break_flag) {
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_bfl_gensys(Bi, Factor, rank, bit_prec) [4]");
			//
			// Check lll - condition for parameter y_par
			//
			tempmz0.assign(y_par);
			LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
			LiDIA::square(tempmz1, my[k-1][k]);
			LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
			LiDIA::add(tempmz1, tempmz1, vect[k]);

			if (tempmz0.compare(tempmz1) > 0) {
				//
				// ::swap vectors k and k-1
				//
				++swapping_steps;
				bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
				bfl_swap_bfl(lattice[k], lattice[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
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
	// restore precision and
	// generate bigfloats
	//
	bigfloat::set_precision(old_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > old_bit_prec) {
				bit_diff = bi_bit_len-old_bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(value[i][j], tempmz0, bit_len-bit_diff);
		}

	//
	// Free allocated storage
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}



void bigfloat_lattice_gensys::Tr_lll_trans_var_bfl_gensys(math_matrix< bigint > & Bi,
							  math_matrix< bigint > &Tr,
							  sdigit bit_len,
							  lidia_size_t& rank, sdigit bit_prec)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, "
		      " Factor, bit_len, rank, bit_prec) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool vector_is_zero;
	bool repeat_loop;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	sdigit short_prec;
	lidia_size_t k = 1;
	lidia_size_t mom_rows = rows;
	bigint tempbiz;
	bigint *tempvect1;
	bigint **BiAddr;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat bound = std::exp(std::log(2.0)*bit_prec/2);
	bigfloat invbound;
	bigfloat *mydel;
	bigfloat **my;
	bigfloat **lattice;
	bigfloat *vect;
	bigfloat *vect_sq;
	bigfloat halb(0.5);

	LiDIA::divide(invbound, 1, bound);

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	short_prec = static_cast<sdigit>(static_cast<double>(bit_prec)*std::log(2.0)/std::log(10.0));
	BiAddr = Bi.get_data_address();

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new bigfloat*[rows*2];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, Factor, bit_prec) :: "
		       "not enough memory !");
	lattice = &my[rows];
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, Factor, bit_prec) :: "
		       "Not enough memory !");


	//
	// think again
	//
	mydel = new bigfloat[rows*(rows+columns+2)+columns];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, Factor, bit_prec) :: "
                       "not enough memory !");
	Tdel = new bigint[rows*rows+((rows > columns)?rows:columns)];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, Factor, bit_prec) :: "
		       "Not enough memory !");

	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		T[i] = &Tdel[i*rows];
		T[i][i].assign_one();
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	//
	// think again Part II
	//

	tempvect1 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	bigfloat::set_precision(short_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > bit_prec) {
				bit_diff = bi_bit_len-bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(lattice[i][j], tempmz0, bit_len-bit_diff);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	bfl_assign_zero_bfl(vect);
	bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Factor, bit_prec) [2]");
			bfl_scalprod_bfl(vect[k], lattice[k], lattice[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_bfl(tempmz0, lattice[k], lattice[j]);
				LiDIA::multiply(tempmz1, invbound, vect_sq[k]);
				LiDIA::multiply(tempmz1, tempmz1, vect_sq[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (abs(tempmz0).compare(tempmz1) < 0) {
					//
					// Compute the correct scalar product
					//
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > bit_prec) {
						bit_diff = bi_bit_len-bit_prec;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(my[j][k], tempmz0, 2*bit_len-bit_diff);
				}
				else
					my[j][k].assign(tempmz0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempmz1, my[i][j], my[i][k]);
					LiDIA::multiply(tempmz1, tempmz1, vect[i]);
					LiDIA::subtract(my[j][k], my[j][k], tempmz1);
				}

				tempmz0.assign(my[j][k]);
				LiDIA::divide(my[j][k], my[j][k], vect[j]);
				LiDIA::multiply(tempmz1, my[j][k], tempmz0);
				LiDIA::subtract(vect[k], vect[k], tempmz1);
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, Factor, bit_prec) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (abs(my[j][k]).compare(halb) > 0) {
						++reduction_steps;
						korr = true;
						round (tempmz0, my[j][k]);
						if (abs(tempmz0).compare(bound) > 0)
							Fc = true;

						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect1, tempbiz, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;

						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempmz1, tempmz0, my[i][j]);
							LiDIA::subtract(my[i][k], my[i][k], tempmz1);
						}
						LiDIA::subtract(my[j][k], my[j][k], tempmz0);

					}
			}

			//
			// Correction Step II
			// lattice=(bigfloat )value
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					if (!(BiAddr[k][j].is_zero())) {
						bi_bit_len = BiAddr[k][j].bit_length();
						if (bi_bit_len > bit_prec) {
							bit_diff = bi_bit_len-bit_prec;
							shift_right(tempbiz, BiAddr[k][j], bit_diff);
							tempmz0.assign(tempbiz);
						}
						else {
							bit_diff = 0;
							tempmz0.assign(BiAddr[k][j]);
						}
						shift_right(lattice[k][j], tempmz0, bit_len-bit_diff);
					}
					else
						lattice[k][j].assign_zero();
				bfl_scalquad_bfl(vect_sq[k], lattice[k]);
				sqrt(vect_sq[k], vect_sq[k]);
				korr = false;
			}
			//
			// Check if it`s 0 - vector
			// delete then
			//
			vector_is_zero = true;
			for  (lidia_size_t j = 0; j < columns; j++)
				if (!(BiAddr[k][j].is_zero())) {
					vector_is_zero = false;
					break;
				}

			if (vector_is_zero) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bfl_swap_bfl(lattice[j], lattice[j+1]);
					bfl_swap_bin(BiAddr[j], BiAddr[j+1]);
					bfl_swap_bin(T[j], T[j+1]);
				}
				bfl_assign_zero_bin(BiAddr[--mom_rows]);
				bfl_assign_zero_bfl(lattice[mom_rows]);
				vectsize = rows;
				bfl_assign_zero_bin(T[mom_rows]);
				vectsize = columns;
				k = 1;
				bfl_assign_zero_bfl(vect);
				bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
				sqrt(vect_sq[0], vect[0]);
				vector_is_zero = false;
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

		if (!break_flag) {
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_bfl_gensys(Bi, Tr, Factor, bit_prec) [4]");
			//
			// Check lll - condition for parameter y_par
			//
			tempmz0.assign(y_par);
			LiDIA::multiply(tempmz0, tempmz0, vect[k-1]);
			LiDIA::square(tempmz1, my[k-1][k]);
			LiDIA::multiply(tempmz1, tempmz1, vect[k-1]);
			LiDIA::add(tempmz1, tempmz1, vect[k]);

			if (tempmz0.compare(tempmz1) > 0) {
				//
				// ::swap vectors k and k-1
				//
				++swapping_steps;
				bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
				bfl_swap_bfl(lattice[k], lattice[k-1]);
				bfl_swap_bin(T[k], T[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_bfl(vect[0], lattice[0], lattice[0]);
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
	// restore precision and
	// generate bigfloats
	//
	bigfloat::set_precision(old_prec);
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > old_bit_prec) {
				bit_diff = bi_bit_len-old_bit_prec;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(value[i][j], tempmz0, bit_len-bit_diff);
		}

	//
	// Copy into math_matrix< bigint > Tr
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
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
