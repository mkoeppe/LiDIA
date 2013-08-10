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
//	Author	:Werner Backes (WB), Thorsten Lauer (TL)
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
// using doubles for approximation
//
//
// Schnorr - Euchner
// * Tr_lll_dbl()
// * Tr_lll_var_dbl();
// * Tr_lll_trans_dbl(T)
// * Tr_lll_trans_var_dbl(T)
// * Tr_lll_deep_insert_dbl()  (not yet implemented)
// * Tr_lll_dbl_gensys()
// * Tr_lll_var_dbl_gensys()
// * Tr_lll_trans_dbl_gensys(T)
// * Tr_lll_trans_var_dbl_gensys(T)
//
// Modified lll
// * Tr_mlll_dbl(v)  (not yet implemented, possible ?)
//


//
// Schnorr - Euchner version of lll optimized for bigfloats using doubles
// result : lll - reduced lattice for parameter y
//
void bigfloat_lattice_gensys::Tr_lll_dbl()
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl() [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	sdigit r_prec;
	sdigit n_prec;
	sdigit old_prec;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint x_factor;
	bigfloat conv;
	bigfloat *tempvect1;

	vectsize = columns;

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= columns; // columns >= rows, basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;
	bigfloat::set_precision(n_prec);

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_dbl() :: "
		       "not enough memory !");
	lattice = &my[rows];
	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_dbl() :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_dbl() :: "
		       "not enough memory !");
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			value[i][j].doublify(lattice[i][j]);

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl() [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Check Precision, correct if nessesary
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl() [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							my[i][k] -= tempdbl0*my[i][j];
						}

						my[j][k] -= tempdbl0;

					}
			}
			//
			// correction (not in paper)
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					value[k][j].doublify(lattice[k][j]);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl() [4]");
		if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
			++swapping_steps;
			bfl_swap_bfl(value[k], value[k-1]);
			bfl_swap_dbl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}

//
// Schnorr - Euchner version of lll optimized for bigints using doubles
// result : transformation lattice and reduced lattice
//
void bigfloat_lattice_gensys::Tr_lll_trans_dbl(math_matrix< bigint > & Tr)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	sdigit r_prec;
	sdigit n_prec;
	sdigit old_prec;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint conv;
	bigint x_factor;
	bigint *tempvect2;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat *tempvect1;

	vectsize = columns;

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= columns; // columns >= rows, basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;
	bigfloat::set_precision(n_prec);

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "not enough memory !");
	lattice = &my[rows];
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "Not enough memory !");

	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
                       "not enough memory !");
	Tdel = new bigint[rows*rows+rows];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
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

	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "not enough memory !");
	tempvect2 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			value[i][j].doublify(lattice[i][j]);

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Check Precision, correct if nessesary
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						tempmz0.bigintify(conv);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect2, conv, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect2);
						vectsize = columns;
						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];

						my[j][k] -= tempdbl0;
					}
			}

			//
			// correction (not in paper)
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					value[k][j].doublify(lattice[k][j]);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) [4]");
		if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
			++swapping_steps;
			bfl_swap_bfl(value[k], value[k-1]);
			bfl_swap_dbl(lattice[k], lattice[k-1]);
			bfl_swap_bin(T[k], T[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
				LiDIA::swap(TrAddr[j][i], T[i][j]);
	//
	// restore precision and
	// free allocated storage
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}

//
// Schnorr - Euchner version of lll optimized for bigfloats using doubles
// result : lll - reduced lattice for parameter y
//
void bigfloat_lattice_gensys::Tr_lll_dbl_gensys(lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool repeat_loop;
	bool vector_is_zero;
	lidia_size_t k;
	sdigit r_prec;
	sdigit n_prec;
	sdigit old_prec;
	lidia_size_t mom_rows = rows;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double u_zero = 0.0;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint x_factor;
	bigfloat u_zero_bfl;
	bigfloat conv;
	bigfloat *tempvect1;

	vectsize = columns;

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	if (r_prec < 200)
		u_zero = 1/std::exp(std::log(10.0)*r_prec);
	else {
		log(u_zero_bfl, 10);
		LiDIA::multiply(u_zero_bfl, u_zero_bfl, r_prec);
		exp(u_zero_bfl, u_zero_bfl);
		LiDIA::divide(u_zero_bfl, 1, u_zero_bfl);
	}
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= (rows > columns)?rows:columns; // columns > rows ? da event. keine basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;
	bigfloat::set_precision(n_prec);

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) :: "
		       "not enough memory !");
	lattice = &my[rows];
	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) :: "
		       "not enough memory !");

	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			value[i][j].doublify(lattice[i][j]);

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute scalar product correct
					//
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];

						my[j][k] -= tempdbl0;
					}
			}

			//
			// Correction Step II,
			// lattice[k]=(double )value[k]
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					value[k][j].doublify(lattice[k][j]);
			}
			//
			// Check if it`s 0 - vector
			// delete then
			//

			if (r_prec < 200) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[k]);
				tempdbl0 = std::sqrt(tempdbl0);
				if (std::fabs(tempdbl0) < u_zero)
					vector_is_zero = true;
				else
					vector_is_zero = false;
			}
			else {
				bigfloat::set_precision(r_prec);
				bfl_scalprod_bfl(tempmz0, value[k], value[k]);
				sqrt(tempmz0, tempmz0);
				if (abs(tempmz0) < u_zero_bfl)
					vector_is_zero = true;
				else
					vector_is_zero = false;
				bigfloat::set_precision(n_prec);
			}

			if (vector_is_zero) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bfl_swap_dbl(lattice[j], lattice[j+1]);
					bfl_swap_bfl(value[j], value[j+1]);
				}
				bfl_assign_zero_bfl(value[--mom_rows]);
				bfl_assign_zero_dbl(lattice[mom_rows]);
				k = 1;
				bfl_assign_zero_dbl(vect);
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
			//
			// Check lll - condition for parameter y_par
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_dbl_gensys(rank) [4]");
			if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
				//
				// ::swap vectors k and k-1
				//
				++swapping_steps;
				bfl_swap_bfl(value[k], value[k-1]);
				bfl_swap_dbl(lattice[k], lattice[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
					vect_sq[0] = std::sqrt(vect[0]);
				}
			}
			else
				++k;
		}
		break_flag = false;
	}
	rank = mom_rows;
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
// Schnorr - Euchner version of lll optimized for bigints using doubles
// result : transformation lattice and reduced lattice
//
void bigfloat_lattice_gensys::Tr_lll_trans_dbl_gensys(math_matrix< bigint > & Tr, lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) [1]");
	bool Fc;
	bool break_flag;
	bool korr;
	bool repeat_loop;
	bool vector_is_zero;
	sdigit r_prec;
	sdigit n_prec;
	sdigit old_prec;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double u_zero = 0.0;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint conv;
	bigint x_factor;
	bigint *tempvect2;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigfloat u_zero_bfl;
	bigfloat *tempvect1;

	vectsize = columns;

	//
	// Compute needed precision for exact arithmetik
	//
	old_prec = bigfloat::get_precision();
	r_prec = compute_read_precision();
	if (r_prec < 200)
		u_zero = 1/std::exp(std::log(10.0)*r_prec);
	else {
		log(u_zero_bfl, 10);
		LiDIA::multiply(u_zero_bfl, u_zero_bfl, r_prec);
		exp(u_zero_bfl, u_zero_bfl);
		LiDIA::divide(u_zero_bfl, 1, u_zero_bfl);
	}
	n_prec = static_cast<sdigit>(r_prec/std::log(10.0))+1;
	n_prec *= (rows > columns)?rows:columns; // columns > rows ? da event. keine basis !!!
	n_prec = static_cast<sdigit>(n_prec*std::log(10.0))+1;
	bigfloat::set_precision(n_prec);

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "not enough memory !");
	lattice = &my[rows];
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "Not enough memory !");

	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
                       "not enough memory !");
	Tdel = new bigint[rows*rows+rows];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
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

	tempvect1 = new bigfloat[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "not enough memory !");
	tempvect2 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			value[i][j].doublify(lattice[i][j]);

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute the correct scalar product
					//
					bfl_scalprod_bfl(tempmz0, value[k], value[j]);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						tempmz0.bigintify(conv);
						bfl_scalmul_bfl(tempvect1, tempmz0, value[j]);
						bfl_subtract_bfl(value[k], value[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect2, conv, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect2);
						vectsize = columns;
						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];

						my[j][k] -= tempdbl0;
					}
			}

			//
			// Corrections Step II
			// lattice[k]=(double )value[k]
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					value[k][j].doublify(lattice[k][j]);
			}
			//
			// Check if it`s 0 - vector
			// then delete it
			//
			if (r_prec < 200) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[k]);
				tempdbl0 = std::sqrt(tempdbl0);
				if (std::fabs(tempdbl0) < u_zero)
					vector_is_zero = true;
				else
					vector_is_zero = false;
			}
			else {
				bigfloat::set_precision(r_prec);
				bfl_scalprod_bfl(tempmz0, value[k], value[k]);
				sqrt(tempmz0, tempmz0);
				if (abs(tempmz0) < u_zero_bfl)
					vector_is_zero = true;
				else
					vector_is_zero = false;
				bigfloat::set_precision(n_prec);
			}

			if (vector_is_zero) {
				for (lidia_size_t j = 0; j < mom_rows-1; j++) {
					bfl_swap_bfl(value[j], value[j+1]);
					bfl_swap_dbl(lattice[j], lattice[j+1]);
					bfl_swap_bin(T[j], T[j+1]);
				}
				bfl_assign_zero_bfl(value[--mom_rows]);
				vectsize = rows;
				bfl_assign_zero_bin(T[mom_rows]);
				vectsize = columns;
				k = 1;
				bfl_assign_zero_dbl(vect);
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
			//
			// Check lll - condition for parameter y_par
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) [4]");
			if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
				//
				// swap vectors k and k-1
				//
				++swapping_steps;
				bfl_swap_bfl(value[k], value[k-1]);
				bfl_swap_dbl(lattice[k], lattice[k-1]);
				bfl_swap_bin(T[k], T[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
					vect_sq[0] = std::sqrt(vect[0]);
				}
			}
			else
				++k;
		}
		break_flag = false;
	}
	rank = mom_rows;
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
	// Free allocated storage and
	// restore precision
	//
	bigfloat::set_precision(old_prec);
	delete[] tempvect1;
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}


//
// Schnorr - Euchner version of lll optimized for bigfloats using doubles
// result : lll - reduced lattice for parameter y
//
void bigfloat_lattice_gensys::Tr_lll_var_dbl(math_matrix< bigint > & Bi, sdigit bit_len)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi, bit_len) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	sdigit old_prec;
	sdigit bit_diff = 0;
	sdigit bi_bit_len;
	sdigit old_bit_prec;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint tempbiz;
	bigint **BiAddr;
	bigint *tempvect1;

	vectsize = columns;
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	bigfloat::set_precision(20);
	BiAddr = Bi.get_data_address();
	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi) :: "
		       "not enough memory !");
	lattice = &my[rows];
	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	tempvect1 = new bigint[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi) :: "
		       "not enough memory !");
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > 8*SIZEOF_DOUBLE) {
				bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(tempmz0, tempmz0, bit_len-bit_diff);
			tempmz0.doublify(lattice[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi, Factor) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Check Precision, correct if nessesary
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > 8*SIZEOF_DOUBLE) {
						bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(tempmz0, tempmz0, 2*bit_len-bit_diff);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//

			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi, Factor) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];

						my[j][k] -= tempdbl0;
					}
			}

			//
			// correction (not in paper)
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = BiAddr[k][j].bit_length();
					if (bi_bit_len > 8*SIZEOF_DOUBLE) {
						bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
						shift_right(tempbiz, BiAddr[k][j], bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(BiAddr[k][j]);
					}
					shift_right(tempmz0, tempmz0, bit_len-bit_diff);
					tempmz0.doublify(lattice[k][j]);
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


		debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl(Bi, Factor) [4]");
		if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
			++swapping_steps;
			bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
			bfl_swap_dbl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
	// and restore precision
	//
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}


//
// Schnorr - Euchner version of lll optimized for bigfloats using doubles
// result : lll - reduced lattice for parameter y
//
void bigfloat_lattice_gensys::Tr_lll_trans_var_dbl(math_matrix< bigint > & Bi,
						   math_matrix< bigint > & Tr,
						   sdigit bit_len)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl(Bi, Tr, Factor, bit_len) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint tempbiz;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigint **BiAddr;
	bigint *tempvect1;

	vectsize = columns;
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	bigfloat::set_precision(20);
	BiAddr = Bi.get_data_address();
	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl(Bi, Tr, Factor) :: "
		       "not enough memory !");
	lattice = &my[rows];
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "Not enough memory !");

	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl(Bi, Tr, Factor) :: "
                       "not enough memory !");
	Tdel = new bigint[rows*rows+columns];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
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

	tempvect1 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > 8*SIZEOF_DOUBLE) {
				bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(tempmz0, tempmz0, bit_len-bit_diff);
			tempmz0.doublify(lattice[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	Fc = false;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl(Bi, Tr, Factor) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Check Precision, correct if nessesary
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > 8*SIZEOF_DOUBLE) {
						bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(tempmz0, tempmz0, 2*bit_len-bit_diff);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl(Bi, Tr, Factor) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect1, tempbiz, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;

						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];

						my[j][k] -= tempdbl0;
					}
			}

			//
			// correction (not in paper)
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					bi_bit_len = BiAddr[k][j].bit_length();
					if (bi_bit_len > 8*SIZEOF_DOUBLE) {
						bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
						shift_right(tempbiz, BiAddr[k][j], bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(BiAddr[k][j]);
					}
					shift_right(tempmz0, tempmz0, bit_len-bit_diff);
					tempmz0.doublify(lattice[k][j]);
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

		debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl(Bi, Tr, Factor) [4]");
		if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
			++swapping_steps;
			bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
			bfl_swap_dbl(lattice[k], lattice[k-1]);
			bfl_swap_bin(T[k], T[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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

//
// Schnorr - Euchner version of lll optimized for bigfloats using doubles
// result : lll - reduced lattice for parameter y
//
void bigfloat_lattice_gensys::Tr_lll_var_dbl_gensys(math_matrix< bigint > & Bi,
						    sdigit bit_len,
						    lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor, bit_len, rank) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool vector_is_zero;
	bool repeat_loop;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint tempbiz;
	bigint **BiAddr;
	bigint *tempvect1;

	vectsize = columns;
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	bigfloat::set_precision(20);
	BiAddr = Bi.get_data_address();
	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor) :: "
		       "not enough memory !");
	lattice = &my[rows];
	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor) :: "
                       "not enough memory !");
	for (lidia_size_t i = 0; i < rows; i++) {
		my[i] = &mydel[i*rows]; // rows x rows
		lattice[i] = &mydel[(rows*rows)+(i*columns)]; // rows x columns
	}

	//
	// Allocating Memory for vector
	//

	tempvect1 = new bigint[columns];
	memory_handler(tempvect1, "bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor) :: "
		       "not enough memory !");
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > 8*SIZEOF_DOUBLE) {
				bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(tempmz0, tempmz0, bit_len-bit_diff);
			tempmz0.doublify(lattice[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	vector_is_zero = false;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Check Precision, correct if nessesary
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > 8*SIZEOF_DOUBLE) {
						bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(tempmz0, tempmz0, 2*bit_len-bit_diff);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							my[i][k] -= tempdbl0*my[i][j];
						}

						my[j][k] -= tempdbl0;

					}
			}

			//
			// correction (not in paper)
			//

			if (korr) {
				++correction_steps;
				for (lidia_size_t j = 0; j < columns; j++) {
					if (!(BiAddr[k][j].is_zero())) {
						bi_bit_len = BiAddr[k][j].bit_length();
						if (bi_bit_len > 8*SIZEOF_DOUBLE) {
							bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
							shift_right(tempbiz, BiAddr[k][j], bit_diff);
							tempmz0.assign(tempbiz);
						}
						else {
							bit_diff = 0;
							tempmz0.assign(BiAddr[k][j]);
						}
						shift_right(tempmz0, tempmz0, bit_len-bit_diff);
						tempmz0.doublify(lattice[k][j]);
					}
					else
						lattice[k][j] = 0.0;
				}
				bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
				vect_sq[k] = std::sqrt(vect[k]);
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
					bfl_swap_dbl(lattice[j], lattice[j+1]);
					bfl_swap_bin(BiAddr[j], BiAddr[j+1]);
				}
				bfl_assign_zero_bin(BiAddr[--mom_rows]);
				bfl_assign_zero_dbl(lattice[mom_rows]);
				k = 1;
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_var_dbl_gensys(Bi, Factor) [4]");
			if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
				++swapping_steps;
				bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
				bfl_swap_dbl(lattice[k], lattice[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
					vect_sq[0] = std::sqrt(vect[0]);
				}
			}
			else
				++k;
		}
		break_flag = false;
	}
	rank = mom_rows;
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
	// and restore precision
	//
	delete[] tempvect1;
	delete[] mydel;
	delete[] my;
}


//
// Schnorr - Euchner version of lll optimized for bigfloats using doubles
// result : lll - reduced lattice for parameter y
//
void bigfloat_lattice_gensys::Tr_lll_trans_var_dbl_gensys(math_matrix< bigint > & Bi,
							  math_matrix< bigint > & Tr,
							  sdigit bit_len,
							  lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, bit_len, rank) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool vector_is_zero;
	bool repeat_loop;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	sdigit old_prec;
	sdigit old_bit_prec;
	sdigit bit_diff;
	sdigit bi_bit_len;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigint tempbiz;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;
	bigint **BiAddr;
	bigint *tempvect1;

	vectsize = columns;
	old_prec = bigfloat::get_precision();
	old_bit_prec = static_cast<sdigit>(static_cast<double>(old_prec)/std::log(2.0)*std::log(10.0)+1);
	bigfloat::set_precision(20);
	BiAddr = Bi.get_data_address();
	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) :: "
		       "not enough memory !");
	lattice = &my[rows];
	T = new bigint*[rows];
	memory_handler(T, "bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) :: "
		       "Not enough memory !");

	mydel = new double[rows*(rows+columns+2)];
	memory_handler(mydel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) :: "
                       "not enough memory !");
	Tdel = new bigint[rows*rows+((rows > columns)?rows:columns)];
	memory_handler(Tdel, "bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) :: "
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

	tempvect1 = &Tdel[rows*rows];
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++) {
			bi_bit_len = BiAddr[i][j].bit_length();
			if (bi_bit_len > 8*SIZEOF_DOUBLE) {
				bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
				shift_right(tempbiz, BiAddr[i][j], bit_diff);
				tempmz0.assign(tempbiz);
			}
			else {
				bit_diff = 0;
				tempmz0.assign(BiAddr[i][j]);
			}
			shift_right(tempmz0, tempmz0, bit_len-bit_diff);
			tempmz0.doublify(lattice[i][j]);
		}

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	break_flag = false;
	vector_is_zero = false;
	bfl_assign_zero_dbl(vect);
	bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) [2]");
			bfl_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bfl_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Check Precision, correct if nessesary
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					bfl_scalprod_bin(tempbiz, BiAddr[k], BiAddr[j]);
					bi_bit_len = tempbiz.bit_length();
					if (bi_bit_len > 8*SIZEOF_DOUBLE) {
						bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
						shift_right(tempbiz, tempbiz, bit_diff);
						tempmz0.assign(tempbiz);
					}
					else {
						bit_diff = 0;
						tempmz0.assign(tempbiz);
					}
					shift_right(tempmz0, tempmz0, 2*bit_len-bit_diff);
					tempmz0.doublify(my[j][k]);
				}
				else
					my[j][k] = tempdbl0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[j][k] -= my[i][j]*my[i][k]*vect[i];

				tempdbl0 = my[j][k];
				my[j][k] /= vect[j];
				vect[k] -= my[j][k]*tempdbl0;
			}
			//
			// end of orthogonalization
			//
			Fc = false;
			korr = false;
			//
			// begin of reduction
			//
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						tempmz0.assign(tempdbl0);
						tempmz0.bigintify(tempbiz);
						bfl_scalmul_bin(tempvect1, tempbiz, BiAddr[j]);
						bfl_subtract_bin(BiAddr[k], BiAddr[k], tempvect1);
						vectsize = rows;
						bfl_scalmul_bin(tempvect1, tempbiz, T[j]);
						bfl_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;

						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];

						my[j][k] -= tempdbl0;
					}
			}

			//
			// correction (not in paper)
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++) {
					if (!(BiAddr[k][j].is_zero())) {
						bi_bit_len = BiAddr[k][j].bit_length();
						if (bi_bit_len > 8*SIZEOF_DOUBLE) {
							bit_diff = bi_bit_len-8*SIZEOF_DOUBLE;
							shift_right(tempbiz, BiAddr[k][j], bit_diff);
							tempmz0.assign(tempbiz);
						}
						else {
							bit_diff = 0;
							tempmz0.assign(BiAddr[k][j]);
						}
						shift_right(tempmz0, tempmz0, bit_len-bit_diff);
						tempmz0.doublify(lattice[k][j]);
					}
					else
						lattice[k][j] = 0.0;
				}
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
					bfl_swap_dbl(lattice[j], lattice[j+1]);
					bfl_swap_bin(BiAddr[j], BiAddr[j+1]);
					bfl_swap_bin(T[j], T[j+1]);
				}
				bfl_assign_zero_bin(BiAddr[--mom_rows]);
				bfl_assign_zero_dbl(lattice[mom_rows]);
				vectsize = rows;
				bfl_assign_zero_bin(T[mom_rows]);
				vectsize = columns;
				k = 1;
				bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[k] = std::sqrt(vect[0]);
				vector_is_zero = false;
				break_flag = true;
				Fc = false;
				korr = false;
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
			debug_handler("bigfloat_lattice_gensys", "Tr_lll_trans_var_dbl_gensys(Bi, Tr, Factor, rank) [4]");
			if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
				++swapping_steps;
				bfl_swap_bin(BiAddr[k], BiAddr[k-1]);
				bfl_swap_dbl(lattice[k], lattice[k-1]);
				bfl_swap_bin(T[k], T[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bfl_scalprod_dbl(vect[0], lattice[0], lattice[0]);
					vect_sq[0] = std::sqrt(vect[0]);
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
