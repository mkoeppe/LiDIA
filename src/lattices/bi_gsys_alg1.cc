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
// using `double` approximation
//
//
// Schnorr - Euchner
// * Tr_lll_dbl()
// * Tr_lll_trans_dbl(T)
// * Tr_lll_deep_insert_dbl() (not working yet)
// * Tr_lll_dbl_gensys()
// * Tr_lll_trans_dbl_gensys(T)
//
// Modified lll
// * Tr_mlll_dbl(v) (not yet implemented, possible ?)
//


//
// Schnorr - Euchner version of lll optimized for bigints using doubles
// result : lll - reduced lattice for parameter y
//
void bigint_lattice_gensys::Tr_lll_dbl()
{
	debug_handler("bigint_lattice_gensys", "Tr_lll_dbl() [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigfloat conv;
	bigint *tempvect1;

	vectsize = columns;
	//
	// Allocating Memory for matrix
	//
	// my[rows][rows], lattice[rows][columns], vect[rows], vect_sq[rows]
	//
	my = new double*[rows*2];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_dbl() :: "
		       "not enough memory !");
	mydel = new double[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_dbl() :: "
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
	memory_handler(tempvect1, "bigint_lattice_gensys", "Tr_lll_dbl() :: "
		       "not enough memory !");
	//
	// Setting pointers
	//
	vect = &mydel[(rows+columns)*rows];
	vect_sq = &vect[rows];

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			lattice[i][j] = dbl(value[i][j]);

	k = 1;
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	vectsize = columns;
	Fc = false;
	bin_assign_zero_dbl(vect);
	bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_dbl() [2]");
			bin_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Correction Step I, only  if nessesary |lattice[k]*lattice[j]| to big
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					my[j][k] = dbl(tempmz0);
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
			debug_handler("bigint_lattice_gensys", "Tr_lll_dbl() [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						conv.assign(tempdbl0);
						conv.bigintify(tempmz0);

						//
						// if tempdbl0 is small enough       ? kommt noch
						// perform conversion over sdigit`s
						//
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i)
							my[i][k] -= tempdbl0*my[i][j];
						my[j][k] -= tempdbl0;
					}
			}

			//
			// Correction Step II
			// needed after Reduction, correct only reduced vector
			// lattice[k] = (double ) value[k];
			//
			if (korr) {
				correction_steps++;
				for  (lidia_size_t j = 0; j < columns; j++)
					lattice[k][j] = dbl(value[k][j]);
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

		debug_handler("bigint_lattice_gensys", "Tr_lll_dbl() [4]");
		//
		// Check lll - conditions for parameter y_par
		//
		if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
			//
			// swap vectors at position k and k-1
			//
			swapping_steps++;
			bin_swap_bin(value[k], value[k-1]);
			bin_swap_dbl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;

			if (k == 1) {
				bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
			}

		}
		else
			++k;
	}
	//
	// Free allocated storage
	//
	delete[] mydel;
	delete[] my;
	delete[] tempvect1;
}



//
// Schnorr - Euchner version of lll optimized for bigints usign doubles
// result : lll - reduced lattice for parameter y
//
void bigint_lattice_gensys::Tr_lll_trans_dbl(math_matrix< bigint > &Tr)
{
	debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) [1]");
	bool Fc;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigfloat conv;
	bigint *tempvect1;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "not enough memory !");
	mydel = new double[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
                       "not enough memory !");
	T = new bigint*[rows];
	memory_handler(T, "bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
		       "not enough memory ! ");
	Tdel = new bigint[rows*rows+((rows > columns)?rows:columns)];
	memory_handler(Tdel, "bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) :: "
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
	// (-> Setting vectors)
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
		for (lidia_size_t j = 0; j < columns; j++)
			lattice[i][j] = dbl(value[i][j]);
	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	vectsize = columns;
	bin_assign_zero_dbl(vect);
	bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) [2]");
			bin_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Correction Step I, correct only if needed, if |lattice[k]*lattice[j]| to big
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute scalar product correct
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					my[j][k] = dbl(tempmz0);
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
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						conv.assign(tempdbl0);
						conv.bigintify(tempmz0);
						//
						// Perform operation on lattice and transformation lattice
						//
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);
						vectsize = rows;
						bin_scalmul_bin(tempvect1, tempmz0, T[j]);
						bin_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;
						for (lidia_size_t i = 0; i < j; ++i) {
							my[i][k] -= tempdbl0*my[i][j];
						}

						my[j][k] -= tempdbl0;


					}
			}

			//
			// Correction Step II
			// lattice[k] = (double )value[k];
			//
			if (korr) {
				correction_steps++;
				for  (lidia_size_t j = 0; j < columns; j++)
					lattice[k][j] = dbl(value[k][j]);
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

		debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl(Tr) [4]");
		//
		// Check lll - condition for parameter y_par
		//
		if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
			//
			// Swap vectors k and k-1, don`t forget the transformation lattice
			//
			swapping_steps++;
			bin_swap_bin(value[k], value[k-1]);
			bin_swap_bin(T[k], T[k-1]);
			bin_swap_dbl(lattice[k], lattice[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);
				vect_sq[0] = std::sqrt(vect[0]);
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
	//
	delete[] mydel;
	delete[] Tdel;
	delete[] my;
	delete[] T;
}



//
// Schnorr - Euchner version of lll for linear generating systems
//
void bigint_lattice_gensys::Tr_lll_dbl_gensys(lidia_size_t& rank)
{
	debug_handler("bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) [1]");
	bool Fc;
	bool korr;
	bool break_flag;
	bool repeat_loop;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigfloat conv;
	bigint *tempvect1;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) :: "
		       "not enough memory !");
	mydel = new double[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) :: "
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
	memory_handler(tempvect1, "bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) :: "
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
			lattice[i][j] = dbl(value[i][j]);

	k = 1;
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
	break_flag = false;
	vectsize = columns;
	bin_assign_zero_dbl(vect);
	bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) [2]");
			bin_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					my[j][k] = dbl(tempmz0);
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
			debug_handler("bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						conv.assign(tempdbl0);
						conv.bigintify(tempmz0);
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);

						for (lidia_size_t i = 0; i < j; ++i) {
							my[i][k] -= tempdbl0*my[i][j];
						}

						my[j][k] -= tempdbl0;
					}

			}
			//
			// Correction Step II
			// lattice[k]=(double )value[k]
			//
			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					lattice[k][j] = dbl(value[k][j]);
				bin_scalquad_dbl(vect[k], lattice[k]);
				vect_sq[k] = std::sqrt(vect[k]);
			}
			//
			// Check if 0 - Vector
			// delete vector k if it`s
			//
			if (vect_sq[k] < 0.5) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bin_swap_dbl(lattice[j], lattice[j+1]);
					bin_swap_bin(value[j], value[j+1]);
				}
				bin_assign_zero_bin(value[--mom_rows]);
				k = 1;
				bin_assign_zero_dbl(vect);
				bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);
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
			debug_handler("bigint_lattice_gensys", "Tr_lll_dbl_gensys(rank) [4]");
			//
			// Check lll - condition for parameter y_par
			//
			if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
				//
				// Swap vectors k and k-1
				//
				++swapping_steps;
				bin_swap_bin(value[k], value[k-1]);
				bin_swap_dbl(lattice[k], lattice[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);
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
	// Free allocated storage
	//
	delete[] mydel;
	delete[] my;
	delete[] tempvect1;
}



//
// Schnorr - Euchner version of lll for linear generating systems
//
void bigint_lattice_gensys::Tr_lll_trans_dbl_gensys(math_matrix< bigint > & Tr, lidia_size_t& rank)
{
	debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) [1]");
	bool Fc;
	bool break_flag;
	bool korr;
	bool repeat_loop;
	lidia_size_t k;
	lidia_size_t mom_rows = rows;
	double bound = std::exp(std::log(2.0)*52/2);
	double invbound = 1/bound;
	double tempdbl0;
	double *mydel;
	double **my;
	double **lattice;
	double *vect;
	double *vect_sq;
	bigfloat conv;
	bigint *tempvect1;
	bigint *Tdel;
	bigint **T;
	bigint **TrAddr;

	vectsize = columns;

	//
	// Allocating Memory for matrix
	//
	my = new double*[2*rows];
	memory_handler(my, "bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) :: "
		       "not enough memory !");
	mydel = new double[rows*((rows+columns)+2)];
	memory_handler(mydel, "bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) :: "
                       "not enough memory !");
	T = new bigint*[rows];
	memory_handler(T, "bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) :: "
		       "not enough memory ! ");
	Tdel = new bigint[rows*rows+((columns > rows)?columns:rows)];
	memory_handler(Tdel, "bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) :: "
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
	tempvect1 = &T[0][rows*rows];
	vect = &my[0][(rows+columns)*rows];
	vect_sq = &vect[rows];

	k = 1;

	//
	// Start
	//
	// Copy into floating point
	//
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			lattice[i][j] = dbl(value[i][j]);

	k = 1;
	reduction_steps = 0;
	correction_steps = 0;
	swapping_steps = 0;
	break_flag = false;
	vectsize = columns;
	bin_assign_zero_dbl(vect);
	bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);

	while (k < mom_rows) {
		Fc = true;
		repeat_loop = false;
		while (Fc) {
			//
			// begin of orthogonalization
			//
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) [2]");
			bin_scalprod_dbl(vect[k], lattice[k], lattice[k]);
			vect_sq[k] = std::sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				bin_scalprod_dbl(tempdbl0, lattice[k], lattice[j]);
				//
				// Correction Step I, only if needed, |lattice[k]*lattice[j]| to big
				//
				if (std::fabs(tempdbl0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute correct scalar product
					//
					bin_scalprod_bin(tempmz0, value[k], value[j]);
					my[j][k] = dbl(tempmz0);
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
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) [3]");
			if (!repeat_loop) {
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (std::fabs(my[j][k]) > 0.5) {
						++reduction_steps;
						korr = true;
						tempdbl0 = rint(my[j][k]);
						if (std::fabs(tempdbl0) > bound)
							Fc = true;

						conv.assign(tempdbl0);
						conv.bigintify(tempmz0);
						bin_scalmul_bin(tempvect1, tempmz0, value[j]);
						bin_subtract_bin(value[k], value[k], tempvect1);
						vectsize = rows;
						bin_scalmul_bin(tempvect1, tempmz0, T[j]);
						bin_subtract_bin(T[k], T[k], tempvect1);
						vectsize = columns;

						for (lidia_size_t i = 0; i < j; ++i) {
							my[i][k] -= tempdbl0*my[i][j];
						}

						my[j][k] -= tempdbl0;
					}
			}

			//
			// Correction Step II
			// lattice[k]=(double )value[k]
			//

			if (korr) {
				++correction_steps;
				for  (lidia_size_t j = 0; j < columns; j++)
					lattice[k][j] = dbl(value[k][j]);
				bin_scalquad_dbl(vect[k], lattice[k]);
				vect_sq[k] = std::sqrt(vect[k]);
			}
			//
			// Check if 0 - Vector
			// delete vector k if it`s
			//
			if (vect_sq[k] < 0.5) {
				for (lidia_size_t j = k; j < mom_rows-1; j++) {
					bin_swap_dbl(lattice[j], lattice[j+1]);
					bin_swap_bin(T[j], T[j+1]);
					bin_swap_bin(value[j], value[j+1]);
				}
				bin_assign_zero_bin(value[--mom_rows]);
				vectsize = rows;
				bin_assign_zero_bin(T[mom_rows]);
				vectsize = columns;
				k = 1;
				bin_assign_zero_dbl(vect);
				bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);
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
			debug_handler("bigint_lattice_gensys", "Tr_lll_trans_dbl_gensys(Tr, rank) [4]");
			//
			// Check lll - condition for parameter y_par
			//
			if (y_par*vect[k-1] > my[k-1][k]*my[k-1][k]*vect[k-1]+vect[k]) {
				//
				// Swap vectors k and k-1
				//
				++swapping_steps;
				bin_swap_bin(value[k], value[k-1]);
				bin_swap_bin(T[k], T[k-1]);
				bin_swap_dbl(lattice[k], lattice[k-1]);

				if (k > 1)
					--k;
				if (k == 1) {
					bin_scalprod_dbl(vect[0], lattice[0], lattice[0]);
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
	//
	delete[] Tdel;
	delete[] T;
	delete[] mydel;
	delete[] my;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
