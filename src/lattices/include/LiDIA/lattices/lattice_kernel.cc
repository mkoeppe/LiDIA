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


#ifndef LIDIA_LATTICE_KERNEL_CC_GUARD_
#define LIDIA_LATTICE_KERNEL_CC_GUARD_


#ifndef LIDIA_LATTICE_KERNEL_H_GUARD_
# include	"LiDIA/lattices/lattice_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif

#ifdef LIDIA_NAMESPACE
using std::sqrt;
using std::floor;
using std::fabs;
#endif



//
// Schnorr-Euchner-lll
// reducing lattice of type E, approximating with type A
//
template< class E, class A, class VectorOp, class AlgParts >
void
lll_kernel_op< E, A, VectorOp, AlgParts >::lll(dense_alg< E > & da,
					       lattice_info& li)
{
	debug_handler("lll_kernel_op", "lll(da, li) [1]");
	bool bigmy;
	bool correct;
	bool repeat_loop;
	lidia_size_t k;
	A bound = std::exp(std::log(2.0)*da.d.bit_prec/2);
	A invbound = 1/bound;
	A tempA0;
	A *mydel;
	A **my;
	A **valueA;
	A *vect;
	A *vect_sq;
	E tempE0;

	//
	// Allocating Memory for matrix
	//
	my = new A*[da.b.rows*2];
	memory_handler(my, "lll_kernel_op", "lll(da, li) :: "
		       "not enough memory !");
	mydel = new A[da.b.rows*((da.b.rows+da.b.columns)+2)];
	memory_handler(mydel, "lll_kernel_op", "lll(da, li) :: "
                       "not enough memory !");
	valueA = &my[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		my[i] = &mydel[i*da.b.rows]; // rows x rows
		valueA[i] = &mydel[(da.b.rows*da.b.rows)+(i*da.b.columns)]; // rows x columns
	}

	//
	// Setting pointers for A vectors
	//
	vect = &mydel[(da.b.rows+da.b.columns)*da.b.rows];
	vect_sq = &vect[da.b.rows];

	//
	// Start
	//
	// Copy into floating point
	//
	alg.prec_startup(da);
	alg.E_convert_lattice_A(da, valueA);
	alg.prec_exact(da);

	//
	// setting startup values
	//
	li.lll.y = da.b.y;
	li.lll.reduction_steps = 0;
	li.lll.correction_steps = 0;
	li.lll.swaps = 0;
	bigmy = false;
	vector.approx.vectsize = 2*da.b.rows;
	vector.approx.assign_zero(vect);
	vector.approx.vectsize = da.b.columns;
	vector.approx.scalprod(vect[0], valueA[0], valueA[0]);
	vector.exact.vectsize = da.b.columns;

	k = 1;
	while (k < da.b.rank) {
		bigmy = true;
		repeat_loop = false;
		while (bigmy) {
			//
			// begin of orthogonalization
			//
			debug_handler("lll_kernel_op", "lll(da, li) [2]");
			vector.approx.scalprod(vect[k], valueA[k], valueA[k]);
			vect_sq[k] = sqrt(vect[k]);

			for (lidia_size_t j = 0; j < k; ++j) {
				vector.approx.scalprod(tempA0, valueA[k], valueA[j]);
				//
				// Correction Step I, only  if nessesary |valueA[k]*valueA[j]| to small
				//
				if (fabs(tempA0) < invbound*vect_sq[k]*vect_sq[j]) {
					//
					// Compute correct scalar product
					//
					alg.prec_correct(da);
					vector.exact.scalprod(tempE0, da.s.value[k],
							      da.s.value[j]);
					da.d.bit_factor <<= 1;
					alg.E_convert_value_A(da, my[k][j], tempE0);
					da.d.bit_factor >>= 1;
					alg.prec_exact(da);
				}
				else
					my[k][j] = tempA0;

				for (lidia_size_t i = 0; i < j; ++i)
					my[k][j] -= my[j][i]*my[k][i]*vect[i];

				tempA0 = my[k][j];
				my[k][j] /= vect[j];
				vect[k] -= my[k][j]*tempA0;
			}
			//
			// end of orthogonalization
			//
			bigmy = false;
			correct = false;
			//
			// begin of reduction
			//
			debug_handler("lll_kernel_op", "lll(da, li) [3]");
			if (!repeat_loop) {
				//
				// trick to allow generating a transformation matrix
				// with the same algorithm
				//
				vector.exact.vectsize = da.b.real_columns;
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (fabs(my[k][j]) > 0.5) {
						++li.lll.reduction_steps;
						correct = true;
						tempA0 = lidia_rint(my[k][j]);
						if (alg.A_convert_value_E_bound(da, tempE0, tempA0, bound))
							bigmy = true;
						vector.exact.scalsub(da.s.value[k], da.s.value[k],
								     tempE0, da.s.value[j]);
						for (lidia_size_t i = 0; i < j; ++i)
							my[k][i] -= tempA0*my[j][i];
						my[k][j] -= tempA0;
					}
				vector.exact.vectsize = da.b.columns;
			}

			//
			// Correction Step II
			// needed after Reduction, correct only reduced vector
			// valueA[k] = (A) value[k];
			//
			if (correct) {
				++li.lll.correction_steps;
				if (alg.E_convert_vector_A(da, valueA, k)) {
					//
					// Found 0 - Vector ->repeat with k = 1;
					//
					bigmy = true;
					repeat_loop = false;
					k = 1;
					vector.approx.scalprod(vect[0], valueA[0], valueA[0]);
					vect_sq[0] = sqrt(vect[0]);
				}
			}

			//
			// end of reduction
			//
			if (bigmy) {
				//
				// Go back if my to big
				//
				if (k > 1)
					--k;
			}
			else {
				//
				// repeat only once
				//
				bigmy = (!repeat_loop) && correct;
				repeat_loop = !repeat_loop;
			}
		}


		debug_handler("lll_kernel_op", "lll(da, li) [4]");
		//
		// Check lll - conditions for parameter y_par
		//
		if (da.b.y*vect[k-1] > my[k][k-1]*my[k][k-1]*vect[k-1]+vect[k]) {
			//
			// swap vectors at position k and k-1
			//
			++li.lll.swaps;
			vector.exact.swap(da.s.value[k], da.s.value[k-1]);
			vector.approx.swap(valueA[k], valueA[k-1]);

			if (k > 1)
				--k;

			if (k == 1) {
				vector.approx.scalprod(vect[0], valueA[0], valueA[0]);
				vect_sq[0] = sqrt(vect[0]);
			}
		}
		else
			++k;
	}

	li.lll.rank = da.b.rank;
	alg.prec_restore(da);
	//
	// Free allocated storage
	//
	delete[] mydel;
	delete[] my;
}



//
// Schnorr-Euchner-lll
// modified version working with gram matrices
// reducing lattice of type E, approximating with type A
//
template< class E, class A, class VectorOp, class AlgParts >
void
lll_kernel_op< E, A, VectorOp, AlgParts >::lll_gram(dense_alg< E > & da,
						    lattice_info& li)
{
	debug_handler("lll_kernel_op", "lll_gram(da, li) [1]");
	bool bigmy;
	bool correct;
	bool repeat_loop;
	lidia_size_t k;
	A bound = std::exp(std::log(2.0)*da.d.bit_prec/2);
	A tempA0;
	A *mydel;
	A **my;
	A **valueA;
	A *vect;
	A *vect_sq;
	E tempE0;
	E tempE1;
	E tempE2;

	//
	// Allocating Memory for matrix
	//
	my = new A*[da.b.rows*2];
	memory_handler(my, "lll_kernel_op", "lll_gram(da, li) :: "
		       "not enough memory !");
	mydel = new A[da.b.rows*((da.b.rows+da.b.columns)+2)];
	memory_handler(mydel, "lll_kernel_op", "lll_gram(da, li) :: "
                       "not enough memory !");
	valueA = &my[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		my[i] = &mydel[i*da.b.rows]; // rows x rows
		valueA[i] = &mydel[(da.b.rows*da.b.rows)+(i*da.b.columns)]; // rows x columns
	}

	//
	// Setting pointers for A vectors
	//
	vect = &mydel[(da.b.rows+da.b.columns)*da.b.rows];
	vect_sq = &vect[da.b.rows];

	//
	// Start
	//
	// Copy into floating point
	//
	alg.prec_startup(da);
	alg.E_convert_lattice_A(da, valueA);
	if (da.b.rank < da.b.rows) {
		for (lidia_size_t i = 0; i < da.b.rows; i++)
			for (lidia_size_t j = i; j < da.b.columns; j++)
				LiDIA::swap(da.s.value[i][j], da.s.value[j][i]);
		da.b.rank = da.b.columns;
		alg.E_convert_lattice_A(da, valueA);
	}
	alg.prec_exact(da);

	//
	// setting startup values
	//
	li.lll.y = da.b.y;
	li.lll.reduction_steps = 0;
	li.lll.correction_steps = 0;
	li.lll.swaps = 0;
	bigmy = false;
	vector.approx.vectsize = 2*da.b.rows;
	vector.approx.assign_zero(vect);
	vect[0] = valueA[0][0];
	vector.approx.vectsize = da.b.columns;
	vector.exact.vectsize = da.b.columns;

	k = 1;
	while (k < da.b.rank) {
		bigmy = true;
		repeat_loop = false;
		while (bigmy) {
			//
			// begin of orthogonalization
			//
			debug_handler("lll_kernel_op", "lll_gram(da, li) [2]");
			vect[k] = valueA[k][k];
			vect_sq[k] = sqrt(vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				my[k][j] = valueA[k][j];
				for (lidia_size_t i = 0; i < j; ++i)
					my[k][j] -= my[j][i]*my[k][i]*vect[i];
				tempA0 = my[k][j];
				my[k][j] /= vect[j];
				vect[k] -= my[k][j]*tempA0;
			}
			//
			// end of orthogonalization
			//
			bigmy = false;
			correct = false;
			//
			// begin of reduction
			//
			debug_handler("lll_kernel_op", "lll_gram(da, li) [3]");
			if (!repeat_loop) {
				//
				// trick to allow generating a transformation matrix
				// with the same algorithm
				//
				vector.exact.vectsize = da.b.real_columns;
				for (lidia_size_t j = k-1; j >= 0; j--)
					if (fabs(my[k][j]) > 0.5) {
						++li.lll.reduction_steps;
						correct = true;
						tempA0 = lidia_rint(my[k][j]);
						if (alg.A_convert_value_E_bound(da, tempE0, tempA0, bound))
							bigmy = true;
						//
						// compute and store g[k][k]
						//
						LiDIA::multiply(tempE1, tempE0, da.s.value[k][j]);
						tempE1.multiply_by_2();
						LiDIA::subtract(tempE1, da.s.value[k][k], tempE1);
						LiDIA::square(tempE2, tempE0);
						LiDIA::multiply(tempE2, tempE2, da.s.value[j][j]);
						LiDIA::add(tempE1, tempE1, tempE2);

						vector.exact.scalsub(da.s.value[k], da.s.value[k],
								     tempE0, da.s.value[j]);
						//
						// set the right value of g[k][k]
						//
						LiDIA::swap(tempE1, da.s.value[k][k]);

						//
						// set also the column
						//
						for (lidia_size_t i = 0; i < da.b.rows; i++)
							da.s.value[i][k].assign(da.s.value[k][i]);

						for (lidia_size_t i = 0; i < j; ++i)
							my[k][i] -= tempA0*my[j][i];
						my[k][j] -= tempA0;
					}
				vector.exact.vectsize = da.b.columns;
			}

			//
			// Correction Step II
			// needed after Reduction, correct only reduced vector
			// valueA[k] = (A) value[k];
			//
			if (correct) {
				++li.lll.correction_steps;
				if (alg.E_convert_vector_A(da, valueA, k)) {
					//
					// Found 0 - Vector ->repeat with k = 1;
					//
					bigmy = true;
					repeat_loop = false;
					for (lidia_size_t i = k; i < da.b.rank; i++)
						for (lidia_size_t j = 0; j < da.b.rows; j++)
							LiDIA::swap(da.s.value[j][i], da.s.value[j][i+1]);
					alg.E_convert_lattice_A(da, valueA);
					k = 1;
				}
				else {
					//
					// Also correct the corresponding column
					//
					for (lidia_size_t i = 0; i < da.b.rows; i++)
						valueA[i][k] = valueA[k][i];
				}
			}

			//
			// end of reduction
			//
			if (bigmy) {
				//
				// Go back if my to big
				//
				if (k > 1)
					--k;
			}
			else {
				//
				// repeat only once
				//
				bigmy = (!repeat_loop) && correct;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("lll_kernel_op", "lll_gram(da, li) [4]");
		//
		// Check lll - conditions for parameter y_par
		//
		if (da.b.y*vect[k-1] > my[k][k-1]*my[k][k-1]*vect[k-1]+vect[k]) {
			//
			// swap vectors at position k and k-1
			//
			++li.lll.swaps;
			vector.exact.swap(da.s.value[k], da.s.value[k-1]);
			vector.approx.swap(valueA[k], valueA[k-1]);
			//
			// swap also the column
			//
			for (lidia_size_t i = 0; i < da.b.rows; i++) {
				LiDIA::swap(da.s.value[i][k], da.s.value[i][k-1]);
				tempA0 = valueA[i][k];
				valueA[i][k] = valueA[i][k-1];
				valueA[i][k-1] = tempA0;
			}

			if (k > 1)
				--k;

			if (k == 1) {
				vect[0] = valueA[0][0];
				vect_sq[0] = sqrt(vect[0]);
			}
		}
		else
			++k;
	}

	li.lll.rank = da.b.rank;
	alg.prec_restore(da);
	//
	// Free allocated storage
	//
	delete[] mydel;
	delete[] my;
}



//
// Kernel implemented for Functions with enhanced functionality.
// These algorithms do more correction steps, increase precision
// if needed.
//
//
// Schnorr-Euchner-lll
// reducing lattice of type E, approximating with type A
//
template< class E, class A, class VectorOp, class AlgParts >
void
lll_kernel_fu< E, A, VectorOp, AlgParts >::lll(dense_alg< E > & da,
					       lattice_info& li)
{
	debug_handler("lll_kernel_fu", "lll(da, li) [1]");
	bool bigmy;
	bool correct;
	bool repeat_loop;
	lidia_size_t k;
	sdigit old_bit_prec;
	E tempE0;
	E tempE1;
	A bound;
	A invbound;
	A tempA0;
	A tempA1;
	A *mydel;
	A **my;
	A **valueA;
	A *vect;
	A *vect_sq;
	A halb(0.5);

	//
	// Allocating Memory for matrix
	//
	my = new A*[da.b.rows*2];
	memory_handler(my, "lll_kernel_fu", "lll(da, li) :: "
		       "not enough memory !");
	mydel = new A[da.b.rows*((da.b.rows+da.b.columns)+2)];
	memory_handler(mydel, "lll_kernel_fu", "lll(da, li) :: "
                       "not enough memory !");
	valueA = &my[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		my[i] = &mydel[i*da.b.rows]; // rowsxrows
		valueA[i] = &mydel[(da.b.rows*da.b.rows)+(i*da.b.columns)]; // rowsxcolumns
	}

	//
	// Setting pointers
	//
	vect = &mydel[(da.b.rows+da.b.columns)*da.b.rows];
	vect_sq = &vect[da.b.rows];

	//
	// Start
	//
	// Compute digit precision and bound for lll - reduction
	// also set precision
	//
	alg.prec_startup(da);
	alg.prec_update(da);
	bound.assign(std::exp(std::log(2.0)*da.d.bit_prec/2));
	LiDIA::divide(invbound, 1, bound);

	//
	// Copy into small A
	//
	alg.E_convert_lattice_A(da, valueA);
	alg.prec_approx(da);

	li.lll.y = da.b.y;
	li.lll.reduction_steps = 0;
	li.lll.correction_steps = 0;
	li.lll.swaps = 0;
	old_bit_prec = da.d.bit_prec;
	vector.approx.vectsize = 2*da.b.rows;
	vector.approx.assign_zero(vect);
	vector.approx.vectsize = da.b.columns;
	vector.exact.vectsize = da.b.columns;
	vector.approx.scalprod(vect[0], valueA[0], valueA[0]);

	k = 1;
	while (k < da.b.rank) {
		bigmy = true;
		repeat_loop = false;
		while (bigmy) {
			//
			// begin of orthogonalization
			//
			debug_handler("lll_kernel_fu", "lll(da, li) [2]");
			vector.approx.scalprod(vect[k], valueA[k], valueA[k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				vector.approx.scalprod(tempA0, valueA[k], valueA[j]);
				LiDIA::multiply(tempA1, invbound, vect_sq[k]);
				LiDIA::multiply(tempA1, tempA1, vect_sq[j]);
				//
				// Correction Step I, only  if nessesary |valueA[k]*valueA[j]| to big
				//
				if (abs(tempA0).compare(tempA1) < 0) {
					//
					// Compute correct scalar product
					//
					alg.prec_correct(da);
					vector.exact.scalprod(tempE0, da.s.value[k],
							      da.s.value[j]);
					da.d.bit_factor <<= 1;
					alg.E_convert_value_A(da, my[k][j], tempE0);
					alg.prec_approx(da);
					da.d.bit_factor >>= 1;
				}
				else
					my[k][j].assign(tempA0);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempA1, my[j][i], my[k][i]);
					LiDIA::multiply(tempA1, tempA1, vect[i]);
					LiDIA::subtract(my[k][j], my[k][j], tempA1);
				}

				tempA0.assign(my[k][j]);
				LiDIA::divide(my[k][j], my[k][j], vect[j]);
				LiDIA::multiply(tempA1, my[k][j], tempA0);
				LiDIA::subtract(vect[k], vect[k], tempA1);
			}

			//
			// end of orthogonalization
			//
			bigmy = false;
			correct = false;
			//
			// begin of reduction
			//
			debug_handler("lll_kernel_fu", "lll(da, li) [3]");
			if (!repeat_loop) {
				//
				// trick to allow generating a transformation matrix
				// with the same algorithm
				//
				vector.exact.vectsize = da.b.real_columns;
				for (lidia_size_t j = k-1; j >= 0; j--) {
					LiDIA::subtract(tempA0, halb, my[k][j]);
					if ((abs(my[k][j]).compare(halb) > 0) &&
					    (!(tempA0.is_approx_zero()))) {
						++li.lll.reduction_steps;
						correct = true;
						tempA0.assign(round(A(my[k][j])));
						//
						// Bug, I don't know why
						// ::round(tempA0, (A)my[k][j]);
						//
						if (alg.A_convert_value_E_bound(da, tempE0,
										tempA0, bound))
							bigmy = true;
						alg.prec_exact(da);
						vector.exact.scalsub(da.s.value[k], da.s.value[k],
								     tempE0, da.s.value[j]);
						alg.prec_approx(da);
						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempA1, tempA0, my[j][i]);
							LiDIA::subtract(my[k][i], my[k][i], tempA1);
						}
						LiDIA::subtract(my[k][j], my[k][j], tempA0);
					}
				}
				vector.exact.vectsize = da.b.columns;
			}

			//
			// Correction Step II
			// needed after Reduction, correct only reduced vector
			// valueA[k] = (A) value[k];
			//
			if (correct) {
				++li.lll.correction_steps;
				if (alg.E_convert_vector_A(da, valueA, k)) {
					bigmy = true;
					repeat_loop = false;
					k = 1;
					vector.approx.scalprod(vect[0], valueA[0], valueA[0]);
					sqrt(vect_sq[0], vect[0]);
				}
			}

			//
			// end of reduction
			//
			if (bigmy) {
				if (k > 1)
					--k;
			}
			else {
				//
				// repeat only once
				//
				bigmy = (!repeat_loop) && correct;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("lll_kernel_fu", "lll(da, li) [4]");
		tempA0.assign(da.b.y);
		LiDIA::multiply(tempA0, tempA0, vect[k-1]);
		LiDIA::square(tempA1, my[k][k-1]);
		LiDIA::multiply(tempA1, tempA1, vect[k-1]);
		LiDIA::add(tempA1, tempA1, vect[k]);

		//
		// Check lll - conditions for parameter y_par
		//
		if (tempA0.compare(tempA1) > 0) {
			//
			// swap vectors at position k and k-1
			//
			++li.lll.swaps;
			vector.exact.swap(da.s.value[k], da.s.value[k-1]);
			vector.approx.swap(valueA[k], valueA[k-1]);

			if (k > 1)
				--k;
			if (k == 1) {
				vector.approx.scalprod(vect[0], valueA[0], valueA[0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}

	li.lll.rank = da.b.rank;
	alg.prec_restore(da);
	//
	// Free allocated storage
	//
	delete[] mydel;
	delete[] my;
}



//
// Schnorr-Euchner-lll
// modified version working with gram matrices
// reducing lattice of type E, approximating with type A
//
template< class E, class A, class VectorOp, class AlgParts >
void
lll_kernel_fu< E, A, VectorOp, AlgParts >::lll_gram(dense_alg< E > & da,
						    lattice_info& li)
{
	debug_handler("lll_kernel_fu", "lll_gram(da, li) [1]");
	bool bigmy;
	bool correct;
	bool repeat_loop;
	lidia_size_t k;
	sdigit old_bit_prec;
	E tempE0;
	E tempE1;
	E tempE2;
	A bound;
	A tempA0;
	A tempA1;
	A *mydel;
	A **my;
	A **valueA;
	A *vect;
	A *vect_sq;
	A halb(0.5);

	//
	// Allocating Memory for matrix
	//
	my = new A*[da.b.rows*2];
	memory_handler(my, "lll_kernel_fu", "lll_gram(da, li) :: "
		       "not enough memory !");
	mydel = new A[da.b.rows*((da.b.rows+da.b.columns)+2)];
	memory_handler(mydel, "lll_kernel_fu", "lll_gram(da, li) :: "
                       "not enough memory !");
	valueA = &my[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		my[i] = &mydel[i*da.b.rows]; // rows x rows
		valueA[i] = &mydel[(da.b.rows*da.b.rows)+(i*da.b.columns)]; // rows x columns
	}

	//
	// Setting pointers
	//
	vect = &mydel[(da.b.rows+da.b.columns)*da.b.rows];
	vect_sq = &vect[da.b.rows];

	//
	// Start
	//
	// Compute digit precision and bound for lll - reduction
	// also set precision
	//
	alg.prec_startup(da);
	alg.prec_update(da);
	bound.assign(std::exp(std::log(2.0)*da.d.bit_prec/2));

	//
	// Copy into small A
	//
	alg.E_convert_lattice_A(da, valueA);
	if (da.b.rank < da.b.rows) {
		for (lidia_size_t i = 0; i < da.b.rows; i++)
			for (lidia_size_t j = i; j < da.b.columns; j++)
				LiDIA::swap(da.s.value[i][j], da.s.value[j][i]);
		da.b.rank = da.b.columns;
		alg.E_convert_lattice_A(da, valueA);
	}
	alg.prec_approx(da);

	li.lll.y = da.b.y;
	li.lll.reduction_steps = 0;
	li.lll.correction_steps = 0;
	li.lll.swaps = 0;
	old_bit_prec = da.d.bit_prec;
	vector.approx.vectsize = 2*da.b.rows;
	vector.approx.assign_zero(vect);
	vect[0].assign(valueA[0][0]);
	vector.approx.vectsize = da.b.columns;
	vector.exact.vectsize = da.b.columns;

	k = 1;
	while (k < da.b.rank) {
		bigmy = true;
		repeat_loop = false;
		while (bigmy) {
			//
			// begin of orthogonalization
			//
			debug_handler("lll_kernel_fu", "lll_gram(da, li) [2]");
			vect[k].assign(valueA[k][k]);
			sqrt(vect_sq[k], vect[k]);
			for (lidia_size_t j = 0; j < k; ++j) {
				my[k][j].assign(valueA[k][j]);

				for (lidia_size_t i = 0; i < j; ++i) {
					LiDIA::multiply(tempA1, my[j][i], my[k][i]);
					LiDIA::multiply(tempA1, tempA1, vect[i]);
					LiDIA::subtract(my[k][j], my[k][j], tempA1);
				}

				tempA0.assign(my[k][j]);
				LiDIA::divide(my[k][j], my[k][j], vect[j]);
				LiDIA::multiply(tempA1, my[k][j], tempA0);
				LiDIA::subtract(vect[k], vect[k], tempA1);
			}

			//
			// end of orthogonalization
			//
			bigmy = false;
			correct = false;
			//
			// begin of reduction
			//
			debug_handler("lll_kernel_fu", "lll_gram(da, li) [3]");
			if (!repeat_loop) {
				//
				// trick to allow generating a transformation matrix
				// with the same algorithm
				//
				vector.exact.vectsize = da.b.real_columns;
				for (lidia_size_t j = k-1; j >= 0; j--) {
					LiDIA::subtract(tempA0, halb, my[k][j]);
					if ((abs(my[k][j]).compare(halb) > 0) &&
					    (!(tempA0.is_approx_zero()))) {
						++li.lll.reduction_steps;
						correct = true;
						tempA0.assign(round(A(my[k][j])));
						//
						// Bug, I don't know why
						// ::round(tempA0, (A)my[k][j]);
						//
						if (alg.A_convert_value_E_bound(da, tempE0,
										tempA0, bound))
							bigmy = true;
						//
						// compute and store g[k][k]
						//
						alg.prec_exact(da);
						LiDIA::multiply(tempE1, tempE0, da.s.value[k][j]);
						tempE1.multiply_by_2();
						LiDIA::subtract(tempE1, da.s.value[k][k], tempE1);
						LiDIA::square(tempE2, tempE0);
						LiDIA::multiply(tempE2, tempE2, da.s.value[j][j]);
						LiDIA::add(tempE1, tempE1, tempE2);

						vector.exact.scalsub(da.s.value[k], da.s.value[k],
								     tempE0, da.s.value[j]);
						//
						// set the right value of g[k][k]
						//
						LiDIA::swap(tempE1, da.s.value[k][k]);

						//
						// set also the column
						//
						for (lidia_size_t i = 0; i < da.b.rows; i++)
							da.s.value[i][k].assign(da.s.value[k][i]);

						alg.prec_approx(da);
						for (lidia_size_t i = 0; i < j; ++i) {
							LiDIA::multiply(tempA1, tempA0, my[j][i]);
							LiDIA::subtract(my[k][i], my[k][i], tempA1);
						}
						LiDIA::subtract(my[k][j], my[k][j], tempA0);
					}
				}
				vector.exact.vectsize = da.b.columns;
			}

			//
			// Correction Step II
			// needed after Reduction, correct only reduced vector
			// valueA[k] = (A) value[k];
			//
			if (correct) {
				++li.lll.correction_steps;
				if (alg.E_convert_vector_A(da, valueA, k)) {
					bigmy = true;
					repeat_loop = false;
					for (lidia_size_t i = k; i < da.b.rank; i++)
						for (lidia_size_t j = 0; j < da.b.rows; j++)
							LiDIA::swap(da.s.value[j][i], da.s.value[j][i+1]);
					alg.E_convert_lattice_A(da, valueA);
					k = 1;
				}
				else {
					//
					// Also correct the corresponding column
					//
					//
					for (lidia_size_t i = 0; i < da.b.rows; i++)
						valueA[i][k].assign(valueA[k][i]);
				}
			}

			//
			// end of reduction
			//
			if (bigmy) {
				if (k > 1)
					--k;
			}
			else {
				//
				// repeat only once
				//
				bigmy = (!repeat_loop) && correct;
				repeat_loop = !repeat_loop;
			}
		}

		debug_handler("lll_kernel_fu", "lll_gram(da, li) [4]");
		tempA0.assign(da.b.y);
		LiDIA::multiply(tempA0, tempA0, vect[k-1]);
		LiDIA::square(tempA1, my[k][k-1]);
		LiDIA::multiply(tempA1, tempA1, vect[k-1]);
		LiDIA::add(tempA1, tempA1, vect[k]);

		//
		// Check lll - conditions for parameter y_par
		//
		if (tempA0.compare(tempA1) > 0) {
			//
			// swap vectors at position k and k-1
			//
			++li.lll.swaps;
			vector.exact.swap(da.s.value[k], da.s.value[k-1]);
			vector.approx.swap(valueA[k], valueA[k-1]);

			//
			// swap also the column
			//
			for (lidia_size_t i = 0; i < da.b.rows; i++) {
				LiDIA::swap(da.s.value[i][k], da.s.value[i][k-1]);
				LiDIA::swap(valueA[i][k], valueA[i][k-1]);
			}

			if (k > 1)
				--k;
			if (k == 1) {
				vect[0].assign(valueA[0][0]);
				sqrt(vect_sq[0], vect[0]);
			}
		}
		else
			++k;
	}

	li.lll.rank = da.b.rank;
	alg.prec_restore(da);
	//
	// Free allocated storage
	//
	delete[] mydel;
	delete[] my;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_LATTICE_KERNEL_CC_GUARD_
