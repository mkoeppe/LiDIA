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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


// functions realize search for eigenvalue with division polynomials
// ev is used to return eigenvalue, fC is Elkies polynomial
// klist[1, ..., klist[0]] is list of possible candidates for ev
// test_all : controls the search for alpha in
//            Phi (X, Y) == alpha * (X, Y)
//
//             0 --> stop searching after first match
//             1 --> test all possible values for alpha
//                   or stop after second match
//



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/weco2_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool eco_gf2n::schoofpart (lidia_size_t & ev, const ff_pol & fC,
			   lidia_size_t* klist, bool test_all)
{
	lidia_size_t   *alpha = NULL;
	lidia_size_t   alpha_num = 1;
	lidia_size_t   k, d;
	lidia_size_t   top;
	lidia_size_t   k_ix;
	ff_element   tmp;
	PsiPowers  *psi_pow; // for storing division polynomials
	ff_polmod ffC; // build poly_modulus for Elkies polynomial
	ff_pol  Q1, Q2; // the right side of the equation phi (P) = k*P

	ffC.build (fC);

#ifdef TIMING
	timer t;
	t.set_print_mode();
#endif


	if (test_all) {
		// allocating memory for possible solutions
		alpha = new lidia_size_t[l+1];
		memory_handler(alpha, "eco_gf2n::schoofpart()", "Allocating alpha");
		alpha[0]  = l+1;
		alpha_num = 1;
	}

	d = static_cast<lidia_size_t>(klist[klist[0]]); // last possibility for ev

	k_ix = 1;

	//****************************** Test Y coordinate ****************

	if (test_y_coordinate()) {
		ff_pol  yq0, yq1; // Y^q = yq0 + Y * yq1 mod fC
		ff_pol Q3;

#ifdef TIMING
		t.start_timer();
#endif

		Ytop_f (yq0, yq1, ffC);

#ifdef TIMING
		t.stop_timer();
		if (info)
			std::cout << "\nComputation of Y^q mod f_C(X) needs time : " << t << "\n" << std::flush;
#endif

		if (info) {
			std::cout << "\nTesting Y-coordinate of " << klist[0];
			if (klist[0] == 1)
				std::cout << " possibility using division polynomials  :\n\n";
			else
				std::cout << " possibilities using division polynomials  :\n\n";
		}

		// Test with k == +-1

		if (klist[1] == 1) {
			if (info)
				std::cout << "1 " << std::flush;

			if (yq1.is_one() && yq0.is_zero()) {
				// match for k = +1.
				if (!test_all) {
					ev = 1;
					if (info)
						std::cout << "\n\nEigenvalue of Frobenius is 1" << std::flush;
					return true;
				}
				alpha[alpha_num++] = 1;
			}
			else
				if (yq1.is_one() && yq0.is_x())   // k == -1 ??
				{
					if (!test_all) {
						ev = -1;
						if (info)
							std::cout << "\n\nEigenvalue of Frobenius is -1" << std::flush;
						return true;
					}
					alpha[alpha_num++] = -1;
				}
			k_ix ++;
		}

		//*********** Tests for abs(k) > 1 ******************

		if (d > 1) {
			int kk, deleted_up_to = -1;

			// initialize division polynomials

			init_psi_comp (psi_pow, top, d+2, ffC);

			// now test successively

			for (k = 2; k <= klist[klist[0]]; k++) {
				if (k+1 > top) {
					next_psi (psi_pow, top, ffC);

					if (k & 1)
						kk = ((k+2) >> 1)-2;
					else
						kk = ((k+2) >> 1)-3;

					if (kk > deleted_up_to) {
						free_psi(psi_pow, deleted_up_to+1, kk);
						deleted_up_to = kk;
					}
				}

				if (k == klist[k_ix]) {
					if (info)
						std::cout << k << " " << std::flush;

					k_ix ++;

					// first we test the y1-part, since it is cheaper

					shift_left(Q1, *(psi_pow[k]).pow2, ffC, 1);
					multiply(Q2, *(psi_pow[k-1]).pow1, *(psi_pow[k+1]).pow1,
						 ffC);
					add (Q2, Q1, Q2);
					multiply(Q1, yq1, Q1, ffC);


					if (Q1 == Q2) {
						// match on y1, test y0 now
						multiply(Q1, *(psi_pow[k-2]).pow1, *(psi_pow[k+1]).pow2,
							 ffC);

						multiply(Q2, *(psi_pow[k]).pow1, *(psi_pow[k-1]).pow1,
							 ffC);
						multiply(Q3, Q2, *(psi_pow[k+1]).pow1, ffC);
						shift_left(Q2, Q3, ffC, 1);
						add(Q1, Q1, Q2);

						shift_left(Q2, Q2, ffC, 1);
						add(Q1, Q1, Q2);

						shift_left(Q2, *(psi_pow[k]).pow3, ffC, 2);
						add(Q1, Q1, Q2);

						shift_left(Q2, *(psi_pow[k]).pow3, ffC, 1);
						multiply(Q2, Q2, yq0, ffC);

						if (Q1 == Q2) {
							if (!test_all) {
								ev = k;

								free_psi (psi_pow, deleted_up_to+1, d+5);
								delete[] psi_pow;
								if (info) {
									std::cout << "\n\nEigenvalue of Frobenius is ";
									std::cout << k << std::flush;
								}
								return true;
							}
							else
								alpha[alpha_num++] = k;
						}
						else {
							shift_left(Q3, Q3, ffC, 1);
							add(Q1, Q1, Q3);

							shift_left(Q3, *(psi_pow[k]).pow3, ffC, 2);
							add(Q1, Q1, Q3);

							if (Q1 == Q2) {
								if (!test_all) {
									ev = -k;

									free_psi (psi_pow, deleted_up_to+1, d+5);
									delete[] psi_pow;
									if (info) {
										std::cout << "\n\nEigenvalue of Frobenius is ";
										std::cout << -k << std::flush;
									}
									return true;
								}
								else
									alpha[alpha_num++] = -k;
							}
						}
					}
				} // end for
			}

			free_psi (psi_pow, deleted_up_to+1, d+5);
			delete[] psi_pow;
		} // end if d > 1

		if (alpha_num == 2)
			ev = alpha[1];
		delete[] alpha;

		if (alpha_num == 1) {
			lidia_error_handler("eco_gf2n", "schoof_part::no eigenvalue found");
			return false;
		}

		if (alpha_num > 2) {
			lidia_error_handler("eco_gf2n", "schoof_part::no unique "
					    "eigenvalue found");
			return false;
		}
		return true;
	} // end of "if (test_y_coordinate())"


	//**************************** Test X coordinate *********************

	if (test_x_coordinate()) {
		ff_pol  xtop; // x^q
		ff_pol tmp_pol;

#ifdef TIMING
		t.start_timer();
#endif

		Xq (xtop, ffC, A6.relative_degree()); // fast exponentiation for x^q

#ifdef TIMING
		t.stop_timer();
		if (info)
			std::cout << "\nComputation of X^q mod f_C(X) needs time : " << t << "\n" << std::flush;
#endif


		if (info) {
			std::cout << "\nTesting X-coordinate of " << klist[0];
			if (klist[0] == 1)
				std::cout << " possibility using division polynomials  :\n\n";
			else
				std::cout << " possibilities using division polynomials  :\n\n";
		}

		// Test with k == +-1

		if (klist[1] == 1) {
			if (info)
				std::cout << "1 " << std::flush;

			if (xtop.is_x()) {
				// match for k = +/- 1.
				if (!test_all) {
					ev = 1;
					if (info)
						std::cout << "\n\nEigenvalue of Frobenius is +/- 1" << std::flush;
					return true;
				}
				alpha[alpha_num++] = 1;
			}
			k_ix ++;
		}

		//*********** Tests for abs(k) > 1 ******************

		if (d > 1) {
			// initialize division polynomials

			init_psi_comp (psi_pow, top, d+1, ffC);

			// Test for k == +-2

			if (klist[k_ix] == 2) {
				k_ix ++;
				shift_left(Q2, *(psi_pow[2]).pow2, ffC);
				add(Q2, Q2, *(psi_pow[3]).pow1);
				multiply(Q1, xtop, *(psi_pow[2]).pow2, ffC);

				if (info)
					std::cout << "2 " << std::flush;

				if (Q1 == Q2) {
					if (!test_all) {
						ev = 2;
						if (info)
							std::cout << "\n\nEigenvalue of Frobenius is +/- 2" << std::flush;
						free_psi(psi_pow, 0, d+4);
						delete[] psi_pow;
						return true;
					}
					else
						alpha[alpha_num++] = 2;
				}
			}

			// now all possibilities 3 <= k <= d

			int kk, deleted_up_to = -1;

			for (k = 3; k <= klist[klist[0]] && (alpha_num == 1 || test_all);
			     k++) {
				if (k+1 > top) {
					next_psi (psi_pow, top, ffC);

					if (k & 1)
						kk = ((k+2) >> 1)-2;
					else
						kk = ((k+2) >> 1)-3;

					if (kk > deleted_up_to) {
						free_psi(psi_pow, deleted_up_to+1, kk);
						deleted_up_to = kk;
					}
				}

				if (k == klist[k_ix]) {
					if (info)
						std::cout << k << " " << std::flush;

					k_ix ++;

					shift_left (Q2, *(psi_pow[k].pow2), ffC);
					multiply (Q1, *(psi_pow[k+1].pow1), *(psi_pow[k-1].pow1),
						  ffC);
					add (Q2, Q2, Q1);
					multiply(Q1, xtop, *(psi_pow[k].pow2), ffC);

					if (Q1 == Q2) {
						if (!test_all) {
							ev = k;
							free_psi(psi_pow, deleted_up_to+1, d+4);
							delete[] psi_pow;
							if (info)
								std::cout << "\n\nEigenvalue of Frobenius is +/- " << k << std::flush;
							return true;
						}
						else
							alpha[alpha_num++] = k;
					}
				} // endif k == klist[]
			} // end for

			// free_psi (psi_pow, deleted_up_to+1, d+3);
			free_psi(psi_pow, deleted_up_to+1, d+4);
			delete[] psi_pow;
		} // end if d > 1

		if (alpha_num == 2)
			ev = alpha[1];
		delete[] alpha;

		if (alpha_num == 1) {
			lidia_error_handler("eco_gf2n", "schoof_part::no eigenvalue found");
			return false;
		}

		if (alpha_num > 2) {
			lidia_error_handler("eco_gf2n", "schoof_part::no unique eigenvalue found");
			return false;
		}
		return true;
	} // end of "if (test_x_coordinate())"

	lidia_error_handler("eco_gf2n", "invalid strategy");
	return false;
}



//------------------------------------------------------------------------
// now we implement the same functionality, but using rational functions
// instead of division polynomials
// Here we have 2 different situations: test x or y coordinate
//                                      with or without table
//
// For table usage: we use three tables:
//
//   i*P is still used number_of_table_usage[i] times, if zero, deallocate
//       memory
//   k*P is stored in table[j], if index_table[j] = k
//   table[...] stores points

bool eco_gf2n::schoofpart_rat_function (lidia_size_t & ev,
					const ff_pol & fC,
					lidia_size_t* klist, bool test_all)
{
	lidia_size_t   *alpha = NULL;
	lidia_size_t   alpha_num = 1;
	lidia_size_t   k, d, i = 0;
	lidia_size_t   k_ix;
	ff_pol  y0, y1; // Y^q = y0 + y1 *Y mod fC for Y-test
	//  y0 = X^q mod fC for X-test
	ff_polmod ffC;
	int last_m = 0;
	bool table_used = false;

	const int size_table = 100;
	long *number_of_table_usage = NULL;
	weco2_rat_function *table = NULL;
	long *index_table = NULL;

#ifdef TIMING
	timer t;
	t.set_print_mode();
#endif

	ffC.build (fC); // first build ff_polmod

#ifdef TIMING
	t.start_timer();
#endif


	if (test_x_coordinate())   // fast exponentiation
		Xq(y0, ffC, A6.relative_degree());
	else
		Ytop_f (y0, y1, ffC);

#ifdef TIMING
	t.stop_timer();
	if (info)
		if (test_x_coordinate())
			std::cout << "\nComputation of X^q mod f_C(X) needs time : " << t << "\n" << std::flush;
		else  std::cout << "\nComputation of Y^q mod f_C(X) needs time : " << t << "\n" << std::flush;
#endif



	if (test_all) {
		alpha = new lidia_size_t[l+1];
		memory_handler(alpha, "eco_gf2n::schoofpart_rat_function",
			       "Allocating alpha");
		alpha[0]  = l+1;
		alpha_num = 1;
	}

	d = static_cast<lidia_size_t>(klist[klist[0]]);
	k_ix = 1;

	if (info) {
		std::cout << "\nTesting " << klist[0];
		if (klist[0] == 1)
			std::cout << " possibility using rational functions :\n\n";
		else
			std::cout << " possibilities using rational functions :\n\n";
	}


	if (klist[1] == 1) {
		// Test with k == +- 1
		if (info)
			std::cout << "1 " << std::flush;

		if ((y1.is_one() && y0.is_zero() && test_y_coordinate()) ||
		    (y0.is_x() && test_x_coordinate())) {
			if (!test_all) {
				ev = 1;
				if (info)
					if (test_y_coordinate())
						std::cout << "\n\nEigenvalue of Frobenius is 1" << std::flush;
					else
						std::cout << "\n\nEigenvalue of Frobenius is +/- 1" << std::flush;
				return true;
			}
			alpha[alpha_num++] = 1;
		}

		if (test_y_coordinate() && y1.is_one() && y0.is_x()) {
			if (!test_all) {
				ev = -1;
				if (info)
					std::cout << "\n\nEigenvalue of Frobenius is -1" << std::flush;
				return true;
			}
			alpha[alpha_num++] = -1;
		}
		k_ix ++;
	}

	//** now rational functions become necessary ***************

	if (klist[0] > 1 || (klist[1] != 1 && klist[0] == 1)) {
		weco2_rat_function::initialize(eco_gf2n::A6, ffC);
		weco2_rat_function P, mP, dP; // P = (X, Y)
		weco2_rat_function::ff_rat help;

		P.assign_xy();
		mP.assign_zero();

		if ((ev_strategy & 0x0f) == EV_RATIONAL_FUNCTION_TABLE) {
			table_used = true;
			number_of_table_usage = new long[d+1];

			for (i = 0; i <= d; i++)
				// negative value == table not
				number_of_table_usage[i] = 0; // yet initialized

			for (i = 2; i <= klist[0]; i++)
				number_of_table_usage[klist[i] - klist[i-1]] --;

			table = new weco2_rat_function [size_table];
			index_table = new long[size_table];

			for (i = 0; i < size_table; i++)
				index_table[i] = 0;
		}

		//*********** Tests for abs(k) > 1 ******************

		for (k = 2; k <= d; k++) {
			if (k == klist[k_ix]) {
				if (info)
					std::cout << k << " " << std::flush;

				k_ix ++;

				if (table_used) {
					// find point in table first
					// determine also the maximal element smaller than k - last_m

					int max_smaller = 0, i_max = -1;

					i = 0;
					while(index_table[i] != k-last_m && index_table[i] != 0
					      && i < size_table) {
						if (index_table[i] < k - last_m
						    && index_table[i] > max_smaller) {
							i_max = i; max_smaller = index_table[i];
						}
						i++;
					}

					if (index_table[i] == 0 || i >= size_table) {
						if (i_max > -1) {
							weco2_rat_function HP;
							multiply(HP, (k - last_m) / max_smaller, table[i_max]);
							multiply(dP, (k - last_m) % max_smaller, P);
							add(dP, HP, dP);
						}
						else
							multiply(dP, (k - last_m), P);
					}
					else
						dP.assign(table[i]);
				}
				else
					multiply(dP, k-last_m, P);

				if (table_used && index_table[i] == 0 && i < size_table) {
					i = 0;
					while(index_table[i] != 0 && i < size_table)
						i++;

					if (i < size_table) {
						table[i].assign(dP);
						index_table[i] = k-last_m;
					}
				}

				last_m = k;
				add(mP, mP, dP);

				if ((equal(y0, mP.get_y0(), ffC) &&
				     equal(y1, mP.get_y1(), ffC) && test_y_coordinate()) ||
				    (equal(y0, mP.get_x(), ffC) && test_x_coordinate())) {
					if (!test_all) {
						ev = k;
						if (info)
							if (test_y_coordinate())
								std::cout << "\n\nEigenvalue of Frobenius is " << k << std::flush;
							else
								std::cout << "\n\nEigenvalue of Frobenius is +/- " << k << std::flush;
						if (table_used) {
							delete[] index_table;
							delete[] number_of_table_usage;
							delete[] table;
						}
						return true;
					}
					else
						alpha[alpha_num++] = k;
				}

				if (test_y_coordinate()) {
					add(help, mP.get_y0(), mP.get_x(), ffC);

					if (equal(y1, mP.get_y1(), ffC)
					    && equal(y0, help, ffC)) {
						if (!test_all) {
							ev = -k;
							if (info)
								std::cout << "\n\nEigenvalue of Frobenius is " << -k << std::flush;
							if (table_used) {
								delete[] index_table;
								delete[] number_of_table_usage;
								delete[] table;
							}
							return true;
						}
						else
							alpha[alpha_num++] = -k;
					}
				}
			} // endif k == klist[]
		} // end for
	}

	if (table_used) {
		delete[] index_table;
		delete[] number_of_table_usage;
		delete[] table;
	}

	if (alpha_num == 2)
		ev = alpha[1];
	delete[] alpha;

	if (alpha_num == 1) {
		lidia_error_handler("eco_gf2n", "schoof_part::no eigenvalue found");
		return false;
	}

	if (alpha_num > 2) {
		lidia_error_handler("eco_gf2n", "schoof_part::no unique eigenvalue found");
		return false;
	}

	return true;
}



void eco_gf2n::compute_sign (lidia_size_t & alpha,
			     const ff_pol & fC)
{
	lidia_size_t klist[2];
	char ev = ev_strategy;

	klist[0] = 1;
	if (alpha < 0)
		klist[1] = l - alpha;
	else
		klist[1] = alpha;

	ev_strategy = eco_gf2n::EV_RATIONAL_FUNCTION;

	schoofpart_rat_function(alpha, fC, klist, false);

	if (alpha < 0)
		alpha += l;

	ev_strategy = ev;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
