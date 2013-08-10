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
//	Author	: Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================

//
// functions realize search for eigenvalue with division polynomials
// ev is used to return eigenvalue, fC is Elkies polynomial
// klist[1, ..., klist[0]] is list of possible candidates for ev, in the
// special format described in compute_lists.c
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
#include	"LiDIA/eco_prime.h"
#include	"LiDIA/wep_rat_function.h"

#include <iostream>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool eco_prime::schoofpart (lidia_size_t & ev, const ff_pol & fC,
			    lidia_size_t* klist, bool test_all)
{
  lidia_size_t   *alpha = NULL;
  lidia_size_t   alpha_num = 1;
  lidia_size_t   k, d;
  lidia_size_t   top;
  lidia_size_t   k_ix;
  ff_element   tmp;
  PsiPowers  *psi_pow; // for storing division polynomials
  ff_pol  sqrYY; // y^4 = (x^3+ax+b)^2
  ff_polmod ffC; // build ff_polmod for Elkies polynomial
  ff_pol  Q1, Q2; // the right side of the equation phi (P) = k*P
  ff_pol  Mone; // Mone = -1 * x^0
  ff_element inv2;
  
  ffC.build (fC);
  Q1.set_modulus(p);
  Q2.set_modulus(p);
  Mone.set_modulus(p); Mone.set_coefficient(-1, 0);
  
#ifdef TIMING
  timer t;
  t.set_print_mode();
#endif
  
  if (test_all) {
    alpha = new lidia_size_t[l+1];
    memory_handler(alpha, "eco_prime::schoofpart()", "Allocating alpha");
    alpha[0]  = l+1;
    alpha_num = 1;
  }
  
  d = static_cast<lidia_size_t>(klist[klist[0]]); // last possibility for ev
  
  k_ix = 1;
  invert (inv2, ff_element(2));
  
  // determine sqrYY = Y^4 = (X^3 + A*X + B)^2
  
  sqrYY.set_modulus(p);
  CurveEqn (sqrYY, p, A, B, ffC);
  square   (sqrYY, sqrYY, ffC);
  
  //******************************* Test Y coordinate ****************
  
  if (test_y_coordinate()) {
    ff_pol  ytop; // 4 * y^(q-1) = 4 * (x^3+ax+b)^((q-1)/2)
    ff_pol  ytop1; // 4 * y^(q-1) * y^4 = 4 * (x^3+ax+b)^((q+3)/2)
    
    ytop.set_modulus(p);
    ytop1.set_modulus(p);
    
#ifdef TIMING
    t.start_timer();
#endif

    Ytop_f (ytop, ffC); // fast exponentiation for Y^(q-1)
    
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
      
      if (ytop.is_one()) {
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
	if (ytop == Mone) {
	  // k == -1 ??
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
    
    //************ Tests for abs(k) > 1 ******************
    
    if (d > 1 && (alpha_num == 1 || test_all)) {
      // determine 4 * Y^(q-1) =  ytop
      
      tmp = 4;
      multiply (ytop, ytop, tmp.mantissa());
      
      // determine ytop1 = 4 * Y^4 * Y^(q-1)
      
      multiply (ytop1, sqrYY, ytop, ffC);
      
      // initialize division polynomials
      
      init_psi_comp (psi_pow, top, d+2, ffC);
      
      // Test for k == + 2
      
      if (klist[k_ix] == 2) {
	k_ix ++; k = 2;
	multiply (Q1, *(psi_pow[k-2].pow1), *(psi_pow[k+1].pow2), ffC);
	multiply (Q2, *(psi_pow[k+2].pow1), *(psi_pow[k-1].pow2), ffC);
	subtract (Q2, Q2, Q1);
	
	multiply(Q1, ytop1, *(psi_pow[k].pow3), ffC);
	
	if (info)
	  std::cout << "2 " << std::flush;
	
	if (Q1 == Q2) {
	  if (!test_all) {
	    ev = 2;
	    if (info)
	      std::cout << "\n\nEigenvalue of Frobenius is 2" << std::flush;
	    free_psi(psi_pow, 0, d+4);
	    delete[] psi_pow;
	    return true;
	  }
	  else
	    alpha[alpha_num++] = 2;
	}
	else
	  if (Q1 == -Q2) {
	    if (!test_all) {
	      ev = -2;
	      if (info)
		std::cout << "\n\nEigenvalue of Frobenius is -2" << std::flush;
	      free_psi(psi_pow, 0, d+4);
	      delete[] psi_pow;
	      return true;
	    }
	    else
	      alpha[alpha_num++] = -2;
	  }
      }
      
      // now all possibilities 3 <= k <= d
      
      int kk, deleted_up_to = -1;
      
      for (k = 3; k <= klist[klist[0]] && (alpha_num == 1 || test_all); k++) {
	next_psi (psi_pow, top, ffC, sqrYY, inv2);
	
	if (k & 1)
	  kk = ((k+3) >> 1)-3;
	else
	  kk = ((k+2) >> 1)-2;
	
	if (kk > deleted_up_to) {
	  free_psi(psi_pow, deleted_up_to+1, kk);
	  deleted_up_to = kk;
	}
	
	if (k == klist[k_ix]) {
	  if (info)
	    std::cout << k << " " << std::flush;
	  
	  k_ix ++;
	  
	  multiply (Q1, *(psi_pow[k-2].pow1), *(psi_pow[k+1].pow2),
		    ffC);
	  multiply (Q2, *(psi_pow[k+2].pow1), *(psi_pow[k-1].pow2),
		    ffC);
	  subtract (Q2, Q2, Q1);
	  
	  if (k & 1) {
	    multiply(Q1, ytop, *(psi_pow[k].pow3), ffC);
	    if (Q1 == Q2) {
	      if (!test_all) {
		ev = k;
		free_psi (psi_pow, deleted_up_to+1, d+4);
		delete[] psi_pow;
		if (info)
		  std::cout << "\n\nEigenvalue of Frobenius is " << k << std::flush;
		
		return true;
	      }
	      else
		alpha[alpha_num++] = k;
	    }
	    else
	      if (Q1 == -Q2) {
		if (!test_all) {
		  ev = -k;
		  free_psi (psi_pow, deleted_up_to+1, d+4);
		  delete[] psi_pow;
		  if (info)
		    std::cout << "\n\nEigenvalue of Frobenius is " << -k << std::flush;
		  return true;
		}
		else
		  alpha[alpha_num++] = -k;
	      }
	  }
	  else {
	    multiply(Q1, ytop1, *(psi_pow[k].pow3), ffC);
	    if (Q1 == Q2) {
	      if (!test_all) {
		ev = k;
		free_psi (psi_pow, deleted_up_to+1, d+4);
		delete[] psi_pow;
		if (info)
		  std::cout << "\n\nEigenvalue of Frobenius is " << k << std::flush;
		return true;
	      }
	      else
		alpha[alpha_num++] = k;
	    }
	    else
	      if (Q1 == -Q2) {
		if (!test_all) {
		  ev = -k;
		  free_psi (psi_pow, deleted_up_to+1, d+4);
		  delete[] psi_pow;
		  if (info)
		    std::cout << "\n\nEigenvalue of Frobenius is " << -k << std::flush;
		  return true;
		}
		else
		  alpha[alpha_num++] = -k;
	      }
	  }
	} // endif k == klist[]
      } // end for
      
      free_psi (psi_pow, deleted_up_to+1, d+4);
      delete[] psi_pow;
    } // end if d > 1
    
    if (alpha_num == 2)
      ev = alpha[1];
    delete[] alpha;
    
    if (alpha_num == 1) {
      lidia_error_handler("eco_prime", "schoof_part::no eigenvalue found");
      return false;
    }
    
    if (alpha_num > 2) {
      lidia_error_handler("eco_prime", "schoof_part::no unique eigenvalue found");
      return false;
    }
    return true;
  } // end of "if (test_y_coordinate())"
  
  
  //***************************** Test X coordinate *********************
  
  if (test_x_coordinate()) {
    ff_pol  xtop; // x^q
    ff_pol right_side; // y^2 = x^3 + A*x + B
    ff_pol tmp_pol;
    
    xtop.set_modulus(p);
    tmp_pol.set_modulus(p);
    
#ifdef TIMING
		t.start_timer();
#endif

		power_x (xtop, pn, ffC); // fast exponentiation for x^q

#ifdef TIMING
		t.stop_timer();
		if (info)
			std::cout << "\nComputation of X^q mod f_C(X) needs time : " << t << "\n" << std::flush;
#endif

		right_side.set_modulus(p);
		CurveEqn (right_side, p, A, B, ffC);

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

		//************ Tests for abs(k) > 1 ******************

		if (d > 1 && (alpha_num == 1 || test_all)) {
			// initialize division polynomials

			init_psi_comp (psi_pow, top, d+1, ffC);

			// Test for k == +-2

			if (klist[k_ix] == 2) {
				k_ix ++;
				multiply(Q1, *(psi_pow[2].pow2), right_side, ffC);
				multiply_by_x_mod (Q2, Q1, fC);
				subtract (Q2, Q2, *(psi_pow[3].pow1));
				multiply(Q1, xtop, Q1, ffC);

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

			for (k = 3; k <= klist[klist[0]] && (alpha_num == 1 || test_all); k++) {
				next_psi (psi_pow, top, ffC, sqrYY, inv2);

				if (k & 1)
					kk = ((k+1) >> 1)-2;
				else
					kk = ((k+2) >> 1)-3;

				if (kk > deleted_up_to) {
					free_psi(psi_pow, deleted_up_to+1, kk);
					deleted_up_to = kk;
				}

				if (k == klist[k_ix]) {
					if (info)
						std::cout << k << " " << std::flush;

					k_ix ++;

					if (k & 1) {
						// odd
						multiply_by_x_mod (Q2, *(psi_pow[k].pow2), fC);
						multiply (Q1, *(psi_pow[k+1].pow1), *(psi_pow[k-1].pow1),
							  ffC);
						multiply(Q1, Q1, right_side, ffC);
						subtract (Q2, Q2, Q1);

						multiply(Q1, xtop, *(psi_pow[k].pow2), ffC);
					}
					else {
						multiply(tmp_pol, *(psi_pow[k].pow2), right_side, ffC);
						multiply_by_x_mod (Q2, tmp_pol, fC);
						multiply (Q1, *(psi_pow[k+1].pow1),
							  *(psi_pow[k-1].pow1), ffC);
						subtract (Q2, Q2, Q1);

						multiply(Q1, xtop, tmp_pol, ffC);
					}

					if (Q1 == Q2) {
						if (!test_all) {
							ev = k;
							free_psi (psi_pow, deleted_up_to+1, d+4);
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

			free_psi (psi_pow, deleted_up_to+1, d+4);
			delete[] psi_pow;
		} // end if d > 1

		if (alpha_num == 2)
			ev = alpha[1];
		delete[] alpha;

		if (alpha_num == 1) {
			lidia_error_handler("eco_prime", "schoof_part::no eigenvalue found");
			return false;
		}

		if (alpha_num > 2) {
			lidia_error_handler("eco_prime", "schoof_part::no unique eigenvalue found");
			return false;
		}
		return true;
	} // end of "if (test_x_coordinate())"

	lidia_error_handler("eco_prime", "invalid strategy");
	return false;
}



//------------------------------------------------------------------------
// now we implement the same functionality, but using rational functions
// instead of division polynomials
// Here we have 2 different situations: test x or y coordinate
//                                      with or without table
//


bool eco_prime::schoofpart_rat_function (lidia_size_t & ev,
					 const ff_pol & fC,
					 lidia_size_t* klist, bool test_all)
{
	lidia_size_t   *alpha = NULL;
	lidia_size_t   alpha_num = 1;
	lidia_size_t   k, d;
	lidia_size_t   k_ix;
	ff_element   inv2;
	ff_element   tmp;
	ff_pol  ytop; // 4 * y^(q-1) = 4 * (x^3+ax+b)^((q-1)/2) for y-test
	//  X^q for x-test
	ff_polmod ffC;
	int last_m = 0;
	bool table_used = false;
	long *number_of_table_usage = NULL;
	wep_rat_function *table = NULL;

#ifdef TIMING
	timer t;
	t.set_print_mode();
#endif

	ffC.build (fC); // first build ff_polmod

	ytop.set_modulus(p); // compute yq

#ifdef TIMING
	t.start_timer();
#endif

	if (test_x_coordinate())
		power_x (ytop, pn, ffC);
	else
		Ytop_f (ytop, ffC);

#ifdef TIMING
	t.stop_timer();
	if (info)
		if (test_x_coordinate())
			std::cout << "\nComputation of X^q mod f_C(X) needs time : " << t << "\n" << std::flush;
		else  std::cout << "\nComputation of Y^q mod f_C(X) needs time : " << t << "\n" << std::flush;
#endif


	if (test_all) {
		alpha = new lidia_size_t[l+1];
		memory_handler(alpha, "eco_prime::schoofpart_rat_function", "Allocating alpha");

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

		if ((ytop.is_one() && test_y_coordinate()) ||
		    (ytop.is_x() && test_x_coordinate())) {
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

		if (test_y_coordinate()) {
			negate(ytop, ytop);
			if (ytop.is_one()) {
				// k == -1 ??
				if (!test_all) {
					ev = -1;
					if (info)
						std::cout << "\n\nEigenvalue of Frobenius is -1" << std::flush;
					return true;
				}
				alpha[alpha_num++] = -1;
			}
			negate(ytop, ytop);
		}
		k_ix ++;
	}

	//*** now rational functions become necessary ***************

	if (klist[0] > 1 || (klist[1] != 1 && klist[0] == 1)) {
		wep_rat_function::initialize(eco_prime::A, eco_prime::B, ffC);
		wep_rat_function P, mP, dP; // P = (X, Y)

		P.assign_xy();
		mP.assign_zero();

		if ((ev_strategy & 0x0f) == EV_RATIONAL_FUNCTION_TABLE) {
			number_of_table_usage = new long[d+1];

			int i;
			for (i = 0; i <= d; i++)           // negative value == table not
				number_of_table_usage[i] = 0; // yet initialized

			for (i = 2; i <= klist[0]; i++)
				number_of_table_usage[klist[i] - klist[i-1]] --;

			table = new wep_rat_function[d+1];
			table_used = true;
		}

		//************ Tests for abs(k) > 1 ******************

		if (alpha_num == 1 || test_all) {
			for (k = 2; k <= d && (alpha_num == 1 || test_all); k++) {
				if (k == klist[k_ix]) {
					if (info)
						std::cout << k << " " << std::flush;

					k_ix ++;

					if (table_used && number_of_table_usage[k-last_m] > 0) {
						dP.assign(table[k-last_m]);
						number_of_table_usage[k-last_m] --;
						if (number_of_table_usage[k-last_m] == 0)  // save the space
							table[k-last_m].assign_zero();
					}
					else
						multiply(dP, k-last_m, P);

					if (table_used && (number_of_table_usage[k-last_m] < -1)) {
						number_of_table_usage[k-last_m] = -number_of_table_usage[k-last_m] - 1;
						table[k-last_m].assign(dP);
					}

					last_m = k;
					add(mP, mP, dP);

					if ((equal(ytop, mP.get_y(), ffC) && test_y_coordinate()) ||
					    (equal(ytop, mP.get_x(), ffC) && test_x_coordinate())) {
						if (!test_all) {
							ev = k;
							if (info)
								if  (test_y_coordinate())
									std::cout << "\n\nEigenvalue of Frobenius is " << k << std::flush;
								else
									std::cout << "\n\nEigenvalue of Frobenius is +/- " << k << std::flush;
							if (table_used) {
								delete[] number_of_table_usage;
								delete[] table;
							}
							return true;
						}
						else
							alpha[alpha_num++] = k;
					}

					if (equal(ytop, -mP.get_y(), ffC) && test_y_coordinate()) {
						if (!test_all) {
							ev = -k;
							if (info)
								std::cout << "\n\nEigenvalue of Frobenius is " << -k << std::flush;

							if (table_used) {
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

	if (alpha_num == 2)
		ev = alpha[1];
	delete[] alpha;

	if (alpha_num == 1) {
		lidia_error_handler("eco_prime", "schoof_part::no eigenvalue found");
		return false;
	}

	if (alpha_num > 2) {
		lidia_error_handler("eco_prime", "schoof_part::no unique eigenvalue found");
		return false;
	}

	if (table_used) {
		delete[] number_of_table_usage;
		delete[] table;
	}
	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
