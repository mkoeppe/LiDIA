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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/quadratic_number_standard.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main (int argc, char*argv[])
{
	bigint                        Delta, h;
	base_vector< bigint > Deltas;
	base_vector< quadratic_number_standard > units;
	xbigfloat        R, L;
	quadratic_order  O;
	quadratic_number_standard u, v, w;
	random_generator rg;

	base_vector< quadratic_number_standard > q;
	base_vector< long > exp;
	base_vector< bigint > bigexp;
	matrix< bigint > M;
	lidia_size_t      i, j;
	lidia_size_t      rows, max_rows;
	lidia_size_t      cols, max_cols;
	long              m, max_exp;
	bigint            g, g1, max_coeff;

	bool             quiet;
	bool             selftest;
	bool             byhand;
	int 		    info;
	lidia_size_t     num_of_tests, k, kk;

	long bR;

	base_vector< quadratic_number_standard > base_u;
	base_vector< bigint > exp_u;
	base_vector< xbigfloat > log_base_u;
	base_vector< long > prec_log_base_u;
	xbigfloat log_u;

	base_vector< xbigfloat > log_q;
	base_vector< long > prec_log_q;

	bigint quot;
	bigint exp_sum;


	//
	//
	// test configuration
	//
	//
	if (argc == 2)
		if (!strcmp(argv[1], "--quiet"))
			quiet = true;
		else
			quiet = false;
	else
		quiet = false;

	byhand = false;
	selftest = true;

	if (quiet)
		info = 0;
	else
		info = 1;

	info = 0;

	max_rows = 10;
	max_cols = 10;
	max_exp = 10;
	max_coeff = 10;
	num_of_tests = 8; // if modified, adjust Deltas and units below

	if (!quiet) {
		std::cout << "Do you want to run the self test ? If not, you" << std::endl;
		std::cout << "have to enter the number of tests and for each" << std::endl;
		std::cout << "test, you must enter a discriminant and, in case" << std::endl;
		std::cout << "of real quadratic fields, a fundamental unit." << std::endl;
		std::cout << std::endl;

		std::cout << "self test (0 = no / 1 = yes) ? ";
		std::cin >> h;
		if (h == 0)
			selftest = false;

		if (!selftest) {
			std::cout << "Please, enter number of tests = " << std::endl;
			std::cin >> num_of_tests;

			std::cout << "Do you want to enter every possible parameter" << std::endl;
			std::cout << "by hand (0 = no / 1 = yes) ?  ";
			std::cin >> h;
			if (h == 1)
				byhand = true;
		}
	}
	else
		selftest = true;

	if (selftest) {
		Deltas.set_capacity(num_of_tests);

		Deltas[0] = -3;
		Deltas[1] = -4;
		Deltas[2] = -500;
		Deltas[3] = -65535;
		Deltas[4] = 5;
		Deltas[5] = 73;
		Deltas[6] = 329;
		Deltas[7] = 493;
	}


	//
	//
	//   the tests
	//
	//

	for (k = 0; k < num_of_tests; k++) {
		if (!selftest) {
			std::cout << "Please, enter quadratic discriminant Delta = ";
			std::cin >> Delta;
			std::cout << std::endl;

			O.assign(Delta);
			u.assign_order(O);
			std::cout << "Please, enter a real quadratic unit u = ";
			std::cin >> u;
		}
		else {
			Delta = Deltas[k];
			O.assign(Deltas[k]);

			if (Deltas[k] > 0) {
				u.assign_order(O);

				switch (k) {
				case 4: u.assign(1, 1, 2);
					break;
				case 5: u.assign(2136, 250, 2);
					break;
				case 6: u.assign(4752830, 262032, 2);
					break;
				case 7: u.assign(111, 5, 2);
					break;
				};
			}
		}


		//
		// check discriminant and unit
		//
		if (Delta > 0)
			assert(abs(norm(u)) == 1);

		//
		// Determine class number
		//
		if (!quiet) {
			std::cout << "Quadratic order O constructed." << std::endl;
			std::cout << "O = " << O << std::endl;
		}

		h = O.class_number();
		if (!quiet) {
			std::cout << "quadratic_order found class number " << h;
			std::cout << std::endl << std::endl;
		}

		//
		// imaginary case
		//
		if (Delta < 0) {
			//
			//
			//  Testing verify_h
			//
			//

			if (!quiet) {
				std::cout << "Testing functions for imaginary quadratic fields";
				std::cout << std::endl;
				std::cout << "================================================";
				std::cout << std::endl << std::endl;
				std::cout << "discriminant = " << Delta << std::endl;
			}
			assert(O.verify_h(Delta, h) == true);
			if (!quiet)
				std::cout << "verify_h() returned true, correct." << std::endl;

			assert(O.verify_h(Delta, 2*h) == false);
			if (!quiet)
				std::cout << "verify_h() returned false, correct." << std::endl;
		}

		//
		// real case
		//
		else {

			for (kk = 0; kk < 10*num_of_tests; kk++) {
				if (!quiet) {
					std::cout << "Testing functions for real quadratic fields";
					std::cout << std::endl;
					std::cout << "===========================================";
					std::cout << std::endl << std::endl;
					std::cout << "discriminant = " << Delta << std::endl;
					std::cout << "unit u = " << u << std::endl;
				}

				//
				//
				//  Testing find_generating_unit
				//
				//

				if (!selftest && !byhand) {
					std::cout << "max. number of rows = ";
					std::cin >> max_rows;
					std::cout << "max. number of cols = ";
					std::cin >> max_cols;
					std::cout << "max. exponent = ";
					std::cin >> max_exp;
					std::cout << "max. coefficient = ";
					std::cin >> max_coeff;
				}

				//
				// generate "rows" powers of u
				//
				if (byhand) {
					std::cout << "enter number of rows = ";
					std::cin >> rows;
				}
				else {
					rg >> rows;
					rows %= max_rows;
					if (rows < 0) rows = -rows;
					if (rows == 0) rows = 1;
				}

				if (!quiet)
					std::cout << "generating " << rows << " powers" << std::endl;
				q.set_capacity(rows);
				exp.set_capacity(rows);

				for (i = 0; i < rows; i++) {
					if (byhand) {
						std::cout << "enter exp[" << i << "] = ";
						std::cin >> exp[i];
					}
					else {
						rg >> exp[i];
						exp[i] %= max_exp;
					}
					power (q[i], u, exp[i]);
					if (!quiet) {
						std::cout << "exponent = " << exp[i] << ", unit = " << q[i];
						std::cout << std::endl;
					}
				}

				//
				// generate a "rows" by "cols" matrix
				// entries bounded by max_coeff
				//
				if (byhand) {
					std::cout << "enter number of cols = ";
					std::cin >> cols;
				}
				else {
					rg >> cols;
					cols %= max_cols;
					if (cols < 0) cols = -cols;
					if (cols == 0) cols = 1;
				}

				M.set_no_of_rows(rows);
				M.set_no_of_columns(cols);

				if (byhand) {
					std::cout << "enter coefficient matrix" << std::endl;
					std::cin >> M;
				}
				else
					M.randomize(max_coeff);

				if (!quiet) {
					std::cout << "generated a " << rows << " by " << cols << " matrix";
					std::cout << std::endl << M << std::endl;
				}

				//
				// determine a lower regulator bound
				//
				if (kk & 1)
					m = O.lower_regulator_bound(Delta);
				else
					m = O.lower_regulator_bound(Delta, h);

				if (!quiet) {
					std::cout << "found lower regulator bound 2^" << m << std::endl;
					if (kk & 1)
						std::cout << "using only Delta.";
					else
						std::cout << "using Delta and class number information.";
					std::cout << std::endl;
				}

				//
				// determine power of generating unit
				//           = gcd( sum_{i=0}^{rows-1} exp_i * M[i][j],
				//                  j = 0,...,cols-1 )
				//
				g = 0;
				for (j = 0; j < cols; j++) {
					g1 = 0;
					for (i = 0; i < rows; i++)
						g1 += exp[i] * M.member(i, j);

					g = gcd(g, g1);
				}
				if (!quiet)
					std::cout << "exponent of generating unit is " << g << std::endl;

				//
				// determine generating unit
				//
				if (!quiet) {
					std::cout << "determing generating unit" << std::endl;
					O.find_generating_unit(M, q, m, 1);
					std::cout << "exponent matrix " << std::endl;
					std::cout << M << std::endl;
				}
				else
					O.find_generating_unit(M, q, m, 0);

				//
				// verifying sum_{i=0}^{rows-1} exp_i * M[i][0] = g
				//
				g1 = 0;
				for (i = 0; i < rows; i++)
					g1 += exp[i] * M.member(i, 0);

				if (!quiet)
					std::cout << "exponent of generating unit (find_generating_unit) = " << g1 << std::endl;

				assert(g == g1);
				if (!quiet)
					std::cout << "find_generating_unit() is correct." << std::endl;

				//
				//
				// Testing reduce_modulo_regulator
				//
				//

				log_base_u.set_capacity(0);
				prec_log_base_u.set_capacity(0);
				base_u.set_capacity(1);
				exp_u.set_capacity(1);

				base_u[0] = u;
				exp_u[0].assign_one();

				std::cout << "fundamental unit u = " << base_u << "^" << exp_u << std::endl;

				O.absolute_Ln_approximation (log_u, log_base_u, prec_log_base_u,
							     exp_u, base_u, -m+1);
				bR = log_u.b_value();

				std::cout << "b(R) +- 1 = " << bR << std::endl;
				std::cout << "R = " << log_u << std::endl;

				log_q.set_capacity(0);
				prec_log_q.set_capacity(0);

				bigexp.set_capacity(exp.get_size());
				v.assign_order(O);
				v.assign_one();
				for (i = 0; i < rows; i++) {
					bigexp[i] = exp[i];
					power(w, q[i], exp[i]);
					multiply (v, v, w);
				}

				O.reduce_modulo_regulator (quot, q, bigexp, log_q, prec_log_q,
							   base_u, exp_u, log_base_u, prec_log_base_u, bR);

				exp_sum = bigexp[0] * bigexp[0];
				for (i = 1; i < bigexp.get_size(); i++)
					add (exp_sum, exp_sum, bigexp[i]*bigexp[i]);

				std::cout << "v = " << v << std::endl;
				std::cout << "sum of exponents = " << exp_sum << std::endl;
				std::cout << "quotient = " << quot << std::endl;

				assert(quot == exp_sum || quot == exp_sum-1 || quot == exp_sum+1 || (quot == 0 && exp_sum < 32));
			}

			//
			//
			//   Testing verify_hR
			//
			//

			// fundamental unit as implicit power product
			q.set_capacity(1);
			q[0] = u;
			bigexp.set_capacity(1);
			bigexp[0].assign_one();

			// lower regulator bound
			L.assign_zero();
			m = O.lower_regulator_bound(L, Delta, h);

			if (!quiet) {
				std::cout << "found lower regulator bound 2^" << m << std::endl;
				std::cout << "L = " << L << std::endl;
			}

			// relative 7-approximation to the regulator
			O.relative_Ln_approximation(R, q, bigexp, 7, m);

			// class number verification
			assert(O.verify_hR(Delta, h, R, info) == true);
			if (!quiet)
				std::cout << "verify_hR(D, h, R) returned true, correct." << std::endl;

			assert(O.verify_hR(L, Delta, h, R, info) == true);
			if (!quiet)
				std::cout << "verify_hR(L, D, h, R) returned true, correct." << std::endl;

			assert(O.verify_hR(Delta, 2*h, R, info) == false);
			if (!quiet)
				std::cout << "verify_hR(D, h, R) returned true, correct." << std::endl;

			assert(O.verify_hR(L, Delta, 2*h, R, info) == false);
			if (!quiet)
				std::cout << "verify_hR() returned false, correct." << std::endl;
		}
	}

	return 0;
}
