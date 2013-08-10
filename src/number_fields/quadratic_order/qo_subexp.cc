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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/lidia_signal.h"
#include	"LiDIA/xbigfloat.h"
#include	"LiDIA/base/b_value.h"
#include	"LiDIA/qi_class.h"
#include	"LiDIA/qi_class_real.h"
#include	"LiDIA/number_fields/qo_list.h"
#include	<cstdio>
#include	<unistd.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef __EMX__
#define lidia_error_handler_sub(f, m)\
{ fcloseall(); \
  char command[50]; \
  sprintf (command, "rm -f /tmp/%ld*QO*", static_cast<long>(getpid())); \
  std::system(command); \
  sprintf (command, "rm -f %ld*QO*", static_cast<long>(getpid())); \
  std::system(command); \
  lidia_error_handler(f, m); }
#else
#define lidia_error_handler_sub(f, m)\
{ char command[50]; \
  sprintf (command, "rm -f /tmp/%ldQO*", static_cast<long>(getpid())); \
  std::system(command); \
  sprintf (command, "rm -f %ldQO*", static_cast<long>(getpid())); \
  std::system(command); \
  lidia_error_handler(f, m); }
#endif


//
// Parameter list
//

const float qo_params[73][7] =
{
	{15, 0.8, 1200, 60, 3, 14, 2},
	{16, 0.8, 1400, 70, 3, 14, 2},
	{17, 0.8, 3000, 80, 3, 14, 2},
	{18, 0.8, 3000, 90, 3, 14, 2},
	{19, 0.8, 3600, 100, 3, 15, 2},
	{20, 0.8, 4000, 110, 3, 15, 2},
	{21, 0.8, 4250, 120, 3, 15, 2},
	{22, 0.8, 4500, 130, 3, 15, 3},
	{23, 0.8, 4750, 140, 4, 16, 3},
	{24, 0.8, 5000, 150, 4, 16, 4},
	{25, 0.8, 5000, 160, 4, 16, 4},
	{26, 0.9, 6000, 170, 5, 16, 4},
	{27, 1.17, 6000, 180, 5, 16, 5},
	{28, 1.17, 6500, 190, 5, 16, 5},
	{29, 1.17, 6500, 200, 5, 16, 5},
	{30, 1.36, 7000, 210, 5, 16, 5},
	{31, 1.36, 7000, 220, 5, 16, 5},
	{32, 1.36, 7500, 230, 5, 16, 5},
	{33, 1.43, 7500, 240, 6, 16, 6},
	{34, 1.43, 7500, 250, 6, 16, 6},
	{35, 1.43, 7500, 260, 6, 18, 6},
	{36, 1.43, 8000, 270, 6, 18, 6},
	{37, 1.92, 9000, 290, 6, 18, 7},
	{38, 1.92, 10000, 300, 6, 18, 7},
	{39, 1.92, 11000, 310, 6, 18, 7},
	{40, 1.72, 14000, 340, 6, 18, 7},
	{41, 1.72, 15500, 360, 6, 18, 8},
	{42, 1.69, 12000, 390, 6, 18, 8},
	{43, 1.69, 12000, 430, 6, 18, 8},
	{44, 1.69, 15000, 490, 7, 18, 9},
	{45, 1.69, 15000, 520, 7, 18, 9},
	{46, 1.69, 15000, 550, 7, 18, 9},
	{47, 1.69, 25000, 590, 7, 18, 10},
	{48, 1.69, 25000, 620, 7, 18, 10},
	{49, 1.69, 30000, 660, 7, 18, 10},
	{50, 1.75, 50000, 730, 7, 19, 10},
	{51, 1.75, 50000, 750, 7, 19, 10},
	{52, 1.75, 50000, 900, 7, 19, 11},
	{53, 1.85, 50000, 950, 7, 19, 11},
	{54, 1.95, 80000, 1000, 7, 20, 11},
	{55, 1.95, 80000, 1100, 7, 20, 11},
	{56, 1.95, 80000, 1200, 7, 20, 11},
	{57, 1.95, 80000, 1300, 8, 20, 12},
	{58, 1.95, 80000, 1400, 8, 20, 12},
	{59, 2.15, 90000, 1500, 8, 20, 12},
	{60, 2.15, 90000, 1600, 8, 20, 12},
	{61, 2.15, 100000, 1700, 8, 21, 13},
	{62, 2.15, 100000, 1900, 8, 21, 13},
	{63, 2.35, 110000, 2100, 8, 21, 13},
	{64, 2.35, 110000, 2300, 8, 21, 13},
	{65, 2.35, 120000, 2500, 8, 21, 13},
	{66, 2.35, 120000, 2700, 8, 21, 13},
	{67, 2.4, 130000, 3000, 8, 21, 13},
	{68, 2.4, 130000, 3300, 8, 22, 13},
	{69, 2.4, 140000, 3600, 8, 22, 13},
	{70, 2.45, 140000, 4000, 8, 22, 13},
	{71, 2.45, 150000, 4200, 9, 22, 13},
	{72, 2.4, 150000, 4400, 9, 23, 13},
	{73, 2.4, 160000, 4600, 9, 23, 13},
	{74, 2.4, 160000, 4800, 9, 24, 13},
	{75, 2.6, 170000, 5000, 9, 24, 13},
	{76, 2.6, 170000, 5200, 9, 25, 13},
	{77, 2.6, 180000, 5400, 9, 25, 13},
	{78, 2.6, 200000, 5600, 9, 26, 13},
	{79, 2.65, 220000, 5800, 9, 26, 13},
	{80, 2.65, 250000, 6000, 9, 27, 13}};



LIDIA_SIGNAL_FUNCTION(qo_stop_it)
{
	lidia_error_handler_sub("quadratic_order", "class_group - the "
				"program was interrupted");
}



// FOR TABLE GENERATION
bool quadratic_order::special = false;
bool quadratic_order::noCL = false;
char quadratic_order::matfile[100];
bool firstcont = true;
bigint hmult;
timer trel, thnf, thnf2, treg, tinit, ttot;
long T1, T2, T3, T4, T5, T6;
int P1 = 0, P2 = 0, P3 = 0, P4 = 0, P5 = 0, P6 = 0, P7 = 0, P8 = 0;



          /////////////////////////////////////////////////////////
          //                                                     //
          // GENERAL FUNCTIONS FOR ALL SUBEXPONENTIAL ALGORITHMS //
          //                                                     //
          /////////////////////////////////////////////////////////

//
// quadratic_order::rBJT_bound
//
// Task:
//      computes R using the BJT algorithm, if R < bound
//

void
quadratic_order::rBJT_bound(const bigfloat & bound)
{
	debug_handler("quadratic_order", "rBJT_bound");

	bigfloat y, usqr, u, temp;
	long upper;
	bigint oa, temp2;
	qi_class_real A, C, D, *F, UU;

	UU.assign_one();

	// get hash table size
	multiply(temp, bound, bigfloat(1.33));
	temp.longify(upper);
	prin_list.initialize(static_cast<lidia_size_t>(upper));
	prin_list.set_key_function(&qi_class_real_key);

	// compute initial step-width = ceil((log(Delta)+3)/4)
	add(temp, log(bigfloat(Delta)), 3);
	shift_right(temp, temp, 2);
	ceil(temp, temp);
	temp.bigintify(temp2);
	if (temp2.is_odd())
		inc(temp2);
	u.assign(temp2);

	A.assign_one();

	D.assign_one();
	y.assign(u);

	// oa = (D - b^2) / 4a
	oa.assign(abs(A.get_c()));

	while ((u <= bound) && (R.is_zero())) {
		// compute more baby steps
		R = A.make_list(oa, u, prin_list);
		if (R.compare(min_log) < 0)
			R.assign_zero();

		// compute giant steps to u^2
		if (R.is_zero()) {
			inverse(C, A);
			C.adjust_abs(u);
			square(usqr, u);
		}
		while ((y.compare(usqr) < 0) && (R.is_zero())) {
			multiply_real(D, D, C);
			D.adjust_abs(y);

			F = prin_list.search(D);
			if (F) {
				// found D in list:  R = y-r
				subtract(R, D.get_distance(), F->get_distance());
				R.absolute_value();
				if (R.compare(min_log) < 0)
					R.assign_zero();
			}

			add(y, y, u);
		}

		u.multiply_by_2();
	}
}



//
// quadratic_order::qo_get_fb
//
// Task:
//      compute factor base for sieving versions of subexponential algorithms
//

void
quadratic_order::qo_get_fb(int & size_FB)
{
	debug_handler("quadratic_order", "qo_get_fb");

	qi_class PID;
	bigfloat temp;
	int bach, upper;
	register int P, n = 0;

	// bound = 0.7 log^2(Delta)
	log(temp, bigfloat(abs(Delta)));
	square(temp, temp);
	if (Delta.is_lt_zero())
		multiply(temp, temp, bigfloat(0.6));
	else
		multiply(temp, temp, bigfloat(0.7));
	floor(temp, temp);
	temp.intify(bach);


	// if prime ideal over 2 exists, add it to the factor base
	if (generate_prime_ideal(PID, bigint(2)))
		fact_base[n++] = PID;


	// add primes until we have size_FB of them and the largest has norm greater
	// then 0.6 log^2 Delta

	P = 3;
	upper = 3+10*(size_FB-n);
	if (upper > 1000000)
		upper = 1000000;
	if ((static_cast<int>(PL.get_upper_bound()) < upper) || (PL.get_lower_bound() != 2))
		PL.resize(2, upper);

	P = PL.get_first_prime(P);
	while ((n < size_FB) || (P < bach)) {
		if (generate_prime_ideal(PID, P)) {
			if (PID.get_a() == bigint(P))
				fact_base[n++] = PID;
			else
				break;
		}

		P = PL.get_next_prime();
		if (P == 0) {
			PL.resize(PL.get_upper_bound() + 1, PL.get_upper_bound() + 1000000);
			P = PL.get_first_prime();
		}
	}

	size_FB = fact_base.size();
}



//
// generate_free_relations
//
// Task:
//      generates relations corresponding to the divisors of Delta in the
//      factor base, and returns the number that were found
//

void
quadratic_order::generate_free_relations(lidia_size_t & curr_A,
					 sort_vector< lidia_size_t > & S)
{
	debug_handler("quadratic_order", "generate_free_relations");

	register lidia_size_t i, j, size_FB;
	int p;

	size_FB = fact_base.size();
	for (i = 0; i < size_FB; ++i) {
		fact_base[i].get_a().intify(p);
		if (!remainder(Delta, p)) {
			AMAT.sto(i, curr_A, 2);
			minima_qn[curr_A].assign(1, 0, p);

			// update index set
			for (j = i; j < S.size()-1; ++j)
				S[j] = S[j+1];
			S.set_size(S.size()-1);

			++curr_A;
		}
	}
}



//
// qo_hnf_and_det
//
// Task:
//      computes the Hermite normal form of the relation matrix, its
//      determinant, and an approximation of the regulator (both bigfloat
//      and xbigfloat)
//

void
quadratic_order::qo_hnf_and_det(matrix< bigint > & A1,
                                sort_vector< lidia_size_t > & S,
                                xbigfloat & LL,
                                bool input_reg)
{
	debug_handler("quadratic_order", "qo_hnf_and_det");

	register lidia_size_t i, j;
	lidia_size_t A1cols;
	int strat = 0, *bad_cols = NULL;
	long lowregbound;
	bool real_order = is_real();
	bigint temp;
	timer t;
	bool firsttime = false;
	static matrix< bigint > kern;
	static bool BSGS_R, do_mod;

	A1cols = A1.get_no_of_columns();

	// select transformation matrix strategy
	if (real_order) {
		bigfloat::set_precision(prec);

		if (A1cols)
			i = A1.get_no_of_rows();
		else
			i = AMAT.get_no_of_rows();

		if (i >= 2000)
			TR.set_mode(1);
		if (i >= 2500)
			TR.set_mode(2);
		if (i > 3000)
			strat = 1;
	}

	if (real_order && (A1cols == 0)) {
		// attempt to compute regulator using Babystep Giantstep

		BSGS_R = false;

		if (info > 1) {
			std::cout << "\nAttempting to compute R with BS-GS..." << std::endl;
			t.start_timer();
		}
		if (special)
			treg.start_timer();

		if (decimal_length(Delta) <= 30)
			rBJT_bound(bigfloat(50.0));
		else if (decimal_length(Delta) <= 60)
			rBJT_bound(bigfloat(100.0));
		else
			rBJT_bound(bigfloat(1000.0));

		if (input_reg) {
			std::cout << "Input the regulator..." << std::endl;
			std::cin >> R;
		}

		if (!R.is_zero()) {
			BSGS_R = true;
			R_xbig.assign(R);
			if (info > 1)
				std::cout << "Found R!  R = " << R << std::endl;
			if (input_reg)
				prin_list.empty();
		}
		else {
			prin_list.empty();
		}

		if (info > 1) {
			t.stop_timer();
			std::cout << "\nBS-GS time:  ";
			MyTime(t.user_time());
			std::cout << "\n\n";
		}
		if (special)
			treg.stop_timer();
	}


	//
	// compute HNF
	//
	// store the HNF in RHNF and
	// the transformation matrix in TR.
	//

///////////// INFO /////////////////////
	// FOR TABLE GENERATION
	if (special) {
		// open matrix file and output, if necessary
		if ((A1cols == 0) && (strlen(matfile) > 0)) {
			std::ofstream out2;
			out2.open(matfile);
			AMAT.set_print_mode(LIDIA_MODE);
			out2 << AMAT << std::endl;
			out2.close();
		}

		if (noCL) {
			ttot.stop_timer();

			hmult.assign_zero();

			T1 = tinit.user_time();
			T2 = trel.user_time();
			T3 = thnf.user_time();
			T4 = thnf2.user_time();
			T5 = treg.user_time();
			T6 = ttot.user_time();

			std::cout << AMAT.get_no_of_rows() << " 0 0 0 0";
			if (is_real())
				std::cout << " 0" << std::endl;
			else
				std::cout << std::endl;

			std::cout << hmult << std::endl;


			sieve.dense_bound(P2);
			std::cout << P1 << "  " << P2 << "  " << P3 << "  " << P4 << "  ";
			std::cout << P5 << "  " << P6 << "  " << P7 << "  " << P8 << std::endl;

			std::cout << T1 << "  " << T2 << "  " << T3 << "  " << T4 << "  ";
			std::cout << T5 << "  " << T6 << std::endl;

			return;
		}


		firsttime = false;
		if (A1cols == 0) {
			thnf.start_timer();
			firstcont = true;
		}
		else {
			if (do_mod) {
				thnf.cont_timer();
				firsttime = true;
			}
			else if (firstcont) {
				firstcont = false;
				thnf2.start_timer();
			}
			else
				thnf2.cont_timer();
		}
	}


	if (info > 1) {
		if (A1cols == 0) {
			std::cout << "\nComputing HNF of " << AMAT.get_no_of_rows() << " x ";
			std::cout << AMAT.get_no_of_columns() << " matrix:" << std::endl;
		}
		else {
			std::cout << "\nComputing HNF of " << A1.get_no_of_rows() << " x ";
			std::cout << A1.get_no_of_columns() << " matrix:" << std::endl;
		}
		t.start_timer();
	}
///////////// INFO /////////////////////

#if 0
	std::cout << "A = " << std::endl;
	AMAT.set_print_mode(LIDIA_MODE);
	std::cout << AMAT << std::endl;
#endif

	TR.set_no_of_columns(AMAT.get_no_of_columns());
	TR.touch_files();

	if (A1cols == 0) {
#ifdef LIDIA_DEBUG
		if (Delta.is_lt_zero()) {
			qi_class TT, PP;
			for (j = 0; j < AMAT.get_no_of_columns(); ++j) {
				TT.assign_one();
				for (i = 0; i < fact_base.size(); ++i) {
					if (AMAT.member(i, j) != 0) {
						power(PP, fact_base[i], AMAT.member(i, j));
						multiply(TT, TT, PP);
					}
				}
				if (!TT.is_one())
					std::cout << "BAD REL (j = " << j << "):  TT = " << TT << std::endl;
			}
		}
		else {
			quadratic_ideal TT, PP, TT2;
			for (j = 0; j < AMAT.get_no_of_columns(); ++j) {
				TT.assign_one(*this);
				for (i = 0; i < fact_base.size(); ++i) {
					if (AMAT.member(i, j) != 0) {
						power(PP, quadratic_ideal(fact_base[i]), AMAT.member(i, j));
						multiply(TT, TT, PP);
					}
				}
//        TT2.assign_principal(minima_qn[j]);
				TT2.assign_one(*this);
				multiply(TT2, TT2, minima_qn[j]);

				if (TT != TT2) {
					std::cout << "\nBAD REL (j = " << j << "):  TT = " << TT;
					std::cout << ", TT2 = " << TT2;
					std::cout << ", minima_qn[j] = " << minima_qn[j] << std::endl;
					std::cout << "TT1.norm = " << TT.norm();
					std::cout << ", TT2.norm = " << TT2.norm();
					std::cout << ", min.norm = " << minima_qn[j].norm() << std::endl;
				}
			}
		}
#endif

		kern.reset();
		if (real_order && !BSGS_R) {
			bad_cols = RHNF.hnf_cg2(AMAT, TR, 5, DET, strat, do_mod);
		}
		else {
			bad_cols = RHNF.hnf_cg2(AMAT, 5, do_mod);
			if (real_order && special)
				std::cout << " 0" << std::endl;
		}
	}
	else {
		if (real_order && !BSGS_R) {
			if (do_mod)
				bad_cols = RHNF.hnf_cg3_mod(A1, TR, strat, do_mod);
			else
				bad_cols = RHNF.hnf_cg3(A1, TR);
		}
		else {
			if (do_mod)
				bad_cols = RHNF.hnf_cg3_mod(A1, do_mod);
			else
				bad_cols = RHNF.hnf_cg3(A1);
		}
	}

#if 0
	std::cout << "HNF(A) = " << std::endl;
	RHNF.set_print_mode(LIDIA_MODE);
	std::cout << "\n\n" << RHNF << std::endl;
	exit(1);
#endif


	// store singularities
	if (bad_cols) {
		S.set_size(bad_cols[0]);
		for (i = 0; i < bad_cols[0]; ++i)
			S[i] = bad_cols[i+1];
		delete [] bad_cols;
	}


	// remove zero columns
	bad_cols = new int[RHNF.get_no_of_columns()];
	bool zero = true;
	bad_cols[0] = 0;
	for (i = 0; i < RHNF.get_no_of_columns() && zero; ++i) {
		for (j = 0; j < RHNF.get_no_of_rows() && zero; ++j) {
			if (RHNF.member(j, i) != 0)
				zero = false;
		}
		if (zero) {
			bad_cols[i+1] = i;
			++bad_cols[0];
		}
	}
	RHNF.remove_columns(bad_cols);


	// compute determinant
	if (RHNF.get_no_of_rows() == RHNF.get_no_of_columns()) {
		h.assign_one();
		for (i = 0; i < RHNF.get_no_of_rows(); ++i)
			multiply(h, h, RHNF.member(i, i));
	}
	else
		h.assign_zero();


///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\nHNF done..." << std::endl;
		std::cout << "HNF time:  ";
		MyTime(t.user_time());
		std::cout << std::endl;
		std::cout << "\nMatrix after HNF:  " << RHNF.get_no_of_rows() << " x ";
		std::cout << RHNF.get_no_of_columns() << std::endl;
		std::cout << "h = " << h << std::endl;
	}

	// FOR TABLE GENERATION
	if (special) {
		if ((A1cols == 0) || firsttime)
			thnf.stop_timer();
		else
			thnf2.stop_timer();
	}
///////////// INFO /////////////////////

	// compute kernel elements and regulator, if necessary
	if (real_order && !do_mod) {
		// FOR TABLE GENERATION
		if (special) {
			treg.cont_timer();
		}

		if (!do_mod && !BSGS_R && (kern.get_no_of_columns() == 0)) {
			if (info > 1) {
				t.start_timer();
				std::cout << "\nComputing kernel elements" << std::endl;
				std::cout << TR.get_size() << " transformation matrices!!!" << std::endl;
			}

			TR.set_no_of_columns(minima_qn.size());

			i = bad_cols[0];
			if ((decimal_length(Delta) > 24) && (i > 30))
				bad_cols[0] = 25;
			TR.get_submatrix(kern, bad_cols);
			bad_cols[0] = i;

			// factor out gcd's
			if (DET != 1) {
				if (info > 3) {
					std::cout << "DET = " << DET << std::endl;
				}
				bigint G;
				for (j = 0; j < kern.get_no_of_columns(); ++j) {
					G = kern.member(0, j);
					for (i = 1; i < kern.get_no_of_rows(); ++i)
						G = gcd(G, kern.member(i, j));
					if (DET.is_lt_zero())
						G.negate();
					if (abs(G) > 1) {
						for (i = 0; i < kern.get_no_of_rows(); ++i)
							kern.sto(i, j, kern.member(i, j)/G);
					}

					if (info > 3) {
						std::cout << "col " << j << ", G = " << G << std::endl;
					}
				}
			}

			if (info > 1) {
				t.stop_timer();
				std::cout << "Kernel done..." << std::endl;
				std::cout << "Kernel time:  ";
				MyTime(t.user_time());
				std::cout << std::endl;
			}
		}


#if 0
		// remove kernel elements from transformation matrix
		if (!BSGS_R && !do_mod)
			TR.remove_columns(bad_cols);
#endif

		// if h is non-zero, compute regulator (if not already computed!)
		if (true) {
			//if (!h.is_zero() && R.is_zero() && !BSGS_R && !do_mod) {
			if (info > 1) {
				t.start_timer();
			}

			kern.set_no_of_rows(minima_qn.size());

			if (info > 1) {
				std::cout << "\nMAX in kernel = " << decimal_length(kern.max());
				std::cout << " decimal digits" << std::endl;
				std::cout << "Determing lower regulator bound ... " << std::endl;
			}

			//lowregbound = lower_regulator_bound(LL,Delta,h,info);

// MJJ
			if (h.is_zero()) {
				h.assign_one();
				lowregbound = lower_regulator_bound(LL, Delta, h, info);
				h.assign_zero();
			}
			else
				lowregbound = lower_regulator_bound(LL, Delta, h, info);
// MJJ



			if (info > 1) {
				t.stop_timer();
				std::cout << "done" << std::endl;
				std::cout << "lower regulator bound is 2^" << lowregbound << std::endl;
				std::cout << std::endl;
				std::cout << "elapsed time:  ";
				MyTime(t.user_time());
				std::cout << std::endl;
				std::cout << "Determing generating unit ... ";
				t.cont_timer();
			}

			kern.set_print_mode(1);
			fund_unit.generating_unit (kern, minima_qn, lowregbound);

			// <MM> FIX ME
			// We should work with a compact representation of the fundamental
			// unit, because the shorter representation allows a faster
			// approximation of the regulator. But then the basis would not be
			// the vector of minima_qn anymore, which would imply a malfunction
			// of is_in_lattice at the moment.
			//
			// fund_unit.compact_representation_of_unit();
			// </MM>


			if (info > 1) {
				t.stop_timer();
				std::cout << "done" << std::endl;
				std::cout << "elapsed time:  ";
				MyTime(t.user_time());
				std::cout << std::endl;
				t.cont_timer();
			}

			// <MM> FIX ME
			// In the following, an approximation R to the regulator is
			// refined until nearest finds the order again.
			// This loop could be replaced by a call of
			//
			// R_xbig = fund_unit.get_absolute_Ln_approximation(4);
			//
			// because an absolute 4 approximation is sufficient to
			// compute the fundamental unit.
			// </MM>


			// determine regulator approximation (should be removed)
			long e_tmp;

			bigfloat::set_precision(prec);
			fund_unit.get_relative_Ln_approximation(R_xbig, xprec, lowregbound);

			R_xbig.absolute_value();
			R.assign(R_xbig.get_mantissa());
			e_tmp = R_xbig.get_exponent()-R_xbig.get_mantissa().bit_length();
			if (e_tmp > 0)
				shift_left(R, R, e_tmp);
			else
				shift_right(R, R, -e_tmp);

			if (info > 3) {
				std::cout << "mantissa:  " << R_xbig.get_mantissa() << "\n";
				std::cout << "exp:       " << R_xbig.get_exponent() << "\n";
				std::cout << "e_tmp:     " << e_tmp << std::endl;
			}

			if (info > 1) {
				std::cout << "\nprecision:  " << bigfloat::get_precision() << std::endl;
				std::cout << "relative prec:  " << xprec << std::endl;
				std::cout << "regulator approx. (xbigfloat) = " << R_xbig << std::endl;
				std::cout << "regulator approx. (bigfloat) = " << R << std::endl;
				std::cout << std::endl;
			}

			//
			// Verify class number by computing relative
			// 7-approximation to presumable regulator and
			// calling verify_hR. Approximation L to L(1,chi_D)
			// computed for lower regulator bound is reused.
			// Generating unit is represented as an implicit power
			// product. Base elements are stored in minima_qn and
			// exponent vector is first column of tr_lat.
			//

			if (info > 1) {
				std::cout << "Computing relative 7-approximation to presumable ";
				std::cout << "regulator ...";
			}

			//<MM>
			fund_unit.get_relative_Ln_approximation(R_xbig, 7, lowregbound);
			//</MM>
			R_xbig.absolute_value();

			if (info > 1) {
				t.stop_timer();
				std::cout << "done" << std::endl << "(" << R_xbig << ")" << std::endl;
				std::cout << "Regulator time:  ";
				MyTime(t.user_time());
				std::cout << std::endl;
			}

			// factor out multiples of the regulator, if necessary
			if (!verify_hR(LL, Delta, h, R_xbig)) {
				if (info > 1) {
					t.start_timer();
					std::cout << "Attempting to factor out small factor..." << std::endl;
				}

				// <MM> FIX ME
				// This call replaces the code for factoring out multiples
				// of the regulator. But the result of the function is
				// a compact representation of the fundamental unit, of which
				// basis is not the vector of minima_qn anymore. Hence,
				// at the moment the function is_in_lattice would not work
				// properly.
				//
				// fund_unit.fundamental_unit( fund_unit, h );
				// </MM>

				qi_class_real Freal;
				bigint fac;
				bigfloat R2, hR, hR2;

				bigfloat hstar = estimate_L1(get_optimal_Q());
				multiply(hstar, hstar, sqrt(bigfloat(Delta/2)));


				multiply(hR, R, bigfloat(h));
				hR2 = ceil(hR/(hstar/2));
				hR2.bigintify(fac);
				while ((hR/fac) > hstar/2)
					++fac;
				--fac;

				while (((R/fac)< hstar) && (fac > 1)) {
					divide(R2, R, bigfloat(fac));
					if (info > 2)
						std::cout << "fac = " << fac << ", R/fac = " << R2 << std::endl;

					Freal.assign_one();
					Freal.assign(nearest(Freal, R2));

					if (info > 3)
						std::cout << "Freal = " << Freal << std::endl;

					if ((Freal.is_one()) && (Freal.get_distance() > 0)) {
						R.assign(Freal.get_distance());
						R_xbig.assign(R);
						fu_mult.assign(fac);

						break;
					}

					--fac;
				}

				if (info > 1) {
					t.stop_timer();
					std::cout << "R = " << R << std::endl;
					std::cout << "Regulator factor time:  ";
					MyTime(t.user_time());
					std::cout << std::endl;
				}
			}

			// FOR TABLE GENERATION
			if (special) {
				treg.stop_timer();
			}

			//<MM>
			// R_xbig = fund_unit.get_absolute_Ln_approximation(4);
			// R_xbig_prec = 4;
			//</MM>

			// if failed to compute the regulator, compute new kernel vectors
			if (R.is_zero())
				kern.reset();
		}
	}

	delete [] bad_cols;
}




/////////////////////////////////
//                             //
// SUBEXPONENTIAL ALGORITHM #1 //
//   - Duellmann's algorithm   //
//                             //
/////////////////////////////////


//
// qo_read_parameters
//
// Task:
//      read the parameters for the non-sieving versions of the subexponential
//      class group algorithm from a table.
//

void
quadratic_order::qo_read_parameters(int decimal_digits, int &size_FB)
{
	debug_handler("quadratic_order", "qo_read_parameters(non-sieving)");

	int i = -1, r8;

	if (decimal_digits < 15)
		i = 0;
	else if (decimal_digits > 80)
		i = 65;
	else
		i = decimal_digits - 15;

	size_FB = static_cast<int>(qo_params[i][3]); // size of factor base

	// smaller factor bases for real orders
	if (is_real() && decimal_digits <= 50)
		size_FB -= 500;


	// use large factor bases for really small class numbers
	bigfloat L1 = estimate_L1(1000);


	// use small factor bases for really big class numbers
	r8 = remainder(Delta, 8);
	if (r8 < 0)
		r8 += 8;
	if (r8 != 1)
		L1.multiply_by_2();
	if (L1 > 9.0) {
		size_FB >>= 1;
		if (is_real() && (decimal_digits > 65))
			size_FB += 750;
		else if (is_real() && (decimal_digits > 50))
			size_FB += 250;
	}
	else if (L1 < 0.4) {
		if (is_real() && decimal_digits <= 50)
			size_FB += 500;

		if ((size_FB > 300) && (size_FB < 500))
			size_FB += size_FB + (size_FB >> 1);
		else
			size_FB <<= 1;
	}
	else if (is_real() && (size_FB > 3500))
		size_FB -= 500;


	// max FB of 7000
	if (is_real()) {
		if (size_FB > 4000)
			size_FB = 4000;
	}
	else {
		if (size_FB > 7000)
			size_FB = 7000;
	}
}




//
// quadratic_order::qo_init_powers()
//
// Task:
//      initializes the table of small powers of factor base elements needed
//      for Duellmann's algorithm.
//

void
quadratic_order::qo_init_powers(int & k0, qi_class ***powers,
                                quadratic_number_standard ***mins,
                                base_vector< lidia_size_t > & power_index)
{
	debug_handler("quadratic_order", "qo_init_powers");

	register lidia_size_t i, j, curr;
	bigfloat lDelta;
	bigint maxP;
	qi_class C, B;
	quadratic_number_standard alpha, beta;
	lidia_size_t Bound = 30;
	bool real_order = is_real();

	power_index.set_mode(EXPAND);

	// compute k0
	lDelta = ceil(log(bigfloat(abs(Delta))));
	lDelta.bigintify(maxP);
	for (k0 = 0; (k0 < fact_base.size() && (fact_base[k0].get_a() < maxP)); ++k0);
	--k0;
	if (k0 < 4)
		k0 = 4;

	// allocate storage
	(*powers) = new qi_class *[k0];
	if (real_order)
		(*mins) = new quadratic_number_standard *[k0];
	for (i = 0; i < k0; ++i) {
		(*powers)[i] = new qi_class[Bound+1];
		if (real_order)
			(*mins)[i] = new quadratic_number_standard[Bound+1];
	}

	// initialize the array
	curr = 0;
	for (i = 0; (curr < k0) && (i < fact_base.size()); ++i) {
		C = fact_base[i];
		if ((Delta % C.get_a()) != 0) {
			power_index[curr] = i;
			(*powers)[curr][1].assign(C);
			if (real_order) {
				alpha.assign_one();
				beta.assign_one();
				(*mins)[curr][1].assign(alpha);
			}
			B.assign(C);
			for (j = 2; j <= Bound; ++j) {
				if (real_order) {
					multiply_real(B, B, C, alpha, alpha, beta);
					(*powers)[curr][j].assign(B);
					(*mins)[curr][j].assign(alpha);
				}
				else {
					multiply(B, B, C);
					(*powers)[curr][j].assign(B);
				}
			}
			++curr;
		}
	}

	k0 = curr-1;
}



//
// quadratic_order::factor_over_FB()
//
// Task:
//      attempts to compute a vector v such that A = FB^v, i.e., factor
//      A over the factor base.
//

bool
quadratic_order::factor_over_FB(qi_class & AA,
                                math_vector< lidia_size_t > &v)
{
	debug_handler("quadratic_order", "factor_over_FB");

	register lidia_size_t i;
	bigint a, b, temp;
	long p, p2, rem, B2p, Bpr2p;
	lidia_size_t ex;

	for (i = 0; i < fact_base.size(); ++i)
		v[i] = 0;

	a.assign(AA.get_a());
	i = 0;
	while ((a > 1) && (i < fact_base.size())) {
		ex = 0;
		fact_base[i].get_a().longify(p);

		div_rem(temp, rem, a, static_cast<long>(p));
		while (rem == 0) {
			++ex;
			a.assign(temp);
			div_rem(temp, rem, a, static_cast<long>(p));
		}

		// compute sign of exponent and add to v
		if (ex > 0) {
			p2 = (p << 1);
			B2p = remainder(AA.get_b(), p2);
			if (B2p < 0)
				B2p += p2;

			b.assign(fact_base[i].get_b());
			Bpr2p = remainder(b, p2);
			if (Bpr2p < 0)
				Bpr2p += p2;

			if (B2p != Bpr2p)
				ex = -ex;

			// if p divides Delta, take mod 2
			if (remainder(Delta, p) == 0) {
				if (ex & 1)
					ex = 1;
				else
					ex = 0;
			}

			v[i] = ex;
		}

		++i;
	}

	return (a == 1);
}




//
// quadratic_order::qo_relations_randexp()
//
// Task:
//      compute n new relations with Duellmann's method, ensuring that each
//      row of the relation matrix corresponding to the indices in S contains
//      a non-zero entry
//

void
quadratic_order::qo_relations_randexp(int size_FB, int n, qi_class **powers,
                                      quadratic_number_standard **mins,
                                      base_vector< lidia_size_t > & power_index,
                                      sort_vector< lidia_size_t > & S,
                                      matrix< bigint > & A1,
                                      lidia_size_t & curr_A,
                                      lidia_size_t & curr_A1)
{
	debug_handler("quadratic_order", "qo_relations_randexp");

	register lidia_size_t i, j;
	lidia_size_t Bound = 30, Acols, A1cols, newrels, ei;
	qi_class B, BO;
	quadratic_number_standard alpha, beta;
	bigint oa;
	math_vector< lidia_size_t > e(size_FB, size_FB), v(size_FB, size_FB);
	math_vector< lidia_size_t > ee(size_FB, size_FB);
	random_generator rg;
	bool real_order = is_real();
	unsigned int percent = 0;
	int count = 0;
	static lidia_size_t numrels = 0, numvec = 0;
	static unsigned int per_level = 10;
	static int times_called = 0;

	Acols = AMAT.get_no_of_columns();
	A1cols = A1.get_no_of_columns();

	newrels = 0;
	if (A1cols == 0) {
		per_level = 10;
		numvec = 0;
		newrels = curr_A;
		numrels = 0;
		times_called = 0;
	}

	++times_called;
	if ((times_called > 5) && (S.size() == 0)) {
		rg >> S[0];
		S[0] %= size_FB;
		if (S[0] < 0)
			S[0] += size_FB;
		if (S[0]< (size_FB >> 1))
			S[0] += (size_FB >> 1);
		S[0] %= size_FB;
	}


///////////// INFO /////////////////////
	timer t;
	if (info > 1)
		t.start_timer();

	// FOR TABLE GENERATION
	if (special) {
		if (A1cols == 0)
			trel.start_timer();
		else
			trel.cont_timer();
	}
///////////// INFO /////////////////////

	// continue to generate relations until S is empty and i>=n
	B.assign_one();

	while ((newrels< n) || (S.size() > 0)) {
		// compute random exponent vector and reduced ideal equivalent to its
		// corresponding power product over the factor base

		++numvec;

		if (real_order && (count)) {
			// if order is real, use reduction steps to generate new ideal
			e = ee;
		}
		else {
			for (i = 0; i < size_FB; ++i)
				e[i] = 0;

			// first dense part
			B.assign_one();
			alpha.assign_one();
			for (i = 0; i < power_index.size(); ++i) {
				rg >> ei;
				ei = (ei % Bound) + 1;

				e[power_index[i]] = ei;
				if (real_order)
					multiply_real(B, B, powers[i][ei], alpha, alpha, mins[i][ei]);
				else
					multiply_imag(B, B, powers[i][ei]);
			}

			// select random index from S
			if (S.size() > 0) {
				rg >> ei;
				ei %= S.size();
				if (e[S[ei]] == 0) {
					e[S[ei]] = 1;
					if (real_order) {
						beta.assign_one();
						multiply_real(B, B, fact_base[S[ei]], alpha, alpha, beta);
					}
					else
						multiply_imag(B, B, fact_base[S[ei]]);
				}
			}

			if (real_order) {
				oa.assign(abs(B.get_c()));
				ee = e;
			}
		}

		// try to factor B over the factor base.  If successful, add relation to
		// the relation matrix
		if (factor_over_FB(B, v)) {
			subtract(e, v, e);

#ifdef LIDIA_DEBUG
			qi_class TT, PP;
			TT.assign_one();
			for (i = 0; i < fact_base.size(); ++i) {
				if (e[i] != 0) {
					power(PP, fact_base[i], e[i]);
					multiply(TT, TT, PP);
				}
			}
			if (TT.is_one())
				std::cout << "REL OK!!!" << std::endl;
			else
				std::cout << "BAD REL:  " << e << ", TT = " << TT << std::endl;
#endif

			if (curr_A == Acols) {
				++Acols;
				AMAT.set_no_of_columns(Acols);
			}
			for (i = 0; i < size_FB; ++i)
				if (e[i] != 0)
					AMAT.sto(i, curr_A, static_cast<long>(e[i]));

			if (real_order)
				minima_qn[curr_A] = alpha;

			++curr_A;

			if (A1cols) {
				if (curr_A1 == A1cols) {
					++A1cols;
					A1.set_no_of_columns(A1cols);
				}
				for (i = 0; i < size_FB; ++i)
					if (e[i] != 0)
						A1.sto(i, curr_A1, bigint(e[i]));

				++curr_A1;
			}

			// update index set
			for (i = S.size()-1; i >= 0; --i) {
				if (e[S[i]] != 0) {
					for (j = i; j < S.size()-1; ++j)
						S[j] = S[j+1];
					S.set_size(S.size()-1);
					break;
				}
			}

			++newrels;

			percent = static_cast<unsigned int>((static_cast<double>(numrels + newrels) /
							     static_cast<double>(numrels + n)) * 100);

			if (percent >= per_level) {
///////////// INFO /////////////////////
				if (info > 1) {
					std::cout << numrels + newrels;
					std::cout << " (" << percent << "%) relations, ";
					std::cout << S.size() << " primes left, ";
					std::cout << numvec << " vectors." << std::endl;
				}
///////////// INFO /////////////////////

				if (per_level > 90 && decimal_length(Delta) > 50)
					per_level += 5;
				else
					per_level += 10;
			}

			if (real_order)
				count = 0;
		}
		else {
			if (real_order) {
				++count;
				if (count == 10)
					count = 0;
				else
					B.rho(oa, alpha);
			}
		}
	}


	percent = static_cast<unsigned int>((static_cast<double>(numrels + newrels) /
					     static_cast<double>(numrels + n)) * 100);

	numrels += newrels;


///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "done finding relations.\n";
		std::cout << numrels << " (" << percent << "%) relations, ";
		std::cout << S.size() << " primes left, ";
		std::cout << numvec << " vectors." << std::endl;

		std::cout << "\nRelation generation time:  ";
		MyTime(t.user_time());
		std::cout << std::endl;
	}

	// FOR TABLE GENERATION
	if (special) {
		trel.stop_timer();

		P6 = numrels;
		P7 = numvec;
		P8 = 0;
	}
///////////// INFO /////////////////////
}




//
// quadratic_order::cg_randexp
//
// Task:
//      computes the class group using Dullmann's algorithm
//

void
quadratic_order::cg_randexp()
{
	debug_handler("quadratic_order", "cg_randexp");

	register lidia_size_t i, j, k, n, c, ii;
	lidia_size_t curr_A, curr_A1;
	int size_FB, k0;
	qi_class **powers;
	quadratic_number_standard **mins;
	base_vector< lidia_size_t > power_index;
	sort_vector< lidia_size_t > S;
	matrix< bigint > A1;
	bool full_lattice;
	xbigfloat LL;

///////////// INFO /////////////////////
	timer t;

	// FOR TABLE GENERATION
	if (special) {
		ttot.start_timer();
		tinit.start_timer();
	}

	if (info > 1) {
		t.start_timer();
		std::cout << "\n\nquadratic order subexponential algorithm (random exponents)\n";
		std::cout << "===========================================================\n\n";
		std::cout << "discriminant of order: " << Delta << " (";
		std::cout << decimal_length(Delta) << ")\n";
		std::cout.flush();
	}
///////////// INFO /////////////////////

	lidia_signal sigterm(LIDIA_SIGTERM, qo_stop_it);
	lidia_signal sigint(LIDIA_SIGINT, qo_stop_it);
	lidia_signal sighup(LIDIA_SIGHUP, qo_stop_it);
	lidia_signal sigquit(LIDIA_SIGQUIT, qo_stop_it);
	lidia_signal sigbus(LIDIA_SIGBUS, qo_stop_it);
	lidia_signal sigseg(LIDIA_SIGSEGV, qo_stop_it);

	// initializations
	A1.set_representation(matrix_flags::sparse_representation);

	if (is_real()) {
		if (decimal_length(Delta) > 51)
			c = 20;
		else
			c = 10;
		fu_mult = 1;
	}
	else
		c = 5;


	// read factor base size from list
	qo_read_parameters(decimal_length(Delta), size_FB);

	// get factor base
	qo_get_fb(size_FB);

	// initialize table of powers
	k0 = 0;
	qo_init_powers(k0, &powers, &mins, power_index);

	// initialize relation matrix and index vector
	S.set_capacity(size_FB);
	S.set_size(size_FB);
	for (i = 0; i < size_FB; ++i)
		S[i] = i;

	// determine number of relations to find
	n = size_FB + c;
	AMAT.resize(size_FB, n);
	curr_A = curr_A1 = 0;


	// generate free relations (divisors of Delta)
	generate_free_relations(curr_A, S);


///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\ninitializtion time:  ";
		MyTime(t.user_time());
		std::cout << "\n\n";
		std::cout << "\nUsed parameters:\n";
		std::cout << "  factor base size " << size_FB << "\n";
		std::cout << "  relations required " << n << "\n";
		std::cout << "  max prime in factor base ";
		std::cout << fact_base[size_FB-1].get_a() << "\n";
		std::cout << "\nstarting relation generation:\n" << std::endl;
		std::cout.flush();
	}

	// FOR TABLE GENERATION - DATA COLLECTION
	if (special) {
		tinit.stop_timer();
		T1 = tinit.user_time();

		P1 = size_FB;
		P2 = 0;
		P3 = 0;
		fact_base[size_FB-1].get_a().intify(P4);
		P5 = 0;
		hmult.assign_zero();
	}
///////////// INFO /////////////////////


	//
	// generate relations until we have a complete generating system of the
	// relation lattice
	//

	full_lattice = false;
	while (!full_lattice) {
		// generate n new relations, ensuring that each row of the relation matrix
		// corresponding to the indices in S contains a non-zero entry
		qo_relations_randexp(size_FB, n, powers, mins, power_index, S, A1,
				     curr_A, curr_A1);

		// compute the determinant of the relations and the regulator
		qo_hnf_and_det(A1, S, LL);

////// REMOVE AFTER FIX!!!
		if (is_real() && !h.is_zero() && R.is_zero())
			return;
//////////////////////////

		// test whether the relation lattice is full
		if (is_imaginary())
			full_lattice = (!h.is_zero() && verify_h(LL, Delta, h));
		else
			full_lattice = (!h.is_zero() && !R.is_zero() &&
					verify_hR(LL, Delta, h, R_xbig));

		// if not, compute the row dependencies and generate more relations
		if (!full_lattice) {
			A1.assign(RHNF);
			if (A1.get_no_of_rows() > A1.get_no_of_columns())
				n = A1.get_no_of_rows() - A1.get_no_of_columns() + c;
			else
				n = c;

			AMAT.set_no_of_columns(AMAT.get_no_of_columns() + n);
			A1.set_no_of_columns(A1.get_no_of_columns() + n);

			curr_A1 = RHNF.get_no_of_columns();

///////////// INFO /////////////////////
			if (info > 1) {
				if (h.is_zero()) {
					std::cout << "\nmatrix does not have full rank:  need more relations\n";
					std::cout << "S = " << S << "\n" << std::endl;
				}
				else {
					std::cout << "\ndeterminant = " << h << std::endl;
					if (is_real())
						std::cout << "hR = " << bigfloat(h)*bigfloat(R) << std::endl;
					std::cout << "determinant too large:  need more relations\n" << std::endl;
				}
			}
///////////// INFO /////////////////////

		}

		// FOR TABLE GENERATION
		if (special && hmult.is_zero())
			hmult.assign(h);
	}

	//
	// compute structure of CL
	//

	if (info > 1)
		t.start_timer();

	if (h == 1)
		CL[0] = 1;
	else {
		// delete all rows and columns with diagonal 1
		contributors.reset();
		i = 0;
		for (j = 0; j < size_FB; ++j)
			if (!(RHNF.member(j, j).is_one()))
				contributors[i++] = j;

		// reduce matrix
		matrix< bigint > Amat, junk;
		k = i;
		Amat.resize(k, k);
		for (i = 0; i < k; ++i)
			for (j = i; j < k; ++j)
				Amat.sto(i, j, RHNF.member(contributors[i], contributors[j]));

		bigint Bjj;
		Amat.snf_havas(U, junk);
		ii = 0;
		for (j = 0; j < Amat.get_no_of_columns(); ++j) {
			Bjj.assign(Amat.member(j, j));
			if (!Bjj.is_one())
				CL[ii++] = Bjj;
		}
	}


///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\nSNF time:  ";
		MyTime(t.user_time());
		std::cout << std::endl;
	}

	// FOR TABLE GENERATION
	if (special) {
		ttot.stop_timer();

		hmult /= h;

		T1 = tinit.user_time();
		T2 = trel.user_time();
		T3 = thnf.user_time();
		T4 = thnf2.user_time();
		T5 = treg.user_time();
		T6 = ttot.user_time();

		std::cout << hmult << std::endl;

		std::cout << P1 << "  " << P2 << "  " << P3 << "  " << P4 << "  " << P5 << "  ";
		std::cout << P6 << "  " << P7 << "  " << P8 << std::endl;

		std::cout << T1 << "  " << T2 << "  " << T3 << "  " << T4 << "  " << T5 << "  ";
		std::cout << T6 << std::endl;
	}
///////////// INFO /////////////////////


	// delete table of powers
	bool real_order = is_real();
	for (i = 0; i < power_index.size(); ++i) {
		delete [] powers[i];
		if (real_order)
			delete [] mins[i];
	}
	delete [] powers;
	if (real_order)
		delete [] mins;
}





/////////////////////////////////
//                             //
// SUBEXPONENTIAL ALGORITHM #2 //
//   - MPQS based sieving      //
//                             //
/////////////////////////////////

//
// qo_read_parameters
//
// Task:
//      read the parameters for the sieving versions of the subexponential
//      class group algorithm from a table.
//

void
quadratic_order::qo_read_parameters(int decimal_digits, double &T, int &M,
                                    int &size_FB, int &P_ONCE, int &P_TOTAL,
                                    int &smallstart, int &POLY)
{
	debug_handler("quadratic_order", "qo_read_parameters(sieving)");

	int  i = -1, r8;

	if (decimal_digits < 15)
		i = 0;
	else if (decimal_digits > 80)
		i = 65;
	else
		i = decimal_digits - 15;

	T = qo_params[i][1]; // tolerance value
	M = static_cast<int>(qo_params[i][2]); // sieve radius
	size_FB = static_cast<int>(qo_params[i][3]); // size of factor base

	P_ONCE = static_cast<int>(qo_params[i][4]); // number of prime factors of
	// coefficient A

	P_TOTAL = static_cast<int>(qo_params[i][5]); // total number of primes
	// the P_ONCEs are taken from

	smallstart = static_cast<int>(qo_params[i][6]); // first prime for sieving


	// smaller factor bases for real orders
	if (is_real() && decimal_digits <= 50)
		size_FB -= 500;


	// use large factor bases for really small class numbers
	bigfloat L1 = estimate_L1(1000);

	// use small factor bases for really big class numbers
	r8 = remainder(Delta, 8);
	if (r8 < 0)
		r8 += 8;
	if (r8 != 1)
		L1.multiply_by_2();
	if (L1 > 9.0) {
		size_FB >>= 1;
		if (is_real() && (decimal_digits > 65))
			size_FB += 750;
		else if (is_real() && (decimal_digits > 50))
			size_FB += 250;
	}
	else if (L1 < 0.4) {
		if (is_real() && decimal_digits <= 50)
			size_FB += 500;

		if ((size_FB > 300) && (size_FB < 500))
			size_FB += size_FB + (size_FB >> 1);
		else
			size_FB <<= 1;
	}
	else if (is_real() && (size_FB > 3500))
		size_FB -= 500;


	// max FB of 7000
	if (is_real()) {
		if (size_FB > 4000)
			size_FB = 4000;
	}
	else {
		if (size_FB > 7000)
			size_FB = 7000;
	}


	if (decimal_digits > 80) {
		i = decimal_digits - 80;
		++P_ONCE;
		P_TOTAL += i;
	}

	// number of coefficients B that are used for one A
	POLY = static_cast<int>(std::pow (2.0, static_cast<double>(P_ONCE - 1)));
}




//
// quadratic_order::qo_relations_mpqs(bool self_init)
//
// Task:
//      compute n new relations with a variation of the MPQS factoring
//      algorithm, ensuring that each row of the relation matrix corresponding
//      to the indices in S contains a non-zero entry.  If self_init is true,
//      self-initialization will be used.
//

void
quadratic_order::qo_relations_mpqs(int size_FB,
                                   int n,
                                   sort_vector< lidia_size_t > & S,
                                   matrix< bigint > & A1,
                                   lidia_size_t & curr_A,
                                   lidia_size_t & curr_A1,
                                   bool self_init)
{
	debug_handler("quadratic_order", "qo_relations_mpqs");

	register lidia_size_t i, j;
	lidia_size_t Acols, A1cols, newrels, new_lp_hits, temp;
	quadratic_number_standard alpha;
	math_vector< lidia_size_t > v(size_FB, size_FB);
	bool real_order = is_real();
	unsigned int percent;
	static lidia_size_t numrels = 0, lp_hits = 0, numpoly = 0;
	static unsigned int per_level = 10;
	static int times_called = 0;

	Acols = AMAT.get_no_of_columns();
	A1cols = A1.get_no_of_columns();

	newrels = new_lp_hits = 0;
	if (A1cols == 0) {
		per_level = 10;
		lp_hits = numpoly = 0;
		newrels = curr_A;
		numrels = 0;
		times_called = 0;
	}

	++times_called;
	if (times_called > 5) {
		sieve.no_self_init();
		n += 20;
	}

///////////// INFO /////////////////////
	timer t;
	if (info > 1)
		t.start_timer();

	// FOR TABLE GENERATION
	if (special) {
		if (A1cols == 0)
			trel.start_timer();
		else
			trel.cont_timer();
	}
///////////// INFO /////////////////////


	// continue to generate relations until S is empty and i>=n
	while (newrels < n) {
		// generate next polynomial and sieve over it
		++numpoly;
		sieve.next_polynomial();
		sieve.sieve();

		// add new relations to the relation matrix
		while ((newrels < n) && sieve.get_relation(v, alpha)) {
#ifdef LIDIA_DEBUG
			if (Delta.is_lt_zero()) {
				qi_class TT, PP;
				TT.assign_one();
				for (i = 0; i < fact_base.size(); ++i) {
					if (v[i] != 0) {
						power(PP, fact_base[i], v[i]);
						multiply(TT, TT, PP);
					}
				}
				if (!TT.is_one())
					std::cout << "BAD REL:  " << v << ", TT = " << TT << std::endl;
			}
			else {
				quadratic_ideal TT, PP, TT2;
				TT.assign_one(*this);
				for (i = 0; i < fact_base.size(); ++i) {
					if (v[i] != 0) {
						power(PP, quadratic_ideal(fact_base[i]), v[i]);
						multiply(TT, TT, PP);
					}
				}
//        TT2.assign_principal(alpha);
				TT2.assign_one(*this);
				multiply(TT2, TT2, alpha);

				if ((TT.get_a() != TT2.get_a()) || (TT.get_b() != TT2.get_b())) {
					std::cout << "BAD REL:  TT = " << TT;
					std::cout << ", TT2 = " << TT2;
					std::cout << ", alpha = " << alpha << std::endl;
				}
			}
#endif

			if (curr_A == Acols) {
				++Acols;
				AMAT.set_no_of_columns(Acols);
			}
			for (i = 0; i < size_FB; ++i)
				if (v[i] != 0)
					AMAT.sto(i, curr_A, static_cast<long>(v[i]));

			if (real_order)
				minima_qn[curr_A] = alpha;

			++curr_A;

			if (A1cols) {
				if (curr_A1 == A1cols) {
					++A1cols;
					A1.set_no_of_columns(A1cols);
				}
				for (i = 0; i < size_FB; ++i)
					if (v[i] != 0)
						A1.sto(i, curr_A1, bigint(v[i]));
				++curr_A1;
			}

			// update index set
			for (i = S.size()-1; i >= 0; --i) {
				if (v[S[i]] != 0) {
					for (j = i; j < S.size()-1; ++j)
						S[j] = S[j+1];
					S.set_size(S.size()-1);
				}
			}

			++newrels;

			percent = static_cast<unsigned int>((static_cast<double>(numrels + newrels + lp_hits +
										 new_lp_hits) /
							     static_cast<double>(numrels + n)) * 100);

			if (percent >= per_level) {
///////////// INFO /////////////////////
				if (info > 1) {
					std::cout << numrels + newrels + new_lp_hits;
					std::cout << " (" << percent << "%) relations ";
					std::cout << "(" << numrels + newrels - lp_hits << " fulls, ";
					std::cout << lp_hits + new_lp_hits << " partials)  ";
					std::cout << S.size() << " primes left.  ";
					std::cout << numpoly << " polys." << std::endl;
				}
///////////// INFO /////////////////////

				if (per_level > 90 && decimal_length(Delta) > 50)
					per_level += 5;
				else
					per_level += 10;

				new_lp_hits = sieve.count_lp_relations();
			}

			// compute large prime relations, if finished
			if (((newrels + new_lp_hits) > n) && (new_lp_hits)) {
				if (info > 1) {
					std::cout << "\ncombining large prime relations..." << std::endl;
				}

				temp = sieve.get_lp_relations();
				lp_hits += temp;

				if (info > 1) {
					std::cout << temp << " large prime relations found (";
					std::cout << new_lp_hits - temp << " zero relations)" << std::endl;
				}

				new_lp_hits = 0;
			}
		}
	}

	// add relations for remaining indices in S

	if ((times_called > 5) && (S.size() == 0)) {
		random_generator rg;
		for (j = 0; j < 5; ++j) {
			rg >> S[j];
			S[j] %= size_FB;
			if (S[j] < 0)
				S[j] += size_FB;
		}
	}

	if ((info > 1) && (S.size() > 0)) {
		std::cout << "\ncomputing relations for unrepresented primes (";
		std::cout << S.size() << " primes)...\n" << std::endl;
	}

	for (j = 0; j < S.size(); ++j) {
		do {
			// generate sieve polynomial with desired factor base element represented
			++numpoly;
			sieve.next_polynomial(S[j]);

			// search for one relation with index S[j] non-zero
		} while (!sieve.sieve_one(S[j], v, alpha));

#ifdef LIDIA_DEBUG
		if (Delta.is_lt_zero()) {
			qi_class TT, PP;
			TT.assign_one();
			for (i = 0; i < fact_base.size(); ++i) {
				if (v[i] != 0) {
					power(PP, fact_base[i], v[i]);
					multiply(TT, TT, PP);
				}
			}
			if (!TT.is_one())
				std::cout << "BAD REL:  " << v << ", TT = " << TT << std::endl;
		}
		else {
			quadratic_ideal TT, PP, TT2;
			TT.assign_one(*this);
			for (i = 0; i < fact_base.size(); ++i) {
				if (v[i] != 0) {
					power(PP, quadratic_ideal(fact_base[i]), v[i]);
					multiply(TT, TT, PP);
				}
			}
//      TT2.assign_principal(alpha);
			TT2.assign_one(*this);
			multiply(TT2, TT2, alpha);

			if ((TT.get_a() != TT2.get_a()) || (TT.get_b() != TT2.get_b())) {
				std::cout << "BAD REL:  TT = " << TT;
				std::cout << ", TT2 = " << TT2;
				std::cout << ", alpha = " << alpha << std::endl;
			}
		}
#endif

		if (curr_A == Acols) {
			++Acols;
			AMAT.set_no_of_columns(Acols);
		}
		for (i = 0; i < size_FB; ++i)
			if (v[i] != 0)
				AMAT.sto(i, curr_A, static_cast<long>(v[i]));

		if (real_order)
			minima_qn[curr_A] = alpha;

		++curr_A;

		if (A1cols) {
			if (curr_A1 == A1cols) {
				++A1cols;
				A1.set_no_of_columns(A1cols);
			}
			for (i = 0; i < size_FB; ++i)
				if (v[i] != 0)
					A1.sto(i, curr_A1, bigint(v[i]));

			++curr_A1;
		}

		++newrels;
	}

	S.set_size(0);


	percent = static_cast<unsigned int>((static_cast<double>(numrels + newrels + lp_hits) /
					     static_cast<double>(numrels + n)) * 100);

	numrels += newrels;

///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "done finding relations.\n";
		std::cout << numrels << " (" << percent << "%) relations ";
		std::cout << "(" << numrels-lp_hits << " fulls, " << lp_hits << " partials)  ";
		std::cout << S.size() << " primes left.  ";
		std::cout << numpoly << " polys." << std::endl;

		std::cout << "\nRelation generation time:  ";
		MyTime(t.user_time());
		std::cout << std::endl;
	}

	// FOR TABLE GENERATION
	if (special) {
		trel.stop_timer();

		P6 = numrels-lp_hits;
		P7 = numpoly;
		P8 = lp_hits;
	}
///////////// INFO /////////////////////
}



//
// quadratic_order::cg_mpqs(bool self_init, bool large_primes)
//
// Task:
//      computes the class group using a relation generation strategy based on
//      the MPQS factoring algorithm.  If self_init is true, then
//      self-initialization is used.  If large_primes is true, the large prime
//      variant is also used.
//

void
quadratic_order::cg_mpqs(bool self_init, bool large_primes)
{
	debug_handler("quadratic_order", "cg_mpqs");

	register lidia_size_t i, j, k, n, c, ii;
	lidia_size_t curr_A, curr_A1;
	int size_FB;
	sort_vector< lidia_size_t > S;
	matrix< bigint > A1;
	bool full_lattice;
	xbigfloat LL;

///////////// INFO /////////////////////
	timer t;

	// FOR TABLE GENERATION
	if (special) {
		ttot.start_timer();
		tinit.start_timer();
	}

	if (info > 1) {
		t.start_timer();
		if (self_init) {
			std::cout << "\n\nquadratic order subexponential algorithm (SIQS)\n";
			std::cout << "===============================================\n\n";
		}
		else {
			std::cout << "\n\nquadratic order subexponential algorithm (MPQS)\n";
			std::cout << "===============================================\n\n";
		}
		std::cout << "discriminant of order: " << Delta << " (";
		std::cout << decimal_length(Delta) << ")\n";
		std::cout.flush();
	}
///////////// INFO /////////////////////

	lidia_signal sigterm(LIDIA_SIGTERM, qo_stop_it);
	lidia_signal sigint(LIDIA_SIGINT, qo_stop_it);
	lidia_signal sighup(LIDIA_SIGHUP, qo_stop_it);
	lidia_signal sigquit(LIDIA_SIGQUIT, qo_stop_it);
	lidia_signal sigbus(LIDIA_SIGBUS, qo_stop_it);
	lidia_signal sigseg(LIDIA_SIGSEGV, qo_stop_it);

	// initializations
	A1.set_representation(matrix_flags::sparse_representation);

	if (is_real()) {
		if (decimal_length(Delta) >= 20)
			c = 50;
		else
			c = 40;
		fu_mult = 1;
	}
	else
		c = 20;

	// read paramters from list and initialize sieve
	int M, P_ONCE, P_TOTAL;
	{
		double Tval;
		int smallstart, POLY;
		long lpbound;

		qo_read_parameters(decimal_length(Delta), Tval, M, size_FB, P_ONCE,
				   P_TOTAL, smallstart, POLY);

		// get factor base
		qo_get_fb(size_FB);

		// initialize the sieve
		if (large_primes)
			lpbound = 10000000;
		else
			lpbound = 0;
		sieve.init(Delta, fact_base, Tval, M, P_ONCE, P_TOTAL, smallstart, POLY,
			   lpbound);

		if (self_init)
			sieve.init_self_initialization();
	}

	// initialize relation matrix and index vector
	S.set_capacity(size_FB);
	S.set_size(size_FB);
	for (i = 0; i < size_FB; ++i)
		S[i] = i;

	// determine number of relations to find
	n = size_FB + c;
	AMAT.resize(size_FB, n);
	curr_A = curr_A1 = 0;



///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\ninitializtion time:  ";
		MyTime(t.user_time());
		std::cout << "\n\n";
		std::cout << "\nUsed parameters:\n";
		std::cout << "  factor base size " << size_FB << "\n";
		std::cout << "  sieve interval [-" << M << ", " << M << "]\n";
		std::cout << "  relations required " << n << "\n";
		std::cout << "  max prime in factor base ";
		std::cout << fact_base[size_FB-1].get_a() << "\n";
		std::cout << "\nstarting relation generation (sieving):\n" << std::endl;
		std::cout.flush();
	}

	// FOR TABLE GENERATION - DATA COLLECTION
	if (special) {
		tinit.stop_timer();
		T1 = tinit.user_time();

		P1 = size_FB;
		P2 = P_TOTAL;
		P3 = P_ONCE;
		fact_base[size_FB-1].get_a().intify(P4);
		P5 = M;
		hmult.assign_zero();
	}
///////////// INFO /////////////////////


	// generate free relations (divisors of Delta)
	generate_free_relations(curr_A, S);


	//
	// generate relations until we have a complete generating system of the
	// relation lattice
	//

	full_lattice = false;
	while (!full_lattice) {
		// generate n new relations, ensuring that each row of the relation matrix
		// corresponding to the indices in S contains a non-zero entry

		qo_relations_mpqs(size_FB, n, S, A1, curr_A, curr_A1, self_init);

		// compute the determinant of the relations and the regulator
		qo_hnf_and_det(A1, S, LL);

		if (special && noCL)
			return;

////// REMOVE AFTER FIX!!!
		if (is_real() && !h.is_zero() && R.is_zero())
			return;
//////////////////////////


		// test whether the relation lattice is full
		if (is_imaginary())
			full_lattice = (!h.is_zero() && verify_h(LL, Delta, h));
		else
			full_lattice = (!h.is_zero() && !R.is_zero() &&
					verify_hR(LL, Delta, h, R_xbig));

		// if not, compute the row dependencies and generate more relations
		if (!full_lattice) {
			A1.assign(RHNF);
			if (A1.get_no_of_rows() > A1.get_no_of_columns())
				n = A1.get_no_of_rows() - A1.get_no_of_columns() + c;
			else
				n = c;

			AMAT.set_no_of_columns(AMAT.get_no_of_columns() + n);
			A1.set_no_of_columns(A1.get_no_of_columns() + n);

			curr_A1 = RHNF.get_no_of_columns();

///////////// INFO /////////////////////
			if (info > 1) {
				if (h.is_zero()) {
					std::cout << "\nmatrix does not have full rank:  need more relations\n";
					std::cout << "S = " << S << "\n" << std::endl;
				}
				else {
					std::cout << "\ndeterminant = " << h << std::endl;
					if (is_real())
						std::cout << "hR = " << bigfloat(h)*bigfloat(R) << std::endl;
					std::cout << "determinant too large:  need more relations\n" << std::endl;
				}
			}
///////////// INFO /////////////////////

		}

		// FOR TABLE GENERATION
		if (special && hmult.is_zero())
			hmult.assign(h);
	}

	//
	// compute structure of CL
	//

	if (info > 1)
		t.start_timer();

	if (h == 1)
		CL[0] = 1;
	else {
		// delete all rows and columns with diagonal 1
		contributors.reset();
		i = 0;
		for (j = 0; j < size_FB; ++j)
			if (!(RHNF.member(j, j).is_one()))
				contributors[i++] = j;

		// reduce matrix
		matrix< bigint > Amat, junk;
		k = i;
		Amat.resize(k, k);
		for (i = 0; i < k; ++i)
			for (j = i; j < k; ++j)
				Amat.sto(i, j, RHNF.member(contributors[i], contributors[j]));

		bigint Bjj;
		Amat.snf_havas(U, junk);
		ii = 0;
		for (j = 0; j < Amat.get_no_of_columns(); ++j) {
			Bjj.assign(Amat.member(j, j));
			if (!Bjj.is_one())
				CL[ii++] = Bjj;
		}
	}

///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\nSNF time:  ";
		MyTime(t.user_time());
		std::cout << std::endl;
	}

	// FOR TABLE GENERATION
	if (special) {
		ttot.stop_timer();

		hmult /= h;

		T1 = tinit.user_time();
		T2 = trel.user_time();
		T3 = thnf.user_time();
		T4 = thnf2.user_time();
		T5 = treg.user_time();
		T6 = ttot.user_time();

		std::cout << hmult << std::endl;


		sieve.dense_bound(P2);
		std::cout << P1 << "  " << P2 << "  " << P3 << "  " << P4 << "  " << P5 << "  ";
		std::cout << P6 << "  " << P7 << "  " << P8 << std::endl;

		std::cout << T1 << "  " << T2 << "  " << T3 << "  " << T4 << "  " << T5 << "  ";
		std::cout << T6 << std::endl;
	}
///////////// INFO /////////////////////

}




//
// quadratic_order::represent_over_generators(math_vector)
//
// Task:
//      computes the representation of the ideal (given by its representation
//      over the factor base) over the system of generators of the class
//      group.
//

math_vector< bigint >
quadratic_order::represent_over_generators(math_vector< bigint > & vec)
{
	debug_handler("quadratic_order", "represent_over_generators(math_vector)");

	lidia_size_t i, numFB, numCL, numCont;
	math_vector< bigint > vec2, nvec, res, temp;

	// rows of U are the exponent vectors for the generators which correspond
	// to each prime ideal in the factor base
	numCL = CL.size();
	numFB = fact_base.size();
	numCont = contributors.size();

	temp.set_capacity(numFB);
	res.set_capacity(numCL);
	vec2.set_capacity(numFB);
	nvec.set_capacity(numCont);
	temp.set_size(numFB);
	res.set_size(numCL);
	vec2.set_size(numFB);
	nvec.set_size(numCont);

	// convert vector to exponents only using factor base elements which
	// contributed to the class group
	vec2.assign(vec);
	for (i = numFB-1; i >= 0; --i) {
		if (RHNF.member(i, i).is_one()) {
			temp.assign((math_vector< bigint > ) RHNF.get_column_vector(i));
			multiply(temp, temp, vec[i]);
			subtract(vec2, vec2, temp);
		}
	}
	for (i = 0; i < numCont; ++i)
		nvec[i].assign(vec2[contributors[i]]);

	// initialize result vector
	for (i = 0; i < numCL; ++i)
		res[i].assign_zero();

	// compute representation over generators
	temp.set_capacity(numCL);
	temp.set_size(numCL);
	for (i = 0; i < numCont; ++i) {
		temp.assign((math_vector< bigint > ) U.get_column_vector(i));
		multiply(temp, temp, nvec[i]);
		add(res, res, temp);
	}

	// reduce modulo orders of generators
	for (i = 0; i < numCL; ++i) {
		remainder(res[i], res[i], CL[i]);
		if (res[i].is_lt_zero())
			add(res[i], res[i], CL[i]);
	}

	return res;
}



/////////////////////////////////
//                             //
// FACTORING AN IDEAL OVER FB  //
//   - random exponents        //
//                             //
/////////////////////////////////

//
// quadratic_order::represent_over_FB
//
// Task:
//      computes the representation of the ideal over the factor base.
//

math_vector< bigint >
quadratic_order::represent_over_FB(const qi_class & AA,
                                   quadratic_number_standard & q)
{
	debug_handler("quadratic_order", "represent_over_FB(qi_class, "
		      "quadratic_number_standard");

	register lidia_size_t i;
	lidia_size_t Bound = 30, ei;
	qi_class B, C;
	random_generator rg;
	bool real_order = is_real();
	int size_FB;
	math_vector< bigint > vec;

	quadratic_order *QO = qi_class::get_current_order_ptr();
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	vec.set_mode(EXPAND);

	qo_read_parameters(decimal_length(Delta), size_FB);

	if (fact_base.size() == 0)
		qo_get_fb(size_FB);
	else
		size_FB = fact_base.size();

	math_vector< lidia_size_t > e(size_FB, size_FB), v(size_FB, size_FB);

	// select random exponent vectors until successful
	B.assign(AA);
	while (!factor_over_FB(B, v)) {
		for (i = 0; i < size_FB; ++i)
			e[i] = 0;

		// first dense part
		B.assign(AA);
		for (i = 0; i < size_FB; ++i) {
			if ((Delta % fact_base[i].get_a()) != 0) {
				rg >> ei;
				ei = (ei % Bound) + 1;

				e[i] = ei;
				if (real_order) {
					power_real(C, fact_base[i], ei);
					multiply_real(B, B, C);
				}
				else {
					power_imag(C, fact_base[i], ei);
					multiply_imag(B, B, C);
				}
			}
		}
	}

	subtract(v, v, e);
	for (i = 0; i < size_FB; ++i)
		vec[i].assign(v[i]);

#ifdef LIDIA_DEBUG
	qi_class TT, PP;
	TT.assign_one();
	for (i = 0; i < fact_base.size(); ++i) {
		if (vec[i] != 0) {
			power(PP, fact_base[i], vec[i]);
			multiply(TT, TT, PP);
		}
	}
	if (TT == AA)
		std::cout << "REL OK!!!" << std::endl;
	else
		std::cout << "BAD REL:  " << vec << ", TT = " << TT << ", A = " << AA << std::endl;
#endif

	qi_class::set_current_order(*QO);

	return vec;
}




//
// quadratic_order::represent_over_generators(qi_class)
//
// Task:
//      computes the representation of the ideal over the system of generators
//      of the class group.
//

math_vector< bigint >
quadratic_order::represent_over_generators(qi_class & AA)
{
	debug_handler("quadratic_order", "represent_over_generators(qi_class)");

	quadratic_number_standard q;
	math_vector< bigint > vec, nvec;
	quadratic_order *QO = qi_class::get_current_order_ptr();

	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	generators();

	vec.set_capacity(fact_base.size());
	vec.set_size(fact_base.size());
	vec = represent_over_FB(AA, q);

	nvec.set_mode(EXPAND);
	nvec = represent_over_generators(vec);

	qi_class::set_current_order(*QO);

	return nvec;
}





/////////////////////////////////
//                             //
// FACTORING AN IDEAL OVER FB  //
//   - MPQS based sieving      //
//                             //
/////////////////////////////////

//
// quadratic_order::represent_over_FB_sieve
//
// Task:
//      computes the representation of the ideal over the factor base.
//

math_vector< bigint >
quadratic_order::represent_over_FB_sieve(const qi_class & AA,
                                         quadratic_number_standard & q)
{
	debug_handler("quadratic_order", "represent_over_FB_sieve()");

	register lidia_size_t i, count;
	qi_class B;
	int size_FB;

	lidia_signal sigterm(LIDIA_SIGTERM, qo_stop_it);
	lidia_signal sigint(LIDIA_SIGINT, qo_stop_it);
	lidia_signal sighup(LIDIA_SIGHUP, qo_stop_it);
	lidia_signal sigquit(LIDIA_SIGQUIT, qo_stop_it);
	lidia_signal sigbus(LIDIA_SIGBUS, qo_stop_it);
	lidia_signal sigseg(LIDIA_SIGSEGV, qo_stop_it);


	qo_read_parameters(decimal_length(Delta), size_FB);

	if (fact_base.size() == 0)
		qo_get_fb(size_FB);

	size_FB = fact_base.size();
	math_vector< bigint > vec(size_FB, size_FB);
	math_vector< lidia_size_t > lvec(size_FB, size_FB);

	// test if AA itself factors over the factor base
	B.assign(AA);
	if (!factor_over_FB(B, lvec)) {
		// initialize the sieve, if necessary
		if (!sieve.is_initialized() || (sieve.discriminant_used() != Delta) ||
		    (size_FB != sieve.FB_size())) {
			double Tval;
			int M, P_ONCE, P_TOTAL, smallstart, POLY;

			for (i = 0; size_FB < qo_params[i][3]; ++i);

			Tval = qo_params[i][1];
			M = static_cast<int>(qo_params[i][2]);
			P_ONCE = static_cast<int>(qo_params[i][4]);
			P_TOTAL = static_cast<int>(qo_params[i][5]);
			smallstart = static_cast<int>(qo_params[i][6]);
			POLY = static_cast<int>(std::pow (2.0, static_cast<double>(P_ONCE - 1)));


			// initialize the sieve
			sieve.init(Delta, fact_base, Tval, M, P_ONCE, P_TOTAL, smallstart, POLY, 0);
		}
		else
			sieve.restart();

		// compute representation
		count = 0;
		do {
			++count;
			if ((count % 2000) == 0)
				std::cout << "\ncount = " << count << ", P = " << AA << std::endl;
			sieve.next_polynomial(AA);
		} while (!sieve.sieve_one(AA, lvec, q));
	}

	for (i = 0; i < size_FB; ++i)
		vec[i].assign(lvec[i]);

#ifdef LIDIA_DEBUG
	qi_class TT, PP;
	TT.assign_one();
	for (i = 0; i < fact_base.size(); ++i) {
		if (vec[i] != 0) {
			power(PP, fact_base[i], vec[i]);
			multiply(TT, TT, PP);
		}
	}
	if (TT == AA)
		std::cout << "REL OK!!!" << std::endl;
	else
		std::cout << "BAD REL:  " << vec << ", TT = " << TT << ", A = " << AA << std::endl;
#endif

	return vec;
}




//
// quadratic_order::represent_over_FB_sieve
//
// Task:
//      computes the representation of the ideal over the factor base.
//

bool
quadratic_order::represent_over_FB_sieve(const qi_class & AA)
{
	debug_handler("quadratic_order", "represent_over_FB(sieve)");

	register lidia_size_t i;
	qi_class B;
	int size_FB;

	lidia_signal sigterm(LIDIA_SIGTERM, qo_stop_it);
	lidia_signal sigint(LIDIA_SIGINT, qo_stop_it);
	lidia_signal sighup(LIDIA_SIGHUP, qo_stop_it);
	lidia_signal sigquit(LIDIA_SIGQUIT, qo_stop_it);
	lidia_signal sigbus(LIDIA_SIGBUS, qo_stop_it);
	lidia_signal sigseg(LIDIA_SIGSEGV, qo_stop_it);

	qo_read_parameters(decimal_length(Delta), size_FB);

	if (fact_base.size() == 0)
		qo_get_fb(size_FB);

	size_FB = fact_base.size();
	math_vector< lidia_size_t > vec(size_FB, size_FB);

	// test if AA itself factors over the factor base
	B.assign(AA);
	if (factor_over_FB(B, vec)) {
		return true;
	}


	// initialize the sieve, if necessary
	if (!sieve.is_initialized() || (sieve.discriminant_used() != Delta) ||
	    (size_FB != sieve.FB_size())) {
		double Tval;
		int M, P_ONCE, P_TOTAL, smallstart, POLY;

		for (i = 0; size_FB < qo_params[i][3]; ++i);

		Tval = qo_params[i][1];
		M = static_cast<int>(qo_params[i][2]);
		P_ONCE = static_cast<int>(qo_params[i][4]);
		P_TOTAL = static_cast<int>(qo_params[i][5]);
		smallstart = static_cast<int>(qo_params[i][6]);
		POLY = static_cast<int>(std::pow (2.0, static_cast<double>(P_ONCE - 1)));


		// initialize the sieve
		sieve.init(Delta, fact_base, Tval, M, P_ONCE, P_TOTAL, smallstart, POLY,
			   20000000);
	}
	else
		sieve.restart();

	// compute representation
	do {
		sieve.next_polynomial(AA);
	} while (!sieve.sieve_one(AA));

	return true;
}




//
// quadratic_order::represent_over_generators_sieve(qi_class)
//
// Task:
//      computes the representation of the ideal over the system of generators
//      of the class group, using sieving.
//

math_vector< bigint >
quadratic_order::represent_over_generators_sieve(qi_class & AA)
{
	debug_handler("quadratic_order", "represent_over_generators_sieve(qi_class)");

	quadratic_number_standard q;
	math_vector< bigint > vec, nvec;
	quadratic_order *QO = qi_class::get_current_order_ptr();

	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	generators();

	vec.set_capacity(fact_base.size());
	vec.set_size(fact_base.size());

	vec = represent_over_FB_sieve(AA, q);

	nvec.set_mode(EXPAND);
	nvec = represent_over_generators(vec);

	qi_class::set_current_order(*QO);

	return nvec;
}



//
// quadratic_order::cg_rel_test(bool self_init)
//
// Task:
//      returns the time required to generate n relations
//

long
quadratic_order::cg_rel_test(lidia_size_t n,
                             int size_FB_in,
                             int M_in,
                             int & size_FB_used,
                             int & M_used,
                             int & PTOTAL_used)
{
	debug_handler("quadratic_order", "cg_rel_test");

	lidia_size_t curr_A, curr_A1;
	int size_FB;
	sort_vector< lidia_size_t > S;
	matrix< bigint > A1;
	timer t;

	t.start_timer();

	// initializations
	qi_class::set_current_order(*this);

	AMAT.set_zero_element(0);
	AMAT.reset();
	AMAT.set_representation(matrix_flags::sparse_representation);
	A1.set_representation(matrix_flags::sparse_representation);

	// read paramters from list and initialize sieve
	int M, P_ONCE, P_TOTAL;
	{
		double Tval;
		int smallstart, POLY;
		long lpbound;

		qo_read_parameters(decimal_length(Delta), Tval, M, size_FB, P_ONCE,
				   P_TOTAL, smallstart, POLY);

		if (size_FB_in > 0)
			size_FB = size_FB_in;
		if (M_in > 0)
			M = M_in;
		size_FB_used = size_FB;
		M_used = M;
		PTOTAL_used = P_TOTAL;

		// get factor base
		qo_get_fb(size_FB);

		// initialize the sieve
		lpbound = 0;
		sieve.init(Delta, fact_base, Tval, M, P_ONCE, P_TOTAL, smallstart, POLY,
			   lpbound);

		sieve.init_self_initialization();
	}

	// initialize relation matrix and index vector
	AMAT.reset();
	S.set_capacity(size_FB);
	S.set_size(0);

	// determine number of relations to find
	AMAT.resize(size_FB, n);
	curr_A = curr_A1 = 0;

	qo_relations_mpqs(size_FB, n, S, A1, curr_A, curr_A1, true);

	t.stop_timer();

	AMAT.reset();
	A1.reset();
	minima_qn.reset();

	return t.user_time();
}





//
// quadratic_order::class_group_siqs(int new_size_FB)
//
// Task:
//      computes the class group using a relation generation strategy based on
//      the MPQS factoring algorithm.  If self_init is true, then
//      self-initialization is used.
//

void
quadratic_order::class_group_siqs(int new_size_FB, int new_M, lidia_size_t c, bool in_reg)
{
	debug_handler("quadratic_order", "class_group_siqs(int)");

	register lidia_size_t i, j, k, n, ii;
	lidia_size_t curr_A, curr_A1;
	int size_FB;
	sort_vector< lidia_size_t > S;
	matrix< bigint > A1;
	bool full_lattice;
	xbigfloat LL;
	timer tt;

	bool self_init = true;

	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	if (info)
		tt.start_timer();

	qi_class::set_current_order(*this);

///////////// INFO /////////////////////
	timer t;

	// FOR TABLE GENERATION
	if (special) {
		ttot.start_timer();
		tinit.start_timer();
	}

	if (info > 1) {
		t.start_timer();
		if (self_init) {
			std::cout << "\n\nquadratic order subexponential algorithm (SIQS)\n";
			std::cout << "===============================================\n\n";
		}
		else {
			std::cout << "\n\nquadratic order subexponential algorithm (MPQS)\n";
			std::cout << "===============================================\n\n";
		}
		std::cout << "discriminant of order: " << Delta << " (";
		std::cout << decimal_length(Delta) << ")\n";
		std::cout.flush();
	}
///////////// INFO /////////////////////

	lidia_signal sigterm(LIDIA_SIGTERM, qo_stop_it);
	lidia_signal sigint(LIDIA_SIGINT, qo_stop_it);
	lidia_signal sighup(LIDIA_SIGHUP, qo_stop_it);
	lidia_signal sigquit(LIDIA_SIGQUIT, qo_stop_it);
	lidia_signal sigbus(LIDIA_SIGBUS, qo_stop_it);
	lidia_signal sigseg(LIDIA_SIGSEGV, qo_stop_it);

	// initializations
	AMAT.set_zero_element(0);
	AMAT.reset();
	AMAT.set_representation(matrix_flags::sparse_representation);
	A1.set_representation(matrix_flags::sparse_representation);

	if (c == 0) {
		if (is_real()) {
			if (decimal_length(Delta) > 51)
				c = 50;
			else
				c = 40;
			fu_mult = 1;
		}
		else
			c = 20;
	}

	// read paramters from list and initialize sieve
	int M, P_ONCE, P_TOTAL;
	{
		double Tval;
		int smallstart, POLY;
		long lpbound;

		qo_read_parameters(decimal_length(Delta), Tval, M, size_FB, P_ONCE,
				   P_TOTAL, smallstart, POLY);

		size_FB = new_size_FB;

		if (new_M > 0)
			M = new_M;

		// get factor base
		qo_get_fb(size_FB);

		// initialize the sieve
		lpbound = 0;
		sieve.init(Delta, fact_base, Tval, M, P_ONCE, P_TOTAL, smallstart, POLY,
			   lpbound);

		if (self_init)
			sieve.init_self_initialization();
	}


	// initialize relation matrix and index vector
	AMAT.reset();
	S.set_capacity(size_FB);
	S.set_size(size_FB);
	for (i = 0; i < size_FB; ++i)
		S[i] = i;

	// determine number of relations to find
	n = size_FB + c;
	AMAT.resize(size_FB, n);
	curr_A = curr_A1 = 0;



///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\ninitializtion time:  ";
		MyTime(t.user_time());
		std::cout << "\n\n";
		std::cout << "\nUsed parameters:\n";
		std::cout << "  factor base size " << size_FB << "\n";
		std::cout << "  relations required " << n << "\n";
		std::cout << "  max prime in factor base ";
		std::cout << fact_base[size_FB-1].get_a() << "\n";
		std::cout << "\nstarting relation generation (sieving):\n\n" << std::flush;
		std::cout.flush();
	}

	// FOR TABLE GENERATION - DATA COLLECTION
	if (special) {
		tinit.stop_timer();
		T1 = tinit.user_time();

		P1 = size_FB;
		P2 = P_TOTAL;
		P3 = P_ONCE;
		fact_base[size_FB-1].get_a().intify(P4);
		P5 = M;
		hmult.assign_zero();
	}
///////////// INFO /////////////////////


	// generate free relations (divisors of Delta)
	generate_free_relations(curr_A, S);


	//
	// generate relations until we have a complete generating system of the
	// relation lattice
	//

	full_lattice = false;
	while (!full_lattice) {
		// generate n new relations, ensuring that each row of the relation matrix
		// corresponding to the indices in S contains a non-zero entry

		qo_relations_mpqs(size_FB, n, S, A1, curr_A, curr_A1, self_init);

		// compute the determinant of the relations and the regulator
		qo_hnf_and_det(A1, S, LL, in_reg);

		// test whether the relation lattice is full
		if (is_imaginary())
			full_lattice = (!h.is_zero() && verify_h(LL, Delta, h));
		else
			full_lattice = (!h.is_zero() && !R.is_zero() &&
					verify_hR(LL, Delta, h, R_xbig));

		// if not, compute the row dependencies and generate more relations
		if (!full_lattice) {
			A1.assign(RHNF);
			if (A1.get_no_of_rows() > A1.get_no_of_columns())
				n = A1.get_no_of_rows() - A1.get_no_of_columns() + c;
			else
				n = c;

			AMAT.set_no_of_columns(AMAT.get_no_of_columns() + n);
			A1.set_no_of_columns(A1.get_no_of_columns() + n);

			curr_A1 = RHNF.get_no_of_columns();

///////////// INFO /////////////////////
			if (info > 1) {
				if (h.is_zero()) {
					std::cout << "\nmatrix does not have full rank:  need more relations\n";
					std::cout << "S = " << S << "\n\n" << std::flush;
				}
				else {
					std::cout << "\ndeterminant = " << h << std::endl;
					if (is_real())
						std::cout << "hR = " << bigfloat(h)*bigfloat(R) << std::endl;
					std::cout << "determinant too large:  need more relations\n\n" << std::flush;
				}
			}
///////////// INFO /////////////////////

		}

		// FOR TABLE GENERATION
		if (special && hmult.is_zero())
			hmult.assign(h);
	}

	//
	// compute structure of CL
	//

	if (info > 1)
		t.start_timer();

	if (h == 1)
		CL[0] = 1;
	else {
		// delete all rows and columns with diagonal 1
		contributors.reset();
		i = 0;
		for (j = 0; j < size_FB; ++j)
			if (!(RHNF.member(j, j).is_one()))
				contributors[i++] = j;

		// reduce matrix
		matrix< bigint > Amat, junk;
		k = i;
		Amat.resize(k, k);
		for (i = 0; i < k; ++i)
			for (j = i; j < k; ++j)
				Amat.sto(i, j, RHNF.member(contributors[i], contributors[j]));

		bigint Bjj;
		Amat.snf_havas(U, junk);
		ii = 0;
		for (j = 0; j < Amat.get_no_of_columns(); ++j) {
			Bjj.assign(Amat.member(j, j));
			if (!Bjj.is_one())
				CL[ii++] = Bjj;
		}
	}

///////////// INFO /////////////////////
	if (info > 1) {
		t.stop_timer();
		std::cout << "\nSNF time:  ";
		MyTime(t.user_time());
		std::cout << std::endl;
	}

	// FOR TABLE GENERATION
	if (special) {
		ttot.stop_timer();

		hmult /= h;

		T1 = tinit.user_time();
		T2 = trel.user_time();
		T3 = thnf.user_time();
		T4 = thnf2.user_time();
		T5 = treg.user_time();
		T6 = ttot.user_time();

		std::cout << hmult << std::endl;


		sieve.dense_bound(P2);
		std::cout << P1 << "  " << P2 << "  " << P3 << "  " << P4 << "  " << P5 << "  ";
		std::cout << P6 << "  " << P7 << "  " << P8 << std::endl;

		std::cout << T1 << "  " << T2 << "  " << T3 << "  " << T4 << "  " << T5 << "  ";
		std::cout << T6 << std::endl;
	}
///////////// INFO /////////////////////

	subused = true;


	if (info) {
		tt.stop_timer();
		std::cout << "\nelapsed CPU time (class_group_siqs) = ";
		MyTime(tt.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_regulator()) {
			if (info)
				std::cout << "Regulator is correct." << std::endl;
		}
		else
			std::cout << "Regulator is not correct!" << std::endl;

		if (verify_class_group()) {
			if (info)
				std::cout << "Class group is correct." << std::endl;
		}
		else
			std::cout << "Class group is not correct!" << std::endl;
	}

	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);
}



// computing HNF and R

void
quadratic_order::compute_HNF_and_R(char *inname, char *outname)
{
	debug_handler("quadratic_order", "compute_HNF_and_R");

	int size_FB, P_ONCE, P_TOTAL, M, bin_index, max_bin, curr_A;
	sort_vector< lidia_size_t > S;
	matrix< bigint > A1;
	xbigfloat LL;
	std::ofstream out;
	std::ifstream in;

	S.set_mode(EXPAND);
	AMAT.set_zero_element(0);
	AMAT.reset();
	AMAT.set_representation(matrix_flags::sparse_representation);
	A1.set_representation(matrix_flags::sparse_representation);

	std::cout << "OPEN " << inname << std::endl;
	// read data
	in.open(inname);
	in >> Delta;
	std::cout << "Delta = " << Delta << std::endl;
	assign(Delta);
	in >> size_FB;
	std::cout << "size_FB = " << size_FB << std::endl;
	in >> M;
	std::cout << "M = " << M << std::endl;
	in >> bin_index;
	std::cout << "bin_index = " << bin_index << std::endl;
	in >> max_bin;
	std::cout << "max_bin = " << max_bin << std::endl;
	in >> P_ONCE;
	std::cout << "P_ONCE = " << P_ONCE << std::endl;
	in >> P_TOTAL;
	std::cout << "P_TOTAL = " << P_TOTAL << std::endl;
	in >> curr_A;
	std::cout << "curr_A = " << curr_A << std::endl;
	AMAT.set_print_mode(LIDIA_MODE);
	in >> AMAT;
	std::cout << "got AMAT" << std::endl;
	if (is_real())
		in >> minima_qn;
	in >> S;
	std::cout << "got S = " << S << std::endl;
	in.close();

	qi_class::set_current_order(*this);
	qo_hnf_and_det(A1, S, LL);

	// output results
	out.open(outname);
	out << size_FB << std::endl;
	out << M << std::endl;
	out << bin_index << std::endl;
	out << max_bin << std::endl;
	out << P_ONCE << std::endl;
	out << P_TOTAL << std::endl;
	out << R << std::endl;
	out << R_xbig << std::endl;
	out << h << std::endl;
	out << LL << std::endl;
	out << curr_A << std::endl;
	AMAT.set_print_mode(LIDIA_MODE);
	out << AMAT;
	RHNF.set_print_mode(LIDIA_MODE);
	out << RHNF;
	out << S << std::endl;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
