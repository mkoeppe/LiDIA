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


#include	"LiDIA/qi_class_real.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	quadratic_order QO;
	qi_class_real A, B, C, F, U;
	base_vector< qi_class_real > Gvec;
	base_vector< bigint > S;
	bigint p, x, y, dlog, D;
	bigfloat dist;
	int i;
	bool testbool;

	Gvec.set_mode(EXPAND);

	// test whether qi_class_real's are valid before an order is initialized
	std::cout << "No order initialized:  D = " << A.discriminant() << "\n" << std::flush;
	std::cout << "Current order:\n" << A.get_current_order() << "\n" << std::flush;

	// perform tests on 3 selected discriminants
	for (i = 0; i < 3; ++i) {
		std::cout << "\n==================================\n" << std::flush;

		if (i == 0)
			QO.assign(241);
		else if (i == 1)
			QO.assign(14401);
		else {
			string_to_bigint("12312312312312301", D);
			QO.assign(D);
		}

		qi_class_real::set_current_order(QO);

		// test whether the current order is correctly set
		std::cout << "Discriminant = " << U.discriminant() << "\n" << std::flush;
		std::cout << "Current order:\n" << A.get_current_order() << "\n" << std::flush;
		if (A.discriminant() != U.discriminant()) {
			std::cout << "ERROR:  discriminants not equal!\n" << std::flush;
			return 1;
		}

		// test assign one function
		U.assign_one();
		std::cout << "(1) = " << U << "\n" << std::flush;

		// generate a prime ideal and test conversions
		p = 3;
		while (!generate_prime_ideal(A, p))
			p = next_prime(p);
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "A as qi_class - " << qi_class(A) << "\n" << std::flush;
		std::cout << "A as quadratic_ideal - " << quadratic_ideal(A) << "\n" << std::flush;
		std::cout << "A as quadratic_form - " << quadratic_form(A) << "\n\n" << std::flush;

		// generate another prime ideal and test conversions
		p = next_prime(p);
		while (!generate_prime_ideal(B, p))
			p = next_prime(p);
		std::cout << "B = " << B << "\n" << std::flush;
		std::cout << "B.a = " << B.get_a() << "\n" << std::flush;
		std::cout << "B.b = " << B.get_b() << "\n" << std::flush;
		std::cout << "B.c = " << B.get_c() << "\n" << std::flush;
		std::cout << "B.dist = " << B.get_distance() << "\n" << std::flush;
		std::cout << "inverse = " << -B << "\n" << std::flush;
		std::cout << "B*inverse(B) = " << B*inverse(B) << "\n" << std::flush;
		std::cout << "B as qi_class - " << qi_class(B) << "\n" << std::flush;
		std::cout << "B as quadratic_ideal - " << quadratic_ideal(B) << "\n" << std::flush;
		std::cout << "B as quadratic_form - " << quadratic_form(B) << "\n" << std::flush;
		std::cout << "B (using constructor) - " << qi_class_real(B.get_a(), B.get_b()) << "\n\n" << std::flush;

		// test arithmetic operations
		std::cout << "Testing multiplication...\n" << std::flush;
		std::cout << "A*B = " << A * B << "\n" << std::flush;
		std::cout << "Testing square...\n" << std::flush;
		square(C, A);
		std::cout << "A*A = " << C << "\n" << std::flush;
		std::cout << "Testing division...\n" << std::flush;
		std::cout << "A/B = " << A / B << "\n\n" << std::flush;

		// test rho functions
		if (i < 2) {
			std::cout << "Testing rho() and inverse_rho()...\n" << std::flush;
			std::cout << U << "\n" << std::flush;
			apply_rho(C, U);
			std::cout << "rho(U) - " << C << "\n" << std::flush;
			C.inverse_rho();
			std::cout << "U.inverse_rho() - " << C << "\n" << std::flush;
			do {
				C.rho();
				std::cout << C << "\n" << std::flush;
			} while (!C.is_one());
			std::cout << "inverse_rho(C) - " << apply_inverse_rho(C) << "\n\n" << std::flush;

			std::cout << "Testing rho() and inverse_rho() (non-principal ideal)...\n" << std::flush;
			std::cout << A << "\n" << std::flush;
			apply_rho(C, A);
			std::cout << "rho(A) - " << C << "\n" << std::flush;
			C.inverse_rho();
			std::cout << "inverse_rho() - " << C << "\n" << std::flush;
			do {
				C.rho();
				std::cout << C << "\n" << std::flush;
			} while (!C.is_equal(A));
			std::cout << "inverse_rho(C) - " << apply_inverse_rho(C) << "\n\n" << std::flush;
		}

		// test principality and equivalence functions
		std::cout << "testing principality test:\n" << std::flush;
		C.assign_principal(70000, 22000);
		std::cout << "P = " << C << "\n" << std::flush;
		testbool = C.is_principal(dist);
		std::cout << "principal? ";
		if (testbool)
			std::cout << "true, dist = " << dist << std::flush;
		else
			std::cout << "false, dist = " << dist << std::flush;
		if (C.verify_principal(dist, testbool))
			std::cout << " --->correct!\n\n" << std::flush;
		else {
			std::cout << " --->WRONG!\n\n" << std::flush;
			return 1;
		}

		std::cout << "A = " << A << "\n" << std::flush;
		testbool = A.is_principal(dist);
		std::cout << "principal? ";
		if (testbool)
			std::cout << "true, dist = " << dist << std::flush;
		else
			std::cout << "false, dist = " << dist << std::flush;
		if (A.verify_principal(dist, testbool))
			std::cout << " --->correct!\n\n" << std::flush;
		else {
			std::cout << " --->WRONG!\n\n" << std::flush;
			return 1;
		}

		std::cout << "testing equivalence test:\n" << std::flush;
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "B = " << B << "\n" << std::flush;
		std::cout << "C = " << apply_rho(A) << "\n" << std::flush;
		testbool = B.is_equivalent(A, dist);
		std::cout << "B equivalent to A? ";
		if (testbool)
			std::cout << "true, dist = " << dist << std::flush;
		else
			std::cout << "false, dist = " << dist << std::flush;
		if (B.verify_equivalent(A, dist, testbool))
			std::cout << " --->correct!\n" << std::flush;
		else {
			std::cout << " --->WRONG!\n" << std::flush;
			return 1;
		}

		C = apply_rho(A);
		testbool = A.is_equivalent(C, dist);
		std::cout << "A equivalent to C? ";
		if (testbool)
			std::cout << "true, dist = " << dist << std::flush;
		else
			std::cout << "false, dist = " << dist << std::flush;
		if (A.verify_equivalent(C, dist, testbool))
			std::cout << " --->correct!\n" << std::flush;
		else {
			std::cout << " --->WRONG!\n" << std::flush;
			return 1;
		}

		C = apply_rho(A);
		testbool = C.is_equivalent(A, dist);
		std::cout << "C equivalent to A? ";
		if (testbool)
			std::cout << "true, dist = " << dist << std::flush;
		else
			std::cout << "false, dist = " << dist << std::flush;
		if (C.verify_equivalent(A, dist, testbool))
			std::cout << " --->correct!\n\n" << std::flush;
		else {
			std::cout << " --->WRONG!\n\n" << std::flush;
			return 1;
		}

		// test order, DL, and subgroup functions
		std::cout << "\nTesting power...\n" << std::flush;
		x.assign(41);
		power(C, A, x);
		std::cout << "C = A^" << x << " = " << C << "\n\n" << std::flush;

		std::cout << "testing order, DL, and subgroup algorithms...\n" << std::flush;
		std::cout << "order of " << A << " = " << A.order_in_CL() << "\n" << std::flush;
		std::cout << "order of " << B << " = " << B.order_in_CL() << "\n" << std::flush;
		if (C.DL(A, dlog, dist))
			std::cout << "log_A C = " << dlog << "\n" << std::flush;
		else
			std::cout << "no log_A C:  order(A) = " << dlog << "\n" << std::flush;
		if (C.DL(B, dlog, dist))
			std::cout << "log_B C = " << dlog << "\n" << std::flush;
		else
			std::cout << "no log_B C:  order(B) = " << dlog << "\n" << std::flush;
		if (i < 2) {
			Gvec[0] = A;
			Gvec[1] = B;
			Gvec[2] = A*B;
			S = subgroup(Gvec);
			std::cout << "< " << Gvec << " >  = " << S << "\n\n" << std::flush;
		}
		else
			std::cout << "\n" << std::flush;

		// repeat, using knowledge of h
		std::cout << "testing order algorithm (with knowledge of h)...\n" << std::flush;
		QO.class_number();
		std::cout << "order of " << A << " = " << A.order_in_CL() << "\n" << std::flush;
		std::cout << "order of " << B << " = " << B.order_in_CL() << "\n" << std::flush;
	}

	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
