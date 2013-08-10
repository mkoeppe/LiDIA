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


#include	"LiDIA/qi_class.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	quadratic_order QO;
	qi_class A, B, C, U;
	base_vector< qi_class > Gvec;
	base_vector< bigint > S;
	bigint p, x, y, dlog, D;
	int i;
	bool testbool;

	Gvec.set_mode(EXPAND);

	// test whether qi_class's are valid before an order is initialized
	std::cout << "No order initialized:  D = " << A.discriminant() << "\n" << std::flush;
	std::cout << "Current order:\n" << A.get_current_order() << "\n" << std::flush;

	// perform tests on 2 selected discriminants
	for (i = 0; i < 2; ++i) {
		std::cout << "\n==================================\n" << std::flush;

		if (i == 0)
			QO.assign(14401);
		else {
			string_to_bigint("-9828323860172600203", D);
			QO.assign(D);
		}

		qi_class::set_current_order(QO);

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
		std::cout << "A as qi_class_real - " << qi_class_real(A) << "\n" << std::flush;
		std::cout << "A as quadratic_ideal - " << quadratic_ideal(A) << "\n" << std::flush;
		std::cout << "A as quadratic_form - " << quadratic_form(A) << "\n\n" << std::flush;

		// generate another prime ideal and test conversions and inversion
		p = next_prime(p);
		while (!generate_prime_ideal(B, p))
			p = next_prime(p);
		std::cout << "B = " << B << "\n" << std::flush;
		std::cout << "B.a = " << B.get_a() << "\n" << std::flush;
		std::cout << "B.b = " << B.get_b() << "\n" << std::flush;
		std::cout << "B.c = " << B.get_c() << "\n" << std::flush;
		std::cout << "inverse = " << -B << "\n" << std::flush;
		std::cout << "B*inverse(B) = " << B*inverse(B) << "\n" << std::flush;
		if (!U.is_equal(B*inverse(B))) {
			std::cout << "ERROR:  " << U << " <>" << B*inverse(B) << "!\n" << std::flush;
			return 1;
		}
		std::cout << "B as qi_class_real - " << qi_class_real(B) << "\n" << std::flush;
		std::cout << "B as quadratic_ideal - " << quadratic_ideal(B) << "\n" << std::flush;
		std::cout << "B as quadratic_form - " << quadratic_form(B) << "\n" << std::flush;
		std::cout << "B (using constructor) - " << qi_class(B.get_a(), B.get_b()) << "\n\n" << std::flush;

		// test arithmetic operations
		std::cout << "Testing multiplication...\n" << std::flush;
		std::cout << "A*B = " << A * B << "\n" << std::flush;
		std::cout << "Testing square...\n" << std::flush;
		square(C, A);
		std::cout << "A*A = " << C << "\n" << std::flush;
		if (QO.is_imaginary()) {
			std::cout << "Testing nucomp...\n" << std::flush;
			nucomp(C, A, B);
			std::cout << "A*B = " << C << "\n" << std::flush;
			if (!C.is_equal(A*B)) {
				std::cout << "ERROR:  " << C << " <>" << A*B << "!\n" << std::flush;
				return 1;
			}
			std::cout << "Testing nudupl...\n" << std::flush;
			nudupl(C, A);
			std::cout << "A*A = " << C << "\n" << std::flush;
			if (!C.is_equal(A*A)) {
				std::cout << "ERROR:  " << C << " <>" << A*A << "!\n" << std::flush;
				return 1;
			}
		}
		std::cout << "Testing division...\n" << std::flush;
		std::cout << "A/B = " << A / B << "\n\n" << std::flush;
		if (!U.is_equal(A/A)) {
			std::cout << "ERROR:  " << U << " <>" << A/A << "!\n" << std::flush;
			return 1;
		}

		// for real orders, test rho functions
		if (U.discriminant() > 0) {
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
		C.assign_principal(7, 22);
		std::cout << "P = " << C << "\n" << std::flush;
		if (C.is_principal())
			std::cout << "principal? true\n\n" << std::flush;
		else
			std::cout << "principal? false\n\n" << std::flush;

		std::cout << "A = " << A << "\n" << std::flush;
		if (A.is_principal())
			std::cout << "principal? true\n\n" << std::flush;
		else
			std::cout << "principal? false\n\n" << std::flush;

		std::cout << "testing equivalence test:\n" << std::flush;
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "B = " << B << "\n" << std::flush;
		if (QO.is_real()) {
			C = apply_rho(A);
			std::cout << "C = " << apply_rho(A) << "\n" << std::flush;
		}
		if (B.is_equivalent(A))
			std::cout << "B equivalent to A? true\n" << std::flush;
		else
			std::cout << "B equivalent to A? false\n" << std::flush;

		if (QO.is_real()) {
			if (A.is_equivalent(C))
				std::cout << "A equivalent to C? true\n" << std::flush;
			else
				std::cout << "A equivalent to C? false\n" << std::flush;
			if (C.is_equivalent(A))
				std::cout << "C equivalent to A? true\n" << std::flush;
			else
				std::cout << "C equivalent to A? false\n" << std::flush;
		}
		std::cout << "\n" << std::flush;

		// test order, DL, and subgroup functions
		std::cout << "Testing power...\n" << std::flush;
		x.assign(41);
		power(C, A, x);
		std::cout << "C = A^" << x << " = " << C << "\n\n" << std::flush;

		std::cout << "testing order, DL, and subgroup algorithms...\n" << std::flush;
		x = A.order_in_CL();
		std::cout << "order of " << A << " = " << x << std::flush;
		if (A.verify_order(x))
			std::cout << " -->correct!\n" << std::flush;
		else {
			std::cout << " -->WRONG!\n" << std::flush;
			return 1;
		}

		x = B.order_in_CL();
		std::cout << "order of " << B << " = " << x << std::flush;
		if (B.verify_order(x))
			std::cout << " -->correct!\n" << std::flush;
		else {
			std::cout << " -->WRONG!\n" << std::flush;
			return 1;
		}

		testbool = C.DL(A, dlog);
		if (testbool)
			std::cout << "log_A C = " << dlog << std::flush;
		else
			std::cout << "no log_A C:  order(A) = " << dlog << std::flush;
		if (C.verify_DL(A, dlog, testbool))
			std::cout << " -->correct!\n" << std::flush;
		else {
			std::cout << " -->WRONG!\n" << std::flush;
			return 1;
		}

		testbool = C.DL(B, dlog);
		if (testbool)
			std::cout << "log_B C = " << dlog << std::flush;
		else
			std::cout << "no log_B C:  order(B) = " << dlog << std::flush;
		if (C.verify_DL(B, dlog, testbool))
			std::cout << " -->correct!\n" << std::flush;
		else {
			std::cout << " -->WRONG!\n" << std::flush;
			return 1;
		}

		Gvec[0] = A;
		Gvec[1] = B;
		Gvec[2] = A*B;
		S = subgroup(Gvec);
		std::cout << "< " << Gvec << " >  = " << S << "\n\n" << std::flush;

		// repeat, using knowledge of h
		std::cout << "\ntesting order algorithm (with knowledge of h)...\n" << std::flush;
		QO.class_number();

		x = A.order_in_CL();
		std::cout << "order of " << A << " = " << x << std::flush;
		if (A.verify_order(x))
			std::cout << " -->correct!\n" << std::flush;
		else {
			std::cout << " -->WRONG!\n" << std::flush;
			return 1;
		}

		x = B.order_in_CL();
		std::cout << "order of " << B << " = " << x << std::flush;
		if (B.verify_order(x))
			std::cout << " -->correct!\n\n" << std::flush;
		else {
			std::cout << " -->WRONG!\n\n" << std::flush;
			return 1;
		}
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
