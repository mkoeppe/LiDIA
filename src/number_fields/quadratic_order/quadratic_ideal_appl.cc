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


#include	"LiDIA/quadratic_order.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	quadratic_order QO, QO1(-3299), QO2(-82475), QO3(241), QO4(2169), QO5, rm;
	bigint D;
	quadratic_ideal A, B, C, E, F, U;
	base_vector< quadratic_ideal > Gvec;
	base_vector< bigint > S;
	bigint p, x, y, dlog;
	int i;

	Gvec.set_mode(EXPAND);

	// test whether ideals are valid before an order is initialized
	std::cout << "No order initialized:  D = " << A.discriminant() << "\n" << std::flush;
	std::cout << "Which order:\n" << A.which_order() << "\n" << std::flush;

	// perform tests on 4 selected discriminants
	for (i = 0; i < 4; ++i) {
		std::cout << "\n==================================\n" << std::flush;

		if (i == 0)
			QO = QO1;
		else if (i == 1)
			QO = QO2;
		else if (i == 2)
			QO = QO3;
		else
			QO = QO4;

		// test assign one function
		if (i % 2)
			U.assign_one();
		else
			U.assign_one(QO);
		std::cout << "(1) = " << U << "\n" << std::flush;

		// test whether the quadratic order is correctly set
		std::cout << "Discriminant = " << U.discriminant() << "\n" << std::flush;
		std::cout << "Which order:\n" << which_order(U) << "\n" << std::flush;
		if (U.discriminant() != QO.discriminant()) {
			std::cout << "ERROR:  discriminants not equal!\n" << std::flush;
			return 1;
		}

		// generate a prime ideal and test conversions
		p = 3;
		while (!generate_prime_ideal(A, p))
			p = next_prime(p);
		std::cout << "A = " << A << ", norm = " << A.norm() << "\n" << std::flush;
		std::cout << "A as qi_class - " << qi_class(A) << "\n" << std::flush;
		std::cout << "A as qi_class_real - " << qi_class_real(A) << "\n" << std::flush;
		std::cout << "A as quadratic_form - " << quadratic_form(A) << "\n" << std::flush;

		A.ring_of_multipliers(rm);
		std::cout << "ring of multipliers:\n" << std::flush;
		std::cout << rm << "\n\n" << std::flush;

		// generate another prime ideal and test conversions
		p = next_prime(p);
		while (!generate_prime_ideal(B, p, QO))
			p = next_prime(p);
		std::cout << "B = " << B << "\n" << std::flush;
		std::cout << "B.a = " << B.get_a() << "\n" << std::flush;
		std::cout << "B.b = " << B.get_b() << "\n" << std::flush;
		std::cout << "B.c = " << B.get_c() << "\n" << std::flush;
		std::cout << "B.q = " << B.get_q() << "\n" << std::flush;
		std::cout << "B.discriminant = " << B.discriminant() << "\n" << std::flush;
		std::cout << "B.smallest_rational = " << B.smallest_rational() << "\n" << std::flush;
		std::cout << "B.norm = " << B.norm() << "\n" << std::flush;
		if (B.is_integral())
			std::cout << "integral? true\n" << std::flush;
		else
			std::cout << "integral? false\n" << std::flush;
		std::cout << "B.conductor = " << B.conductor() << "\n" << std::flush;
		if (B.is_invertible())
			std::cout << "invertible? true\n" << std::flush;
		else
			std::cout << "invertible? false\n" << std::flush;
		std::cout << "conjugate = " << get_conjugate(B) << "\n" << std::flush;
		std::cout << "B*conjugate(B) = " << B*get_conjugate(B) << "\n" << std::flush;
		std::cout << "inverse = " << -B << "\n" << std::flush;
		C = inverse(B);
		std::cout << "norm(B^-1) = " << C.norm() << "\n" << std::flush;
		std::cout << "B*inverse(B) = " << B*inverse(B) << "\n" << std::flush;
		if (B.is_reduced())
			std::cout << "reduced? true\n" << std::flush;
		else
			std::cout << "reduced? false\n" << std::flush;
		std::cout << "B as qi_class - " << qi_class(B) << "\n" << std::flush;
		std::cout << "B as qi_class_real - " << qi_class_real(B) << "\n" << std::flush;
		std::cout << "B as quadratic_form - " << quadratic_form(B) << "\n" << std::flush;
		std::cout << "B (using constructor) - " << quadratic_ideal(B.get_a(), B.get_b(), B.get_q()) << "\n" << std::flush;
		std::cout << "ring of multipliers:\n" << std::flush;
		B.ring_of_multipliers(QO5);
		std::cout << QO5 << "\n\n" << std::flush;

		// test arithmetic operations
		std::cout << "Testing multiplication...\n" << std::flush;
		std::cout << "A*B = " << A * B << "\n" << std::flush;
		C = A*B;
		C.reduce();
		std::cout << "reduced = " << C << "\n" << std::flush;
		std::cout << "Testing square...\n" << std::flush;
		square(C, A);
		std::cout << "A*A = " << C << "\n" << std::flush;
		C.reduce();
		std::cout << "reduced = " << C << "\n" << std::flush;
		std::cout << "Testing division...\n" << std::flush;
		std::cout << "A/B = " << A / B << "\n" << std::flush;
		C = A/B;
		C.reduce();
		std::cout << "reduced = " << C << "\n\n" << std::flush;

		// test rho functions
		if (U.discriminant() > 0) {
			std::cout << "Testing rho() and inverse_rho()...\n" << std::flush;
			std::cout << U << "\n" << std::flush;
			apply_rho(C, U);
			std::cout << "C = rho(U) - " << C << "\n" << std::flush;
			C.inverse_rho();
			std::cout << "C.inverse_rho() - " << C << "\n" << std::flush;
			do {
				C.rho();
				std::cout << C << "\n" << std::flush;
			} while (!C.is_one());
			std::cout << "inverse_rho(C) - " << apply_inverse_rho(C) << "\n\n" << std::flush;
		}

		std::cout << "Testing rho() and inverse_rho() (non-reduced ideal)...\n" << std::flush;
		power(C, A, 41);
		std::cout << "A^41 = " << C << "\n" << std::flush;
		C.reduce();
		std::cout << "C reduced = " << C << "\n\n" << std::flush;
		power(C, A, 41);
		std::cout << "A^41 = " << C << "\n" << std::flush;
		while (!C.is_reduced()) {
			C.rho();
			std::cout << C << "\n" << std::flush;
		}
		std::cout << "\n" << std::flush;

		power(C, A, 41);
		std::cout << "A^41 = " << C << "\n" << std::flush;
		while (!C.is_reduced()) {
			C.inverse_rho();
			std::cout << C << "\n" << std::flush;
		}
		std::cout << "\n" << std::flush;

		// test principality and equivalence functions
		std::cout << "testing principality test:\n" << std::flush;
		C.assign_principal(70000, 22000);
		std::cout << "P = " << C << "\n" << std::flush;
		if (C.is_principal())
			std::cout << "principal? true\n" << std::flush;
		else
			std::cout << "principal? false\n" << std::flush;
		C.reduce();
		std::cout << "reduced = " << C << "\n\n" << std::flush;

		std::cout << "A = " << A << "\n" << std::flush;
		if (A.is_principal())
			std::cout << "principal? true\n\n" << std::flush;
		else
			std::cout << "principal? false\n\n" << std::flush;

		std::cout << "testing equivalence test:\n" << std::flush;
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "B = " << B << "\n" << std::flush;
		if (QO.is_real())
			C = apply_rho(A);
		else
			power(C, A, A.class_number());
		std::cout << "C = " << C << "\n" << std::flush;
		if (B.is_equivalent(A))
			std::cout << "B equivalent to A? true\n" << std::flush;
		else
			std::cout << "B equivalent to A? false\n" << std::flush;
		if (A.is_equivalent(C))
			std::cout << "A equivalent to C? true\n" << std::flush;
		else
			std::cout << "A equivalent to C? false\n" << std::flush;
		if (C.is_equivalent(A))
			std::cout << "C equivalent to A? true\n\n" << std::flush;
		else
			std::cout << "C equivalent to A? false\n\n" << std::flush;

		// test order, DL, and subgroup functions
		std::cout << "Testing power...\n" << std::flush;
		x.assign(41);
		power(C, A, x);
		std::cout << "C = A^" << x << " = " << C << "\n" << std::flush;
		E.assign(C);
		E.reduce();
		std::cout << "reduced = " << E << "\n\n" << std::flush;

		std::cout << "testing order, DL, and subgroup algorithms...\n" << std::flush;
		std::cout << "order of " << A << " = " << A.order_in_CL() << "\n" << std::flush;
		std::cout << "order of " << B << " = " << B.order_in_CL() << "\n" << std::flush;
		if (C.DL(A, dlog))
			std::cout << "log_A C = " << dlog << "\n" << std::flush;
		else
			std::cout << "no log_A C:  order(A) = " << dlog << "\n" << std::flush;
		if (C.DL(B, dlog))
			std::cout << "log_B C = " << dlog << "\n" << std::flush;
		else
			std::cout << "no log_B C:  order(B) = " << dlog << "\n" << std::flush;
		Gvec[0] = A;
		Gvec[1] = B;
		Gvec[2] = A*B;
		S = subgroup(Gvec);
		std::cout << "< " << Gvec << " >  = " << S << "\n\n" << std::flush;

		// testing regulator, class number and class group
		std::cout << "testing regulator, class number, class group...\n" << std::flush;
		if (QO.is_real())
			std::cout << "R = " << A.regulator() << "\n" << std::flush;
		std::cout << "h = " << A.class_number() << "\n" << std::flush;
		std::cout << "CL = " << A.class_group() << "\n\n" << std::flush;
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
