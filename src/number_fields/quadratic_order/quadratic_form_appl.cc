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


#include	"LiDIA/quadratic_form.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/qi_class.h"
#include	"LiDIA/qi_class_real.h"
#include	"LiDIA/quadratic_ideal.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	bigint D;
	quadratic_form A, B, C, F, U;
	base_vector< quadratic_form > Gvec;
	base_vector< bigint > S;
	bigint p, x, y, dlog, N;
	matrix_GL2Z M1, NMAT(1, 0, 0, -1);
	int i;
	sort_vector< pair < bigint, bigint > > Reps;

	Gvec.set_mode(EXPAND);

	// test whether forms are valid before an order is initialized
	std::cout << "No order initialized:  D = " << A.discriminant() << "\n" << std::flush;
	std::cout << "Which order:\n" << A.which_order() << "\n" << std::flush;

	// perform tests on 5 selected discriminants
	for (i = 0; i < 5; ++i) {
		std::cout << "\n==================================\n" << std::flush;

		if (i == 0)
			D = -3299;
		else if (i == 1)
			D = -82475;
		else if (i == 2)
			D = 241;
		else if (i == 3)
			D = 2169;
		else
			D = 169;

		// test assign one function
		U.assign_one(D);
		std::cout << "Unit form = " << U << "\n" << std::flush;

		// test whether the quadratic order is correctly set
		std::cout << "Discriminant = " << U.discriminant() << "\n" << std::flush;
		std::cout << "Which order:\n" << which_order(U) << "\n" << std::flush;
		if (U.discriminant() != D) {
			std::cout << "ERROR:  discriminants not equal!\n" << std::flush;
			return 1;
		}

		// generate a prime ideal and test conversions
		p = 3;
		while (!generate_prime_form(A, p, D))
			p = next_prime(p);
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "A as qi_class - " << qi_class(A) << "\n" << std::flush;
		std::cout << "A as qi_class_real - " << qi_class_real(A) << "\n" << std::flush;
		std::cout << "A as quadratic_ideal - " << quadratic_ideal(A) << "\n" << std::flush;
		std::cout << "definiteness - " << A.definiteness() << "\n" << std::flush;
		std::cout << "A.content = " << A.content() << "\n" << std::flush;
		if (A.is_primitive())
			std::cout << "primitive? true\n" << std::flush;
		else
			std::cout << "primitive? false\n" << std::flush;
		if (A.is_normal())
			std::cout << "normal? true\n" << std::flush;
		else
			std::cout << "normal? false\n" << std::flush;
		C = A;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		C.normalize(M1);
		std::cout << "A normalized = " << C << "\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		C = A;
		C.transform(M1);
		std::cout << "A*M = " << C << "\n" << std::flush;
		if (A.is_reduced())
			std::cout << "reduced? true\n" << std::flush;
		else
			std::cout << "reduced? false\n" << std::flush;
		C = A;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		C.reduce(M1);
		std::cout << "A reduced = " << C << "\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		C = A;
		C.transform(M1);
		std::cout << "A*M = " << C << "\n\n" << std::flush;

		// generate another prime ideal and test conversions
		p = next_prime(p);
		while (!generate_prime_form(B, p, D))
			p = next_prime(p);
		std::cout << "B = " << B << "\n" << std::flush;
		std::cout << "B.a = " << B.get_a() << "\n" << std::flush;
		std::cout << "B.b = " << B.get_b() << "\n" << std::flush;
		std::cout << "B.c = " << B.get_c() << "\n" << std::flush;
		if (A < B)
			std::cout << "A < B? true\n" << std::flush;
		else
			std::cout << "A < B? false\n" << std::flush;
		if (A <= B)
			std::cout << "A <= B? true\n" << std::flush;
		else
			std::cout << "A <= B? false\n" << std::flush;
		if (A > B)
			std::cout << "A > B? true\n" << std::flush;
		else
			std::cout << "A > B? false\n" << std::flush;
		if (A >= B)
			std::cout << "A >= B? true\n" << std::flush;
		else
			std::cout << "A >= B? false\n" << std::flush;
		std::cout << "conjugate = " << get_conjugate(B) << "\n" << std::flush;
		std::cout << "B*conjugate(B) = " << B*get_conjugate(B) << "\n" << std::flush;
		C = B;
		C.reduce();
		std::cout << "B reduced = " << C << "\n" << std::flush;
		std::cout << "B as qi_class - " << qi_class(B) << "\n" << std::flush;
		std::cout << "B as qi_class_real - " << qi_class_real(B) << "\n" << std::flush;
		std::cout << "B as quadratic_ideal - " << quadratic_ideal(B) << "\n" << std::flush;
		std::cout << "B (using constructor) - " << quadratic_form(B.get_a(), B.get_b(), B.get_c()) << "\n\n" << std::flush;

		// test arithmetic operations
		std::cout << "Testing composition...\n" << std::flush;
		std::cout << "A*B = " << A * B << "\n" << std::flush;
		C = A*B;
		C.reduce();
		std::cout << "reduced = " << C << "\n" << std::flush;
		std::cout << "Testing square...\n" << std::flush;
		square(C, A);
		std::cout << "A*A = " << C << "\n" << std::flush;
		C = A*A;
		C.reduce();
		std::cout << "reduced = " << C << "\n" << std::flush;
		if (!A.is_indefinite() && A.is_regular()) {
			std::cout << "Testing nucomp...\n" << std::flush;
			nucomp(C, A, B);
			std::cout << "A*B = " << C << "\n" << std::flush;
			F = A*B;
			F.reduce();
			if (!C.is_equal(F)) {
				std::cout << "ERROR:  " << C << " <>" << F << "!\n" << std::flush;
				return 1;
			}
			std::cout << "Testing nudupl...\n" << std::flush;
			nudupl(C, A);
			std::cout << "A*A = " << C << "\n" << std::flush;
			F = A*A;
			F.reduce();
			if (!C.is_equal(F)) {
				std::cout << "ERROR:  " << C << " <>" << A*A << "!\n" << std::flush;
				return 1;
			}
		}
		std::cout << "Testing division...\n" << std::flush;
		std::cout << "A/B = " << A / B << "\n\n" << std::flush;

		// test rho functions
		if (U.is_indefinite()) {
			std::cout << "Testing rho() and inverse_rho()...\n" << std::flush;
			std::cout << U << "\n" << std::flush;
			C.assign(U);
			C.rho();
			std::cout << "C = rho(U) - " << C << "\n" << std::flush;
			C.inverse_rho();
			std::cout << "C.inverse_rho() - " << C << "\n" << std::flush;
			do {
				C.rho();
				std::cout << C << "\n" << std::flush;
			} while (!C.is_one());
			C.inverse_rho();
			std::cout << "inverse_rho(C) - " << C << "\n\n" << std::flush;

			std::cout << "Testing rho() and inverse_rho() (non-principal ideal)...\n" << std::flush;
			std::cout << A << "\n" << std::flush;
			C = A;
			C.normalize();
			std::cout << C << " (A normalized) \n" << std::flush;
			F = C;
			do {
				C.rho();
				std::cout << C << "\n" << std::flush;
			} while (!C.is_equal(F));
		}
		std::cout << "\n" << std::flush;

		std::cout << "Testing rho() and inverse_rho()...\n" << std::flush;
		power(C, A, 41);
		std::cout << "A^41 = " << C << "\n" << std::flush;
		C.reduce();
		std::cout << "reduced = " << C << "\n" << std::flush;
		power(C, A, 41);
		M1 = matrix_GL2Z(1, 0, 0, 1);
		std::cout << "A^41 = " << C << "\n" << std::flush;
		while (!C.is_reduced()) {
			C.rho(M1);
			std::cout << C << "\n" << std::flush;
		}
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		power(C, A, 41);
		C.transform(M1);
		std::cout << "A^41*M = " << C << "\n\n" << std::flush;

		power(C, A, 41);
		M1 = matrix_GL2Z(1, 0, 0, 1);
		std::cout << "A^41 = " << C << "\n" << std::flush;
		while (!C.is_reduced()) {
			C.inverse_rho(M1);
			std::cout << C << "\n" << std::flush;
		}
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		power(C, A, 41);
		C.transform(M1);
		std::cout << "A^41*M = " << C << "\n\n" << std::flush;

		// test principality and equivalence functions
		std::cout << "testing principality test:\n" << std::flush;
		power(C, A, A.class_number());
		std::cout << "P = " << C << "\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (C.is_principal(M1))
			std::cout << "principal? true\n" << std::flush;
		else
			std::cout << "principal? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		C = U;
		C.transform(M1);
		std::cout << "U*M = " << C << "\n\n" << std::flush;

		M1 = matrix_GL2Z(1, 0, 0, 1);
		std::cout << "A = " << A << "\n" << std::flush;
		if (A.is_principal(M1))
			std::cout << "principal? true\n" << std::flush;
		else
			std::cout << "principal? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		C = U;
		C.transform(M1);
		std::cout << "U*M = " << C << "\n\n" << std::flush;

		std::cout << "testing equivalence test:\n" << std::flush;
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "B = " << B << "\n" << std::flush;
		C = A;
		if (D > 0)
			C.rho();
		else
			power(C, C, A.class_number());
		std::cout << "C = " << C << "\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (B.is_equivalent(A, M1))
			std::cout << "B equivalent to A? true\n" << std::flush;
		else
			std::cout << "B equivalent to A? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		F = A;
		F.transform(M1);
		std::cout << "A*M = " << F << "\n\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (A.is_equivalent(C, M1))
			std::cout << "A equivalent to C? true\n" << std::flush;
		else
			std::cout << "A equivalent to C? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		F = C;
		F.transform(M1);
		std::cout << "C*M = " << F << "\n\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (C.is_equivalent(A))
			std::cout << "C equivalent to A? true\n" << std::flush;
		else
			std::cout << "C equivalent to A? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		F = A;
		F.transform(M1);
		std::cout << "A*M = " << F << "\n\n" << std::flush;

		std::cout << "testing proper equivalence test:\n" << std::flush;
		std::cout << "A = " << A << "\n" << std::flush;
		std::cout << "B = " << B << "\n" << std::flush;
		C = A;
		if (D > 0)
			C.rho();
		else
			power(C, C, A.class_number());
		std::cout << "C = " << C << "\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (A.is_prop_equivalent(B, M1))
			std::cout << "A properly equivalent to B? true\n" << std::flush;
		else
			std::cout << "A properly equivalent to B? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		F = B;
		F.transform(M1);
		std::cout << "B*M = " << F << "\n\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (A.is_prop_equivalent(C, M1))
			std::cout << "A properly equivalent to C? true\n" << std::flush;
		else
			std::cout << "A properly equivalent to C? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		F = C;
		F.transform(M1);
		std::cout << "C*M = " << F << "\n\n" << std::flush;
		M1 = matrix_GL2Z(1, 0, 0, 1);
		if (C.is_prop_equivalent(A, M1))
			std::cout << "C properly equivalent to A? true\n" << std::flush;
		else
			std::cout << "C properly equivalent to A? false\n" << std::flush;
		std::cout << "transformation = \n" << M1 << "\n" << std::flush;
		F = A;
		F.transform(M1);
		std::cout << "A*M = " << F << "\n\n" << std::flush;

		// test representations
		std::cout << "Testing representations...\n" << std::flush;
		N.assign(A.eval(-101, 33));
		A.representations(Reps, N);
		std::cout << "Representations of " << N << " by " << A << ":" << std::endl;
		std::cout << Reps << std::endl << std::endl;

		// test order, DL, and subgroup functions
		std::cout << "Testing power...\n" << std::flush;
		x.assign(41);
		power(C, A, x);
		std::cout << "C = A^" << x << " = " << C << "\n" << std::flush;
		C.reduce();
		std::cout << "reduced = " << C << "\n\n" << std::flush;
		power(C, A, x);

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
		if (A.is_indefinite())
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
