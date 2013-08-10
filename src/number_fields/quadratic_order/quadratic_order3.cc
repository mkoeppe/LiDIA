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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/prime_list.h"
#include	"LiDIA/quadratic_order.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//#define DEBUG
//#define DEBUG_WITH_CHECK


//
//
//  Approximation of l(n, Delta)
//
//


void quadratic_order
::bachs_L1_error_triplet(const bigint & Q,
			 xbigfloat & A,
			 xbigfloat & B) const
{
	debug_handler("mm_quadratic_order", "bachs_L1_error_triplet");

	if (Q >= 1000000) {
		A = xbigfloat(6.246);
		B = xbigfloat(16.271);
	}
	else if (Q >= 500000) {
		A = xbigfloat(6.269);
		B = xbigfloat(16.409);
	}
	else if (Q >= 100000) {
		A = xbigfloat(6.338);
		B = xbigfloat(17.031);
	}
	else if (Q >= 50000) {
		A = xbigfloat(6.378);
		B = xbigfloat(17.397);
	}
	else if (Q >= 10000) {
		A = xbigfloat(6.510);
		B = xbigfloat(18.606);
	}
	else if (Q >= 5000) {
		A = xbigfloat(6.593);
		B = xbigfloat(19.321);
	}
	else if (Q >= 1000) {
		A = xbigfloat(6.897);
		B = xbigfloat(21.528);
	}
	else if (Q >= 500) {
		A = xbigfloat(7.106);
		B = xbigfloat(22.845);
	}
	else if (Q >= 100) {
		A = xbigfloat(7.962);
		B = xbigfloat(27.145);
	}
	else if (Q >= 50) {
		A = xbigfloat(8.628);
		B = xbigfloat(29.587);
	}
	else if (Q >= 10) {
		A = xbigfloat(12.170);
		B = xbigfloat(38.831);
	}
	else {
		A = xbigfloat(16.397);
		B = xbigfloat(47.183);
	}
}


//
// number_of_terms:
//
//  returns n such that C(n) <= 2^{-k}, where
//  C(n) = A + B * log|Delta| / (sqrt(n) * log n)
//
//  see Buchmann/Maurer: Approximate evaluation of L(1,chi).
//

bigint quadratic_order
::number_of_terms (const bigint & Delta, long k, int info) const
{
	xbigfloat d, A, B, C, S, L;
	bigint    l, m, n;
	long	     bm;

	// d = log(abs(Delta), k+7);
	// log(d, abs(Delta));
	log(d, abs(Delta), k+7);

	// A = Amax, B = Bmax
	bachs_L1_error_triplet (0, A, B);

	if (info > 2) {
		std::cout << "Amax = " << A << std::endl;
		std::cout << "Bmax = " << B << std::endl;
	}

	// C = 2^k * (A + B * d) + 0.5;
	C = A + B * d;
	shift_left (C, C, k);
	C += xbigfloat(0.5);

	// n = ceil(C)^2;
	ceil(n, C);
	square (n, n);

	l = 2;

	if (info > 2)
		std::cout << "n0 = " << n << std::endl;

	do {
		m = l + (n-l)/2;

		if (info > 2) {
			std::cout << "============================" << std::endl;
			std::cout << "Entering loop with (l, m, n) = ";
			std::cout << "(" << l << ", ";
			std::cout << m << ", ";
			std::cout << n << ")" << std::endl;
			std::cout << std::endl;
		}

		// S = sqrt(m, 3+ ceil(b(m)/2) + b(b(m)));
		// sqrt(S,m);
		bm = b_value(m);
		sqrt(S, xbigfloat(m), 3+ (bm+1) >> 1 + b_value(bm));

		// L = log (m, 3+ ceil(b(m)/2));
		//log(L,m);
		log(L, xbigfloat(m), 3+ (bm+1) >> 1);


		// C = 2^k * (A(m) + B(m) * d);
		bachs_L1_error_triplet (m, A, B);

		if (info > 2) {
			std::cout << "A(m) = " << A << std::endl;
			std::cout << "B(m) = " << B << std::endl;
		}

		C = A + B * d;
		shift_left (C, C, k);

		if ((S*L - 1) >= C)
			n = m;
		else
			l = m+1;
	} while (l < n);

#if 0
	// Relative 1-approximation to C(n)
	log(d, abs(Delta), 7);
	bachs_L1_error_triplet (n, A, B);
	C = A + B * d;

	sqrt(S, bigfloat(n), 7);
	log (L, bigfloat(n), 7);
	S *= L;
	divide (C, C, S, 7);

	std::cout << "C = " << C << std::endl;

	d.assign_one();
	shift_right(d, d, k);
	d -= C;

	std::cout << "lower bound on the difference is " << d << std::endl;
	std::cout << "b(lower bound) = " << b_value(d) << std::endl;
#endif

	return n;
}


//
// ell
//
//  returns absolute k-approximation xl to l(n,Delta), where
//
// l(n,Delta) = \sum_{p < 2n-1} wpn * log(p/ p-chi_{Delta}(p))
// and
//         1, if p < n
// wpn =   \sum_{j=p-n+1}^{n-1} a_j(n), if n \le p < 2n-1
//
// and
//
//  a_i(n) = (n+i) log(n+i) / \sum_{j=0}^{n-1} (n+j) log(n+j),
//         0 \leq i \leq n-1.
//
//  see Buchmann/Maurer: Approximate evaluation of L(1,chi).
//

void quadratic_order
::ell (xbigfloat & xl,
       unsigned long n,
       const bigint & Delta,
       long k,
       int info) const
{
	lidia_size_t  i, N;
	long		 b, bN, bln, blp, s, t, t1;
	unsigned long p, q, j;
	bigfloat      l, w_num, w_den, w, oneover2, h, x, nj;
	int		 old_rnd_mode, kron, kbN9;
	int           use_double_for_w, use_double_for_l;
	bigint        big_t;
	double        dnj, dh, dw_den = 0.0, dw_num, dw = 0.0, dx, dl;

	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode(MP_EXACT);

	//
	//   initialize prime table
	//

	prime_list P(2*n-1);
	N = P.get_number_of_primes();
	bN = b_value(N);

	if (info > 2) {
		std::cout << "Absolute " << k << " approximation required." << std::endl;
		std::cout << "Primes computed up to 2n-1 = " << 2*n-1 << std::endl;
		std::cout << "Number of primes N = " << N << std::endl;
		std::cout << "b(N) = " << bN << std::endl;
	}

	//
	//  Approximate denominator of the weights. The precision is chosen such
	//  that each approximation of the denominator of w_p(n) for a
	//  single p can be deduced from that approximation by truncation.
	//

	//
	// We need a relative k+b(N)+bl(n)+3 approximation
	// w_den to the denominator of w_p(n).
	// First, we check whether we can use doubles for
	// the whole approximation, i.e., we test whether
	// max{0,k+b(N)+bl(n)}+b(n-1)+10 <= 53
	//

	// use_double_for_w = 0; don't use double approximation
	// use_double_for_w = 1; use double approximation for products
	//                       in w_num and w_den, but add to xbigfloat
	// use_double_for_w = 2; compute w_num and w_den as double

	if (n == 2)
		bln = 0;
	else
		bln = -b_value(n-2)+1;

	big_t = k;
	big_t += bN;
	big_t += bln;

	if (big_t < 0)
		big_t = 0;

	big_t += b_value(n-1);
	big_t += 10;

	if (big_t <= 53) {
		if (info > 2) {
			std::cout << "Computing w_den as double (rounding to " << std::endl;
			std::cout << big_t << " bits) required for each step." << std::endl;
		}
		use_double_for_w = 2;
		dw_den = 0.0;

		for (j = 0; j < n; j++) {
			// no error by adding and casting
			dnj = static_cast<double>(n+j);

			// relative error <= 2^-53 for log
			dh = std::log(dnj);

			// relative error <= 2^-53 for
			// multiplication and addition
			dw_den += (dnj * dh);
		}
		// no error by casting
		w_den = bigfloat(dw_den);

	}
	else // big_t <= 53
	{
		if (info > 2) {
			std::cout << "Computing w_den as double not possible (rounding" << std::endl;
			std::cout << "to " << big_t << " bits) required for each step." << std::endl;
		}

		// b = b(log(n,1)-0.5)
		log_absolute (h, n, 1);
		oneover2.assign_one();
		shift_right(oneover2, oneover2, 1);
		b = b_value(h - oneover2);

		// t = max{0,k+b(N)+bl(n)+8-b}
		big_t = k;
		big_t += bN;
		big_t += 8;
		big_t += bln;

		if (big_t.is_lt_zero())
			big_t.assign_zero();

		//
		// We need absolute t-approximations to (n+j) * log(n+j).
		// Check for overflow.
		//

		if (big_t.longify(t)) {
			lidia_error_handler ("quadratic_order3.c::Ell",
					     "Precision overflow.");
			return;
		}

		//
		// Check whether double precision is sufficient. i.e.,
		//   t+t1+2 <= 53 (for IEEE double), where
		// t1 = b(floor(7 (b(n)+1)/10) + 1).
		//

		// w_den = sum(j=0,n) (n+j) * log(n+j,max{0,k+b(N)+bl(n)+8-b})
		w_den.assign_zero();

		t1 = b_value(((7* (b_value(n)+1)) / 10) + 1);
		big_t += t1;
		big_t += 2;

		if (big_t <= 53) {
			if (info > 2) {
				std::cout << "Each summand w_den computed as double" << std::endl;
				std::cout << "(rounding to " << big_t << " bits required)." << std::endl;
			}

			use_double_for_w = 1;
			for (j = 0; j < n; j++) {
				// no error by adding and casting
				dnj = static_cast<double>(n+j);

				// relative error <= 2^-53 for log
				dh = std::log(dnj);

				// relative error <= 2^-53 for multiplication
				// and no error by casting and adding
				w_den += bigfloat(dnj * dh);
			}
		}
		else {
			// big_t <= 53 (2)
			if (info > 2) {
				std::cout << "Using bigfloats only for w_den" << std::endl;
				std::cout << "(rounding to " << big_t << " bits required)." << std::endl;
				std::cout << "Computing each log(n+j) with absolute precision t = ";
				std::cout << t << std::endl;
			}

			use_double_for_w = 0;
			for (j = 0; j < n; j++) {
				// no error by adding and casting
				nj = bigfloat(n+j);

				// absolut error < 2^-t for log
				log_absolute(h, nj, t);

				// no error for multiplication and adding
				w_den += (nj * h);
			}
		}
	}

	if (info > 2)
		std::cout << "w_den = " << w_den << std::endl;


	//
	// ell(n, Delta): part for n <= p < 2n-1
	//
	// -) approximate the numerator w_num for each p,
	// -) approximate the denominator w_den for each p using the
	//    w_den approximation computed above,
	// -) compute the weight w for each p as quotient of w_num and w_den,
	// -) approximate log(p/(p-chi(p))),
	// -) multiply the logarithm h by the weight w and
	// -) add the product to the sum l.
	//


	// use_double_for_l = 0:  no usage of doubles
	// use_double_for_l = 1:  use doubles for the whole computation of l

	// check whether k+bN+9 <= 53.

	big_t = k;
	big_t += bN;
	big_t += 9;

	if (big_t <= 53)
		kbN9 = 1;
	else
		kbN9 = 0;


	// Decide whether the complete sum w(p)*log(p/(p-chi)) can be computed
	// as double. To be able to do this we require the following:
	//
	// 1) w(p) must be computable as double (use_double_for_w = 2)
	// 2) Each log(p/(p-chi)) must be computable as double
	// 3) k+2b(N)+4 <= 53

	// Check 1)
	if (use_double_for_w == 2)
		use_double_for_l = 1;
	else
		use_double_for_l = 0;

	// Check 2), i.e.,
	// b(p_max+2) + 6, k+bN+9 <= 53

	if (!kbN9 || b_value(P[N-1]+2)+6 > 53)
		use_double_for_l = 0;

	// Check 3)
	big_t = k;
	big_t += bN;
	big_t += bN;
	big_t += 4;

	if (big_t > 53)
		use_double_for_l = 0;


	l.assign_zero();
	dl = 0.0;
	w_num.assign_zero();
	dw_num = 0.0;

	q = n;
	i = N-1;
	s = k + bN;

	if (info > 2) {
		std::cout << "s = " << s << std::endl;
		if (use_double_for_l) {
			std::cout << "Using doubles for complete sum." << std::endl;
			std::cout << "k+b(N)+9 = " << k+bN+9 << std::endl;
			std::cout << "b_value(pmax+2)+6 = " << b_value(P[N-1]+2)+6 << std::endl;
			std::cout << "k+2b(N)+4 = " << k+2*bN+4 << std::endl;
		}
	}

	for (p = P[i]; p >= n; p = P[i]) {
		// inefficient, we should treat the case n = 2 separately.
		if (p == 2)
			blp = 0;
		else
			blp = -b_value(p-2)+1;

		//
		// compute the numerator w_num for p
		//
		if (use_double_for_w == 2) {
			for (j = q-1; j >= p-n+1; j--) {
				// no error by adding and casting
				dnj = static_cast<double>(n+j);

				// relative error <= 2^-53 for log
				dh = std::log(dnj);

				// relative error <= 2^-53 for
				// multiplication and addition
				dw_num += (dnj * dh);
			}
		}
		else if (use_double_for_w == 1) {
			for (j = q-1; j >= p-n+1; j--) {
				// no error by adding and casting
				dnj = static_cast<double>(n+j);

				// relative error <= 2^-53 for log
				dh = std::log(dnj);

				// relative error <= 2^-53 for multiplication
				// and no error by casting and adding
				add(w_num, w_num, bigfloat(dnj * dh));
			}
		}
		else {
			for (j = q-1; j >= p-n+1; j--) {
				// no error by adding and casting
				nj = bigfloat(n+j);

				// absolut error < 2^-t for log
				log_absolute(h, nj, t);

				// no error for multiplication and adding
				add(w_num, w_num, nj * h);
			}
		}

		//
		// compute the quotient w = w_num / w_den for p
		//
		if (use_double_for_w == 2)
			dw = dw_num / dw_den;
		else {
			if (s+7 < -blp)
				divide(w, w_num, w_den, 1);
			else {
				divide(w, w_num, w_den, (s+7)+blp);
				w.cut(w, (s+6)+blp);
			}
		}

		//
		// approximate log(p/(p-chi)), multiply by w and add to sum.
		//
		kron = kronecker(Delta, bigint(p));
		if (kron == 1 || kron == -1) {
			// check whether log(p/(p-chi)) can be computed
			// as double
			if (kbN9 && b_value(p+2)+6 <= 53) {
				dx = static_cast<double>(p) / static_cast<double>(p-kron);
				dh = std::log(dx);

				// check how to multiply and how to add to sum:
				if (use_double_for_w == 2)
					if (use_double_for_l)      // mutiply with and add to double
						dl += (dw * dh);
					else
						add(l, l, bigfloat(dw*dh)); // multiply with double, but
				// add to bigfloat
				else
					add(l, l, w*bigfloat(dh)); // multiply with and add
				// to bigfloat
			}
			else {
				// compute log as bigfloat. Note, that in this case
				// the conditions above guarantee use_double_for_l == 0.
				divide(x, bigfloat(p), bigfloat(p-kron), s+5);
				log_absolute(h, x, s+4);

				if (use_double_for_w == 2)
					add(l, l, bigfloat(w)*h);
				else
					add(l, l, w * h);
			}
		}

		q = p-n+1;
		i--;
	}


	if (info > 2)
		std::cout << "PART 2" << std::endl;

	//
	//
	//  ell(n, \Delta) : part for 2 <= p < n
	//
	//
	do {
		p = P[i];
		kron = kronecker(Delta, bigint(p));
		if (kron == 1 || kron == -1) {
			// check whether log(p/(p-chi)) can be computed
			// as double
			if (kbN9 && b_value(p+2)+6 <= 53) {
				dx = static_cast<double>(p) / static_cast<double>(p-kron);
				dh = std::log(dx);

				if (use_double_for_l)
					dl += dh;
				else
					add(l, l, bigfloat(dh));
			}
			else {
				divide(x, bigfloat(p), bigfloat(p-kron), s+3);
				log_absolute(h, x, s+2);
				add(l, l, h);
			}
		}
		i--;
	} while (i >= 0);

	if (use_double_for_l)
		l = bigfloat(dl);

	l.cut(l.exponent()+b_value(l.mantissa())+k+1);
	xl = l;
	bigfloat::set_mode (old_rnd_mode);
}



//
//
//   Approximation of L(1, chi)
//
//

 //
 // relative_L1chi_approximation
 //
 //   returns a relative k-approximation L to L(1,chi_Delta)
 //
 // see Buchmann/Maurer: Approximate evaluation of L(1,chi).
 //

xbigfloat quadratic_order
::relative_L1chi_approximation (

	long k,
	int info) const
{
	xbigfloat L;
	this->relative_L1chi_approximation(L, k, info);
	return L;
}

void quadratic_order
::relative_L1chi_approximation (

	xbigfloat & L,
	long k,
	int info) const
{
	bigint   N;
	long     n; // Would like to use unsigned long, but there is
	// no ulongify for bigints.
	xbigfloat l;
//long t;


	// check precision
	if (k < 0) {
		warning_handler("relative_L1chi_approximation",
				"precision negative; computing an absolute 0 approx.");
		k = 0;
	}

	// compute the number of terms
	N = number_of_terms (Delta, k+2);

	if (info > 2) {
		std::cout << "number of terms = " << N;
		std::cout << " (" << N.bit_length() << " bits) " << std::endl;
	}

	if (N.longify(n)) {
		lidia_error_handler ("relative_L1chi_approximation",
				     "number of terms to large for long.");
		return;
	}

	if (info > 2)
		std::cout << "number of terms (as long) n = " << n << std::endl;

#if 0
	std::cout << "Increase k by->";
	std::cin >> t;
	k += t;

	std::cout << "New k = " << k << std::endl;
#endif

	// Absolute k+3 approximation to l(n,Delta)
	ell (l, n, Delta, k+3, 1);

	if (info > 2) {
		std::cout << "absolute k+4-approximation l to l(n, Delta) is ";
		std::cout << l << std::endl;
	}

	// Relative k+2 approximation L to exp(l)
	exp (L, l, k+3);

	if (info > 2)
		std::cout << "relative k+2-approximation to exp(l) is " << L << std::endl;

	// Truncate mantissa to k+3 bits to obtain a
	// relative k-approximation to L(1,chi_Delta).
	L.truncate(k+3);

	if (info > 2) {
		std::cout << "relative k-approximation L to L(1, chi_Delta) is ";
		std::cout << L << std::endl;
	}
}


//
// Determine approximation L to L(1,chi_Delta) with
//
// |L/exp(l(n,Delta))-1| < 2^{-7}
// C(n) <= 1/4
//
// (used for class number verification and computing
//  lower regulator bounds, see below).
//

void quadratic_order
::approximate_exp_l_n_delta(
	xbigfloat & L,
	const bigint & Delta,
	int info) const
{
	long n;
	bigint N;
	xbigfloat l;
	timer t;

	if (info > 2)
		t.start_timer();

	N = number_of_terms (Delta, 2);

	if (info > 2) {
		std::cout << "number of terms = " << N;
		std::cout << " (" << N.bit_length() << " bits) " << std::endl;
	}
	if (info > 2) {
		t.stop_timer();
		std::cout << "(Timing: ";
		t.print();
		std::cout << ")" << std::endl;
	}

	if (N.longify(n)) {
		lidia_error_handler ("approximate_sqrt_2_exp_l_n_delta",
				     "number of terms to large for long.");
		return;
	}

	if (info > 2) {
		std::cout << "number of terms (as long) n = " << n << std::endl;
	}
	if (info > 2) t.start_timer();
	ell(l, n, Delta, 8, info);

	if (info > 2) {
		t.stop_timer();
		std::cout << "(Timing for ell: ";
		t.print();
		std::cout << ")" << std::endl;
	}

	exp(L, l, 9);
}



//
//
//

void quadratic_order
::test_L1chi ()
{
	long	     k;
	bigint    Delta;
	xbigfloat l;

	quadratic_order O;

	std::cout << "discriminant Delta = ";
	std::cin >> Delta;
	std::cout << std::endl;

	O.assign(Delta);

	std::cout << "precision k = ";
	std::cin >> k;
	std::cout << std::endl;

	O.relative_L1chi_approximation(l, k, 2);

	std::cout << "l = " << l << std::endl;

#if 0
	lidia_size_t n;

	std::cout << "n = ";
	std::cin >> n;

	prime_list P(n);
	lidia_size_t i;
	lidia_size_t N;

	N = P.get_number_of_primes();

	for (i = 0; i < N; i++) {
		std::cout << P[i] << " ";
		std::cout << static_cast<int>(P.diff[i]) << " ";
	}

	std::cout << "===================" << std::endl;

	for (i = N-1; i >= 0; i--) {
		std::cout << P[i] << " ";
	}
#endif

#if 0
	lidia_size_t j, n;
	double nj, w_den;

	w_den = 0;
	n = 5184;

	for (j = 0; j < n; j++) {
		w_den += (static_cast<double>(n+j) * log(n+j));
	}

	std::cout << "w_den = " << w_den << std::endl;
#endif

}




//
//
//  Verification of the class number
//
//


 //
 // imaginary case: verify_h
 //
 //  Delta is discriminant,
 //  h is integer multiple of the class number,
 //
 //  return true iff h is the class number of O(Delta).
 //

 // for convenience
bool quadratic_order
::verify_h(const bigint & Delta, const bigint & h, int info)
{
	xbigfloat L;
	return (verify_h(L, Delta, h, info));
}

bool quadratic_order
::verify_h(
	xbigfloat & L,
	const bigint & Delta,
	const bigint & h,
	int info)
{
	long w;
	bigint N;
	xbigfloat l, z, d, G, H;
	timer t;

	//
	// If L is non-zero, than
	//
	// |L/exp(l(n,Delta))-1| < 2^{-7}
	// C(n) <= 1/4
	//
	// otherwise it is computed and returned.
	//
	if (L.is_zero())
		approximate_exp_l_n_delta(L, Delta, info);

	//
	// Relative 6-approximation H to sqrt(2) * exp(l(n,Delta))
	//

	multiply(H, L, 1.41421);
	H.truncate(9);

	if (info > 2) {
		std::cout << "relative 6-approx. H to sqrt(2) * exp(l(n, Delta))";
		std::cout << ", H = " << H << std::endl;
	}

	//
	// w = number of roots of unity
	//
	if (Delta == -3)
		w = 6;
	else if (Delta == -4)
		w = 4;
	else
		w = 2;

	//
	// relative 6-approximation G to h * kappa(Delta)
	//
	// G = divide(truncate(h, 12), (w/2) * d, 11);
	// G = truncate(G * Pi(7), 9);
	sqrt(d, -Delta, 11);
	d.divide_by_2();
	multiply(G, w, d);
	d = h;
	d.truncate(12);
	divide (G, d, G, 11);
	G *= 3.14159;
	G.truncate(9);
	t.stop_timer();

	if (info > 2) {
		std::cout << "relative 6-approximation G to h * kappa(Delta) ";
		std::cout << "G = " << G << std::endl;
		if (info > 2) {
			std::cout << "(Timing for H and G: ";
			t.print();
			std::cout << ")" << std::endl;
		}
	}

	//
	// Check  h * kappa(Delta) < sqrt(2) * exp(l(n,Delta))
	//
	if (G < H)
		return true;
	else
		return false;
}


//
// real case: verify_hR
//
//  Delta is discriminant
//  r is relative 7-approx. to integer multiple R of regulator R(\D)
//  h is integer multiple of class number h(\D)
//
//  returns true iff h*R == h(\D)*R(\D)
//

 // for convenience
bool quadratic_order
::verify_hR(
	const bigint & Delta,
	const bigint & h,
	const xbigfloat & r,
	int info)
{
	xbigfloat L;
	return (verify_hR(L, Delta, h, r, info));
}


bool quadratic_order
::verify_hR(xbigfloat & L,
	    const bigint & Delta,
	    const bigint & h,
	    const xbigfloat & r,
	    int info)
{
	xbigfloat d, H, G;

	//
	// If L is non-zero, than
	//
	// |L/exp(l(n,Delta))-1| < 2^{-7}
	// C(n) <= 1/4
	//
	// otherwise it is computed and returned.
	//
	if (L.is_zero())
		approximate_exp_l_n_delta(L, Delta, info);


	//
	// Relative 6-approximation H to sqrt(2) * exp(l(n,Delta))
	// C(n) <= 1/4
	//
	multiply(H, L, 1.41421);
	H.truncate(9);

	if (info > 2) {
		std::cout << "relative 6-approx. H to sqrt(2) * exp(l(n, Delta))";
		std::cout << ", H = " << H << std::endl;
	}

	//
	// relative 6-approximation G to h * kappa(Delta)
	//

	// G = divide(truncate(h, 12), d/2, 11);
	sqrt(d, Delta, 11);
	d.divide_by_2();
	G = h;
	G.truncate(12);
	divide(G, G, d, 11);

	// G = truncate(G * r, 9);
	G *= r;
	G.truncate(9);

	if (info > 2) {
		std::cout << "relative 6-approximation G to h * kappa(Delta) ";
		std::cout << "G = " << G << std::endl;
	}

	//
	// Check h * kappa(Delta) < sqrt(2) * exp(l(n,Delta))
	//
	if (G < H)
		return true;
	else
		return false;
}




//
//
//  Lower bound on the regulator
//
//

 //
 // lower_regulator_bound
 //
 // returns m with 2^{m} <= R(Delta)
 //

long quadratic_order
::lower_regulator_bound(const bigint & Delta)
{
	if (Delta == 5)
		return -2;

	if (Delta == 8)
		return -1;

	if (Delta <= 29)
		return 0;

	//
	// We use RD >= log(1/2 (\sqrt{D-4} + \sqrt{D}))
	//		 >= log(\sqrt{D-4})
	//
	// (see Jacobson, Lukes, Williams, An investigation of Bounds
	//   for the Regulator of Quadratic Fields, Exp. Math., Vol. 4,
	//   1995, No. 3)
	//
	//  Determine d = \sqrt{D-4} (1+e), |e| < 2^{-1-b(\sqrt{D-4})}.
	//  This implies \sqrt{D-4} > d -1/2.
	//
	//  Determine |l - log(d-1/2)| < 1/2.
	//  This implies l -1/2 = |l| -1/2 < log(d-1/2) < log(\sqrt{D-4}).
	//
	//  Note Delta >= 29 implies that all values are positive.
	//

	xbigfloat d, l;

	//
	// |\sqrt{D-4}| < \sqrt{D} < 2^{b(D)/2 + 1}
	//

	sqrt(d, Delta-4, b_value(Delta)/2 + 2);
	d -= xbigfloat(0.5);

	//
	// l >= 2^{b(l)-1}
	//

	log(l, d, 1);
	l -= xbigfloat(0.5);

	return (l.b_value()-1);
}



//
// integer multiple of class number known
//

// for convenience
long quadratic_order
::lower_regulator_bound(
	const bigint & Delta,
	const bigint & h,
	int info)
{
	xbigfloat L;
	return (lower_regulator_bound(L, Delta, h, info));
}



long quadratic_order
::lower_regulator_bound(xbigfloat & L,
			const bigint & Delta,
			const bigint & h,
			int info)
{
	xbigfloat l, d;

	//
	// If L is non-zero, than
	//
	// |L/exp(l(n,Delta))-1| < 2^{-7}
	// C(n) <= 1/4
	//
	// otherwise it is computed and returned.
	//
	if (L.is_zero())
		approximate_exp_l_n_delta(L, Delta, info);

	sqrt(d, Delta, 2);

	multiply(d, d, L);
	divide (l, d, h, 2);

	//
	// Now we have L = L(1,chi_D) (1+e3), |e3| < 1,
	//	           d = sqrt{D}    (1+e4), |e4| < 2^{-2}.
	//
	// We use RD = L(1,chi_D) \sqrt{D} / (2 hD)
	//           > L/2 * d/(1+1/4) / (2h)
	//		 > l / 8
	//  to determine the lower bound for RD.
	//

	return(l.b_value() - 4);
}



//
//
// Approximations to Ln of power product
// of quadratic numbers
//
//

//
// Determines an absolute k-approximation l to
//
//   L = Ln \prod_{i=0}^{n-1} q[i]^M[i][j]
//
// where 0 <= j < M.get_no_of_columns and n = M.get_no_of_rows()
// = q.get_size()
//
// If an absolute approximation log_q[i] is not accurate enough
// for computing l, a more accurate approximation is
// computed and stored in log_q[i] and the precision
// is stored in prec_log_q[i]; otherwise
// log_q[i] is truncated to the precision necessary for
// computing l.
//

void quadratic_order
::absolute_Ln_approximation (
	xbigfloat & l,
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	const matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	lidia_size_t j,
	long k,
	int info)
{
	lidia_size_t n, i;
	bool no_precomp_yet;
	xbigfloat h;
	long t, kbn2;

	n = q.get_size();
	no_precomp_yet = false;
	l.assign_zero();

	if (info > 2) {
		std::cout << "quadratic_order3.c::";
		std::cout << "absolute_Ln_approximation (...matrix...)";
		std::cout << std::endl;
	}


	//
	// Verify preconditions
	//

	if (log_q.get_size() == 0) {
		log_q.set_capacity(n);
		prec_log_q.set_capacity(n);
		no_precomp_yet = true;
	}

	//
	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(n)+2
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(n))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	//
	// Determine the approximation
	//

	for (i = 0; i < n; i++) {
		if (info > 2)
			std::cout << "row " << i << " : ";

		if (!M.member(i, j).is_zero()) {
			// absolute precision t for Ln(q[i])
			// t = kbn2 + b_value(M[i][j]);

			if (check_overflow(t, kbn2, b_value(M.member(i, j)))) {
				lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
						     "precision overflow.");
				return;
			}

			if (info > 2) {
				std::cout << " need precision " << t;
				std::cout << " (" << prec_log_q[i] << ") for Ln(q[" << i << "])." << std::endl;
			}

			// approximate M[i][j] * Ln(q[i])
			if (no_precomp_yet || t > prec_log_q[i]) {
				if (t < -1) t = -1;
				q[i].get_absolute_Ln_approximation_old(log_q[i], t);
				prec_log_q[i] = t;
				multiply (h, log_q[i], M.member(i, j));
			}
			else if (t == prec_log_q[i]) {
				multiply (h, log_q[i], M.member(i, j));
			}
			else {
				// t += 1+b_value(log_q[i])

				if (check_overflow(t, t, 1)) {
					lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
							     "precision overflow.");
					return;
				}

				if (check_overflow(t, t, log_q[i].b_value())) {
					lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
							     "precision overflow.");
					return;
				}

				truncate (h, log_q[i], t);
				multiply (h, h, M.member(i, j));
			}

			// truncate and add h to the sum;
			// t = kbn2 + b_value(h)

			if (check_overflow(t, kbn2, b_value(h))) {
				lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
						     "precision overflow.");
				return;
			}

			truncate (h, h, t);
			add (l, l, h);
		}
		else {
			if (info > 2)
				std::cout << "zero exponent." << std::endl;

			if (no_precomp_yet) {
				t = 10;
				q[i].get_absolute_Ln_approximation_old(log_q[i], t);
				prec_log_q[i] = t;
			}
		}
	}

	// truncate l
	// t = k + 1 + b_value(l)

	if (check_overflow(t, k+1, b_value(l))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}
	truncate(l, l, t);
}



//
// For debugging only
//

// log |alpha| with bigfloats

void quadratic_order
::absolute_log_approximation (
	bigfloat & l,
	const matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	lidia_size_t j,
	int info)
{
	bigint 			max_exp;
	bigfloat 			h, delta;
	long 		 	old_prec, new_prec;
	base_vector< bigfloat > log_q;
	lidia_size_t			i, n;


	// initialize

	n = q.get_size();
	l.assign_zero();

	if (n == 0)
		return;


	// determine largest exponent
	max_exp = abs(M.member(0, j));

	for (i = 1; i < n; i++)
		if (abs(M.member(i, j)) > max_exp)
			max_exp = abs(M.member(i, j));

	std::cout << "Maximal exponent is " << max_exp << std::endl;


	// set bigfloat precision
	old_prec = bigfloat::get_precision();

	new_prec = 2*(b_value(q[0].get_discriminant()) + b_value(n) + 2 * b_value(max_exp));

	std::cout << "Setting bigfloat precision to " << new_prec << std::endl;
	bigfloat::set_precision(new_prec);


	// precompute log_q

	delta = sqrt(bigfloat(q[0].get_discriminant()));
	log_q.set_capacity(n);

	for (i = 0; i < n; i++) {
		multiply (h, q[i].get_b(), delta);
		add      (h, h, q[i].get_a());
		divide   (h, h, q[i].get_d());
		h.absolute_value();
		log (log_q[i], h);
	}


	// approximate log|alpha|

	for (i = 0; i < n; i++) {
		if (!M.member(i, j).is_zero()) {
			multiply (h, log_q[i], M.member(i, j));
			add (l, l, h);
		}
	}

	// restore bigfloat precision
	//bigfloat::set_precision(old_prec);
}



//
// For debugging only
//

// Ln(alpha) with bigfloats

void quadratic_order
::absolute_Ln_approximation (
	bigfloat & l,
	const matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	lidia_size_t j,
	int info)
{
	bigint 			max_exp;
	bigfloat 			h, delta;
	long 		 	old_prec, new_prec;
	base_vector< bigfloat > log_q;
	lidia_size_t			i, n;


	// initialize

	n = q.get_size();
	l.assign_zero();

	if (n == 0)
		return;


	// determine largest exponent
	max_exp = abs(M.member(0, j));

	for (i = 1; i < n; i++)
		if (abs(M.member(i, j)) > max_exp)
			max_exp = abs(M.member(i, j));

	std::cout << "Maximal exponent is " << max_exp << std::endl;


	// set bigfloat precision
	old_prec = bigfloat::get_precision();

	new_prec = 2*(b_value(q[0].get_discriminant()) + b_value(n) + 2 * b_value(max_exp));

	std::cout << "Setting bigfloat precision to " << new_prec << std::endl;
	bigfloat::set_precision(new_prec);


	// precompute Ln(q)

	delta = sqrt(bigfloat(q[0].get_discriminant()));
	log_q.set_capacity(n);

	for (i = 0; i < n; i++) {
		multiply (h, q[i].get_b(), delta);
		add      (h, h, q[i].get_a());
		divide   (h, h, q[i].get_d());
		h.absolute_value();
		log (log_q[i], h);

		h = bigfloat(norm(q[i]));
		h.absolute_value();
		log (h, h);
		h.divide_by_2();

		subtract(log_q[i], log_q[i], h);
	}


	// approximate Ln(alpha)

	for (i = 0; i < n; i++) {
		if (!M.member(i, j).is_zero()) {
			multiply (h, log_q[i], M.member(i, j));
			add (l, l, h);
		}
	}

	// restore bigfloat precision
	//bigfloat::set_precision(old_prec);
}




// for debugging only

void quadratic_order
::absolute_Ln_approximation_with_check (
	xbigfloat & l,
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	const matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	lidia_size_t j,
	long k,
	int info)
{

	//
	// bigfloat precomputation
	//

	bigint 			max_exp;
	bigfloat 			l_big, h_big, delta;
	long 		 	old_prec, new_prec;
	base_vector< bigfloat > log_q_big;
	lidia_size_t			i, n;


	// initialize

	n = q.get_size();
	l_big.assign_zero();

	if (n == 0)
		return;


	// determine largest exponent
	max_exp = abs(M.member(0, j));

	for (i = 1; i < n; i++)
		if (abs(M.member(i, j)) > max_exp)
			max_exp = abs(M.member(i, j));

	std::cout << "Maximal exponent is " << max_exp << std::endl;


	// set bigfloat precision
	old_prec = bigfloat::get_precision();

	new_prec = 2*(b_value(q[0].get_discriminant()) + b_value(n) + 2 * b_value(max_exp));

	std::cout << "Setting bigfloat precision to " << new_prec << std::endl;
	bigfloat::set_precision(new_prec);


	// precompute Ln(q)

	delta = sqrt(bigfloat(q[0].get_discriminant()));
	log_q_big.set_capacity(n);

	for (i = 0; i < n; i++) {
		multiply (h_big, q[i].get_b(), delta);
		add      (h_big, h_big, q[i].get_a());
		divide   (h_big, h_big, q[i].get_d());
		h_big.absolute_value();
		log (log_q_big[i], h_big);

		h_big = bigfloat(norm(q[i]));
		h_big.absolute_value();
		log (h_big, h_big);
		h_big.divide_by_2();

		subtract(log_q_big[i], log_q_big[i], h_big);
	}


	//
	// now the xbigfloat computation
	//

	bool no_precomp_yet;
	xbigfloat h;
	long t, tt, kbn2;

	n = q.get_size();
	no_precomp_yet = false;
	l.assign_zero();

	std::cout << "quadratic_order3.c::absolute_Ln_approximation_with_check (...matrix...)";
	std::cout << std::endl;
	std::cout << "number of terms n = " << n << std::endl;
	std::cout << "b(n) = " << b_value(n) << std::endl;
	std::cout << "Looking for absolute k = " << k << " approximation" << std::endl;

	//
	// Verify preconditions
	//

	std::cout << "log_q.get_size() == " << log_q.get_size() << std::endl;

	if (log_q.get_size() == 0) {
		log_q.set_capacity(n);
		prec_log_q.set_capacity(n);
		no_precomp_yet = true;
	}

	if (no_precomp_yet)
		std::cout << "no precomputation yet" << std::endl;
	else
		std::cout << "precomputed logs exist" << std::endl;

	//
	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(n)+2
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(n))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	//
	// Determine the approximation
	//

	for (i = 0; i < n; i++) {
		std::cout << "row " << i << std::endl;

		if (!M.member(i, j).is_zero()) {
			// absolute precision t for Ln(q[i])
			// t = kbn2 + b_value(M[i][j]);

			if (check_overflow(t, kbn2, b_value(M.member(i, j)))) {
				lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
						     "precision overflow.");
				return;
			}

			// approximate M[i][j] * Ln(q[i])
			if (no_precomp_yet || t > prec_log_q[i]) {
				if (t < -1) t = -1;
				q[i].get_absolute_Ln_approximation_with_check(log_q[i], t);
				prec_log_q[i] = t;

				// check Ln(q[i])
				if (!(b_value(log_q[i]-log_q_big[i]) <= -t)) {
					std::cout << "ERROR: report follows:" << std::endl;
					std::cout << "if (no_precomp_yet || t > prec_log_q[i])" << std::endl;
					std::cout << "absolute Ln precision t = " << t << " for " << std::endl;
					std::cout << "relative generator with index i = " << i << std::endl;
					std::cout << "q[i] = " << q[i] << std::endl;
					std::cout << "prec_log_q[i] = " << prec_log_q[i] << std::endl;
					std::cout << "Ln(q[i]) as xbigfloat = " << log_q[i] << std::endl;
					std::cout << "Ln(q[i]) as bigfloat = " << log_q_big[i] << std::endl;
					std::cout << "M.member(i, j) = " << M.member(i, j) << std::endl;
					std::cout << "b(M.member(i, j)) = " << b_value(M.member(i, j)) << std::endl;
					exit(1);
				}

				multiply (h, log_q[i], M.member(i, j));
				multiply (h_big, log_q_big[i], M.member(i, j));
			}
			else if (t == prec_log_q[i]) {
				// check Ln(q[i])
				if (!(b_value(log_q[i]-log_q_big[i]) <= -t)) {
					std::cout << "ERROR: report follows:" << std::endl;
					std::cout << "else if (t == prec_log_q[i])" << std::endl;
					std::cout << "absolute Ln precision t = " << t << " for " << std::endl;
					std::cout << "relative generator with index i = " << i << std::endl;
					std::cout << "q[i] = " << q[i] << std::endl;
					std::cout << "prec_log_q[i] = " << prec_log_q[i] << std::endl;
					std::cout << "Ln(q[i]) as xbigfloat = " << log_q[i] << std::endl;
					std::cout << "Ln(q[i]) as bigfloat = " << log_q_big[i] << std::endl;
					std::cout << "M.member(i, j) = " << M.member(i, j) << std::endl;
					std::cout << "b(M.member(i, j)) = " << b_value(M.member(i, j)) << std::endl;
					exit(1);
				}

				multiply (h, log_q[i], M.member(i, j));
				multiply (h_big, log_q_big[i], M.member(i, j));
			}
			else {
				// t += 1+b_value(log_q[i])

				if (check_overflow(tt, t, 1)) {
					lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
							     "precision overflow.");
					return;
				}

				if (check_overflow(tt, tt, log_q[i].b_value())) {
					lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
							     "precision overflow.");
					return;
				}

				truncate (h, log_q[i], tt);

				// check Ln(q[i])
				if (!(b_value(h-log_q_big[i]) <= -t)) {
					std::cout << "else (t < prec_log_q[i])" << std::endl;
					std::cout << "ERROR: report follows:" << std::endl;
					std::cout << "absolute Ln precision t = " << t << " for " << std::endl;
					std::cout << "relative generator with index i = " << i << std::endl;
					std::cout << "q[i] = " << q[i] << std::endl;
					std::cout << "prec_log_q[i] = " << prec_log_q[i] << std::endl;
					std::cout << "Ln(q[i]) as xbigfloat = " << h << std::endl;
					std::cout << "Ln(q[i]) as bigfloat = " << log_q_big[i] << std::endl;
					std::cout << "M.member(i, j) = " << M.member(i, j) << std::endl;
					std::cout << "b(M.member(i, j)) = " << b_value(M.member(i, j)) << std::endl;
					exit(1);
				}

				multiply (h, h, M.member(i, j));
				multiply (h_big, log_q_big[i], M.member(i, j));
			}

			// truncate and add h to the sum;
			// t = kbn2 + b_value(h)

			if (check_overflow(t, kbn2, b_value(h))) {
				lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
						     "precision overflow.");
				return;
			}

			if (info > 1)
				std::cout << "h before truncating (" << t << ") = " << h << std::endl;

			truncate (h, h, t);
			add (l, l, h);
			add (l_big, l_big, h_big);

			if (info > 1) {
				std::cout << "added h = " << h << std::endl;
				std::cout << "new l = " << l << std::endl;
			}
		}
		else {
			if (no_precomp_yet) {
				t = 10;
				q[i].get_absolute_Ln_approximation_old(log_q[i], t);
				prec_log_q[i] = t;
			}

			if (info > 1)
				std::cout << "Skipping zero exponent." << std::endl;
		}
	}

	// truncate l
	// t = k + 1 + b_value(l)

	if (check_overflow(t, k+1, b_value(l))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (info > 1)
		std::cout << "final l before truncating = " << l << std::endl;

	truncate(l, l, t);

	// check l
	if (!(b_value(l-l_big) <= -k)) {
		std::cout << "ERROR: report follows:" << std::endl;
		std::cout << "absolute Ln precision k = " << k;
		std::cout << "l as xbigfloat = " << l << std::endl;
		std::cout << "l as bigfloat = " << l_big << std::endl;
		exit(1);
	}
}





void quadratic_order
::absolute_Ln_approximation (
	xbigfloat & l,
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	const base_vector< bigint > & exponents,
	const base_vector< quadratic_number_standard > & q,
	long k,
	int info)
{
	lidia_size_t n, i;
	bool no_precomp_yet;
	xbigfloat h;
	long t, kbn2;

	n = q.get_size();
	no_precomp_yet = false;
	l.assign_zero();

	//
	// Verify preconditions
	//

	if (log_q.get_size() == 0) {
		log_q.set_capacity(n);
		prec_log_q.set_capacity(n);
		no_precomp_yet = true;
	}

	//
	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(n)+2
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(n))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	//
	// Determine the approximation
	//

	for (i = 0; i < n; i++) {
		if (!exponents[i].is_zero()) {
			// absolute precision t for Ln(q[i])
			// t = kbn2 + b_value(exponents[i]);

			if (check_overflow(t, kbn2, b_value(exponents[i]))) {
				lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
						     "precision overflow.");
				return;
			}

			// approximate exponents[i] * Ln(q[i])
			if (no_precomp_yet || t > prec_log_q[i]) {
				if (t < -1) t = -1;
				q[i].get_absolute_Ln_approximation_old(log_q[i], t);
				prec_log_q[i] = t;
				multiply (h, log_q[i], exponents[i]);

				if (info > 1) {
					std::cout << "absolute_Ln_approximation(pp< qn > )::" << std::endl;
					std::cout << "(Re-)computing abs. " << t << "approx. ";
					std::cout << "to q(" << i << ") = " << log_q[i] << std::endl;
				}
			}
			else if (t == prec_log_q[i]) {
				multiply (h, log_q[i], exponents[i]);
			}
			else {
				// t += 1+b_value(log_q[i])

				if (check_overflow(t, t, 1)) {
					lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
							     "precision overflow.");
					return;
				}

				if (check_overflow(t, t, log_q[i].b_value())) {
					lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
							     "precision overflow.");
					return;
				}

				truncate (h, log_q[i], t);
				multiply (h, h, exponents[i]);
			}

			// truncate and add h to the sum
			// t = kbn2 + b_value(h)

			if (check_overflow(t, kbn2, b_value(h))) {
				lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
						     "precision overflow.");
				return;
			}

			if (info > 1)
				std::cout << "h before truncating = " << h << std::endl;

			truncate (h, h, t);
			add (l, l, h);

			if (info > 1) {
				std::cout << "added h = " << h << std::endl;
				std::cout << "new l = " << l << std::endl;
			}
		}
		else {
			if (no_precomp_yet) {
				t = 10;
				q[i].get_absolute_Ln_approximation_old(log_q[i], t);
				prec_log_q[i] = t;
			}

			if (info > 1)
				std::cout << "Skipping zero exponent." << std::endl;
		}
	}

	// truncate l
	// t = k + 1 + b_value(l)

	if (check_overflow(t, k+1, b_value(l))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (info > 1)
		std::cout << "final l before truncating = " << l << std::endl;

	truncate(l, l, t);
}





void quadratic_order
::absolute_Ln_approximation (
	xbigfloat & l,
	const base_vector< quadratic_number_standard > & q,
	const base_vector< bigint > & e,
	long k,
	int info)
{
	lidia_size_t n, i;
	xbigfloat h;
	long t, kbn2;

	(void)info;
	n = q.get_size();
	l.assign_zero();

	//
	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(n)+2
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(n))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	//
	// Determine the approximation
	//

	for (i = 0; i < n; i++) {
		// absolute precision t for Ln(q[i])
		// t = kbn2 + b_value(e[i]);

		if (check_overflow(t, kbn2, b_value(e[i]))) {
			lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
					     "precision overflow.");
			return;
		}

		// approximate e[i] * Ln(q[i])

		if (t < -1) t = -1;
		q[i].get_absolute_Ln_approximation_old(h, t);
		multiply (h, h, e[i]);

		// truncate and add h to the sum
		// t = kbn2 + b_value(h)

		if (check_overflow(t, kbn2, b_value(h))) {
			lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
					     "precision overflow.");
			return;
		}

		truncate (h, h, t);
		add (l, l, h);
	}

	// truncate l
	// t = k + 1 + b_value(l)

	if (check_overflow(t, k+1, b_value(l))) {
		lidia_error_handler ("absolute_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	truncate(l, l, t);
}


void quadratic_order
::relative_Ln_approximation (
	xbigfloat & l,
	const base_vector< quadratic_number_standard > & q,
	const base_vector< bigint > & e,
	long k,
	long m,
	int info)
{
	long t;

	// t = k+1-m
	if (check_overflow(t, k, 1)) {
		lidia_error_handler ("relative_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}
	if (check_overflow(t, t, -m)) {
		lidia_error_handler ("relative_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	absolute_Ln_approximation (l, q, e, t, info);

	// t = k+1+b_value(l)
	k++;
	if (check_overflow(t, k, b_value(l))) {
		lidia_error_handler ("relative_Ln_approximation(pp< qn > )",
				     "precision overflow.");
		return;
	}

	l.truncate (t);
}





//
//
//  Detecting and removing +-1 units.
//
//

//
// Removes all columns j from M with
//
//  |alpha_j| = \prod_{i=0}^{n-1} q[i]^M[i][j] = 1
//
// 0 <= j < M.get_no_of_columns() by computing an absolute
// -m+1 approximation l to Ln(alpha_j) and verifying
// |l| < 2^{m-1}. Here, 2^{m} is a lower bound for the regulator
// and n = M.get_no_of_rows() = q.get_size().
//

//
// Condition: All alpha_j must be positive units,
// 	      0 <= j < M.nof_colums().
//

void quadratic_order
::remove_one_from_unit_array (
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	long m,
	int info)
{
	lidia_size_t c, i;
	lidia_size_t nzi; // non-zero index
	xbigfloat l;
	long mmp1; // minus m plus 1
	long mm1; // m minus 1

	c = M.get_no_of_columns();
	nzi = 0;

	if (check_overflow(mmp1, -m, 1)) {
		lidia_error_handler ("remove_one_from_unit_array",
				     "precision overflow.");
		return;
	}

	mm1 = -mmp1;

	if (info > 1) {
		std::cout << "-m+1 = " << mmp1 << std::endl;
		std::cout << "m-1 = " << mm1 << std::endl;
	}

	for (i = 0; i < c; i++) {
		if (info > 1)
			std::cout << "verifying column " << i << ": " << std::flush;

		// absolute -m+1 approx. to Ln(alpha_i)
#ifndef DEBUG_WITH_CHECK
		absolute_Ln_approximation (
			l, log_q, prec_log_q, M, q, i, mmp1);
#else
		absolute_Ln_approximation_with_check (
			l, log_q, prec_log_q, M, q, i, mmp1);
#endif

		// alpha_i == +-1 iff |l| < 2^{m-1} iff (l == 0 || b(l) <= m-1)
		if (l.is_zero() || b_value(l) <= mm1) {
			if (info > 1) {
				std::cout << "quadratic_order3.c::remove_one_from_unit_array::" << std::endl;
				std::cout << "+-1 found (column " << i << ")" << std::endl;
			}
		}
		else {
			// store column i in column nzi
			if (i > nzi)
				M.swap_columns(nzi, i);

			nzi++;
		}
	}

	if (nzi > 0)
		M.set_no_of_columns(nzi);
	else {
		c = M.get_no_of_rows();
		M.set_no_of_columns(1);
		for (i = 0; i < c; i++)
			M.sto(i, 0, 0);
	}
}



//
// Returns true, if alpha = prod_{i=0}^{n-1} q[i]^tmp_exponents[i] == +-1.
// Otherwise, returns false.
//
// alpha must be a unit.
//
// 2^m <= regulator.
//

bool quadratic_order
::is_one(base_vector< xbigfloat > & log_q,
	 base_vector< long > & prec_log_q,
	 const base_vector< quadratic_number_standard > & q,
	 const base_vector< bigint > & tmp_exponents,
	 lidia_size_t m,
	 int info)
{
	xbigfloat l;
	long mmp1; // minus m plus 1
	long mm1; // m minus 1

	if (check_overflow(mmp1, -m, 1)) {
		lidia_error_handler ("quadratic_order3.c::is_one",
				     "precision overflow.");
		return false;
	}

	mm1 = -mmp1;

	if (info > 1) {
		std::cout << "-m+1 = " << mmp1 << std::endl;
		std::cout << "m-1 = " << mm1 << std::endl;
	}


	// absolute -m+1 approx. to Ln(alpha)
	absolute_Ln_approximation (l, log_q, prec_log_q, tmp_exponents, q, mmp1, info);

	// alpha == +-1 iff |l| < 2^{m-1} iff (l == 0 || b(l) <= m-1)
	if (l.is_zero() || b_value(l) <= mm1) {
		if (info > 1)
			std::cout << "is equal to one." << std::endl;

		return true;
	}
	else {
		if (info > 1)
			std::cout << "is not equal to one." << std::endl;

		return false;
	}
}




//
//  Verifies whether
//
//   u = (alpha_i)^x * (alpha_j)^y
//
//  generates the subgroup generated by (alpha_i,alpha_j)
//  by verifying that
//
//  u^M1 / alpha_i == 1 and u^M2 / alpha_j == 1.
//

int quadratic_order
::is_generator (
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	const matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	lidia_size_t i,
	lidia_size_t j,
	const bigint & x,
	const bigint & y,
	const bigint & M1,
	const bigint & M2,
	long m,
	int info)
{
	base_vector< bigint > u_exponents;
	base_vector< bigint > tmp_exponents;
	lidia_size_t n, k;
	bigint tmp, tmp2;

	bigfloat	s_bigfloat, t_bigfloat;


	n = M.get_no_of_rows();
	u_exponents.set_capacity(n);
	tmp_exponents.set_capacity(n);


	// determine the exponents for u
	for (k = 0; k < n; ++k) {
		multiply(tmp, M.member(k, i), x);
		multiply(tmp2, M.member(k, j), y);
		add(u_exponents[k], tmp, tmp2);
	}


	// determine the exponents for u^M1 / alpha_i

	for (k = 0; k < n; ++k) {
		multiply(tmp_exponents[k], u_exponents[k], M1);
		subtract(tmp_exponents[k], tmp_exponents[k], M.member(k, i));
	}

	// check for 1

	if (!is_one(log_q, prec_log_q, q, tmp_exponents, m)) {
		return 1;
	}


	// determine the exponents for u^M2 / alpha_j

	for (k = 0; k < n; ++k) {
		multiply(tmp_exponents[k], u_exponents[k], M2);
		subtract(tmp_exponents[k], tmp_exponents[k], M.member(k, j));
	}

	// check for 1

	if (!is_one(log_q, prec_log_q, q, tmp_exponents, m)) {
		return 2;
	}

	return 0;
}



//
//
// Real gcd algorithm for two power products of units
//
//
//


//
// Function: cfrac_expansion_bounded_den
//
//  Computes the last convergent conv_num/conv_den
//  of the continued fraction expansion of num/den,
//  of which absolute value of the denominator is
//  less or equal to S, where S is positive.
//

void quadratic_order
::cfrac_expansion_bounded_den (
	bigint & conv_num,
	bigint & conv_den,
	bigint num,
	bigint den,
	const bigint & S,
	int info)
{

	bigint 	 r; // contains the rest
	// in the euclidian step

	bigint 	 q; // contains the quotient
	// in the euclidean step
	int 	 sign;

	bigint 	prev1_conv_num, prev1_conv_den;
	bigint 	prev2_conv_num, prev2_conv_den;

	//
	// sign * (num/den) = n/d,
	// num, den are positive
	//

	sign = 1;

	if (num.is_negative()) {
		num.negate();
		sign *= -1;
	}

	if (den.is_negative()) {
		den.negate();
		sign *= -1;
	}


	// compute the partial quotients
	// and the convergents

	// first convergent: q[0]/1

	if (S.is_negative())
		lidia_error_handler("quadratic_order3.c::cfrac_expansion_bounded_den",
				    "Negative upper bound.");

	if (S.is_zero())
		lidia_error_handler("quadratic_order3.c::cfrac_expansion_bounded_den",
				    "Upper bound is zero.");


	// first quotient q
	// first convergent: conv_num / d_n
	// first convergent: q / 1

	div_rem (q, r, num, den);
	num = den;
	den = r;

	conv_num = q;
	conv_den.assign_one();

	if (info > 1) {
		std::cout << "cfrac::first convergent ";
		std::cout << conv_num << " / " << conv_den << std::endl;
	}

	if (!den.is_zero()) {
		// store previous convergent

		prev1_conv_num = conv_num;
		prev1_conv_den = conv_den;

		conv_num = q;

		// second quotient q

		div_rem (q, r, num, den);
		num = den;
		den = r;

		// second convergent conv_num / conv_den
		// second convergent: (q[0]q[1]+1)/q[1]

		multiply (conv_num, conv_num, q);
		inc (conv_num);

		conv_den = q;

		if (info > 1) {
			std::cout << "cfrac::second convergent ";
			std::cout << conv_num << " / " << conv_den << std::endl;
		}

		// stop ?

		if (conv_den > S) {
			conv_num = prev1_conv_num;
			conv_den = prev1_conv_den;
		}
		else if (!den.is_zero()) {
			do {
				prev2_conv_num = prev1_conv_num;
				prev2_conv_den = prev1_conv_den;

				prev1_conv_num = conv_num;
				prev1_conv_den = conv_den;

				// next quotient q

				div_rem (q, r, num, den);
				num = den;
				den = r;

				// next convergent conv_num / conv_den
				// conv_num(k) = q(k) * conv_num(k-1) + conv_num(k-2)
				// conv_den(k) = q(k) * conv_den(k-1) + conv_den(k-2)

				multiply (conv_num, q, prev1_conv_num);
				multiply (conv_den, q, prev1_conv_den);

				add (conv_num, conv_num, prev2_conv_num);
				add (conv_den, conv_den, prev2_conv_den);

				if (info > 1) {
					std::cout << "cfrac::convergent ";
					std::cout << conv_num << " / " << conv_den << std::endl;
				}
			} while (!den.is_zero() && conv_den <= S);

			if (conv_den > S) {
				conv_num = prev1_conv_num;
				conv_den = prev1_conv_den;
			}
		}
	}

	if (sign == -1)
		conv_num.negate();
}


//
// Function: rgcd
//
// Determines a pair (x,y) with
//
// x Ln alpha_j1 + y Ln alpha_j2 = rgcd(Ln alpha_j1, Ln alpha_j2)
//
// where alpha_jk = \prod_{i=0}_{n-1} q[i]^M[i][jk],  k = 1,2
//
// and n = M.get_no_of_rows() = q.get_size()
// and 0 \leq j1 <= j2 < M.get_no_of_columns().
//
// 2^{m} is a lower regulator bound.
//
// for log_q, prec_log_q see above.
//
// Condition: alpha_j1, alpha_j2 != 1
//

void quadratic_order
::rgcd (
	bigint & x,
	bigint & y,
	bigint & M1,
	bigint & M2,
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	const matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	lidia_size_t j1,
	lidia_size_t j2,
	long m,
	int info)
{
	xbigfloat s, t, h;
	bigint    S, num, den, gcdM1M2;
	bigint    big_k, big_ks, big_kt;
	bigint    bs, bt;
	long      k, ks, kt;
	bool      s_is_negative, t_is_negative;

#ifdef DEBUG
	bigfloat  s_bigfloat, t_bigfloat;
	xbigfloat diff_xbigfloat;
	int 	     M_index;
	if (info < 2) info = 2;
#endif

	//
	// s = relative 1-approximation to Ln(alpha_j1)
	// t = relative 1-approximation to Ln(alpha_j2)
	// |b(s) - b(Ln(alpha_j1))| <= 1
	// |b(t) - b(Ln(alpha_j2))| <= 1
	//

	big_k = -m;
	big_k++;

	if (big_k.longify(k)) {
		lidia_error_handler ("rgcd(pp< qn > )",
				     "precision overflow.");
		return;
	}
	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "COMPUTING absolute " << k << " approximations";
		std::cout << std::endl;
	}

#ifndef DEBUG_WITH_CHECK
	absolute_Ln_approximation (
		s, log_q, prec_log_q, M, q, j1, k, info);

	absolute_Ln_approximation (
		t, log_q, prec_log_q, M, q, j2, k, info);
#else
	absolute_Ln_approximation_with_check (
		s, log_q, prec_log_q, M, q, j1, k, 0);

	absolute_Ln_approximation_with_check (
		t, log_q, prec_log_q, M, q, j2, k, 0);
#endif

	bs = b_value(s);
	bt = b_value(t);

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "FOUND s = " << s << std::endl;
		std::cout << "FOUND t = " << t << std::endl;
		std::cout << "FOUND b_value(s) = " << bs << std::endl;
		std::cout << "FOUND b_value(t) = " << bt << std::endl;
	}

	//
	// upper bound S on |M2|
	// S = ceil( 2*|t| / 2^m );
	//
	h.assign(t);
	h.absolute_value();
	if (k >= 0)
		shift_left (h, h, k);
	else
		shift_right (h, h, -k);

	ceil(S, h);

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "FOUND upper bound h on |M2|, h = " << h << std::endl;
		std::cout << "S = ceil(h), S = " << S << std::endl;
	}

	//
	// precisions to guarantee |R1/R2 - s/t| < 1/(2S^2)
	// ks = 2*b_value(S)
	// kt = b_value(t)-b_value(s)
	// k = max {ks, kt} +1
	//

	big_ks = 2 * b_value(S);
	big_kt = bt-bs;

	if (big_ks > big_kt)
		big_k = big_ks;
	else
		big_k = big_kt;

	big_k++;

	//
	// s = relative ks approx. to Ln(alpha_j1)
	// t = relative kt approx. to Ln(alpha_j2)
	//
	// ks = b_value(s)-b_value(t)+4+k;
	// kt = b_value(s)-b_value(t)+6+k;
	//

	big_ks = bs - bt + 4 + big_k;
	big_kt = bs - bt + 6 + big_k;

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "COMPUTING relative " << big_ks << " approximation s" << std::endl;
		std::cout << "COMPUTING relative " << big_kt << " approximation t" << std::endl;
	}

	//
	// ks += -b_value(s)+2
	// kt += -b_value(t)+2
	//

	big_ks -= (bs - 2);
	big_kt -= (bt - 2);

	if (big_ks.longify(ks)) {
		lidia_error_handler ("rgcd(pp< qn > )",
				     "precision overflow.");
		return;
	}
	if (big_kt.longify(kt)) {
		lidia_error_handler ("rgcd(pp< qn > )",
				     "precision overflow.");
		return;
	}

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "COMPUTING absolute " << ks << " approximation s" << std::endl;
		std::cout << "COMPUTING absolute " << kt << " approximation t" << std::endl;
	}


#ifndef DEBUG_WITH_CHECK
	absolute_Ln_approximation (
		s, log_q, prec_log_q, M, q, j1, ks, info);

	absolute_Ln_approximation (
		t, log_q, prec_log_q, M, q, j2, kt, info);
#else
	absolute_Ln_approximation_with_check (
		s, log_q, prec_log_q, M, q, j1, ks);

	absolute_Ln_approximation_with_check (
		t, log_q, prec_log_q, M, q, j2, kt);
#endif

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "FOUND s = " << s << std::endl;
		std::cout << "FOUND t = " << t << std::endl;
	}

	//
	// Finding the convergent
	//
	//
	// k =   s.get_exponent() - b_value(s.get_mantissa());
	// k -= (t.get_exponent() - b_value(t.get_mantissa()));
	//

	big_k = s.get_exponent();
	big_k -= b_value(s.get_mantissa());
	big_k -= t.get_exponent();
	big_k += b_value(t.get_mantissa());

	if (big_k.longify(k)) {
		lidia_error_handler ("rgcd(pp< qn > )",
				     "precision overflow.");
		return;
	}

	num = s.get_mantissa();
	den = t.get_mantissa();

	if (k > 0)
		shift_left (num, num, k);
	else
		shift_left (den, den, -k);

	if (num.is_negative()) {
		s_is_negative = true;
		num.negate();
	}
	else
		s_is_negative = false;

	if (den.is_negative()) {
		t_is_negative = true;
		den.negate();
	}
	else
		t_is_negative = false;

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "numerator = " << num << std::endl;
		std::cout << "denominator = " << den << std::endl;
		std::cout << "denominator bound = " << S << std::endl;
	}

	cfrac_expansion_bounded_den (M1, M2, num, den, S, 0);

	if (s_is_negative)
		M1.negate();

	if (t_is_negative)
		M2.negate();

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "Found convergent " << M1 << " / " << M2 << std::endl;
	}

	//
	// compute the x,y with 1 = x M1 + y M2
	//

	gcdM1M2 = xgcd (x, y, M1, M2);

	if (!gcdM1M2.is_one()) {
		lidia_error_handler ("real_gcd(pp< qn > )",
				     "gcd(M1, M2) != 1");
	}

	if (!(x*M1+y*M2 == 1)) {
		lidia_error_handler ("real_gcd(pp< qn > )",
				     "x*M1 + y*M2) != 1");
	}

	if (info > 2) {
		std::cout << "rgcd(pp< qn > as matrix< bigint > )::" << std::endl;
		std::cout << "FOUND x = " << x << std::endl;
		std::cout << "FOUND y = " << y << std::endl;
	}

#ifdef DEBUG
	M_index = is_generator (log_q, prec_log_q, M, q, j1, j2, x, y, M1, M2, m, info);

	if (M_index) {
		std::cout << "ERROR: report follows" << std::endl;

		absolute_Ln_approximation(s_bigfloat, M, q, j1, info);
		absolute_Ln_approximation(t_bigfloat, M, q, j2, info);

		std::cout << "Ln(alpha_j1) as bigfloat = " << s_bigfloat << std::endl;
		std::cout << "Ln(alpha_j2) as bigfloat = " << t_bigfloat << std::endl;

		if (b_value(s-xbigfloat(s_bigfloat)) <= -ks)
			std::cout << "b_value(s-xbigfloat(s_bigfloat)) <= -ks :-)" << std::endl;
		else
			std::cout << "b_value(s-xbigfloat(s_bigfloat)) > -ks :-(" << std::endl;

		if (b_value(t-xbigfloat(t_bigfloat)) <= -kt)
			std::cout << "b_value(t-xbigfloat(t_bigfloat)) <= -kt :-)" << std::endl;
		else
			std::cout << "b_value(t-xbigfloat(t_bigfloat)) > -kt :-(" << std::endl;

		absolute_log_approximation(s_bigfloat, M, q, j1, info);
		absolute_log_approximation(t_bigfloat, M, q, j2, info);

		std::cout << "log|alpha_j1| as bigfloat = " << s_bigfloat << std::endl;
		std::cout << "log|alpha_j2| as bigfloat = " << t_bigfloat << std::endl;

		std::cout << "j1 = " << j1 << std::endl;
		std::cout << "j2 = " << j2 << std::endl;
		std::cout << "M1 = " << M1 << std::endl;
		std::cout << "M2 = " << M2 << std::endl;
		std::cout << "x = " << x << std::endl;
		std::cout << "y = " << y << std::endl;

		if (M_index == 1)
			lidia_error_handler("quadratic_order3.c::rgcd",
					    "u^M1 != alpha_i");
		else if (M_index == 2)
			lidia_error_handler("quadratic_order3.c::rgcd",
					    "u^M2 != alpha_j");
		else
			lidia_error_handler("quadratic_order3.c::rgcd",
					    "unknown M_index");
	}
	else {
		std::cout << "quadratic_order3.c::find_generating_unit::Verified generator" << std::endl;
	}
#endif

}




//
//
//  Determing a generating unit from a given
//  array of units
//
//

//
// Determines a unit
//
//   g = \prod_{j=0}^{n-1} q[i]^x_i
//
// that generates the subgroup of units
// U = <alpha_1,...,alpha_k> (modulo <-1>), where
//
//  alpha_j = \prod_{i=0}_{n-1} q[i]^M[i][j], 0 <= j < k
//  k = M.get_no_of_columns()
//  n = M.get_no_of_rows() = q.get_size()
//  and
//  2^{m} is a lower regulator bound.
//
// Method:
//
//  1) Set all columns of M that form a unit +-1
//     to zero. If all colums are zero, stop.
//     Otherwise
//
//  2) Use a binary tree method, where in each step
//     two units are combined using a rgcd.
//     Manipulate only the columns of M until
//     the first column is the only non-zero one.
//     This column contains x_0,...,x_n-1.
//

void quadratic_order
::find_generating_unit (
	matrix< bigint > & M,
	const base_vector< quadratic_number_standard > & q,
	long m,
	int info)

{
	base_vector< xbigfloat > log_q;
	base_vector< long > prec_log_q;
	bigint used_memory;
	long abs_Ln_q_bound;

#ifdef DEBUG
	if (info < 3) info = 3;
#endif

	if (info > 1) {
		std::cout << "quadratic_order3.c::find_generating_unit::" << std::endl;
		std::cout << M.get_no_of_columns() << " units given." << std::endl;
		//std::cout << "MATRIX" << std::endl;
		//std::cout << M << std::endl;
		//std::cout << "RELATIVE GENERATORS" << std::endl;
		//std::cout << q << std::endl;
	}

	//
	// Verify preconditions
	//

	if (q.get_size() != M.get_no_of_rows()) {
		lidia_error_handler ("find_generating_unit(pp< qn > )",
				     "Number of exponents does not match "
				     "number of quadratic numbers.");
		return;
	}

	//
	// Precompute approximations to logarithms of quadratic numbers.
	//
	// Note: Based on a heuristic upper bound on the absolute values
	//       of the Ln of the quadratic numbers, a heuristic upper bound
	//       on the largest required absolute precision is computed and used
	//       to compute approximations to the Ln's of the quadratic numbers.
	//
	//       Nevertheless, during the computation all sufficient accuracies
	//       are determined and if necessary, the approximations to the Ln's
	//       are recomputed to a higher precision. Hence, the algorithm
	//       really yields the generating unit.
	//

	abs_Ln_q_bound = 10000;
	precompute_logarithm_bases (log_q, prec_log_q, M , q,
				    estimate_pp_accuracy (M, m, abs_Ln_q_bound));

	//
	// remove the roots of unity
	//
	remove_one_from_unit_array (log_q, prec_log_q, M, q, m, info-1);
	if (info > 1) {
		std::cout << "quadratic_order3.c::find_generating_unit::" << std::endl;
		std::cout << M.get_no_of_columns() << " remain after removing +-1." << std::endl;
	}

	if (M.get_no_of_columns() <= 1) {
		if (info > 1) {
			std::cout << "quadratic_order3.c::find_generating_unit::" << std::endl;
			std::cout << M.get_no_of_columns() << " units remain." << std::endl;
			//std::cout << "MATRIX" << std::endl;
			//std::cout << M << std::endl;
		}
		return;
	}
	//
	// binary tree method with application of rgcd
	//
	bigint x, y;
	bigint M1, M2;
	bigint tmp, tmp2;
	lidia_size_t i, j, k;
	lidia_size_t step_size, neighbour;
	lidia_size_t n = M.get_no_of_rows();
	lidia_size_t c = M.get_no_of_columns();
	lidia_size_t iterations = b_value(static_cast<long>(c));

	step_size = 2;
	neighbour = 1;

	for (i = 0; i < iterations; i++) {
		if (info > 1) {
			std::cout << "quadratic_order3.c::find_generating_unit::" << std::endl;
			std::cout << "STEP " << i << std::endl;
		}

		for (j = 0; j+neighbour < c; j += step_size) {
			if (info > 1) {
				std::cout << "combining ";
				std::cout << "(" << j << ", " << j+neighbour << ")" << std::endl;
			}

			// representation of real gcd
			rgcd(x, y, M1, M2, log_q, prec_log_q, M, q, j, j+neighbour, m, info-1);

			if (info > 2) {
				std::cout << "(x, y) = (" << x << ", " << y << ")" << std::endl;
			}

			// build power product by combining columns of M

			for (k = 0; k < n; ++k) {
				multiply(tmp, M.member(k, j), x);
				multiply(tmp2, M.member(k, j+neighbour), y);
				add(tmp, tmp, tmp2);
				M.sto(k, j, tmp);
			}
		}

		step_size <<= 1;
		neighbour <<= 1;
	}

	M.set_no_of_columns(1);

	//
	// compute size of memory used.
	//

	if (info > 1) {
		used_memory = 0;
		for (i = 0; i < q.get_size(); i++)
			add(used_memory, used_memory, bigint(log_q[i].get_mantissa().bit_length())/8);

		std::cout << "quadratic_order3.c::rgcd::" << std::endl;
		std::cout << "number of minima = " << log_q.get_size() << std::endl;
		std::cout << used_memory << " bytes for logarithms." << std::endl;

		used_memory = 0;
		for (i = 0; i < M.get_no_of_rows(); i++)
			add(used_memory, used_memory, bigint(M.member(i, 0).bit_length())/8);

		std::cout << used_memory << " bytes for exponents of fundamental unit." << std::endl;
	}
}



#if 0
void floor_quotient (bigint & Q,
		     const base_vector< quadratic_number_standard > q,
		     const base_vector< bigint > n1,
		     const base_vector< bigint > n2)
{

}
#endif



long quadratic_order
::estimate_pp_accuracy (const matrix< bigint > & M,
			long m,
			long c)

	//
	// estimate accuracy of power product
	//
	// Returns k such that k heuristically is an upper bound on the
	// largest absolute accuracy needed in the rgcd algorithm.
	// 2^m should be a lower bound on the regulator and
	// c an upper bound on the absolute values of Ln(q), where q
	// runs over all principal ideal generators, which form the
	// bases in the power products of the units.
	//
	//  k = b(c) + max(j) b(sum |e_ij|) - 2m + 9
	//

{
	bigint col_sum, max_col_sum;
	bigint big_k;
	long   k;

	lidia_size_t i, j;
	lidia_size_t nofr, nofc;


	// initialize

	nofr = M.get_no_of_rows();
	nofc = M.get_no_of_columns();
	max_col_sum.assign_zero();


	// determine column sum norm

	for (j = 0; j < nofc; j++) {
		col_sum.assign_zero();
		for (i = 0; i < nofr; i++) {
			// add absolute value
			if (M.member(i, j).is_negative())
				subtract(col_sum, col_sum, M.member(i, j));
			else
				add(col_sum, col_sum, M.member(i, j));
		}

		// compare max
		if (col_sum > max_col_sum)
			max_col_sum.assign(col_sum);
	}


	// determine precision

	big_k = bigint(b_value(c)) + bigint(b_value(max_col_sum));
	big_k += 9 - 2 * bigint(m);

	if (big_k.longify(k)) {
		lidia_error_handler ("quadratic_order3.c::"
				     "estimate_pp_accuray",
				     "Precision overflow.");
	}

	return k;
}



void quadratic_order
::precompute_logarithm_bases (
	base_vector< xbigfloat > & log_q,
	base_vector< long > & prec_log_q,
	const bigint_matrix & M,
	const base_vector< quadratic_number_standard > & q,
	long k)

	//
	// For each i, determine an absolute prec_log_q[i] approximation
	// log_q[i] to Ln(q[i]), such that sum e_ij log_q[i] is an absolute
	// k approximation to sum e_ij Ln(q[i]) for each j, i.e.,
	//
	//  prec_log_q[i] = k + 2 + b(nofr) + max(j) b(e_ij).
	//
	// nofr = M.get_no_of_rows
	// nofc = M.get_no_of_colums
	//
	// 0 <= i < nofr,  0 <= j < nofc
	//

{
	lidia_size_t nofr, nofc, i, j;
	long t, kbn2;
	bigint max_row_entry;

	// initialize

	nofr = M.get_no_of_rows();
	nofc = M.get_no_of_columns();

	if (nofr != q.get_size()) {
		lidia_error_handler ("quadratic_order3.c::precompute_logarithm_bases",
				     "Number of exponents does not match number of "
				     "quadratic numbers.");
	}

	log_q.set_capacity(nofr);
	prec_log_q.set_capacity(nofr);


	//
	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(nofr)+2
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler ("quadratic_order3.c::precompute_logarithm_bases",
				     "precision overflow.");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(nofr))) {
		lidia_error_handler ("quadratic_order3.c::precompute_logarithm_bases",
				     "precision overflow.");
		return;
	}


	//
	// Determine the approximation for each quadratic number q[i].
	//

	for (i = 0; i < nofr; i++) {

		// find abs max row entry
		max_row_entry.assign_zero();

		for (j = 0; j < nofc; j++)
			if (max_row_entry.abs_compare(M.member(i, j)) < 0)
				max_row_entry.assign(M.member(i, j));

		// determine precision
		// t = kbn2 + b_value(max_row_entry);

		if (check_overflow(t, kbn2, b_value(max_row_entry))) {
			lidia_error_handler ("quadratic_order3.c::precompute_logarithm_bases",
					     "precision overflow.");
			return;
		}

		// approximate Ln(q[i])

		if (t < -1) t = -1;
		q[i].get_absolute_Ln_approximation_old(log_q[i], t);
		prec_log_q[i] = t;
	}
}


void quadratic_order
::reduce_modulo_regulator (
	bigint & q,
	const base_vector< quadratic_number_standard > & base_w,
	const base_vector< bigint > & exp_w,
	base_vector< xbigfloat > & log_base_w,
	base_vector< long > & prec_log_base_w,
	const base_vector< quadratic_number_standard > & base_r,
	const base_vector< bigint > & exp_r,
	base_vector< xbigfloat > & log_base_r,
	base_vector< long > & prec_log_base_r,
	long br)

	//
	// Reduction of Ln(w) modulo the regulator Ln(r).
	//
	//
	//  The function computes a quotient q such that
	//
	//   | floor(Ln(w) / Ln(r)) - q | <= 1,
	//
	//  if |Ln(w)| is larger than |Ln(r)| such that the condition
	//  below is satisfied. Otherwise it sets q = 0. Note that
	//  the condition is satisfied, if |Ln(w)| >= 32 |Ln(r)|, but
	//  the constant will usually be smaller than 32 in practice.
	//
	//  w and r are given as power product where base_... are the
	//  base elements and exp_... are the corresponding exponents.
	//  log_base_w[i] is an absolute prec_log_base_w[i] approximation
	//  to the Ln of the base_w[i], if the size of the log vector
	//  is the same as the size of the vector for the base elements.
	//  Same holds for r. If larger precisions are required during the
	//  computation of q, the logarithms are recomputed.
	//
	//  It is assumed that |b(Ln(r)) - br| <= 1.
	//
	//  Ln(w) must not be equal to zero.
	//

{
	xbigfloat log_w, log_r;
	long k, bw;
	bigint big_k;

	// compute absolute -br+6 approx. to Ln(w).
	absolute_Ln_approximation (log_w, log_base_w, prec_log_base_w, exp_w, base_w, -br+6);

	// If |log_w| < 2^(br+2), q = 0.
	bw = log_w.b_value();

	if (bw <= (br+2)) {
		q.assign_zero();
		return;
	}

	// Otherwise, compute absolute -2*br+bw+8 approx. to Ln(r)
	big_k = bigint(8);
	big_k += bigint(bw);
	big_k -= bigint(2) * bigint(br);

	if (big_k.longify(k)) {
		lidia_error_handler ("quadratic_order3.c::reduce_modulo_regulator",
				     "Precision overflow.");
		return;
	}

	absolute_Ln_approximation (log_r, log_base_r, prec_log_base_r, exp_r, base_r, k);

	// determine quotient
	divide (q, log_w, log_r);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
