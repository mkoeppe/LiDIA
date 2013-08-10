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
//	$Id: qo_sieve.cc,v 2.9 2004/06/15 10:20:00 lidiaadm Exp $
//
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/number_fields/qo_sieve.h"
#include	"LiDIA/udigit.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/matrix_GL2Z.h"
#include	"LiDIA/random_generator.h"
#include	<cstdio>
#include	<unistd.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void concat (char* , char* , char*);
void sort_file_line(char*, char*);

//
// utility functions
//


//
// qo_reverse
//
// Task:
//      reverses s
//

inline void
qo_reverse(char *s)
{
	debug_handler("qo_sieve", "qo_reverse");

	int c, i, j;

	for (i = 0, j = strlen(s) - 1; i < j; ++i, --j) {
		c = *(s + i);
		*(s + i) = *(s + j);
		*(s + j) = c;
	}
}



//
// qo_itoa
//
// Task:
//      converts n to a string
//

inline void
qo_itoa(unsigned long n, char *s, int base)
{
	debug_handler("qo_sieve", "qo_itoa");

	unsigned long   i, digit, dig16;

	if (base > 16) {
		lidia_error_handler("qo_itoa", "itoa - maxbase = 16");
		return;
	}

	i = 0;
	do {
		digit = n % base + '0';
		if (base == 16) {
			dig16 = digit - '0';
			switch (dig16) {
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
				s[i++] = (char) 'a' + (char) dig16 - 10;
				break;
			default:
				s[i++] = (char) digit;
				break;
			}
		} else
			s[i++] = (char) digit;
	} while ((n /= base) > 0);
	s[i] = '\0';
	qo_reverse(s);
}



//
// qo_get_name
//
// Task:
//      computes the name of a temporary file
//

void
qo_sieve::qo_get_name(char *s)
{
	debug_handler("qo_sieve", "qo_get_name");

	FILE *fp;
	static int first_one = 1;
	static char *hilf;
	char *physical;

	physical = new char[100];

	if (first_one) {
		first_one = 0;
		hilf = new char[100];

		if ((fp = fopen("/tmp/FAKT", "w")) == NULL) {
			strcpy(hilf, "./");
			qo_itoa(getpid(), physical, 10);
			strcat(hilf, physical);
		}
		else {
			fclose(fp); std::remove("/tmp/FAKT");
			strcpy(hilf, "/tmp/");

			qo_itoa(getpid(), physical, 10);
			strcat(hilf, physical);
		}
	}

	if (!strcmp(s, "PART")) {
		strcpy(physical, hilf);
		strcat(physical, "QO_PART");
		strcpy(PART, physical);
	}
	if (!strcmp(s, "TEMP")) {
		strcpy(physical, hilf);
		strcat(physical, "QO_TEMP");
		strcpy(TEMP, physical);
	}
	if (!strcmp(s, "ZUSA")) {
		strcpy(physical, hilf);
		strcat(physical, "QO_ZUSA");
		strcpy(ZUSA, physical);
	}

	delete [] physical;
}



//
// qo_count_ones
//
// Task:
//      determine the number of ones in binary representation of bi
//

inline int
qo_count_ones(unsigned int bi)
{
	int s = 0;

	while (bi > 0) {
		if (bi & 1)
			++s;
		bi >>= 1;
	}

	return s;
}



//
// new_polynomial_si()
//
// Task:
//      computes a new polynomial for sieving using self-initialization
//

void
qo_sieve::new_polynomial_si()
{
	debug_handler("qo_sieve", "new_polynomial_si");

	int SIEBS, tmp, tmp1, tmp2, p, p2, B2p, Bpr2p;
	lidia_size_t j, nu_2, i, i1, i2;
	bigint reserve, TMP, A, B, Delta2, temp;
	int SIEBSTART = M;

	// use Delta/2 if Delta = 0 (mod 4)
	if (Delta.is_even())
		shift_right(Delta2, Delta, 2);
	else
		Delta2.assign(Delta);

	if (index_i == 0) {
		++bin_index;

		while (qo_count_ones(bin_index) != P_ONCE)
			++bin_index;

		// reset exponent vector
		if (Q_prime_si.size() > 0) {
			for (i = 0; i < Q_prime_si.size(); ++i)
				e_si[Q_prime_si[i]] = 0;
		}

		i = 0;
		for (j = 0; j < P_TOTAL; ++j) {  // determine primes used for A
			// in this iteration
			if (bin_index & (1 << j))
				Q_prime_si[i++] = Q_prime_glob[j];
		}

		A.assign (FB[Q_prime_si[0]]); // compute coefficient A
		for (i = 1; i < P_ONCE; ++i)
			multiply (A, A, FB[Q_prime_si[i]]);

		// compute BG[0] to BG[P_ONCE-1]
		// each B is of the form
		// v_0*BG[0]+v_1*BG[1]+...+v_(P_ONCE-1)*BG[(P_ONCE-1)],
		// where each v_j is +1 or -1

		// this has to be done only once (for index_i==0) for the
		// coefficient A; if index_i > 0 then there is a linear
		// recursion for B

		for (i = 0; i < P_ONCE; ++i) {
			p = FB[Q_prime_si[i]];
			divide (reserve, A, p);
			tmp = static_cast<int>(remainder(reserve, p));
			if (tmp < 0)
				tmp += p;
			multiply (reserve, reserve, invert(tmp, p));

			if (Delta2.is_gt_zero())
				tmp = ressol (static_cast<int>(remainder(Delta2, p)), p);
			else
				tmp = ressol (p + static_cast<int>(remainder(Delta2, p)), p);

			if (tmp < 0)
				tmp += p;
			multiply (reserve, reserve, tmp);
			remainder (BG[i], reserve, A);

			add(B, B, BG[i]); // compute actual B coefficient
			if (Delta.is_even())
				BG[i].multiply_by_2();
		}

		if ((B.is_even()) && (Delta.is_odd()))          // assure   B = 1 mod 4
			add(B, (A.least_significant_digit() & 3)*A, B);
		if (Delta.is_even())
			B.multiply_by_2();

		square(reserve, B); // compute coefficient -C
		subtract(reserve, reserve, Delta);
		shift_left(TMP, A, 2);
		divide (TMP, reserve, TMP);
		F_si.assign(A, B, TMP);

		// a_inv[i] = 1/(2*A) mod p_i
		for (i = 0; i < size_FB; ++i)
			a_inv[i] = invert (static_cast<int>(remainder(A << 1, static_cast<long>(FB[i]))), FB[i]);

		// vorb[i][j] = 1/A * B[i] mod p_j
		for (i = 0; i < P_ONCE; ++i) {
			for (j = 0; j < size_FB; ++j) {
				p = FB[j];
				multiply (reserve, BG[i], a_inv[j] << 1);
				if ((tmp = static_cast<int>(remainder(reserve, p))) < 0)
					tmp += p;
				vorb[i][j] = tmp;
			}
		}

		for (j = 0; j < size_FB; ++j) {
			p = FB[j];
			SIEBS = SIEBSTART % p;
			tmp = static_cast<int>(remainder(-B, p));

			if ((D_primes.size() > 0) &&
			    ((D_primes.bin_search(p, i1)) || (D_primes.bin_search(-p, i2)))) {
				// p divides Delta - root = -b * (2a)^-1 (mod p)

				if (tmp < 0)
					tmp += p;

				if (p == 2)
					if (F_si.get_c().is_odd())
						tmp2 = 1;
					else
						tmp2 = 0;
				else
					tmp2 = static_cast<int>(multiply_mod(
						static_cast<udigit>(tmp), udigit(a_inv[j]), static_cast<udigit>(p)));

				tmp2 = (tmp2 + SIEBS) % p;
				if (tmp2 < 0)
					START1[j] = tmp2 + p;
				else
					START1[j] = tmp2;
				START2[j] = START1[j];


			}
			else {
				tmp1 = (tmp - SQRTkN[j]) % p;
				if (tmp1 < 0)
					tmp1 += p;
				tmp = (tmp + SQRTkN[j]) % p;
				if (tmp < 0)
					tmp += p;

				tmp2 = static_cast<int>(multiply_mod(
					static_cast<udigit>(tmp), udigit(a_inv[j]), static_cast<udigit>(p)));

				tmp2 = (tmp2 + SIEBS) % p;
				if (tmp2 < 0)
					START1[j] = tmp2 + p;
				else
					START1[j] = tmp2;

				tmp2 = static_cast<int>(multiply_mod(
					static_cast<udigit>(tmp1), udigit(a_inv[j]), static_cast<udigit>(p)));
				tmp2 = (tmp2 + SIEBS) % p;
				if (tmp2 < 0)
					START2[j] = tmp2 + p;
				else
					START2[j] = tmp2;

			}
		}

		// compute signs of exponents
		i = 0;
		for (j = 0; j < P_ONCE; ++j) {
			p = FB[Q_prime_si[j]];

			p2 = (p << 1);
			B2p = static_cast<int>(remainder(F_si.get_b(), p2));
			if (B2p < 0)
				B2p += p2;

			while (FB[i] != p)
				++i;
			temp = FB_b[i];
			Bpr2p = static_cast<int>(remainder(temp, p2));
			if (Bpr2p < 0)
				Bpr2p += p2;
			if (B2p == Bpr2p)
				e_si[Q_prime_si[j]] = 1;
			else
				e_si[Q_prime_si[j]] = -1;
		}
	}
	else {  // no "real" computation -- use recursive formula
		// first:   update of B, compute B[index_i], index_i > 0

		A.assign(F_si.get_a());
		B.assign(F_si.get_b());

		nu_2 = 0; // nu_2 = nu_2(index_i)
		j = index_i;
		while ((j & 1) == 0) {
			++nu_2;
			j >>= 1;
		}

		// negate appropriate exponent
		e_si[Q_prime_si[nu_2]] = -e_si[Q_prime_si[nu_2]];

		shift_left(TMP, BG[nu_2], 1);

		if ((((j+1)/2) & 1) == 1) {
			i = -1;
			subtract(B, B, TMP);
		}
		else {
			i = 1;
			add(B, B, TMP);
		}

		square(reserve, B); // compute coefficient -C
		subtract(reserve, reserve, Delta);
		shift_left(TMP, A, 2);
		divide (TMP, reserve, TMP);
		F_si.assign(A, B, TMP);

		// determine new starting positions for sieving
		if (i == -1) {
			for (j = 0; j < size_FB; ++j) {
				p = FB[j];
				START1[j] += vorb[nu_2][j];
				if (START1[j] >= p)
					START1[j] -= p;
				START2[j] += vorb[nu_2][j];
				if (START2[j] >= p)
					START2[j] -= p;
			}
		}
		else {
			for (j = 0; j < size_FB; ++j) {
				p = FB[j];
				START1[j] -= vorb[nu_2][j];
				if (START1[j] < 0)
					START1[j] += p;
				START2[j] -= vorb[nu_2][j];
				if (START2[j] < 0)
					START2[j] += p;
			}
		}

	}

// now compute zeros of polynomials that have only one zero mod p
// because p divides coefficient A

	TMP.assign(-F_si.get_c());

	for (j = 0; j < P_ONCE; ++j) {
		p = FB[Q_prime_si[j]];

		tmp = invert(static_cast<int>(remainder(B, p)), p);
		if (tmp < 0)
			tmp += p;
		tmp2 = static_cast<int>(remainder(TMP, p));
		if (tmp2 < 0)
			tmp2 += p;

		tmp = static_cast<int>(multiply_mod(static_cast<udigit>(tmp2), static_cast<udigit>(tmp),
							   static_cast<udigit>(p)));
		START1[Q_prime_si[j]] = START2[Q_prime_si[j]] = (tmp + SIEBSTART) % p;
	}

#ifdef LIDIA_DEBUG
	if (F_si.discriminant() != Delta) {
		std::cerr << "Wrong coefficients!  F = " << F_si << std::endl;
		std::cerr << "Delta = " << F_si.discriminant() << std::endl;
		std::cerr << "Q_prime = " << Q_prime_si << std::endl;
		std::cerr << " = [ ";
		for (i = 0; i < Q_prime_si.size(); ++i)
			std::cerr << FB[Q_prime_si[i]] << " ";
		std::cerr << "]" << std::endl;
		std::cerr << "e = " << e_si << std::endl;
		exit(1);
	}

	quadratic_form FF, FP;
	FF.assign_one(Delta);
	for (j = 0; j < Q_prime_si.size(); ++j) {
		generate_prime_form(FP, FB[Q_prime_si[j]], Delta);
		if (e_si[Q_prime_si[j]] < 0)
			compose_regular(FF, FF, get_conjugate(FP));
		else
			compose_regular(FF, FF, FP);
	}
	FP.assign(F_si);
	FP.normalize_regular();
	FF.normalize_regular();
	if (FP != FF) {
		std::cerr << "Wrong polynomial!  F = " << F_si << std::endl;
		std::cerr << "            computed = " << FF << std::endl;
		std::cerr << "Q_prime = " << Q_prime_si << std::endl;
		std::cerr << " = [ ";
		for (i = 0; i < Q_prime_si.size(); ++i)
			std::cerr << FB[Q_prime_si[i]] << " ";
		std::cerr << "]" << std::endl;
		std::cerr << "e = " << e_si << std::endl;
		exit(1);
	}

	for (j = 0; j < size_FB; ++j) {
		p = FB[j];
		SIEBS = SIEBSTART % p;
		if (F_si.eval(START1[j]-SIEBS, 1) % p != 0) {
			std::cerr << "\n new_polynomial_si: found wrong polynomial (1).\n\n";
			std::cerr << "p = " << p << ", F = " << F_si << std::endl;
			std::cerr << "Q_prime = " << Q_prime_si << std::endl;
			std::cerr << " = [ ";
			for (i = 0; i < Q_prime_si.size(); ++i)
				std::cerr << FB[Q_prime_si[i]] << " ";
			std::cerr << "]" << std::endl;
			std::cerr << "e = " << e_si << std::endl;
			std::cerr << "START1[j] = " << START1[j] << ", SIEBS = " << SIEBS << std::endl;
			std::cerr << "F[r] = " << F_si.eval(START1[j]-SIEBS, 1) % p << std::endl;
			exit(1);
		}
		if (F_si.eval(START2[j]-SIEBS, 1) % p != 0) {
			std::cerr << "\n new_polynomial_si: found wrong polynomial (2).\n\n";
			std::cerr << "p = " << p << ", F = " << F_si << std::endl;
			std::cerr << "Q_prime = " << Q_prime_si << std::endl;
			std::cerr << " = [ ";
			for (i = 0; i < Q_prime_si.size(); ++i)
				std::cerr << FB[Q_prime_si[i]] << " ";
			std::cerr << "]" << std::endl;
			std::cerr << "e = " << e_si << std::endl;
			std::cerr << "START2[j] = " << START2[j] << ", SIEBS = " << SIEBS << std::endl;
			std::cerr << "F[r] = " << F_si.eval(START2[j]-SIEBS, 1) % p << std::endl;
			exit(1);
		}
	}

#endif
}



//
// new_polynomial
//
// Task:
//      computes a new polynomial for sieving by selectin a random exponent
//      vector
//

void
qo_sieve::new_polynomial(lidia_size_t force_index, int B)
{
	debug_handler("qo_sieve", "new_polynomial");

	lidia_size_t i, j;
	quadratic_form FF;
	bigint *oldpoly;
	bigint poly_pair;
	int p, e, i1, i2, ex, rand_size, attempts;
	bigint temp;
	bigint dbl_tmp;
	random_generator rg;

	for (i = 0; i < size_FB; ++i)
		e_ra[i] = 0;

	oldpoly = (bigint *) true;
	Q_prime_ra.reset();
	j = 0;
	if (force_index >= 0) {
		e_ra[force_index] = -1;
		Q_prime_ra[j] = force_index;
		F_ra.assign(qi_class(bigint(FB[force_index]), FB_b[force_index]));
		F_ra.conjugate();
		F_ra.normalize_regular();
		++j;
	}
	else {
		qi_class A;
		A.assign_one();
		F_ra.assign(A);
	}
	dbl_tmp.assign(size_A);
	if (B <= 0)
		rand_size = size_FB;
	else
		rand_size = B;
	attempts = 0;

	while (oldpoly) {
		do {
			rg >> i;
			i %= rand_size;
			p = FB[i];
			FF.assign(qi_class(bigint(p), FB_b[i]));

			if ((!Q_prime_ra.linear_search(i, i1)) && (!D_primes.linear_search(p, i2))) {
				rg >> e;
				e %= 2;
				if (e) {
					compose_regular(F_ra, F_ra, FF);
					ex = 1;
				}
				else {
					compose_regular(F_ra, F_ra, get_conjugate(FF));
					ex = -1;
				}

				e_ra[i] = ex;
				Q_prime_ra[j] = i;
				++j;

				F_ra.normalize_regular();
			}
			else {
				++attempts;
				if (attempts > 5) {
					rand_size += 5;
					if (rand_size > size_FB)
						rand_size = size_FB;
					attempts = 0;
				}
			}

		} while (F_ra.get_a() < dbl_tmp);

		shift_left(poly_pair, F_ra.get_a(), 1);
		add(poly_pair, poly_pair, F_ra.get_b());
		oldpoly = used_polys.search(poly_pair);

		++attempts;
		if (attempts > 5) {
			rand_size += 5;
			if (rand_size > size_FB)
				rand_size = size_FB;
			attempts = 0;
		}
	}

	if (B == pBound)
		pBound = rand_size;

	Q_prime_ra.sort();
}



//
// compute_roots
//
// Task:
//      computes the roots of the polynomial given by F, as well as the
//      start values for sieving.
//

void
qo_sieve::compute_roots()
{
	debug_handler("qo_sieve", "compute_roots");

	bigint A, B, TMP;
	int p, SIEBS, i1, i2, tmp, tmp1, tmp2, AINV, divA;
	lidia_size_t j;
	int SIEBSTART = M;

	A.assign(F_ra.get_a());
	B.assign(F_ra.get_b());
	for (j = 0; j < size_FB; ++j) {
		p = FB[j];
		SIEBS = SIEBSTART % p;
		divA = static_cast<int>(remainder(A, p));

		if ((D_primes.size() > 0) &&
		    ((D_primes.bin_search(p, i1)) || (D_primes.bin_search(-p, i2)))) {
			// p divides Delta - root = -b * (2a)^-1 (mod p)

			if (!divA) {
				// no roots mod p
				START1_2[j] = START2_2[j] = -1;
			}
			else {
				AINV = invert (static_cast<int>(remainder(A << 1, static_cast<long>(p))), p);
				tmp = static_cast<int>(remainder(-B, p));

				if (tmp < 0)
					tmp += p;

				if (p == 2)
					if (F_ra.get_c().is_odd())
						tmp2 = 1;
					else
						tmp2 = 0;
				else
					tmp2 = static_cast<int>(multiply_mod(
						static_cast<udigit>(tmp), udigit(AINV), static_cast<udigit>(p)));

				tmp2 = (tmp2 + SIEBS) % p;
				if (tmp2 < 0)
					START1_2[j] = tmp2 + p;
				else
					START1_2[j] = tmp2;
				START2_2[j] = START1_2[j];
			}
		}

		else if (!divA) {
			TMP.assign(-F_ra.get_c());

			tmp = invert (static_cast<int>(remainder(B, p)), p);
			if (tmp < 0)
				tmp += p;
			tmp2 = static_cast<int>(remainder (TMP, p));
			if (tmp2 < 0)
				tmp2 += p;

			tmp = static_cast<int>(multiply_mod(static_cast<udigit>(tmp2),
								   static_cast<udigit>(tmp),
								   static_cast<udigit>(p)));
			START1_2[j] = START2_2[j] = (tmp + SIEBSTART) % p;
		}


		else {
			AINV = invert (static_cast<int>(remainder(A << 1, static_cast<long>(p))), p);
			tmp = static_cast<int>(remainder(-B, p));

			tmp1 = (tmp - SQRTkN[j]) % p;
			if (tmp1 < 0)
				tmp1 += p;
			tmp = (tmp + SQRTkN[j]) % p;
			if (tmp < 0)
				tmp += p;

			tmp2 = static_cast<int>(multiply_mod(
				static_cast<udigit>(tmp), udigit(AINV), static_cast<udigit>(p)));

			tmp2 = (tmp2 + SIEBS) % p;
			if (tmp2 < 0)
				START1_2[j] = tmp2 + p;
			else
				START1_2[j] = tmp2;

			tmp2 = static_cast<int>(multiply_mod(
				static_cast<udigit>(tmp1), udigit(AINV), static_cast<udigit>(p)));
			tmp2 = (tmp2 + SIEBS) % p;
			if (tmp2 < 0)
				START2_2[j] = tmp2 + p;
			else
				START2_2[j] = tmp2;

		}
	}

#ifdef LIDIA_DEBUG
	for (j = 0; j < size_FB; ++j) {
		p = FB[j];
		SIEBS = SIEBSTART % p;
		if (START1_2[j] >= 0) {
			if (F_ra.eval(START1_2[j]-SIEBS, 1) % p != 0) {
				std::cerr << "\n compute_roots: found wrong polynomial (1).\n\n";
				std::cerr << "p = " << p << ", F = " << F_ra << std::endl;
				std::cerr << "START1[j] = " << START1_2[j] << ", SIEBS = ";
				std::cerr << SIEBS << std::endl;
				exit(1);
			}
			if (F_ra.eval(START2_2[j]-SIEBS, 1) % p != 0) {
				std::cerr << "\n compute_roots: found wrong polynomial (2).\n\n";
				std::cerr << "p = " << p << ", F = " << F_ra << std::endl;
				std::cerr << "START2[j] = " << START2_2[j] << ", SIEBS = ";
				std::cerr << SIEBS << std::endl;
				exit(1);
			}
		}
		else {
			if (remainder(A, p) || remainder(B, p) || !remainder(F_ra.get_c(), p)) {
				std::cerr << "\n compute_roots: found wrong polynomial (-1).\n\n";
				std::cerr << "p = " << p << ", F = " << F_ra << std::endl;
				std::cerr << "START1[j] = " << START1_2[j] << std::endl;
				exit(1);
			}
		}
	}

#if 0
	quadratic_form f1, f2;
	f1.assign_one(Delta);
	for (j = 0; j < Q_prime_ra.size(); ++j) {
		generate_prime_form(f2, FB[Q_prime_ra[j]], Delta);
		if (e_ra[Q_prime_ra[j]] < 0)
			compose_regular(f1, f1, get_conjugate(f2));
		else
			compose_regular(f1, f1, f2);
	}
	f1.normalize_regular();
	f2.assign(F_ra);
	f2.normalize_regular();
	if (f1 != f2) {
		std::cerr << "ERROR:  f2 <> f1!!!" << std::endl;
		std::cerr << "f1 = " << f1 << std::endl;
		std::cerr << "f2 = " << f2 << std::endl;
		std::cerr << "Q_prime = " << Q_prime_ra << std::endl;
		std::cerr << "e = " << e_ra << std::endl;
		exit(1);
	}
#endif
#endif
}


// FIXME: not portable!
#ifdef WORDS_BIGENDIAN
// big endian
# define INT32_OCTET_0_0x80	0x80000000U
# define INT32_OCTET_1_0x80	0x00800000U
# define INT32_OCTET_2_0x80	0x00008000U
# define INT32_OCTET_3_0x80	0x00000080U
#else
// little endian
# define INT32_OCTET_0_0x80	0x00000080U
# define INT32_OCTET_1_0x80	0x00008000U
# define INT32_OCTET_2_0x80	0x00800000U
# define INT32_OCTET_3_0x80	0x80000000U
#endif

//
// qo_sieve_interval
//
// Task:
//      sieves the interval given the start values, roots, and factor base
//

void
qo_sieve::qo_sieve_interval()
{
	debug_handler("qo_sieve", "qo_sieve_interval");

	register int p, l, *fbp, *lsieb = (int *)sieb;	// FIXME: SIEBTYP and int are incompatible types
	register SIEBTYP logp;
	register SIEBTYP *begin;
	register int x, counter = 0, M_2 = M << 1;
	register int oldstart1;
	register int *ST1, *ST2;

	if (next_is_si) {
		ST1 = START1;
		ST2 = START2;
	}
	else {
		ST1 = START1_2;
		ST2 = START2_2;
	}

	memset(sieb, '\0', (M_2) * sizeof(SIEBTYP));
	fbp = &FB[smallstart];
	l = smallstart;

	while ((p = *fbp++) != 0) {
		if (ST1[l] >= 0) {
			logp = LOGP[l];
			begin = sieb + ST1[l];
			oldstart1 = ST1[l];

			for (;;) {
				if (begin <= ende) {
					(*begin) += logp;
					begin += p;
				}
				else
					break;
			}

			if (oldstart1 != ST2[l]) {
				begin = sieb + ST2[l];
				for (;;) {
					if (begin <= ende) {
						(*begin) += logp;
						begin += p;
					}
					else
						break;
				}
			}
		}

		++l;
	}

	l = 0;
	while (l < M_2) {
		// check whether at least one of four consequent sieve entries
		// is a candidate
		if ((p = (*lsieb))& (0x80808080)) {
			x = l;

			if (p & (INT32_OCTET_0_0x80 | INT32_OCTET_1_0x80)) {
				if (p & INT32_OCTET_0_0x80) {
					CANDIDATE[counter++] = x;
					if (p & INT32_OCTET_1_0x80)
						CANDIDATE[counter++] = x+1;
				}
				else
					CANDIDATE[counter++] = x+1;

				if (p & INT32_OCTET_2_0x80)
					CANDIDATE[counter++] = x+2;

				if (p & INT32_OCTET_3_0x80)
					CANDIDATE[counter++] = x+3;
			}
			else {
				if (p & INT32_OCTET_2_0x80)
					CANDIDATE[counter++] = x+2;
				else
					CANDIDATE[counter++] = x+3;
			}
		}
		++lsieb;
		l += 4;
	}
	CANDIDATE[counter] = 0;
}



//
// test_candidates
//
// Task:
//      tests if the candidates found by sieving are relations.  Any relations
//      found are appended to the end of the relation array.
//

bool
qo_sieve::test_candidates(int needed, lidia_size_t force_index, bool get_vec) {
	debug_handler("qo_sieve", "test_candidates");

	int *ST1, *ST2;
	matrix_GL2Z UT;
	quadratic_form F, G, PF;
	bigint u, v, H, Qx, temp, big_p2, big_B2p, big_Bpr2p;
	int div2, fak_i, vorber, x, p, divides, counter, p2, B2p, Bpr2p, rest_i;
	long rest;
	short small_value = 0;
	lidia_size_t i, k, ra_cols = relation_array.get_no_of_columns();
	int i1, i2, exp_i;
	bool found_force, is_divisor, good_one, one_entry;
	math_vector< int > e(size_FB, size_FB), vec(size_FB, size_FB);
	sort_vector< int > Q_prime;
	partial_relation part_rel;

	quadratic_number_standard pig; // principal ideal generator
	bigint pig_h, pig_d;

	if (force_index >= 0) {
		found_force = false;
		needed = 1000;
	}
	else
		found_force = true;

	if (next_is_si) {
		F = F_si;
		e = e_si;
		ST1 = START1;
		ST2 = START2;
		Q_prime = Q_prime_si;
	}
	else {
		F = F_ra;
		e = e_ra;
		ST1 = START1_2;
		ST2 = START2_2;
		Q_prime = Q_prime_ra;
	}


				// while there are candidates to test
	counter = 0;
	while (((x = CANDIDATE[counter++]) != 0) &&
	       ((small_value < needed) || (needed < 0) || (!found_force))) {
		x -= M;
		Qx.assign(F.eval(bigint(x), bigint(1)));
		Qx.absolute_value();
		multiply(H, F.get_a(), bigint(x));
		H.multiply_by_2();
		add(H, H, F.get_b());
		remainder(H, H, abs(Delta));
		if (H.is_lt_zero())
			add(H, H, abs(Delta));

		xgcd(u, v, bigint(x), bigint(1));
		UT = matrix_GL2Z(bigint(x), -v, bigint(1), u);
		G.assign(F);
		G.transform(UT);
		G.conjugate();

				// compute powers for primes p
		fak_i = 0;
		one_entry = false;
		k = 0;

				// initialize relation vector
		for (i = 0; i < size_FB; ++i)
			vec[i] = 0;


		pig_d.assign_one();
		while ((p = FB[fak_i]) != 0) {
			div2 = divides = 0;

			// check if p divides norm of ideal -- adjust exponent accordingly
			if (Q_prime.linear_search(fak_i, i1))
				divides += e[Q_prime[i1]];

			if (ST1[fak_i] >= 0) {
				vorber = (M + x) %p;

				if (p == 2)
					is_divisor = (Qx.is_even());
				else
					is_divisor = ((vorber == ST1[fak_i]) || (vorber == ST2[fak_i]));

				// if p divides Qx
				if (is_divisor) {
					div2 = 0;
					do {
						div_rem(temp, rest, Qx, static_cast<long>(p));
						if (rest == 0) {
							++div2;
							Qx.assign(temp);
						}
					} while (rest == 0);

					if (get_vec) {
						// compute sign of exponent
						p2 = (p << 1);
						B2p = static_cast<int>(remainder(G.get_b(), p2));
						if (B2p < 0)
							B2p += p2;

						temp = FB_b[fak_i];
						Bpr2p = static_cast<int>(remainder(temp, p2));
						if (Bpr2p < 0)
							Bpr2p += p2;
						if (B2p != Bpr2p)
							div2 = -div2;
					}
					divides += div2;
				}

				if ((D_primes.size() > 0) &&
				    ((D_primes.bin_search(p, i1)) || (D_primes.bin_search(-p, i2)))) {
					// p divides Delta -- take exponent mod 2
					if (divides & 1)
						divides = 1;
					else
						divides = 0;
				}
			}

			if (divides) {
				if (fak_i == force_index)
					found_force = true;
				one_entry = true;
				vec[fak_i] = divides;

				if (divides > 0) {
					power(temp, p, divides);
					multiply(pig_d, pig_d, temp);
				}
			}

			++fak_i;
		}

				// to avoid duplication of relations (which might cause problems in the
				// linear equation system) use hash function on H(x)

		if ((Qx.is_one()) && (one_entry)) {
			// exit if only existance of a relation is required
			if (!get_vec)
				return true;

			if (Delta.is_gt_zero()) {
				// add minimum to the list
				// pig = ((2*x*a + b) + sqrt(Delta) / 2)
				multiply (pig_h, F.get_a(), x);
				pig_h.multiply_by_2();
				add (pig_h, pig_h, F.get_b());
				pig_d.multiply_by_2();
				pig.assign(pig_h, bigint(1), pig_d);
			}
			else
				pig.assign_one();

			if (((force_index < 0) || found_force) && !htable.search(H)) {
				htable.hash(H);
				relation_array.set_no_of_columns(ra_cols+1);
				relation_array.sto_column_vector(vec, size_FB, ra_cols);
				minima_array[ra_cols] = pig;
				++ra_cols;
				++small_value;
			}
		}
		else if (((Qx.intify(rest_i)) == 0) &&
			 ((rest_i < BachBound) || (rest_i < lpbound))) {
			if (!get_vec) {
				if (lp_table.search(rest_i))
					return true;
			}

			if (is_prime(bigint(rest_i), 1)) {
				if (!lp_table.search(rest_i))
					lp_table.hash(rest_i);

				if ((rest_i < lpbound) && get_vec) {
					// compute sign of exponent
					exp_i = 0;
					if (generate_prime_form(PF, bigint(rest_i), Delta)) {
						big_p2.assign(rest_i);
						big_p2.multiply_by_2();
						remainder(big_B2p, G.get_b(), big_p2);
						if (big_B2p.is_lt_zero())
							big_B2p += big_p2;

						temp = PF.get_b();
						remainder(big_Bpr2p, temp, big_p2);
						if (big_Bpr2p.is_lt_zero())
							big_Bpr2p += big_p2;
						if (big_B2p != big_Bpr2p)
							exp_i = 1; // negative exponent
						else
							exp_i = 0; // positive exponent

						if (get_vec)
							good_one = true;
						else
							good_one = false;
					}
					else
						good_one = false;

					if (good_one && !htable.search(H)) {
						htable.hash(H);
						if (Delta.is_gt_zero()) {
							// compute minima
							// pig = ((2*x*a + b) + sqrt(Delta)) / 2
							multiply (pig_h, F.get_a(), x);
							pig_h.multiply_by_2();
							add (pig_h, pig_h, F.get_b());
							pig.assign(pig_h, bigint(1), pig_d);
						}
						else
							pig.assign_one();

						memset(faktor, '\0', LINE_LEN);
						for (i = 0; i < size_FB; ++i) {
							if (vec[i])
								sprintf(faktor, "%s %d %d", faktor, vec[i], i);
						}

						part_rel.assign(rest_i, exp_i, H, faktor, pig);
						lp_file << part_rel;
					}
				}
			}
		}
	}

	return (small_value > 0);
}




//
// test_single
//
// Task:
//      tests if one of the candidates found by sieving is a relation
//

bool
qo_sieve::test_single(math_vector< int > &vec, quadratic_number_standard & pig) {
	debug_handler("qo_sieve", "test_single");

	bool found;

	found = test_candidates(1, -1, true);

	if (found) {
				// get relation from arrays
		lidia_size_t ra_cols = relation_array.get_no_of_columns();
		relation_array.get_column_vector(vec, ra_cols-1);
		relation_array.set_no_of_columns(ra_cols-1);
		pig = minima_array[ra_cols-1];
		minima_array.set_size(ra_cols-1);
	}

	return found;
}




//
// constructor
//	- initializes parameters to zero
//

qo_sieve::qo_sieve() {
	debug_handler("qo_sieve", "qo_sieve()");

	Delta.assign_zero();
	is_init = is_si = false;

	size_FB = M = P_ONCE = P_TOTAL = smallstart = POLY = 0;
	lpbound = 0;

	Tval = 0.0;

	FB = static_cast<int *>(NULL);
	a_inv = static_cast<int *>(NULL);
	BG = static_cast<bigint *>(NULL);

	D_primes.set_mode(EXPAND);

	relation_array.set_representation(matrix_flags::sparse_representation);
	relation_array.set_orientation(matrix_flags::column_oriented);
	relation_array.set_zero_element(0);

	minima_array.set_mode(EXPAND);

	Q_prime_si.set_mode(EXPAND);
	Q_prime_ra.set_mode(EXPAND);
	Q_prime_glob.set_mode(EXPAND);

	ende = static_cast<SIEBTYP *>(0);

	strcpy(PART, "\n");
	strcpy(ZUSA, "\n");
	strcpy(TEMP, "\n");
}


//
// destructor
//	- frees dynamically allocated memory
//

qo_sieve::~qo_sieve() {
	debug_handler("qo_sieve", "~qo_sieve()");

	reset();

	if (FB) {
		delete [] FB; delete [] FB_b; delete [] LOGP; delete [] SQRTkN;
		delete [] START1; delete [] START2; delete [] sieb; delete [] CANDIDATE;
		delete [] START1_2; delete [] START2_2;
	}

	if (PART && (strlen(PART) > 0))
		std::remove(PART);
	if (TEMP && (strlen(TEMP) > 0))
		std::remove(TEMP);
	if (ZUSA && (strlen(ZUSA) > 0))
		std::remove(ZUSA);
}



//
// init
//
// Task:
//      initializes the sieve (allocate memory, set parameters)
//

void
qo_sieve::init(bigint & Delta_in,
	       base_vector< qi_class > & fact_base,
	       double Tval_in,
	       int M_in,
	       int P_ONCE_in,
	       int P_TOTAL_in,
	       int smallstart_in,
	       int POLY_in,
	       long lpbound_in) {
	debug_handler("qo_sieve", "init");

	register lidia_size_t i, j;
	int p, oldFB;
	bigfloat temp;
	bool use_old = false;

	if (!fact_base.size())
		return;

	if (Delta_in == Delta) {
		use_old = true;
		oldFB = size_FB;
	}
	else {
		oldFB = 0;
		D_primes.reset();
		lp_table.empty();
	}


				// allocate arrays if necessary
	if (!FB) {
		FB = new int[MAX_FB];
		FB_b = new bigint[MAX_FB];
		LOGP = new SIEBTYP[MAX_FB];
		SQRTkN = new int[MAX_FB];
		START1 = new int[MAX_FB];
		START2 = new int[MAX_FB];
		sieb = new SIEBTYP[MAX_SIEB];
		CANDIDATE = new int[CAND_NUMBER];
		START1_2 = new int[MAX_FB];
		START2_2 = new int[MAX_FB];
	}


				// destroy old initialization
	if (is_init) {
				// reset sieving arrays
		if (vorb) {
			for (i = 0; i < static_cast<int>(P_TOTAL); ++i)
				if (vorb[i])
					delete[] vorb[i];
			delete[] vorb;
			vorb = static_cast<int **>(NULL);
		}

				// reset hash tables
		htable.empty();
		used_polys.empty();

				// reset exponent vectors
		e_si.set_size(0);
		e_ra.set_size(0);

		Q_prime_glob.reset();
		Q_prime_si.reset();
		Q_prime_ra.reset();

				// reset parameters

		Delta.assign_zero();
		is_init = is_si = false;

		size_FB = M = P_ONCE = P_TOTAL = smallstart = POLY = 0;
		Tval = 0.0;

		relation_array.reset();
		minima_array.set_size(0);

				// delete large prime files
		if (lpbound) {
			lpbound = 0;

			lp_file.close();

			std::remove(PART);
			std::remove(TEMP);
			std::remove(ZUSA);
		}
	}
	is_si = false;

	is_init = true;


				// set parameters
	Delta.assign(Delta_in);
	size_FB = fact_base.size();
	Tval = Tval_in;
	M = M_in;
	P_ONCE = P_ONCE_in;
	P_TOTAL = P_TOTAL_in;
	smallstart = smallstart_in;
	POLY = POLY_in;
	lpbound = lpbound_in;

	square(temp, log(bigfloat(abs(Delta))));
	ceil(temp*bigfloat(12.0)).intify(BachBound);

	relation_array.set_no_of_rows(fact_base.size());
	minima_array.set_size(0);

	memset(sieb, '\0', (M << 1)*sizeof(SIEBTYP));
	memset(START1, '\0', (size_FB+1)*sizeof(int));
	memset(START2, '\0', (size_FB+1)*sizeof(int));
	memset(CANDIDATE, '\0', (CAND_NUMBER)*sizeof(int));
	memset(START1_2, '\0', (size_FB+1)*sizeof(int));
	memset(START2_2, '\0', (size_FB+1)*sizeof(int));
	if (oldFB == 0)
		memset(SQRTkN, '\0', (size_FB+1)*sizeof(int));
	if (oldFB != size_FB)
		memset(LOGP, '\0', (size_FB+1)*sizeof(SIEBTYP));


	ende = sieb + (M << 1) - 1;


				// compute optimal size of leading coefficient - should be about
				// sqrt(kN)/M
	d_wurz = std::sqrt(dbl(abs(Delta)));
	size_A = static_cast<double>(d_wurz) / (std::sqrt(static_cast<double>(2.0)) * M);


	if (size_FB != oldFB) {
		if (size_FB < oldFB) {
			FB[size_FB] = 0;
			FB_b[size_FB] = 0;
			SQRTkN[size_FB] = 0;
			LOGP[size_FB] = 0;
			for (i = 0; i < D_primes.size(); ++i) {
				if (D_primes[i] > FB[size_FB-1]) {
					D_primes.set_size(i);
					break;
				}
			}
		}
		else if (size_FB > oldFB) {
			// compute array of factor base ideal norms
			// also compute index bound for non-self-initialization
			int *fbb, isizeA = static_cast<int>(std::ceil(size_A));
			bigint *fbb_b, TA, maxp;
			bigfloat lDelta;

			lDelta = ceil(log(bigfloat(abs(Delta))));
			lDelta.bigintify(maxp);

			if (!use_old) {
				fbb = FB;
				fbb_b = FB_b;

				j = 0;
				D_primes.reset();
				pBound = 0;
			}
			else {
				fbb = &FB[oldFB];
				fbb_b = &FB_b[oldFB];

				j = D_primes.size();
			}

			TA = 1;
			for (i = oldFB; i < size_FB; ++i) {
				fact_base[i].get_a().intify(p);
				*fbb++ = p;
				*fbb_b++ = fact_base[i].get_b();
				if (remainder(Delta, p) == 0)
					D_primes[j++] = p;
				else {
					if (!use_old && (pBound == 0)) {
						multiply(TA, TA, bigint(p));
						if ((TA > isizeA) && (p > maxp)) {
							pBound = i;
							if (pBound < 4)  pBound = 4;
						}
					}
				}
			}
			*fbb = 0;
			*fbb_b = 0;
		}


				// compute logarithms and square roots of Delta
		double LOGMUL;

		LOGMUL = (2 * static_cast<SIEBTYP>(0.5 * LiDIA::log2 (dbl(abs(Delta))) +
						   LiDIA::log2 (static_cast<double>(M)) - Tval *
						   LiDIA::log2 (static_cast<double>(FB[size_FB - 1]))));
		LOGMUL = 127.0 / LOGMUL;

				// compute log(p) and sqrt(Delta) mod p for all p in the factor base
		for (i = 0; i < size_FB; ++i) {
			p = FB[i];
			LOGP[i] = static_cast<SIEBTYP>(LOGMUL*LiDIA::log2(static_cast<double>(p))*2);
			if (i >= oldFB) {
				if (Delta.is_gt_zero()) {
					if ((SQRTkN[i] = ressol(static_cast<int>(remainder(Delta, p)), p)) < 0)
						SQRTkN[i] += p;
				}
				else {
					if ((SQRTkN[i] = ressol(p + static_cast<int>(remainder(Delta, p)), p)) < 0)
						SQRTkN[i] += p;
				}
			}
		}

		SQRTkN[size_FB] = 0;
	}


				// initialize hash tables
	htable.initialize(10*size_FB);
	htable.set_key_function(&bigint_key);

	used_polys.initialize(10*size_FB);
	used_polys.set_key_function(&bigint_key);


				// initialize exponent vectors
	e_si.set_capacity(size_FB);
	e_si.set_size(size_FB);

	e_ra.set_capacity(size_FB);
	e_ra.set_size(size_FB);


	if (!use_old) {
		lp_table.initialize(100000);
		lp_table.set_key_function(&int_key);
	}


				// get large prime file names
	if (lpbound) {
		memset(faktor, '\0', (LINE_LEN)*sizeof(char));

		qo_get_name("PART");
		qo_get_name("TEMP");
		qo_get_name("ZUSA");

		lp_file.open(PART);
		if (lp_file.fail()) {
			std::remove(PART);
			std::remove(TEMP);
			std::remove(ZUSA);
			lidia_error_handler("qo_sieve", "init() - can't open PARTIALRELATIONS");
			return;
		}
	}
}



//
// init_self_initialization
//
// Task:
//      initializes the sieve for self_initialization
//

void
qo_sieve::init_self_initialization() {
	debug_handler("qo_sieve", "init_self_initialization");

	register lidia_size_t i, j;
	int p;
	lidia_size_t i1, i2;

	if (!is_init)
		return;

	is_si = true;

	// recompute size of A coefficients
	size_A = static_cast<double>(d_wurz)/M;

	if (!a_inv) {
		a_inv = new int[MAX_FB];
		BG = new bigint[MAX_FB];

		vorb = new int *[P_TOTAL];
		for (i = 0; i < static_cast<int>(P_TOTAL); ++i) {
			if (!(vorb[i] = new int[size_FB + 1])) {
				std::remove(PART);
				std::remove(TEMP);
				std::remove(ZUSA);
				lidia_error_handler("qo_sieve", "init_self_initialization:can't "
						    "allocate vorb");
				return;
			}
		}
	}



	// the size of coefficient A should be approximately
	// sqrt(kN)/M, so the size of the primes p dividing
	// A should be approximately (sqrt(kN/M))^(1/P_ONCE)

	double T = d_wurz / M;
	T = std::pow(2.0, LiDIA::log2(T) / P_ONCE);

	if (T > FB[size_FB-1]) {
		std::remove(PART);
		std::remove(TEMP);
		std::remove(ZUSA);

		lidia_error_handler ("qo_sieve", "init_self_initialization: P_ONCE too "
				     "small");
		return;
	}

	i = 0;
	while (FB[i] < T)
		++i;
	start_fb = static_cast<int>(i); // P_TOTAL consecutive primes p[start_fb], ...,
				// are chosen from the factor basis

	start_fb -= (P_ONCE >> 1);
	if (start_fb < 0)  start_fb = 0;

	j = start_fb+1;
	for (i = 0; i < static_cast<int>(P_TOTAL); ++i) {
		p = FB[j];
				// don't take divisors of the discriminant!
		while ((D_primes.size() > 0) &&
		       ((D_primes.bin_search(p, i1)) || (D_primes.bin_search(-p, i2))))
			p = FB[++j];
		Q_prime_glob[i] = j; // collect prime numbers which
				// will build the A coefficents
		++j;
	}

				// which product of small primes we are currently using
	bin_index = (1 << P_ONCE) - 2;
	max_bin = 0;
	for (i = 1; i <= P_ONCE; ++i)
		max_bin += (1 << (P_TOTAL-i));

				// first B value will have index 0
	index_i = -1;
}


//
// reset
//
// Task:
//      resets the sieve to pre-initialization state
//

void
qo_sieve::reset() {
	debug_handler("qo_sieve", "reset");

	register lidia_size_t i;

	if (is_init) {
		// reset sieving arrays
		if (a_inv) {
			delete [] a_inv;
			delete [] BG;
			a_inv = static_cast<int *>(NULL);
			BG = static_cast<bigint *>(NULL);

			for (i = 0; i < static_cast<int>(P_TOTAL); ++i)
				if (vorb[i])
					delete[] vorb[i];
			delete[] vorb;
			vorb = static_cast<int **>(NULL);
		}

		// reset hash tables
		htable.empty();
		used_polys.empty();
		lp_table.empty();

		// reset sieve polynomials
		F_si.assign_zero();
		F_ra.assign_zero();

		// reset exponent vectors
		e_si.set_size(0);
		e_ra.set_size(0);

		Q_prime_glob.reset();
		Q_prime_si.reset();
		Q_prime_ra.reset();

		// reset parameters
		Delta.assign_zero();
		D_primes.reset();

		size_FB = M = P_ONCE = P_TOTAL = smallstart = POLY = 0;
		Tval = d_wurz = size_A = 0.0;
		pBound = BachBound = 0;

		// reset si-parameters
		start_fb = bin_index = max_bin = index_i = 0;

		// reset relations arrays
		relation_array.reset();
		minima_array.set_size(0);

		// delete large prime files
		if (lpbound) {
			lpbound = 0;

			lp_file.close();

			std::remove(PART);
			std::remove(TEMP);
			std::remove(ZUSA);
		}
	}

	is_init = is_si = false;
}



//
// restart
//
// Task:
//      resets the sieve to the state immediately after initialization
//

void
qo_sieve::restart() {
	debug_handler("qo_sieve", "restart");

	register lidia_size_t i;

	relation_array.reset();
	relation_array.set_no_of_rows(size_FB);
	minima_array.set_size(0);

				// reset hash tables
	htable.empty();
	used_polys.empty();

	if (lpbound) {
		lp_file.close();

		std::remove(PART);
		std::remove(TEMP);
		std::remove(ZUSA);
	}

	if (a_inv) {
		delete [] a_inv;
		delete [] BG;
		a_inv = static_cast<int *>(NULL);
		BG = static_cast<bigint *>(NULL);

		for (i = 0; i < static_cast<int>(P_TOTAL); ++i)
			if (vorb[i])
				delete[] vorb[i];
		delete[] vorb;
		vorb = static_cast<int **>(NULL);
	}

	is_si = false;
}




//
// is_initialized
//
// Task:
//      returns true if the sieve has been initialized
//

bool
qo_sieve::is_initialized() {
	debug_handler("qo_sieve", "is_initialized");

	return is_init;
}



//
// is_si_initialized
//
// Task:
//      returns true if the sieve has been initialized for self-initialization
//

bool
qo_sieve::is_si_initialized() {
	debug_handler("qo_sieve", "is_si_initialized");

	return is_si;
}



//
// no_self_init
//
// Task:
//      turns off self-initialization
//

void
qo_sieve::no_self_init() {
	debug_handler("qo_sieve", "no_self_init");

	if (is_si)
		bin_index = max_bin;

	pBound = size_FB;
	size_A *= 2;
}



//
// discriminant
//
// Task:
//	returns the discriminant of the sieveing polynomials
//

bigint
qo_sieve::discriminant_used() {
	debug_handler("qo_sieve", "discriminant_used");

	return Delta;
}



//
// FB_size
//
// Task:
//	returns the size of the factor base used
//

int
qo_sieve::FB_size() {
	debug_handler("qo_sieve", "FB_size");

	return size_FB;
}



//
// dense_bound
//
// Task:
//	returns the index of the largest prime used for non-self-initialization
//

void
qo_sieve::dense_bound(int & val) {
	debug_handler("qo_sieve", "dense_bound");

	if (!is_si)
		val = pBound;
}




//
// next_polynomial
//
// Task:
//      generates the next polynomial to be used for sieving, either with
//      self-initialization or randomly.  In both cases, the roots of the
//      polynomial modulo the factor base primes are also computed.
//

void
qo_sieve::next_polynomial() {
	debug_handler("qo_sieve", "next_polynomial");

	if (is_si && (bin_index < max_bin)) {
				// generate poly using self-initialization

		if (index_i == static_cast<int>(POLY)-1)
			index_i = 0;
		else
			++index_i;

		new_polynomial_si();

		next_is_si = true;
	}
	else {
		bigint poly_pair;

				// generate randomly
		new_polynomial(-1, pBound);

				// process polynomial
		compute_roots();

		shift_left(poly_pair, F_ra.get_a(), 1);
		add(poly_pair, poly_pair, F_ra.get_b());
		used_polys.hash(poly_pair);
		next_is_si = false;
	}

}



//
// next_polynomial(lidia_size_t)
//
// Task:
//      generates the next polynomial to be used for sieving, forcing its
//      exponent vector to have a non-zero entry in the given position.
//      These polynomials are always generated randomly, even when
//      self-initialization is being used.  The roots of the
//      polynomial modulo the factor base primes are also computed.
//

void
qo_sieve::next_polynomial(lidia_size_t s) {
	debug_handler("qo_sieve", "next_polynomial(lidia_size_t)");

	bigint poly_pair;

				// generate randomly
	new_polynomial(s, 0);

				// process polynomial
	compute_roots();

	shift_left(poly_pair, F_ra.get_a(), 1);
	add(poly_pair, poly_pair, F_ra.get_b());
	used_polys.hash(poly_pair);
	next_is_si = false;
}



//
// next_polynomial(qi_class & A)
//
// Task:
//      generates the next polynomial to be used for sieving.  The ideal
//      used is A multiplied by some random power product of factor base
//      elements.  These polynomials are always generated randomly, even when
//      self-initialization is being used.  The roots of the
//      polynomial modulo the factor base primes are also computed.
//

void
qo_sieve::next_polynomial(const qi_class & A) {
	debug_handler("qo_sieve", "next_polynomial(qi_class)");

	bigint poly_pair;
	static qi_class last_A;
	static int numrho = 0;

				// generate randomly
	if (Delta.is_gt_zero() && last_A.is_equal(A) && (numrho < 100)) {
		F_ra.rho_indef();
		++numrho;

		shift_left(poly_pair, F_ra.get_a(), 1);
		add(poly_pair, poly_pair, F_ra.get_b());
		if (used_polys.search(poly_pair)) {
			new_polynomial(-1, 0);
			compose_regular(F_ra, F_ra, quadratic_form(inverse(A)));
			F_ra.reduce_regular();
			numrho = 0;
			shift_left(poly_pair, F_ra.get_a(), 1);
			add(poly_pair, poly_pair, F_ra.get_b());
		}
	}
	else {
		new_polynomial(-1, 0);
		compose_regular(F_ra, F_ra, quadratic_form(inverse(A)));
		F_ra.reduce_regular();
		last_A = A;
		numrho = 0;
		shift_left(poly_pair, F_ra.get_a(), 1);
		add(poly_pair, poly_pair, F_ra.get_b());
	}

				// process polynomial
	compute_roots();

	used_polys.hash(poly_pair);
	next_is_si = false;
}




//
// sieve
//
// Task:
//	sieves the current polynomial and appends any relations found to the
//      storage matrix.
//

void
qo_sieve::sieve() {
	debug_handler("qo_sieve", "sieve");

				// sieve using one polynomial
	qo_sieve_interval();

				// output relations (all of them)
	test_candidates(-1, -1, true);
}




//
// sieve_one(lidia_size_t)
//
// Task:
//	sieves the current polynomial and returns the first relation found
//      in which the coefficient of the given index is non-zero
//

bool
qo_sieve::sieve_one(lidia_size_t s, math_vector< lidia_size_t > & v,
		    quadratic_number_standard & alpha) {
	debug_handler("qo_sieve", "sieve_one(lidia_size_t)");

	bool found;

				// sieve using one polynomial
	qo_sieve_interval();

				// output relations (only one)
	found = test_candidates(1, s, true);

	if (found) {
		lidia_size_t ra_cols = relation_array.get_no_of_columns();
		--ra_cols;
		relation_array.get_column_vector(v, ra_cols);
		alpha.assign(minima_array[ra_cols]);

		relation_array.set_no_of_columns(ra_cols);
		minima_array.set_size(ra_cols);

		if (ra_cols == 0) {
			relation_array.reset();
			relation_array.set_no_of_rows(size_FB);
		}
	}

	return found;
}




//
// sieve_one(qi_class)
//
// Task:
//	sieves the current polynomial and returns the first representation of
//      A over FB that is found.
//

bool
qo_sieve::sieve_one(const qi_class & A,
		    math_vector< lidia_size_t > & v,
		    quadratic_number_standard & alpha) {
	debug_handler("qo_sieve", "sieve_one(qi_class, vec, qn)");

	bool found;

				// sieve using one polynomial
	qo_sieve_interval();

				// output relations
	found = test_single(v, alpha);

	return found;
}




//
// sieve_one(qi_class)
//
// Task:
//	sieves the current polynomial and returns true if A can be represented
//      over FB.
//

bool
qo_sieve::sieve_one(const qi_class & A) {
	debug_handler("qo_sieve", "sieve_one(qi_class)");

	bool found;
	int p;

	if (A.get_a().intify(p) == 0)
		if (lp_table.search(p))
			return true;

				// sieve using one polynomial
	qo_sieve_interval();

				// output relations
	found = test_candidates(1, -1, false);

	return found;
}




//
// count_lp_relations()
//
// Task:
//      counts the number of large prime relations that can be produced from
//      the current set of partial relations
//

lidia_size_t
qo_sieve::count_lp_relations() {
	debug_handler("qo_sieve", "count_lp_relations");

	FILE *fp;
	char line[LINE_LEN];

	if (!lpbound)
		return 0;

	lp_file.close();

	if ((fp = fopen(PART, "r")) != NULL) {
		fclose(fp);
		rename(PART, TEMP);
		concat(TEMP, ZUSA, PART);

		std::remove(ZUSA);
		std::remove(TEMP);
		sort_file_line(PART, ZUSA);
		std::remove(PART);
	}


	fp = fopen (ZUSA, "r");

	long counter = 0, large, large_old = -1, i = 1;

	while (fgets (line, LINE_LEN, fp)) {
		if (strlen(line) > 20) {
			large = atol(line);

			if (large == large_old) {
				counter += i;
				++i;
			}
			else {
				large_old = large;
				i = 1;
			}
		}
	}
	fclose(fp);

	lp_file.open(PART);
	if (lp_file.fail()) {
		std::remove(PART);
		std::remove(TEMP);
		std::remove(ZUSA);
		lidia_error_handler("qo_sieve", "count_lp_relations() - "
				    "can't open PARTIALRELATIONS");
		return (-1);
	}

	return (counter);
}




//
// get_lp_relations()
//
// Task:
//      generates all the large prime relations that can be produced from
//      the current set of partial relations and stores them in the storage
//      matrix
//

lidia_size_t
qo_sieve::get_lp_relations() {
	debug_handler("qo_sieve", "get_lp_relations");

	lidia_size_t i, j, pos, ra_cols;
	int small_value, large_old, zeros;
	partial_relation partial, old_partial, new_partial;
	base_vector< partial_relation > rvec;
	std::ifstream in;
	std::ofstream out;

	if (!lpbound)
		return 0;

	lp_file.close();

				// sort the file of partial relations
	rename(ZUSA, TEMP);

				// combine partials into full relations
	in.open(TEMP);
	if (in.fail()) {
		std::remove(PART);
		std::remove(TEMP);
		std::remove(ZUSA);
		lidia_error_handler("qo_sieve", "get_lp_relations() - can't "
				    "open PARTIALRELATIONS");
		return(-1);
	}
	out.open(ZUSA);

	ra_cols = relation_array.get_no_of_columns();
	rvec.set_mode(EXPAND);
	small_value = zeros = 0;
	large_old = -1;

	in >> rvec[0];
	large_old = rvec[0].get_large();
	if (!lp_table.search(large_old))
		lp_table.hash(large_old);

	while (!in.eof()) {
		in >> partial;

		if (!partial.is_zero() && (partial.get_large() == large_old)) {
			pos = 0;
			while (partial.get_large() == large_old && !in.eof()) {
				rvec[++pos] = partial;
				in >> partial;
			}
			++pos;

			for (i = 0; i < pos-1; ++i) {
				for (j = i+1; j < pos; ++j) {
					if (!make_full(new_partial, rvec[i], rvec[j], Delta, size_FB)) {
						std::remove(PART);
						std::remove(TEMP);
						std::remove(ZUSA);
						lidia_error_handler("qo_sieve", "get_lp_relations() - "
								    "error combining partial relations");
						return(-1);
					}

					if (!new_partial.is_zero()) {
						htable.hash(new_partial.get_index());
						relation_array.set_no_of_columns(ra_cols+1);
						new_partial.add_to_matrix(relation_array, minima_array, ra_cols);

						++ra_cols;
						++small_value;
					}
					else
						++zeros;
				}
			}
		}
		else
			out << rvec[0] << std::flush;

		rvec[0] = partial;
		large_old = rvec[0].get_large();

		if (!lp_table.search(large_old))
			lp_table.hash(large_old);
	}

	in.close();
	out.close();

	std::remove(TEMP);

	lp_file.open(PART);
	if (lp_file.fail()) {
		std::remove(PART);
		std::remove(ZUSA);
		lidia_error_handler("qo_sieve", "init() - can't open PARTIALRELATIONS");
		return (-1);
	}

	return small_value;
}



//
// get_relation()
//
// Task:
//      returns one relation.  If no relations are currently in the storage
//      matrix, false is returned.
//

bool
qo_sieve::get_relation(math_vector< lidia_size_t > & v, quadratic_number_standard &alpha) {
	debug_handler("qo_sieve", "get_relation");

	lidia_size_t ra_cols = relation_array.get_no_of_columns();

	if (ra_cols > 0) {
		--ra_cols;
		relation_array.get_column_vector(v, ra_cols);
		alpha = minima_array[ra_cols];
		relation_array.set_no_of_columns(ra_cols);
		minima_array.set_size(ra_cols);
		if (ra_cols == 0) {
			relation_array.reset();
			relation_array.set_no_of_rows(size_FB);
		}

		return true;
	}
	else
		return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
