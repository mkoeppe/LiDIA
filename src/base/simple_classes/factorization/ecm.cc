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
//	Author	: Andreas Mueller (AM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================

// Description : see the diploma thesis of the first author


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const unsigned int rational_factorization::ecm_params[29][3] = {
	{151, 1009, 3},
	{211, 2003, 4},
	{283, 3203, 5},
	{349, 4297, 5},
	{411, 8861, 6},
	{659, 17981, 8},
	{997, 31489, 10},
	{1439, 54617, 13},
	{2111, 89501, 17},
	{3067, 144961, 21},
	{4391, 245719, 27},
	{6287, 382919, 33},
	{8839, 577721, 42},
	{12241, 845951, 54},
	{16339, 1285633, 68},
	{21977, 1788863, 88},
	{29860, 2508269, 111},
	{41011, 3550163, 137},
	{53214, 5139913, 172},
	{71537, 6996393, 215},
	{97499, 9650917, 261},
	{125653, 13767227, 325},
	{168273, 18372001, 398},
	{223559, 24108067, 489},
	{280891, 33036147, 612},
	{376513, 43533099, 736},
	{489249, 54671451, 909},
	{631461, 77351879, 1088},
	{830049, 97901685, 1315}
};



unsigned int rational_factorization::
ecm_read_max(int stell)               // returns the limit for prime numbers
{                                     // to execute the given strategy
	return (ecm_params[stell-6][1]);
}



int rational_factorization::
ecm_job_planing(int strat[30], int buf[30][5])
{                                     // read the parameters for ECM for
	// a given strategy
	int i = 0, j = strat[0];
	unsigned int job = 1;

	if (info) {
		std::cout << "Strategy for ECM : \n\n";
		std::cout << "   #    |             |     prime limits    | #curves \n";
		std::cout << " digits |   curves    |  step 1  |  cont    |  total \n";
		std::cout << "-------------------------------------------------------\n";
		std::cout.flush();
	}

	while (j) {
		j -= 6;

		buf[i][0] = ecm_params[j][0];
		buf[i][1] = buf[i][0];
		buf[i][2] = ecm_params[j][1];
		buf[i][3] = ecm_params[j][2];
		buf[i][4] = buf[i][3];

		if (info)  printf("   %2d   | %4d -", strat[i], job);

		job += buf[i][3];

		if (info)
			printf(" %4d |%9d | %8d | %4d\n", job-1, buf[i][0], buf[i][2], buf[i][3]);

		j = strat[++i];
	}

	if (info) {
		fflush(stdout);
		std::cout << "========================================================\n";
		std::cout.flush();
	}

	return(i);
}



void mecfind(ec_curve & me, ec_point_M & mp, bigint & factor, unsigned long B, unsigned long C, ecm_primes & prim)
{
	// first step of ECM:
	// compute  mp  <--  K * mp
	// K is the product of all primes powers up to B
	unsigned long p;

	bigmod a_2, m1, inv_4(4L);

	long k, count = 0;
	unsigned long sqrt_C = static_cast<long>(std::sqrt (static_cast<double>(C)));

	inv_4.invert(0); // precompute the constant a_2 = (me.a+2)/4
	add(m1, me.A(), 2L); // for the given curve me
	multiply(a_2, inv_4, m1);

	for (k = 0; k < 31; k++)
		multiply_by_2(mp, mp, a_2); // mp = 2^30 * mp

                             //  successive computation of    mp = p^k * mp
                             //  for all primes p = 3, ..., B   and
                             //  k integer with p^k <= C && p^(k+1) > C
	while ((p = prim.getprimes()) <= sqrt_C) {
		// loop: exponent k > 1
		k = static_cast<long>(std::log(static_cast<double>(C))/std::log(static_cast<double>(p)));
		multiply(mp, mp, a_2, p, power(p, k));
	}

	do {
		// loop: exponent k = 1
		multiply(mp, mp, a_2, p, p);

		if (count++ > 1000) {
			m1.assign(mp.Z());
			factor = m1.invert(1);
			if (!factor.is_one())
				return;

			count = 0;
			multiply(m1, m1, mp.X());
			mp.assign(m1, bigmod(1L));
		}
	} while ((p = prim.getprimes()) < B);

	m1.assign(mp.Z());
	factor = m1.invert(1);
	if (!factor.is_one())
		return;

	multiply(m1, m1, mp.X());
	mp.assign(m1, bigmod(1L));
}



void cont(ec_point_W & wp, bigint & factor, unsigned long B, unsigned long D, ecm_primes & prim)
{
	register long w, u, v, count;

	unsigned long p;

	bigmod * x, h, c(1L);

	ec_point_W wq;

	// compute the number of babysteps w
	// and an initial value v for then giantsteps
	u = static_cast<long>(std::ceil(std::sqrt(static_cast<double>(D))));

	w = (u / 30) * 30;
	if (w + 30 - u < u - w)
		w += 30;
	D = w * w;
	v = static_cast<long>(B) / w;

	x = new bigmod[w+2];

	// precompute babysteps:
	// k * wp for k = 1, .., w

	for (count = 1; count <= w; count += 2) {
		add(wq, wp, wq, factor);
		if (!factor.is_one()) {
			delete [] x;
			return;
		}

		if (gcd(w, count) == 1)
			x[count] = wq.X();

		add(wq, wp, wq, factor);
		if (!factor.is_one()) {
			delete [] x;
			return;
		}
	}

	// computed: wq = w * wp
	// initialization for giantsteps:
	//     wp = v * wq = v * w * wp

	multiply(wp, wq, v, factor);
	if (!factor.is_one()) {
		delete [] x;
		return;
	}

	count = 1;
	p = prim.getprimes();
	// giantsteps
	while (p < D && p > 1) {
		if ((u = v*w - p) < 1) {
			v++;
			add(wp, wp, wq, factor);
			if (!factor.is_one()) {
				delete [] x;
				return;
			}

			u = v*w - p;
			if (u < 0)
				u += w;
		}

		subtract(h, wp.X(), x[u]);

		multiply(c, c, h);

		if (count++ >= 500) {
			count = 1;
			factor = c.invert(1);
			if (!factor.is_one()) {
				delete [] x;
				return;
			}
		}
		p = prim.getprimes();
	}

	factor = c.invert(1);
	delete [] x;
}



void rational_factorization::
ecm(lidia_size_t index, int ecm_restlength, int jobs_avail, int job_buffer[30][5], ecm_primes & prim)
{
	bigint N (factors[index].single_base);
	int exp = factors[index].single_exponent;
	lidia_size_t n = no_of_comp();

	unsigned int B, C, D;
	int uc, current_length, count = 0, m, job_nr = 0, job_grp = 0, qs_2;

	char *n_string = new char[(N.bit_length() / 3) + 10];

	random_generator rg;

	current_length = bigint_to_string(N, n_string);
	uc = ((current_length-30) / 3) << 1;

	if (current_length <= ecm_restlength && uc <= 0) {
		if (info) {
			std::cout << "\nUsing MPQS for further factorization ....\n";
			std::cout.flush();
		}

		delete [] n_string;
		return;
	}

	if ((factors[index].factor_state >> 2) > 0) {
		int params_for_computed = ecm_params[(factors[index].factor_state >> 2)-6][0];

		while (job_buffer[job_grp][0] < params_for_computed)
			job_grp ++;

		jobs_avail -= job_grp;
	}

	B = job_buffer[job_grp][0]; // first job
	C = job_buffer[job_grp][1];
	D = job_buffer[job_grp][2];

	bigint old_modul = bigmod::modulus();
	bigmod::set_modulus(N);

	rf_single_factor fact;
	bigint factor, Q, R;
	factor.assign_one();

	ec_curve ell_curve;
	ec_point_M mp;
	ec_point_W wp;

	timer ti;
	qs_2 = static_cast<int>(zeitqs(current_length)) >> 1;

	ti.start_timer();

	while (!N.is_one() && jobs_avail > 0) {
		while (jobs_avail > 0) {
			prim.resetprimes(1);

			if (job_buffer[job_grp][3] < 1) {
				if (--jobs_avail <= 0)
					break;

				job_grp++;
				B = job_buffer[job_grp][0]; // next job
				C = job_buffer[job_grp][1];
				D = job_buffer[job_grp][2];
			}

			job_buffer[job_grp][3]--;

			if (info) {
				std::cout << "\n" << ++job_nr << ".curve ";
				std::cout.flush();
			}

			rg >> m;

			mecgen16(ell_curve, mp, factor, m);
			if (!factor.is_one())                    break;

			mecfind(ell_curve, mp, factor, B, C, prim);
			if (!factor.is_one())                    break;

			factor = trans(ell_curve, wp, mp);
			if (!factor.is_one())                    break;

			cont(wp, factor, B, D, prim);
			if (!factor.is_one())                    break;

			ti.stop_timer();
			ti.cont_timer();

			if ((current_length <= ecm_restlength && --uc <= 0) ||
			    ((ti.user_time()+ti.sys_time())/100 > qs_2)) {
				fact.factor_state = rf_single_factor::not_prime;

				fact.single_base = N;
				fact.single_exponent = exp;
				factors[index] = fact;

				delete [] n_string;
				if (!old_modul.is_zero()) bigmod::set_modulus(old_modul);
				return;
			}
		}

		if (jobs_avail > 0) {
			if (!(factor.is_one() || factor.is_zero() || factor == N)) {
				div_rem(Q, R, N, factor);

				while (R.is_zero()) {
					count++;
					N.assign(Q);
					div_rem(Q, R, N, factor);
				}

				if (count) {
					if (info) {
						if (count > 1)
							std::cout << "\nfactor " << factor << " ^ " << count;
						else
							std::cout << "\nfactor " << factor;
						std::cout.flush();
					}

					if (is_prime(factor, 8))
						fact.factor_state = rf_single_factor::prime;
					else
						fact.factor_state = rf_single_factor::not_prime;

					fact.single_base = factor;
					fact.single_exponent = count * exp;

					count = 0; // ts

					if (N.is_one()) {
						factors[index] = fact;
						break;
					}
					else
						factors[n++] = fact;

					if (is_prime(N, 8)) {
						fact.factor_state = rf_single_factor::prime;

						fact.single_base = N;
						fact.single_exponent = exp;
						factors[index] = fact;

						if (info) {
							std::cout << "\nprime number " << N << "\n";
							std::cout.flush();
						}

						N.assign_one();
						break;
					}

					current_length = bigint_to_string(N, n_string);
					uc = ((current_length-30) / 3) << 1;

					qs_2 = static_cast<int>(zeitqs(current_length)) >> 1;
					if ((current_length <= ecm_restlength && uc <= 0) ||
					    ((ti.user_time()+ti.sys_time())/100 > qs_2)) {
						fact.factor_state = rf_single_factor::not_prime;
						fact.single_base = N;
						fact.single_exponent = exp;
						factors[index] = fact;

						delete [] n_string;
						if (!old_modul.is_zero()) bigmod::set_modulus(old_modul);
						return;
					}

					bigmod::set_modulus(N);

				}
			}
		}
	}


	if (!N.is_one()) {
		int i = 0;

		unsigned int table = job_buffer[job_grp][0];

		while (ecm_params[i][0] != table)
			i++;

		i += 6;

		fact.single_base = N;
		fact.single_exponent = exp;
		fact.factor_state = (rf_single_factor::not_prime) | (i << 2);
		factors[index] = fact;
	}

	delete [] n_string;
	if (!old_modul.is_zero()) bigmod::set_modulus(old_modul);
}



void rational_factorization::
ecm(lidia_size_t index, int jobs_avail, int job_buffer[30][5], ecm_primes & prim)
{
	bigint N (factors[index].single_base);
	int exp = factors[index].single_exponent;
	lidia_size_t n = no_of_comp();

	unsigned int B, C, D;
	int count = 0, m, job_nr = 0, job_grp = 0;

	random_generator rg;

	if ((factors[index].factor_state >> 2) > 0) {
		int params_for_computed = ecm_params[(factors[index].factor_state >> 2)-6][0];

		while (job_buffer[job_grp][0] < params_for_computed)
			job_grp ++;

		jobs_avail -= job_grp;
	}

	B = job_buffer[job_grp][0]; // first job
	C = job_buffer[job_grp][1];
	D = job_buffer[job_grp][2];

	bigint old_modul = bigmod::modulus();
	bigmod::set_modulus(N);

	rf_single_factor fact;
	bigint factor, Q, R;
	factor.assign_one();

	ec_curve ell_curve;
	ec_point_M mp;
	ec_point_W wp;

	while (!N.is_one() && jobs_avail > 0) {
		while (jobs_avail > 0) {
			prim.resetprimes(1);

			if (info) {
				std::cout << "\n" << ++job_nr << ".curve ";
				std::cout.flush();
			}

			if (job_buffer[job_grp][3] <= 0) {
				if (--jobs_avail <= 0)
					break;

				job_grp++;
				B = job_buffer[job_grp][0]; // next job
				C = job_buffer[job_grp][1];
				D = job_buffer[job_grp][2];
			}
			else
				job_buffer[job_grp][3]--;

			rg >> m;

			mecgen16(ell_curve, mp, factor, m);
			if (!factor.is_one())		    break;

			mecfind(ell_curve, mp, factor, B, C, prim);
			if (!factor.is_one())		    break;

			factor = trans(ell_curve, wp, mp);
			if (!factor.is_one())               break;

			cont(wp, factor, B, D, prim);
			if (!factor.is_one())	            break;
		}

		if (jobs_avail > 0) {
			if (!(factor.is_one() || factor.is_zero() || factor == N)) {
				div_rem(Q, R, N, factor);

				while (R.is_zero()) {
					count++;
					N.assign(Q);
					div_rem(Q, R, N, factor);
				}

				if (count) {
					if (info) {
						if (count > 1)
							std::cout << "\nfactor " << factor << " ^ " << count;
						else
							std::cout << "\nfactor " << factor;
						std::cout.flush();
					}

					if (is_prime(factor, 8))
						fact.factor_state = rf_single_factor::prime;
					else
						fact.factor_state = rf_single_factor::not_prime;

					fact.single_base = factor;
					fact.single_exponent = count * exp;

					count = 0; // ts

					if (N.is_one()) {
						factors[index] = fact;
						break;
					}
					else
						factors[n++] = fact;

					if (is_prime(N, 8)) {
						fact.factor_state = rf_single_factor::prime;

						fact.single_base = N;
						fact.single_exponent = exp;
						factors[index] = fact;

						if (info) {
							std::cout << "\nprime number " << N << "\n";
							std::cout.flush();
						}

						N.assign(1);
						break;
					}

					count = 0;
					bigmod::set_modulus(N);
				}
			}
		}
	}

	if (!N.is_one()) {
		int i = 0;

		unsigned int table = job_buffer[job_grp][0];

		while (ecm_params[i][0] != table)
			i++;

		i += 6;

		fact.single_base = N;
		fact.single_exponent = exp;
		fact.factor_state = (rf_single_factor::not_prime) | (i << 2);
		factors[index] = fact;
	}

	if (!old_modul.is_zero()) bigmod::set_modulus(old_modul);
}



void rational_factorization::
ecm_comp(lidia_size_t index, int upper_bound, int lower_bound, int step)
{
	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "ecm_comp::index out of range");
		return;
	}

	if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime) {
		if (info)
			std::cout << "index " << index << " :  prime number " << factors[index].single_base << "\n";
		return;
	}

	if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::dont_know) {
		if (is_prime(factors[index].single_base, 8)) {
			if (info)
				std::cout << "index " << index << " :  prime number " << factors[index].single_base << "\n";

			factors[index].factor_state = rf_single_factor::prime;
			return;
		}
		else
			factors[index].factor_state = rf_single_factor::not_prime;
	}

	if (step <= 0) {
		lidia_warning_handler("rational_factorization", "ecm_comp::step <= 0");
		step = 3;
	}

	if (upper_bound > 34) {
		lidia_warning_handler("rational_factorization", "ecm_comp::upper_bound > 34");
		upper_bound = 34;
	}

	if (lower_bound < 6) {
		lidia_warning_handler("rational_factorization", "ecm_comp::lower_bound < 6");
		lower_bound = 6;
	}

	if (lower_bound > upper_bound) {
		lidia_warning_handler("rational_factorization", "ecm_comp::lower_bound > upper_bound");
		upper_bound = lower_bound;
	}

	int k, D, jobs_avail, job_buffer[30][5], strategy[30];
	int n = ((factors[index].single_base).bit_length() / 3) + 10;

	if (upper_bound == 34) {
		char *n_string;
		n_string = new char[n];
		upper_bound = (bigint_to_string(factors[index].single_base, n_string) >> 1) + 1;
		delete [] n_string;
		if (upper_bound > 34)  upper_bound = 34;
	}

	strategy[0] = lower_bound;
	n = lower_bound + step;
	k = 1;

	while (n < upper_bound) {
		strategy[k++] = n;
		n += step;
	}

	if (lower_bound < upper_bound)
		strategy[k++] = upper_bound;

	strategy[k] = 0;

	D = ecm_read_max(upper_bound);

	ecm_primes prim(1, static_cast<unsigned int>(D+200), 200000);

	jobs_avail = ecm_job_planing(strategy, job_buffer);

	ecm(index, jobs_avail, job_buffer, prim);

	compose();
}



void rational_factorization::
ecm(int upper_bound, int lower_bound, int step)
{
	if (is_prime_factorization())
		return;
	lidia_size_t index, len = no_of_comp();

	int D, k, jobs_avail, job_buffer[30][5], strategy[30];
	int n = ((factors[len-1].single_base).bit_length() / 3) + 10;

	if (step <= 0) {
		lidia_warning_handler("rational_factorization", "ecm::step < 0");
		step = 3;
	}

	if (upper_bound > 34) {
		lidia_warning_handler("rational_factorization", "ecm::upper_bound > 34");
		upper_bound = 34;
	}

	if (lower_bound < 6) {
		lidia_warning_handler("rational_factorization", "ecm::lower_bound < 6");
		lower_bound = 6;
	}

	if (lower_bound > upper_bound) {
		lidia_warning_handler("rational_factorization", "ecm::lower_bound > upper_bound");
		upper_bound = lower_bound;
	}

	if (upper_bound == 34) {
		char *n_string;
		n_string = new char[n];
		upper_bound = (bigint_to_string(factors[len-1].single_base, n_string) >> 1) + 1;
		delete [] n_string;
		if (upper_bound > 34)  upper_bound = 34;
		if (upper_bound < 6)   upper_bound = 6;
	}

	D = ecm_read_max(upper_bound);

	ecm_primes prim(1, D+200, 200000);

	strategy[0] = lower_bound; k = 1;
	n = lower_bound + step;

	while (n < upper_bound) {
		strategy[k++] = n;
		n += step;
	}

	if (lower_bound < upper_bound) {
		strategy[k] = upper_bound;
		strategy[k+1] = 0;
	}
	else strategy[k] = 0;

	if (info) {
		std::cout << "\n ECM with same strategy for each number";
		std::cout << "\n --------------------------------------\n\n";
	}

	jobs_avail = ecm_job_planing(strategy, job_buffer);

	for (index = 0; index < len; index++) {
		if (info) {
			std::cout << "\nindex " << index;
			std::cout << "\n--------\n";
		}

		if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime) {
			if (info)
				std::cout << "index " << index << " :  prime number " << factors[index].single_base << "\n";
			continue;
		}

		if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::dont_know) {
			if (is_prime(factors[index].single_base, 8)) {
				if (info)
					std::cout << "index " << index << " :  prime number " << factors[index].single_base << "\n";

				factors[index].factor_state = rf_single_factor::prime;

				continue;
			}
			else
				factors[index].factor_state = rf_single_factor::not_prime;
		}

		ecm(index, jobs_avail, job_buffer, prim);

		for (D = 0; D < jobs_avail; D++)
			job_buffer[D][3] = job_buffer[D][4];
	}

	compose();
}



rational_factorization ecm (const bigint &N, int upper_bound, int lower_bound, int step)
{
	rational_factorization f(N);

	f.ecm(upper_bound, lower_bound, step);

	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
