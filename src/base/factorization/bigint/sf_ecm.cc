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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/base/ecm_primes.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const unsigned int single_factor< bigint >::ecm_params[29][3] = {
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



unsigned int single_factor< bigint >::
ecm_read_max(int stell) const         // returns the limit for prime numbers
{                                     // to execute the given strategy
	return (ecm_params[stell-6][1]);
}



int single_factor< bigint >::
ecm_job_planing(int strat[30], int buf[30][5])
{                                     // read the parameters for ecm for
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



void single_factor< bigint >::ECM(factorization< bigint > & tmp,
				   int jobs_avail,
				   int job_buffer[30][5],
				   ecm_primes & prim, bool jump_QS)
{
	unsigned int B, C, D;
	int count, m, job_nr = 0, job_grp = 0;

	B = job_buffer[job_grp][0];
	C = job_buffer[job_grp][1];
	D = job_buffer[job_grp][2];

	bigint old_modul = bigmod::modulus();
	bigmod::set_modulus(rep);

	bigint factor(1), Q, R;
	ec_curve ell_curve;
	ec_point_M mp;
	ec_point_W wp;
	random_generator rg;
	timer ti;
	int qs_2;

	// the following statements are only used in combination
	// ECM and QS --> determine when to jump to QS

	qs_2 = static_cast<int>(zeitqs(decimal_length(rep))) >> 1;
	ti.start_timer();

	// start the next ECM job

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
			std::cout << "\ncurve " << ++job_nr;
			std::cout.flush();
		}

		rg >> m;

		mecgen16(ell_curve, mp, factor, m);
		if (!factor.is_one())
			break;

		mecfind(ell_curve, mp, factor, B, C, prim);
		if (!factor.is_one())
			break;

		factor = trans(ell_curve, wp, mp);
		if (!factor.is_one())
			break;

		cont(wp, factor, B, D, prim);
		if (!factor.is_one())
			break;

		// check whether we should jump to QS

		ti.stop_timer();
		ti.cont_timer();
		if ((ti.user_time()+ti.sys_time())/100 > qs_2 && jump_QS) {
			if (info)
				std::cout << "\n\nUsing MPQS for further factorization\n" << std::flush;

			if (decimal_length(rep) > 65)
				prim.initprimes(2, 500000, 100000);
			else
				prim.initprimes(2, 200000, 100000);

			MPQS(tmp, prim);
			return;
		}
	}

	if (!(factor.is_one() || factor.is_zero() || factor == rep))  // factor found
	{
		count = 0;
		div_rem(Q, R, rep, factor);

		while (R.is_zero()) {
			count++;
			rep.assign(Q);
			div_rem(Q, R, rep, factor);
		}

		if (count) {
			if (info) {
				if (count > 1)
					std::cout << "\nfactor: " << factor << " ^ " << count;
				else
					std::cout << "\nfactor: " << factor;
				std::cout.flush();
			}

			single_factor< bigint > fact(factor);

			if (is_prime(factor, 8))
				fact.set_prime_flag(prime);
			else
				fact.set_prime_flag(not_prime);

			tmp.append(fact, count);
			tmp.append(rep);
			rep.assign_one();
			return;
		}
	}

	if (!old_modul.is_zero())
		bigmod::set_modulus(old_modul);
}



factorization< bigint > single_factor< bigint >::
ECM(int upper_bound, int lower_bound, int step, bool jump_to_QS)
{
	factorization< bigint > h;

	if (lower_bound < 6) {
		lidia_warning_handler("single_factor< bigint >",
				      "ecm::lower_bound < 6");
		lower_bound = 6;
	}

	if (upper_bound > 34) {
		lidia_warning_handler("single_factor< bigint >",
				      "ecm::upper_bound > 34");
		upper_bound = 34;
	}

	if (lower_bound > upper_bound) {
		lidia_warning_handler("single_factor< bigint >",
				      "ecm::lower_bound > upper_bound");
		upper_bound = lower_bound;
	}

	if (step <= 0) {
		lidia_warning_handler("single_factor< bigint >",
				      "ecm::step <= 0");
		step = 3;
	}

	if (rep.is_negative()) {
		h.append(single_factor< bigint > (-1));
		rep.negate();
	}

	if (is_prime_factor(1)) {
		h.append(*this);
		return h;
	}

	int k, D = 0, jobs_avail = 0, job_buffer[30][5], strategy[30];
	int n = (rep.bit_length() / 3) + 10;

	if (upper_bound == 34) {
		upper_bound = (decimal_length (rep) >> 1) + 1;
		if (upper_bound > 34)
			upper_bound = 34;
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
	jobs_avail = ecm_job_planing(strategy, job_buffer);

	ecm_primes prim(1, static_cast<unsigned int>(D+200), 200000);

	TrialDiv(h, 2, D+200, prim);
	if (rep.is_one())
		return h;

	if (is_prime_factor(1)) {
		h.append(*this);
		rep.assign_one();
		return h;
	}

	ECM(h, jobs_avail, job_buffer, prim, jump_to_QS);

	if (!rep.is_one())
		h.append(*this);

	return h;
}



factorization< bigint > ECM(const bigint & x,
			    int upper_bound,
			    int lower_bound,
			    int step)
{
	single_factor< bigint > a(x);
	factorization< bigint > f;
	f = a.ECM(upper_bound, lower_bound, step);

	if (!a.is_one())
		f.append(a);
	return f;
}

factorization< bigint > ECM(const bigint & x,
			    int upper_bound,
			    int lower_bound)
{
	single_factor< bigint > a(x);
	factorization< bigint > f;
	f = a.ECM(upper_bound, lower_bound);

	if (!a.is_one())
		f.append(a);
	return f;
}

factorization< bigint > ECM(const bigint & x,
			    int upper_bound)
{
	single_factor< bigint > a(x);
	factorization< bigint > f;
	f = a.ECM(upper_bound);

	if (!a.is_one())
		f.append(a);
	return f;
}

factorization< bigint > ECM(const bigint & x)
{
	single_factor< bigint > a(x);
	factorization< bigint > f;
	f = a.ECM();

	if (!a.is_one())
		f.append(a);
	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
