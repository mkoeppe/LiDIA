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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/alg_number.h"
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef LIDIA_IMPLICIT_CAST_EXPLICIT
#define alg_ideal_cast(O) alg_ideal(O)
#else
#define alg_ideal_cast(O) O
#endif

bool single_factor< alg_ideal >::verbose_flag = false;


//********************************************************************
//                              class single_factor< alg_ideal >
//********************************************************************

single_factor< alg_ideal >::single_factor()
{
	debug_handler("single_factor< alg_ideal >", "single_factor< alg_ideal > ()");

	// DEFAULT VALUE MUST BE '1', i.e. the identity element of multiplication !!!
	// we have a problem here : we must assign '1' to an ideal without
	// really knowing the base the ideal is defined over.
	// However, this is what our default nf_base is for !!!
	// Note that this differs from what the polynomials are doing quite a bit!

	set_prime_flag(decomposable_object::not_prime);
	rep_a.assign_one();
}



single_factor< alg_ideal >::
single_factor(const single_factor< alg_ideal > & x):
	rep_p(x.rep_p), rep_a(x.rep_a)


{
	debug_handler("single_factor< alg_ideal >",
		      "single_factor< alg_ideal > (single_factor< alg_ideal > &)");
	set_prime_flag(x.prime_flag());
}



single_factor< alg_ideal >::
single_factor(const alg_ideal & x):
	rep_a(x)
	// everything else is initialized to zero by default !!
{
	debug_handler("single_factor< alg_ideal >",
		      "single_factor< alg_ideal > (alg_ideal&)");
	set_prime_flag(decomposable_object::unknown);
}



single_factor< alg_ideal >::
single_factor(const prime_ideal & x):
	rep_p(x)
{
	debug_handler("single_factor< alg_ideal >",
		      "single_factor< alg_ideal > (alg_ideal&)");
	set_prime_flag(decomposable_object::prime);
}



void swap(single_factor< alg_ideal > & a, single_factor< alg_ideal > &b)
{
	debug_handler("single_factor< alg_ideal >",
		      "swap (single_factor< alg_ideal > &, single_factor< alg_ideal > &)");

	swap(a.rep_p, b.rep_p);
	swap(a.rep_a, b.rep_a);

	single_factor< alg_ideal >::decomp_state tmp = a.prime_flag();
	a.set_prime_flag(b.prime_flag());
	b.set_prime_flag(tmp);
}



single_factor< alg_ideal > &
single_factor< alg_ideal >::
operator = (const single_factor< alg_ideal > & x)
{
	debug_handler("single_factor< alg_ideal >",
		      "operator = (single_factor< alg_ideal > &)");

	assign(x);
	return *this;
}



single_factor< alg_ideal > &
single_factor< alg_ideal >::
operator = (const alg_ideal & x)
{
	debug_handler("single_factor< alg_ideal >", "operator = (alg_ideal&)");

	assign(x);
	return *this;
}



single_factor< alg_ideal > &
single_factor< alg_ideal >::
operator = (const prime_ideal & x)
{
	debug_handler("single_factor< alg_ideal >", "operator = (alg_ideal&)");

	assign(x);
	return *this;
}



void
single_factor< alg_ideal >::
assign(const single_factor< alg_ideal > & x)
{
	debug_handler("single_factor< alg_ideal >",
		      "assign(single_factor< alg_ideal > &)");

	set_prime_flag(x.prime_flag());
	rep_p = x.rep_p;
	rep_a.assign(x.rep_a);
}



void
single_factor< alg_ideal >::
assign(const alg_ideal & x)
{
	debug_handler("single_factor< alg_ideal >", "assign(alg_ideal&)");

	rep_a.assign(x);
	set_prime_flag(decomposable_object::unknown);
}



void
single_factor< alg_ideal >::
assign(const prime_ideal & x)
{
	debug_handler("single_factor< alg_ideal >", "assign(alg_ideal&)");

	rep_p = x;
	set_prime_flag(decomposable_object::prime);
}



alg_ideal
single_factor< alg_ideal >::
extract_unit()
{
	debug_handler("single_factor< alg_ideal >", "extract_unit()");
	if (is_prime_factor())
		return (order(rep_p.second_generator().which_base())).operator alg_ideal();
	else
		return (order(rep_a.which_base())).operator alg_ideal();
}



//
// arithmetic operations
//

void multiply(single_factor< alg_ideal > & c,
	      const single_factor< alg_ideal > & a,
	      const single_factor< alg_ideal > & b)
{
	//c = a*b
	if (a.is_prime_factor()) {
		if (b.is_prime_factor())
			multiply(c.rep_a, a.rep_p, b.rep_p);
		else
			multiply(c.rep_a, a.rep_p, b.rep_a);
	}
	else {
		if (b.is_prime_factor())
			multiply(c.rep_a, a.rep_a, b.rep_p);
		else
			multiply(c.rep_a, a.rep_a, b.rep_a);
	}
	c.set_prime_flag(decomposable_object::unknown);
}



void divide(single_factor< alg_ideal > & c,
	    const single_factor< alg_ideal > & a,
	    const single_factor< alg_ideal > & b)
{

#ifdef __xlC__
	single_factor< alg_ideal > s = c.rep_a;

	if (a.is_prime_factor()) {
		if (b.is_prime_factor())
			divide(s, a.rep_p, b.rep_p);
		else
			divide(s, a.rep_p, b.rep_a);
	}
	else {
		if (b.is_prime_factor())
			divide(s, a.rep_a, b.rep_p);
		else
			divide(s, a.rep_a, b.rep_a);
	}
	c.rep_a = s.rep_a;
#else
	if (a.is_prime_factor()) {
		if (b.is_prime_factor())
			divide(c.rep_a, a.rep_p, b.rep_p);
		else
			divide(c.rep_a, a.rep_p, b.rep_a);
	}
	else {
		if (b.is_prime_factor())
			divide(c.rep_a, a.rep_a, b.rep_p);
		else
			divide(c.rep_a, a.rep_a, b.rep_a);
	}
#endif

	c.set_prime_flag(decomposable_object::unknown);
	if (!(c.rep_a.denominator().is_one()))
		lidia_error_handler("single_factor< alg_ideal >",
				    "divide(...)::quotient is not integral!");

}



lidia_size_t ord_divide(const single_factor< alg_ideal > &a,
                        single_factor< alg_ideal > &b)
{
	if (a.is_one())
		lidia_error_handler("single_factor< alg_ideal >",
				    "ord_divide::1st argument mustn't be 1");

	lidia_size_t expo = 0;
	if (a.is_prime_factor()) {
		if (b.is_prime_factor()) {
			if ((a.rep_p.first_generator() == b.rep_p.first_generator()) &&
			    (a.rep_p.second_generator() == b.rep_p.second_generator())) {
				expo = 1;
				b.rep_a.assign_one();
				b.set_prime_flag(decomposable_object::not_prime);
			}
		}
		else {
			const prime_ideal & aa = a.rep_p;
			expo = static_cast<lidia_size_t>(ord(aa, b.rep_a));
			if (expo) {
				alg_ideal divisor;
				power(divisor, aa, expo);
				divide(b.rep_a, b.rep_a, divisor);
			}
		}
	}
	else {
		alg_ideal q;
		if (b.is_prime_factor())
			divide(q, b.rep_p, a.rep_a);
		else
			divide(q, b.rep_a, a.rep_a);
		while (q.denominator().is_one()) {
			expo++;
			swap(q, b.rep_a);
			divide(q, b.rep_a, a.rep_a);
		}
	}
	return expo;
}



bool
single_factor< alg_ideal >::
is_prime_factor(int test)
{
	lidia_debug_handler("single_factor< alg_ideal >", "is_prime_factor(int)");

	if (prime_flag() == decomposable_object::prime)
		return true;

	if (test == 0)              // =  > no explicit primality test
		return false;

	lidia_error_handler("single_factor< alg_ideal >",
			    "Doing a hard primality test of an "
			    "ideal is not supported\n "
			    "(And it woudn't make any sense anyway. "
			    "Just factor the ideal!)");
	return false;
}



factorization< alg_ideal >
single_factor< alg_ideal >::factor(int upper_bound) const
{
	factorization< alg_ideal > fact;
	factor(fact, upper_bound);
	return fact;
}



void
single_factor< alg_ideal >::factor(factorization< alg_ideal > &fact,
				    int upper_bound) const
{
	if (is_prime_factor()) {
		fact.assign(*this);
		return;
	}
	rational_factorization rf;

	bigrational expo = exponent(rep_a);
	rf.assign(expo*rep_a.denominator());
	rf.verbose(verbose_flag);
	rf.factor(upper_bound);
	if (verbose_flag)
		std::cout << "Factorization of exponent is" << rf << std::endl;
	finish(fact, rf);
	if (!rep_a.denominator().is_one()) {
		factorization< alg_ideal > fact2;
		rf.assign(rep_a.denominator());
		rf.verbose(verbose_flag);
		rf.factor(upper_bound);
		if (verbose_flag)
			std::cout << "Factorization of denominator is" << rf << std::endl;
		alg_ideal tmp1(rep_a.denominator(), alg_number(0, rep_a.which_base()));
		(static_cast< single_factor< alg_ideal > >(tmp1)).finish(fact2, rf);
		divide(fact, fact, fact2);
	}
}



factorization< alg_ideal >
factor(const alg_ideal & a, int upper_bound)
{
	factorization< alg_ideal > fact;
	(static_cast< single_factor< alg_ideal > >(a)).factor(fact, upper_bound);
	return fact;
}



void factor(factorization< alg_ideal > &fact, const alg_ideal & a,
            int upper_bound)
{
	(static_cast< single_factor< alg_ideal > >(a)).factor(fact, upper_bound);
}



void decompose_prime(prime_ideal * &factor, lidia_size_t & num,
		     const bigint & p, const order & O)
{
	register lidia_size_t j;

	if (!(O.base_denominator() % p).is_zero()) {
		if (single_factor< alg_ideal >::verbose_flag) {
			std::cout << "\n   Using factorization of polynomials to compute PI over ";
			std::cout << p << ":\n";
		}
		debug_handler_c("single_factor< alg_ideal >", "finish(...)::", 3,
				std::cout << "We compute in " << O);
		Fp_polynomial g (O.which_polynomial(), p);
		if (single_factor< alg_ideal >::verbose_flag) {
			std::cout << "Factoring polynomial " << g << std::endl;
		}
		factorization< Fp_polynomial > poly_factor;
		LiDIA::factor(poly_factor, g);
		// We assume, that factor computes a
		// decomposition into irreducible polynomials.
		poly_factor.sort();

		// Convert (p, poly_factor) to prime_ideals;
		num = poly_factor.no_of_prime_components();

		factor = new prime_ideal[num];
		if (num == 1 && poly_factor.prime_exponent(0) == 1)
			factor[0] = prime_ideal(p, O);
		else {
			if (single_factor< alg_ideal >::verbose_flag) {
				std::cout << "Contructing prime ideals" << std::endl;
			}
			bigint * tmp;
			for (j = 0; j < num; j++) {
				bigint_matrix A(g.degree(), 1);
				bigint deno(1);

				bigint_matrix B(O.base_numerator());
				lidia_size_t i;
				const Fp_polynomial & fact = poly_factor.prime_base(j).base();
				for (i = 0; i <= fact.degree(); i++)
					A.sto(i, 0, fact[i]);
				if (B.get_no_of_columns() > 1) {
					A.reginvimage(B, A);
					deno = A.member(g.degree(), 0);
					A.set_no_of_rows(g.degree());
					A *= O.base_denominator();
				}
				alg_number x(tmp = A.get_column(0), deno, O);
				delete[] tmp;

				factor[j] = prime_ideal(p, x, poly_factor.prime_exponent(j),
							fact.degree());
				if (single_factor< alg_ideal >::verbose_flag) {
					std::cout << " Constructed prime_ideal(" << p << ", " << x << ")" << std::endl;
				}
			}
		}
	}
	else {
		if (single_factor< alg_ideal >::verbose_flag) {
			std::cout << "\n   Using index divisor algorithm to compute PI over ";
			std::cout << p << ":\n";
		}
		factor_p(p, O, num, factor);
	}
	if (single_factor< alg_ideal >::verbose_flag) {
		std::cout << "Found the following candidates :\n";
		for (j = 0; j < num; j++)
			std::cout << j << ": " << factor[j] << std::endl;
	}
}



factorization< alg_ideal >
single_factor< alg_ideal >::finish(rational_factorization &rf) const
{
	factorization< alg_ideal > fact;
	finish(fact, rf);
	return fact;
}



factorization< alg_ideal > finish(const alg_ideal &a,
				  rational_factorization &rf)
{
	factorization< alg_ideal > fact;
	(static_cast< single_factor< alg_ideal > >(a)).finish(fact, rf);
	return fact;
}



void
single_factor< alg_ideal >::finish(factorization< alg_ideal > &fact,
				    rational_factorization &rf) const
{
	if (is_prime_factor()) {
		fact.assign(*this);
		return;
	}

	if (rf.base(0).is_one()) {
		fact.assign(factorization< alg_ideal > ());
		return;
	}

	prime_ideal *factor;
	bigint p;
	lidia_size_t num;
	lidia_size_t expo;
	int max_inertia;
	order O (rep_a.which_base());
	alg_ideal tmp(rep_a.which_base());
	alg_ideal tmp_this(rep_a);

//  alg_ideal divisor(O);

	register lidia_size_t j;

	fact.reset();
	timer T1;
	if (verbose_flag)
		std::cout << "Finish " << (*this) << " with exponent " << rf << std::endl;
	for (long i = 0; i < rf.no_of_comp(); i++) {
		p = rf.base(i);
		if (!rf.is_prime_factor(i)) {
			add(tmp, rep_a, alg_ideal(alg_number(p, rep_a.which_base())));
			while (!tmp.is_whole_order()) {
				alg_ideal non_prime_factor(tmp);
				lidia_size_t expo = 0;
				while (tmp == non_prime_factor) {
					divide(tmp_this, tmp_this, tmp);
					expo++;
					add(tmp, non_prime_factor,
					    alg_ideal(alg_number(p, rep_a.which_base())));
				}
				if (expo)
					fact.append(non_prime_factor, expo);
			}
			continue;
		}
		// We know, we are working on a prime number:
		max_inertia = rf.exponent(i) * O.degree();
		decompose_prime(factor, num, p, O);

		for (j = 0; j < num; j++) {
			int current_inertia = static_cast<int>(factor[j].degree_of_inertia());
			bool divides = (current_inertia <= max_inertia);
			if (verbose_flag) {
				std::cout << "Computing exponent for " << j << "-th candidate(";
				std::cout << factor[j] << ")." << std::endl;
			}
			if (divides) {
				T1.start_timer();
				expo = static_cast<lidia_size_t>(ord(factor[j], rep_a));
				T1.stop_timer();
				if (verbose_flag) {
					std::cout << "Computing ord needed " << T1 << std::endl;
				}
				if (expo) {
					if (verbose_flag) {
						std::cout << "Found factor: " << factor[j] << " ^ " << expo << std::endl;
					}
					single_factor< alg_ideal > tmp_factor(factor[j]);
					fact.swapping_append(tmp_factor, expo);
					max_inertia -= expo * current_inertia;
					//    power(tmp, factor[j], expo);
					//    multiply(divisor,divisor,tmp);
				}
			}
		}
		delete[] factor;
	}
}



void finish(factorization< alg_ideal > &fact,
            const alg_ideal &a, rational_factorization &rf)
{
	(static_cast < single_factor< alg_ideal > >(a)).finish(fact, rf);
}



factorization< alg_ideal >
single_factor< alg_ideal >::trialdiv(unsigned int upper_bound,
				      unsigned int lower_bound) const
{
	factorization< alg_ideal > fact;
	trialdiv(fact, upper_bound, lower_bound);
	return fact;
}



void single_factor< alg_ideal >::trialdiv(factorization< alg_ideal > &fact,
					   unsigned int upper_bound,
					   unsigned int lower_bound) const
{
	if (is_prime_factor()) {
		fact.assign(*this);
		return;
	}
	rational_factorization rf;
	rf.verbose(verbose_flag);

	bigrational expo = exponent(rep_a);
	rf.assign(expo*rep_a.denominator());
	rf.verbose(verbose_flag);
	rf.trialdiv(upper_bound, lower_bound);
	if (verbose_flag)
		std::cout << "Factorization of exponent is" << rf << std::endl;
	finish(fact, rf);
	if (!rep_a.denominator().is_one()) {
		factorization< alg_ideal > fact2;
		rf.assign(rep_a.denominator());
		rf.verbose(verbose_flag);
		rf.trialdiv(upper_bound, lower_bound);
		if (verbose_flag)
			std::cout << "Factorization of denominator is" << rf << std::endl;
		alg_ideal tmp1(rep_a.denominator(), alg_number(0, rep_a.which_base()));
		(static_cast< single_factor< alg_ideal > >(tmp1)).finish(fact2, rf);
		divide(fact, fact, fact2);
	}
}



factorization< alg_ideal > trialdiv(const alg_ideal &a,
				    unsigned int upper_bound,
				    unsigned int lower_bound)
{
	factorization< alg_ideal > fact;
	(static_cast< single_factor< alg_ideal > >(a)).trialdiv(fact, upper_bound, lower_bound);
	return fact;
}



void trialdiv(factorization< alg_ideal > &fact, const alg_ideal &a,
              unsigned int upper_bound,
              unsigned int lower_bound)
{
	(static_cast< single_factor< alg_ideal > >(a)).trialdiv(fact, upper_bound, lower_bound);
}



factorization< alg_ideal >
single_factor< alg_ideal >::ecm(int upper_bound,
				 int lower_bound, int step) const
{
	factorization< alg_ideal > fact;
	ecm(fact, upper_bound, lower_bound, step);
	return fact;
}



factorization< alg_ideal > ecm(const alg_ideal &a,
			       int upper_bound,
			       int lower_bound, int step)
{
	factorization< alg_ideal > fact;
	(static_cast< single_factor< alg_ideal > >(a)).ecm(fact, upper_bound, lower_bound, step);
	return fact;
}



void single_factor< alg_ideal >::ecm(factorization< alg_ideal > &fact,
				      int upper_bound, int lower_bound,
				      int step) const
{
	if (is_prime_factor()) {
		fact.assign(*this);
		return;
	}
	rational_factorization rf;
	rf.verbose(verbose_flag);

	bigrational expo = exponent(rep_a);
	rf.assign(expo*rep_a.denominator());
	rf.verbose(verbose_flag);
	rf.ecm(upper_bound, lower_bound, step);
	if (verbose_flag)
		std::cout << "Factorization of exponent is" << rf << std::endl;
	finish(fact, rf);
	if (!rep_a.denominator().is_one()) {
		factorization< alg_ideal > fact2;
		rf.assign(rep_a.denominator());
		rf.verbose(verbose_flag);
		rf.ecm(upper_bound, lower_bound, step);
		if (verbose_flag)
			std::cout << "Factorization of denominator is" << rf << std::endl;
		alg_ideal tmp1(rep_a.denominator(), alg_number(0, rep_a.which_base()));
		(static_cast< single_factor< alg_ideal > >(tmp1)).finish(fact2, rf);
		divide(fact, fact, fact2);
	}
}



void ecm(factorization< alg_ideal > &fact, const alg_ideal &a,
         int upper_bound, int lower_bound, int step)
{
	(static_cast< single_factor< alg_ideal > >(a)).ecm(fact, upper_bound, lower_bound, step);
}



factorization< alg_ideal > single_factor< alg_ideal >::mpqs() const
{
	factorization< alg_ideal > fact;
	mpqs(fact);
	return fact;
}



factorization< alg_ideal > mpqs(const alg_ideal &a)
{
	factorization< alg_ideal > fact;
	(static_cast< single_factor< alg_ideal > >(a)).mpqs(fact);
	return fact;
}



void single_factor< alg_ideal >::mpqs(factorization< alg_ideal > &fact) const
{
	if (is_prime_factor()) {
		fact.assign(*this);
		return;
	}
	rational_factorization rf;
	rf.verbose(verbose_flag);

	bigrational expo = exponent(rep_a);
	rf.assign(expo*rep_a.denominator());
	rf.verbose(verbose_flag);
	rf.mpqs_comp(0);
	if (verbose_flag)
		std::cout << "Factorization of exponent is" << rf << std::endl;
	finish(fact, rf);
	if (!rep_a.denominator().is_one()) {
		factorization< alg_ideal > fact2;
		rf.assign(rep_a.denominator());
		rf.verbose(verbose_flag);
		rf.mpqs_comp(0);
		if (verbose_flag)
			std::cout << "Factorization of denominator is" << rf << std::endl;
		alg_ideal tmp1(rep_a.denominator(), alg_number(0, rep_a.which_base()));
		(static_cast< single_factor< alg_ideal > >(tmp1)).finish(fact2, rf);
		divide(fact, fact, fact2);
	}
}



void mpqs(factorization< alg_ideal > &fact, const alg_ideal &a)
{
	(static_cast< single_factor< alg_ideal > >(a)).mpqs(fact);
}



// we offer a faster comparisons using additional information for prime_ideals
bool operator == (const single_factor< alg_ideal > & a,
                  const single_factor< alg_ideal > & b)
{
	debug_handler("single_factor< alg_ideal >",
		      "operator == (single_factor< alg_ideal > &)");
	if (a.prime_flag() == decomposable_object::prime) {
		if (b.prime_flag() == decomposable_object::prime)
			return ((a.rep_p.first_generator() == b.rep_p.first_generator()) &&
				(a.rep_p.second_generator() == b.rep_p.second_generator()));
		else
			return (static_cast<alg_ideal>(a.rep_p) == a.rep_a);
	}
	else {
		if (b.prime_flag() == decomposable_object::prime)
			return (a.rep_a == static_cast<alg_ideal>(b.rep_p));
		else
			return (a.rep_a == b.rep_a);
	}
}



void gcd(single_factor< alg_ideal > &c,
         const single_factor< alg_ideal > &a, const single_factor< alg_ideal > &b)
{
	debug_handler("single_factor< alg_ideal >", "gcd(...)");

	if (a.prime_flag() == decomposable_object::prime)
		if ((b.prime_flag() == decomposable_object::prime
		     && a.rep_p.first_generator() == b.rep_p.first_generator()
		     && a.rep_p.second_generator() == b.rep_p.second_generator())
		    || (ord(a.rep_p, b.rep_a)))
			c.assign(a);
		else {
			c.rep_a.assign_one();
			c.set_prime_flag(decomposable_object::not_prime);
		}
	else if (b.prime_flag() == decomposable_object::prime)
		if (ord(b.rep_p, a.rep_a))
			c.assign(b);
		else {
			c.rep_a.assign_one();
			c.set_prime_flag(decomposable_object::not_prime);
		}
	else {
		c.set_prime_flag(decomposable_object::unknown);
		add(c.rep_a, a.rep_a, b.rep_a);
	}
	return;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
