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
//	Author	: Emre Binisik, Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include "LiDIA/bigint.h"
#include "LiDIA/single_factor.h"
#include "LiDIA/factorization.h"
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool single_factor< bigint >::info = false;


//*********************************************************************
// class single_factor< bigint >
//*********************************************************************

single_factor< bigint >::
single_factor() : base_factor< bigint > ()
{
	debug_handler("single_factor< bigint >",
		      "single_factor< bigint > ()");

	rep.assign_one();
	algo_info = 0;
}



single_factor< bigint >::
single_factor(const single_factor< bigint > & b) :
	base_factor< bigint > (b)
{
	debug_handler("single_factor< bigint >",
		      "single_factor< bigint > (const single_factor< bigint > &)");

	rep.assign(b.rep);
	algo_info = 0;
	set_prime_flag(b.prime_flag());
}



single_factor< bigint >::
single_factor(const bigint & b):
	base_factor< bigint > ()
{
	debug_handler("single_factor< bigint >",
		      "single_factor< bigint > (const bigint & b)");

	rep.assign(b);
	algo_info = 0;
	set_prime_flag(unknown);
}



single_factor< bigint >::
~single_factor()
{
	debug_handler("single_factor< bigint >", "~single_factor< bigint > ()");
}



void
single_factor< bigint >::
swap(single_factor< bigint > & a)
{
	debug_handler("single_factor< bigint >", "swap(single_factor< bigint > &)");
	LiDIA::swap(rep, a.rep);
	decomp_state tmp = prime_flag();
	set_prime_flag(a.prime_flag());
	a.set_prime_flag(tmp);

	unsigned long h = algo_info;
	algo_info = a.algo_info;
	a.algo_info = h;
}



void
single_factor< bigint >::
assign(const single_factor< bigint > & x)
{
	debug_handler("single_factor< bigint >", "assign(single_factor< bigint > &)");
	rep.assign(x.rep);
	algo_info = x.algo_info;
	set_prime_flag(x.prime_flag());
}



void
single_factor< bigint >::
assign(const bigint & x)
{
	debug_handler("single_factor< bigint >", "assign(bigint &)");
	rep.assign(x);
	set_prime_flag(unknown);
	algo_info = 0;
}



single_factor< bigint > &
single_factor< bigint >::
operator = (const single_factor< bigint > & x)
{
	debug_handler("single_factor< bigint >", "operator = (bigint &)");
	assign(x);
	return *this;
}



single_factor< bigint > &
single_factor< bigint >::
operator = (const bigint & x)
{
	debug_handler("single_factor< bigint >", "operator = (bigint&)");
	assign(x);
	return *this;
}



bigint
single_factor< bigint >:: extract_unit()
{
	debug_handler("single_factor< bigint >", "extract_unit()");
	if (rep.is_negative()) {
		rep.negate();
		return bigint(-1);
	}
	else
		return bigint(1);
}



bool
single_factor< bigint >::is_prime_factor(int i)
{
	if (i == 0)
		return (prime_flag() == decomposable_object::prime);
	else
		if (prime_flag() == decomposable_object::unknown) {
			bool w = is_prime(rep, 8);
			if (w == true)
				set_prime_flag(prime);
			else
				set_prime_flag(not_prime);
			return w;
		}
		else
			if (prime_flag() == decomposable_object::prime)
				return true;
			else
				return false;
}



lidia_size_t ord_divide(const single_factor< bigint > & a,
			const single_factor< bigint > & b)
{
	if (abs(a.rep) == 1) {
		precondition_error_handler(a.rep, "a", "|a.rep| != 1",
					   "lidia_size_t "
					   "ord_divide("
					   "const single_factor<bigint>& a, "
					   "const single_factor<bigint>& b)",
					   "single_factor<bigint>",
					   "ord_divide::1st argument "
					   "mustn't be 1");
		return 0;
	}
	bigint hb(b.rep), r;
	lidia_size_t e = 0;

	div_rem(hb, r, hb, a.rep);

	while (r.is_zero()) {
		e++;
		div_rem(hb, r, hb, a.rep);
	}
	return e;
}



void gcd(single_factor< bigint > & c, const single_factor< bigint > & a,
         const single_factor< bigint > & b)
{
	debug_handler("single_factor< bigint >", "gcd(...)");
	c.set_prime_flag(decomposable_object::unknown);
	c.rep = gcd(a.rep, b.rep);
}



factorization< bigint > sf_factor(const bigint & N, int size)
{
	single_factor< bigint > f(N);
	return f.factor(size);
}

factorization< bigint > sf_factor(const bigint & N)
{
	single_factor< bigint > f(N);
	return f.factor();
}


factorization< bigint > single_factor< bigint >::
completely_factor() const
{
	factorization< bigint > f(*this);
	f.factor_all_components();
	return f;
}



factorization< bigint > completely_factor(const bigint & N)
{
	factorization< bigint > f(N);
	f.factor_all_components();
	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
