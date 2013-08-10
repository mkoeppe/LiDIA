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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/Fp_polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool single_factor< Fp_polynomial >::verbose_flag = false;



//*********************************************************************
// class single_factor< Fp_polynomial >
//*********************************************************************

single_factor< Fp_polynomial >::
single_factor() : base_factor< Fp_polynomial > ()
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< Fp_polynomial >", "single_factor< Fp_polynomial > ()");

	//DEFAULT VALUE MUST BE '1', i.e. the identity element of multiplication !!!
	//we have a problem here : we must assign '1' to a polynomial without
	//knowing the field the polynomial is defined over
	//therefore we use the additional flag 'know_modulus', and leave 'rep' uninitialized

	//rep.assign_one();		//rep = 1;
	know_modulus = false;
}



single_factor< Fp_polynomial >::
single_factor(const single_factor< Fp_polynomial > & x) :
	base_factor< Fp_polynomial > (x),
	know_modulus(x.know_modulus)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< Fp_polynomial >", "single_factor< Fp_polynomial > (single_factor< Fp_polynomial > &)");
}



single_factor< Fp_polynomial >::
single_factor(const Fp_polynomial & x)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< Fp_polynomial >", "single_factor< Fp_polynomial > (Fp_polynomial&)");

	rep.assign(x);
	set_prime_flag(unknown);
	know_modulus = true;
}



single_factor< Fp_polynomial >::
~single_factor()
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASSS

	debug_handler("single_factor< Fp_polynomial >", "~single_factor< Fp_polynomial > ()");
}



void
swap(single_factor< Fp_polynomial > & a, single_factor< Fp_polynomial > & b)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< Fp_polynomial >", "swap(single_factor< Fp_polynomial > &, single_factor< Fp_polynomial > &)");

	LiDIA::swap(b.rep, a.rep);

	single_factor< Fp_polynomial >::decomp_state tmp = b.prime_flag();
	b.set_prime_flag(a.prime_flag());
	a.set_prime_flag(tmp);

	bool bb = b.know_modulus;
	b.know_modulus = a.know_modulus;
	a.know_modulus = bb;
}



single_factor< Fp_polynomial > &
single_factor< Fp_polynomial >::
operator = (const single_factor< Fp_polynomial > & x)
{
	debug_handler("single_factor< Fp_polynomial >", "operator = (single_factor< Fp_polynomial > &)");

	assign(x);
	return *this;
}



single_factor< Fp_polynomial > &
single_factor< Fp_polynomial >::
operator = (const Fp_polynomial & x)
{
	debug_handler("single_factor< Fp_polynomial >", "operator = (Fp_polynomial&)");

	assign(x);
	return *this;
}



void
single_factor< Fp_polynomial >::
assign(const single_factor< Fp_polynomial > & x)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< Fp_polynomial >", "assign(single_factor< Fp_polynomial > &)");

	if (x.know_modulus == false) {
		know_modulus = false;
		return;
	}
	rep.assign(x.rep);
	set_prime_flag(x.prime_flag());
	know_modulus = true;
}



void
single_factor< Fp_polynomial >::
assign(const Fp_polynomial & x)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< Fp_polynomial >", "assign(Fp_polynomial&)");

	rep.assign(x);
	set_prime_flag(unknown);
	know_modulus = true;
}



Fp_polynomial
single_factor< Fp_polynomial >::
extract_unit()
{
	debug_handler("single_factor< Fp_polynomial >", "extract_unit()");

	if (!know_modulus)
		lidia_error_handler("single_factor< Fp_polynomial >", "extract_unit(void)::single_factor is not initialized - modulus is unknown");

	bigint lc = rep.lead_coeff();
	if (!lc.is_one()) {
		bigint lc_inv;
		InvMod(lc_inv, lc, rep.modulus());
		LiDIA::multiply(rep, rep, lc_inv);
	}
	Fp_polynomial tmp;
	tmp.set_modulus(rep.modulus());
	tmp.set_coefficient(lc, 0);
	return tmp;

	//
	// this function 'normalizes' rep (i.e. rep is always made monic),
	// and returns the former leading coefficient of rep
	//
}



lidia_size_t ord_divide(const single_factor< Fp_polynomial > &a,
			single_factor< Fp_polynomial > &b)
{
	if (!a.know_modulus || !b.know_modulus)
		lidia_error_handler("single_factor< Fp_polynomial >", "ord_divide::modulus must be set");
	if (a.is_one())
		lidia_error_handler("single_factor< Fp_polynomial >", "ord_divide::1st argument mustn't be 1");

	Fp_polynomial q, r;
	lidia_size_t e = 0;
	div_rem(q, r, b.rep, a.rep);
	while (r.is_zero()) {
		e++;
		swap(q, b.rep);
		div_rem(q, r, b.rep, a.rep);
	}
	return e;
}



bool
single_factor< Fp_polynomial >::
is_prime_factor(int test)
{
	debug_handler("single_factor< Fp_polynomial >", "is_prime_factor(int)");

	if (prime_flag() == prime)
		return true;

	if (test == 0)		//=  > no explicit primality test
		return false;

	//********************************************************
	// IMPLEMENT YOUR OWN PRIMALITY TEST HERE
	// BEGIN
	//********************************************************

	if (prob_irred_test(rep, 4)) {
		set_prime_flag(prime);
		return true;
	}
	else {
		set_prime_flag(not_prime);
		return false;
	}

	//********************************************************
	// END OF YOUR OWN PRIMALITY TEST
	//********************************************************
}



factorization< Fp_polynomial >
factor(const single_factor< Fp_polynomial > & f)
{
	return f.factor();
}



//we need new comparisons because of the flag 'know_modulus'
//remember: if know_modulus == false, then rep == '1' (we don't know over
//which finite field - but this shouldn't matter)
bool operator == (const single_factor< Fp_polynomial > & a, const single_factor< Fp_polynomial > & b)
{
	debug_handler ("single_factor< Fp_polynomial >", "operator == (single_factor< Fp_polynomial > &)");

	if (!a.know_modulus) {
		//i.e. if a == '1'
		if (!b.know_modulus) return true; //because a == b == '1'
		if (b.rep.is_one())  return true;
		else                 return false;
	}

	if (!b.know_modulus) {
		//i.e. if b == '1'
		if (a.rep.is_one())  return true;
		else                 return false;
	}

	return (a.rep == b.rep);
}



void multiply(single_factor< Fp_polynomial > &c, const single_factor< Fp_polynomial > &a, const single_factor< Fp_polynomial > &b)
{
	debug_handler("single_factor< Fp_polynomial >", "multiply(...)");

	c.set_prime_flag(decomposable_object::unknown);
	if (a.know_modulus && b.know_modulus) {
		multiply(c.rep, a.rep, b.rep);
		c.know_modulus = true;
		return;
	}
	if (!a.know_modulus && !b.know_modulus) {
		c.rep.kill();
		c.know_modulus = false;
		return;
	}

	if (a.know_modulus)
		c.assign(a);
	else	// if (b.know_modulus)  
		c.assign(b);
}



void divide(single_factor< Fp_polynomial > &c, const single_factor< Fp_polynomial > &a, const single_factor< Fp_polynomial > &b)
{
	debug_handler("single_factor< Fp_polynomial >", "divide(...)");

	c.set_prime_flag(decomposable_object::unknown);
	if (a.know_modulus && b.know_modulus) {
		Fp_polynomial qwe;
		div_rem(c.rep, qwe, a.rep, b.rep);
		if (!qwe.is_zero())
			lidia_error_handler("single_factor< Fp_polynomial >", "divide(...)::remainder is not zero !");
		//	instead of: divide(c.rep, a.rep, b.rep);

		c.know_modulus = true;
		return;
	}
	if (!a.know_modulus && !b.know_modulus) {
		c.rep.kill();
		c.know_modulus = false;
		return;
	}

	if (a.know_modulus)
		c.assign(a);
	else {
		// if (b.know_modulus)  
		Fp_polynomial tmp(b.rep); //copy b in order to set tmp.modulus()
		if (tmp.degree() != 0)
			lidia_error_handler("single_factor< Fp_polynomial >", "divide(...)::deg != 0 !");

		tmp.assign_one();
		divide(c.rep, tmp, b.rep);
		c.know_modulus = true;
	}
}



void gcd(single_factor< Fp_polynomial > &c, const single_factor< Fp_polynomial > &a, const single_factor< Fp_polynomial > &b)
{
	debug_handler("single_factor< Fp_polynomial >", "gcd(...)");

	c.set_prime_flag(decomposable_object::unknown);
	if (a.know_modulus && b.know_modulus) {
		gcd(c.rep, a.rep, b.rep);
		c.know_modulus = true;
		return;
	}
	if (!a.know_modulus && !b.know_modulus) {
		c.rep.kill();
		c.know_modulus = false;
		return;
	}

	if (a.know_modulus)
		c.rep.set_modulus(a.rep.modulus());
	else	// if (b.know_modulus)  
		c.rep.set_modulus(b.rep.modulus());
	c.rep.assign_one();
	c.know_modulus = true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
