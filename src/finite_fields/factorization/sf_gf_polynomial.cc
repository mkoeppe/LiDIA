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
#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool single_factor< gf_polynomial >::verbose_flag = false;



//*********************************************************************
//	class single_factor< gf_polynomial >
//*********************************************************************

single_factor< gf_polynomial >::single_factor ()
	: base_factor< gf_polynomial > ()
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< gf_polynomial >", "single_factor< gf_polynomial > ()");

	//DEFAULT VALUE MUST BE '1', i.e. the identity element of multiplication !!!
	//we have a problem here : we must assign '1' to a polynomial without
	//knowing the field the polynomial is defined over
	//therefore we use the additional flag 'know_field', and leave 'rep' uninitialized

	//rep.assign_one();		//rep = 1;
	know_field = false;
}



single_factor< gf_polynomial >::single_factor(const single_factor< gf_polynomial > & x)
	: base_factor< gf_polynomial > (x),
	know_field(x.know_field)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< gf_polynomial >",
		      "single_factor< gf_polynomial > (single_factor< gf_polynomial > &)");
}



single_factor< gf_polynomial >::single_factor (const gf_polynomial & x)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< gf_polynomial >", "single_factor< gf_polynomial > (gf_polynomial&)");

	rep.assign(x);
	set_prime_flag(unknown);
	know_field = true;
}



single_factor< gf_polynomial >::~single_factor ()
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASSS

	debug_handler("single_factor< gf_polynomial >", "~single_factor< gf_polynomial > ()");
}



void
swap (single_factor< gf_polynomial > & a, single_factor< gf_polynomial > & b)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< gf_polynomial >", "swap(single_factor< gf_polynomial > &)");

	LiDIA::swap(b.rep, a.rep);

	single_factor< gf_polynomial >::decomp_state tmp = b.prime_flag();
	b.set_prime_flag(a.prime_flag());
	a.set_prime_flag(tmp);

	bool bb = b.know_field;
	b.know_field = a.know_field;
	a.know_field = bb;
}



single_factor< gf_polynomial > &
single_factor< gf_polynomial >::operator = (const single_factor< gf_polynomial > & x)
{
	debug_handler("single_factor< gf_polynomial >", "operator = (single_factor< gf_polynomial > &)");

	assign(x);
	return *this;
}



single_factor< gf_polynomial > &
single_factor< gf_polynomial >::operator = (const gf_polynomial & x)
{
	debug_handler("single_factor< gf_polynomial >", "operator = (gf_polynomial&)");

	assign(x);
	return *this;
}



void
single_factor< gf_polynomial >::assign (const single_factor< gf_polynomial > & x)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< gf_polynomial >", "assign(single_factor< gf_polynomial > &)");

	if (x.know_field == false) {
		know_field = false;
		return;
	}
	rep.assign(x.rep);
	set_prime_flag(x.prime_flag());
	know_field = true;
}



void
single_factor< gf_polynomial >::assign (const gf_polynomial & x)
{
	//WARNING : MUST BE REWRITTEN IF YOU DECLARE NEW DATA ELEMENTS IN THIS CLASS

	debug_handler("single_factor< gf_polynomial >", "assign(gf_polynomial&)");

	rep.assign(x);
	set_prime_flag(unknown);
	know_field = true;
}



gf_polynomial
single_factor< gf_polynomial >::extract_unit ()
{
	debug_handler("single_factor< gf_polynomial >", "extract_unit()");

	if (!know_field)
		lidia_error_handler("single_factor< gf_polynomial >",
				    "extract_unit(void)::single_factor is not initialized - field is unknown");

	gf_element lc = rep.lead_coeff();
	if (!lc.is_one()) {
		gf_element lc_inv;
		invert(lc_inv, lc);
		LiDIA::multiply(rep, rep, lc_inv);
	}
	gf_polynomial tmp(lc.get_field());
	tmp.set_degree(0);
	tmp[0] = lc;
	return tmp;

	//
	// this function 'normalizes' rep (i.e. rep is always made monic),
	// and returns the former leading coefficient of rep
	//
}



lidia_size_t
ord_divide (const single_factor< gf_polynomial > &a,
	    single_factor< gf_polynomial > &b)
{
	if (!a.know_field || !b.know_field)
		lidia_error_handler("single_factor< gf_polynomial >", "ord_divide::modulus must be set");
	if (a.is_one())
		lidia_error_handler("single_factor< gf_polynomial >", "ord_divide::1st argument mustn't be 1");

	gf_polynomial q, r;
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
single_factor< gf_polynomial >::is_prime_factor (int test)
{
	debug_handler("single_factor< gf_polynomial >", "is_prime_factor(int)");

	if (prime_flag() == prime)
		return true;

	if (test == 0)		//=  > no explicit primality test
		return false;

	//********************************************************
	// IMPLEMENT YOUR OWN PRIMALITY TEST HERE
	// BEGIN
	//********************************************************

	if (det_irred_test(rep)) {
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



factorization< gf_polynomial >
factor (const single_factor< gf_polynomial > & f)
{
	return f.factor();
}



//we need new comparisons because of the flag 'know_field'
//remember: if know_field == false, then rep == '1' (we don't know over
//which finite field - but this shouldn't matter)
bool
operator == (const single_factor< gf_polynomial > & a, const single_factor< gf_polynomial > & b)
{
	debug_handler ("single_factor< gf_polynomial >", "operator == (single_factor< gf_polynomial > &)");

	if (!a.know_field) {
		//i.e. if a == '1'
		if (!b.know_field) return true; //because a == b == '1'
		if (b.rep.is_one())  return true;
		else                 return false;
	}

	if (!b.know_field) {
		//i.e. if b == '1'
		if (a.rep.is_one())  return true;
		else                 return false;
	}

	return (a.rep == b.rep);
}



void
multiply (single_factor< gf_polynomial > &c,
	  const single_factor< gf_polynomial > &a,
	  const single_factor< gf_polynomial > &b)
{
	debug_handler("single_factor< gf_polynomial >", "multiply(...)");

	c.set_prime_flag(decomposable_object::unknown);
	if (a.know_field && b.know_field) {
		multiply(c.rep, a.rep, b.rep);
		c.know_field = true;
		return;
	}
	if (!a.know_field && !b.know_field) {
		//c.rep.set_degree(-1);
		c.know_field = false;
		return;
	}

	if (a.know_field)
		c.assign(a);
	else	// if (b.know_field)  
		c.assign(b);
}



void
divide (single_factor< gf_polynomial > &c,
	const single_factor< gf_polynomial > &a,
	const single_factor< gf_polynomial > &b)
{
	debug_handler("single_factor< gf_polynomial >", "divide(...)");

	c.set_prime_flag(decomposable_object::unknown);
	if (a.know_field && b.know_field) {
		gf_polynomial qwe;
		div_rem(c.rep, qwe, a.rep, b.rep);
		if (!qwe.is_zero())
			lidia_error_handler("single_factor< gf_polynomial >", "divide(...)::remainder is not zero !");
		//divide(c.rep, a.rep, b.rep);

		c.know_field = true;
		return;
	}
	if (!a.know_field && !b.know_field) {
		c.rep.set_degree(-1);
		c.know_field = false;
		return;
	}

	if (a.know_field)
		c.assign(a);
	else {
		// if (b.know_field)  
		if (b.rep.degree() != 0)
			lidia_error_handler("single_factor< gf_polynomial >", "divide(...)::deg != 0 !");
		gf_polynomial tmp;
		tmp.assign_one(b.rep.get_field()); //must set tmp.base()
		divide(c.rep, tmp, b.rep);
		c.know_field = true;
	}
}



void
gcd (single_factor< gf_polynomial > &c,
     const single_factor< gf_polynomial > &a,
     const single_factor< gf_polynomial > &b)
{
	debug_handler("single_factor< gf_polynomial >", "gcd(...)");

	c.set_prime_flag(decomposable_object::unknown);
	if (a.know_field && b.know_field) {
		gcd(c.rep, a.rep, b.rep);
		c.know_field = true;
		return;
	}
	if (!a.know_field && !b.know_field) {
		c.rep.set_degree(-1);
		c.know_field = false;
		return;
	}

	if (a.know_field)
		c.rep.assign_one(a.rep.get_field());
	else	// if (b.know_field)  
		c.rep.assign_one(b.rep.get_field());
	c.know_field = true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
