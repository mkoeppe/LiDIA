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
//	Author	: Detlef Anton (DA), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/galois_field.h"
#include	"LiDIA/finite_fields/galois_field_rep.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************************
//*			class galois_field				*
//***********************************************************************

//
// constructors and destructor
//

galois_field::galois_field ()
{
	debug_handler("galois_field", "galois_field()");
	rep = galois_field_rep::create(0, 0);
}



galois_field::galois_field (const galois_field_rep& internal_rep)
{
	debug_handler("galois_field", "galois_field(galois_field_rep&)");
	rep = &internal_rep;
	rep->rc.inc_ref_counter();
}


galois_field::galois_field (const galois_field& K)
{
	debug_handler("galois_field", "galois_field(galois_field&)");
	rep = K.rep;
	K.rep->rc.inc_ref_counter();
}



galois_field::galois_field (const bigint & characteristic, unsigned int deg)
{
	debug_handler("galois_field", "galois_field(bigint&, unsigned int)");
	if (!is_prime(characteristic))
		lidia_error_handler("galois_field", "characteristic no prime "
				    "number");
	rep = galois_field_rep::create(characteristic, deg);
}



galois_field::galois_field (const bigint & characteristic, unsigned int deg,
			    const rational_factorization &fact)
{
	debug_handler("galois_field",
		      "galois_field(bigint&, unsigned int, rational_factorization&)");
	rep = galois_field_rep::create(characteristic, deg, fact);
}



galois_field::galois_field (const Fp_polynomial &pol)
{
	debug_handler("galois_field", "galois_field(Fp_polynomial&)");
	rep = galois_field_rep::create(pol);
}



galois_field::galois_field (const Fp_polynomial &pol,
			    const rational_factorization &fact)
{
	debug_handler("galois_field",
		      "galois_field(Fp_polynomial&, rational_factorization&)");
	rep = galois_field_rep::create(pol, fact);
}



galois_field::~galois_field()
{
	debug_handler("galois_field", "destructor");
	rep->rc.dec_ref_counter();
	if (rep->rc.get_ref_counter() == 0)
		delete rep;
}



//
// iterator access
//

galois_field_iterator
galois_field::begin() const {
  return galois_field_iterator(*this);
}


galois_field_iterator
galois_field::end() const {
  galois_field_iterator result(*this);
  result.move_past_the_end();
  return result;
}


//
// access functions
///

const bigint &
galois_field::characteristic () const
{
	return rep->characteristic();
}



unsigned int
galois_field::degree() const
{
	return rep->degree();
}



const bigint &
galois_field::number_of_elements () const
{
	return rep->number_of_elements();
}



const rational_factorization &
galois_field::factorization_of_mult_order () const
{
	return rep->factorization_of_mult_order();
}



Fp_polynomial
galois_field::irred_polynomial () const
{
	return rep->irred_polynomial();
}


gf_element const&
galois_field::generator() const {
  return rep->generator();
}


//
// assignment
//

void
galois_field::assign (const galois_field & K)
{
	if (rep != K.rep) {
		rep->rc.dec_ref_counter();
		if (rep->rc.get_ref_counter() == 0)
			delete rep;
		rep = K.rep;
		K.rep->rc.inc_ref_counter();
	}
}



galois_field &
galois_field::operator = (const galois_field& K)
{
	assign(K);
	return *this;
}



void
swap (galois_field& a, galois_field& b)
{
	galois_field_rep const* p = a.rep;
	a.rep = b.rep;
	b.rep = p;
}



//
// comparisons
//

bool
galois_field::operator < (const galois_field& K) const
{
	if (characteristic() != K.characteristic())
		return false;
	if (degree() == 1 && K.degree() > 1)
		return true;

	//if (K.degree() % degree() != 0)
	//if ...
	return false;
}



bool
galois_field::operator > (const galois_field& K) const
{
	return (K < *this);
}



//
// input / output
//
// I/O format:  (1)  "[p, n]"	   where p is the char. and n is the degree
//              (2)  "f(x) mod p"  where f is an irreducible polynomial mod p
//
// Input: format (1) or (2)
// Output: always format (2)
//

std::istream &
operator >> (std::istream & in, galois_field & K)
{
	bigint p;
	lidia_size_t n;
	Fp_polynomial pol;

	char c;
	in >> c;
	while (c == ' ') in >> c;
	if (c == '[') {
		// format (1)
		in >> p;
		in >> c;
		while (c == ' ') in >> c;
		if (c != ',')
			lidia_error_handler("galois_field", "operator >> : ', ' expected");
		in >> n;
		in >> c;
		while (c == ' ') in >> c;
		if (c != ']')
			lidia_error_handler("galois_field", "operator >> : ')' expected");

		galois_field K2(p, n);
		K.assign(K2);
	}
	else {
		// format (2)
		in.putback(c); // bigfix by Francois Dionne
		Fp_polynomial pol;
		in >> pol;
		galois_field K2(pol);
		K.assign(K2);
	}
	return in;
}



std::ostream &
operator << (std::ostream & out, const galois_field & K)
{
	// simple format:
	//out << "[" << K.characteristic() << ", " << K.degree() << "]";

	if (K.degree() > 0)
		out << K.irred_polynomial();
	else
		out << Fp_polynomial();
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
