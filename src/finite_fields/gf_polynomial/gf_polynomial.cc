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



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

template<>
gf_element base_polynomial< gf_element >::
operator() (gf_element const& value) const {
	debug_handler_l(DM_BP,
			"in base_polynomial< gf_element >::operator "
			"() (gf_element const&)", DV_BP + 5);

	if (deg == 0)
		return coeff[0];

	gf_element result(value.get_field()); // initializes result = 0

	if (deg < 0)
		return result;

	result = coeff[deg];
	for (lidia_size_t i = deg - 1; i >= 0; i--) {
		LiDIA::multiply(result, result, value);
		LiDIA::add(result, result, coeff[i]);
	}
	return result;
}



static bool inline is_undefined(const galois_field &K)
{
	return K.degree() == 0;
}



const galois_field* polynomial< gf_element >::FIELD = 0;



const galois_field&
polynomial< gf_element >::
common_field(const galois_field &a, const galois_field &b)
// returns a common superfield of a and b
// in this version, an error is raised if a!=b
{
	if (is_undefined(a) && is_undefined(b))
		std::cout << "common_field::WARNING (both fields are uninitialized)" << std::endl;

	if (a != b)
		lidia_error_handler("polynomial< gf_element >",
				    "common_field::arguments are different fields");
	return a;
}



void polynomial< gf_element >::check_coefficients()
{
	lidia_size_t i, d = pol.degree();
//    if (ffield.is_undefined())
//	lidia_error_handler("gf_element","check_coefficients()::don't know over which field");

	if (d < 0)
		return;

	// look for field which is initialized
	galois_field K;
	for (i = d; i >= 0; i--) {
		if (!is_undefined(pol[i].get_field())) {
			K = pol[i].get_field();
			break;
		}
	}

	if (is_undefined(K))
		K = ffield;
	else
		ffield = K;

	if (!is_undefined(K)) {
		gf_element null(K);
		for (i = 0; i <= d; i++) {
			if (is_undefined(pol[i].get_field()))
				add(pol[i], pol[i], null); // set field
			else if (pol[i].get_field() != K)
				lidia_error_handler("gf_polynomial", "check_coefficients()::"
						    "coefficients lie in different fields");
		}
	}
}



void polynomial< gf_element >::set_degree(lidia_size_t d)
{
	build_frame(ffield);
	pol.set_degree(d);
	delete_frame();
}



gf_element polynomial< gf_element >::lead_coeff() const
{
//    if (ffield.is_undefined())
//	lidia_error_handler("polynomial<gf_element>","lead_coeff::no field given for this polynomial - don't know which 'zero' to return" );

	if (pol.degree() < 0)
		return gf_element(ffield);
	else
		return pol[pol.degree()];
}



gf_element polynomial< gf_element >::const_term() const
{
//    if (ffield.is_undefined())
//	lidia_error_handler("polynomial<gf_element>","const_term::deg < 0 and no field given for this polynomial - don't know which 'zero' to return" );

	if (pol.degree() < 0)
		return gf_element(ffield);
	else
		return pol[0];
}



void swap(polynomial< gf_element > &a, polynomial< gf_element > &b)
{
	swap(a.pol, b.pol);
	swap(a.ffield, b.ffield);
}



void polynomial< gf_element >::remove_leading_zeros()
{
	lidia_size_t d = degree();
	while (d >= 0 && pol[d].is_zero())
		d--;
	if (d != degree()) {
		build_frame(ffield);
		set_degree(d);
		delete_frame();
	}
}



void polynomial< gf_element >::set_field(const galois_field &K)
{
	this->assign_zero(K);
}



bool operator == (const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b)
{
	return (a.ffield == b.ffield && a.pol == b.pol);
}



bool operator != (const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b)
{
	return (a.ffield != b.ffield || a.pol != b.pol);
}



bool polynomial< gf_element >::is_monic() const
{
	if (is_zero()) return false;
	return (*this)[degree()].is_one();
}



void polynomial< gf_element >::assign(const polynomial< gf_element > &a)
{
	build_frame(a.ffield);
	ffield = a.ffield;
	pol.assign(a.pol);
	delete_frame();
}



void polynomial< gf_element >::assign(const gf_element &a)
{
	build_frame(a.get_field());
	ffield = a.get_field();
	pol.assign(a);
	delete_frame();
}



polynomial< gf_element > &
polynomial< gf_element >::operator = (const gf_element &a)
{
	this->assign(a);
	return *this;
}



void polynomial< gf_element >::assign_zero()
{
	if (is_undefined(ffield)) lidia_error_handler("polynomial< gf_element >", "assign_zero() :: field unknown");
	pol.assign_zero();
}



void polynomial< gf_element >::assign_one()
{
	if (is_undefined(ffield)) lidia_error_handler("polynomial< gf_element >", "assign_one() :: field unknown");
	build_frame(ffield);
	pol.set_degree(0);
	delete_frame();
	pol[0].assign_one(ffield);
}



void polynomial< gf_element >::assign_x()
{
	if (is_undefined(ffield)) lidia_error_handler("polynomial< gf_element >", "assign_x() :: field unknown");
	build_frame(ffield);
	pol.set_degree(1);
	delete_frame();
	pol[0].assign_zero(ffield);
	pol[1].assign_one(ffield);
}



void polynomial< gf_element >::assign_zero(const galois_field &K)
{
	ffield = K;
	pol.assign_zero();
}



void polynomial< gf_element >::assign_one(const galois_field &K)
{
	ffield = K;
	build_frame(K);
	pol.set_degree(0);
	delete_frame();
	pol[0].assign_one(K);
}



void polynomial< gf_element >::assign_x(const galois_field &K)
{
	ffield = K;
	build_frame(K);
	pol.set_degree(1);
	delete_frame();
	pol[0].assign_zero(K);
	pol[1].assign_one(K);
}



//for class single_factor< gf_polynomial >
bool operator< (const polynomial < gf_element > &a,
		const polynomial< gf_element > &b)
{
	lidia_size_t i, n = a.degree(), m = b.degree();
	if (n != m)
		return (n < m);
	for (i = n; i >= 0; i--) {
		if (a[i].polynomial_rep() < b[i].polynomial_rep())  return true;
		if (!(a[i].polynomial_rep() <= b[i].polynomial_rep()))  return false;
	}
	return false;
}



bool operator <= (const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b)
{
	lidia_size_t i, n = a.degree(), m = b.degree();
	if (n != m)
		return (n < m);
	for (i = n; i >= 0; i--) {
		if (a[i].polynomial_rep() < b[i].polynomial_rep())  return true;
		if (!(a[i].polynomial_rep() <= b[i].polynomial_rep()))  return false;
	}
	return true;
}



void negate(polynomial< gf_element > & c,
	    const polynomial< gf_element > &a)
{
	polynomial< gf_element >::build_frame(a.ffield);
	c.ffield = a.ffield;
	negate(c.pol, a.pol);
	polynomial< gf_element >::delete_frame();
}



void add(polynomial< gf_element > & c,
	 const polynomial< gf_element > & a,
	 const polynomial< gf_element > & b)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.ffield);
	polynomial< gf_element >::build_frame(c.ffield);
	add(c.pol, a.pol, b.pol);
	polynomial< gf_element >::delete_frame();
}



void add(polynomial< gf_element > & c,
	 const polynomial< gf_element > & a, const gf_element & b)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	add(c.pol, a.pol, b);
	polynomial< gf_element >::delete_frame();
}



void add(polynomial< gf_element > & c,
	 const gf_element & b, const polynomial< gf_element > & a)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	add(c.pol, a.pol, b);
	polynomial< gf_element >::delete_frame();
}



void subtract(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a,
	      const polynomial< gf_element > & b)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.ffield);
	polynomial< gf_element >::build_frame(c.ffield);
	subtract(c.pol, a.pol, b.pol);
	polynomial< gf_element >::delete_frame();
}



void subtract(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a, const gf_element & b)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	subtract(c.pol, a.pol, b);
	polynomial< gf_element >::delete_frame();
}



void subtract(polynomial< gf_element > & c,
	      const gf_element & b, const polynomial< gf_element > & a)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	subtract(c.pol, b, a.pol);
	polynomial< gf_element >::delete_frame();
}



void multiply(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a, const gf_element & b)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	multiply(c.pol, a.pol, b);
	polynomial< gf_element >::delete_frame();
}



void multiply(polynomial< gf_element > & c,
	      const gf_element & b, const polynomial< gf_element > & a)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	multiply(c.pol, a.pol, b);
	c.delete_frame();
	polynomial< gf_element >::delete_frame();
}



void divide(polynomial< gf_element > & c,
	    const polynomial< gf_element > & a, const gf_element & b)
{
	c.ffield = polynomial< gf_element >::common_field(a.ffield, b.get_field());
	polynomial< gf_element >::build_frame(c.ffield);
	divide(c.pol, a.pol, b);
	c.delete_frame();
	polynomial< gf_element >::delete_frame();
}



gf_element
polynomial< gf_element >::operator() (const gf_element & value) const
{
	if (pol.degree() < 0)
		return gf_element(ffield);
	else
		return pol(value);
}



polynomial< gf_element > operator - (const polynomial< gf_element > &a)
{
	polynomial< gf_element > c;
	negate(c, a);
	return c;
}



polynomial< gf_element > operator + (const polynomial< gf_element > &a,
				     const polynomial< gf_element > &b)
{
	polynomial< gf_element > c;
	add(c, a, b);
	return c;
}



polynomial< gf_element > operator + (const polynomial< gf_element > & a,
				     const gf_element& b)
{
	polynomial< gf_element > c;
	add(c, a, b);
	return c;
}



polynomial< gf_element > operator + (const gf_element & b,
				     const polynomial< gf_element > & a)
{
	polynomial< gf_element > c;
	add(c, a, b);
	return c;
}



polynomial< gf_element > operator - (const polynomial< gf_element > &a,
				     const polynomial< gf_element > &b)
{
	polynomial< gf_element > c;
	subtract(c, a, b);
	return c;
}



polynomial< gf_element > operator - (const polynomial< gf_element > &a,
				     const gf_element &b)
{
	polynomial< gf_element > c;
	subtract(c, a, b);
	return c;
}



polynomial< gf_element > operator - (const gf_element &a,
				     const polynomial< gf_element > &b)
{
	polynomial< gf_element > c;
	subtract(c, a, b);
	return c;
}



polynomial< gf_element > operator * (const polynomial< gf_element > &a,
				     const polynomial< gf_element > &b)
{
	polynomial< gf_element > c;
	multiply(c, a, b);
	return c;
}



polynomial< gf_element > operator * (const polynomial< gf_element > &a,
				     const gf_element &b)
{
	polynomial< gf_element > c;
	multiply(c, a, b);
	return c;
}



polynomial< gf_element > operator * (const gf_element &b,
				     const polynomial< gf_element > &a)
{
	polynomial< gf_element > c;
	multiply(c, a, b);
	return c;
}



polynomial< gf_element > &
polynomial< gf_element >::operator += (const polynomial< gf_element > &a)
{
	add(*this, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator += (const gf_element &a)
{
	add(*this, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator -= (const polynomial< gf_element > &a)
{
	subtract(*this, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator -= (const gf_element &a)
{
	subtract(*this, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator *= (const polynomial< gf_element > &a)
{
	multiply(*this, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator *= (const gf_element &a)
{
	multiply(*this, *this, a);
	return *this;
}



void derivative(polynomial< gf_element > &c,
		const polynomial< gf_element > &a)
{
	c.build_frame(a.ffield);
	c.ffield = a.ffield;
	derivative(c.pol, a.pol);
	c.delete_frame();
}



polynomial< gf_element > derivative(const polynomial< gf_element > & a)
{
	polynomial< gf_element > c;
	derivative(c, a);
	return c;
}



void integral(polynomial< gf_element > &c,
	      const polynomial< gf_element > &a)
{
	c.build_frame(a.ffield);
	c.ffield = a.ffield;
	integral(c.pol, a.pol);
	c.delete_frame();
}



polynomial< gf_element > randomize(const galois_field &K, lidia_size_t n)
{
	polynomial< gf_element > f;
	f.ffield = K;
	f.set_degree(n);
	for (lidia_size_t i = 0; i < n; i++) {
		f[i].assign_zero(K); // initialize with field K
		f[i].randomize();
	}
	f[n].assign_zero(K); // initialize with field K
	do
		f[n].randomize();
	while (f[n].is_zero());
	return f;
}



std::istream & polynomial< gf_element >::read_verbose(std::istream &is)
{
	if (is_undefined(ffield))
		lidia_error_handler("polynomial< gf_element >", "read_verbose(...)::polynomial must be assigned to a field before any input");
	build_frame(ffield);
	pol.read_verbose(is);
	delete_frame();
	check_coefficients();
	remove_leading_zeros();
	return is;
}



polynomial< gf_element > &
polynomial< gf_element >::operator /= (const polynomial< gf_element > &a)
{
	polynomial< gf_element > r;
	div_rem(*this, r, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator /= (const gf_element &a)
{
	divide(*this, *this, a);
	return *this;
}



polynomial< gf_element > &
polynomial< gf_element >::operator %= (const polynomial< gf_element > &a)
{
	polynomial< gf_element > q;
	div_rem(q, *this, *this, a);
	return *this;
}



#if 0
polynomial< gf_element > gcd(const polynomial< gf_element > &aa,
			     const polynomial< gf_element > &bb)
{
	polynomial< gf_element > tmp;
	gcd(tmp, aa, bb);
	return tmp;
}



polynomial< gf_element > xgcd(polynomial< gf_element > &x,
			      polynomial< gf_element > &y,
			      const polynomial< gf_element > &aa,
			      const polynomial< gf_element > &bb)
{
	polynomial< gf_element > tmp;
	xgcd(tmp, x, y, aa, bb);
	return tmp;
}
#endif



void xgcd(polynomial< gf_element > &d, polynomial< gf_element > &x,
	  polynomial< gf_element > &y, const polynomial< gf_element > &aa,
	  const polynomial< gf_element > &bb)
{
	x.ffield = polynomial< gf_element >::common_field(aa.ffield, bb.ffield);
	polynomial< gf_element >::build_frame(x.ffield);
	xgcd(d.pol, x.pol, y.pol, aa.pol, bb.pol);
	polynomial< gf_element >::delete_frame();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
