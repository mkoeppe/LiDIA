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


#ifndef LIDIA_POLY_INTERN_CC_GUARD_
#define LIDIA_POLY_INTERN_CC_GUARD_



#ifndef LIDIA_POLY_INTERN_H_GUARD_
# include	"LiDIA/base/poly_intern.h"
#endif
#include	<cctype>



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define DV_BP LDBL_UNIPOL
#define DM_BP "base_polynomial"
#define DV_P LDBL_UNIPOL+12
#define DM_P "polynomial< T >"
#define LP_ERROR poly_error_msg

// BASE_POLYNOMIAL FUNCTIONS

//
// constructors and destructor
//

template< class T >
base_polynomial< T >::base_polynomial()
{
	debug_handler_l(DM_BP, "in constructor "
			"base_polynomial()", DV_BP);
	this->deg = -1;
	this->coeff = NULL;
}



template< class T >
base_polynomial< T >::base_polynomial(const T & x)
{
	debug_handler_l(DM_BP, "in constructor "
			"base_polynomial(const T &)", DV_BP);
	if (x == 0) {
		this->deg = -1;
		this->coeff = NULL;
	}
	else {
		this->deg = 0;
		this->coeff = new T[1];
		memory_handler(this->coeff, DM_BP, "base_polynomial(const T &) :: "
			       "Error in memory allocation (coeff)");
		this->coeff[0] = x;
	}
}



template< class T >
base_polynomial< T >::base_polynomial(const T * v, lidia_size_t d)
{
	debug_handler_l(DM_BP, "in constructor "
			"base_polynomial(const T *, lidia_size_t)", DV_BP);
	if (d < 0 || v == NULL)
		lidia_error_handler_para(d, "d", "d >= 0",
					 PRT, "v", "v != NULL",
					 "base_polynomial< T >::"
					 "base_polynomial(const T * v, lidia_size_t d)",
					 DM_BP, LP_ERROR[0]);
	this->deg = d;
	this->coeff = new T[d + 1];
	memory_handler(this->coeff, DM_BP, "base_polynomial(const T *, lidia_size_t)"
		       " :: Error in memory allocation (coeff)");
	this->copy_data(this->coeff, v, this->deg);
	this->remove_leading_zeros();
}



template< class T >
base_polynomial< T >::base_polynomial(const base_vector< T > v)
{
	debug_handler_l(DM_BP, "in constructor "
			"base_polynomial(const base_vector< T > )", DV_BP);
	this->deg = v.size()-1;
	this->coeff = v.get_data();   // Note: the ownership of the array
                            	// returned by v.get_data() is taken 
                                // by this->coeff!
	this->remove_leading_zeros();
}



template< class T >
base_polynomial< T >::base_polynomial(const base_polynomial< T > &p)
{
	debug_handler_l(DM_BP, "in constructor "
			"base_polynomial(const base_polynomial< T > )", DV_BP);
	this->deg = p.deg;
	if (this->deg < 0)
		this->coeff = NULL;
	else {
		this->coeff = new T[this->deg + 1];
		memory_handler(this->coeff, DM_BP, "base_polynomial(const T *, lidia_size_t)"
			       " :: Error in memory allocation (coeff)");
		this->copy_data(coeff, p.coeff, this->deg);
	}
}



template< class T >
base_polynomial< T >::~base_polynomial()
{
	debug_handler_l(DM_BP, "in destructor "
			"~base_polynomial()", DV_BP);
	if (this->deg >= 0) {
		delete[] this->coeff;
	}
}



//
// comparisons
//

template< class T >
bool base_polynomial< T >::equal(const base_polynomial< T > &b) const
{
	debug_handler_l(DM_BP, "in member - function "
			"bool equal(const base_polynomial< T > &)", DV_BP + 4);
	const T *ap, *bp;
	lidia_size_t i;
	if (this->deg != b.deg)
		return false;
	for (i = this->deg + 1, ap = this->coeff, bp = b.coeff; i; i--, ap++, bp++)
		if (*ap != *bp)
			return false;
	return true;
}



//
// member functions
//

template< class T >
int base_polynomial< T >::set_data (const T * d, lidia_size_t l)
{
	debug_handler_l(DM_BP, "in member - function "
			"set_data (const T *, lidia_size_t)" , DV_BP + 2);
	if (l <= 0 || d == NULL)
		lidia_error_handler_para(l, "l", "l > 0",
					 PRT, "d", "d != NULL",
					 "base_polynomial< T >::"
					 "set_data(const T * d, lidia_size_t l)",
					 DM_BP, LP_ERROR[0]);
	if (this->deg >= 0)
		delete[] this->coeff;
	this->coeff = new T[l];
	memory_handler(this->coeff, DM_BP, "set_data(const T *, lidia_size_t) :: "
		       "Error in memory allocation (this->coeff)");

	this->deg = l-1;

	for (register lidia_size_t i = 0; i < l; i++)
		this->coeff[i] = d[i];
	this->remove_leading_zeros();
	return 0;
}



template< class T >
T* base_polynomial< T >::get_data () const
{
	debug_handler_l(DM_BP, "in member - function "
			"get_data ()" , DV_BP + 2);

	T * d;

	if (this->deg < 0) {
		d = new T[1];
		memory_handler(d, DM_BP, "get_data () :: "
			       "Error in memory allocation (d)");
		d[0] = 0;
		return d;
	}

	d = new T [this->deg+1];
	memory_handler(d, DM_BP, "get_data () :: "
		       "Error in memory allocation (d)");

	for (register lidia_size_t i = 0; i <= this->deg; i++)
		d[i] = this->coeff[i];
	return d;
}



template< class T >
void base_polynomial< T >::copy_data(T * d, const T * vd, lidia_size_t al)
{
	debug_handler_l(DM_BP, "in member - function "
			"copy_data(T *, const T *, lidia_size_t)", DV_BP + 1);
	for (lidia_size_t i = al +1; i; i--, d++, vd++)
		(*d) = (*vd);
}



template< class T >
void base_polynomial< T >::remove_leading_zeros()
{
	debug_handler_c(DM_BP, "in member - function remove_leading_zeros ()",
			DV_BP + 1, std::cout << "original degree is " << deg << std::endl << std::flush);
	T c, *np;
	lidia_size_t d, i;

	d = this->deg;
	np = this->coeff + d;
#ifdef LIDIA_T_IS_BUILTIN
	while (d >= 0 && (*np) == 0)
#else
		while (d >= 0 && np->is_zero())
#endif
			d--, np--;

	if (d < 0) {
		this->deg = d;
		delete[] this->coeff;
		this->coeff = NULL;
	}
	else if (d != this->deg) {
		debug_handler_c(DM_BP, "in member - function remove_leading_zeros()",
				DV_BP + 1, std::cout << "new degree is " << d << std::endl << std::flush);
		this->deg = d;
		np = new T[d + 1];
		memory_handler(np, DM_BP, "remove_leading_zeros() :: "
			       "Error in memory allocation (np)");
		for (i = 0; i <= d; i++)
			np[i] = this->coeff[i];
		delete[] this->coeff;
		this->coeff = np;
	}
}



template< class T >
void base_polynomial< T >::set_degree(lidia_size_t d)
{
	debug_handler_l(DM_BP, "in member - function "
			"set_degree(lidia_size_t)", DV_BP + 2);

	if (d < -1)
		lidia_error_handler_para(d, "d", "d >= -1",
					 "void base_polynomial< T >::"
					 "set_degree(lidia_size_t d)",
					 DM_BP, LP_ERROR[2]);

	if (d == this->deg)
		return;

	if (d < 0) {
		delete[] this->coeff; // Note: coeff != NULL, since otherwise d == deg !
		this->coeff = NULL;
		this->deg = d;
		return;
	}

	T *tmp = this->coeff;
	this->coeff = new T [d + 1];
	memory_handler(this->coeff, DM_BP, "set_degree(lidia_size_t d) :: "
		       "Error in memory allocation (this->coeff)");

	register lidia_size_t minimum = (d < this->deg)? d : this->deg;

	for (register lidia_size_t i = 0; i <= minimum; i++)
		this->coeff[i] = tmp[i];

	if (tmp != NULL)
		delete[] tmp;
	this->deg = d;
}



//
// assignment
//

template< class T >
void base_polynomial< T >::assign(const T & a)
{
	debug_handler_l(DM_BP, "in member - function "
			"assign (const T &)", DV_BP + 3);
	if (a == 0) {
		if (this->deg >= 0) {
			delete[] this->coeff;
			this->coeff = NULL;
			this->deg = -1;
		}
		return;
	}
	if (this->deg > 0)
		delete[] this->coeff;
	if (this->deg != 0) {
		this->deg = 0;
		this->coeff = new T[1];
		memory_handler(this->coeff, DM_BP, "assign(const T &) :: "
			       "Error in memory allocation (coeff)");
	}
	this->coeff[0] = a;
}



template< class T >
void base_polynomial< T >::assign(const base_polynomial< T > &a)
{
	debug_handler_l(DM_BP, "in member - function "
			"assign (const base_polynomial< T > &)", DV_BP + 3);
	if (this->deg != a.deg) {
		if (this->deg >= 0)
			delete[] this->coeff;
		this->deg = a.deg;
		if (this->deg >= 0) {
			this->coeff = new T[this->deg + 1];
			memory_handler(this->coeff, DM_BP, "assign(const base_polynomial< T > &) :: "
				       "Error in memory allocation (coeff)");
		}
		else this->coeff = NULL;
	}
	for (register lidia_size_t i = 0; i <= this->deg; i++)
		this->coeff[i] = a.coeff[i];
}



template< class T >
void base_polynomial< T >::assign_zero()
{
	debug_handler_l(DM_BP, "in member - function "
			"assign_zero ()" , DV_BP + 3);
	if (this->deg >= 0) {
		this->deg = -1;
		delete[] this->coeff;
		this->coeff = NULL;
	}
}



template< class T >
void base_polynomial< T >::assign_one()
{
	debug_handler_l(DM_BP, "in member - function "
			"assign_one ()" , DV_BP + 3);
	if (this->deg > 0)
		delete[] this->coeff;
	if (this->deg != 0) {
		this->deg = 0;
		this->coeff = new T[1];
		memory_handler(this->coeff, DM_BP, "base_polynomial(const T &) :: "
			       "Error in memory allocation (coeff)");
	}
	this->coeff[0] = 1;
}



template< class T >
void base_polynomial< T >::assign_x()
{
	debug_handler_l(DM_BP, "in member - function "
			"assign_x ()" , DV_BP + 3);
	if (this->deg > 1 || this->deg == 0)
		delete[] this->coeff;
	if (this->deg != 1) {
		this->deg = 1;
		this->coeff = new T[2];
		memory_handler(this->coeff, DM_BP, "base_polynomial(const T &) :: "
			       "Error in memory allocation (coeff)");
	}
	this->coeff[0] = 0;
	this->coeff[1] = 1;
}



template< class T >
void base_polynomial< T >::swap(base_polynomial< T > &b)
{
	debug_handler_l(DM_BP, "in member - function "
			"swap(const base_polynomial< T > &)" , DV_BP + 3);
	T* tmp = this->coeff;
	this->coeff = b.coeff;
	b.coeff = tmp;
	lidia_size_t deg_tmp = this->deg;
	this->deg = b.deg;
	b.deg = deg_tmp;
}



//
// operator overloading
//

template< class T >
T base_polynomial< T >::operator() (const T & value) const
{
	debug_handler_l(DM_BP, "in operator "
			"() (const T &)", DV_BP + 5);

	if (this->deg == 0)
		return this->coeff[0];

	T result;
	result = 0;

	if (this->deg < 0)
		return result;

	result = this->coeff[deg];
	for (lidia_size_t i = this->deg - 1; i >= 0; i--) {
		LiDIA::multiply(result, result, value);
		LiDIA::add(result, result, this->coeff[i]);
	}
	return result;
}



//
// arithmetic procedures
//

template< class T >
void base_polynomial< T >::negate (const base_polynomial< T > &a)
{
	debug_handler_l(DM_BP, "in member - function "
			"negate (const base_polynomial< T > &)", DV_BP + 5);
	register lidia_size_t d = a.deg;
	this->set_degree(d);

	const T * ap = a.coeff;
	T * cp = this->coeff;

	for (d++; d; d--, ap++, cp++)
		LiDIA::negate(*cp, *ap);
}



template< class T >
void base_polynomial< T >::add(const base_polynomial< T > & a,
				const base_polynomial< T > & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"add (const base_polynomial< T > &, "
			"const base_polynomial< T > &)", DV_BP + 5);
	const T *ap, *bp;
	T *cp;
	lidia_size_t deg_a = a.deg, deg_b = b.deg;
	lidia_size_t i, min_deg_ab, max_deg_ab;

	if (deg_a > deg_b) {
		max_deg_ab = deg_a;
		min_deg_ab = deg_b;
	}
	else {
		max_deg_ab = deg_b;
		min_deg_ab = deg_a;
	}

	this->set_degree(max_deg_ab);
	if (max_deg_ab < 0) return;

	ap = a.coeff;
	bp = b.coeff;
	cp = this->coeff;

	for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
		LiDIA::add(*cp, *ap, *bp);

	if (deg_a > min_deg_ab)
		for (i = deg_a - min_deg_ab; i; i--, cp++, ap++)
			*cp = *ap;
	else if (deg_b > min_deg_ab)
		for (i = deg_b - min_deg_ab; i; i--, cp++, bp++)
			*cp = *bp;
	else
		this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::add(const base_polynomial< T > & a, const T & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"add (base_polynomial< T > &, "
			"const T &)", DV_BP + 5);
	if (a.deg < 0) {
		if (b != 0) {
			this->set_degree(0);
			this->coeff[0] = b;
			return;
		}
		this->set_degree(-1);
		return;
	}
	this->set_degree(a.deg);

	const T *ap = a.coeff;
	T *cp = this->coeff;

	LiDIA::add(*cp , *ap, b);
	if (a.deg > 0 && this != &a)
		for (register lidia_size_t i = a.deg; i; i--)
			*(++cp) = *(++ap);
	else
		this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::add(const T & b,
				const base_polynomial< T > & a)
{
	debug_handler_l(DM_BP, "in member - function "
			"add (const T &, "
			"const base_polynomial< T > &)", DV_BP + 5);
	if (a.deg < 0) {
		if (b != 0) {
			this->set_degree(0);
			this->coeff[0] = b;
			return;
		}
		this->set_degree(-1);
		return;
	}
	this->set_degree(a.deg);

	const T *ap = a.coeff;
	T *cp = this->coeff;

	LiDIA::add(*cp, *ap, b);
	if (a.deg > 0 && &a != this)
		for (register lidia_size_t i = a.deg; i; i--)
			*(++cp) = *(++ap);
	else
		this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::subtract(const base_polynomial< T > & a,
				     const base_polynomial< T > & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"subtract (const base_polynomial< T > &, "
			"const base_polynomial< T > &)", DV_BP + 5);
	const T *ap, *bp;
	T *cp;
	lidia_size_t deg_a = a.deg, deg_b = b.deg;
	lidia_size_t i, min_deg_ab, max_deg_ab;

	if (deg_a > deg_b) {
		max_deg_ab = deg_a;
		min_deg_ab = deg_b;
	}
	else {
		max_deg_ab = deg_b;
		min_deg_ab = deg_a;
	}

	this->set_degree(max_deg_ab);
	if (max_deg_ab < 0) return;

	ap = a.coeff;
	bp = b.coeff;
	cp = this->coeff;
	for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
		LiDIA::subtract(*cp, *ap, *bp);

	if (deg_a > min_deg_ab && this != &a)
		for (i = deg_a - min_deg_ab; i; i--, cp++, ap++)
			*cp = *ap;
	else if (deg_b > min_deg_ab)
		for (i = deg_b - min_deg_ab; i; i--, cp++, bp++)
			LiDIA::negate(*cp, *bp);
	else
		this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::subtract(const base_polynomial< T > & a,
				     const T & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"subtract (const base_polynomial< T > &, "
			"const T &)", DV_BP + 5);
	if (a.deg < 0) {
		if (b != 0) {
			this->set_degree(0);
			this->coeff[0] = - b;
			return;
		}
		this->set_degree(-1);
		return;
	}
	this->set_degree(a.deg);

	const T *ap = a.coeff;
	T *cp = this->coeff;

	LiDIA::subtract(*cp, *ap, b);
	if (a.deg > 0 && this != &a)
		for (register lidia_size_t i = a.deg; i; i--)
			*(++cp) = *(++ap);
	else
		this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::subtract(const T & b,
				     const base_polynomial< T > & a)
{
	debug_handler_l(DM_BP, "in member - function "
			"subtract (const T &, "
			"const base_polynomial< T > &)", DV_BP + 5);
	if (a.deg < 0) {
		if (b != 0) {
			this->set_degree(0);
			this->coeff[0] = b;
			return;
		}
		this->set_degree(-1);
		return;
	}
	this->set_degree(a.deg);

	const T *ap = a.coeff;
	T *cp = this->coeff;

	LiDIA::subtract(*cp, b, *ap);
	if (a.deg > 0)
		for (register lidia_size_t i = a.deg; i; i--)
			LiDIA::negate(*(++cp), *(++ap));
	else
		this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::multiply(const base_polynomial< T > & a,
				     const base_polynomial< T > & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"multiply (const base_polynomial< T > &, "
			"const base_polynomial< T > &)", DV_BP + 5);
	const T *ap, *bp;
	T * cp, temp;
	lidia_size_t deg_a = a.deg, deg_b = b.deg;

	if (deg_a < 0 || deg_b < 0) {
		this->set_degree(-1);
		return;
	}

	lidia_size_t i, j, deg_ab = deg_a + deg_b;

	if (this->coeff == a.coeff || this->coeff == b.coeff) {
		base_polynomial< T > c_temp;
		c_temp.set_degree(deg_ab);

		for (i = deg_ab + 1, cp = c_temp.coeff; i; i--, cp++)
			(*cp) = 0;

		for (i = 0, ap = a.coeff; i <= deg_a; i++, ap++)
			for (j = deg_b + 1, bp = b.coeff, cp = c_temp.coeff + i;
			     j; j--, bp++, cp++) {
				LiDIA::multiply(temp, *ap, *bp);
				LiDIA::add(*cp, *cp, temp);
			}
		c_temp.remove_leading_zeros();
		assign(c_temp);
	}
	else {
		this->set_degree(deg_ab);

		for (i = deg_ab + 1, cp = this->coeff; i; i--, cp++)
			(*cp) = 0;

		for (i = 0, ap = a.coeff; i <= deg_a; i++, ap++)
			for (j = deg_b + 1, bp = b.coeff, cp = this->coeff + i;
			     j; j--, bp++, cp++) {
				LiDIA::multiply(temp, *ap, *bp);
				LiDIA::add(*cp, *cp, temp);
			}
		this->remove_leading_zeros();
	}
}



template< class T >
void base_polynomial< T >::multiply(const base_polynomial< T > & a,
				     const T & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"multiply (const base_polynomial< T > &, "
			"const T &)", DV_BP + 5);
	if (b == 0) {
		this->set_degree(-1);
		return;
	}
	const T *ap;
	T *cp;
	lidia_size_t deg_a = a.deg;

	this->set_degree(deg_a);

	register lidia_size_t i = 0;
	for (ap = a.coeff, cp = this->coeff; i <= deg_a; i++, ap++, cp++)
		LiDIA::multiply(*cp, *ap, b);
	this->remove_leading_zeros(); // necessary if characteristic != 0
}



template< class T >
void base_polynomial< T >::multiply(const T & b,
				     const base_polynomial< T > & a)
{
	debug_handler_l(DM_BP, "in member - function "
			"multiply (const T &, "
			"const base_polynomial< T > &)", DV_BP + 5);
	if (b == 0) {
		this->set_degree(-1);
		return;
	}
	const T *ap;
	T *cp;
	lidia_size_t deg_a = a.deg;

	this->set_degree(deg_a);

	register lidia_size_t i = 0;
	for (ap = a.coeff, cp = this->coeff; i <= deg_a; i++, ap++, cp++)
		LiDIA::multiply(*cp, b, *ap);
	this->remove_leading_zeros(); // necessary if characteristic != 0
}



template< class T >
void base_polynomial< T >::power(const base_polynomial< T > & a,
				  const bigint & b)
{
	debug_handler_l(DM_BP, "in member - function "
			"power (const base_polynomial< T > &, "
			"const bigint &)", DV_BP + 5);
	bigint exponent;
	base_polynomial< T > multiplier;
	if (b.is_negative())
		assign_zero();
	else if (b.is_zero() || a.is_one())
		assign_one();
	else {
		exponent.assign(b);
		multiplier.assign(a);
		assign_one();
		while (exponent.is_gt_zero()) {
			if (!exponent.is_even())
				LiDIA::multiply(*this, *this, multiplier);
			LiDIA::multiply(multiplier, multiplier, multiplier);
			exponent.divide_by_2();
		}
	}
}



//
// functions
//

template< class T >
void base_polynomial< T >::derivative(const base_polynomial< T > & a)
{
	debug_handler_l(DM_BP, "in member - function "
			"derivative (const base_polynomial< T > &)", DV_BP + 5);

	lidia_size_t d = a.deg;

	if (d <= 0) {
		this->set_degree(-1);
		return;
	}

	this->set_degree(d - 1);
	const T *ap = a.coeff + 1;
	T *cp = this->coeff;
	T temp;
	for (lidia_size_t i = 1; i <= d; i++, cp++, ap++) {
		temp = i; // necessary, since bigcomplex does not
                                // support automatic cast !!
		LiDIA::multiply(*cp, *ap, temp);
	}
	this->remove_leading_zeros(); // necessary if characteristic != 0
}



//
// input / output
//

template< class T >
void base_polynomial< T >::read(std::istream &s)
{
	debug_handler_l(DM_BP, "in member - function "
			"read(std::istream &)", DV_BP + 6);
	char c;
	s >> std::ws >> c;
	s.putback(c);

	if (c != '[')
		read_verbose(s);
	else {
		base_vector< T > temp;
		s >> temp;
		temp.reverse();

		this->set_degree(temp.size() - 1);
		if (this->deg >= 0) {
		        // Note: the ownership of the array
                        // returned by temp.get_data() is taken 
                        // by this->coeff!
		        this->coeff = temp.get_data();
		}
	}
	this->remove_leading_zeros();
}



template< class T >
void base_polynomial< T >::read_verbose(std::istream & s)
{
	// This function reads a univariate polynomial in any variable.
	// input format : a_n*x^n+ ... + a_1*x + a_0
	// e.g. :   3*x^2 + 4*x - 1
	// Monomials need not be sorted, and powers of x may even appear
	// repeated, '*' may be omitted and coefficients may follow the variable:
	//        -4 + 8x^5 + 2 - x^2 3 + x^5 + x^5*17
	// Note however, that the routine will work faster, if the leading monomial
	// is read first.

	debug_handler_l(DM_BP, "in member-function "
			"void read_verbose(std::istream &)", DV_BP + 6);

	register lidia_size_t n = 0;
	lidia_size_t sz;
	char c;

	this->set_degree(8);
	for (; n <= this->deg; n++)
		this->coeff[n] = 0;

	char variable = 0;
	T coeff_tmp;
	coeff_tmp = 1;
	T tmp;

  // Read a monomial, i.e. "x^k" or "- x^k"
  // or "a*x^k" or "a x^k" or "x^k*a" or "x^k a"

	do {
		c = s.get();
	} while (isspace(c) && c != '\n');
	while (c != '\n' && c != EOF && s.good()) {
		sz = 0; // Assume we read coeffizient of x^0;
		if (c == '+') {
			coeff_tmp = 1;
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
		}
		if (c == '-') {
			coeff_tmp = -1;
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
		}
#ifdef POLYREAD_DEBUG
		std::cout << "\n 1) Now looking at " << c;
#endif
		if (c >= '0' && c <= '9' || c == '(') {
			s.putback(c);
			s >> tmp;
			coeff_tmp *= tmp;
			do {
				c = s.get();
#ifdef POLYREAD_DEBUG
				std::cout << ", looking at " << c;
#endif
			} while (isspace(c) && c != '\n');
			if (c == '*')
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
#ifdef POLYREAD_DEBUG
			std::cout << "\n coeff_tmp is now " << coeff_tmp;
			std::cout << ", looking at " << c;
#endif
		}
#ifdef POLYREAD_DEBUG
		std::cout << "\n 2) Now looking at " << c;
#endif
		if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
			if (variable == 0)
				variable = c;
			else if (variable != c)
				lidia_error_handler_c("base_polynomial", "member function "
						      "read_verbose: The given string is not "
						      "recognized to be a univariate polynomial",
						      std::cout << "Variable name seemed to be " << variable;
						      std::cout << " and now you used " << c << "." << std::endl);
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
#ifdef POLYREAD_DEBUG
			std::cout << "\n 3) Now looking at " << c;
#endif

			if (c != '^') sz = 1;
			else {
				s >> sz;
#ifdef POLYREAD_DEBUG
				std::cout << "sz ist " << sz;
#endif
				do {
					c = s.get();
#ifdef POLYREAD_DEBUG
					std::cout << "\n4') Looking at " << c;
#endif
				} while (isspace(c) && c != '\n');
#ifdef POLYREAD_DEBUG
				std::cout << "\n 5) Now looking at " << c;
#endif
			}
			if (c == '*') {
				s >> tmp;
				coeff_tmp *= tmp;
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
#ifdef POLYREAD_DEBUG
			std::cout << "\n 6) Now looking at " << c;
#endif

			if (c >= '0' && c <= '9' || c == '(') {
				s.putback(c);
				s >> tmp;
#ifdef POLYREAD_DEBUG
				std::cout << "\n Old coeff_tmp: " << coeff_tmp;
#endif
				coeff_tmp *= tmp;
#ifdef POLYREAD_DEBUG
				std::cout << "\n New coeff_tmp: " << coeff_tmp;
#endif
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
		}

		if (c != '+' && c != '-') {
			// No next monomial, so assume end of input is reached
			s.putback(c);
			c = '\n'; // set c to end--marker
		}
		if (sz >= n) {
			this->set_degree(sz);
			for (; n <= this->deg; n++)
				this->coeff[n] = 0;
		}
		this->coeff[sz] += coeff_tmp;

#ifdef POLYREAD_DEBUG
		std::cout << "\nSuccessfully read next monomial; Poly now is " << (*this);
		std::cout << "\n Now looking at " << c;
#endif
	}
	this->remove_leading_zeros();
}



// print polynomial
template< class T >
void base_polynomial< T >::print_verbose(std::ostream &s, char x) const
{
	debug_handler_l(DM_BP, "in member - function "
			"print_verbose(std::ostream &, char)", DV_BP + 6);
	lidia_size_t d = this->deg;

	if (d < 0)
		s << 0;
	else if (d == 0)
		s << this->coeff[0];
	else if (d == 1) {
		if (this->coeff[1] == 1)
			s << x;
		else
			s << this->coeff[1] << " * " << x;
		if (this->coeff[0] != 0)
			s << "+ " << this->coeff[0];
	}
	else {
		if (this->coeff[d] == 1)
			s << x << "^" << d;
		else
			s << this->coeff[d] << " * " << x << "^" << d;
		for (register lidia_size_t i = d - 1; i > 1; i--)
			if (this->coeff[i] == 1)
				s << " + " << x << "^" << i;
			else if (this->coeff[i] != 0)
				s << " + " << this->coeff[i] << " * " << x << "^" << i;
		if (this->coeff[1] == 1)
			s << " + " << x;
		else if (this->coeff[1] != 0)
			s << " + " << this->coeff[1] << " * " << x;
		if (this->coeff[0] != 0)
			s << " + " << this->coeff[0];
	}
}



#undef DV_BP
#undef DM_BP
#undef DV_P
#undef DM_P
#undef LP_ERROR



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_POLY_INTERN_CC_GUARD_
