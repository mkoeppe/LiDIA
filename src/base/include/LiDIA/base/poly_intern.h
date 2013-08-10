// -*- C++ -*-
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

//
// This  include  file  supports basic polynomial operations over
// all types
//
// Most   of   this   code (+, -, *, /, gcd, xgcd)  was
// originally  written by Victor Shoup (finite fields); Thomas Papanikolaou
// changed his implementation to use templates and the names of the
// functions to our standard. Then I extended this code to be able to
// treat polynomials over other types, especially over bigints
//
//                                             Stefan Neis
//


#ifndef LIDIA_POLY_INTERN_H_GUARD_
#define LIDIA_POLY_INTERN_H_GUARD_


#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_XERROR_H_GUARD_
# include	"LiDIA/xerror.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint;


//
// The first class in this file (base_polynomial)
// is for internal use only !!!
//
// The user is supposed to always use "polynomial< T >".
//

//
// debug defines / error defines
//
extern const char *PRT;
extern const char *poly_error_msg[];

#define DV_BP LDBL_UNIPOL
#define DM_BP "base_polynomial"
#define LP_ERROR poly_error_msg

//
// debug level
//
//   0 : constructors / destructor
//   1 : internal functions
//   2 : access functions
//   3 : assignments
//   4 : comparison
//   5 : arithmetic functions
//   6 : input / output
//

template <class T >
class base_polynomial;

template <class T >
std::istream& operator>>(std::istream&, base_polynomial<T>&);
template <class T >
std::ostream& operator<<(std::ostream&, base_polynomial<T> const &);

template< class T > class base_polynomial
// Type T must at least supply the functions add(c, a, b)
// subtract(c, a, b), negate(c, a) and multiply (c, a, b).
{
protected:

	//
	// base_polynomial< T > a = c[0] + c[1]*X + ... + c[deg]*X^deg
	//

	T * coeff;
	lidia_size_t deg;

	void copy_data(T * d, const T * vd, lidia_size_t al);

public:

	//
	// constructors and destructor
	//

	base_polynomial();
	base_polynomial(const T & x);
	base_polynomial(const T * v, lidia_size_t d);
	base_polynomial(const base_vector< T > v);
	base_polynomial(const base_polynomial< T > & p);
	~base_polynomial();


	// comparisons

	bool equal(const base_polynomial< T > &) const;


	//
	// access
	//

	T lead_coeff() const;
	T const_term() const;


	// no range checking --- use at your own risk
	T & operator[] (lidia_size_t i)
	{
		return coeff[i];
	}

	const T & operator[] (lidia_size_t i) const
	{
		return coeff[i];
	}


	T & at(lidia_size_t i)
	{
		debug_handler_l(DM_BP, "in member "
				"at (lidia_size_t)", DV_BP + 2);
		if ((i< 0) || (i > deg)) {
			lidia_error_handler_para(i, "i", "0 <= i <= degree",
						 "T & operator[] (lidia_size_t i) const",
						 DM_BP, LP_ERROR[2]);
#ifdef LIDIA_NO_EXIT_ON_ERROR
			return T();
#endif
		}
		return coeff[i];
	}

	const T & at(lidia_size_t i) const
	{
		debug_handler_l(DM_BP, "in member "
				"at (lidia_size_t)", DV_BP + 2);
		if ((i< 0) || (i > deg)) {
			lidia_error_handler_para(i, "i", "0 <= i <= degree",
						 "T & operator[] (lidia_size_t i) const",
						 DM_BP, LP_ERROR[2]);
#ifdef LIDIA_NO_EXIT_ON_ERROR
			return T();
#endif
		}
		return coeff[i];
	}



	//
	// member functions
	//

	void remove_leading_zeros();

	int set_data (const T * d, lidia_size_t l);
	T* get_data () const;

	lidia_size_t degree() const
	{
		debug_handler_l(DM_BP, "in member - function "
				"degree ()" , DV_BP + 2);
		return deg;
	}

	void set_degree(lidia_size_t);

	bool is_zero() const
	{
		return (deg < 0);
	}

	bool is_one() const
	{
		return (deg == 0 && coeff[0] == 1);
	}

	bool is_x() const
	{
		return (deg == 1 && coeff[1] == 1 && coeff[0] == 0);
	}

	//
	// assignment
	//

	void assign(const  T &);

	base_polynomial< T > & operator = (const T & a);

	void assign(const base_polynomial< T > &);

	base_polynomial< T > & operator = (const base_polynomial< T > & a);

	void assign_zero();
	void assign_one();
	void assign_x();

	void swap(base_polynomial< T > &);

	//
	// operator overloading
	//

	T operator() (const T & value) const;


	//
	// arithmetic procedures
	//

	void negate(const base_polynomial< T > &);
	void add(const base_polynomial< T > &, const base_polynomial< T > &);
	void add(const base_polynomial< T > &, const T & );
	void add(const T &, const base_polynomial< T > &); //mmmh ??
	void subtract(const base_polynomial< T > &, const base_polynomial< T > &);
	void subtract(const base_polynomial< T > &, const T &);
	void subtract(const T &, const base_polynomial< T > &);
	void multiply(const base_polynomial< T > &, const base_polynomial< T > &);
	void multiply(const base_polynomial< T > &, const T&);
	void multiply(const T &, const base_polynomial< T > &);
	void power(const base_polynomial< T > &, const bigint&);


public:
	//
	// functions
	//

	void derivative(const base_polynomial< T > &);


	//
	// input / output
	//

protected:
	void read(std::istream &);

public:
	void read_verbose(std::istream &);
	void print_verbose(std::ostream &, char v = 'x') const;
	friend std::istream & operator >> < T >(std::istream & s, base_polynomial< T > & c);
	friend std::ostream & operator << < T >(std::ostream & s, const base_polynomial< T > & a);


};



//
// Now the general class for polynomials:
//

template< class T >
class polynomial : public base_polynomial< T >
{
public:
	//
	// constructors and destructor
	//

	polynomial(): base_polynomial< T > ()
	{ }

	polynomial(T x): base_polynomial< T > (x)
	{ }

	polynomial(const T * v, lidia_size_t d): base_polynomial< T > (v, d)
	{ }

	polynomial(const base_vector< T > v): base_polynomial< T > (v)
	{ }

	polynomial(const base_polynomial< T > & p): base_polynomial< T > (p)
	{ }

	~polynomial()
	{ }

	polynomial< T > & operator = (const base_polynomial< T > & a)
	{
		base_polynomial< T >::assign(a);
		return *this;
	}
};




#undef DV_BP
#undef DM_BP
#undef LP_ERROR



template< class T >
inline base_polynomial< T > &
base_polynomial< T >::operator = (const base_polynomial< T > &a)
{
	assign(a);
	return *this;
}



template< class T >
inline base_polynomial< T > &
base_polynomial< T >::operator = (const T & a)
{
	assign(a);
	return *this;
}



//
// comparators
//

template< class T >
inline bool
operator == (const base_polynomial< T > & a, const base_polynomial< T > & b)
{
	return a.equal(b);
}



template< class T >
inline bool
operator != (const base_polynomial< T > & a, const base_polynomial< T > & b)
{
	return !a.equal(b);
}



template< class T >
inline void
swap(base_polynomial< T > & a, base_polynomial< T > & b)
{
	a.swap(b);
}



//
// accessors
//

template< class T >
inline T
base_polynomial< T >::lead_coeff() const
{
	if (deg < 0)
		return T(0);
	else
		return coeff[deg];
}



template< class T >
inline T
base_polynomial< T >::const_term() const
{
	if (deg < 0)
		return T(0);
	else
		return coeff[0];
}



template< class T >
inline T
lead_coeff(const base_polynomial< T > & a)
{
	return a.lead_coeff();
}



template< class T >
inline T
const_term(const base_polynomial< T > & a)
{
	return a.const_term();
}



template< class T >
inline void
negate(base_polynomial< T > & c, const base_polynomial< T > & a)
{
	c.negate(a);
}



template< class T >
inline void
add(base_polynomial< T > & c,
    const base_polynomial< T > & a,
    const base_polynomial< T > & b)
{
	c.add(a, b);
}



template< class T >
inline void
add(base_polynomial< T > & c,
    const base_polynomial< T > & a,
    const T & b)
{
	c.add(a, b);
}



template< class T >
inline void
add(base_polynomial< T > & c,
    const T & b,
    const base_polynomial< T > & a)
{
	c.add(a, b);
}



template< class T >
inline void
subtract(base_polynomial< T > & c,
	 const base_polynomial< T > &a,
	 const base_polynomial< T > & b)
{
	c.subtract(a, b);
}



template< class T >
inline void
subtract(base_polynomial< T > & c,
	 const base_polynomial< T > & a,
	 const T & b)
{
	c.subtract(a, b);
}



template< class T >
inline void
subtract(base_polynomial< T > & c,
	 const T & b,
	 const base_polynomial< T > & a)
{
	c.subtract(b, a);
}



template< class T >
inline void
multiply(base_polynomial< T > & c,
	 const base_polynomial< T > & a,
	 const base_polynomial< T > & b)
{
	c.multiply(a, b);
}



template< class T >
inline void
multiply(base_polynomial< T > & c,
	 const base_polynomial< T > & a,
	 const T & b)
{
	c.multiply(a, b);
}



template< class T >
inline void
multiply(base_polynomial< T > & c,
	 const T & b,
	 const base_polynomial< T > & a)
{
	c.multiply(b, a);
}



template< class T >
inline void
power(base_polynomial< T > & c,
      const base_polynomial< T > & a,
      const bigint & b)
{
	c.power(a, b);
}



template< class T >
inline base_polynomial< T >
operator - (const base_polynomial< T > & a)
{
	base_polynomial< T > c;

	c.negate(a);
	return c;
}



template< class T >
inline base_polynomial< T >
operator + (const base_polynomial< T > & a,
	    const base_polynomial< T > & b)
{
	base_polynomial< T > c;

	add(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator + (const base_polynomial< T > & a,
	    const T & b)
{
	base_polynomial< T > c;

	add(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator + (const T & b,
	    const base_polynomial< T > & a)
{
	base_polynomial< T > c;

	add(c, b, a);
	return c;
}



template< class T >
inline base_polynomial< T >
operator - (const base_polynomial< T > & a,
	    const base_polynomial< T > & b)
{
	base_polynomial< T > c;

	subtract(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator - (const base_polynomial< T > & a,
	    const T & b)
{
	base_polynomial< T > c;

	subtract(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator - (const T & a,
	    const base_polynomial< T > & b)
{
	base_polynomial< T > c;

	subtract(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator * (const base_polynomial< T > & a,
	    const base_polynomial< T > & b)
{
	base_polynomial< T > c;

	multiply(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator * (const base_polynomial< T > & a,
	    const T & b)
{
	base_polynomial< T > c;

	multiply(c, a, b);
	return c;
}



template< class T >
inline base_polynomial< T >
operator * (const T & b,
	    const base_polynomial< T > & a)
{
	base_polynomial< T > c;

	multiply(c, b, a);
	return c;
}



template< class T >
inline base_polynomial< T > &
operator += (base_polynomial< T > & a, const base_polynomial< T > & b)
{
	add(a, a, b);
	return a;
}



template< class T >
inline base_polynomial< T > &
operator += (base_polynomial< T > & a, const T & b)
{
	add(a, a, b);
	return a;
}



template< class T >
inline base_polynomial< T > &
operator -= (base_polynomial< T > & a, const base_polynomial< T > & b)
{
	subtract(a, a, b);
	return a;
}



template< class T >
base_polynomial< T > &
operator -= (base_polynomial< T > & a, const T & b)
{
	subtract(a, a, b);
	return a;
}



template< class T >
inline base_polynomial< T > &
operator *= (base_polynomial< T > & a, const base_polynomial< T > & b)
{
	multiply(a, a, b);
	return a;
}



template< class T >
inline base_polynomial< T > &
operator *= (base_polynomial< T > & a, const T & b)
{
	multiply(a, a, b);
	return a;
}



template< class T >
inline void
derivative(base_polynomial< T > & c,
	   const base_polynomial< T > & a)
{
	c.derivative(a);
}



template< class T >
inline base_polynomial< T >
derivative(const base_polynomial< T > & a)
{
	base_polynomial< T > c;

	c.derivative(a);
	return c;
}



template< class T >
inline std::istream &
operator >> (std::istream & s, base_polynomial< T > & c)
{
	c.read(s);
	return s;
}



template< class T >
inline std::ostream &
operator << (std::ostream & s, const base_polynomial< T > & a)
{
	a.print_verbose(s);
	return s;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/base/poly_intern.cc"
#endif



#endif	// LIDIA_POLY_INTERN_H_GUARD_
