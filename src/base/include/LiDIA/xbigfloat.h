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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_XBIGFLOAT_H_GUARD_
#define LIDIA_XBIGFLOAT_H_GUARD_



#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class xdouble;
class bigint;



class xbigfloat
{
private:

	bigfloat x;



public:


	//
	// c'tors and d'tor
	//

	xbigfloat();
	xbigfloat(int i);
	xbigfloat(long l);
	xbigfloat(unsigned long ul);
	xbigfloat(double d);
	xbigfloat(const xdouble & xd);
	xbigfloat(const bigint & bi);
	xbigfloat(const bigfloat & bf);
	xbigfloat(const xbigfloat & xbf);
	~xbigfloat();



	//
	// accessors
	//

	const bigfloat & get_x() const;

	long get_exponent() const;
	long exponent() const;
	const bigint & get_mantissa() const;
	const bigint & mantissa() const;
	long b_value() const;



	//
	// assigners
	//

	void assign_zero();
	void assign_one();

	void assign(int i);
	void assign(long l);
	void assign(unsigned long ul);
	void assign(double d);
	void assign(const xdouble & xd);
	void assign(const bigint & bi);
	void assign(const bigfloat & bf);
	void assign(const xbigfloat & bf);

	xbigfloat & operator = (int i);
	xbigfloat & operator = (long l);
	xbigfloat & operator = (unsigned long ul);
	xbigfloat & operator = (double d);
	xbigfloat & operator = (const xdouble & xd);
	xbigfloat & operator = (const bigint & bi);
	xbigfloat & operator = (const bigfloat & bf);
	xbigfloat & operator = (const xbigfloat & xbf);



	//
	// comparators
	//

	bool is_zero() const;
	bool is_one() const;

	bool is_negative() const;
	bool is_positive() const;

	int get_sign() const;

	int sign() const;

	int compare(const xbigfloat & x) const;



	//
	// modifiers
	//

	void inc();
	void dec();

	void multiply_by_2();
	void divide_by_2();

	void negate();
	void absolute_value();
	void truncate(lidia_size_t k);

	void randomize(const bigint & n, long f);

	void swap (xbigfloat & a);



	//
	// arithmetic procedures
	//

	friend void add(xbigfloat & c, const xbigfloat & a, const xbigfloat & b);
	friend void subtract(xbigfloat & c, const xbigfloat & a, const xbigfloat & b);
	friend void multiply(xbigfloat & c, const xbigfloat & a, const xbigfloat & b);
	friend void square(xbigfloat & c, const xbigfloat & a);
	friend void divide(xbigfloat & c, const xbigfloat & a, const xbigfloat & b, long k);
	friend void divide (bigint & q, const xbigfloat & a, const xbigfloat & b);

	friend void shift_left(xbigfloat & c, const xbigfloat & a, long i);
	friend void shift_right(xbigfloat & c, const xbigfloat & a, long i);



	//
	// procedures
	//

	friend void sqrt(xbigfloat & c, const xbigfloat & a, long k);
	friend void exp(xbigfloat & c, const xbigfloat & a, long k);
	friend void log(xbigfloat & c, const xbigfloat & a, long k);

	friend void truncate(xbigfloat & c, const xbigfloat & a, lidia_size_t k);
	friend long b_value(const xbigfloat & a);
	friend void ceil(bigint & b, const xbigfloat & a);
	friend void floor(bigint & b, const xbigfloat & a);


	//
	// Verification of approximations
	//

	static bool check_relative_error(
		const xbigfloat & x,
		const xbigfloat & y,
		long k,
		long c);

	static bool check_absolute_error(
		const xbigfloat & x,
		const xbigfloat & y,
		long k,
		long c);

	//
	// I/O
	//

	friend std::ostream & operator << (std::ostream & out, const xbigfloat & a);
	friend std::istream & operator >> (std::istream & in , xbigfloat & a);

	void print_as_bigfloat() const;
	void print_as_bigfloat(std::ostream & out) const;
};



//
// c'tors and d'tor
//

inline
xbigfloat::xbigfloat ()
{
	assign_zero();
}



inline
xbigfloat::xbigfloat (int i)
{
	assign(i);
}



inline
xbigfloat::xbigfloat (long l)
{
	assign(l);
}



inline
xbigfloat::xbigfloat (unsigned long ul)
{
	assign(ul);
}



inline
xbigfloat::xbigfloat (double d)
{
	assign(d);
}



inline
xbigfloat::xbigfloat (const xdouble & xd)
{
	assign(xd);
}



inline
xbigfloat::xbigfloat (const bigint & bi)
{
	assign(bi);
}



inline
xbigfloat::xbigfloat (const bigfloat & bf)
{
	assign(bf);
}



inline
xbigfloat::xbigfloat (const xbigfloat & xbf)
{
	assign(xbf);
}



inline
xbigfloat::~xbigfloat()
{
	// nothing to do
}



//
// accessors
//

inline const bigfloat &
xbigfloat::get_x () const
{
	return x;
}



inline long
xbigfloat::get_exponent () const
{
	return (x.exponent() + x.mantissa().bit_length());
}



inline long
xbigfloat::exponent () const
{
	return (x.exponent() + x.mantissa().bit_length());
}



inline const bigint &
xbigfloat::get_mantissa () const
{
	return x.mantissa();
}



inline const bigint &
xbigfloat::mantissa () const
{
	return x.mantissa();
}



inline long
xbigfloat::b_value () const
{
	return (x.exponent() + x.mantissa().bit_length());
}



//
// assigners
//

inline xbigfloat &
xbigfloat::operator = (int i)
{
	assign(i);
	return *this;
}



inline xbigfloat &
xbigfloat::operator = (long l)
{
	assign(l);
	return *this;
}



inline xbigfloat &
xbigfloat::operator = (unsigned long ul)
{
	assign(ul);
	return *this;
}



inline xbigfloat &
xbigfloat::operator = (double d)
{
	assign(d);
	return *this;
}



inline xbigfloat &
xbigfloat::operator = (const xdouble & xd)
{
	assign(xd);
	return *this;
}



inline xbigfloat &
xbigfloat::operator = (const bigint & bi)
{
	assign(bi);
	return *this;
}



#if 0
inline xbigfloat &
xbigfloat::operator = (const bigrational & br)
{
	assign(br);
	return *this;
}
#endif



inline xbigfloat &
xbigfloat::operator = (const bigfloat & bf)
{
	assign(bf);
	return *this;
}



inline xbigfloat &
xbigfloat::operator = (const xbigfloat & xbf)
{
	if (&xbf != this) {
		assign(xbf);
	}
	return *this;
}



//
// comparators
//

inline bool
xbigfloat::is_zero () const
{
	return x.mantissa().is_zero();
}



inline bool
xbigfloat::is_positive () const
{
	return x.mantissa().is_positive();
}



inline bool
xbigfloat::is_negative () const
{
	return x.mantissa().is_negative();
}



inline int
xbigfloat::get_sign () const
{
	return x.sign();
}



inline int
xbigfloat::sign () const
{
	return x.sign();
}



inline bool
operator == (const xbigfloat & a, const xbigfloat & b)
{
	return (a.compare(b) == 0);
}



inline bool
operator != (const xbigfloat & a, const xbigfloat & b)
{
	return (a.compare(b) != 0);
}



inline bool
operator < (const xbigfloat & a, const xbigfloat & b)
{
	return (a.compare(b) < 0);
}



inline bool
operator <= (const xbigfloat & a, const xbigfloat & b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator > (const xbigfloat & a, const xbigfloat & b)
{
	return (a.compare(b) > 0);
}



inline bool
operator >= (const xbigfloat & a, const xbigfloat & b)
{
	return (a.compare(b) >= 0);
}



//
// modifiers
//

inline void
xbigfloat::multiply_by_2 ()
{
	x.multiply_by_2();
}



inline void
xbigfloat::divide_by_2 ()
{
	x.divide_by_2();
}



inline void
xbigfloat::negate ()
{
	x.negate();
}



inline void
xbigfloat::absolute_value ()
{
	x.absolute_value();
}



inline void
xbigfloat::truncate (lidia_size_t k)
{
	x.cut(k);
}



inline void
xbigfloat::randomize (const bigint & n, long f)
{
	x.randomize(n, f);
}



inline void
xbigfloat::swap (xbigfloat & a)
{
	x.swap(a.x);
}



//
// arithmetic procedures
//

inline void
swap (xbigfloat & a, xbigfloat & b)
{
	a.swap(b);
}



inline void
sqrt (xbigfloat & c, const xbigfloat & a, long k)
{
	sqrt(c.x, a.x, k);
}



inline void
exp (xbigfloat & c, const xbigfloat & a, long k)
{
	exp(c.x, a.x, k);
}



inline void
log (xbigfloat & c, const xbigfloat & a, long k)
{
	log_absolute(c.x, a.x, k);
}



inline void
truncate (xbigfloat & c, const xbigfloat & a, lidia_size_t k)
{
	c.x.cut(a.x, k);
}



inline long
b_value (const xbigfloat & a)
{
	return a.b_value();
}



//
// arithmetic procedures
//

inline void
inc (xbigfloat & a)
{
	a.inc();
}



inline void
dec (xbigfloat & a)
{
	a.dec();
}



void add(xbigfloat & c, const xbigfloat & a, const xbigfloat & b);
void subtract(xbigfloat & c, const xbigfloat & a, const xbigfloat & b);
void multiply(xbigfloat & c, const xbigfloat & a, const xbigfloat & b);
void square(xbigfloat & c, const xbigfloat & a);
void divide(xbigfloat & c, const xbigfloat & a, const xbigfloat & b, long k);
void divide (bigint & q, const xbigfloat & a, const xbigfloat & b);

void shift_left(xbigfloat & c, const xbigfloat & a, long i);
void shift_right(xbigfloat & c, const xbigfloat & a, long i);



//
// arithmetic operators
//

inline xbigfloat
operator - (const xbigfloat & a)
{
	xbigfloat b = a;

	b.negate();
	return b;
}



inline xbigfloat
operator + (const xbigfloat & a, const xbigfloat & b)
{
	xbigfloat c;

	add(c, a, b);
	return c;
}



inline xbigfloat
operator - (const xbigfloat & a, const xbigfloat & b)
{
	xbigfloat c;

	subtract(c, a, b);
	return c;
}



inline xbigfloat
operator * (const xbigfloat & a, const xbigfloat & b)
{
	xbigfloat c;

	multiply(c, a, b);
	return c;
}



#if 0
inline xbigfloat
operator / (const xbigfloat & a, const xbigfloat & b)
{
	xbigfloat c;

	divide(c, a, b);
	return c;
}
#endif



inline xbigfloat
operator << (const xbigfloat & a, long u)
{
	xbigfloat c;

	shift_left(c, a, u);
	return c;
}



inline xbigfloat
operator >> (const xbigfloat & a, long u)
{
	xbigfloat c;

	shift_right(c, a, u);
	return c;
}



inline xbigfloat &
operator += (xbigfloat & a, const xbigfloat & b)
{
	add(a, a, b);
	return a;
}



inline xbigfloat &
operator -= (xbigfloat & a, const xbigfloat & b)
{
	subtract(a, a, b);
	return a;
}



inline xbigfloat &
operator *= (xbigfloat & a, const xbigfloat & b)
{
	multiply(a, a, b);
	return a;
}



#if 0
inline xbigfloat &
operator /= (xbigfloat & a, const xbigfloat & b)
{
	divide(a, a, b);
	return a;
}
#endif



inline xbigfloat &
operator <<= (xbigfloat & a, long ui)
{
	shift_left(a, a, ui);
	return a;
}



inline xbigfloat &
operator >>= (xbigfloat & a, long ui)
{
	shift_right(a, a, ui);
	return a;
}



//
//
//


//
// I/O
//

std::ostream & operator << (std::ostream & out, const xbigfloat & a);
std::istream & operator >> (std::istream & in , xbigfloat & a);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_XBIGFLOAT_H_GUARD_
