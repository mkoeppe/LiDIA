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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGFLOAT_H_GUARD_
#define LIDIA_BIGFLOAT_H_GUARD_



//#define BIGFLOAT_DEBUG

#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#include	"LiDIA/bigfloat_config.h"
#ifndef LIDIA_XDOUBLE_H_GUARD_
# include	"LiDIA/xdouble.h"
#endif
#ifndef LIDIA_ARITH_INL_GUARD_
# include	"LiDIA/arith.inl"
#endif
#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
# include	"LiDIA/random_generator.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Rounding Modes

enum rnd_mode_t {
	MP_TRUNC =	0,	// round to zero
	MP_RND =	1,	// round to nearest
	MP_RND_UP =	2,	// round to +infinity
	MP_RND_DOWN =	3,	// round to -infinity
	MP_EXACT =	4	// do not round
};



//
// A bigfloat x is represented by a pair (bigint m, long e)
// This pair denotes the number
//
//               x = m * pow(2, e)
//
// where m is a bigint of size prec bits. In order to
// be able to work with exceptional values, each bigfloat
// contains a flag indicating that it is a not exact number,
// exact number (integer), +0, -0, +oo, -oo, overflow or underflow.
//
//

enum bigfloat_value {
	NotExact = 0,
	Exact = 1,
	PlusZero = 2,
	MinusZero = 4,
	PlusInf = 8,
	MinusInf = 16,
	Nan = 32,
	Overflow = 64,
	Underflow = 128
};



class bigfloat_flag
{
	bigfloat_value v;

public:

	bigfloat_flag()
		: v(NotExact)
	{
	}
	bigfloat_flag(const bigfloat_flag & x)
		: v(x.v)
	{
	}

	~bigfloat_flag()
	{
	}

#ifndef HEADBANGER

	int is_exact() const
	{
		return (v & Exact);
	}
	int is_not_exact() const
	{
		return !(v & Exact);
	}
	int is_plus_zero() const
	{
		return (v & PlusZero);
	}
	int is_minus_zero() const
	{
		return (v & MinusZero);
	}
	int is_pm_zero() const
	{
		return ((v & PlusZero) || (v & MinusZero));
	}
	int is_plus_infinity() const
	{
		return (v & PlusInf);
	}
	int is_minus_infinity() const
	{
		return (v & MinusInf);
	}
	int is_nan() const
	{
		return (v & Nan);
	}
	int is_overflow() const
	{
		return (v & Overflow);
	}
	int is_underflow() const
	{
		return (v & Underflow);
	}

	const bigfloat_value & flag() const
	{
		return v;
	}

	bigfloat_value & flag()
	{
		return v;
	}

#endif	// HEADBANGER

};



class bigfloat : public bigfloat_flag
{
private:

	bigint m;
	long e;
	long prec;



	//
	// These are the global variables for every bigfloat
	//

	static bigfloat C_Pi;
	static bigfloat C_E;
	static bigfloat C_Euler;
	static bigfloat C_Catalan;

	static long binary_precision;
	static long decimal_precision;
	static long digit_precision;
	static int  rounding_mode;

	static int  old_rnd_mode;
	static long old_bin_prec;


	//
	// Private functions used for normalisation
	//

	int leading_zeros() const;
	void base_digit_normalize();
	void normalize();

public:

	void cut(const bigfloat &x, long l);
	void cut(long l);

public:

	//
	// c'tors and d'tor
	//

	bigfloat();
	bigfloat(int i);
	bigfloat(long l);
	bigfloat(unsigned long ul);
	bigfloat(double d);
	bigfloat(const xdouble & xd);
	bigfloat(const bigint & bi);
	bigfloat(const bigint & n, const bigint & d);
	bigfloat(const bigrational & br);
	bigfloat(const bigfloat & bf);
	~bigfloat();



#ifndef HEADBANGER

	//
	// accessors
	//

	const bigint& mantissa() const;
	long exponent() const;
	long precision () const;



	//
	// assigners
	//

	void assign_zero();
	void assign_exact_zero();
	void assign_one();
	void assign(int i);
	void assign(long l);
	void assign(unsigned long ul);
	void assign(double d);
	void assign(const xdouble & d);
	void assign(const bigint & bi);
	void assign(const bigint & n, const bigint & d);
	void assign(const bigfloat & bf);

	bigfloat & operator = (int i);
	bigfloat & operator = (long l);
	bigfloat & operator = (unsigned long ul);
	bigfloat & operator = (double d);
	bigfloat & operator = (const xdouble & d);
	bigfloat & operator = (const bigint & bi);
	bigfloat & operator = (const bigrational & br);
	bigfloat & operator = (const bigfloat & bf);



	//
	// comparators
	//

	int abs_compare(const bigfloat & b) const;
	int compare(const bigfloat & b) const;

	bool is_zero() const;
	bool is_exact_zero() const;
	bool is_approx_zero() const;
	bool is_gt_zero() const;
	bool is_lt_zero() const;
	bool is_ge_zero() const;
	bool is_le_zero() const;
	bool is_one() const;

	int sign() const;



#endif	// HEADERBANGER

	//
	// predicates
	//

	bool is_char() const;
	bool is_uchar() const;
	bool is_short() const;
	bool is_ushort() const;
	bool is_int() const;
	bool is_uint() const;
	bool is_long() const;
	bool is_ulong() const;
	bool is_double() const;



	//
	// properties
	//

	lidia_size_t length() const;
	lidia_size_t bit_length() const;



	//
	// converters
	//

	void bigintify(bigint & bi) const;
	bool doublify(double & d) const;
	bool xdoublify(xdouble & xd) const;
	bool intify(int & i) const;
	bool longify(long & i) const;



	//
	// modifiers
	//

	void absolute_value();

#ifndef HEADBANGER

	void negate();
	void invert();

	void inc();
	void dec();

	void multiply_by_2();
	void divide_by_2();

	void swap(bigfloat & x);

#endif	// HEADBANGER



	void randomize(const bigint & n, long f);



#ifndef HEADBANGER

	//
	// procedural arithmetic
	//

	friend void add(bigfloat & c, const bigfloat & a, const bigfloat & b);
	friend void subtract(bigfloat & c, const bigfloat & a, const bigfloat & b);

	friend void multiply(bigfloat & c, const bigfloat & a, const bigfloat & b);
	friend void multiply(bigfloat & c, const bigfloat & a, long i);
	friend void multiply(bigfloat & c, long i, const bigfloat & b);

	friend void square(bigfloat & c, const bigfloat & a);

	friend void divide(bigfloat & c, const bigfloat & a, const bigfloat & b);
	friend void divide(bigfloat & c, const bigfloat & a, long i);
	friend void divide(bigfloat & c, long i, const bigfloat & a);
	friend void divide(bigfloat & c, long i, long j);

	friend void shift_left(bigfloat & c, const bigfloat & a, long i);
	friend void shift_right(bigfloat & c, const bigfloat & a, long i);

	friend void truncate(bigint & b, const bigfloat & a);
	friend void truncate(bigfloat & b, const bigfloat & a);



	//
	// Transcedental functions(procedural)
	//

	friend void sqrt(bigfloat & c, const bigfloat & a);

	friend void exp(bigfloat & c, const bigfloat & a);
	friend void log(bigfloat & c, const bigfloat & a);

	friend void power(bigfloat & c, const bigfloat & a, long y);
	friend void power(bigfloat & c, const bigfloat & a, const bigfloat & b);

	friend void sin(bigfloat & c, const bigfloat & a);
	friend void cos(bigfloat & c, const bigfloat & a);
	friend void tan(bigfloat & c, const bigfloat & a);
	friend void cot(bigfloat & c, const bigfloat & a);
	friend void asin(bigfloat & c, const bigfloat & a);
	friend void acos(bigfloat & c, const bigfloat & a);
	friend void atan(bigfloat & c, const bigfloat & a);
	friend void atan2(bigfloat & c, const bigfloat & a, const bigfloat & b);
	friend void acot(bigfloat & c, const bigfloat & a);

	friend void sinh(bigfloat & c, const bigfloat & a);
	friend void cosh(bigfloat & c, const bigfloat & a);
	friend void tanh(bigfloat & c, const bigfloat & a);
	friend void coth(bigfloat & c, const bigfloat & a);
	friend void asinh(bigfloat & c, const bigfloat & a);
	friend void acosh(bigfloat & c, const bigfloat & a);
	friend void atanh(bigfloat & c, const bigfloat & a);
	friend void acoth(bigfloat & c, const bigfloat & a);

	friend void besselj(bigfloat & c, int n, const bigfloat & a);
	
#endif	// HEADBANGER



	//
	// constants
	//

	friend void constant_E(bigfloat & x);
	friend void constant_Pi(bigfloat & x);
	friend void constant_Euler(bigfloat & x);
	friend void constant_Catalan(bigfloat & x);

	friend bigfloat E();
	friend bigfloat Pi();
	friend bigfloat Euler();
	friend bigfloat Catalan();


	//
	// input/output
	//

	friend std::istream & operator >> (std::istream & in, bigfloat & a);
	friend std::ostream & operator << (std::ostream & out, const bigfloat & a);

	friend int string_to_bigfloat(char *s, bigfloat & a);
	friend int bigfloat_to_string(const bigfloat & a, char *s);

#ifdef C_STDIO
	//
	// Using fscanf/fprintf
	//

	int scan_from_file(FILE * fp);
	int print_to_file(FILE * fp);
#endif	// C_STDIO

	void print(char *s);



	//
	//  Arithmetic with error bounds, MM
	//

	friend void add      (bigfloat & c, const bigfloat & a,
			      const bigfloat & b,
			      long t);

	friend void subtract (bigfloat & c, const bigfloat & a,
			      const bigfloat & b,
			      long t);

	friend void multiply (bigfloat & c, const bigfloat & a,
			      const bigfloat & b,
			      long t);

	friend void square   (bigfloat & c, const bigfloat & a,
			      long t);

	friend void divide   (bigfloat & c, const bigfloat & a,
			      const bigfloat & b,
			      long t);

	friend void invert   (bigfloat & c, const bigfloat & a,
			      long t);

	friend void sqrt (bigfloat & y, const bigfloat & x,
			  long t);

	friend void log(bigfloat & y, const bigfloat & x,
			long t);

	friend void log_absolute(bigfloat & y, const bigfloat & x,
				 long t);

	friend void exp(bigfloat & y, const bigfloat & x,
			long t);


	void print_binary (std::ostream & out) const;
	void cut_absolute (long k);



	//
	// class modifiers and accessors
	//

	// Precision and rounding mode setting
	//static void precision(long l); // obsolete, use set_precision
	static void set_precision(long l);
	static long get_precision();
	//static void mode(int i); // obsolete, use set_mode
	static void set_mode(int i);
	static void set_mode_no_check(int mode);
	static int get_mode();



	//
	// misc. functions
	//

	static bigint characteristic();

};


void add(bigfloat & c, const bigfloat & a, const bigfloat & b);
void subtract(bigfloat & c, const bigfloat & a, const bigfloat & b);

void multiply(bigfloat & c, const bigfloat & a, const bigfloat & b);
void multiply(bigfloat & c, const bigfloat & a, long i);
void multiply(bigfloat & c, long i, const bigfloat & b);

void square(bigfloat & c, const bigfloat & a);

void divide(bigfloat & c, const bigfloat & a, const bigfloat & b);
void divide(bigfloat & c, const bigfloat & a, long i);
void divide(bigfloat & c, long i, const bigfloat & a);
void divide(bigfloat & c, long i, long j);

void shift_left(bigfloat & c, const bigfloat & a, long i);
void shift_right(bigfloat & c, const bigfloat & a, long i);

void truncate(bigint & b, const bigfloat & a);
void truncate(bigfloat & b, const bigfloat & a);



//
// Transcedental functions(procedural)
//

void sqrt(bigfloat & c, const bigfloat & a);

void exp(bigfloat & c, const bigfloat & a);
void log(bigfloat & c, const bigfloat & a);

void power(bigfloat & c, const bigfloat & a, long y);
void power(bigfloat & c, const bigfloat & a, const bigfloat & b);

void sin(bigfloat & c, const bigfloat & a);
void cos(bigfloat & c, const bigfloat & a);
void tan(bigfloat & c, const bigfloat & a);
void cot(bigfloat & c, const bigfloat & a);
void asin(bigfloat & c, const bigfloat & a);
void acos(bigfloat & c, const bigfloat & a);
void atan(bigfloat & c, const bigfloat & a);
void atan2(bigfloat & c, const bigfloat & a, const bigfloat & b);
void acot(bigfloat & c, const bigfloat & a);

void sinh(bigfloat & c, const bigfloat & a);
void cosh(bigfloat & c, const bigfloat & a);
void tanh(bigfloat & c, const bigfloat & a);
void coth(bigfloat & c, const bigfloat & a);
void asinh(bigfloat & c, const bigfloat & a);
void acosh(bigfloat & c, const bigfloat & a);
void atanh(bigfloat & c, const bigfloat & a);
void acoth(bigfloat & c, const bigfloat & a);

void besselj(bigfloat & c, int n, const bigfloat & a);

//
// constants
//

void constant_E(bigfloat & x);
void constant_Pi(bigfloat & x);
void constant_Euler(bigfloat & x);
void constant_Catalan(bigfloat & x);

bigfloat E();
bigfloat Pi();
bigfloat Euler();
bigfloat Catalan();


//
// input/output
//

std::istream & operator >> (std::istream & in, bigfloat & a);
std::ostream & operator << (std::ostream & out, const bigfloat & a);

int string_to_bigfloat(char *s, bigfloat & a);
int bigfloat_to_string(const bigfloat & a, char *s);

#ifdef C_STDIO
//
// Using fscanf/fprintf
//

int scan_from_file(FILE * fp);
int print_to_file(FILE * fp);
#endif	// C_STDIO

void print(char *s);


//
//  Arithmetic with error bounds, MM
//

void add      (bigfloat & c, const bigfloat & a,
	       const bigfloat & b,
	       long t);

void subtract (bigfloat & c, const bigfloat & a,
	       const bigfloat & b,
	       long t);

void multiply (bigfloat & c, const bigfloat & a,
	       const bigfloat & b,
	       long t);

void square   (bigfloat & c, const bigfloat & a,
	       long t);

void divide   (bigfloat & c, const bigfloat & a,
	       const bigfloat & b,
	       long t);

void invert   (bigfloat & c, const bigfloat & a,
	       long t);

void sqrt (bigfloat & y, const bigfloat & x,
	   long t);

void log(bigfloat & y, const bigfloat & x,
	 long t);

void log_absolute(bigfloat & y, const bigfloat & x,
		  long t);

void exp(bigfloat & y, const bigfloat & x,
	 long t);


#ifndef HEADBANGER
int check_overflow(long &z, long x, long y);
#endif



//
// c'tors and d'tor
//

inline
bigfloat::bigfloat ()
{
	this->assign_zero();
}



inline
bigfloat::bigfloat (int i)
{
	this->assign(i);
}



inline
bigfloat::bigfloat (long l)
{
	this->assign(l);
}



inline
bigfloat::bigfloat (unsigned long ul)
{
	this->assign(ul);
}



inline
bigfloat::bigfloat (double d)
{
	this->assign(d);
}



inline
bigfloat::bigfloat (const xdouble & xd)
{
	this->assign(xd);
}



inline
bigfloat::bigfloat (const bigint & bi)
{
	this->assign(bi);
}



inline
bigfloat::bigfloat (const bigrational & br)
{
	this->assign(br.numerator(), br.denominator());
}



inline
bigfloat::bigfloat (const bigfloat & bf)
{
	this->assign(bf);
}



inline
bigfloat::bigfloat (const bigint & n, const bigint & d)
{
	this->assign(n, d);
}



inline
bigfloat::~bigfloat()
{
	// nothing to do
}



//
// accessors
//

inline long
bigfloat::exponent () const
{
	return e;
}



inline const bigint &
bigfloat::mantissa () const
{
	return m;
}



inline long
bigfloat::precision () const
{
	return prec;
}



//
// assigners
//

inline bigfloat &
bigfloat::operator = (int i)
{
	this->assign(i);
	return *this;
}



inline bigfloat &
bigfloat::operator = (long l)
{
	this->assign(l);
	return *this;
}



inline bigfloat &
bigfloat::operator = (unsigned long ul)
{
	this->assign(ul);
	return *this;
}



inline bigfloat &
bigfloat::operator = (double d)
{
	this->assign(d);
	return *this;
}



inline bigfloat &
bigfloat::operator = (const xdouble & xd)
{
	this->assign(xd);
	return *this;
}



inline bigfloat &
bigfloat::operator = (const bigint & bi)
{
	this->assign(bi);
	return *this;
}



inline bigfloat &
bigfloat::operator = (const bigrational & br)
{
	this->assign(br);
	return *this;
}



inline bigfloat &
bigfloat::operator = (const bigfloat & bf)
{
	this->assign(bf);
	return *this;
}



//
// comparators
//

inline bool
bigfloat::is_zero () const
{
	return (m.is_zero() || this->is_pm_zero());
}



inline bool
bigfloat::is_exact_zero () const
{
	return (m.is_zero() && this->is_exact());
}



inline bool
bigfloat::is_approx_zero () const
{
	return (m.is_zero() || ((e + prec) < (-bigfloat::binary_precision + rounding_bits)));
}



inline bool
bigfloat::is_gt_zero () const
{
	return m.is_gt_zero();
}



inline bool
bigfloat::is_lt_zero () const
{
	return m.is_lt_zero();
}



inline bool
bigfloat::is_ge_zero () const
{
	return m.is_ge_zero();
}



inline bool
bigfloat::is_le_zero () const
{
	return m.is_le_zero();
}



inline bool
bigfloat::is_one () const
{
	bigfloat tmp(1UL);

	subtract(tmp, tmp, *this);
	return tmp.is_zero();
}



inline bool
operator == (const bigfloat & x, const bigfloat & y)
{
	return (x.compare(y) == 0);
}



inline bool
operator != (const bigfloat & x, const bigfloat & y)
{
	return (x.compare(y) != 0);
}



inline bool
operator < (const bigfloat & x, const bigfloat & y)
{
	return (x.compare(y) < 0);
}



inline bool
operator <= (const bigfloat & x, const bigfloat & y)
{
	return (x.compare(y) <= 0);
}



inline bool
operator > (const bigfloat & x, const bigfloat & y)
{
	return (x.compare(y) > 0);
}



inline bool
operator >= (const bigfloat & x, const bigfloat & y)
{
	return (x.compare(y) >= 0);
}



inline int
compare (const bigfloat & x, const bigfloat & y)
{
	return x.compare(y);
}



inline int
bigfloat::sign () const
{
	return m.sign();
}



//
// properties
//

inline lidia_size_t
bigfloat::length () const
{
	return m.length();
}



inline lidia_size_t
bigfloat::bit_length () const
{
	return prec;
}



//
// predicates
//

inline bool
is_char(const bigfloat & a)
{
	return a.is_char();
}



inline bool
is_uchar(const bigfloat & a)
{
	return a.is_uchar();
}



inline bool
is_short(const bigfloat & a)
{
	return a.is_short();
}



inline bool
is_ushort(const bigfloat & a)
{
	return a.is_ushort();
}



inline bool
is_int(const bigfloat & a)
{
	return a.is_int();
}



inline bool
is_uint(const bigfloat & a)
{
	return a.is_uint();
}



inline bool
is_long(const bigfloat & a)
{
	return a.is_long();
}



inline bool
is_ulong(const bigfloat & a)
{
	return a.is_ulong();
}



inline bool
is_double(const bigfloat & a)
{
	return a.is_double();
}



//
// modifiers
//

inline void
swap (bigfloat & x, bigfloat & y)
{
	x.swap(y);
}



inline void
bigfloat::absolute_value ()
{
	m.absolute_value();
}



inline void
bigfloat::negate ()
{
	m.negate();
}



inline void
bigfloat::inc ()
{
	add(*this, *this, 1UL);
}



inline void
bigfloat::dec ()
{
	subtract(*this, *this, 1UL);
}



inline void
bigfloat::multiply_by_2 ()
{
	if (!m.is_zero())
		e++;
}



inline void
bigfloat::divide_by_2 ()
{
	if (!m.is_zero())
		e--;
}



//
// procedural arithmetic
//

inline void
negate (bigfloat & x, const bigfloat & y)
{
	x.assign(y);
	x.negate();
}



void add(bigfloat & c, const bigfloat & a, const bigfloat & b);
void subtract(bigfloat & c, const bigfloat & a, const bigfloat & b);

void multiply(bigfloat & c, const bigfloat & a, const bigfloat & b);
void multiply(bigfloat & c, const bigfloat & a, long i);
void multiply(bigfloat & c, long i, const bigfloat & b);

void square(bigfloat & c, const bigfloat & a);

void divide(bigfloat & c, const bigfloat & a, const bigfloat & b);
void divide(bigfloat & c, const bigfloat & a, long i);
void divide(bigfloat & c, long i, const bigfloat & a);
void divide(bigfloat & c, long i, long j);



inline void
shift_left (bigfloat & res, const bigfloat & x, long i)
{
	if (&res != &x)
		res.assign(x);
	res.e = x.e + i;
}



inline void
shift_right (bigfloat & res, const bigfloat & x, long i)
{
	if (&res != &x)
		res.assign(x);
	res.e = x.e - i;
}



inline void
invert (bigfloat & x, const bigfloat & y)
{
	x.assign(y);
	x.invert();
}



inline void
inc (bigfloat & x)
{
	x.inc();
}



inline void
dec (bigfloat & x)
{
	x.dec();
}



//
// operators
//

inline bigfloat
operator - (const bigfloat & a)
{
	bigfloat c(a);

	c.negate();
	return c;
}



inline bigfloat
operator + (const bigfloat & a, const bigfloat & b)
{
	bigfloat c;

	add(c, a, b);
	return c;
}



inline bigfloat
operator - (const bigfloat & a, const bigfloat & b)
{
	bigfloat c;

	subtract(c, a, b);
	return c;
}



inline bigfloat
operator * (const bigfloat & a, const bigfloat & b)
{
	bigfloat c;

	multiply(c, a, b);
	return c;
}



inline bigfloat
operator / (const bigfloat & a, const bigfloat & b)
{
	bigfloat c;

	divide(c, a, b);
	return c;
}



inline bigfloat &
operator += (bigfloat & a, const bigfloat & b)
{
	add(a, a, b);
	return a;
}



inline bigfloat &
operator -= (bigfloat & a, const bigfloat & b)
{
	subtract(a, a, b);
	return a;
}



inline bigfloat &
operator *= (bigfloat & a, const bigfloat & b)
{
	multiply(a, a, b);
	return a;
}



inline bigfloat &
operator /= (bigfloat & a, const bigfloat & b)
{
	divide(a, a, b);
	return a;
}



inline bigfloat &
operator ++ (bigfloat & a)
{
	a.inc();
	return a;
}



inline bigfloat
operator ++ (bigfloat & a, int)
{
	bigfloat c(a);

	a.inc();
	return c;
}



inline bigfloat &
operator -- (bigfloat & a)
{
	a.dec();
	return a;
}



inline bigfloat
operator -- (bigfloat & a, int)
{
	bigfloat c(a);

	a.dec();
	return c;
}



inline bigfloat
operator << (const bigfloat & a, unsigned long ui)
{
	bigfloat c;

	shift_left(c, a, ui);
	return c;
}



inline bigfloat
operator >> (const bigfloat & a, unsigned long ui)
{
	bigfloat c;

	shift_right(c, a, ui);
	return c;
}



inline bigfloat &
operator <<= (bigfloat & a, unsigned long ui)
{
	shift_left(a, a, ui);
	return a;
}



inline bigfloat &
operator >>= (bigfloat & a, unsigned long ui)
{
	shift_right(a, a, ui);
	return a;
}



inline bool
operator ! (const bigfloat & a)
{
	return a.is_zero();
}



//
// functions
//

inline int
sign (const bigfloat & x)
{
	return x.sign();
}



inline bigfloat
abs (const bigfloat & x)
{
	bigfloat y(x);

	y.absolute_value();
	return y;
}



inline bigfloat
inverse (const bigfloat & y)
{
	bigfloat x(y);

	x.invert();
	return x;
}



//
// functional accessors
//

inline long
exponent (const bigfloat & x)
{
	return x.exponent();
}



inline const bigint &
mantissa (const bigfloat & x)
{
	return x.mantissa();
}


inline long
precision (const bigfloat & x)
{
	return x.precision();
}



//
// rounding functions
//

void truncate(bigint & b, const bigfloat & a);
void round(bigint & b, const bigfloat & a);
void floor(bigint & b, const bigfloat & a);
void ceil(bigint & b, const bigfloat & a);

void truncate(bigfloat & b, const bigfloat & a);
void round(bigfloat & b, const bigfloat & a);
void floor(bigfloat & b, const bigfloat & a);
void ceil(bigfloat & b, const bigfloat & a);



inline bigfloat
ceil (const bigfloat & x)
{
	bigfloat y;

	ceil(y, x);
	return y;
}



inline bigfloat
floor (const bigfloat & x)
{
	bigfloat y;

	floor(y, x);
	return y;
}



inline bigfloat
truncate (const bigfloat & x)
{
	bigfloat y(x);

	truncate(y, x);
	return y;
}



inline bigfloat
round (const bigfloat & x)
{
	bigfloat y;

	round(y, x);
	return y;
}



//
// exp, log, power, and related
//

void sqrt(bigfloat & c, const bigfloat & a);

void exp(bigfloat & c, const bigfloat & a);
void log(bigfloat & c, const bigfloat & a);

void power(bigfloat & c, const bigfloat & a, long y);
void power(bigfloat & c, const bigfloat & a, const bigfloat & b);



inline bigfloat
sqrt (const bigfloat & x)
{
	bigfloat y;

	sqrt(y, x);
	return y;
}



inline bigfloat
square (const bigfloat & a)
{
	bigfloat c;

	square(c, a);
	return c;
}



inline bigfloat
exp (const bigfloat & x)
{
	bigfloat y;

	exp(y, x);
	return y;
}



inline bigfloat
log (const bigfloat & x)
{
	bigfloat y;

	log(y, x);
	return y;
}



inline bigfloat
power (const bigfloat & x, const bigfloat & y)
{
	bigfloat z;

	power(z, x, y);
	return z;
}



inline bigfloat
power (const bigfloat & x, long y)
{
	bigfloat z;

	power(z, x, y);
	return z;
}



//
// (Inverse) trigonometrical functions(procedural)
//

void sin(bigfloat & c, const bigfloat & a);
void cos(bigfloat & c, const bigfloat & a);
void tan(bigfloat & c, const bigfloat & a);
void cot(bigfloat & c, const bigfloat & a);
void asin(bigfloat & c, const bigfloat & a);
void acos(bigfloat & c, const bigfloat & a);
void atan(bigfloat & c, const bigfloat & a);
void atan2(bigfloat & c, const bigfloat & a, const bigfloat & b);
void acot(bigfloat & c, const bigfloat & a);



inline bigfloat
sin (const bigfloat & x)
{
	bigfloat y;
	sin(y, x);
	return y;
}



inline bigfloat
cos (const bigfloat & x)
{
	bigfloat y;

	cos(y, x);
	return y;
}



inline bigfloat
tan (const bigfloat & x)
{
	bigfloat y;

	tan(y, x);
	return y;
}



inline bigfloat
cot (const bigfloat & x)
{
	bigfloat y;

	cot(y, x);
	return y;
}



inline bigfloat
asin (const bigfloat & x)
{
	bigfloat y;

	asin(y, x);
	return y;
}



inline bigfloat
acos (const bigfloat & x)
{
	bigfloat y;

	acos(y, x);
	return y;
}



inline bigfloat
atan (const bigfloat & x)
{
	bigfloat y;

	atan(y, x);
	return y;
}



inline bigfloat
atan2 (const bigfloat & y, const bigfloat & x)
{
	bigfloat z;

	atan2(z, y, x);
	return z;
}



inline bigfloat
acot (const bigfloat & x)
{
	bigfloat y;

	acot(y, x);
	return y;
}



//
// (Inverse) Hyperbolic trigonometrical functions(procedural)
//

void sinh(bigfloat & c, const bigfloat & a);
void cosh(bigfloat & c, const bigfloat & a);
void tanh(bigfloat & c, const bigfloat & a);
void coth(bigfloat & c, const bigfloat & a);
void asinh(bigfloat & c, const bigfloat & a);
void acosh(bigfloat & c, const bigfloat & a);
void atanh(bigfloat & c, const bigfloat & a);
void acoth(bigfloat & c, const bigfloat & a);



inline bigfloat
sinh (const bigfloat & x)
{
	bigfloat y;

	sinh(y, x);
	return y;
}



inline bigfloat
cosh (const bigfloat & x)
{
	bigfloat y;

	cosh(y, x);
	return y;
}



inline bigfloat
tanh (const bigfloat & x)
{
	bigfloat y;

	tanh(y, x);
	return y;
}



inline bigfloat
coth (const bigfloat & x)
{
	bigfloat y;

	coth(y, x);
	return y;
}



inline bigfloat
asinh (const bigfloat & x)
{
	bigfloat y;

	asinh(y, x);
	return y;
}



inline bigfloat
acosh (const bigfloat & x)
{
	bigfloat y;

	acosh(y, x);
	return y;
}



inline bigfloat
acoth (const bigfloat & x)
{
	bigfloat y;

	acoth(y, x);
	return y;
}



inline bigfloat
atanh (const bigfloat & x)
{
	bigfloat y;

	atanh(y, x);
	return y;
}



//
// bessel function
//

void besselj(bigfloat & c, int n, const bigfloat & a);



inline bigfloat
besselj (int n, const bigfloat & x)
{
	bigfloat y;

	besselj(y, n, x);
	return y;
}



//
// constants
//

void constant_E(bigfloat & x);
void constant_Pi(bigfloat & x);
void constant_Euler(bigfloat & x);
void constant_Catalan(bigfloat & x);



inline bigfloat
Catalan ()
{
	bigfloat x;

	constant_Catalan(x);
	return x;
}



inline bigfloat
E ()
{
	bigfloat x;

	constant_E(x);
	return x;
}



inline bigfloat
Euler ()
{
	bigfloat x;

	constant_Euler(x);
	return x;
}



inline bigfloat
Pi ()
{
	bigfloat x;

	constant_Pi(x);
	return x;
}



//
// misc
//

inline void
bigfloat::set_mode_no_check(int mode)
{
	bigfloat::rounding_mode = mode;
}



inline int
bigfloat::get_mode()
{
	return bigfloat::rounding_mode;
}



inline long
bigfloat::get_precision ()
{
	return bigfloat::decimal_precision;
}



inline bigint
bigfloat::characteristic ()
{
	return 0UL;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_BIGFLOAT

#include	"LiDIA/specialization/bigfloat.special"



#endif	// LIDIA_BIGFLOAT_H_GUARD_
