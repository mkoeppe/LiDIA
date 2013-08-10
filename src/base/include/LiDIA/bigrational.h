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


#ifndef LIDIA_BIGRATIONAL_H_GUARD_
#define LIDIA_BIGRATIONAL_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigrational
{
private:

	//
	// the C++ type we use to represent a bigrational
	//

	bigint num, den;


	void normalize();



public:

	//
	// c'tors and d'tor
	//

	bigrational();
	bigrational(int n);
	bigrational(long n);
	bigrational(unsigned long n);
	bigrational(double d);
	bigrational(const bigint & n);
	bigrational(const bigint & n, const bigint & d);
	bigrational(const bigrational & a);
	~bigrational();



#ifndef HEADBANGER

	//
	// accessors
	//

	const bigint & numerator() const;
	const bigint & denominator() const;
	


	//
	// assigners
	//

	void assign_zero();
	void assign_one();

	void assign(int n);
	void assign(long n);
	void assign(unsigned long n);
	void assign(double d);
	void assign(const bigint & n);
	void assign(const bigint & n, const bigint & d);
	void assign(const bigrational & a);
	
	void assign_numerator(const bigint & n);
	void assign_denominator(const bigint & d);

	bigrational & operator = (int n);
	bigrational & operator = (long n);
	bigrational & operator = (unsigned long n);
	bigrational & operator = (double d);
	bigrational & operator = (const bigint & a);
	bigrational & operator = (const bigrational & a);



	//
	// comparators
	//

	int abs_compare(const bigrational & a) const;
	int abs_compare(const bigint & a) const;
	int abs_compare(unsigned long a) const;
	int compare(const bigrational & a) const;
	int compare(const bigint & a) const;
	int compare(const long a) const;
	int compare(const unsigned long a) const;

	int sign() const;

	bool is_positive() const;
	bool is_negative() const;
	bool is_zero() const;
	bool is_gt_zero() const;
	bool is_ge_zero() const;
	bool is_lt_zero() const;
	bool is_le_zero() const;
	bool is_one() const;



	//
	// converters
	//

	bool intify(int & i) const;
	bool longify(long & l) const;



	//
	// modifiers
	//

	void absolute_value();
	void negate();
	void invert();

	void inc();
	void dec();

	void multiply_by_denominator();

	void multiply_by_2();
	void divide_by_2();

	void swap(bigrational & a);



	//
	// procedural arithmetic
	//

	friend void add(bigrational & c, const bigrational & a, const bigrational & b);
	friend void add(bigrational & c, const bigrational & a, const bigint & b);
	friend void add(bigrational & c, const bigrational & a, long b);
	friend void add(bigrational & c, const bigrational & a, unsigned long b);
	friend void subtract(bigrational & c, const bigrational & a, const bigrational & b);
	friend void subtract(bigrational & c, const bigrational & a, const bigint & b);
	friend void subtract(bigrational & c, const bigrational & a, long b);
	friend void subtract(bigrational & c, const bigrational & a, unsigned long b);

	friend void multiply(bigrational & c, const bigrational & a, const bigrational & b);
	friend void multiply(bigrational & c, const bigrational & a, const bigint & b);
	friend void multiply(bigrational & c, const bigrational & a, long b);
	friend void multiply(bigrational & c, const bigrational & a, unsigned long b);
	friend void square(bigrational & a, const bigrational & b);
	friend void square(bigrational & a, const bigint & b);

	friend void divide(bigrational & c, const bigrational & a, const bigrational & b);
	friend void divide(bigrational & c, const bigrational & a, const bigint & b);
	friend void divide(bigrational & c, const bigrational & a, long b);
	friend void divide(bigrational & c, const bigrational & a, unsigned long b);

	friend void shift_left(bigrational & c, const bigrational & a, long ui);
	friend void shift_right(bigrational & c, const bigrational & a, long ui);


#endif

	//
	// input / output
	//

	friend std::istream & operator >> (std::istream & in, bigrational & a);
	friend std::ostream & operator << (std::ostream & out, const bigrational & a);

	friend int string_to_bigrational(char *s, char *t, bigrational & a);
	friend int bigrational_to_string(const bigrational & a, char *s, char *t);

#ifdef C_STDIO
	//
	// using fread/fwrite
	//

	void read_from_file(FILE * fp);
	void write_to_file(FILE * fp);

	//
	// using fscanf/fprintf
	//

	void scan_from_file(FILE * fp);
	void print_to_file(FILE * fp);
#endif	// C_STDIO

	// for use in elliptic curve code
	static bigint characteristic();

};

void add(bigrational & c, const bigrational & a, const bigrational & b);
void add(bigrational & c, const bigrational & a, const bigint & b);
void add(bigrational & c, const bigrational & a, long b);
void add(bigrational & c, const bigrational & a, unsigned long b);
void subtract(bigrational & c, const bigrational & a, const bigrational & b);
void subtract(bigrational & c, const bigrational & a, const bigint & b);
void subtract(bigrational & c, const bigrational & a, long b);
void subtract(bigrational & c, const bigrational & a, unsigned long b);

void multiply(bigrational & c, const bigrational & a, const bigrational & b);
void multiply(bigrational & c, const bigrational & a, const bigint & b);
void multiply(bigrational & c, const bigrational & a, long b);
void multiply(bigrational & c, const bigrational & a, unsigned long b);
void square(bigrational & a, const bigrational & b);
void square(bigrational & a, const bigint & b);

void divide(bigrational & c, const bigrational & a, const bigrational & b);
void divide(bigrational & c, const bigrational & a, const bigint & b);
void divide(bigrational & c, const bigrational & a, long b);
void divide(bigrational & c, const bigrational & a, unsigned long b);

void shift_left(bigrational & c, const bigrational & a, long ui);
void shift_right(bigrational & c, const bigrational & a, long ui);

//
// input / output
//

std::istream & operator >> (std::istream & in, bigrational & a);
std::ostream & operator << (std::ostream & out, const bigrational & a);

int string_to_bigrational(char *s, char *t, bigrational & a);
int bigrational_to_string(const bigrational & a, char *s, char *t);


//
// c'tors and d'tor
//

inline
bigrational::bigrational ()
	: num(0UL),
	  den(1UL)
{
	// nothing to do
}



inline
bigrational::bigrational (int n)
	: num(n),
	  den(1UL)
{
	// nothing to do
}



inline
bigrational::bigrational (long n)
	: num(n),
	  den(1UL)
{
	// nothing to do
}



inline
bigrational::bigrational (unsigned long n)
	: num(n),
	  den(1UL)
{
	// nothing to do
}



inline
bigrational::bigrational (double d)
{
	assign(d);
}



inline
bigrational::bigrational (const bigint & n)
	: num(n),
	  den(1UL)
{
	// nothing to do
}



inline
bigrational::bigrational (const bigint & n, const bigint & d)
	: num(n),
	  den(d)
{
	if (d.is_zero()) {
		lidia_error_handler("bigrational", "constructor(n, d)::division by zero.");
		// actually call assign_zero(), but this is declared inline later
		num.assign_zero();
		den.assign_zero();
	}
	else {
		normalize();
	}
}



inline
bigrational::bigrational (const bigrational & a)
	: num(a.num),
	  den(a.den)
{
	// nothing to do
}



inline
bigrational::~bigrational()
{
}



//
// accessors
//

inline const bigint &
bigrational::numerator () const
{
	return num;
}



inline const bigint &
bigrational::denominator () const
{
	return den;
}



inline const bigint &
numerator (const bigrational & a)
{
	return a.numerator();
}



inline const bigint &
denominator (const bigrational & a)
{
	return a.denominator();
}



//
// assigners
//

inline void
bigrational::assign_zero ()
{
	num.assign_zero();
	den.assign_one();
}



inline void
bigrational::assign_one ()
{
	num.assign_one();
	den.assign_one();
}



inline void
bigrational::assign(int n)
{
	num.assign(n);
	den.assign_one();
}



inline void
bigrational::assign(long n)
{
	num.assign(n);
	den.assign_one();
}



inline void
bigrational::assign(unsigned long n)
{
	num.assign(n);
	den.assign_one();
}



inline void
bigrational::assign(const bigint & n)
{
	num.assign(n);
	den.assign_one();
}



inline void
bigrational::assign (const bigint & n, const bigint & d)
{
	if (d.is_zero()) {
		lidia_error_handler("bigrational", "assign(n, d)::division by zero.");
		assign_zero();
		return;
	}
	num.assign(n);
	den.assign(d);
	normalize();
}



inline void
bigrational::assign (const bigrational & a)
{
	if (&a != this) {
		num.assign(a.num);
		den.assign(a.den);
	}
}



inline void
bigrational::assign_numerator (const bigint & n)
{
	if (&n != &num) {
		num.assign(n);
	}
}



inline void
bigrational::assign_denominator (const bigint & d)
{
	if (d.is_zero()) {
		lidia_error_handler("bigrational", "assign(n, d)::division by zero.");
		return;
	}
	den.assign(d);
	if (den.is_negative()) {
		num.negate();
		den.negate();
	}
}



inline bigrational &
bigrational::operator = (int n)
{
	assign(n);
	return *this;
}



inline bigrational &
bigrational::operator = (long n)
{
	assign(n);
	return *this;
}



inline bigrational &
bigrational::operator = (unsigned long n)
{
	assign(n);
	return *this;
}



inline bigrational &
bigrational::operator = (double d)
{
	assign(d);
	return *this;
}



inline bigrational &
bigrational::operator = (const bigint & n)
{
	assign(n);
	return *this;
}



inline bigrational &
bigrational::operator = (const bigrational & a)
{
	assign(a);
	return *this;
}



//
// comparators
//

inline int
bigrational::abs_compare (const bigrational & a) const
{
	if (&a == this) {
		return 0;
	}

	bigint r1, r2;

	multiply(r1, num, a.den);
	multiply(r2, den, a.num);
	return r1.abs_compare(r2);
}



inline int
bigrational::abs_compare (const bigint & a) const
{
	if (den.is_one()) {
		return num.compare(a);
	}

	bigint r;

	multiply(r, den, a);
	return num.abs_compare(r);
}



inline int
bigrational::abs_compare (unsigned long a) const
{
	if (den.is_one()) {
		return num.compare(a);
	}

	bigint r;

	multiply(r, den, a);
	return num.abs_compare(r);
}



inline int
abs_compare (const bigrational & a, const bigrational & b)
{
	return a.abs_compare(b);
}



inline int
abs_compare (const bigrational & a, const bigint & b)
{
	return a.abs_compare(b);
}



inline int
abs_compare (const bigrational & a, const unsigned long b)
{
	return a.abs_compare(b);
}



inline int
bigrational::compare (const bigrational & a) const
{
	if (&a == this) {
		return 0;
	}

	bigint r1, r2;

	multiply(r1, num, a.den);
	multiply(r2, den, a.num);
	return r1.compare(r2);
}



inline int
bigrational::compare (const bigint & a) const
{
	if (den.is_one()) {
		return num.compare(a);
	}

	bigint r;

	multiply(r, den, a);
	return num.compare(r);
}



inline int
bigrational::compare (unsigned long a) const
{
	if (den.is_one()) {
		return num.compare(a);
	}

	bigint r;

	multiply(r, den, a);
	return num.compare(r);
}



inline int
bigrational::compare (long a) const
{
	if (den.is_one()) {
		return num.compare(a);
	}

	bigint r;

	multiply(r, den, a);
	return num.compare(r);
}



inline int
compare (const bigrational & a, const bigrational & b)
{
	return a.compare(b);
}



inline int
compare (const bigrational & a, const bigint & b)
{
	return a.compare(b);
}



inline int
compare (const bigrational & a, unsigned long b)
{
	return a.compare(b);
}



inline int
compare (const bigrational & a, long b)
{
	return a.compare(b);
}



inline bool
operator == (const bigrational & a, const bigrational & b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigrational & a, const bigint & b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigrational & a, unsigned long b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigrational & a, long b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigrational & a, int b)
{
	return (a.compare(static_cast<long>(b)) == 0);
}



inline bool
operator == (const bigint & a, const bigrational & b)
{
	return (b.compare(a) == 0);
}



inline bool
operator == (const unsigned long a, const bigrational & b)
{
	return (b.compare(a) == 0);
}



inline bool
operator == (const long a, const bigrational & b)
{
	return (b.compare(a) == 0);
}



inline bool
operator == (const int a, const bigrational & b)
{
	return (b.compare(static_cast<long>(a)) == 0);
}



inline bool
operator != (const bigrational & a, const bigrational & b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigrational & a, const bigint & b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigrational & a, unsigned long b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigrational & a, long b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigrational & a, int b)
{
	return (a.compare(static_cast<long>(b)) != 0);
}



inline bool
operator != (const bigint & a, const bigrational & b)
{
	return (b.compare(a) != 0);
}



inline bool
operator != (long a, const bigrational & b)
{
	return (b.compare(a) != 0);
}



inline bool
operator != (unsigned long a, const bigrational & b)
{
	return (b.compare(a) != 0);
}



inline bool
operator != (int a, const bigrational & b)
{
	return (b.compare(static_cast<long>(a)) != 0);
}



inline bool
operator > (const bigrational & a, const bigrational & b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigrational & a, const bigint & b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigrational & a, long b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigrational & a, unsigned long b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigrational & a, int b)
{
	return (a.compare(static_cast<long>(b)) > 0);
}



inline bool
operator > (const bigint & a, const bigrational & b)
{
	return (b.compare(a) < 0);
}



inline bool
operator > (long a, const bigrational & b)
{
	return (b.compare(a) < 0);
}



inline bool
operator > (unsigned long a, const bigrational & b)
{
	return (b.compare(a) < 0);
}



inline bool
operator > (int a, const bigrational & b)
{
	return (b.compare(static_cast<long>(a)) < 0);
}



inline bool
operator >= (const bigrational & a, const bigrational & b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigrational & a, const bigint & b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigrational & a, long b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigrational & a, unsigned long b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigrational & a, int b)
{
	return (a.compare(static_cast<long>(b)) >= 0);
}



inline bool
operator >= (const bigint & a, const bigrational & b)
{
	return (b.compare(a) <= 0);
}



inline bool
operator >= (long a, const bigrational & b)
{
	return (b.compare(a) <= 0);
}



inline bool
operator >= (unsigned long a, const bigrational & b)
{
	return (b.compare(a) <= 0);
}



inline bool
operator >= (int a, const bigrational & b)
{
	return (b.compare(static_cast<long>(a)) <= 0);
}



inline bool
operator < (const bigrational & a, const bigrational & b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigrational & a, const bigint & b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigrational & a, long b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigrational & a, unsigned long b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigrational & a, int b)
{
	return (a.compare(static_cast<long>(b)) < 0);
}



inline bool
operator < (const bigint & a, const bigrational & b)
{
	return (b.compare(a) > 0);
}



inline bool
operator < (long a, const bigrational & b)
{
	return (b.compare(a) > 0);
}



inline bool
operator < (unsigned long a, const bigrational & b)
{
	return (b.compare(a) > 0);
}



inline bool
operator < (int a, const bigrational & b)
{
	return (b.compare(static_cast<long>(a)) > 0);
}



inline bool
operator <= (const bigrational & a, const bigrational & b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigrational & a, const bigint & b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigrational & a, long b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigrational & a, unsigned long b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigrational & a, int b)
{
	return (a.compare(static_cast<long>(b)) <= 0);
}



inline bool
operator <= (const bigint & a, const bigrational & b)
{
	return (b.compare(a) >= 0);
}



inline bool
operator <= (long a, const bigrational & b)
{
	return (b.compare(a) >= 0);
}



inline bool
operator <= (unsigned long a, const bigrational & b)
{
	return (b.compare(a) >= 0);
}



inline bool
operator <= (int a, const bigrational & b)
{
	return (b.compare(static_cast<long>(a)) >= 0);
}



inline int
bigrational::sign () const
{
	return num.sign();
}



inline bool
bigrational::is_positive () const
{
	return (num.sign() > 0);
}



inline bool
bigrational::is_negative () const
{
	return (num.sign() < 0);
}



inline bool
bigrational::is_zero () const
{
	return (num.sign() == 0);
}



inline bool
bigrational::is_gt_zero () const
{
	return (num.sign() > 0);
}



inline bool
bigrational::is_ge_zero () const
{
	return (num.sign() >= 0);
}



inline bool
bigrational::is_lt_zero () const
{
	return (num.sign() < 0);
}



inline bool
bigrational::is_le_zero () const
{
	return (num.sign() <= 0);
}



inline bool
bigrational::is_one () const
{
	return (num.is_one() && den.is_one());
}



//
// converters
//

inline bool
bigrational::intify (int & i) const
{
	bigint I;

	divide(I, num, den);
	return I.intify(i);
}



inline bool
bigrational::longify (long & i) const
{
	bigint I;

	divide(I, num, den);
	return I.longify(i);
}



double dbl (const bigrational & a);



//
// predicates
//

inline bool
is_bigint (const bigrational & a)
{
	return a.denominator().is_one();
}



//
// modifiers
//

inline void
bigrational::absolute_value ()
{
	num.absolute_value();
}



inline void
bigrational::negate ()
{
	num.negate();
}



inline void
bigrational::multiply_by_denominator ()
{
	den.assign_one();
}



inline void
bigrational::inc ()
{
	add(num, num, den);
}



inline void
bigrational::dec ()
{
	subtract(num, num, den);
}



inline void
bigrational::multiply_by_2 ()
{
	if (den.is_even()) {
		den.divide_by_2();
	}
	else {
		num.multiply_by_2();
	}
}



inline void
bigrational::divide_by_2 ()
{
	if (num.is_even()) {
		num.divide_by_2();
	}
	else {
		den.multiply_by_2();
	}
}



inline void
bigrational::swap (bigrational & a)
{
	num.swap(a.num);
	den.swap(a.den);
}



inline void
swap (bigrational & a, bigrational & b)
{
	a.swap(b);
}



//
// procedural arithmetic
//

void add(bigrational & c, const bigrational & a, const bigrational & b);
void subtract(bigrational & c, const bigrational & a, const bigrational & b);
void multiply(bigrational & c, const bigrational & a, const bigrational & b);
void divide(bigrational & c, const bigrational & a, const bigrational & b);

void shift_left(bigrational & c, const bigrational & a, long ui);
void shift_right(bigrational & c, const bigrational & a, long ui);



inline void
negate (bigrational & a, const bigrational & b)
{
	a.assign(b);
	a.negate();
}



inline void
add (bigrational & c, const bigrational & a, const bigint & b)
{
	bigint tmp;

	multiply(tmp, a.den, b);
	add(c.num, a.num, tmp);
	c.den.assign(a.den);
}



inline void
add (bigrational & c, const bigrational & a, long b)
{
	bigint tmp;

	multiply(tmp, a.den, b);
	add(c.num, a.num, tmp);
	c.den.assign(a.den);
}



inline void
add (bigrational & c, const bigrational & a, unsigned long b)
{
	bigint tmp;

	multiply(tmp, a.den, b);
	add(c.num, a.num, tmp);
	c.den.assign(a.den);
}



inline void
add (bigrational & c, const bigrational & a, int b)
{
	add(c, a, static_cast<long>(b));
}



inline void
subtract (bigrational & c, const bigrational & a, const bigint & b)
{
	bigint tmp;

	multiply(tmp, a.den, b);
	subtract(c.num, a.num, tmp);
	c.den.assign(a.den);
}



inline void
subtract (bigrational & c, const bigrational & a, long b)
{
	bigint tmp;

	multiply(tmp, a.den, b);
	subtract(c.num, a.num, tmp);
	c.den.assign(a.den);
}



inline void
subtract (bigrational & c, const bigrational & a, unsigned long b)
{
	bigint tmp;

	multiply(tmp, a.den, b);
	subtract(c.num, a.num, tmp);
	c.den.assign(a.den);
}



inline void
subtract (bigrational & c, const bigrational & a, int b)
{
	subtract(c, a, static_cast<long>(b));
}



inline void
multiply (bigrational & c, const bigrational & a, const bigint & b)
{
	multiply(c.num, a.num, b);
	c.den.assign(a.den);
	c.normalize();
}



inline void
multiply (bigrational & c, const bigrational & a, long b)
{
	multiply(c.num, a.num, b);
	c.den.assign(a.den);
	c.normalize();
}



inline void
multiply (bigrational & c, const bigrational & a, unsigned long b)
{
	multiply(c.num, a.num, b);
	c.den.assign(a.den);
	c.normalize();
}



inline void
multiply (bigrational & c, const bigrational & a, int b)
{
	multiply(c, a, static_cast<long>(b));
}



inline void
square (bigrational & c, const bigrational & a)
{
	square(c.num, a.num);
	square(c.den, a.den);
}



inline void
square (bigrational & c, const bigint & a)
{
	square(c.num, a);
	c.den.assign_one();
}



inline void
invert (bigrational & c, const bigrational & a)
{
	c.assign(a);
	c.invert();
}



inline void
invert (bigrational & c, const bigint & a)
{
	c.assign(1UL, a);
}



inline void
divide (bigrational & c, const bigrational & a, const bigint & b)
{
	c.num.assign(a.num);
	multiply(c.den, a.den, b);
	c.normalize();
}



inline void
divide (bigrational & c, const bigrational & a, long b)
{
	c.num.assign(a.num);
	multiply(c.den, a.den, b);
	c.normalize();
}



inline void
divide (bigrational & c, const bigrational & a, unsigned long b)
{
	c.num.assign(a.num);
	multiply(c.den, a.den, b);
	c.normalize();
}



inline void
divide (bigrational & c, const bigrational & a, int b)
{
	divide(c, a, static_cast<long>(b));
}



void power(bigrational & c, const bigrational & a, const bigint & b);
void power(bigrational & c, const bigrational & a, long b);

bigint round(const bigrational & a);
bigint floor(const bigrational & a);
bigint ceiling(const bigrational & a);
bigint truncate(const bigrational & a);

bool square_root(bigrational & root, const bigrational & s);
bool cube_root(bigrational & root, const bigrational & s);



//
// functions
//

inline bigrational
abs (const bigrational & a)
{
	bigrational c(a);

	c.absolute_value();
	return c;
}



inline bigrational
inverse (const bigrational & a)
{
	bigrational c(a);

	c.invert();
	return c;
}



//
// operators
//

inline bigrational
operator - (const bigrational & a)
{
	bigrational c(a);

	c.negate();
	return c;
}



inline bigrational
operator + (const bigrational & a, const bigrational & b)
{
	bigrational c;

	add(c, a, b);
	return c;
}



inline bigrational
operator + (const bigrational & a, int b)
{
	bigrational c;

	add(c, a, static_cast<long>(b));
	return c;
}



inline bigrational
operator + (const bigrational & a, long b)
{
	bigrational c;

	add(c, a, b);
	return c;
}



inline bigrational
operator + (const bigrational & a, const unsigned long b)
{
	bigrational c;

	add(c, a, b);
	return c;
}



inline bigrational
operator + (const bigrational & a, const bigint & b)
{
	bigrational c;

	add(c, a, b);
	return c;
}



inline bigrational
operator + (int a, const bigrational & b)
{
	bigrational c;

	add(c, b, static_cast<long>(a));
	return c;
}



inline bigrational
operator + (long a, const bigrational & b)
{
	bigrational c;

	add(c, b, a);
	return c;
}



inline bigrational
operator + (unsigned long a, const bigrational & b)
{
	bigrational c;

	add(c, b, a);
	return c;
}



inline bigrational
operator + (const bigint & a, const bigrational & b)
{
	bigrational c;

	add(c, b, a);
	return c;
}



inline bigrational
operator - (const bigrational & a, const bigrational & b)
{
	bigrational c;

	subtract(c, a, b);
	return c;
}



inline bigrational
operator - (const bigrational & a, int b)
{
	bigrational c;

	subtract(c, a, static_cast<long>(b));
	return c;
}



inline bigrational
operator - (const bigrational & a, long b)
{
	bigrational c;

	subtract(c, a, b);
	return c;
}



inline bigrational
operator - (const bigrational & a, unsigned long b)
{
	bigrational c;

	subtract(c, a, b);
	return c;
}



inline bigrational
operator - (const bigrational & a, const bigint & b)
{
	bigrational c;

	subtract(c, a, b);
	return c;
}



inline bigrational
operator - (int a, const bigrational & b)
{
	bigrational c;

	subtract(c, b, static_cast<long>(a));
	c.negate();
	return c;
}



inline bigrational
operator - (long a, const bigrational & b)
{
	bigrational c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigrational
operator - (unsigned long a, const bigrational & b)
{
	bigrational c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigrational
operator - (const bigint & a, const bigrational & b)
{
	bigrational c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigrational
operator * (const bigrational & a, const bigrational & b)
{
	bigrational c;

	multiply(c, a, b);
	return c;
}



inline bigrational
operator * (const bigrational & a, int b)
{
	bigrational c;

	multiply(c, a, static_cast<long>(b));
	return c;
}



inline bigrational
operator * (const bigrational & a, long b)
{
	bigrational c;

	multiply(c, a, b);
	return c;
}



inline bigrational
operator * (const bigrational & a, unsigned long b)
{
	bigrational c;

	multiply(c, a, b);
	return c;
}



inline bigrational
operator * (const bigrational & a, const bigint & b)
{
	bigrational c;

	multiply(c, a, b);
	return c;
}



inline bigrational
operator * (int a, const bigrational & b)
{
	bigrational c;

	multiply(c, b, static_cast<long>(a));
	return c;
}



inline bigrational
operator * (long a, const bigrational & b)
{
	bigrational c;

	multiply(c, b, a);
	return c;
}



inline bigrational
operator * (unsigned long a, const bigrational & b)
{
	bigrational c;

	multiply(c, b, a);
	return c;
}



inline bigrational
operator * (const bigint & a, const bigrational & b)
{
	bigrational c;

	multiply(c, b, a);
	return c;
}



inline bigrational
operator / (const bigrational & a, const bigrational & b)
{
	bigrational c;

	divide(c, a, b);
	return c;
}



inline bigrational
operator / (const bigrational & a, int b)
{
	bigrational c;

	divide(c, a, static_cast<long>(b));
	return c;
}



inline bigrational
operator / (const bigrational & a, long b)
{
	bigrational c;

	divide(c, a, b);
	return c;
}



inline bigrational
operator / (const bigrational & a, unsigned long b)
{
	bigrational c;

	divide(c, a, b);
	return c;
}



inline bigrational
operator / (const bigrational & a, const bigint & b)
{
	bigrational c;

	divide(c, a, b);
	return c;
}



inline bigrational
operator / (int a, const bigrational & b)
{
	bigrational c;

        if(a == 0) {
          if(b.is_zero()) {
            lidia_error_handler("bigrational", "operator/()::division by zero.");
          } 
          c.assign_zero();
        }
        else {
          divide(c, b, static_cast<long>(a));
          c.invert();
        }
	return c;
}



inline bigrational
operator / (long a, const bigrational & b)
{
	bigrational c;

        if(a == 0) {
          if(b.is_zero()) {
            lidia_error_handler("bigrational", "operator/()::division by zero.");
          } 
          c.assign_zero();
        }
        else {
          divide(c, b, a);
          c.invert();
        }
	return c;
}



inline bigrational
operator / (unsigned long a, const bigrational & b)
{
	bigrational c;

        if(a == 0) {
          if(b.is_zero()) {
            lidia_error_handler("bigrational", "operator/()::division by zero.");
          } 
          c.assign_zero();
        }
        else {          
          divide(c, b, a);
          c.invert();
        }
	return c;
}



inline bigrational
operator / (const bigint & a, const bigrational & b)
{
	bigrational c;

        if(a == 0) {
          if(b.is_zero()) {
            lidia_error_handler("bigrational", "operator/()::division by zero.");
          } 
          c.assign_zero();
        }
        else {          
          divide(c, b, a);
          c.invert();
        }
	return c;
}



inline bigrational
operator << (const bigrational & a, unsigned long ui)
{
	bigrational c;

	shift_left(c, a, ui);
	return c;
}



inline bigrational
operator >> (const bigrational & a, unsigned long ui)
{
	bigrational c;

	shift_right(c, a, ui);
	return c;
}



inline bigrational &
operator += (bigrational & a, const bigrational & b)
{
	add(a, a, b);
	return a;
}



inline bigrational &
operator += (bigrational & a, int b)
{
	add(a, a, static_cast<long>(b));
	return a;
}



inline bigrational &
operator += (bigrational & a, long b)
{
	add(a, a, b);
	return a;
}



inline bigrational &
operator += (bigrational & a, unsigned long b)
{
	add(a, a, b);
	return a;
}



inline bigrational &
operator += (bigrational & a, const bigint & b)
{
	add(a, a, b);
	return a;
}



inline bigrational &
operator -= (bigrational & a, const bigrational & b)
{
	subtract(a, a, b);
	return a;
}



inline bigrational &
operator -= (bigrational & a, int b)
{
	subtract(a, a, static_cast<long>(b));
	return a;
}



inline bigrational &
operator -= (bigrational & a, long b)
{
	subtract(a, a, b);
	return a;
}



inline bigrational &
operator -= (bigrational & a, unsigned long b)
{
	subtract(a, a, b);
	return a;
}



inline bigrational &
operator -= (bigrational & a, const bigint & b)
{
	subtract(a, a, b);
	return a;
}



inline bigrational &
operator *= (bigrational & a, const bigrational & b)
{
	multiply(a, a, b);
	return a;
}



inline bigrational &
operator *= (bigrational & a, int b)
{
	multiply(a, a, static_cast<long>(b));
	return a;
}



inline bigrational &
operator *= (bigrational & a, long b)
{
	multiply(a, a, b);
	return a;
}



inline bigrational &
operator *= (bigrational & a, unsigned long b)
{
	multiply(a, a, b);
	return a;
}



inline bigrational &
operator *= (bigrational & a, const bigint & b)
{
	multiply(a, a, b);
	return a;
}



inline bigrational &
operator /= (bigrational & a, const bigrational & b)
{
	divide(a, a, b);
	return a;
}



inline bigrational &
operator /= (bigrational & a, int b)
{
	divide(a, a, static_cast<long>(b));
	return a;
}



inline bigrational &
operator /= (bigrational & a, long b)
{
	divide(a, a, b);
	return a;
}



inline bigrational &
operator /= (bigrational & a, unsigned long b)
{
	divide(a, a, b);
	return a;
}



inline bigrational &
operator /= (bigrational & a, const bigint & b)
{
	divide(a, a, b);
	return a;
}



inline bigrational &
operator <<= (bigrational & a, unsigned long ui)
{
	shift_left(a, a, ui);
	return a;
}



inline bigrational &
operator >>= (bigrational & a, unsigned long ui)
{
	shift_right(a, a, ui);
	return a;
}



inline bigrational &
operator ++ (bigrational & a)
{
	a.inc();
	return a;
}



inline bigrational
operator ++ (bigrational & a, int)
{
	bigrational c(a);

	a.inc();
	return c;
}



inline bigrational &
operator -- (bigrational & a)
{
	a.dec();
	return a;
}



inline bigrational
operator -- (bigrational & a, int)
{
	bigrational c(a);

	a.dec();
	return c;
}



inline bool
operator ! (const bigrational & a)
{
	return a.is_zero();
}



//
// procedural modifiers
//

inline void
multiply_by_two (bigrational & c, const bigrational & a)
{
	c.assign(a);
	c.multiply_by_2();
}



inline void
divide_by_two (bigrational & c, const bigrational & a)
{
	c.assign(a);
	c.divide_by_2();
}



inline void
inc (bigrational & c)
{
	c.inc();
}



inline void
dec (bigrational & c)
{
	c.dec();
}



//
// I/O and converters
//

std::istream & operator >> (std::istream & in, bigrational & a);
std::ostream & operator << (std::ostream & out, const bigrational & a);

int string_to_bigrational(char *s, char *t, bigrational & a);
int bigrational_to_string(const bigrational & a, char *s, char *t);



//
// some properties of Q
//

inline bigint
bigrational::characteristic ()
{
	return 0UL;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_BIGRATIONAL

#include	"LiDIA/specialization/bigrational.special"



#endif	// LIDIA_BIGRATIONAL_H_GUARD_
