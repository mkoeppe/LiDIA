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


#ifndef LIDIA_BIGINT_H_GUARD_
#define LIDIA_BIGINT_H_GUARD_

#include <cstring>

#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif

#include	"LiDIA/kernel/bigint_def.h"

#ifndef LIDIA_XDOUBLE_H_GUARD_
# include	"LiDIA/xdouble.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint
{
private:

	//
	// the C type we use to represent a bigint
	//

	bigint_rep_t	I;


public:

	// The next 2 functions are implemented in bigint_share.cc
	// and exist only once for all interfaces.

	static void set_chars_per_line (int cpl);
	static int  get_chars_per_line ();


	//
	// c'tors and d'tor
	//

	bigint();
	bigint(int i);
	bigint(long l);
	bigint(unsigned long ul);
	bigint(double d);
	explicit bigint(const bigint_rep_t & a);
	bigint(const bigint & a);
	~bigint();



	//
	// accessors
	//

	const bigint_rep_t &	bigint_rep() const;

	int bit(unsigned long i) const;
	unsigned long most_significant_digit() const;
	unsigned long least_significant_digit() const;



	//
	// assigners
	//

	void assign_zero();
	void assign_one();

	void assign(int i);
	void assign(long l);
	void assign(unsigned long ul);
	void assign(double d);
	void assign(const bigint_rep_t & a);
	void assign(const bigint & a);

	bigint & operator = (int i);
	bigint & operator = (long l);
	bigint & operator = (unsigned long ul);
	bigint & operator = (double d);
	bigint & operator = (const bigint_rep_t & a);
	bigint & operator = (const bigint & a);


	//
	// comparators
	//

	int abs_compare(const bigint & a) const;
	int abs_compare(unsigned long a) const;
	int compare(const bigint & a) const;
	int compare(unsigned long a) const;
	int compare(long a) const;
	
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
	// properties
	//

	lidia_size_t length() const;
	lidia_size_t bit_length() const;
	lidia_size_t decimal_length() const;



	//
	// predicates
	//

	bool is_odd() const;
	bool is_even() const;

	bool is_char() const;
	bool is_uchar() const;
	bool is_short() const;
	bool is_ushort() const;
	bool is_int() const;
	bool is_uint() const;
	bool is_long() const;
	bool is_ulong() const;

	bool is_prime(int k = 10) const;
	bool is_square() const;
	bool is_square(bigint & r) const;
	long is_power(bigint & b) const;



#ifndef HEADBANGER

	//
	// converters
	//

	bool intify(int & i) const;
	bool longify(long & l) const;

	double dbl() const;
	xdouble xdbl() const;



	//
	// functions
	//

	bigint next_prime() const;
	bigint previous_prime() const;



	//
	// modifiers
	//

	void absolute_value();
	void abs();
	void negate();

	void inc();
	void dec();

	void multiply_by_2();
	void divide_by_2();

	void swap(bigint & a);



	//
	// procedural arithmetics
	//

	friend void negate(bigint & a, const bigint & b);
	friend void add(bigint & c, const bigint & a, const bigint & b);
	friend void add(bigint & c, const bigint & a, long b);
	friend void add(bigint & c, const bigint & a, unsigned long b);
	friend void subtract(bigint & c, const bigint & a, const bigint & b);
	friend void subtract(bigint & c, const bigint & a, long b);
	friend void subtract(bigint & c, const bigint & a, unsigned long b);

	friend void multiply(bigint & c, const bigint & a, const bigint & b);
	friend void multiply(bigint & c, const bigint & a, long b);
	friend void multiply(bigint & c, const bigint & a, unsigned long b);
	friend void square(bigint & a, const bigint & b);
	friend void invert(bigint & a, const bigint & b);
	friend void divide(bigint & c, const bigint & a, const bigint & b);
	friend void divide(bigint & c, const bigint & a, long b);
	friend void divide(bigint & c, const bigint & a, unsigned long b);
	friend void remainder(bigint & c, const bigint & a, const bigint & b);
	friend void remainder(bigint & c, const bigint & a, long b);	
	friend void remainder(bigint & c, const bigint & a, unsigned long b);	
	friend void remainder(long & c, const bigint & a, long b);
	friend void remainder(unsigned long & c, const bigint & a, unsigned long b);
	friend long remainder(const bigint & a, long b);	
	friend void div_rem(bigint & q, bigint & r, const bigint & a, const bigint & b);
	friend void div_rem(bigint & q, bigint & r, const bigint & a, long b);
	friend void div_rem(bigint & q, bigint & r, const bigint & a, unsigned long b);
	friend void div_rem(bigint & q, long & r, const bigint & a, long b);

	friend void power(bigint & c, const bigint & a, const bigint & b);
	friend void power(bigint & c, const bigint & a, long i);

	friend void sqrt(bigint & a, const bigint & b);

	friend void shift_left(bigint & c, const bigint & a, long ui);
	friend void shift_right(bigint & c, const bigint & a, long ui);

	friend void bitwise_and(bigint & c, const bigint & a, const bigint & b);
	friend void bitwise_or(bigint & c, const bigint & a, const bigint & b);
	friend void bitwise_xor(bigint & c, const bigint & a, const bigint & b);
	friend void bitwise_not(bigint & c, const bigint & a);

#endif	// HEADBANGER



	//
	// gcd's
	//

	friend bigint gcd(const bigint & a, const bigint & b);
	friend bigint bgcd(const bigint & a, const bigint & b);
	friend bigint dgcd(const bigint & a, const bigint & b);
	friend bigint xgcd(bigint & u, bigint & v, const bigint & a, const bigint & b);
	friend bigint xgcd_left(bigint & u, const bigint & a, const bigint & b);
	friend bigint xgcd_right(bigint & v, const bigint & a, const bigint & b);

	friend bigint lcm(const bigint& a, const bigint & b);



	//
	// other arithmetical functions
	//

	friend bigint abs(const bigint & a);



	//
	// random numbers
	//

	static bool is_seeded;

	static void seed(const bigint & a);
	static void seed();

	void randomize(const bigint & a);



	//
	// input / output
	//

	friend std::istream & operator >> (std::istream & in, bigint & a);
	friend std::ostream & operator << (std::ostream & out, const bigint & a);

	friend int string_to_bigint(const char *s, bigint & a);
	friend int bigint_to_string(const bigint & a, char *s);



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



	//
	// misc. functions
	//

	static bigint characteristic ();



	//
	// the next two definitions are needed by bigfloat
	//

	//static const unsigned long radix();
	static double radix();
	static int bits_per_digit();


private:

	//
	//  input / output  utilities, implemented in bigint_share.cc
	//

	static int chars_per_line; // number of characters in one line
                                // when printing a bigint


	// The next 5 functions are implemented in bigint_share.cc
	// and exist only once for all interfaces.

	static void allocate (char * &s, int old_size, int new_size);
	static void append_char (char * &s, int & sz, int pos, char c);
	static int  skip_backslash_new_line (std::istream & in);

	void scan  (std::istream & in);
	void print (std::ostream & out, char const *s) const;

};



#if INLINE_INTERFACE
# include	"LiDIA/kernel/bigint_interface.h"
#endif

//
// procedural arithmetics
//

void negate(bigint & a, const bigint & b);
void add(bigint & c, const bigint & a, const bigint & b);
void add(bigint & c, const bigint & a, long b);
void add(bigint & c, const bigint & a, unsigned long b);
void subtract(bigint & c, const bigint & a, const bigint & b);
void subtract(bigint & c, const bigint & a, long b);
void subtract(bigint & c, const bigint & a, unsigned long b);

void multiply(bigint & c, const bigint & a, const bigint & b);
void multiply(bigint & c, const bigint & a, long b);
void multiply(bigint & c, const bigint & a, unsigned long b);
void square(bigint & a, const bigint & b);
void invert(bigint & a, const bigint & b);
void divide(bigint & c, const bigint & a, const bigint & b);
void divide(bigint & c, const bigint & a, long b);
void divide(bigint & c, const bigint & a, unsigned long b);
void remainder(bigint & c, const bigint & a, const bigint & b);
void remainder(bigint & c, const bigint & a, long b);	
void remainder(bigint & c, const bigint & a, unsigned long b);	
void remainder(long & c, const bigint & a, long b);
void remainder(unsigned long & c, const bigint & a, unsigned long b);
long remainder(const bigint & a, long b);	
void div_rem(bigint & q, bigint & r, const bigint & a, const bigint & b);
void div_rem(bigint & q, bigint & r, const bigint & a, long b);
void div_rem(bigint & q, bigint & r, const bigint & a, unsigned long b);
void div_rem(bigint & q, long & r, const bigint & a, long b);

void power(bigint & c, const bigint & a, const bigint & b);
void power(bigint & c, const bigint & a, long i);

void sqrt(bigint & a, const bigint & b);

void shift_left(bigint & c, const bigint & a, long ui);
void shift_right(bigint & c, const bigint & a, long ui);

void bitwise_and(bigint & c, const bigint & a, const bigint & b);
void bitwise_or(bigint & c, const bigint & a, const bigint & b);
void bitwise_xor(bigint & c, const bigint & a, const bigint & b);
void bitwise_not(bigint & c, const bigint & a);

//
// gcd's
//

bigint gcd(const bigint & a, const bigint & b);
bigint bgcd(const bigint & a, const bigint & b);
bigint dgcd(const bigint & a, const bigint & b);
bigint xgcd(bigint & u, bigint & v, const bigint & a, const bigint & b);
bigint xgcd_left(bigint & u, const bigint & a, const bigint & b);
bigint xgcd_right(bigint & v, const bigint & a, const bigint & b);

bigint lcm(const bigint& a, const bigint & b);


//
// other arithmetical functions
//

bigint abs(const bigint & a);


//
// input / output
//

std::istream & operator >> (std::istream & in, bigint & a);
std::ostream & operator << (std::ostream & out, const bigint & a);

int string_to_bigint(const char *s, bigint & a);
int bigint_to_string(const bigint & a, char *s);


//
// accessors
//

inline const bigint_rep_t &
bigint::bigint_rep () const
{
	return I;
}



//
// assigners
//

inline bigint &
bigint::operator = (int i)
{
	assign(i);
	return *this;
}



inline bigint &
bigint::operator = (long l)
{
	assign(l);
	return *this;
}



inline bigint &
bigint::operator = (unsigned long ul)
{
	assign(ul);
	return *this;;
}



inline bigint &
bigint::operator = (double d)
{
	assign(d);
	return *this;
}



inline bigint &
bigint::operator = (const bigint_rep_t & a)
{
	assign(a);
	return *this;
}



inline bigint &
bigint::operator = (const bigint & a)
{
	assign(a);
	return *this;
}




//
// comparators
//

inline bool
bigint::is_positive () const
{
	return (sign() > 0);
}



inline bool
bigint::is_negative () const
{
	return (sign() < 0);
}



inline bool
bigint::is_zero () const
{
	return (sign() == 0);
}



inline bool
bigint::is_gt_zero () const
{
	return (sign() > 0);
}



inline bool
bigint::is_ge_zero () const
{
	return (sign() >= 0);
}



inline bool
bigint::is_lt_zero () const
{
	return (sign() < 0);
}



inline bool
bigint::is_le_zero () const
{
	return (sign() <= 0);
}



inline bool
bigint::is_one () const
{
	return (compare(1L) == 0);
}



//
// operational comparators
//

inline bool
operator == (const bigint & a, const bigint & b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigint & a, unsigned long b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigint & a, long b)
{
	return (a.compare(b) == 0);
}



inline bool
operator == (const bigint & a, int b)
{
	return (a.compare(static_cast<long>(b)) == 0);
}



inline bool
operator == (unsigned long a, const bigint & b)
{
	return (b.compare(a) == 0);
}



inline bool
operator == (long a, const bigint & b)
{
	return (b.compare(a) == 0);
}



inline bool
operator == (int a, const bigint & b)
{
	return (b.compare(static_cast<long>(a)) == 0);
}



inline bool
operator != (const bigint & a, const bigint & b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigint & a, unsigned long b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigint & a, long b)
{
	return (a.compare(b) != 0);
}



inline bool
operator != (const bigint & a, int b)
{
	return (a.compare(static_cast<long>(b)) != 0);
}



inline bool
operator != (unsigned long a, const bigint & b)
{
	return (b.compare(a) != 0);
}



inline bool
operator != (long a, const bigint & b)
{
	return (b.compare(a) != 0);
}



inline bool
operator != (int a, const bigint & b)
{
	return (b.compare(static_cast<long>(a)) != 0);
}



inline bool
operator > (const bigint & a, const bigint & b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigint & a, unsigned long b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigint & a, long b)
{
	return (a.compare(b) > 0);
}



inline bool
operator > (const bigint & a, int b)
{
	return (a.compare(static_cast<long>(b)) > 0);
}



inline bool
operator > (unsigned long a, const bigint & b)
{
	return (b.compare(a) < 0);
}



inline bool
operator > (long a, const bigint & b)
{
	return (b.compare(a) < 0);
}



inline bool
operator > (int a, const bigint & b)
{
	return (b.compare(static_cast<long>(a)) < 0);
}



inline bool
operator >= (const bigint & a, const bigint & b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigint & a, unsigned long b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigint & a, long b)
{
	return (a.compare(b) >= 0);
}



inline bool
operator >= (const bigint & a, int b)
{
	return (a.compare(static_cast<long>(b)) >= 0);
}



inline bool
operator >= (unsigned long a, const bigint & b)
{
	return (b.compare(a) <= 0);
}



inline bool
operator >= (long a, const bigint & b)
{
	return (b.compare(a) <= 0);
}



inline bool
operator >= (int a, const bigint & b)
{
	return (b.compare(static_cast<long>(a)) <= 0);
}



inline bool
operator < (const bigint & a, const bigint & b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigint & a, unsigned long b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigint & a, long b)
{
	return (a.compare(b) < 0);
}



inline bool
operator < (const bigint & a, int b)
{
	return (a.compare(static_cast<long>(b)) < 0);
}



inline bool
operator < (unsigned long a, const bigint & b)
{
	return (b.compare(a) > 0);
}



inline bool
operator < (long a, const bigint & b)
{
	return (b.compare(a) > 0);
}



inline bool
operator < (int a, const bigint & b)
{
	return (b.compare(static_cast<long>(a)) > 0);
}



inline bool
operator <= (const bigint & a, const bigint & b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigint & a, unsigned long b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigint & a, long b)
{
	return (a.compare(b) <= 0);
}



inline bool
operator <= (const bigint & a, int b)
{
	return (a.compare(static_cast<long>(b)) <= 0);
}



inline bool
operator <= (unsigned long a, const bigint & b)
{
	return (b.compare(a) >= 0);
}



inline bool
operator <= (long a, const bigint & b)
{
	return (b.compare(a) >= 0);
}



inline bool
operator <= (int a, const bigint & b)
{
	return (b.compare(static_cast<long>(a)) >= 0);
}



//
// functional comparators
//

inline int
abs_compare (const bigint& a, const bigint& b)
{
	return a.abs_compare(b);
}



inline int
abs_compare (const bigint & a, unsigned long b)
{
	return a.abs_compare(b);
}



inline int
compare (const bigint& a, const bigint& b)
{
	return a.compare(b);
}



inline int
compare (const bigint & a, unsigned long b)
{
	return a.compare(b);
}



inline int
compare (const bigint & a, long b)
{
	return a.compare(b);
}



inline bool
is_positive (const bigint & a)
{
	return (a.sign() > 0);
}



inline bool
is_negative (const bigint & a)
{
	return (a.sign() < 0);
}



inline bool
is_zero (const bigint & a)
{
	return (a.sign() == 0);
}



inline bool
is_one (const bigint & a)
{
	return (a.compare(1UL) == 0);
}



//
// functional predicates
//

inline bool
is_even (const bigint & a)
{
	return a.is_even();
}



inline bool
is_odd (const bigint & a)
{
	return a.is_odd();
}



inline bool
is_char (const bigint & a)
{
	return a.is_char();
}



inline bool
is_uchar (const bigint & a)
{
	return a.is_uchar();
}



inline bool
is_short (const bigint & a)
{
	return a.is_short();
}



inline bool
is_ushort (const bigint & a)
{
	return a.is_ushort();
}



inline bool
is_int (const bigint & a)
{
	return a.is_int();
}



inline bool
is_uint (const bigint & a)
{
	return a.is_uint();
}



inline bool
is_long (const bigint & a)
{
	return a.is_long();
}



inline bool
is_ulong (const bigint & a)
{
	return a.is_ulong();
}


inline bool
is_prime (const bigint & a, int k = 10)
{
	return a.is_prime(k);
}



inline bool
is_square (bigint & r, const bigint & a)
{
	return a.is_square(r);
}



inline bool
is_square (const bigint & a)
{
	return a.is_square();
}



inline long
is_power (bigint & b, const bigint & a)
{
	return a.is_power(b);
}



//
// properties
//

inline lidia_size_t
bigint::decimal_length () const
{
  if (is_zero())
    return 1;

  lidia_size_t r = 0;
  bigint a(*this);
 
  a.abs();

  while(!a.is_zero())
    {
      divide(a, a, static_cast<unsigned long>(10));
      r++;
    }

  return r;
}



inline lidia_size_t
decimal_length (const bigint & a)
{
  return a.decimal_length();
}



//
// procedural modifiers
//

inline void
swap (bigint & a, bigint & b)
{
	a.swap(b);
}



inline void
inc (bigint & a)
{
	a.inc();
}



inline void
dec (bigint & a)
{
	a.dec();
}



//
// arithmetic procedures
//

inline void
add (bigint & c, const bigint & a, int b)
{
	add(c, a, static_cast<long>(b));
}



inline void
subtract (bigint & c, const bigint & a, int b)
{
	subtract(c, a, static_cast<long>(b));
}



inline void
multiply (bigint & c, const bigint & a, int b)
{
	multiply(c, a, static_cast<long>(b));
}



inline void
divide (bigint & c, const bigint & a, int b)
{
	divide(c, a, static_cast<long>(b));
}



inline void
remainder (bigint & c, const bigint & a, int b)
{
	remainder(c, a, static_cast<long>(b));
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, int b)
{
	div_rem(q, r, a, static_cast<long>(b));
}



//
// arithmetic operators
//

inline bigint
operator - (const bigint & a)
{
	bigint c;

	negate(c, a);
	return c;
}



inline bigint
operator + (const bigint & a, const bigint & b)
{
	bigint c;

	add(c, a, b);
	return c;
}



inline bigint
operator + (const bigint & a, long b)
{
	bigint c;

	add(c, a, b);
	return c;
}



inline bigint
operator + (const bigint & a, unsigned long b)
{
	bigint c;

	add(c, a, b);
	return c;
}



inline bigint
operator + (const bigint & a, int b)
{
	bigint c;

	add(c, a, static_cast<long>(b));
	return c;
}



inline bigint
operator + (long a, const bigint & b)
{
	bigint c;

	add(c, b, a);
	return c;
}



inline bigint
operator + (unsigned long a, const bigint & b)
{
	bigint c;

	add(c, b, a);
	return c;
}



inline bigint
operator + (int a, const bigint & b)
{
	bigint c;

	add(c, b, static_cast<long>(a));
	return c;
}



inline bigint
operator - (const bigint & a, const bigint & b)
{
	bigint c;

	subtract(c, a, b);
	return c;
}



inline bigint
operator - (const bigint & a, long b)
{
	bigint c;

	subtract(c, a, b);
	return c;
}



inline bigint
operator - (const bigint & a, unsigned long b)
{
	bigint c;

	subtract(c, a, b);
	return c;
}



inline bigint
operator - (const bigint & a, int b)
{
	bigint c;

	subtract(c, a, static_cast<long>(b));
	return c;
}



inline bigint
operator - (long a, const bigint & b)
{
	bigint c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigint
operator - (unsigned long a, const bigint & b)
{
	bigint c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigint
operator - (int a, const bigint & b)
{
	bigint c;

	subtract(c, b, static_cast<long>(a));
	c.negate();
	return c;
}



inline bigint
operator * (const bigint & a, const bigint & b)
{
	bigint c;

	multiply(c, a, b);
	return c;
}



inline bigint
operator * (const bigint & a, long b)
{
	bigint c;

	multiply(c, a, b);
	return c;
}



inline bigint
operator * (const bigint & a, unsigned long b)
{
	bigint c;

	multiply(c, a, b);
	return c;
}



inline bigint
operator * (const bigint & a, int b)
{
	bigint c;

	multiply(c, a, static_cast<long>(b));
	return c;
}



inline bigint
operator * (long a, const bigint & b)
{
	bigint c;

	multiply(c, b, a);
	return c;
}



inline bigint
operator * (unsigned long a, const bigint & b)
{
	bigint c;

	multiply(c, b, a);
	return c;
}



inline bigint
operator * (int a, const bigint & b)
{
	bigint c;

	multiply(c, b, static_cast<long>(a));
	return c;
}



inline bigint
operator / (const bigint & a, const bigint & b)
{
	bigint c;

	divide(c, a, b);
	return c;
}



inline bigint
operator / (const bigint & a, unsigned long b)
{
	bigint c;

	divide(c, a, b);
	return c;
}



inline bigint
operator / (const bigint & a, long b)
{
	bigint c;

	divide(c, a, b);
	return c;
}



inline bigint
operator / (const bigint & a, int b)
{
	bigint c;

	divide(c, a, static_cast<long>(b));
	return c;
}



inline bigint
operator / (long a, const bigint & b)
{
	bigint c;

	divide(c, bigint(a), b);
	return c;
}



inline bigint
operator / (unsigned long a, const bigint & b)
{
	bigint c;

	divide(c, bigint(a), b);
	return c;
}



inline bigint
operator / (int a, const bigint & b)
{
	bigint c;

	divide(c, bigint(static_cast<long>(a)), b);
	return c;
}



inline bigint
operator % (const bigint & a, const bigint & b)
{
	bigint c;

	remainder(c, a, b);
	return c;
}



inline bigint
operator % (const bigint & a, long b)
{
	bigint c;

	remainder(c, a, b);
	return c;
}



inline bigint
operator % (const bigint & a, unsigned long b)
{
	bigint c;

	remainder(c, a, b);
	return c;
}



inline bigint
operator % (const bigint & a, int b)
{
	bigint c;

	remainder(c, a, static_cast<long>(b));
	return c;
}



inline bigint
operator % (long a, const bigint & b)
{
	bigint c;

	remainder(c, bigint(a), b);
	return c;
}



inline bigint
operator % (unsigned long a, const bigint & b)
{
	bigint c;

	remainder(c, bigint(a), b);
	return c;
}



inline bigint
operator % (int a, const bigint & b)
{
	bigint c;

	remainder(c, bigint(static_cast<long>(a)), b);
	return c;
}



inline bigint
operator << (const bigint & a, long ui)
{
	bigint c;

	if (ui < 0) {
		lidia_error_handler("bigint", "operator << ::index is negative.");
	}
	shift_left(c, a, ui);
	return c;
}



inline bigint
operator >> (const bigint & a, long ui)
{
	bigint c;

	if (ui < 0)
		lidia_error_handler("bigint", "operator >>::index is negative.");
	shift_right(c, a, ui);
	return c;
}



inline bigint
operator & (const bigint & a, const bigint & b)
{
	bigint c;

	bitwise_and(c, a, b);
	return c;
}



inline bigint
operator | (const bigint & a, const bigint & b)
{
	bigint c;

	bitwise_or(c, a, b);
	return c;
}



inline bigint
operator ^ (const bigint & a, const bigint & b)
{
	bigint c;

	bitwise_xor(c, a, b);
	return c;
}



inline bigint &
operator += (bigint & a, const bigint & b)
{
	add(a, a, b);
	return a;
}



inline bigint &
operator += (bigint & a, long b)
{
	add(a, a, b);
	return a;
}



inline bigint &
operator += (bigint & a, unsigned long b)
{
	add(a, a, b);
	return a;
}



inline bigint &
operator += (bigint & a, int b)
{
	add(a, a, static_cast<long>(b));
	return a;
}



inline bigint &
operator -= (bigint & a, const bigint & b)
{
	subtract(a, a, b);
	return a;
}



inline bigint &
operator -= (bigint & a, long b)
{
	subtract(a, a, b);
	return a;
}



inline bigint &
operator -= (bigint & a, unsigned long b)
{
	subtract(a, a, b);
	return a;
}



inline bigint &
operator -= (bigint & a, int b)
{
	subtract(a, a, static_cast<long>(b));
	return a;
}



inline bigint &
operator *= (bigint & a, const bigint & b)
{
	multiply(a, a, b);
	return a;
}



inline bigint &
operator *= (bigint & a, long b)
{
	multiply(a, a, b);
	return a;
}



inline bigint &
operator *= (bigint & a, unsigned long b)
{
	multiply(a, a, b);
	return a;
}



inline bigint &
operator *= (bigint & a, int b)
{
	multiply(a, a, static_cast<long>(b));
	return a;
}



inline bigint &
operator /= (bigint & a, const bigint & b)
{
	divide(a, a, b);
	return a;
}



inline bigint &
operator /= (bigint & a, long b)
{
	divide(a, a, b);
	return a;
}



inline bigint &
operator /= (bigint & a, unsigned long b)
{
	divide(a, a, b);
	return a;
}



inline bigint &
operator /= (bigint & a, int b)
{
	divide(a, a, static_cast<long>(b));
	return a;
}



inline bigint &
operator %= (bigint & a, const bigint & b)
{
	if (&a == &b) {
		a.assign_zero();
	}
	else {
		remainder(a, a, b);
	}
	return a;
}



inline bigint &
operator %= (bigint & a, long b)
{
	remainder(a, a, b);
	return a;
}



inline bigint &
operator %= (bigint & a, unsigned long b)
{
	remainder(a, a, b);
	return a;
}



inline bigint &
operator %= (bigint & a, int b)
{
	remainder(a, a, static_cast<long>(b));
	return a;
}



inline bigint &
operator <<= (bigint & a, long ui)
{
	shift_left(a, a, ui);
	return a;
}



inline bigint &
operator >>= (bigint & a, long ui)
{
	shift_right(a, a, ui);
	return a;
}



inline bigint &
operator &= (bigint & a, const bigint & b)
{
	bitwise_and(a, a, b);
	return a;
}



inline bigint &
operator |= (bigint & a, const bigint & b)
{
	bitwise_or(a, a, b);
	return a;
}



inline bigint &
operator ^= (bigint & a, const bigint & b)
{
	bitwise_xor(a, a, b);
	return a;
}



inline bigint &
operator ++ (bigint & a)
{
	a.inc();
	return a;
}



inline bigint &
operator -- (bigint & a)
{
	a.dec();
	return a;
}



inline bigint
operator ++ (bigint & a, int)
{
	bigint c(a);

	a.inc();
	return c;
}



inline bigint
operator -- (bigint & a, int)
{
	bigint c(a);

	a.dec();
	return c;
}



inline bool
operator ! (const bigint & a)
{
	return (a.sign() == 0);
}



inline bigint
operator ~ (const bigint & a)
{
	bigint c;

	bitwise_not(c, a);
	return c;
}



//
// functional converters
//

inline double
dbl(const bigint & a)
{
	return a.dbl();
}



inline xdouble
xdbl(const bigint & a)
{
	return a.xdbl();
}



//
// functions
//

inline bigint
abs (const bigint & a)
{
	bigint c(a);

	c.absolute_value();
	return c;
}



inline bigint
next_prime (const bigint & a)
{
	return a.next_prime();
}



inline bigint
previous_prime (const bigint & a)
{
	return a.previous_prime();
}



//
// misc. functions
//

inline bigint
randomize(const bigint & a)
{
	bigint c;

	c.randomize(a);
	return c;
}



inline bigint
bigint::characteristic ()
{
	return 0UL;
}



//
// class accesors and modfiers
//

inline void
bigint::set_chars_per_line (int cpl)
{
	bigint::chars_per_line = cpl;
}



inline int
bigint ::get_chars_per_line ()
{
	return bigint::chars_per_line;
}



//
// other functions and procedures operating on bigint
//

bool fermat (const bigint & n);

bigint chinese_remainder(const bigint & x1, const bigint & m1,
			 const bigint & x2, const bigint & m2);

int jacobi (const bigint &a1, const bigint &b1);

void nearest (bigint & z, const bigint & u, const bigint & v);
void newton_root (bigint & b, const bigint & a, int n);

void power_mod (bigint & res, const bigint & a, const bigint & n, const bigint & m , int err = 0);
bigint * mgcd (const bigint * a, const int n);

void ressol   (bigint&, const bigint &, const bigint &);
void ressol_p (bigint&, const bigint &, const bigint &);
void ressol_prime_power(bigint & , const bigint & , 
						const bigint & , const int );

bool cornacchia (bigint & x, bigint & y, const bigint & DD, const bigint & p);
bool cornacchia_prime_power( bigint & x, bigint & y, const bigint & DD, 
							 const bigint & p, const int exp );


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_BIGINT

#include	"LiDIA/specialization/bigint.special"



#endif	// LIDIA_BIGINT_H_GUARD_
