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


#ifndef LIDIA_BIGMOD_H_GUARD_
#define LIDIA_BIGMOD_H_GUARD_



#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_BASE_BIGMOD_H_GUARD_
# include	"LiDIA/base/base_bigmod.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class multi_bigmod;



class bigmod : public base_bigmod
{
private :

#ifndef HEADBANGER
	static LiDIA::residue_class< bigint > *Mp;
#endif
	static bigint M;

public :

	//
	// class modifiers and accessors
	//

	static void set_modulus (const bigint & m);
	static const bigint & modulus ();
	static LiDIA::residue_class< bigint > * residue_class();


	//
	// class properties
	//

	static const bigint & characteristic ();



        //
	// c'tors and d'tor
	//

	bigmod ();
	bigmod (int i);
	bigmod (long l);
	bigmod (unsigned long ul);
	bigmod (double d);
	bigmod (const bigint & i);
	bigmod (const bigmod & a);
        virtual ~bigmod ();



#ifndef HEADBANGER

        //
        // assigners
        //

        void assign_zero ();
        void assign_one ();

        void assign (int i);
        void assign (long l);
        void assign (unsigned long ul);
	void assign (double d);
        void assign (const bigint & a);
        void assign (const bigmod & a);

#endif	// HEADBANGER

	bigmod & operator = (int i);
	bigmod & operator = (long l);
	bigmod & operator = (unsigned long ul);
	bigmod & operator = (double d);
	bigmod & operator = (const bigint & a);
	bigmod & operator = (const bigmod & a);



        //
	// comparators
	//

	bool is_equal(const bigmod &) const;
	bool is_equal(const bigint &) const;
	bool is_equal(long) const;
	bool is_equal(unsigned long) const;



        //
	// modifiers
	//

        bigint invert (int verbose = 0);
        void negate ();
	void normalize ();

	void inc();
	void dec();

	void multiply_by_2 ();
	void divide_by_2 ();

        void randomize ();

	void swap(bigmod & a);



	//
	// I/O
	//

        friend std::istream & operator >> (std::istream & in, bigmod & a);
        friend std::ostream & operator << (std::ostream & out, const bigmod & a);

	friend int string_to_bigmod (char * s, bigmod & a);
	friend int bigmod_to_string (const bigmod & a, char * s);

#ifdef C_STDIO

	void read_from_file (FILE * fp);
	void write_to_file  (FILE * fp);
	void scan_from_file (FILE * fp);
	void print_to_file  (FILE * fp);

#endif	// C_STDIO

};


//
// I/O
//

std::istream & operator >> (std::istream & in, bigmod & a);
std::ostream & operator << (std::ostream & out, const bigmod & a);

int string_to_bigmod (char * s, bigmod & a);
int bigmod_to_string (const bigmod & a, char * s);

//
// class accessors
//

inline const bigint &
bigmod::modulus ()
{
	return bigmod::M;
}



inline
LiDIA::residue_class< bigint > *
bigmod::residue_class ()
{
	return bigmod::Mp;
}



//
// normalizer (inline definition must precede c'tors and assigners)
//

inline void
bigmod::normalize ()
{
	base_bigmod::normalize(bigmod::M);
}



//
// c'tors and d'tor
//

inline
bigmod::bigmod ()
	: base_bigmod()
{
	// nothing to do
}



inline
bigmod::bigmod (int i)
	: base_bigmod(i)
{
	normalize();
}



inline
bigmod::bigmod (long l)
	: base_bigmod(l)
{
	normalize();
}



inline
bigmod::bigmod (unsigned long ul)
	: base_bigmod(ul)
{
	normalize();
}



inline
bigmod::bigmod (double d)
	: base_bigmod(d)
{
	normalize();
}



inline
bigmod::bigmod (const bigint & i)
	: base_bigmod(i)
{
	normalize();
}



inline
bigmod::bigmod (const bigmod & a)
	: base_bigmod(a.I)
{
	// nothing to do
}



inline
bigmod::~bigmod ()
{
}



//
// assigners
//

inline void
bigmod::assign_zero ()
{
	I.assign_zero();
}



inline void
bigmod::assign_one ()
{
	I.assign_one();
}



inline void
bigmod::assign (int i)
{
	I.assign(i);
	normalize();
}




inline void
bigmod::assign (long l)
{
	I.assign(l);
	normalize();
}



inline void
bigmod::assign (unsigned long ul)
{
	I.assign(ul);
	normalize();
}



inline void
bigmod::assign (double d)
{
	I.assign(d);
	normalize();
}



inline void
bigmod::assign (const bigint & a)
{
	I.assign(a);
	normalize();
}



inline void
bigmod::assign (const bigmod & a)
{
	if (&a != this) {
		I.assign(a.I);
	}
}



inline bigmod &
bigmod::operator = (int i)
{
	assign(i);
	return *this;
}



inline bigmod &
bigmod::operator = (long l)
{
	assign(l);
	return *this;
}



inline bigmod &
bigmod::operator = (unsigned long ul)
{
	assign(ul);
	return *this;
}



inline bigmod &
bigmod::operator = (double d)
{
	assign(d);
	return *this;
}



inline bigmod &
bigmod::operator = (const bigint & a)
{
	assign(a);
	return *this;
}



inline bigmod &
bigmod::operator = (const bigmod & a)
{
	assign(a);
	return *this;
}



//
// comparators
//

inline bool
bigmod::is_equal (const bigmod & a) const
{
	return (&a == this || I.compare(a.I) == 0);
}



inline bool
bigmod::is_equal (const bigint & a) const
{
	return a.is_ge_zero() ?
	    (I.compare(a) == 0) : (I.compare(M + a) == 0);
}



inline bool
bigmod::is_equal (long a) const
{
	return a >= 0 ?
	    (I.compare(a) == 0) : (I.compare(M + a) == 0);
}



inline bool
bigmod::is_equal (unsigned long a) const
{
	return I.compare(a) == 0;
}



inline bool
operator == (const bigmod & a, const bigmod & b)
{
	return a.is_equal(b);
}



inline bool
operator == (const bigmod & a, const bigint & b)
{
	return a.is_equal(b);
}



inline bool
operator == (const bigmod & a, long b)
{
	return a.is_equal(b);
}



inline bool
operator == (const bigmod & a, unsigned long b)
{
	return a.is_equal(b);
}



inline bool
operator == (const bigmod & a, int b)
{
	return a.is_equal(static_cast<long>(b));
}



inline bool
operator == (const bigint & a, const bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (long a, const bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (unsigned long a, const bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (int a, const bigmod & b)
{
	return b.is_equal(static_cast<long>(a));
}



inline bool
operator != (const bigmod & a, const bigmod & b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const bigmod & a, const bigint & b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const bigmod & a, long b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const bigmod & a, unsigned long b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const bigmod & a, int b)
{
	return !a.is_equal(static_cast<long>(b));
}



inline bool
operator != (const bigint & a, const bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (long a, const bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (unsigned long a, const bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (int a, const bigmod & b)
{
	return !b.is_equal(static_cast<long>(a));
}



//
// modifiers
//

inline bigint
bigmod::invert (int verbose)
{
	return base_bigmod::invert(bigmod::M, verbose);
}



inline void
bigmod::negate ()
{
	base_bigmod::negate(bigmod::M);
}



inline void
bigmod::inc ()
{
	base_bigmod::inc(bigmod::M);
}



inline void
bigmod::dec ()
{
	base_bigmod::dec(bigmod::M);
}



inline void
bigmod::multiply_by_2 ()
{
	base_bigmod::multiply_by_2(bigmod::M);
}



inline void
bigmod::divide_by_2 ()
{
	base_bigmod::divide_by_2(bigmod::M);
}



inline void
bigmod::swap (bigmod & a)
{
	I.swap(a.I);
}



//
// arithmetic procedures
//

inline void
normalize (bigmod & a, const bigmod & b)
{
	a.assign(b);
	a.normalize();
}



inline void
invert (bigmod & a, const bigmod & b)
{
	a.assign(b);
	a.invert();
}



inline void
negate (bigmod & a, const bigmod & b)
{
	a.assign(b);
	a.negate();
}



inline void
add (bigmod & c, const bigmod & a, const bigmod & b)
{
	add(c, a, b, bigmod::modulus());
}



inline void
add (bigmod & c, const bigmod & a, const bigint & b)
{
	add(c, a, b, bigmod::modulus());
}



inline void
add (bigmod & c, const bigmod & a, long b)
{
	add(c, a, b, bigmod::modulus());
}



inline void
add (bigmod & c, const bigmod & a, unsigned long b)
{
	add(c, a, b, bigmod::modulus());
}



inline void
add (bigmod & c, const bigmod & a, int b)
{
	add(c, a, static_cast<long>(b));
}



inline void
subtract (bigmod & c, const bigmod & a, const bigmod & b)
{
	subtract(c, a, b, bigmod::modulus());
}



inline void
subtract (bigmod & c, const bigmod & a, const bigint & b)
{
	subtract(c, a, b, bigmod::modulus());
}



inline void
subtract (bigmod & c, const bigmod & a, long b)
{
	subtract(c, a, b, bigmod::modulus());
}



inline void
subtract (bigmod & c, const bigmod & a, unsigned long b)
{
	subtract(c, a, b, bigmod::modulus());
}



inline void
subtract (bigmod & c, const bigmod & a, int b)
{
	subtract(c, a, static_cast<long>(b));
}



inline void
multiply (bigmod & c, const bigmod & a, const bigmod & b)
{
	multiply(c, a, b, bigmod::modulus());
}



inline void
multiply (bigmod & c, const bigmod & a, const bigint & b)
{
	multiply(c, a, b, bigmod::modulus());
}



inline void
multiply (bigmod & c, const bigmod & a, long b)
{
	multiply(c, a, b, bigmod::modulus());
}



inline void
multiply (bigmod & c, const bigmod & a, unsigned long b)
{
	multiply(c, a, b, bigmod::modulus());
}



inline void
multiply (bigmod & c, const bigmod & a, int b)
{
	multiply(c, a, static_cast<long>(b));
}



inline void
square (bigmod & c, const bigmod & a)
{
	square(c, a, bigmod::modulus());
}



inline void
square (bigmod & c, const bigint & a)
{
	square(c, a, bigmod::modulus());
}



inline void
divide (bigmod & c, const bigmod & a, const bigmod & b)
{
	divide(c, a, b, bigmod::modulus());
}



inline void
divide (bigmod & c, const bigmod & a, const bigint & b)
{
	divide(c, a, b, bigmod::modulus());
}



inline void
divide (bigmod & c, const bigmod & a, long b)
{
	divide(c, a, b, bigmod::modulus());
}



inline void
divide (bigmod & c, const bigmod & a, unsigned long b)
{
	divide(c, a, b, bigmod::modulus());
}



inline void
divide (bigmod & c, const bigmod & a, int b)
{
	divide(c, a, static_cast<long>(b));
}



inline void
power (bigmod & c, const bigmod & a, const bigint & b)
{
	power(c, a, b, bigmod::modulus());
}



inline void
power (bigmod & c, const bigmod & a, long b)
{
	power(c, a, b, bigmod::modulus());
}



inline void
inc (bigmod & c)
{
	c.inc();
}



inline void
dec (bigmod & c)
{
	c.dec();
}



//
// arithmetic operators
//

inline bigmod
operator - (const bigmod & a)
{
	bigmod c(a);

	c.negate();
	return c;
}



inline bigmod
operator + (const bigmod & a, const bigmod & b)
{
	bigmod c;

	add(c, a, b);
	return c;
}



inline bigmod
operator + (const bigmod & a, const bigint & b)
{
	bigmod c;

	add(c, a, b);
	return c;
}



inline bigmod
operator + (const bigmod & a, int b)
{
	bigmod c;

	add(c, a, static_cast<long>(b));
	return c;
}



inline bigmod
operator + (const bigmod & a, long b)
{
	bigmod c;

	add(c, a, b);
	return c;
}



inline bigmod
operator + (const bigmod & a, unsigned long b)
{
	bigmod c;

	add(c, a, b);
	return c;
}



inline bigmod
operator + (const bigint & a, const bigmod & b)
{
	bigmod c;

	add(c, b, a);
	return c;
}



inline bigmod
operator + (int a, const bigmod & b)
{
	bigmod c;

	add(c, b, static_cast<long>(a));
	return c;
}



inline bigmod
operator + (long a, const bigmod & b)
{
	bigmod c;

	add(c, b, a);
	return c;
}



inline bigmod
operator + (unsigned long a, const bigmod & b)
{
	bigmod c;

	add(c, b, a);
	return c;
}



inline bigmod
operator - (const bigmod & a, const bigmod & b)
{
	bigmod c;

	subtract(c, a, b);
	return c;
}



inline bigmod
operator - (const bigmod & a, const bigint & b)
{
	bigmod c;

	subtract(c, a, b);
	return c;
}



inline bigmod
operator - (const bigmod & a, int b)
{
	bigmod c;

	subtract(c, a, static_cast<long>(b));
	return c;
}



inline bigmod
operator - (const bigmod & a, long b)
{
	bigmod c;

	subtract(c, a, b);
	return c;
}



inline bigmod
operator - (const bigmod & a, unsigned long b)
{
	bigmod c;

	subtract(c, a, b);
	return c;
}



inline bigmod
operator - (const bigint & a, const bigmod & b)
{
	bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigmod
operator - (int a, const bigmod & b)
{
	bigmod c;

	subtract(c, b, static_cast<long>(a));
	c.negate();
	return c;
}



inline bigmod
operator - (long a, const bigmod & b)
{
	bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigmod
operator - (unsigned long a, const bigmod & b)
{
	bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline bigmod
operator * (const bigmod & a, const bigmod & b)
{
	bigmod c;

	multiply(c, a, b);
	return c;
}



inline bigmod
operator * (const bigmod & a, const bigint & b)
{
	bigmod c;

	multiply(c, a, b);
	return c;
}



inline bigmod
operator * (const bigmod & a, int b)
{
	bigmod c;

	multiply(c, a, static_cast<long>(b));
	return c;
}



inline bigmod
operator * (const bigmod & a, long b)
{
	bigmod c;

	multiply(c, a, b);
	return c;
}



inline bigmod
operator * (const bigmod & a, unsigned long b)
{
	bigmod c;

	multiply(c, a, b);
	return c;
}



inline bigmod
operator * (const bigint & a, const bigmod & b)
{
	bigmod c;

	multiply(c, b, a);
	return c;
}



inline bigmod
operator * (int a, const bigmod & b)
{
	bigmod c;

	multiply(c, b, static_cast<long>(a));
	return c;
}



inline bigmod
operator * (long a, const bigmod & b)
{
	bigmod c;

	multiply(c, b, a);
	return c;
}



inline bigmod
operator * (unsigned long a, const bigmod & b)
{
	bigmod c;

	multiply(c, b, a);
	return c;
}



inline bigmod
operator / (const bigmod & a, const bigmod & b)
{
	bigmod c;

	divide(c, a, b);
	return c;
}



inline bigmod
operator / (const bigmod & a, const bigint & b)
{
	bigmod c;

	divide(c, a, b);
	return c;
}



inline bigmod
operator / (const bigmod & a, int b)
{
	bigmod c;

	divide(c, a, static_cast<long>(b));
	return c;
}



inline bigmod
operator / (const bigmod & a, long b)
{
	bigmod c;

	divide(c, a, b);
	return c;
}



inline bigmod
operator / (const bigmod & a, unsigned long b)
{
	bigmod c;

	divide(c, a, b);
	return c;
}



inline bigmod
operator / (const bigint & a, const bigmod & b)
{
	bigmod c;

	divide(c, b, a);
	c.invert();
	return c;
}



inline bigmod
operator / (int a, const bigmod & b)
{
	bigmod c;

	divide(c, b, static_cast<long>(a));
	c.invert();
	return c;
}



inline bigmod
operator / (long a, const bigmod & b)
{
	bigmod c;

	divide(c, b, a);
	c.invert();
	return c;
}



inline bigmod
operator / (unsigned long a, const bigmod & b)
{
	bigmod c;

	divide(c, b, a);
	c.invert();
	return c;
}



inline bigmod &
operator += (bigmod & a, const bigmod & b)
{
	add(a, a, b);
	return a;
}



inline bigmod &
operator += (bigmod & a, const bigint & b)
{
	add(a, a, b);
	return a;
}



inline bigmod &
operator += (bigmod & a, int b)
{
	add(a, a, static_cast<long>(b));
	return a;
}



inline bigmod &
operator += (bigmod & a, long b)
{
	add(a, a, b);
	return a;
}



inline bigmod &
operator += (bigmod & a, unsigned long b)
{
	add(a, a, b);
	return a;
}



inline bigmod &
operator -= (bigmod & a, const bigmod & b)
{
	subtract(a, a, b);
	return a;
}



inline bigmod &
operator -= (bigmod & a, const bigint & b)
{
	subtract(a, a, b);
	return a;
}



inline bigmod &
operator -= (bigmod & a, int b)
{
	subtract(a, a, static_cast<long>(b));
	return a;
}



inline bigmod &
operator -= (bigmod & a, long b)
{
	subtract(a, a, b);
	return a;
}



inline bigmod &
operator -= (bigmod & a, unsigned long b)
{
	subtract(a, a, b);
	return a;
}



inline bigmod &
operator *= (bigmod & a, const bigmod & b)
{
	multiply(a, a, b);
	return a;
}



inline bigmod &
operator *= (bigmod & a, const bigint & b)
{
	multiply(a, a, b);
	return a;
}



inline bigmod &
operator *= (bigmod & a, int b)
{
	multiply(a, a, static_cast<long>(b));
	return a;
}



inline bigmod &
operator *= (bigmod & a, long b)
{
	multiply(a, a, b);
	return a;
}



inline bigmod &
operator *= (bigmod & a, unsigned long b)
{
	multiply(a, a, b);
	return a;
}



inline bigmod &
operator /= (bigmod & a, const bigmod & b)
{
	divide(a, a, b);
	return a;
}



inline bigmod &
operator /= (bigmod & a, const bigint & b)
{
	divide(a, a, b);
	return a;
}



inline bigmod &
operator /= (bigmod & a, int b)
{
	divide(a, a, static_cast<long>(b));
	return a;
}



inline bigmod &
operator /= (bigmod & a, long b)
{
	divide(a, a, b);
	return a;
}



inline bigmod &
operator /= (bigmod & a, unsigned long b)
{
	divide(a, a, b);
	return a;
}



inline bigmod &
operator ++ (bigmod & a)
{
	a.inc();
	return a;
}



inline bigmod
operator ++ (bigmod & a, int)
{
	bigmod c(a);

	a.inc();
	return c;
}



inline bigmod &
operator -- (bigmod & a)
{
	a.dec();
	return a;
}



inline bigmod
operator -- (bigmod & a, int)
{
	bigmod c(a);

	a.dec();
	return c;
}



inline bool
operator ! (const bigmod & a)
{
	return a.is_zero();
}



//
// functions
//

inline bigmod
inverse (const bigmod & a)
{
	bigmod c(a);

	c.invert();
	return c;
}



inline bigmod
randomize (const bigmod & a)
{
	bigmod c(randomize(a.mantissa()));

	return c;
}



//
// misc. procedures
//

inline void
swap (bigmod & a, bigmod & b)
{
	a.swap(b);
}



inline const bigint &
bigmod::characteristic ()
{
	return bigmod::M;
}



//
// I/O
//

std::istream & operator >> (std::istream & in, bigmod & a);
std::ostream & operator << (std::ostream & out, const bigmod & a);

int string_to_bigmod (char * s, bigmod & a);
int bigmod_to_string (const bigmod & a, char * s);



//
// bigmod_lib functions
//

//**** relative degree of element over prime field *************

inline unsigned int
get_relative_degree(const bigmod &)
{
	return 1;
}


//**** absolute degree over prime field ************************

inline unsigned int
get_absolute_degree(const bigmod &)
{
	return 1;
}


//**** get generator of finite field ***************************

bigmod get_generator(const bigmod & x);


//**** hash functions ******************************************

udigit hash(const bigmod & x);


//**** is_square ***********************************************

bool is_square(const bigmod & a);


//**** sqrt ****************************************************

bigmod sqrt(const bigmod & a);


//**** solve_quadratic *****************************************

bool solve_quadratic(bigmod & root, const bigmod & a1,
		     const bigmod & a0);


//**** characteristic ******************************************

bigint characteristic(const bigmod & a);


//**** number_of_elements **************************************

bigint number_of_elements(const bigmod & a);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_BIGMOD


#include	"LiDIA/specialization/bigmod.special"



#endif	// LIDIA_BIGMOD_H_GUARD_
