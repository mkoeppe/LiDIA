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


#ifndef LIDIA_MULTI_BIGMOD_H_GUARD_
#define LIDIA_MULTI_BIGMOD_H_GUARD_



#ifndef LIDIA_BASE_BIGMOD_H_GUARD_
# include	"LiDIA/base/base_bigmod.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class multi_bigmod : public base_bigmod
{
private :

#ifndef HEADBANGER

        static residue_class< bigint > *zero;
        residue_class< bigint > *Mp;



        void assign_modulus (residue_class< bigint > * a);

        void initialize();
	void read  (std::istream & in);
        void print (std::ostream & out) const;

#endif	// HEADBANGER

public :

        //
	// c'tors and d'tor
	//

	// (remark: There exists no constructor multi_bigmod(bigmod)
	// to avoid implicit casts (int, long, bigint, ...) ->bigmod ->multi_bigmod.
	// These implicit casts are traps for the programmer, who
	// used the wrong constructor, i.e. multi_bigmod(5), which sets the
        // modulus to bigmod::modulus(), instead of multi_bigmod(5, m),
	// which sets the modulus to m.)

	multi_bigmod ();
	multi_bigmod (int i, const bigint & m);
	multi_bigmod (long l, const bigint & m);
	multi_bigmod (unsigned long ul, const bigint & m);
	multi_bigmod (double d, const bigint & m);
	multi_bigmod (const bigint & i, const bigint & m);
	multi_bigmod (const multi_bigmod & a);
        virtual ~multi_bigmod ();

#ifndef HEADBANGER



        //
	// modulus modifiers and accessors
	//

        void set_modulus (const bigint & m);
        void set_modulus (const multi_bigmod & a);
	const bigint & modulus () const;



        //
        // assigners
        //

        void assign_zero ();
        void assign_one ();

        void assign_zero (const bigint & m);
        void assign_one (const bigint & m);
        void assign (int i, const bigint & m);
        void assign (long i, const bigint & m);
        void assign (unsigned long i, const bigint & m);
        void assign (const bigint & a, const bigint & m);
        void assign (const multi_bigmod & a);
	void assign (const bigmod & a);

        void set_mantissa (int i);
        void set_mantissa (long i);
        void set_mantissa (unsigned long i);
        void set_mantissa (const bigint & a);

	multi_bigmod & operator = (int i);
	multi_bigmod & operator = (long l);
	multi_bigmod & operator = (unsigned long ul);
	multi_bigmod & operator = (double d);
	multi_bigmod & operator = (const bigint & a);
	multi_bigmod & operator = (const multi_bigmod  & a);
	multi_bigmod & operator = (const bigmod & a);



        //
	// comparators
	//

	bool is_equal (const multi_bigmod & a) const;
	bool is_equal (const bigmod & a) const;
	bool is_equal (const bigint & a) const;
	bool is_equal (long a) const;
	bool is_equal (unsigned long a) const;



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

        void randomize (const bigint & m);
        void randomize (const multi_bigmod & m);

	void swap (multi_bigmod & a);



        //
        // arithmetic procedures
	//

	friend void add (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
	friend void add (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
	friend void add (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
	friend void add (multi_bigmod & c, const multi_bigmod & a, long b);
	friend void add (multi_bigmod & c, const multi_bigmod & a, unsigned long b);

	friend void subtract (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
	friend void subtract (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
	friend void subtract (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
	friend void subtract (multi_bigmod & c, const multi_bigmod & a, long b);
	friend void subtract (multi_bigmod & c, const multi_bigmod & a, unsigned long b);

	friend void multiply (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
	friend void multiply (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
	friend void multiply (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
	friend void multiply (multi_bigmod & c, const multi_bigmod & a, long b);
	friend void multiply (multi_bigmod & c, const multi_bigmod & a, unsigned long b);
	friend void square (multi_bigmod & a, const multi_bigmod & b);

	friend void divide (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
	friend void divide (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
	friend void divide (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
	friend void divide (multi_bigmod & c, const multi_bigmod & a, long b);
	friend void divide (multi_bigmod & c, const multi_bigmod & a, unsigned long b);

        friend void power (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
        friend void power (multi_bigmod & c, const multi_bigmod & a, long b);

#endif



	//
	// I/O
	//

	friend int string_to_multi_bigmod (char * s, multi_bigmod & a);
	friend int multi_bigmod_to_string (const multi_bigmod & a, char * &s);

	friend std::istream & operator >> (std::istream & in, multi_bigmod & a);
	friend std::ostream & operator << (std::ostream & out, const multi_bigmod & a);

#ifdef C_STDIO

	void read_from_file (FILE * fp);
	void write_to_file  (FILE * fp);
	void scan_from_file (FILE * fp);
	void print_to_file  (FILE * fp);

#endif	// C_STDIO

};


//
// arithmetic procedures
//

void add (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
void add (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
void add (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
void add (multi_bigmod & c, const multi_bigmod & a, long b);
void add (multi_bigmod & c, const multi_bigmod & a, unsigned long b);

void subtract (multi_bigmod & c, const multi_bigmod & a,
	       const multi_bigmod & b);
void subtract (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
void subtract (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
void subtract (multi_bigmod & c, const multi_bigmod & a, long b);
void subtract (multi_bigmod & c, const multi_bigmod & a, unsigned long b);

void multiply (multi_bigmod & c, const multi_bigmod & a,
	       const multi_bigmod & b);
void multiply (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
void multiply (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
void multiply (multi_bigmod & c, const multi_bigmod & a, long b);
void multiply (multi_bigmod & c, const multi_bigmod & a, unsigned long b);
void square (multi_bigmod & a, const multi_bigmod & b);

void divide (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
void divide (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);
void divide (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
void divide (multi_bigmod & c, const multi_bigmod & a, long b);
void divide (multi_bigmod & c, const multi_bigmod & a, unsigned long b);

void power (multi_bigmod & c, const multi_bigmod & a, const bigint & b);
void power (multi_bigmod & c, const multi_bigmod & a, long b);

//
// I/O
//

int string_to_multi_bigmod (char * s, multi_bigmod & a);
int multi_bigmod_to_string (const multi_bigmod & a, char * &s);

std::istream & operator >> (std::istream & in, multi_bigmod & a);
std::ostream & operator << (std::ostream & out, const multi_bigmod & a);


inline void
multi_bigmod::assign_modulus (residue_class< bigint > * a)
{
	if (Mp != a) {
		base_bigmod::L->clear(Mp);
		Mp = base_bigmod::L->set_to(a);
	}
}



//
// normalizer (inline definition must precede c'tors and assigners)
//

inline void
multi_bigmod::normalize ()
{
	base_bigmod::normalize(Mp->get_mod());
}



//
// c'tors and d'tor
//

inline
multi_bigmod::multi_bigmod ()
	: base_bigmod()
{
	initialize();
	Mp = base_bigmod::L->set_to(zero);
}



inline
multi_bigmod::multi_bigmod (int i, const bigint & m)
	: base_bigmod(i),
	  Mp(NULL)
{
	initialize();
	set_modulus(m);
	normalize();
}



inline
multi_bigmod::multi_bigmod (long l, const bigint & m)
	: base_bigmod(l),
	  Mp(NULL)
{
	initialize();
	set_modulus(m);
	normalize();
}



inline
multi_bigmod::multi_bigmod (unsigned long ul, const bigint & m)
	: base_bigmod(ul),
	  Mp(NULL)
{
	initialize();
	set_modulus(m);
	normalize();
}



inline
multi_bigmod::multi_bigmod (double d, const bigint & m)
	: base_bigmod(d),
	  Mp(NULL)
{
	initialize();
	set_modulus(m);
	normalize();
}



inline
multi_bigmod::multi_bigmod (const bigint & i, const bigint & m)
	: base_bigmod(i),
	  Mp(NULL)
{
	initialize();
	set_modulus(m);
	normalize();
}



inline
multi_bigmod::multi_bigmod (const multi_bigmod & a)
	: base_bigmod(a.I),
	  Mp(base_bigmod::L->set_to(a.Mp))
{
}



inline
multi_bigmod::~multi_bigmod ()
{
	base_bigmod::L->clear(Mp);
}



//
// accessors
//

inline const bigint &
multi_bigmod::modulus () const
{
	return Mp->get_mod();
}



inline void
multi_bigmod::set_modulus (const multi_bigmod & a)
{
	this->assign_modulus(a.Mp);
}



//
// assigners
//

inline void
multi_bigmod::assign_zero ()
{
	I.assign_zero();
}



inline void
multi_bigmod::assign_one ()
{
	I.assign_one();
}



inline void
multi_bigmod::assign_zero (const bigint & m)
{
	I.assign_zero();
	set_modulus(m);
}



inline void
multi_bigmod::assign_one (const bigint & m)
{
	I.assign_one();
	set_modulus(m);
}



inline void
multi_bigmod::assign (int i, const bigint & m)
{
	I.assign (i);
	set_modulus(m);
	normalize();
}



inline void
multi_bigmod::assign (long i, const bigint & m)
{
	I.assign (i);
	set_modulus(m);
	normalize();
}



inline void
multi_bigmod::assign (unsigned long i, const bigint & m)
{
	I.assign (i);
	set_modulus(m);
	normalize();
}



inline void
multi_bigmod::assign (const bigint & a, const bigint & m)
{
	I.assign (a);
	set_modulus(m);
	normalize();
}



inline void
multi_bigmod::assign (const multi_bigmod & a)
{
	if (&a != this) {
		I.assign(a.I);
		assign_modulus(a.Mp);
	}
}



inline void
multi_bigmod::assign (const bigmod & a)
{
	I.assign(a.mantissa());
	base_bigmod::L->clear(Mp);
	Mp = base_bigmod::L->set_to(bigmod::residue_class());
}



inline void
multi_bigmod::set_mantissa (int i)
{
	I.assign (i);
	normalize();
}



inline void
multi_bigmod::set_mantissa (long i)
{
	I.assign (i);
	normalize();
}



inline void
multi_bigmod::set_mantissa (unsigned long i)
{
	I.assign (i);
	normalize();
}



inline void
multi_bigmod::set_mantissa (const bigint & a)
{
	I.assign (a);
	normalize();
}



inline multi_bigmod &
multi_bigmod::operator = (int i)
{
	I.assign(i);
	return *this;
}



inline multi_bigmod &
multi_bigmod::operator = (long l)
{
	I.assign(l);
	return *this;
}



inline multi_bigmod &
multi_bigmod::operator = (unsigned long ul)
{
	I.assign(ul);
	return *this;
}



inline multi_bigmod &
multi_bigmod::operator = (double d)
{
	I.assign(d);
	return *this;
}



inline multi_bigmod &
multi_bigmod::operator = (const bigint & a)
{
	I.assign(a);
	return *this;
}



inline multi_bigmod &
multi_bigmod::operator = (const multi_bigmod & a)
{
	assign(a);
	return *this;
}



inline multi_bigmod &
multi_bigmod::operator = (const bigmod & a)
{
	assign(a);
	return *this;
}



//
// comparators
//

inline bool
multi_bigmod::is_equal (const multi_bigmod & a) const
{
	return (this == &a || (I.compare(a.I) == 0 && Mp == a.Mp));
}



inline bool
multi_bigmod::is_equal (const bigmod & a) const
{
	return (I.compare(a.mantissa()) == 0 && Mp == a.residue_class());
}



inline bool
multi_bigmod::is_equal (const bigint & a) const
{
	return (a.is_ge_zero() ? (I.compare(a) == 0) : I.compare(modulus() + a));
}



inline bool
multi_bigmod::is_equal (long a) const
{
	return ((a >= 0) ? (I.compare(a) == 0) : I.compare(modulus() + a));
}



inline bool
multi_bigmod::is_equal (unsigned long a) const
{
	return (I.compare(a) == 0);
}



inline bool
operator == (const multi_bigmod & a, const multi_bigmod & b)
{
	return a.is_equal(b);
}



inline bool
operator == (const multi_bigmod & a, const bigmod & b)
{
	return a.is_equal(b);
}



inline bool
operator == (const multi_bigmod & a, const bigint & b)
{
	return a.is_equal(b);
}



inline bool
operator == (const multi_bigmod & a, long b)
{
	return a.is_equal(b);
}



inline bool
operator == (const multi_bigmod & a, unsigned long b)
{
	return a.is_equal(b);
}



inline bool
operator == (const multi_bigmod & a, int b)
{
	return a.is_equal(static_cast<long>(b));
}



inline bool
operator == (const bigmod & a, const multi_bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (const bigint & a, const multi_bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (long a, const multi_bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (unsigned long a, const multi_bigmod & b)
{
	return b.is_equal(a);
}



inline bool
operator == (int a, const multi_bigmod & b)
{
	return b.is_equal(static_cast<long>(a));
}



inline bool
operator != (const multi_bigmod & a, const multi_bigmod & b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const multi_bigmod & a, const bigmod & b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const multi_bigmod & a, const bigint & b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const multi_bigmod & a, long b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const multi_bigmod & a, unsigned long b)
{
	return !a.is_equal(b);
}



inline bool
operator != (const multi_bigmod & a, int b)
{
	return !a.is_equal(static_cast<long>(b));
}



inline bool
operator != (const bigmod & a, const multi_bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (const bigint & a, const multi_bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (long a, const multi_bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (unsigned long a, const multi_bigmod & b)
{
	return !b.is_equal(a);
}



inline bool
operator != (int a, const multi_bigmod & b)
{
	return !b.is_equal(static_cast<long>(a));
}



//
// modifiers
//

inline bigint
multi_bigmod::invert (int verbose)
{
	return base_bigmod::invert(Mp->get_mod(), verbose);
}



inline void
multi_bigmod::negate ()
{
	base_bigmod::negate(Mp->get_mod());
}



inline void
multi_bigmod::inc ()
{
	base_bigmod::inc(Mp->get_mod());
}



inline void
multi_bigmod::dec ()
{
	base_bigmod::dec(Mp->get_mod());
}



inline void
multi_bigmod::multiply_by_2 ()
{
	base_bigmod::multiply_by_2(Mp->get_mod());
}



inline void
multi_bigmod::divide_by_2 ()
{
	base_bigmod::divide_by_2(Mp->get_mod());
}



inline void
multi_bigmod::swap (multi_bigmod & a)
{
	I.swap(a.I);

	residue_class< bigint > * h;
	h = Mp;
	Mp = a.Mp;
	a.Mp = h;
}



//
// arithmetic procedures
//

inline void
negate (multi_bigmod & a, const multi_bigmod & b)
{
	a.assign(b);
	a.negate();
}



void add (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
void add (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);



inline void
add (multi_bigmod & c, const multi_bigmod & a, const bigint & b)
{
	add (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
add (multi_bigmod & c, const multi_bigmod & a, long b)
{
	add (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
add (multi_bigmod & c, const multi_bigmod & a, unsigned long b)
{
	add (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
add (multi_bigmod & c, const multi_bigmod & a, int b)
{
	add(c, a, static_cast<long>(b));
}



void subtract (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
void subtract (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);



inline void
subtract (multi_bigmod & c, const multi_bigmod & a, const bigint & b)
{
	subtract (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
subtract (multi_bigmod & c, const multi_bigmod & a, long b)
{
	subtract (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
subtract (multi_bigmod & c, const multi_bigmod & a, unsigned long b)
{
	subtract (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
subtract (multi_bigmod & c, const multi_bigmod & a, int b)
{
	subtract(c, a, static_cast<long>(b));
}



void multiply (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
void multiply (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);



inline void
multiply (multi_bigmod & c, const multi_bigmod & a, const bigint & b)
{
	multiply (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
multiply (multi_bigmod & c, const multi_bigmod & a, long b)
{
	multiply (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
multiply (multi_bigmod & c, const multi_bigmod & a, unsigned long b)
{
	multiply (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
multiply (multi_bigmod & c, const multi_bigmod & a, int b)
{
	multiply(c, a, static_cast<long>(b));
}



void square (multi_bigmod & a, const multi_bigmod & b);



inline void
invert (multi_bigmod & a, const multi_bigmod & b)
{
	a.assign(b);
	a.invert();
}



void divide (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b);
void divide (multi_bigmod & c, const multi_bigmod & a, const bigmod & b);



inline void
divide (multi_bigmod & c, const multi_bigmod & a, const bigint & b)
{
	divide (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
divide (multi_bigmod & c, const multi_bigmod & a, long b)
{
	divide (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
divide (multi_bigmod & c, const multi_bigmod & a, unsigned long b)
{
	divide (c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
divide (multi_bigmod & c, const multi_bigmod & a, int b)
{
	divide(c, a, static_cast<long>(b));
}



inline void
power (multi_bigmod & c, const multi_bigmod & a, const bigint & b)
{
	power(c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
power (multi_bigmod & c, const multi_bigmod & a, long b)
{
	power(c, a, b, a.modulus());
	c.assign_modulus (a.Mp);
}



inline void
inc (multi_bigmod & c)
{
	c.inc();
}



inline void
dec (multi_bigmod & c)
{
	c.dec();
}



//
// functions
//

inline multi_bigmod
inverse (const multi_bigmod & a)
{
	multi_bigmod c;

	invert(c, a);
	return c;
}



inline multi_bigmod
randomize (const multi_bigmod & a)
{
	multi_bigmod c;

	c.randomize(a);
	return c;
}



//
// arithmetic operators
//

inline multi_bigmod
operator - (const multi_bigmod & a)
{
	multi_bigmod c(a);

	c.negate();
	return c;
}



inline multi_bigmod
operator + (const multi_bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	add(c, a, b);
	return c;
}



inline multi_bigmod
operator + (const multi_bigmod & a, const bigmod & b)
{
	multi_bigmod c;

	add(c, a, b);
	return c;
}



inline multi_bigmod
operator + (const multi_bigmod & a, const bigint & b)
{
	multi_bigmod c;

	add(c, a, b);
	return c;
}



inline multi_bigmod
operator + (const multi_bigmod & a, int b)
{
	multi_bigmod c;

	add(c, a, static_cast<long>(b));
	return c;
}



inline multi_bigmod
operator + (const multi_bigmod & a, long b)
{
	multi_bigmod c;

	add(c, a, b);
	return c;
}



inline multi_bigmod
operator + (const multi_bigmod & a, unsigned long b)
{
	multi_bigmod c;

	add(c, a, b);
	return c;
}



inline multi_bigmod
operator + (const bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	add(c, b, a);
	return c;
}



inline multi_bigmod
operator + (const bigint & a, const multi_bigmod & b)
{
	multi_bigmod c;

	add(c, b, a);
	return c;
}



inline multi_bigmod
operator + (int a, const multi_bigmod & b)
{
	multi_bigmod c;

	add(c, b, static_cast<long>(a));
	return c;
}



inline multi_bigmod
operator + (long a, const multi_bigmod & b)
{
	multi_bigmod c;

	add(c, b, a);
	return c;
}



inline multi_bigmod
operator + (unsigned long a, const multi_bigmod & b)
{
	multi_bigmod c;

	add(c, b, a);
	return c;
}



inline multi_bigmod
operator - (const multi_bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	subtract(c, a, b);
	return c;
}



inline multi_bigmod
operator - (const multi_bigmod & a, const bigmod & b)
{
	multi_bigmod c;

	subtract(c, a, b);
	return c;
}



inline multi_bigmod
operator - (const multi_bigmod & a, const bigint & b)
{
	multi_bigmod c;

	subtract(c, a, b);
	return c;
}



inline multi_bigmod
operator - (const multi_bigmod & a, int b)
{
	multi_bigmod c;

	subtract(c, a, static_cast<long>(b));
	return c;
}



inline multi_bigmod
operator - (const multi_bigmod & a, long b)
{
	multi_bigmod c;

	subtract(c, a, b);
	return c;
}



inline multi_bigmod
operator - (const multi_bigmod & a, unsigned long b)
{
	multi_bigmod c;

	subtract(c, a, b);
	return c;
}



inline multi_bigmod
operator - (const bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline multi_bigmod
operator - (const bigint & a, const multi_bigmod & b)
{
	multi_bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline multi_bigmod
operator - (int a, const multi_bigmod & b)
{
	multi_bigmod c;

	subtract(c, b, static_cast<long>(a));
	c.negate();
	return c;
}



inline multi_bigmod
operator - (long a, const multi_bigmod & b)
{
	multi_bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline multi_bigmod
operator - (unsigned long a, const multi_bigmod & b)
{
	multi_bigmod c;

	subtract(c, b, a);
	c.negate();
	return c;
}



inline multi_bigmod
operator * (const multi_bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	multiply(c, a, b);
	return c;
}



inline multi_bigmod
operator * (const multi_bigmod & a, const bigmod & b)
{
	multi_bigmod c;

	multiply(c, a, b);
	return c;
}



inline multi_bigmod
operator * (const multi_bigmod & a, const bigint & b)
{
	multi_bigmod c;

	multiply(c, a, b);
	return c;
}



inline multi_bigmod
operator * (const multi_bigmod & a, int b)
{
	multi_bigmod c;

	multiply(c, a, static_cast<long>(b));
	return c;
}



inline multi_bigmod
operator * (const multi_bigmod & a, long b)
{
	multi_bigmod c;

	multiply(c, a, b);
	return c;
}



inline multi_bigmod
operator * (const multi_bigmod & a, unsigned long b)
{
	multi_bigmod c;

	multiply(c, a, b);
	return c;
}



inline multi_bigmod
operator * (const bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	multiply(c, b, a);
	return c;
}



inline multi_bigmod
operator * (const bigint & a, const multi_bigmod & b)
{
	multi_bigmod c;

	multiply(c, b, a);
	return c;
}



inline multi_bigmod
operator * (int a, const multi_bigmod & b)
{
	multi_bigmod c;

	multiply(c, b, static_cast<long>(a));
	return c;
}



inline multi_bigmod
operator * (long a, const multi_bigmod & b)
{
	multi_bigmod c;

	multiply(c, b, a);
	return c;
}



inline multi_bigmod
operator * (unsigned long a, const multi_bigmod & b)
{
	multi_bigmod c;

	multiply(c, b, a);
	return c;
}



inline multi_bigmod
operator / (const multi_bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c;

	divide(c, a, b);
	return c;
}



inline multi_bigmod
operator / (const multi_bigmod & a, const bigmod & b)
{
	multi_bigmod c;

	divide(c, a, b);
	return c;
}



inline multi_bigmod
operator / (const multi_bigmod & a, const bigint & b)
{
	multi_bigmod c;

	divide(c, a, b);
	return c;
}



inline multi_bigmod
operator / (const multi_bigmod & a, int b)
{
	multi_bigmod c;

	divide(c, a, static_cast<long>(b));
	return c;
}



inline multi_bigmod
operator / (const multi_bigmod & a, long b)
{
	multi_bigmod c;

	divide(c, a, b);
	return c;
}



inline multi_bigmod
operator / (const multi_bigmod & a, unsigned long b)
{
	multi_bigmod c;

	divide(c, a, b);
	return c;
}



inline multi_bigmod
operator / (const bigmod & a, const multi_bigmod & b)
{
	multi_bigmod c(b);

	c.invert();
	multiply(c, c, a);
	return c;
}



inline multi_bigmod
operator / (const bigint & a, const multi_bigmod & b)
{
	multi_bigmod c(b);

	c.invert();
	multiply(c, c, a);
	return c;
}



inline multi_bigmod
operator / (int a, const multi_bigmod & b)
{
	multi_bigmod c(b);

	c.invert();
	multiply(c, c, static_cast<long>(a));
	return c;
}



inline multi_bigmod
operator / (long a, const multi_bigmod & b)
{
	multi_bigmod c(b);

	c.invert();
	multiply(c, c, a);
	return c;
}



inline multi_bigmod
operator / (unsigned long a, const multi_bigmod & b)
{
	multi_bigmod c(b);

	c.invert();
	multiply(c, c, a);
	return c;
}



inline multi_bigmod &
operator += (multi_bigmod & a, const multi_bigmod & b)
{
	add(a, a, b);
	return a;
}



inline multi_bigmod &
operator += (multi_bigmod & a, const bigmod & b)
{
	add(a, a, b);
	return a;
}



inline multi_bigmod &
operator += (multi_bigmod & a, const bigint & b)
{
	add(a, a, b);
	return a;
}



inline multi_bigmod &
operator += (multi_bigmod & a, int b)
{
	add(a, a, static_cast<long>(b));
	return a;
}



inline multi_bigmod &
operator += (multi_bigmod & a, long b)
{
	add(a, a, b);
	return a;
}



inline multi_bigmod &
operator += (multi_bigmod & a, unsigned long b)
{
	add(a, a, b);
	return a;
}



inline multi_bigmod &
operator -= (multi_bigmod & a, const multi_bigmod & b)
{
	subtract(a, a, b);
	return a;
}



inline multi_bigmod &
operator -= (multi_bigmod & a, const bigmod & b)
{
	subtract(a, a, b);
	return a;
}



inline multi_bigmod &
operator -= (multi_bigmod & a, const bigint & b)
{
	subtract(a, a, b);
	return a;
}



inline multi_bigmod &
operator -= (multi_bigmod & a, int b)
{
	subtract(a, a, static_cast<long>(b));
	return a;
}



inline multi_bigmod &
operator -= (multi_bigmod & a, long b)
{
	subtract(a, a, b);
	return a;
}



inline multi_bigmod &
operator -= (multi_bigmod & a, unsigned long b)
{
	subtract(a, a, b);
	return a;
}



inline multi_bigmod &
operator *= (multi_bigmod & a, const multi_bigmod & b)
{
	multiply(a, a, b);
	return a;
}



inline multi_bigmod &
operator *= (multi_bigmod & a, const bigmod & b)
{
	multiply(a, a, b);
	return a;
}



inline multi_bigmod &
operator *= (multi_bigmod & a, const bigint & b)
{
	multiply(a, a, b);
	return a;
}



inline multi_bigmod &
operator *= (multi_bigmod & a, int b)
{
	multiply(a, a, static_cast<long>(b));
	return a;
}


inline multi_bigmod &
operator *= (multi_bigmod & a, long b)
{
	multiply(a, a, b);
	return a;
}


inline multi_bigmod &
operator *= (multi_bigmod & a, unsigned long b)
{
	multiply(a, a, b);
	return a;
}


inline multi_bigmod &
operator /= (multi_bigmod & a, const multi_bigmod & b)
{
	divide(a, a, b);
	return a;
}



inline multi_bigmod &
operator /= (multi_bigmod & a, const bigmod & b)
{
	divide(a, a, b);
	return a;
}



inline multi_bigmod &
operator /= (multi_bigmod & a, const bigint & b)
{
	divide(a, a, b);
	return a;
}



inline multi_bigmod &
operator /= (multi_bigmod & a, int b)
{
	divide(a, a, static_cast<long>(b));
	return a;
}



inline multi_bigmod &
operator /= (multi_bigmod & a, long b)
{
	divide(a, a, b);
	return a;
}



inline multi_bigmod &
operator /= (multi_bigmod & a, unsigned long b)
{
	divide(a, a, b);
	return a;
}



inline multi_bigmod &
operator ++ (multi_bigmod & a)
{
	a.inc();
	return a;
}



inline multi_bigmod
operator ++ (multi_bigmod & a, int)
{
	multi_bigmod c(a);

	a.inc();
	return c;
}



inline multi_bigmod &
operator -- (multi_bigmod & a)
{
	a.dec();
	return a;
}



inline multi_bigmod
operator -- (multi_bigmod & a, int)
{
	multi_bigmod c(a);

	a.dec();
	return c;
}



inline bool
operator ! (const multi_bigmod & a)
{
	return a.is_zero();
}



//
// procedural modifiers
//

inline void
swap (multi_bigmod & a, multi_bigmod & b)
{
	a.swap(b);
}



//
// I/O
//

std::istream & operator >> (std::istream & in, multi_bigmod & a);
std::ostream & operator << (std::ostream & out, const multi_bigmod & a);



int string_to_multi_bigmod (char * s, multi_bigmod & a);
int multi_bigmod_to_string (const multi_bigmod & a, char * &s);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MULTI_BIGMOD_H_GUARD_
