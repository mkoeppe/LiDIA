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


#ifndef LIDIA_BASE_BIGMOD_H_GUARD_
#define LIDIA_BASE_BIGMOD_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_RESIDUE_CLASS_LIST_H_GUARD_
# ifndef HEADBANGER
#  include	"LiDIA/base/residue_class_list.h"
# endif
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class base_bigmod {
protected:

	bigint I;
#ifndef HEADBANGER
	static residue_class_list< bigint > * L;
#endif


	void normalize (const bigint & m);


	//
	// c'tors and dtor
	//

	base_bigmod();
	base_bigmod(int);
	base_bigmod(long);
	base_bigmod(unsigned long);
	base_bigmod(double);
	base_bigmod(const bigint &);
	base_bigmod(const base_bigmod &);
	virtual ~base_bigmod ();



	//
	// assigners
	//

	void assign_zero();
	void assign_one();
	void assign (const base_bigmod & a);

	base_bigmod & operator = (const base_bigmod & a);



	//
	// modifiers
	//

	bigint invert (const bigint & m, int verbose = 0);
	void negate (const bigint & m);

	void inc(const bigint & m);
	void dec(const bigint & m);

	void multiply_by_2(const bigint & m);
	void divide_by_2(const bigint & m);

	void swap(base_bigmod & a);



	//
	// arithmetical procedures
	//

	friend void negate (base_bigmod & c, const base_bigmod & a, const bigint & m);

	friend void add (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m);
	friend void add (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);
	friend void add (base_bigmod & c, const base_bigmod & a, long b, const bigint & m);
	friend void add (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m);
	//friend void add (base_bigmod & c, long i, const base_bigmod & a, const bigint & m);
	//friend void add (base_bigmod & c, const bigint & i, const base_bigmod & a, const bigint & m);

	friend void subtract (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m);
	friend void subtract (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);
	friend void subtract (base_bigmod & c, const base_bigmod & a, long b, const bigint & m);
	friend void subtract (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m);
	//friend void subtract (base_bigmod & c, long i, const base_bigmod & a, const bigint & m);
	//friend void subtract (base_bigmod & c, const bigint & i, const base_bigmod & a, const bigint & m);

	friend void multiply (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m);
	friend void multiply (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);
	friend void multiply (base_bigmod & c, const base_bigmod & a, long b, const bigint & m);
	friend void multiply (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m);
	//friend void multiply (base_bigmod & c, long i, const base_bigmod & a, const bigint & m);
	//friend void multiply (base_bigmod & c, const bigint & i, const base_bigmod & a, const bigint & m);

	friend void square(base_bigmod & a, const base_bigmod & b, const bigint & m);
	//friend void square(base_bigmod & a, int b, const bigint & m);
	//friend void square(base_bigmod & a, long b, const bigint & m);
	//friend void square(base_bigmod & a, unsigned long b, const bigint & m);
	friend void square(base_bigmod & a, const bigint & b, const bigint & m);

	friend void invert(base_bigmod & a, const base_bigmod & b, const bigint & m);
	friend void divide (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m);
	friend void divide (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);
	friend void divide (base_bigmod & c, const base_bigmod & a, long b, const bigint & m);
	friend void divide (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m);
	//friend void divide (base_bigmod & c, long i, const base_bigmod & a, const bigint & m);
	//friend void divide (base_bigmod & c, const bigint & i, const base_bigmod & a, const bigint & m);

	friend void power (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);
	friend void power (base_bigmod & c, const base_bigmod & a, long b, const bigint & m);


	friend base_bigmod inverse(const base_bigmod & a, const bigint & m);


	//
	// ******************************************
	// ***** THE PUBLIC PART OF BASE_BIGMOD *****
	// ******************************************
	//

public :


	//
	// accessors
	//

	const bigint & mantissa () const;



	//
	// comparators
	//

	bool is_zero() const;
	bool is_one() const;



	//
	// properties
	//

	lidia_size_t length() const;
	lidia_size_t bit_length() const;



	//
	// predicates
	//

	bool is_char () const;
	bool is_uchar () const;
	bool is_short () const;
	bool is_ushort () const;
	bool is_int () const;
	bool is_uint () const;
	bool is_long () const;
	bool is_ulong () const;



	//
	// converters
	//

	bool intify(int & i) const;
	bool longify(long & i) const;

};



//
// c'tors and d'tor
//

inline
base_bigmod::base_bigmod ()
{
	debug_handler ("base_bigmod", "base_bigmod()");
}



inline
base_bigmod::base_bigmod (int i)
	: I(i)
{
	debug_handler ("base_bigmod", "base_bigmod(int)");
}



inline
base_bigmod::base_bigmod (long l)
	: I(l)
{
	debug_handler ("base_bigmod", "base_bigmod(long)");
}



inline
base_bigmod::base_bigmod (unsigned long ul)
	: I(ul)
{
	debug_handler ("base_bigmod", "base_bigmod(unsigned long)");
}



inline
base_bigmod::base_bigmod (double d)
	: I(d)
{
	debug_handler ("base_bigmod", "base_bigmod(double)");
}



inline
base_bigmod::base_bigmod (const bigint & bi)
	: I(bi)
{
	debug_handler ("base_bigmod", "base_bigmod(const bigint &)");
}



inline
base_bigmod::base_bigmod (const base_bigmod & a)
	: I(a.I)
{
	debug_handler ("base_bigmod", "base_bigmod(const base_bigmod &)");
}



inline
base_bigmod::~base_bigmod ()
{
	debug_handler ("base_bigmod", "~base_bigmod()");
}



//
// accessors
//

inline const bigint &
base_bigmod::mantissa () const
{
	return I;
}



//
// assigners
//

inline void
base_bigmod::assign_zero()
{
	I.assign_zero();
}



inline void
base_bigmod::assign_one()
{
	I.assign_one();
}



inline void
base_bigmod::assign (const base_bigmod & a)
{
	I.assign(a.I);
}



inline base_bigmod &
base_bigmod::operator = (const base_bigmod & a)
{
	assign(a);
	return *this;
}



//
// predicates
//

inline bool
base_bigmod::is_char () const
{
	return I.is_char();
}



inline bool
base_bigmod::is_uchar () const
{
	return I.is_uchar();
}



inline bool
base_bigmod::is_short () const
{
	return I.is_short();
}



inline bool
base_bigmod::is_ushort () const
{
	return I.is_ushort();
}



inline bool
base_bigmod::is_int () const
{
	return I.is_int();
}



inline bool
base_bigmod::is_uint () const
{
	return I.is_uint();
}



inline bool
base_bigmod::is_long () const
{
	return I.is_long();
}



inline bool
base_bigmod::is_ulong () const
{
	return I.is_ulong();
}



//
// properties
//

inline lidia_size_t
base_bigmod::length() const
{
	return I.length();
}



inline lidia_size_t
base_bigmod::bit_length() const
{
	return I.bit_length();
}



//
// comparators
//

inline bool
base_bigmod::is_zero() const
{
	return I.is_zero();
}



inline bool
base_bigmod::is_one() const
{
	return I.is_one();
}



//
// converters
//

inline bool
base_bigmod::intify(int & i) const
{
	return I.intify(i);
}



inline bool
base_bigmod::longify(long & i) const
{
	return I.longify(i);
}



//
// modifiers
//

inline void
base_bigmod::negate (const bigint & m)
{
	if (!I.is_zero()) {
		I.negate();
		add(I, I, m);
	}
}



inline void
base_bigmod::inc (const bigint & m)
{
	I.inc();
	if (I.abs_compare(m) == 0)
		I.assign_zero();
}



inline void
base_bigmod::dec (const bigint & m)
{
	if (I.is_zero())
		I.assign(m);
	I.dec();
}



//
// accessors
//

inline bigint
mantissa (const base_bigmod & a)
{
	return (a.mantissa());
}



inline double
dbl (const base_bigmod & a)
{
	return dbl(a.mantissa());
}



//
// functional predicates
//

inline bool
is_char (const base_bigmod & a)
{
	return a.is_char();
}



inline bool
is_uchar (const base_bigmod & a)
{
	return a.is_uchar();
}



inline bool
is_short (const base_bigmod & a)
{
	return a.is_short();
}



inline bool
is_ushort (const base_bigmod & a)
{
	return a.is_ushort();
}



inline bool
is_int (const base_bigmod & a)
{
	return a.is_int();
}



inline bool
is_uint (const base_bigmod & a)
{
	return a.is_uint();
}



inline bool
is_long (const base_bigmod & a)
{
	return a.is_long();
}



inline bool
is_ulong (const base_bigmod & a)
{
	return a.is_ulong();
}



//
// arithmetic procedures
//

inline void
negate(base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	a.assign(b);
	a.negate(m);
}



inline void
add (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	add (c.I, a.I, b.I);
	if (c.I.compare(m) >= 0)
		subtract (c.I, c.I, m);
}



inline void
add (base_bigmod & c, const base_bigmod & a, long b, const bigint & m)
{
	add (c.I, a.I, b);
	c.normalize (m);
}



inline void
add (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m)
{
	add (c.I, a.I, b);
	c.normalize (m);
}



inline void
add (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m)
{
	add (c.I, a.I, b);
	c.normalize (m);
}



inline void
subtract (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	subtract (c.I, a.I, b.I);
	if (c.I.is_negative())
		add (c.I, c.I, m);
}



inline void
subtract (base_bigmod & c, const base_bigmod & a, long b, const bigint & m)
{
	subtract (c.I, a.I, b);
	c.normalize (m);
}



inline void
subtract (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m)
{
	subtract (c.I, a.I, b);
	c.normalize (m);
}



inline void
subtract (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m)
{
	subtract (c.I, a.I, b);
	c.normalize (m);
}



inline void
multiply (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	multiply (c.I, a.I, b.I);
	remainder (c.I, c.I, m);
}



inline void
multiply (base_bigmod & c, const base_bigmod & a, long b, const bigint & m)
{
	multiply (c.I, a.I, b);
	c.normalize (m);
}



inline void
multiply (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m)
{
	multiply (c.I, a.I, b);
	c.normalize (m);
}



inline void
multiply (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m)
{
	multiply (c.I, a.I, b);
	c.normalize (m);
}



inline void
square(base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	square(a.I, b.I);
	remainder(a.I, a.I, m);
}



inline void
square(base_bigmod & a, const bigint & b, const bigint & m)
{
	square(a.I, b);
	remainder(a.I, a.I, m);
}



inline void
invert(base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	a.assign(b);
	a.invert(m);
}



void divide (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m);
void divide (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);



inline void
divide (base_bigmod & c, const base_bigmod & a, long b, const bigint & m)
{
	divide(c, a, base_bigmod(b), m);
}



inline void
divide (base_bigmod & c, const base_bigmod & a, unsigned long b, const bigint & m)
{
	divide(c, a, base_bigmod(b), m);
}



inline base_bigmod
inverse(const base_bigmod & a, const bigint & m)
{
	base_bigmod c(a);

	c.invert(m);
	return c;
}



void power (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m);
void power (base_bigmod & c, const base_bigmod & a, long b, const bigint & m);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BASE_BIGMOD_H_GUARD_
