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
//	Author	: Thorsten Rottschaefer (TR)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_UDIGIT_MOD_H_GUARD_
#define LIDIA_UDIGIT_MOD_H_GUARD_



#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class udigit_mod
{
private:

        static udigit prime;
	static udigit h;
	udigit u;



public:

        //
        // modulus handling
        //

        static void set_modulus(udigit modulus);
        static udigit get_modulus();



        //
        // c'tors and d'tor
        //

        udigit_mod();
	udigit_mod(const udigit_mod & n);
        udigit_mod(udigit n);
        ~udigit_mod();




        //
        // accessors
        //

        udigit get_mantissa() const;



	//
	// assigners
	//

        void assign_zero();
        void assign_one();

	void assign(const udigit_mod & n);
	void assign(udigit n);

        udigit_mod & operator = (const udigit_mod & n);
        udigit_mod & operator = (udigit n);



	//
	// comparators
	//

        bool is_zero() const;
        bool is_one() const;
	bool is_equal(const udigit_mod & b) const;



	//
	// modifiers
	//

	void inc();
	void dec();

        void swap(udigit_mod & b);



        //
        // friend functions
        //

        friend void negate(udigit_mod & a, const udigit_mod & b);

        friend void add (udigit_mod & c, const udigit_mod & a, const udigit_mod & b);
        friend void add (udigit_mod & c, const udigit_mod & a, udigit b);

        friend void subtract (udigit_mod & c, const udigit_mod & a, const udigit_mod & b);
        friend void subtract (udigit_mod & c, const udigit_mod & a, udigit b);
        friend void subtract (udigit_mod & c, udigit a, const udigit_mod & b);

	friend void multiply (udigit_mod & c, const udigit_mod & a, const udigit_mod & b);
	friend void multiply (udigit_mod & c, const udigit_mod & a, udigit b);
        friend void square (udigit_mod & c, const udigit_mod & a);

        friend void divide (udigit_mod & c, const udigit_mod & a, const udigit_mod & b);
        friend void divide (udigit_mod & c, const udigit_mod & a, udigit b);
        friend void divide (udigit_mod & c, udigit a, const udigit_mod & b);
        friend void invert(udigit_mod & a, const udigit_mod & b);



	//
	// I/O
	//

	void read(std::istream & in);
	void write(std::ostream & out) const;

};



inline void
udigit_mod::set_modulus(udigit modulus)
{
	udigit_mod::prime = modulus;
}



inline udigit
udigit_mod::get_modulus()
{
	return prime;
}



inline
udigit_mod::udigit_mod()
	: u(0UL)
{
}



inline
udigit_mod::udigit_mod (const udigit_mod & n)
	: u(n.u)
{
}



inline
udigit_mod::udigit_mod (udigit n)
{
	u = divide(h, 0, n, prime);
}



inline
udigit_mod::~udigit_mod ()
{
}



//
// accessors
//

inline udigit
udigit_mod::get_mantissa() const
{
	return u;
}



//
// assigners
//

inline void
udigit_mod::assign_zero()
{
	u = 0UL;
}



inline void
udigit_mod::assign_one()
{
	u = 1UL;
}



inline void
udigit_mod::assign (udigit n)
{
	u = divide(h, 0, n, prime);
}



inline void
udigit_mod::assign (const udigit_mod & n)
{
	u = divide(h, 0, n.u, prime);
}



inline udigit_mod &
udigit_mod::operator = (udigit n)
{
	assign(n);
	return *this;
}



inline udigit_mod &
udigit_mod::operator = (const udigit_mod & n)
{
	assign(n);
	return *this;
}



//
// comparators
//

inline bool
udigit_mod::is_zero() const
{
	return (u == 0);
}



inline bool
udigit_mod::is_one() const
{
	return (u == 1);
}



inline bool
udigit_mod::is_equal (const udigit_mod & b) const
{
	return (u == b.u);
}



inline bool
operator == (const udigit_mod & a, const udigit_mod & b)
{
	return a.is_equal(b);
}



inline bool
operator != (const udigit_mod & a, const udigit_mod & b)
{
	return !a.is_equal(b);
}



//
// modifiers
//

inline void
udigit_mod::inc ()
{
	u++;
	if (u == prime) {
		u = 0;
	}
}



inline void
udigit_mod::dec ()
{
	u--;
	if (u < 0) {
		u = prime + u;
	}
}



inline void
udigit_mod::swap (udigit_mod & b)
{
	udigit t;

	t = u;
	u = b.u;
	b.u = t;
}



//
// arithmetic procedures
//

inline void
negate (udigit_mod & a, const udigit_mod & b)
{
	if (b.u == 0)
		a.u = 0;
	else
		subtract(a.u, udigit_mod::prime, b.u, 0);
}



inline void
add (udigit_mod & c, const udigit_mod & a, const udigit_mod & b)
{
	c.u = add_mod(a.u, b.u, udigit_mod::prime);
}



inline void
add (udigit_mod & c, const udigit_mod & a, udigit b)
{
	c.u = add_mod(a.u, b, udigit_mod::prime);
}



inline void
subtract (udigit_mod & c, const udigit_mod & a, const udigit_mod & b)
{
	c.u = subtract_mod(a.u, b.u, udigit_mod::prime);
}



inline void
subtract (udigit_mod & c, const udigit_mod & a, udigit b)
{
	c.u = subtract_mod(a.u, b, udigit_mod::prime);
}



inline void
subtract (udigit_mod & c, udigit a, const udigit_mod & b)
{
	c.u = subtract_mod(a, b.u, udigit_mod::prime);
}



inline void
multiply (udigit_mod & c, const udigit_mod & a, const udigit_mod & b)
{
	c.u = multiply_mod(a.u, b.u, udigit_mod::prime);
}



inline void
multiply (udigit_mod & c, const udigit_mod & a, udigit b)
{
	c.u = multiply_mod(a.u, b, udigit_mod::prime);
}



inline void
square (udigit_mod & c, const udigit_mod & a)
{
	c.u = multiply_mod(a.u, a.u, udigit_mod::prime);
}



inline void
divide (udigit_mod & c, const udigit_mod & a, const udigit_mod & b)
{
	c.u = divide_mod(a.u, b.u, udigit_mod::prime);
}



inline void
divide (udigit_mod & c, const udigit_mod & a, udigit b)
{
	c.u = divide_mod(a.u, b, udigit_mod::prime);
}



inline void
divide (udigit_mod & c, udigit a, const udigit_mod & b)
{
	c.u = divide_mod(a, b.u, udigit_mod::prime);
}



inline void
invert(udigit_mod & a, const udigit_mod & b)
{
	a.u = invert_mod(b.u, udigit_mod::prime);
}



inline void
swap(udigit_mod & a, udigit_mod & b)
{
	a.swap(b);
}



//
// arithmetic operators
//

inline udigit_mod
operator - (const udigit_mod & a)
{
	udigit_mod c;

	negate(c, a);
	return c;
}



inline udigit_mod
operator + (const udigit_mod & a, const udigit_mod & b)
{
	udigit_mod c;

	add(c, a, b);
	return c;
}



inline udigit_mod
operator + (const udigit_mod & a, udigit b)
{
	udigit_mod c;

	add(c, a, b);
	return c;
}



inline udigit_mod
operator + (udigit a, const udigit_mod & b)
{
	udigit_mod c;

	add(c, b, a);
	return c;
}



inline udigit_mod
operator - (const udigit_mod & a, const udigit_mod & b)
{
	udigit_mod c;

	subtract(c, a, b);
	return c;
}



inline udigit_mod
operator - (const udigit_mod & a, udigit b)
{
	udigit_mod c;

	subtract(c, a, b);
	return c;
}



inline udigit_mod
operator - (udigit a, const udigit_mod & b)
{
	udigit_mod c;

	subtract(c, a, b);
	return c;
}



inline udigit_mod
operator * (const udigit_mod & a, const udigit_mod & b)
{
	udigit_mod c;

	multiply(c, a, b);
	return c;
}



inline udigit_mod
operator * (const udigit_mod & a, udigit b)
{
	udigit_mod c;

	multiply(c, a, b);
	return c;
}



inline udigit_mod
operator * (udigit a, const udigit_mod & b)
{
	udigit_mod c;

	multiply(c, b, a);
	return c;
}



inline udigit_mod
operator / (const udigit_mod & a, const udigit_mod & b)
{
	udigit_mod c;

	divide(c, a, b);
	return c;
}



inline udigit_mod
operator / (const udigit_mod & a, udigit b)
{
	udigit_mod c;

	divide(c, a, b);
	return c;
}



inline udigit_mod
operator / (udigit a, const udigit_mod & b)
{
	udigit_mod c;

	divide(c, a, b);
	return c;
}



inline udigit_mod &
operator += (udigit_mod & a, const udigit_mod & b)
{
	add(a, a, b);
	return a;
}



inline udigit_mod &
operator += (udigit_mod & a, udigit b)
{
	add(a, a, b);
	return a;
}



inline udigit_mod &
operator -= (udigit_mod & a, const udigit_mod & b)
{
	subtract(a, a, b);
	return a;
}



inline udigit_mod &
operator -= (udigit_mod & a, udigit b)
{
	subtract(a, a, b);
	return a;
}



inline udigit_mod &
operator *= (udigit_mod & a, const udigit_mod & b)
{
	multiply(a, a, b);
	return a;
}



inline udigit_mod &
operator *= (udigit_mod & a, udigit b)
{
	multiply(a, a, b);
	return a;
}



inline udigit_mod &
operator /= (udigit_mod & a, const udigit_mod & b)
{
	divide(a, a, b);
	return a;
}



inline udigit_mod &
operator /= (udigit_mod & a, udigit b)
{
	divide(a, a, b);
	return a;
}



inline udigit_mod &
operator ++ (udigit_mod & a)
{
	a.inc();
	return a;
}



inline udigit_mod
operator ++ (udigit_mod & a, int)
{
	udigit_mod c(a);

	a.inc();
	return c;
}



inline udigit_mod &
operator -- (udigit_mod & a)
{
	a.dec();
	return a;
}



inline udigit_mod
operator -- (udigit_mod & a, int)
{
	udigit_mod c(a);

	a.dec();
	return c;
}



inline bool
operator ! (const udigit_mod & a)
{
	return !a.is_zero();
}



//
// I/O
//

inline void
udigit_mod::read (std::istream & in)
{
	in >> u;
}



inline void
udigit_mod::write (std::ostream & out) const
{
	out << u;
}



inline std::istream &
operator >> (std::istream & in, udigit_mod & a)
{
	a.read(in);
	return in;
}



inline std::ostream &
operator << (std::ostream & out, const udigit_mod & a)
{
	a.write(out);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_UDIGIT_MOD

#include	"LiDIA/specialization/udigit_mod.special"



#endif	// LIDIA_UDIGIT_MOD_H_GUARD_
