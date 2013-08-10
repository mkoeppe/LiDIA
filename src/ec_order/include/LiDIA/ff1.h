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
//	Author	: Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


// The class ff1 realizes single precision prime field arithmetic.
// Surely this should be changed in future.


#ifndef LIDIA_FF1_H_GUARD_
#define LIDIA_FF1_H_GUARD_



#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class ff1
{
private:

	udigit x;
	static udigit p;

public:

	ff1 ();
	ff1 (const ff1 &);
	ff1 (const udigit &);
	~ff1 ();

	static void  set_characteristic (udigit);
	static udigit characteristic ()
	{
		return ff1::p;
	}

	ff1 & operator = (const ff1 &);
	bool operator== (const ff1 &) const;
	void assign_one ();
	void assign_zero ();
	udigit as_udigit()
	{
		return x;
	}
	void assign_non_square();
	bool is_one () const;
	bool is_zero () const;
	bool is_square () const;
	udigit multiplicative_order () const;
	void assign (udigit a);

	friend std::istream & operator >> (std::istream &, ff1 &);
	friend std::ostream & operator << (std::ostream &, const ff1 &);

	friend void swap(ff1 &, ff1 &);
	friend void add (ff1 &, const ff1 &, const ff1 &);
	friend void subtract (ff1 &, const ff1 &, const ff1 &);
	friend void multiply (ff1 &, const ff1 &, const ff1 &);
	friend void divide (ff1 &, const ff1 &, const ff1 &);
	friend void invert (ff1 &, const ff1 &);
	friend void negate (ff1 &, const ff1 &);
	friend void power (ff1 &, const ff1 &, long);
	friend void sqrt (ff1 &, const ff1 &);
	friend void nearly_all_of_order (udigit, base_vector< ff1 >&);
};

// friend functions of ff1

std::istream & operator >> (std::istream &, ff1 &);
std::ostream & operator << (std::ostream &, const ff1 &);

void swap(ff1 &, ff1 &);
void add (ff1 &, const ff1 &, const ff1 &);
void subtract (ff1 &, const ff1 &, const ff1 &);
void multiply (ff1 &, const ff1 &, const ff1 &);
void divide (ff1 &, const ff1 &, const ff1 &);
void invert (ff1 &, const ff1 &);
void negate (ff1 &, const ff1 &);
void power (ff1 &, const ff1 &, long);
void sqrt (ff1 &, const ff1 &);
void nearly_all_of_order (udigit, base_vector< ff1 >&);



inline void square (ff1 & res, const ff1 & a)
{
	multiply(res, a, a);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FF1_H_GUARD_
