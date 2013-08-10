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


// The class trace_mod is responsible for handling the possibilities
// for (c mod l) for some prime l.


#ifndef LIDIA_TRACE_MOD_H_GUARD_
#define LIDIA_TRACE_MOD_H_GUARD_



#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_LIDIA_VECTOR_H_GUARD_
# include	"LiDIA/lidia_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class trace_mod
{
public:
	friend int compare(const trace_mod & tm1, const trace_mod & tm2);

	friend std::ostream & operator << (std::ostream &, const trace_mod &);
	friend std::istream & operator >> (std::istream &, trace_mod &);

private:

	udigit    modulus; 	// modulus l
	base_vector< udigit > cpos; // possible traces modulo l
	float  inv_density; // [ l / cpos.size() ]

public:

	trace_mod ();
	~trace_mod ();

	trace_mod & operator = (const trace_mod &);
	void set_first_element (udigit p_l, udigit c);
	void set_vector (udigit p_l, const base_vector<udigit> & p_cpos);

	lidia_size_t get_size_of_trace_mod () const
	{
		return cpos.size();
	}

	udigit get_modulus () const
	{
		return modulus;
	}

	udigit get_first_element () const
	{
		return cpos[0];
	}

	udigit get_element(lidia_size_t num_of_element) const
	{
		return (cpos[num_of_element]);
	}

	friend void swap(trace_mod & t1, trace_mod & t2);
};

// friend functions of class trace_mod

int compare(const trace_mod & tm1, const trace_mod & tm2);

std::ostream & operator << (std::ostream &, const trace_mod &);
std::istream & operator >> (std::istream &, trace_mod &);

void swap(trace_mod & t1, trace_mod & t2);



inline int compare(const trace_mod & tm1, const trace_mod & tm2)
{
	if (&tm1 == &tm2)
		return 0;
	if (tm1.inv_density == tm2.inv_density)
		return 0;
	else
		return (tm1.inv_density < tm2.inv_density ? -1 : 1);
}



inline bool operator < (const trace_mod & tm1, const trace_mod & tm2)
{
	return (compare (tm1, tm2) == -1);
}



inline bool operator <= (const trace_mod & tm1, const trace_mod & tm2)
{
	return (compare (tm1, tm2) <= 0);
}



inline bool operator == (const trace_mod & tm1, const trace_mod & tm2)
{
	return (compare (tm1, tm2) == 0);
}



inline bool operator > (const trace_mod & tm1, const trace_mod & tm2)
{
	return (compare (tm1, tm2) == 1);
}



inline bool operator >= (const trace_mod & tm1, const trace_mod & tm2)
{
	return (compare (tm1, tm2) >= 0);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_TRACE_MOD_H_GUARD_
