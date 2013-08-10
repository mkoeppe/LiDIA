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
//	Author	: Frank Lehmann (FL), Markus Maurer (MM),
//                Thomas Papanikolaou (TP), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *
// *    File        :  comparator.h
// *
// *    Description :  -) definition of template class 'comparator < T > '
// *			  to provide functions for comparing elements
// *		          of any type T
// *
// *		       -) to use this class comparator, T must support
// *                      the operators < , <= , ==
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


#ifndef LIDIA_COMPARATOR_H_GUARD_
#define LIDIA_COMPARATOR_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif

#include	<cstring>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class comparator
{

public:

	//
	// by reference
	//

	static bool less_than(const T & a, const T & b)
	{
		return (a < b);
	}

	static bool less_equal(const T & a, const T & b)
	{
		return (a <= b);
	}

	static bool greater_than(const T & a, const T & b)
	{
		return (b < a);
	}

	static bool greater_equal(const T & a, const T & b)
	{
		return (b <= a);
	}

	static bool equal(const T & a, const T & b)
	{
		return (b == a);
	}

	//
	// by value
	//

	static bool val_less_than(T a, T b)
	{
		return (a < b);
	}

	static bool val_less_equal(T a, T b)
	{
		return (a <= b);
	}

	static bool val_greater_than(T a, T b)
	{
		return (b < a);
	}

	static bool val_greater_equal (T a, T b)
	{
		return (b <= a);
	}

	static bool val_equal(T a, T b)
	{
		return (b == a);
	}

	//
	// min and max
	//

	static T min(const T & a, const T & b)
	{
		return (a <= b ? a : b);
	}

	static T max(const T & a, const T & b)
	{
		return (b <= a ? a : b);
	}
};



template<>
class comparator< char * >
{
public:

	//
	// by reference
	//

	static bool less_than(char *a, char *b)
	{
		return (strcmp(a, b) < 0);
	}

	static bool less_equal(char *a, char *b)
	{
		return (strcmp(a, b) <= 0);
	}

	static bool greater_than(char *a, char *b)
	{
		return (strcmp(a, b) > 0);
	}

	static bool greater_equal(char *a, char *b)
	{
		return (strcmp(a, b) >= 0);
	}

	static bool equal(char *a, char *b)
	{
		return (strcmp(a, b) == 0);
	}

	//
	// by value
	//

	static bool val_less_than(char *a, char *b)
	{
		return (strcmp(a, b) < 0);
	}

	static bool val_less_equal (char *a, char *b)
	{
		return (strcmp(a, b) <= 0);
	}

	static bool val_greater_than(char *a, char *b)
	{
		return (strcmp(a, b) > 0);
	}

	static bool val_greater_equal(char *a, char *b)
	{
		return (strcmp(a, b) >= 0);
	}

	static bool val_equal(char *a, char *b)
	{
		return (strcmp(a, b) == 0);
	}

	//
	// min and max
	//

	static char* min(char *a, char *b)
	{
		return ((strcmp(a, b) <= 0) ? a : b);
	}

	static char* max(char *a, char *b)
	{
		return ((strcmp(a, b) >= 0) ? a : b);
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_COMPARATOR_H_GUARD_
