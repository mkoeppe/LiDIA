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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


// Description  : This file describes an implementation for
//                computing a table of powers, computed with
//                squaring trick.
//                We assume that the template type T supports functions
//                multiply, square and assign_one.


#ifndef LIDIA_POWER_TABLE_H_GUARD_
#define LIDIA_POWER_TABLE_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class power_tab
{

	//
	// this class represents powers of an instance of
	// class T produced by squaring_trick
	//

private:

	bool initialized;
	T * f_powers;
	lidia_size_t N;

public:

	//
	// constructors and destructor
	//

	power_tab();
	~power_tab();

	//
	// member functions
	//

	void initialize(const T & f, lidia_size_t n);
	void get_power(T &f, lidia_size_t i) const;
	T get_power(lidia_size_t i) const;
};



#define LIDIA_CLASS_POWER_TAB



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_POWER_TABLE_H_GUARD_
