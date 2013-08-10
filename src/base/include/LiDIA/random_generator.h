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


#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
#define LIDIA_RANDOM_GENERATOR_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class random_generator
{
private:

	static bool initialized;


public:

	random_generator();
	~random_generator();


private:

	// inhibit:
	random_generator(const random_generator&);
	random_generator& operator = (const random_generator&);


public:

	static void seed ();
	static void seed (unsigned int);

};



//
// constructor / destructor
//

inline
random_generator::random_generator ()
{
	if (!initialized) {
		seed();
	}
}



inline
random_generator::~random_generator ()
{
	// nothing to do
}



random_generator & operator >> (random_generator & rg, int  & i);

random_generator & operator >> (random_generator & rg, unsigned int  & i);

random_generator & operator >> (random_generator & rg, long & l);

random_generator & operator >> (random_generator & rg, unsigned long & l);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_RANDOM_GENERATOR_H_GUARD_

