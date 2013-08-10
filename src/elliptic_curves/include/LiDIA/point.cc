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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_POINT_CC_GUARD_
#define LIDIA_POINT_CC_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include	"LiDIA/point.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructors / destructor
//
template< class T >
point< T >
::point() : base_point< T > ()
{
	debug_handler("point< T >",
		      "point()");
	this->info = 0;
}



template< class T >
point< T >
::point(const base_point< T > & bp) : base_point< T > (bp)
{
	debug_handler("point< T >",
		      "point(const base_point< T > &)");
	this->info = 0;
}



template< class T >
point< T >
::point(const T & xp, const T & yp,
	const elliptic_curve< T > & e) : base_point< T > ()
{
	debug_handler("point< T >",
		      "point(constT&, constT&, "
		      "const elliptic_curve< T > &)");

	this->assign(xp, yp, e);
}



template< class T >
point< T >
::point(const T & xp, const T & yp, const T & zp,
	const elliptic_curve< T > & e) : base_point< T > ()
{
	debug_handler("point< T >",
		      "point(constT&, constT&, constT&"
		      "const elliptic_curve< T > &)");

	this->assign(xp, yp, zp, e);
}



template< class T >
point< T >
::point(const point< T > & P) : base_point< T > ()
{
	this->assign(P);
}



template< class T >
point< T >
::point(const elliptic_curve< T > &e) : base_point< T > ()
{
	debug_handler("point< T >",
		      "point(const elliptic_curve< T > &");

	this->assign_zero(e);
}



template< class T >
point< T >
::~point()
{
	debug_handler("point< T >",
		      "~point()");
}



//
// assignment
//

template< class T >
point< T > & point< T >
::operator = (const point< T > & P)
{
	this->assign(P);
	return *this;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_POINT_CC_GUARD_
