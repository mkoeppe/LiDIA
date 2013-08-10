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


#ifndef LIDIA_BASE_PPAIR_CC_GUARD_
#define LIDIA_BASE_PPAIR_CC_GUARD_


#ifndef LIDIA_BASE_PPAIR_H_GUARD_
# include	"LiDIA/base_ppair.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// swap function
//

template< class T1, class T2 >
void
base_ppair< T1, T2 >::swap (base_ppair< T1, T2 > & q)
{
	T1 *tmp_1 = l;
	l = q.l;
	q.l = tmp_1;

	LiDIA::swap(r, q.r);
}



template< class T1, class T2 >
void
base_ppair< T1, T2 >::read (std::istream & in)
{
	debug_handler("base_ppair", "read(std::istream &)");

	char c = 0;

	in >> c;
	if (c != '(')
		lidia_error_handler("base_ppair", "read()::(expected");
	else {
		in >> *l >> c;
		if (c != ',')
			lidia_error_handler("base_ppair", "read()::, expected");
		in >> r >> c;
		while (c != ')')
			in >> c;
	}
}



template< class T1, class T2 >
void
base_ppair< T1, T2 >::print (std::ostream & out) const
{
	debug_handler("base_ppair", "print(std::ostream &)");
	out << "(" << *l << " , " << r << ")";
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_PPAIR_CC_GUARD_
