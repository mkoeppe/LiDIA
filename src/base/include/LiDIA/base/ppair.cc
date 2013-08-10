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
//	Author	: Thomas Papanikolaou (TP), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_PPAIR_CC_GUARD_
#define LIDIA_PPAIR_CC_GUARD_



#ifndef LIDIA_PPAIR_H_GUARD_
# include	"LiDIA/base/ppair.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// member functions
//

template< class T1, class T2 >
void
ppair< T1, T2 >::read (std::istream & in)
{
	debug_handler("pair", "read()");
	char c = 0;

	in >> c;
	if (c != '(')
		lidia_error_handler("pair", "read()::(expected");
	else {
		in >> *l >> c;
		if (c != ',')
			lidia_error_handler("pair", "read()::, expected");
		in >> r >> c;
		while (c != ')')
			in >> c;
	}
}



template< class T1, class T2 >
void
ppair< T1, T2 >::write (std::ostream & out) const
{
	debug_handler("pair", "print()");
	out << "(" << *l << " , " << r << ")";
}



template< class T1, class T2 >
int
ppair< T1, T2 >::compare (const ppair< T1, T2 > & p) const
{
	int res = 0;

	if (*l > *p.l)
		res = 1;
	if (*l < *p.l)
		res = -1;

	if (res == 0) {
		if (r > p.r)
			res = 1;
		if (r < p.r)
			res = -1;
	}
	return res;
}



template< class T1, class T2 >
void
ppair< T1, T2 >::swap (ppair< T1, T2 > & p)
{
	T1 *tmp_1 = p.l;
	p.l = l;
	l = tmp_1;

	LiDIA::swap(r, p.r);
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_PPAIR_CC_GUARD_
