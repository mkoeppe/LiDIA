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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_PAIR_CC_GUARD_
#define LIDIA_PAIR_CC_GUARD_



#ifndef LIDIA_PAIR_H_GUARD_
# include	"LiDIA/pair.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



template< class T1, class T2 >
void
pair< T1, T2 >::read (std::istream & in)
{
	debug_handler("pair", "read()");
	char c = 0;

	in >> c;
	if (c != '(')
		lidia_error_handler("pair", "read()::(expected");
	else {
		in >> l >> c;
		if (c != ',')
			lidia_error_handler("pair", "read()::, expected");
		in >> r >> c;
		while (c != ')')
			in >> c;
	}
}



template< class T1, class T2 >
void
pair< T1, T2 >::write (std::ostream & out) const
{
	debug_handler("pair", "print()");
	out << "(" << l << " , " << r << ")";
}



template< class T1, class T2 >
int
compare (const pair< T1, T2 > & p, const pair< T1, T2 > & q)
{
	int res = 0;

	if (&p == &q)
		return 0;

	if (p.left() > q.left())
		res = 1;
	if (p.left() < q.left())
		res = -1;

	if (res == 0) {
		if (p.right() > q.right())
			res = 1;
		if (p.right() < q.right())
			res = -1;
	}
	return res;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_PAIR_CC_GUARD_
