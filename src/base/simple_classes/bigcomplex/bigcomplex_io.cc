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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigcomplex.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Input/Output
//

std::istream&
operator >> (std::istream& s, bigcomplex& x)
{
	bigfloat r(0), i(0);
	char ch;

	s.get(ch);
	if (ch == '(') {
		s >> r;
		s.get(ch);
		if (ch == ',') {
			s >> i;
			s.get(ch);
		}
		else
			i = 0;
		if (ch != ')')
			s.clear(std::ios::failbit);
	}
	else {
		s.putback(ch);
		s >> r;
		i = 0;
	}
	x = bigcomplex(r, i);
	return s;
}



std::ostream &
operator << (std::ostream& s, const bigcomplex& x)
{
	return s << "(" << x.real() << "," << x.imag() << ")";
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
