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


#ifndef LIDIA_B_VALUE_H_GUARD_
#define LIDIA_B_VALUE_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



inline
long b_value (long m)
{
	long b = 0;

	if (m < 0) m = -m;

	while (m != 0) {
		b++;
		m >>= 1;
	}
	return b;
}



inline
long b_value (const bigint & m)
{
	return static_cast<long>(m.bit_length());
}



inline
long b_value (const bigfloat & m)
{
	return static_cast<long>(m.exponent() + m.bit_length());
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_B_VALUE_H_GUARD_
