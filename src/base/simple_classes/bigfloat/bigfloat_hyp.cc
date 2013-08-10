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
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/bigfloat_config.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
acosh (bigfloat & y, const bigfloat & x)
{
	bigfloat p1;

	square(p1, x);
	dec(p1);
	if (p1.is_lt_zero())
		lidia_error_handler("bigfloat", "acosh::cannot handle arguments <= 1");
	sqrt(p1, p1);
	add(y, x, p1);
	log(y, y);
}



void
acoth (bigfloat & y, const bigfloat & x)
{
	bigfloat p1(x);

	dec(p1);
	divide(y, 2UL, p1);
	inc(p1);
	if (p1.is_lt_zero())
		lidia_error_handler("bigfloat", "acoth::cannot handle arguments >= -1 and <= 1");
	log(y, y);
	y.e--;
}



void
asinh (bigfloat & y, const bigfloat & x)
{
	bigfloat p1;

	square(p1, x);
	inc(p1);
	sqrt(p1, p1);
	add(y, x, p1);
	log(y, y);
}



void
atanh (bigfloat & y, const bigfloat & x)
{
	bigfloat p1(x);

	p1.negate();
	inc(p1);
	divide(y, 2UL, p1);
	dec(y);
	if (p1.is_lt_zero())
		lidia_error_handler("bigfloat", "atanh::cannot handle arguments <= -1 and >= 1");
	log(y, y);
	y.e--;
}



void
cosh (bigfloat & y, const bigfloat & x)
{
	bigfloat tmp(x);

	exp(y, tmp);
	divide(tmp, 1UL, y);
	add(y, y, tmp);
	y.e--;
}



void
coth (bigfloat & y, const bigfloat & x)
{
	bigfloat tmp1, tmp2;

	exp(y, x);
	divide(tmp1, 1UL, y);
	subtract(tmp2, y, tmp1);
	add(y, y, tmp1);
	divide(y, y, tmp2);
}



void
sinh (bigfloat & y, const bigfloat & x)
{
	bigfloat tmp(x);

	exp(y, tmp);
	divide(tmp, 1UL, y);
	subtract(y, y, tmp);
	y.e--;
}



void
tanh (bigfloat & y, const bigfloat & x)
{
	bigfloat tmp1, tmp2;

	exp(y, x);
	divide(tmp1, 1UL, y);
	add(tmp2, y, tmp1);
	subtract(y, y, tmp1);
	divide(y, y, tmp2);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
