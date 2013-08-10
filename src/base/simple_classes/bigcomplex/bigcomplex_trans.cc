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


#include	"LiDIA/bigcomplex.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
sin (bigcomplex & c, const bigcomplex & x)
{
	bigfloat tmp;

	sin(c.re, x.re);
	cosh(tmp, x.im);
	multiply(c.re, c.re, tmp);
	cos(c.im, x.re);
	sinh(tmp, x.im);
	multiply(c.im, c.im, tmp);
}



void
cos (bigcomplex & c, const bigcomplex & x)
{
	bigfloat tmp;

	cos(c.re, x.re);
	cosh(tmp, x.im);
	multiply(c.re, c.re, tmp);
	sin(c.im, x.re);
	sinh(tmp, x.im);
	multiply(c.im, c.im, tmp);
	c.im.negate();
}



void
sinh (bigcomplex & c, const bigcomplex & x)
{
	bigfloat tmp;

	cos(c.re, x.im);
	sinh(tmp, x.re);
	multiply(c.re, c.re, tmp);
	sin(c.im, x.im);
	cosh(tmp, x.re);
	multiply(c.im, c.im, tmp);
}



void
cosh (bigcomplex & c, const bigcomplex & x)
{
	bigfloat tmp;

	cos(c.re, x.im);
	cosh(tmp, x.re);
	multiply(c.re, c.re, tmp);
	sin(c.im, x.im);
	sinh(tmp, x.re);
	multiply(c.im, c.im, tmp);
}



void
exp (bigcomplex & c, const bigcomplex & x)
{
	bigfloat tmp;

	exp(tmp, x.re);
	cos(c.re, x.im);
	multiply(c.re, c.re, tmp);
	sin(c.im, x.im);
	multiply(c.im, c.im, tmp);
}



void
log (bigcomplex & c, const bigcomplex & x)
{
	bigfloat h(hypot(x));

	if (h.is_le_zero()) {
		lidia_error_handler("bigcomplex", "log of number less equal zero");
		c.assign_zero();
		return;
	}
	log(c.re, h);
	atan2(c.im, x.im, x.re);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
