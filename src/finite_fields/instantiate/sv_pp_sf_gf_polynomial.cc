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
//	Author	: Safuat Hamdy (SH)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/single_factor.h"
#include	"LiDIA/base/ppair.h"


#define TYPE ppair< single_factor< gf_polynomial >, lidia_size_t >

#define SORT_VECTOR


#include	"LiDIA/instantiate/vector.cc"
