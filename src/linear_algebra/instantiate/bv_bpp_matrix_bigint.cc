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


#include	"LiDIA/bigint_matrix.h"
#include	"LiDIA/base_ppair.h"


#define TYPE base_ppair< ring_matrix< bigint >, lidia_size_t >

#define BASE_VECTOR


#include	"LiDIA/instantiate/vector.cc"
