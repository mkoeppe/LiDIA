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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/LiDIA.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifndef LIDIA_POINTER
# define LIDIA_POINTER

const char *PRT = "(pointer)";

#endif



#ifndef LIDIA_VECTOR_ERROR_MSG_CC_GUARD_
#define LIDIA_VECTOR_ERROR_MSG_CC_GUARD_



const char *vector_error_msg[] =
{
	"(0) Parameter has to be greater than zero !",
	"(1) Error in parameter !",
	"(2) Wrong vector format !!",
	"(3) Parameter out of range !!",
	"(4) Error in dimensions of the matrices !",
	"(5) Error in format !",
	"(6) Number of allocated elements has to be greater than zero.",
	"(7) Vector cannot expanded. Please check mode !!",
	"(8) Ratio has to be greater than one.",
	"(9) Vectormode has to be EXPAND or FIXED.",
	"(10) Invalid sort direction !",
	"(11) Invalid cut !",
	"(12) Invalid position or number of elements !",
	"(13) Invalid cut or empty vector !",
	"(14) Different sizes !"
};



#endif	// LIDIA_VECTOR_ERROR_MSG_CC_GUARD_



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
