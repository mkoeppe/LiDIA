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



const char *matrix_error_msg[] =
{
	"(0) Parameters have to be greater than zero !",
	"(1) Error found in parameters !",
	"(2) Unknown matrix format !!",
	"(3) Parameters out of range !!",
	"(4) Error in dimensions of the matrices !",
	"(5) Error in format !",
	"(6) No memory allocated for parameters !",
	"(7) Non square matrix !",
	"(8) No solution !!",
	"(9) File LIDIA_PRIMES not found. Please check path and filename.",
	"(10) Number of rows is not equal rank of matrix !",
	"(11) Matrix is not regular !"
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
