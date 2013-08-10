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


//
//  usage:
//
//	include all headers that are necessary to describe the type TYPE
//
//	define the type TYPE
//
//	include this file
//


#include	"LiDIA/power_functions.h"
#include	"LiDIA/power_functions.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef EXP_LONG

template void lidia_power_left_to_right (TYPE &, const TYPE &, unsigned long);
template void lidia_power_left_to_right (TYPE &, const TYPE &, unsigned long,
					 void (*fast_mul) (TYPE &, const TYPE &, const TYPE &));
template void lidia_power_right_to_left (TYPE &, const TYPE &, unsigned long);

#endif


#ifdef EXP_BIGINT

template void lidia_power_left_to_right (TYPE &, const TYPE &, bigint);
template void lidia_power_left_to_right (TYPE &, const TYPE &, bigint,
					 void (*fast_mul) (TYPE &, const TYPE &, const TYPE &));
template void lidia_power_right_to_left (TYPE &, const TYPE &, bigint);

#endif


#ifdef SLIDING_WINDOW

template void lidia_power_sliding_window (TYPE &, const TYPE &, const bigint &, unsigned long);

#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
