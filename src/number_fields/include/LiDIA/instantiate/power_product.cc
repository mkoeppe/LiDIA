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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


//
//  usage:
//
//	include all headers that are necessary to describe the types TYPE1 and TYPE2
//
//	define the types TYPE1 and TYPE2
//
//	include this file
//


#include	"LiDIA/base_power_product.h"
#include	"LiDIA/base_power_product.cc"
#include	"LiDIA/base_ppair.h"
#include	"LiDIA/base_ppair.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// base_power_product
//

template class base_power_product< TYPE1 , TYPE2 >;



//
// base_pair
//

template class base_ppair< TYPE1 , TYPE2 >;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
