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


#include	"LiDIA/pair.h"
#include	"LiDIA/pair.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



template class pair< TYPE1, TYPE2 >;

template int compare(pair< TYPE1, TYPE2 > const &,
		     pair< TYPE1, TYPE2 > const &);




#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
