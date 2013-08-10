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


// =================== hash_table ===================

#if defined (HASH_TABLE)

# include	"LiDIA/hash_table.h"
# include	"LiDIA/hash_table.cc"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class hentry< TYPE >;
template class hash_table< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// HASH_TABLE



// =================== indexed_hash_table ===================

#if defined (I_HASH_TABLE)

# include	"LiDIA/indexed_hash_table.h"
# include	"LiDIA/indexed_hash_table.cc"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class indexed_hash_table< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// I_HASH_TABLE
