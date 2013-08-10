// -*- C++ -*-
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


#ifndef LIDIA_VECTOR_REPRESENTATION_H_GUARD_
#define LIDIA_VECTOR_REPRESENTATION_H_GUARD_



#ifndef LIDIA_VECTOR_FLAGS_H_GUARD_
# include	"LiDIA/vector_flags.h"
#endif
#ifndef LIDIA_LIDIA_DEFINES_H_GUARD_
# include	"LiDIA/lidia_defines.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class vector_representation
{

protected:

	T *value;

	lidia_size_t *index_array;

	lidia_size_t length;
	lidia_size_t allocated;

	vector_flags bitfield;

	char sort_dir; // sort-direction
	int (*el_cmp) (const T &, const T &); // compare-func. for vec.elements

	//
	// friend classes
	//

	friend class crt;

	vector_representation()
	{
		debug_handler("vector_representation",
			      "vector_representation()");
	}

	~vector_representation()
	{
		debug_handler("vector_representation",
			      "~vector_representation()");
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_VECTOR_REPRESENTATION_H_GUARD_
