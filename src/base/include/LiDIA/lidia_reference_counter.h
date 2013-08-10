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
//	$Id: lidia_reference_counter.h,v 2.3 2002/06/24 09:41:29 lidiaadm Exp $
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_REFERENCE_COUNTER_H_GUARD_
#define LIDIA_REFERENCE_COUNTER_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class lidia_reference_counter
{
	int rc;
public:
	lidia_reference_counter()
	{
		rc = 1;
	}

	void inc_ref_counter()
	{
		rc++;
	}

	void dec_ref_counter()
	{
		rc--;
	}

	int get_ref_counter() const
	{
		return rc;
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_REFERENCE_COUNTER_H_GUARD_
