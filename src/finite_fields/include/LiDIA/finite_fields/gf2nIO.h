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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_GF2NIO_H_GUARD_
#define LIDIA_GF2NIO_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

class gf2n;

class gf2nIO
{
public:

	enum base {
		Dec,
		Hex
	};

private:

	static base iobase; // actual base
	static char *ioprefix; // prefix

	friend class gf2n;
	friend std::ostream & operator << (std::ostream &, const gf2n &);

	// ----------------------------------------------------------------------

public:

	gf2nIO();
	gf2nIO(base, char *);

	static void           setbase      (base new_base);
	static base           showbase     ();
	static void           setprefix    (base actual_base);
	static void           setprefix    (char *prefix);
	static void           noprefix     ();
	static char*          showprefix   ();

	// ----------------------------------------------------------------------
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GF2NIO_H_GUARD_
