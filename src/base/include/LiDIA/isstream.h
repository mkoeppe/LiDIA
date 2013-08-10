// -*-C++-*-

//=============================================================================
//
//     This file is part of LiDIA --- a library for computational number theory
//
//     Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//     See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//-----------------------------------------------------------------------------
//
//     $Id$
//
//     Author  : Christoph Ludwig
//     Changes : See CVS log
//
//=============================================================================

/******************************************************************************
**
** Rational:
** =========
**
** The typedef and the function declared in this header are intended for use
** within LiDIA. User applications that are targeted at one platform only do
** not need to use this typedef / function and adhere the restrictions below.
**
** Unfortunately, some older standard library implementations (notably
** the one shipping with g++ 2.95) lack support for the new string streams
** declared in the standard header sstream.
**
** But we want to discourage explicit usage of the deprecated
** std::istrstream. We therefore create a class isstream that is derived from
** the new std::istringstream if available and from std::strstream else.
**
** You must not call any method of the string stream besides the usual
** extractors. Otherwise LiDIA may show different behaviour on diffent
** platforms (or LiDIA may not be buildable at all on some platforms)!
**
******************************************************************************/

#if !defined(LIDIA_ISSTREAM_H_GUARD_)
#define LIDIA_ISSTREAM_H_GUARD_

#include "LiDIA/LiDIA.h"


#if !defined(OSSTREAM_USE_STRSTREAM)
# include <sstream>   // std::istringstream
#else
# include <strstream> // std::strstream
#endif

# include <string>    // std::string


#if defined(LIDIA_NAMESPACE)
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

#if !defined(OSSTREAM_USE_STRSTREAM)

  class isstream : public std::istringstream {
  public:

      explicit isstream(const std::string& data) : std::istringstream(data) {
	  // Nothing to do
      }
  };

#else
  class isstream : public std::strstream {
  public:
      explicit isstream(const std::string& data) {
	  // Note: We can't initialize strstream with data.c_str(), because
	  // the strstream-buffer takes ownership of the array the
	  // passed pointer refers to.
	  *this << data << std::ends;
      }
  };
#endif	  


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif
