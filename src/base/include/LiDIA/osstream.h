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
** std::ostrstream. We therefore create a typedef osstream that refers to
** the new std::ostringstream if available and to std::ostrstream else.
**
** The two string stream classes differ in who takes ownership of the
** string buffer. Since we do not want to deal with this difference
** via #ifdef's cluttered all over the place, osstream.h declares a
** function for extracting the string written to the stream. It makes sure
** the stream stays responsible for its buffer.
**
** You must not call any method of the string stream besides the usual
** inserters. Otherwise LiDIA may show different behaviour on diffent
** platforms (or LiDIA may not be buildable at all on some platforms)!
**
******************************************************************************/

#if !defined(LIDIA_OSSTREAM_H_GUARD_)
#define LIDIA_OSSTREAM_H_GUARD_

#include "LiDIA/LiDIA.h"


#if !defined(OSSTREAM_USE_STRSTREAM)
# include <sstream>   // std::ostringstream
#else
# include <strstream> // std::ostrstream
#endif

# include <string>    // std::string


#if defined(LIDIA_NAMESPACE)
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

#if !defined(OSSTREAM_USE_STRSTREAM)

  typedef std::ostringstream osstream;

#else

  typedef std::ostrstream osstream;

#endif

  
  // Get a copy of the output string stream's buffer in a std::string.
  // If osstream happens to be std::ostrstream then we return
  // ownership of the buffer to the stream.
  //
  // We cannot declare oss to be a const reference since extracting
  // the string from std::ostrstream is not a const operation.
  //
  // CAVEAT: std::ostrstream offers no way to determine the number
  // of characters in the output buffer. extractString must therefore
  // insert '\0' at the current write position. The original write
  // position is restored afterwards, but if you moved the write position
  // previous to the call of extractString() you may get different results
  // depending on whether your system has sstream or not.
  // Therefore, calls to osstream::seekp() must be avoided!
  std::string extractString(osstream& oss);

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif
